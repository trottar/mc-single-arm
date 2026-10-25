"""Audit a notebook without mounting Drive or executing file-dependent cells.

Only reviewed function definitions and one reviewed inline dataset assignment are
compiled from AST. This is not a general notebook runner.
"""
from pathlib import Path
import argparse
import ast
import contextlib
import csv
from datetime import datetime,timezone
import hashlib
import importlib.metadata
import io
import json
import math
import shutil
import sys
import numpy as np
from scipy.optimize import least_squares
from scipy.stats import chi2
from collect_run import ROOT,sha,write_json

def csvwrite(path,rows):
    with path.open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def jsonable(v):
    if isinstance(v,np.ndarray):return v.tolist()
    if isinstance(v,np.generic):return v.item()
    if isinstance(v,dict):return {str(k):jsonable(x) for k,x in v.items()}
    if isinstance(v,(list,tuple)):return [jsonable(x) for x in v]
    return v

def load_selected(cells):
    env={"np":np,"least_squares":least_squares,"chi2":chi2}
    selected=[]
    for number,name in [(14,"fit_data_c1c2"),(15,"fit_data_c1c2_beta"),(6,"read_data_file")]:
        nodes=[n for n in ast.parse(cells[number]).body if isinstance(n,ast.FunctionDef) and n.name==name]
        if len(nodes)!=1:raise ValueError("expected reviewed function absent: "+name)
        selected.extend(nodes)
    # Reviewed assignment contains only dictionary/list literals and np.array.
    assignment=next(n for n in ast.parse(cells[5]).body if isinstance(n,ast.Assign)
        and any(isinstance(t,ast.Name) and t.id=="data_pvdis6" for t in n.targets))
    selected.append(assignment)
    module=ast.fix_missing_locations(ast.Module(body=selected,type_ignores=[]))
    exec(compile(module,"reviewed_notebook_definitions","exec"),env)
    return env

def findings():
    return [
      {"id":"N01","severity":"high","cells":[14,15],"title":"Joint confidence threshold overwritten",
       "evidence":"delta = chi2.ppf(cl, df=2) is immediately followed by delta = 1.0 in both fit functions.",
       "impact":"The displayed 68.3% ellipse encloses 39.35% for a two-dimensional Gaussian. A joint 68.3% contour needs delta=2.2977 and axes 1.5158 times larger.",
       "recommendation":"Expose delta-chi2/coverage explicitly and retain separate one-parameter errors and joint contours."},
      {"id":"N02","severity":"high","cells":[3,8],"title":"Hall C time scaling is inverted and applied to every setting",
       "evidence":"hallc_stat_loq=sqrt(split/0.5), hallc_stat_hiq=sqrt((1-split)/0.5); cell 8 applies loq to all HMS and SHMS files.",
       "impact":"At split=0.2 the notebook multiplies errors by 0.63246 instead of 1.58114 for reduced low-angle exposure relative to the same-total 50:50 baseline. It does not enforce the user's 30-day total budget.",
       "recommendation":"Use sqrt(reference_time/allocated_time) from explicit setting-pair exposure records; reject zero exposure and duplicate time allocations."},
      {"id":"N03","severity":"high","cells":[4,14,15],"title":"Amplitude correlation fraction used as variance fraction",
       "evidence":"solid11_corr=sqrt(0.004^2+0.002^2)/solid11_syst, but fit covariance uses syst_corr=syst*sqrt(f_corr).",
       "impact":"If the listed 0.4% and 0.2% sources are the intended correlated components, f_corr should be their squared sum divided by total systematic variance. Off-diagonal variance is inflated by about 27.5% for the nominal SoLID input.",
       "recommendation":"Define correlated source amplitudes directly or use the squared amplitude fraction consistently; document the intended decomposition."},
      {"id":"N04","severity":"high","cells":[22,28],"title":"Undefined setting groups prevent a clean sequential run",
       "evidence":"data_sets_hms1, data_sets_shms1, data_sets_hms2, data_sets_shms2 are used without assignment; result_11gev_hallc1 and result_11gev_hallc1_all are also used without definition.",
       "impact":"Supplying the missing files alone would not make the exported notebook execute top to bottom.",
       "recommendation":"Define the setting groups from an explicit schedule and remove dependencies on hidden interactive state."},
      {"id":"N05","severity":"high","cells":[18,19,20,23,30],"title":"Display recentering overwrites fitted parameter records",
       "evidence":"Multiple result['params'] assignments replace the fitted centers; the four-parameter result is replaced by a two-element display vector in cell 23.",
       "impact":"Returned params can disagree with fitted, cov_full, chi2, and residuals. Shifted ellipses illustrate precision, not measured best fits.",
       "recommendation":"Preserve fitted_params and covariances; store display_center separately, with the shift and purpose recorded."},
      {"id":"N06","severity":"medium","cells":[14,15],"title":"Common systematics stop at dataset boundaries",
       "evidence":"Covariance is block diagonal between dataset dictionaries, with one same-sign correlated source inside each block.",
       "impact":"Polarimetry or other common sources across HMS/SHMS settings are not represented unless all relevant rows are put into one explicitly shared source block.",
       "recommendation":"Represent each physical systematic as a source vector spanning every affected measurement; validate signs and overlap."},
      {"id":"N07","severity":"medium","cells":[5,15],"title":"Hadronic templates supplied, but application and authority need review",
       "evidence":"A=(c1*X+c2*Y)*(1+beta_HT/((1-x)^3*Q2)+beta_CSV*x^2); beta applies only where the second coefficient is nonzero.",
       "impact":"beta_HT has units GeV^2 and beta_CSV is dimensionless. Both terms multiply the full asymmetry. The coefficient mask is not an explicit dataset-physics label; no priors or bounds are supplied.",
       "recommendation":"Retain this exact ansatz as notebook provenance; validate its scientific interpretation and use an explicit apply_hadronic flag before the PAC fit."},
      {"id":"N08","severity":"medium","cells":[5,18],"title":"Derived external constraints treated as independent priors",
       "evidence":"Elastic+Cs, the Qweak-derived X constraint, and projected P2 X constraint are separate independent scalar datasets with placeholder x=0.99,Q2=10000.",
       "impact":"Independence of the derived Qweak constraint from the elastic+Cs input is not established by the notebook. P2 is a projection, not a measurement.",
       "recommendation":"Recover the upstream likelihood/provenance and shared information; do not claim an independently combined world fit until this is resolved."},
      {"id":"N09","severity":"medium","cells":[6,8,11,12],"title":"External tables missing and parsers silently skip rows",
       "evidence":"28 Hall C tables, four HRS tables, solid11.dat, and pvdis_11gev_CJ15nlo.csv are referenced but not embedded. Parsers ignore wrong column counts and nonnumeric lines.",
       "impact":"Full projections cannot be reproduced from this notebook alone; malformed numerical rows can disappear without a diagnostic.",
       "recommendation":"Supply exact tables and hashes; use explicit schemas, row rejection logs, unit checks, and finite-value/rate validation."},
      {"id":"N10","severity":"low","cells":[8,30,36],"title":"Labeling and legacy contour geometry need cleanup",
       "evidence":"SHMS datasets are named HMS; the Lite legend can use nominal syst; PCA/max-extent approximates historical contour point sets.",
       "impact":"Labels can misidentify datasets, and an enclosing display ellipse is not an inferred likelihood/covariance contour.",
       "recommendation":"Correct labels and retain historical contour points with their original confidence convention and provenance."}
    ]

def audit(source):
    source=Path(source).resolve();original=source.read_bytes();j=json.loads(original)
    if len(j["cells"])!=39:raise ValueError("this reviewed harness expects the supplied 39-cell notebook")
    cells={i:"".join(c.get("source",[])) for i,c in enumerate(j["cells"])}
    for i,c in enumerate(j["cells"]):
        if c["cell_type"]=="code":ast.parse(cells[i])
    notebook_hash=hashlib.sha256(original).hexdigest();harness_hash=sha(Path(__file__))
    reviews=ROOT/"notebook_reviews";reviews.mkdir(exist_ok=True)
    review_id="notebook_"+notebook_hash[:12]+"_audit_"+harness_hash[:8]
    dest=reviews/review_id
    if dest.exists():raise ValueError("review snapshot already exists; use report mode to rebuild the PDF")
    dest.mkdir()
    (dest/"source.ipynb").write_bytes(original)
    shutil.copyfile(__file__,dest/"audit_notebook.py")
    (dest/"cells.txt").write_text("\n\n".join(f"CELL {i} [{j['cells'][i]['cell_type']}]\n{cells[i]}" for i in cells),encoding="utf-8")
    inventory=[]
    for i,c in enumerate(j["cells"]):
        inventory.append({"cell_index_zero_based":i,"cell_type":c["cell_type"],"cell_id":c.get("id",c.get("metadata",{}).get("id")),
            "source_sha256":hashlib.sha256(cells[i].encode()).hexdigest(),"source_lines":len(cells[i].splitlines()),
            "execution_count":c.get("execution_count"),"outputs":[{"type":o.get("output_type"),"mime_types":list(o.get("data",{})),
                "text_chars":len("".join(o.get("text",[])))} for o in c.get("outputs",[])]})
    write_json(dest/"cell_inventory.json",inventory)
    write_json(dest/"findings.json",{"schema_version":"pvdis_notebook_findings/v1","cell_numbering":"zero-based","findings":findings()})
    external=[]
    for arm in ["hms","shms"]:
        for i in range(1,15):external.append({"file":f"hallc_bin_{arm}1{i:02d}.dat","role":arm.upper()+" projection", "status":"not_supplied"})
    external += [{"file":f"halla_hrs{arm}{i}.dat","role":"HRS rate check","status":"not_supplied"} for arm in ["l","r"] for i in [1,2]]
    external += [{"file":n,"role":"SoLID projection","status":"not_supplied"} for n in ["solid11.dat","pvdis_11gev_CJ15nlo.csv"]]
    csvwrite(dest/"external_inputs.csv",external)
    env=load_selected(cells);log=io.StringIO(); probes=[]
    with contextlib.redirect_stdout(log):
        inline=env["data_pvdis6"]
        baseline=env["fit_data_c1c2"](inline,name="Inline Hall A DIS implementation reproduction")
        fixed=env["fit_data_c1c2_beta"](inline,beta_HT=0.,beta_CSV=0.,name="Fixed-beta implementation cross-check")
        noaxischange=env["fit_data_c1c2"](inline,confidence_levels=(.683,.95),name="Confidence-level override probe")
    if not np.allclose(baseline["params"],fixed["params"][:2],rtol=1e-6,atol=1e-8):raise AssertionError("fixed-beta/linear mismatch")
    if not np.allclose(baseline["cov"],fixed["cov"],rtol=1e-5,atol=1e-8):raise AssertionError("fixed-beta covariance mismatch")
    assert noaxischange["ellipses"][.683]["semi_major"]==noaxischange["ellipses"][.95]["semi_major"]
    A=inline[0]["A"];M=inline[0]["coeff"];stat=inline[0]["stat"];syst=inline[0]["syst"]
    V=np.diag(stat**2)+np.outer(syst,syst)
    info=M.T@np.linalg.solve(V,M);cov=np.linalg.inv(info);params=np.linalg.solve(info,M.T@np.linalg.solve(V,A))
    assert np.allclose(params,baseline["params"],atol=1e-12)
    d68=float(chi2.ppf(.683,2));factor=math.sqrt(d68)
    expected_cached=np.array([-.777397,-.030059])
    assert np.allclose(params,expected_cached,atol=5e-7,rtol=0)
    probes += [{"probe":"two_confidence_levels_same_notebook_axes","status":"confirmed","delta_683":noaxischange["ellipses"][.683]["delta_chi2"],"delta_95":noaxischange["ellipses"][.95]["delta_chi2"]},
        {"probe":"DeltaChi2_1_enclosed_probability_df2","probability":float(chi2.cdf(1,2)),"joint_683_delta":d68,"joint_683_axis_multiplier":factor},
        {"probe":"fixed_zero_betas_match_linear_GLS","status":"PASS"},
        {"probe":"fresh_inline_fit_matches_cached_printed_center","status":"PASS","tolerance":5e-7}]
    # A separately defined X-only row is a numerical probe of the beta mask,
    # not a new physical measurement or projected dataset.
    probe_rows={"name":"mask probe (not physics data)","A":np.array([1.]),"stat":np.array([1.]),
        "syst":np.array([0.]),"coeff":np.array([[1.,0.]]),"f_corr":0.,"x":np.array([.3]),"Q2":np.array([1.])}
    with contextlib.redirect_stdout(log):
        a=env["fit_data_c1c2_beta"](inline+[probe_rows],beta_HT=0.,beta_CSV=0.)
        # Changing fixed betas also affects DIS rows, so do not infer a global
        # fit invariant here. The mask's local behavior is inspected explicitly.
    mask_suppresses=(probe_rows["coeff"][:,1]!=0).tolist()==[False]
    probes.append({"probe":"X_only_row_hadronic_mask","is_numerical_probe":True,"hadronic_correction_disabled":mask_suppresses})
    # Exercise the notebook parser only on an audit-owned temporary input.
    parser_fixture=dest/"parser_probe_input.dat"
    parser_fixture.write_text("11 3 20 .5 2 .3 4 10 -70 2 100 30\n11 3 20 .5 2 .3 4 10 -70 2 100\n")
    parsed=env["read_data_file"](str(parser_fixture),"parser audit probe")
    probes.append({"probe":"malformed_numeric_row_silently_skipped","input_rows":2,"returned_rows":len(parsed["A"]),"confirmed":len(parsed["A"])==1})
    time_rows=[]
    for label,f in [("low_angle",.2),("high_angle",.8)]:
        time_rows.append({"allocation":label,"fraction":f,"notebook_formula":math.sqrt(f/.5),
            "inverse_sqrt_time_factor":math.sqrt(.5/f),"reference":"same-total exposure; 50:50 reference"})
    csvwrite(dest/"time_scaling_probes.csv",time_rows)
    total=.0057;amp=math.hypot(.004,.002)/total;variance=amp*amp
    correlation={"systematic_total":total,"intended_common_amplitude":math.hypot(.004,.002),
        "notebook_f":amp,"variance_fraction_for_listed_sources":variance,
        "notebook_common_variance":total*total*amp,"listed_common_variance":.004**2+.002**2,
        "ratio":amp/variance}
    result={"schema_version":"pvdis_notebook_numeric_audit/v1","scope":"implementation checks, not a new world fit or Hall C projection",
        "inline_data":jsonable(inline),"fresh_linear_fit":jsonable(baseline),"fixed_beta_fit":jsonable(fixed),
        "independent_GLS":{"params":params.tolist(),"cov":cov.tolist()},"correlation_fraction_probe":correlation,"probes":probes}
    write_json(dest/"numeric_diagnostics.json",jsonable(result))
    (dest/"selected_function_stdout.log").write_text(log.getvalue())
    csvwrite(dest/"inline_halla_DIS.csv",[{"point":i+1,"A_ppm":float(A[i]),"stat_ppm":float(stat[i]),"syst_ppm":float(syst[i]),
        "x":float(inline[0]["x"][i]),"Q2_GeV2":float(inline[0]["Q2"][i]),"coeff_C1_ppm":float(M[i,0]),"coeff_C2_ppm":float(M[i,1])} for i in range(2)])
    csvwrite(dest/"ellipse_convention.csv",[{"convention":"notebook_DeltaChi2_1","delta_chi2":1.,"enclosed_probability_df2":float(chi2.cdf(1,2)),
        "semi_major":baseline["ellipses"][.683]["semi_major"],"semi_minor":baseline["ellipses"][.683]["semi_minor"]},
        {"convention":"joint_683_percent","delta_chi2":d68,"enclosed_probability_df2":.683,
        "semi_major":factor*baseline["ellipses"][.683]["semi_major"],"semi_minor":factor*baseline["ellipses"][.683]["semi_minor"]}])
    # One analytic geometric diagnostic; original notebook plots are not claimed
    # reproduced, because their external datasets and hidden state are absent.
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Ellipse
    fig,ax=plt.subplots(figsize=(7,5.5),layout="constrained")
    center=baseline["params"];geometry=baseline["ellipses"][.683]
    for scale,label,color in [(1.,"Notebook: Delta chi2 = 1 (39.35% joint)","#b34b00"),(factor,"Correct 68.3% joint contour","#1761a0")]:
        ax.add_patch(Ellipse(center,2*geometry["semi_major"]*scale,2*geometry["semi_minor"]*scale,
            angle=baseline["angle_deg"],fill=False,color=color,lw=2,label=label))
    ax.plot(*center,"ko",label="Fresh inline-data fit; no recentering")
    ax.set(xlim=(-1.01,-.54),ylim=(-.36,.30),xlabel="X = 2C1u - C1d",ylabel="Y = 2C2u - C2d",
        title="Notebook confidence-convention diagnostic")
    ax.grid(alpha=.2);ax.legend(fontsize=8,loc="lower left")
    fig.savefig(dest/"ellipse_convention.png",dpi=180);plt.close(fig)
    meta={"schema_version":"pvdis_notebook_review/v1","review_id":review_id,"captured_utc":datetime.now(timezone.utc).isoformat(),
        "source_sha256":notebook_hash,"source_bytes":len(original),"source_original_path":str(source),"harness_sha256":harness_hash,
        "cell_count":39,"external_file_count":len(external),"notebook_executed_top_to_bottom":False,
        "executed_cells":"Only AST-selected reviewed function definitions from 6,14,15 and data_pvdis6 assignment from 5.",
        "drive_mounted":False,"external_projection_data_loaded":False,"cached_outputs_used_as_fresh_results":False,
        "templates":{"formula":"A=(c1*X+c2*Y)*(1+beta_HT/((1-x)^3*Q2)+beta_CSV*x^2)","beta_HT_unit":"GeV^2","beta_CSV_unit":"dimensionless","prior_or_bounds":"none in fit function","apply_mask":"second coupling coefficient != 0"},
        "missing_definitions":["data_sets_hms1","data_sets_shms1","data_sets_hms2","data_sets_shms2","result_11gev_hallc1","result_11gev_hallc1_all"],
        "packages":{n:importlib.metadata.version(n) for n in ["numpy","scipy","matplotlib","Pillow"]}}
    write_json(dest/"review.json",meta)
    (dest/"REVIEW.md").write_text("# Notebook review\n\nSource SHA-256: "+notebook_hash+"\n\nCell numbers are zero-based. Full notebook execution and Hall C fits were not performed.\n\n"+
        "\n\n".join("## "+f["id"]+": "+f["title"]+"\n\nCells: "+str(f["cells"])+"\n\n"+f["evidence"]+"\n\n"+f["impact"]+"\n\n"+f["recommendation"] for f in findings()),encoding="utf-8")
    files={p.relative_to(dest).as_posix():{"sha256":sha(p),"bytes":p.stat().st_size} for p in dest.rglob("*") if p.is_file()}
    write_json(dest/"manifest.json",{"schema_version":"pvdis_notebook_manifest/v1","files":files})
    index_file=reviews/"index.json";index=json.loads(index_file.read_text()) if index_file.exists() else {"schema_version":"pvdis_notebook_index/v1","reviews":[]}
    entry={"review_id":review_id,"path":dest.relative_to(ROOT).as_posix(),"source_sha256":notebook_hash,"manifest_sha256":sha(dest/"manifest.json")}
    index["reviews"].append(entry);write_json(index_file,index)
    print("Reviewed",review_id,"with",len(external),"external files and",len(findings()),"findings")
    print("Inline fit:",params.tolist(),"68.3% joint axis multiplier:",factor)
    return entry

if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("source");a=p.parse_args();audit(a.source)
