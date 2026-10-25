"""Fresh resonance/inelastic model curves and explicitly labeled response probes."""
from pathlib import Path
import csv
from datetime import datetime,timezone
import hashlib
import importlib.util
import json
import math
import os
import shutil
import subprocess
import sys
import numpy as np
from collect_run import ROOT,sha,write_json,read_csv

def get_model():
    spec=importlib.util.spec_from_file_location("reference_asymmetry",ROOT/"src/analyze.py")
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod);return mod

def csvwrite(path,rows):
    with path.open("w",newline="",encoding="utf-8") as f:
        out=csv.DictWriter(f,fieldnames=list(rows[0]));out.writeheader();out.writerows(rows)

def run(command,cwd,logbase):
    p=subprocess.run(command,cwd=cwd,capture_output=True,text=True,encoding="utf-8",errors="replace")
    logbase.with_suffix(".stdout.log").write_text(p.stdout,encoding="utf-8")
    logbase.with_suffix(".stderr.log").write_text(p.stderr,encoding="utf-8")
    write_json(logbase.with_suffix(".json"),{"command":command,"cwd":str(cwd),"returncode":p.returncode})
    if p.returncode:raise RuntimeError("model stage failed; see "+str(logbase))

def shell(script,args,cwd,logbase):
    def linux(p):
        p=Path(p).resolve();return "/mnt/"+p.drive[0].lower()+p.as_posix()[2:]
    command=["wsl","--exec","sh",linux(script),*[linux(p) for p in args]] if os.name=="nt" else ["sh",str(script),*map(str,args)]
    run(command,cwd,logbase)

def weak_response(model,e,q,w,f1,f2,axial_scale=1.,longitudinal_scale=1.,vector_ratio=None):
    if axial_scale<0 or longitudinal_scale<0:raise ValueError("negative response multiplier")
    kin=model.kinematics(e,q,w);x,y,r2=kin["x"],kin["y"],kin["r2"]
    r=r2*f2/(2*x*f1)-1
    if r< -1e-10:raise ValueError("negative electromagnetic longitudinal response")
    r=max(0.,r);rgz=longitudinal_scale*r
    y1,y3=model.y_factors(e,q,w,r,rgz)
    sw=model.SW2;gae=-.5;gve=-.5+2*sw
    k1=6/5*(1.5-10/3*sw) if vector_ratio is None else vector_ratio
    k3=18/5*axial_scale
    pref=model.GF*q/(4*math.sqrt(2)*math.pi*model.ALPHA)*1e6
    vector=pref*2*gae*k1*y1;axial=pref*gve*k3*y3
    # Check the general response, including unequal longitudinal ratios, by
    # a direct response contraction rather than only the Y-factor expression.
    f1gz=k1*f1;f2gz=2*x*f1gz*(1+rgz)/r2;f3gz=k3*f1
    l=1-y-model.M*x*y/(2*e);b=y*y*x*f1+l*f2;bgz=y*y*x*f1gz+l*f2gz
    direct=pref*(2*gae*bgz/b+gve*x*y*(2-y)*f3gz/b)
    if not math.isclose(vector+axial,direct,rel_tol=1e-12,abs_tol=1e-10):raise ValueError("general response check failed")
    return dict(kin,R_gamma=r,R_gammaZ=rgz,Y1=y1,Y3=y3,F1gammaZ_over_F1gamma=k1,
        F3gammaZ_over_F1gamma=k3,A_vector_ppm=vector,A_axial_ppm=axial,A_ppm=vector+axial,
        direct_response_residual_ppm=vector+axial-direct)

def main():
    config=json.loads((ROOT/"data/resonance_physics_config.json").read_text())
    model=get_model();stamp=datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    case=ROOT/"resonance_studies"/stamp;case.mkdir(parents=True,exist_ok=False)
    (case/"logs").mkdir();(case/"source").mkdir();(case/"f1f2tables").mkdir()
    # Freeze every numerical input and algorithm used by this new study.
    for n in ["src/analyze.py","src/sf_probe.f90","vendor/F1F2IN21_v1.0.f","tools/build_probe.sh",
              "tools/run_resonance_probe.sh","tools/run_resonance_study.py"]:
        shutil.copyfile(ROOT/n,case/"source"/Path(n).name)
    shutil.copyfile(ROOT/"data/resonance_table_I.csv",case/"published_data.csv")
    write_json(case/"physics_config.json",config)
    identity=json.loads((ROOT/"results/input_receipt.json").read_text())
    if sha(ROOT/"vendor/F1F2IN21_v1.0.f")!=identity["model_source"]["local_sha256"]:raise ValueError("model identity changed")
    for name,meta in identity["tables"].items():
        src=ROOT/"build/f1f2tables"/name
        if sha(src)!=meta["sha256"]:raise ValueError("table identity changed")
        shutil.copyfile(src,case/"f1f2tables"/name)
    write_json(case/"input_provenance.json",identity)
    plans=[];published=read_csv(case/"published_data.csv")
    for i,r in enumerate(published,1):
        plans.append({"point_id":i,"role":"published_mean","setting":r["setting"],"E_GeV":float(r["E_GeV"]),
            "Q2_GeV2":float(r["Q2_GeV2"]),"W_GeV":float(r["W_GeV"]),"requested_theta_deg":None})
    ws=sorted(set(round(1.1+j*.01,6) for j in range(111))|set(config["landmark_W_GeV"]))
    for i,theta in enumerate(config["angles_deg"],1):
        e=config["beam_energy_GeV"];sin2=math.sin(math.radians(theta)/2)**2
        for j,w in enumerate(ws):
            ep=(model.M*model.M+2*model.M*e-w*w)/(2*model.M+4*e*sin2)
            q=4*e*ep*sin2;kin=model.kinematics(e,q,w)
            if not math.isclose(kin["theta_deg"],theta,abs_tol=1e-10):raise ValueError("angle inversion failed")
            if not (0<q<30 and w*w<20 and 0<kin["x"]<1):raise ValueError("outside chosen domain")
            plans.append({"point_id":10000+1000*i+j,"role":"11GeV_scan","setting":"reference_angle",
                "E_GeV":e,"Q2_GeV2":q,"W_GeV":w,"requested_theta_deg":theta})
    csvwrite(case/"planned_kinematics.csv",plans)
    (case/"points.dat").write_text("".join(f'{r["point_id"]} {r["E_GeV"]:.17g} {r["Q2_GeV2"]:.17g} {r["W_GeV"]:.17g}\n' for r in plans))
    shell(ROOT/"tools/build_probe.sh",[],ROOT,case/"logs/build_model")
    shutil.copyfile(ROOT/"build/compiler_version.txt",case/"compiler_version.txt")
    shutil.copyfile(ROOT/"build/build.log",case/"logs/compiler.log")
    binary_sha=sha(ROOT/"build/sf_probe")
    shell(ROOT/"tools/run_resonance_probe.sh",[case],ROOT,case/"logs/run_model")
    if (case/"raw_sf.dat").read_bytes()!=(case/"raw_sf_repeat.dat").read_bytes():raise ValueError("repeat mismatch")
    byid={r["point_id"]:r for r in plans};nominal=[];scenarios=[];perpoint=[];delta=[]
    for line in (case/"raw_sf.dat").read_text().splitlines():
        fields=line.split();pid=int(fields[0]);e,q,w,fp,f2p,fn,f2n,fd,f2d=map(float,fields[1:])
        plan=byid[pid]
        if any(not math.isclose(v,plan[k],rel_tol=1e-14) for v,k in [(e,"E_GeV"),(q,"Q2_GeV2"),(w,"W_GeV")]):raise ValueError("raw/planned point mismatch")
        free=weak_response(model,e,q,w,fp+fn,f2p+f2n)
        closure=weak_response(model,e,q,w,fd,f2d)
        assert math.isclose(free["A_ppm"],model.asymmetry(e,q,w,fp+fn,f2p+f2n)["A_ppm"],abs_tol=1e-10)
        base={"point_id":pid,"role":plan["role"],"E_GeV":e,"Q2_GeV2":q,"W_GeV":w,
            "theta_deg":free["theta_deg"],"Eprime_GeV":free["Eprime_GeV"],"x":free["x"],"y":free["y"],
            "F1_free_sum":fp+fn,"F2_free_sum":f2p+f2n,"F1d_smeared":fd,"F2d_smeared":f2d,
            "R_free":free["R_gamma"],"R_closure":closure["R_gamma"],
            "A_free_ppm":free["A_ppm"],"A_closure_ppm":closure["A_ppm"],
            "A_vector_ppm":closure["A_vector_ppm"],"A_axial_ppm":closure["A_axial_ppm"],
            "axial_fraction_abs":abs(closure["A_axial_ppm"]/closure["A_ppm"]),
            "statistical_error_status":"unavailable_unverified_rates"}
        nominal.append(base)
        pointcases=[]
        for eta in config["axial_response_scales"]:
            for rho in config["RgammaZ_over_Rgamma"]:
                response=weak_response(model,e,q,w,fd,f2d,eta,rho)
                row={"point_id":pid,"axial_scale":eta,"RgammaZ_over_Rgamma":rho,
                    "A_ppm":response["A_ppm"],"A_vector_ppm":response["A_vector_ppm"],"A_axial_ppm":response["A_axial_ppm"],
                    "shift_from_nominal_ppm":response["A_ppm"]-closure["A_ppm"]}
                scenarios.append(row);pointcases.append(dict(row,response=response))
        if math.isclose(w,1.232,abs_tol=1e-9) or pid==1:
            iso=weak_response(model,e,q,w,fd,f2d,axial_scale=0,vector_ratio=2-4*model.SW2)
            delta.append({"point_id":pid,"W_GeV":w,"theta_deg":free["theta_deg"],"Q2_GeV2":q,
                "valence_vector_ratio":closure["F1gammaZ_over_F1gamma"],"pure_I32_vector_ratio":2-4*model.SW2,
                "valence_vector_only_ppm":closure["A_vector_ppm"],"pure_I32_vector_only_ppm":iso["A_ppm"],
                "interpretation":"component_limit_only_not_inclusive_prediction"})
        perpoint.append({"point_id":pid,"planned":plan,"free_reference":free,"common_response_closure":closure,"response_scenarios":pointcases})
    if len(nominal)!=len(plans) or len({r["point_id"] for r in nominal})!=len(plans):raise ValueError("incomplete model output")
    csvwrite(case/"nominal_predictions.csv",nominal);csvwrite(case/"response_scenarios.csv",scenarios)
    csvwrite(case/"isospin_component_limits.csv",delta)
    representative=[r for r in nominal if r["role"]=="11GeV_scan" and r["W_GeV"] in config["landmark_W_GeV"]]
    csvwrite(case/"representative_kinematics.csv",representative)
    comparisons=[]
    for i,data in enumerate(published,1):
        p=next(r for r in nominal if r["point_id"]==i);rows=[r for r in scenarios if r["point_id"]==i]
        get=lambda eta,rho:next(r["A_ppm"] for r in rows if r["axial_scale"]==eta and r["RgammaZ_over_Rgamma"]==rho)
        comparisons.append({"setting":data["setting"],"point_id":i,"data_ppm":float(data["A_phys_ppm"]),
            "total_error_ppm":float(data["total_reported_ppm"]),"nominal_ppm":p["A_closure_ppm"],
            "vector_only_ppm":get(0,1),"axial_half_ppm":get(.5,1),"axial_150_percent_ppm":get(1.5,1),
            "RgammaZ_08R_ppm":get(1,.8),"RgammaZ_12R_ppm":get(1,1.2)})
    csvwrite(case/"benchmark_response_comparison.csv",comparisons)
    checks={"model_source_and_table_hashes":"PASS","fixed_angle_inversion":"PASS","physical_and_fit_domain_masks":"PASS",
        "raw_fresh_process_repeat":"PASS","complete_unique_model_points":"PASS","nominal_formula_parity":"PASS",
        "general_response_vs_direct_contraction":"PASS","scenario_count":len(scenarios),"model_point_count":len(nominal)}
    # Verify exact response endpoints and vector/axial decomposition at all points.
    for p in nominal:
        values=[r for r in scenarios if r["point_id"]==p["point_id"]]
        zero=next(r for r in values if r["axial_scale"]==0 and r["RgammaZ_over_Rgamma"]==1)
        one=next(r for r in values if r["axial_scale"]==1 and r["RgammaZ_over_Rgamma"]==1)
        assert math.isclose(zero["A_ppm"],p["A_vector_ppm"],abs_tol=1e-10)
        assert math.isclose(one["A_ppm"],p["A_closure_ppm"],abs_tol=1e-10)
        assert math.isclose(p["A_vector_ppm"]+p["A_axial_ppm"],p["A_closure_ppm"],abs_tol=1e-10)
    checks["axial_zero_nominal_and_decomposition_limits"]="PASS"
    write_json(case/"diagnostics.json",{"schema_version":"pvdis_resonance_diagnostics/v1","checks":checks,"points":perpoint})
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig,(ax,ay)=plt.subplots(2,1,figsize=(9,8),sharex=True,layout="constrained")
    for angle in [12.5,20.,30.]:
        rows=[r for r in nominal if r["role"]=="11GeV_scan" and math.isclose(r["theta_deg"],angle,abs_tol=1e-8)]
        w=[r["W_GeV"] for r in rows];a=[r["A_closure_ppm"] for r in rows]
        line=ax.plot(w,a,label=f'{angle:g} deg nominal')[0]
        for eta,style in [(0,":"),(1.5,"--")]:
            vals=[next(x["A_ppm"] for x in scenarios if x["point_id"]==r["point_id"] and x["axial_scale"]==eta and x["RgammaZ_over_Rgamma"]==1) for r in rows]
            ax.plot(w,vals,style,color=line.get_color(),alpha=.75)
        ay.plot(w,[r["Q2_GeV2"] for r in rows],color=line.get_color(),label=f'{angle:g} deg')
    ax.set(title="Fresh 11 GeV inelastic/resonance A_PV reference curves",ylabel="Physical A_PV (ppm)")
    ax.legend(fontsize=9);ax.grid(alpha=.2)
    ay.set(xlabel="W (GeV)",ylabel="Q2 (GeV²)");ay.grid(alpha=.2)
    fig.supxlabel("Solid: nominal light-valence closure; dotted: F3=0; dashed: F3 x 1.5.\nControlled input scenarios, not an uncertainty band; no acceptance or radiative folding.",fontsize=9)
    fig.savefig(case/"Apv_11GeV_reference.png",dpi=170);plt.close(fig)
    fig,ax=plt.subplots(figsize=(9,5.5),layout="constrained")
    x=np.arange(4);ax.errorbar(x,[r["data_ppm"] for r in comparisons],yerr=[r["total_error_ppm"] for r in comparisons],fmt="ko",capsize=4,label="Published physics asymmetry (total error)")
    for column,offset,marker,label in [("nominal_ppm",-.12,"D","Nominal closure"),("vector_only_ppm",-.04,"s","F3=0"),("RgammaZ_08R_ppm",.04,"v","RgammaZ=0.8R"),("RgammaZ_12R_ppm",.12,"^","RgammaZ=1.2R")]:
        ax.scatter(x+offset,[r[column] for r in comparisons],marker=marker,label=label)
    ax.set(xticks=x,xticklabels=[r["setting"] for r in comparisons],ylabel="A_PV (ppm)",title="Published resonance means: controlled weak-response tests")
    ax.legend(fontsize=8);ax.grid(axis="y",alpha=.2)
    fig.supxlabel("No parameter fitted; scenario spread is not a model uncertainty estimate.",fontsize=9)
    fig.savefig(case/"benchmark_response_scenarios.png",dpi=170);plt.close(fig)
    meta={"schema_version":"pvdis_resonance_study/v1","study_id":stamp,"created_utc":datetime.now(timezone.utc).isoformat(),
        "primary_goal":config["primary_goal"],"scan_points":len(nominal)-4,"published_means":4,"response_scenarios":len(scenarios),
        "representative_kinematics":len(representative),"angles_deg":config["angles_deg"],"beam_energy_GeV":11.,
        "constants":{"GF":model.GF,"alpha":model.ALPHA,"sin2thetaW":model.SW2,"M_GeV":model.M},
        "compiled_binary_sha256":binary_sha,"raw_sha256":sha(case/"raw_sf.dat"),
        "statistical_reach":"unavailable; rates not evaluated or normalized in this study",
        "interpretation":"fresh conditional reference curves, not a complete weak resonance model or accepted detector prediction",
        "unresolved":["weak axial resonance amplitudes","resonance/background isospin decomposition","weak longitudinal response","acceptance/radiative folding","validated rates and exposure"]}
    write_json(case/"study.json",meta)
    write_json(case/"data_dictionary.json",{"schema_version":"pvdis_resonance_columns/v1",
        "units":{"E_GeV":"GeV","Eprime_GeV":"GeV","W_GeV":"GeV","Q2_GeV2":"GeV^2","theta_deg":"degree","A_*_ppm":"ppm","F1*":"dimensionless per nucleus (p+n reference or deuteron)","R_*":"dimensionless"},
        "scenarios":{"axial_scale":"multiplier of assumed F3gammaZ/F1gamma=18/5","RgammaZ_over_Rgamma":"multiplier of weak longitudinal ratio; transverse vector ratio unchanged"},
        "reference_kinematics":"Angle/W points, not spectrometer acceptance-integrated settings or an optimized schedule."})
    files={p.relative_to(case).as_posix():{"sha256":sha(p),"bytes":p.stat().st_size} for p in case.rglob("*") if p.is_file()}
    write_json(case/"manifest.json",{"schema_version":"pvdis_resonance_manifest/v1","files":files})
    index_file=case.parent/"index.json";index=json.loads(index_file.read_text()) if index_file.exists() else {"schema_version":"pvdis_resonance_index/v1","studies":[]}
    index["studies"].append({"study_id":stamp,"path":case.relative_to(ROOT).as_posix(),"manifest_sha256":sha(case/"manifest.json")})
    write_json(index_file,index)
    print("Executed fresh resonance study:",stamp)
    print("Scan points:",len(nominal)-4,"Published means: 4; Scenarios:",len(scenarios),"Representative kinematics:",len(representative))
    print("Response checks:",checks)
    return case

if __name__=="__main__":main()
