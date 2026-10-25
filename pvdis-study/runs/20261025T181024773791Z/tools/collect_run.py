"""Validate and snapshot an executed resonance benchmark without rerunning science."""
from pathlib import Path
import csv
from datetime import datetime, timezone
import hashlib
import importlib.util
import importlib.metadata
import json
import math
import platform
import shutil
import sys

ROOT=Path(__file__).resolve().parents[1]
RESULTS=ROOT/"results"
SCHEMA="pvdis_resonance_run/v1"
MODEL_INPUTS=["vendor/F1F2IN21_v1.0.f","src/sf_probe.f90","tools/build_probe.sh"]
ANALYSIS_INPUTS=["src/analyze.py","src/plot_results.py","data/resonance_table_I.csv"]
TABLES=["build/f1f2tables/"+n+"-onshell-cdbonn.dat" for n in ["f11","f12","f22"]]

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def blob(b):return hashlib.sha1(b"blob "+str(len(b)).encode()+b"\0"+b).hexdigest()
def read_csv(path):
    with Path(path).open(newline="",encoding="utf-8") as f:return list(csv.DictReader(f))
def write_json(path,obj):Path(path).write_text(json.dumps(obj,indent=2,allow_nan=False)+"\n",encoding="utf-8")
def load_analysis():
    spec=importlib.util.spec_from_file_location("pvdis_analysis",ROOT/"src/analyze.py")
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod);return mod

def bind_generation(mode="retained_raw"):
    """Only call after an executed refresh/full run, or verified legacy bootstrap."""
    inputs=MODEL_INPUTS+ANALYSIS_INPUTS+TABLES+["build/points.dat","build/sf_values.dat"]
    outputs=["results/"+n for n in ["resonance_benchmark.csv","resonance_W_scan.csv",
        "resonance_benchmark.png","resonance_W_scan.png","resonance_benchmark.svg","resonance_W_scan.svg",
        "benchmark_receipt.json","input_receipt.json"]]
    record={"schema_version":"pvdis_generation/v1","mode":mode,
        "captured_utc":datetime.now(timezone.utc).isoformat(),
        "inputs":{n:sha(ROOT/n) for n in inputs},"outputs":{n:sha(ROOT/n) for n in outputs}}
    write_json(RESULTS/"generation_receipt.json",record)
    return record

def validate_generation():
    p=RESULTS/"generation_receipt.json"
    if not p.exists():
        # The pre-collector run has an existing hashed review manifest. Admit it
        # only when all numerical code, input, and output hashes still agree.
        old=json.loads((RESULTS/"artifact_hashes.json").read_text())
        required=MODEL_INPUTS+ANALYSIS_INPUTS+TABLES+["build/points.dat","build/sf_values.dat",
            "results/resonance_benchmark.csv","results/resonance_W_scan.csv",
            "results/resonance_benchmark.png","results/resonance_W_scan.png",
            "results/input_receipt.json","results/benchmark_receipt.json"]
        for n in required:
            if old.get(n)!=sha(ROOT/n):raise ValueError("legacy generation identity mismatch: "+n)
        return bind_generation("legacy_executed_run_with_prior_hash_manifest")
    record=json.loads(p.read_text())
    for group in ["inputs","outputs"]:
        for n,digest in record[group].items():
            if sha(ROOT/n)!=digest:raise ValueError("stale generation: "+n)
    return record

def validate():
    generation=validate_generation()
    receipt=json.loads((RESULTS/"benchmark_receipt.json").read_text())
    inp=json.loads((RESULTS/"input_receipt.json").read_text())
    raw=(ROOT/"build/sf_values.dat").read_bytes()
    if hashlib.sha256(raw).hexdigest()!=receipt["raw_output_sha256"]:
        raise ValueError("raw output differs from analysis receipt")
    if raw!=(ROOT/"build/sf_values_repeat.dat").read_bytes():raise ValueError("raw repeat mismatch")
    if blob((ROOT/"vendor/F1F2IN21_v1.0.f").read_bytes())!=inp["model_source"]["git_blob_sha1"]:
        raise ValueError("pinned model blob mismatch")
    for name,meta in inp["tables"].items():
        b=(ROOT/"build/f1f2tables"/name).read_bytes()
        if blob(b)!=meta["git_blob_sha1"] or hashlib.sha256(b).hexdigest()!=meta["sha256"]:
            raise ValueError("pinned table mismatch: "+name)
    if sha(ROOT/"data/resonance_table_I.csv")!=inp["data_source"]["sha256"]:
        raise ValueError("published data fingerprint mismatch")
    data=read_csv(ROOT/"data/resonance_table_I.csv")
    bench=read_csv(RESULTS/"resonance_benchmark.csv")
    scan=read_csv(RESULTS/"resonance_W_scan.csv")
    if len(data)!=4 or len(bench)!=4 or len(scan)!=804:raise ValueError("unexpected point counts")
    published={i:d for i,d in enumerate(data,1)}
    saved={int(r["point_id"]):r for r in bench+scan}
    if len(saved)!=808:raise ValueError("duplicate point IDs")
    planned={int(p[0]):tuple(map(float,p[1:])) for p in
             (line.split() for line in (ROOT/"build/points.dat").read_text().splitlines())}
    if set(saved)!=set(planned):raise ValueError("planned and saved point IDs differ")
    model=load_analysis(); points=[]; maximum=0.
    for line in raw.decode().splitlines():
        f=line.split();pid=int(f[0]);e,q,w,fp,f2p,fn,f2n,fd,f2d=map(float,f[1:])
        if tuple([e,q,w])!=planned[pid]:raise ValueError("raw/planned kinematics mismatch")
        row=saved[pid]
        if not all(math.isfinite(float(v)) for k,v in row.items() if k!="setting"):
            raise ValueError("nonfinite CSV value")
        free=model.asymmetry(e,q,w,fp+fn,f2p+f2n)
        closure=model.asymmetry(e,q,w,fd,f2d)
        lo=model.y_factors(e,q,w,closure["R_gamma"],.8*closure["R_gamma"])[0]
        hi=model.y_factors(e,q,w,closure["R_gamma"],1.2*closure["R_gamma"])[0]
        c1=-1.5+10/3*model.SW2;c2=-1.5+6*model.SW2
        factor=3*model.GF*q/(10*math.sqrt(2)*math.pi*model.ALPHA)*1e6
        expected={"E_GeV":e,"Q2_GeV2":q,"W_GeV":w,"F1p":fp,"F1n":fn,
            "F1_free_sum":fp+fn,"F1d_smeared":fd,"F2_free_sum":f2p+f2n,"F2d_smeared":f2d,
            "x":free["x"],"y":free["y"],"theta_deg":free["theta_deg"],
            "R_free":free["R_gamma"],"R_smeared":closure["R_gamma"],
            "Y3_free":free["Y3"],"Y3_closure":closure["Y3"],
            "A_free_ppm":free["A_ppm"],"A_closure_ppm":closure["A_ppm"],
            "coeff_C1_closure_ppm":closure["coeff_C1_ppm"],"coeff_C2_closure_ppm":closure["coeff_C2_ppm"],
            "A_closure_SW2_0232_ppm":model.asymmetry(e,q,w,fd,f2d,sw2=.232)["A_ppm"],
            "A_Rgz_08R_ppm":factor*(lo*c1+closure["Y3"]*c2),
            "A_Rgz_12R_ppm":factor*(hi*c1+closure["Y3"]*c2)}
        if pid in published:
            d=published[pid]
            if row["setting"]!=d["setting"]:raise ValueError("setting mismatch")
            expected.update({k:float(v) for k,v in d.items() if k!="setting"})
            expected["residual_closure_ppm"]=closure["A_ppm"]-float(d["A_phys_ppm"])
            expected["residual_over_total"]=expected["residual_closure_ppm"]/float(d["total_reported_ppm"])
            if abs(math.hypot(float(d["stat_ppm"]),float(d["syst_ppm"]))-float(d["total_reported_ppm"]))>.011:
                raise ValueError("published error-column inconsistency")
        for key,v in expected.items():
            err=abs(float(row[key])-v);maximum=max(maximum,err)
            if not math.isclose(float(row[key]),v,rel_tol=1e-12,abs_tol=1e-10):
                raise ValueError(f"CSV disagrees with raw-model reconstruction: {pid} {key}")
        points.append({"point_id":pid,"E_GeV":e,"Q2_GeV2":q,"W_GeV":w,
            "free_response":free,"smeared_common_response_closure":closure,
            "effective_U":6/5*(4*fp-fn),"effective_D":6/5*(4*fn-fp),
            "all_saved_columns_reconstructed":True})
    if len(points)!=808:raise ValueError("raw point count mismatch")
    gates=[{"name":n,"status":"PASS"} for n in ["pinned_model_blob","pinned_table_blobs",
        "published_data_fingerprint","generation_input_and_output_hashes","raw_output_hash",
        "byte_identical_fresh_process_repeat","unique_and_complete_point_ids",
        "all_808_points_reconstructed_from_raw","all_csv_numeric_columns_consistent"]]
    for file,phrase,name in [("invalid.log","outside benchmark domain","Fortran_invalid_input_rejection"),
                             ("missing_tables.log","required smearing table missing","missing_table_rejection")]:
        if phrase not in (ROOT/"build"/file).read_text():raise ValueError("missing guard evidence: "+file)
        gates.append({"name":name,"status":"PASS","evidence":"build/"+file})
    return {"schema_version":"pvdis_diagnostics/v1","gates":gates,
        "maximum_csv_reconstruction_abs_error":maximum,"point_count":len(points),"points":points,
        "limitations":receipt["limitations"],"generation":generation}

def data_dictionary():
    columns=list(read_csv(RESULTS/"resonance_benchmark.csv")[0])
    special={"point_id":("integer","1-4: published means; 1000*s+j: grid, s=1..4, j=0..200"),
        "setting":("label","Published setting I-IV"),"E_GeV":("GeV","Incident beam energy"),
        "Q2_GeV2":("GeV^2","Positive squared four-momentum transfer"),"W_GeV":("GeV","Hadronic mass using nucleon convention"),
        "theta_deg":("degree","Angle implied by mean Q2,W,E; not central spectrometer angle"),
        "residual_over_total":("dimensionless","(closure model - measured asymmetry)/quoted total error; not a fit statistic")}
    out={}
    for c in columns:
        if c in special:u,d=special[c]
        elif c.endswith("_ppm"):u,d="ppm",c.replace("_"," ")
        else:u,d="dimensionless",c.replace("_"," ")
        out[c]={"unit":u,"description":d}
    return {"schema_version":"pvdis_columns/v1","columns":out,
        "notes":["Structure functions are dimensionless; deuteron values are per nucleus; free sum is p+n.",
        "R=longitudinal/transverse electromagnetic response; Y3 is dimensionless.",
        "No rate, exposure, efficiency, or detector-response columns exist in this benchmark.",
        "A_Rgz_* are illustrative input variations, not theory-error endpoints."]}

def collect(label="First resonance benchmark",mode="collect_retained_evidence"):
    diagnostics=validate()
    now=datetime.now(timezone.utc);run_id=now.strftime("%Y%m%dT%H%M%S%fZ")
    runs=ROOT/"runs";runs.mkdir(exist_ok=True);dest=runs/run_id
    dest.mkdir(exist_ok=False)
    files=[ROOT/n for n in ["README.md","DEPENDENCIES.md","ROADMAP.md"]]
    for folder in ["data","src","tools","vendor","results"]:
        files += [p for p in (ROOT/folder).rglob("*") if "__pycache__" not in p.parts and p.is_file()]
    files += [ROOT/"build"/n for n in ["points.dat","sf_values.dat","sf_values_repeat.dat",
        "compiler_version.txt","build.log","invalid.log","missing_tables.log"]]
    files += [ROOT/n for n in TABLES]
    if (ROOT/"build/commands").exists():
        files += [p for p in (ROOT/"build/commands").rglob("*") if p.is_file()]
    for src in files:
        rel=src.relative_to(ROOT);t=dest/rel;t.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(src,t)
    write_json(dest/"diagnostics.json",diagnostics)
    write_json(dest/"data_dictionary.json",data_dictionary())
    meta={"schema_version":SCHEMA,"run_id":run_id,"label":label,"capture_utc":now.isoformat(),
        "collection_mode":mode,"science_scope":"mean-kinematics model benchmark",
        "capture_did_rerun_fortran":mode=="full_fortran_run",
        "python_executable":sys.executable,"python_version":sys.version,"platform":platform.platform(),
        "python_packages":{n:importlib.metadata.version(n) for n in ["matplotlib","numpy","Pillow"]},
        "total_PAC_days":30,"exposure_used_in_benchmark":False,
        "model_point_count":808,"measured_benchmark_count":4,"scan_point_count":804,
        "source_provenance":"results/input_receipt.json","generation_provenance":"results/generation_receipt.json",
        "diagnostics":"diagnostics.json","column_dictionary":"data_dictionary.json"}
    write_json(dest/"run.json",meta)
    manifest={str(p.relative_to(dest)).replace("\\","/"):{"sha256":sha(p),"bytes":p.stat().st_size}
        for p in dest.rglob("*") if p.is_file()}
    write_json(dest/"manifest.json",{"schema_version":"pvdis_snapshot_manifest/v1","files":manifest})
    entry={"run_id":run_id,"label":label,"capture_utc":meta["capture_utc"],"mode":mode,
        "path":"runs/"+run_id,"manifest_sha256":sha(dest/"manifest.json"),"diagnostic_gate_count":len(diagnostics["gates"])}
    index_file=runs/"index.json"
    index=json.loads(index_file.read_text()) if index_file.exists() else {"schema_version":"pvdis_run_index/v1","runs":[]}
    index["runs"].append(entry);write_json(index_file,index)
    print("Collected",run_id,"with",len(manifest),"files and",len(diagnostics["gates"]),"PASS gates")
    return entry

def verify_snapshot(entry):
    root=ROOT/entry["path"]
    if sha(root/"manifest.json")!=entry["manifest_sha256"]:raise ValueError("snapshot manifest changed")
    manifest=json.loads((root/"manifest.json").read_text())["files"]
    actual={str(p.relative_to(root)).replace("\\","/") for p in root.rglob("*") if p.is_file()}
    if actual!=set(manifest)|{"manifest.json"}:raise ValueError("snapshot inventory changed")
    for name,record in manifest.items():
        p=root/name
        if sha(p)!=record["sha256"] or p.stat().st_size!=record["bytes"]:
            raise ValueError("snapshot changed: "+name)
    return root

if __name__=="__main__":collect()
