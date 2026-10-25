"""Verify a sealed resonance study and optionally reconstruct its numeric CSVs."""
import importlib.util
import json
import math
import sys
from collect_run import ROOT,sha,read_csv
from verify_snapshot import verify

def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path);mod=importlib.util.module_from_spec(spec)
    old=sys.dont_write_bytecode
    try:
        sys.dont_write_bytecode=True;spec.loader.exec_module(mod)
    finally:sys.dont_write_bytecode=old
    return mod

def verify_study(entry,recompute=False):
    root=ROOT/entry["path"]
    if sha(root/"manifest.json")!=entry["manifest_sha256"]:raise ValueError("resonance manifest changed")
    verify(root)
    meta=json.loads((root/"study.json").read_text())
    raw=(root/"raw_sf.dat").read_bytes()
    if raw!=(root/"raw_sf_repeat.dat").read_bytes() or sha(root/"raw_sf.dat")!=meta["raw_sha256"]:
        raise ValueError("resonance raw output mismatch")
    identity=json.loads((root/"input_provenance.json").read_text())
    for name,record in identity["tables"].items():
        if sha(root/"f1f2tables"/name)!=record["sha256"]:raise ValueError("table differs from provenance")
    if sha(root/"source/F1F2IN21_v1.0.f")!=identity["model_source"]["local_sha256"]:
        raise ValueError("model differs from provenance")
    if recompute:
        model=load("sealed_reference_model",root/"source/analyze.py")
        calculation=load("sealed_resonance_calculation",root/"source/run_resonance_study.py")
        nominal={int(r["point_id"]):r for r in read_csv(root/"nominal_predictions.csv")}
        scenario={}
        for r in read_csv(root/"response_scenarios.csv"):
            key=(int(r["point_id"]),float(r["axial_scale"]),float(r["RgammaZ_over_Rgamma"]))
            if key in scenario:raise ValueError("duplicate response scenario")
            scenario[key]=r
        largest=0.
        for line in raw.decode().splitlines():
            f=line.split();pid=int(f[0]);e,q,w,fp,f2p,fn,f2n,fd,f2d=map(float,f[1:])
            base=calculation.weak_response(model,e,q,w,fd,f2d)
            free=calculation.weak_response(model,e,q,w,fp+fn,f2p+f2n)
            n=nominal[pid]
            for col,value in [("A_closure_ppm",base["A_ppm"]),("A_free_ppm",free["A_ppm"]),
                    ("A_vector_ppm",base["A_vector_ppm"]),("A_axial_ppm",base["A_axial_ppm"]),
                    ("Eprime_GeV",base["Eprime_GeV"]),("theta_deg",base["theta_deg"])]:
                diff=abs(float(n[col])-value);largest=max(largest,diff)
                if diff>1e-10:raise ValueError("nominal reconstruction mismatch")
            for eta in [0.,.5,1.,1.5]:
                for rho in [.8,1.,1.2]:
                    expected=calculation.weak_response(model,e,q,w,fd,f2d,eta,rho)
                    r=scenario[(pid,eta,rho)]
                    if abs(float(r["A_ppm"])-expected["A_ppm"])>1e-10:raise ValueError("scenario reconstruction mismatch")
        # Negative-response and invalid-kinematics checks use no new physical data.
        for args in [(11.,4.,1.232,1.,1.,-1.,1.),(11.,4.,1.232,1.,1.,1.,-1.)]:
            try:calculation.weak_response(model,*args)
            except ValueError:pass
            else:raise AssertionError("negative response multiplier accepted")
        print("PASS: reconstructed",len(nominal),"nominal points and",len(scenario),"response scenarios; maximum nominal difference",largest)
    return root

if __name__=="__main__":
    index=json.loads((ROOT/"resonance_studies/index.json").read_text())
    for entry in index["studies"]:verify_study(entry,"--recompute" in sys.argv)
