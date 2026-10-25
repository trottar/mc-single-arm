"""Meaningful negative controls for stale or internally resealed evidence."""
from pathlib import Path
import csv
import json
import shutil
import tempfile
import collect_run as c

def main():
    source=c.ROOT
    entry=json.loads((source/"runs/index.json").read_text())["runs"][0]
    root=c.verify_snapshot(entry)
    results=[]
    # Copy only the sealed snapshot's science files, never the original snapshot.
    with tempfile.TemporaryDirectory(prefix="collection_check_",dir=source/"build") as name:
        base=Path(name)
        for scenario in ["stale_csv","resealed_wrong_number","changed_source"]:
            trial=base/scenario;shutil.copytree(root,trial)
            c.ROOT=trial;c.RESULTS=trial/"results"
            if scenario=="changed_source":
                with (trial/"src/analyze.py").open("a") as f:f.write("\n# changed source\n")
            else:
                p=trial/"results/resonance_benchmark.csv"
                rows=c.read_csv(p);rows[0]["A_free_ppm"]=str(float(rows[0]["A_free_ppm"])+1)
                with p.open("w",newline="") as f:
                    w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
                if scenario=="resealed_wrong_number":
                    gp=trial/"results/generation_receipt.json";g=json.loads(gp.read_text())
                    g["outputs"]["results/resonance_benchmark.csv"]=c.sha(p);c.write_json(gp,g)
            try:c.validate()
            except ValueError as error:
                expected="reconstruction" if scenario=="resealed_wrong_number" else "stale generation"
                # The detailed reconstruction error uses "disagrees".
                if scenario=="resealed_wrong_number":assert "CSV disagrees" in str(error)
                else:assert expected in str(error)
                results.append({"test":scenario,"status":"PASS","rejection":str(error)})
            else:raise AssertionError("bad evidence accepted: "+scenario)
    c.ROOT=source;c.RESULTS=source/"results"
    c.write_json(source/"results/collection_negative_controls.json",{
        "schema_version":"pvdis_collection_tests/v1","tests":results,
        "original_snapshots_verified_after_tests":bool(c.verify_snapshot(entry))})
    print("PASS:",len(results),"negative controls; original snapshot unchanged")

if __name__=="__main__":main()
