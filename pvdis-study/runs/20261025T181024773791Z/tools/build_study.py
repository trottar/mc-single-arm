"""One entry point for collection, analysis refresh, full execution, and reporting."""
from pathlib import Path
import argparse
from datetime import datetime,timezone
import importlib.util
import json
import os
import shutil
import subprocess
import sys
import zipfile
from collect_run import ROOT,MODEL_INPUTS,TABLES,sha,write_json,bind_generation,collect,verify_snapshot

def report_python(explicit):
    if explicit:return explicit
    if os.environ.get("PVDIS_REPORT_PYTHON"):return os.environ["PVDIS_REPORT_PYTHON"]
    if importlib.util.find_spec("reportlab"):return sys.executable
    bundled=Path.home()/".cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe"
    if bundled.is_file():return str(bundled)
    raise RuntimeError("reportlab unavailable; install requirements-report.txt or pass --report-python")

def run(command,label):
    folder=ROOT/"build/commands";folder.mkdir(parents=True,exist_ok=True)
    stamp=datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    proc=subprocess.run(command,cwd=ROOT,capture_output=True,text=True,encoding="utf-8",errors="replace")
    (folder/f"{stamp}_{label}.stdout.log").write_text(proc.stdout,encoding="utf-8")
    (folder/f"{stamp}_{label}.stderr.log").write_text(proc.stderr,encoding="utf-8")
    write_json(folder/f"{stamp}_{label}.json",{"schema_version":"pvdis_command/v1","command":command,
        "cwd":str(ROOT),"returncode":proc.returncode,"capture_utc":stamp})
    if proc.returncode:raise RuntimeError(f"{label} failed ({proc.returncode}); inspect {folder}")
    print(label+": PASS")

def shell(script):
    path=ROOT/"tools"/script
    if os.name=="nt":
        if not shutil.which("wsl"):raise RuntimeError("full model execution requires WSL")
        drive=path.drive[0].lower();linux="/mnt/"+drive+path.as_posix()[2:]
        run(["wsl","--exec","sh",linux],script)
    else:run(["sh",str(path)],script)

def refresh(full):
    if full:
        run([sys.executable,str(ROOT/"tools/prepare_inputs.py")],"prepare_inputs")
        shell("build_probe.sh");shell("run_probe.sh")
    else:
        previous=json.loads((ROOT/"results/generation_receipt.json").read_text())
        for n in MODEL_INPUTS+TABLES+["build/points.dat","build/sf_values.dat"]:
            if previous["inputs"].get(n)!=sha(ROOT/n):
                raise ValueError("model-generation input changed; use --mode full: "+n)
    run([sys.executable,str(ROOT/"src/analyze.py")],"analysis")
    run([sys.executable,str(ROOT/"src/plot_results.py")],"plots")
    bind_generation("full_fortran_run" if full else "analysis_refresh_from_retained_raw")

def package():
    out=ROOT/"output/pdf";path=out/"pvdis_pac_study_bundle.zip"
    files=[ROOT/n for n in ["README.md","ROADMAP.md","DEPENDENCIES.md","requirements-analysis.txt","requirements-report.txt"]]
    for d in ["data","src","tools","vendor","results","runs","output/pdf"]:
        files += [p for p in (ROOT/d).rglob("*") if "__pycache__" not in p.parts and p.is_file() and p.suffix!=".zip" and p.name!="bundle_manifest.json"]
    files += [ROOT/"build"/n for n in ["points.dat","sf_values.dat","sf_values_repeat.dat","compiler_version.txt",
        "build.log","invalid.log","missing_tables.log"]]+[ROOT/n for n in TABLES]
    files += [p for p in (ROOT/"build/commands").rglob("*") if p.is_file()]
    manifest={p.relative_to(ROOT).as_posix():{"sha256":sha(p),"bytes":p.stat().st_size} for p in files}
    write_json(out/"bundle_manifest.json",{"schema_version":"pvdis_bundle_manifest/v1","files":manifest})
    files.append(out/"bundle_manifest.json")
    with zipfile.ZipFile(path,"w",zipfile.ZIP_DEFLATED) as z:
        for p in files:z.write(p,p.relative_to(ROOT).as_posix())
    with zipfile.ZipFile(path) as z:
        if z.testzip() is not None:raise ValueError("archive CRC failure")
        import hashlib
        for n,r in manifest.items():
            if hashlib.sha256(z.read(n)).hexdigest()!=r["sha256"]:raise ValueError("archived hash mismatch: "+n)
    print("Review bundle:",path)

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode",choices=["collect","refresh","full","report"],default="collect")
    parser.add_argument("--label",default="Resonance benchmark evidence collection")
    parser.add_argument("--report-python")
    args=parser.parse_args()
    if args.mode in ["refresh","full"]:refresh(args.mode=="full")
    if args.mode!="report":
        collect(args.label,{"collect":"collect_retained_evidence","refresh":"analysis_refresh_from_retained_raw","full":"full_fortran_run"}[args.mode])
    index=json.loads((ROOT/"runs/index.json").read_text())
    for entry in index["runs"]:verify_snapshot(entry)
    run([report_python(args.report_python),str(ROOT/"tools/render_report.py")],"render_report")
    run([report_python(args.report_python),str(ROOT/"tools/verify_report.py")],"verify_report")
    package()

if __name__=="__main__":main()
