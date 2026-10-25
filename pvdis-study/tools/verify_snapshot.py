"""Standalone verification of a copied/expanded run snapshot."""
from pathlib import Path
import hashlib
import json
import sys

def verify(root):
    root=Path(root).resolve()
    manifest=json.loads((root/"manifest.json").read_text())["files"]
    actual={p.relative_to(root).as_posix() for p in root.rglob("*") if p.is_file()}
    if actual!=set(manifest)|{"manifest.json"}:raise ValueError("snapshot inventory mismatch")
    for name,record in manifest.items():
        p=root/name
        if not p.resolve().is_relative_to(root):raise ValueError("manifest path escapes snapshot")
        if hashlib.sha256(p.read_bytes()).hexdigest()!=record["sha256"] or p.stat().st_size!=record["bytes"]:
            raise ValueError("snapshot file changed: "+name)
    print("PASS:",len(manifest),"snapshot file identities")

if __name__=="__main__":
    if len(sys.argv)!=2:raise SystemExit("usage: python verify_snapshot.py <snapshot-directory>")
    verify(sys.argv[1])
