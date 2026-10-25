"""Verify sealed notebook evidence; optionally recompute selected inline fits."""
from pathlib import Path
import contextlib
import importlib.util
import io
import json
import sys
from collect_run import ROOT,sha
from verify_snapshot import verify

def verify_review(entry,recompute=False):
    root=ROOT/entry["path"]
    if sha(root/"manifest.json")!=entry["manifest_sha256"]:raise ValueError("notebook manifest changed")
    verify(root)
    meta=json.loads((root/"review.json").read_text())
    if sha(root/"source.ipynb")!=meta["source_sha256"]:raise ValueError("notebook source changed")
    if recompute:
        import numpy as np
        spec=importlib.util.spec_from_file_location("sealed_notebook_audit",root/"audit_notebook.py")
        harness=importlib.util.module_from_spec(spec)
        previous=sys.dont_write_bytecode
        try:
            sys.dont_write_bytecode=True
            spec.loader.exec_module(harness)
        finally:
            sys.dont_write_bytecode=previous
        j=json.loads((root/"source.ipynb").read_text())
        cells={i:"".join(c.get("source",[])) for i,c in enumerate(j["cells"])}
        env=harness.load_selected(cells)
        with contextlib.redirect_stdout(io.StringIO()):
            fit=env["fit_data_c1c2"](env["data_pvdis6"])
        saved=json.loads((root/"numeric_diagnostics.json").read_text())["fresh_linear_fit"]
        assert np.allclose(fit["params"],saved["params"],atol=1e-12,rtol=0)
        assert np.allclose(fit["cov"],saved["cov"],atol=1e-12,rtol=0)
        print("PASS: selected inline fit reproduced from sealed notebook definitions")
    return root

if __name__=="__main__":
    index=json.loads((ROOT/"notebook_reviews/index.json").read_text())
    for entry in index["reviews"]:verify_review(entry,recompute="--recompute" in sys.argv)
