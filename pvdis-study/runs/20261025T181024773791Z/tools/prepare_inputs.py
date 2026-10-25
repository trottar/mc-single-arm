"""Prepare exact, provenance-checked local inputs. Standard library only."""
from pathlib import Path
import csv
import hashlib
import json
import shutil

ROOT = Path(__file__).resolve().parents[1]
WORKSPACE = ROOT.parent

def git_blob(data):
    return hashlib.sha1(b"blob " + str(len(data)).encode() + b"\0" + data).hexdigest()

def main():
    identity = json.loads((ROOT / "data/table_git_identity.json").read_text())
    assert not identity["truncated"] and len(identity["tables"]) == 3
    receipt = {"table_source": identity, "tables": {}}
    target = ROOT / "build/f1f2tables"
    target.mkdir(parents=True, exist_ok=True)
    for table in identity["tables"]:
        retained = target / Path(table["path"]).name
        src = retained if retained.exists() else WORKSPACE / "mc-single-arm/src" / table["path"]
        data = src.read_bytes()
        if len(data) != table["size"] or git_blob(data) != table["sha"]:
            raise ValueError(f"table differs from pinned GitHub blob: {src}")
        dest = target / src.name
        if src.resolve() != dest.resolve():
            shutil.copyfile(src, dest)
        receipt["tables"][src.name] = {"git_blob_sha1": git_blob(data),
            "sha256": hashlib.sha256(data).hexdigest(), "bytes": len(data)}
    source = ROOT / "vendor/F1F2IN21_v1.0.f"
    source_bytes = source.read_bytes()
    # apply_patch can normalize line endings. Verify text identity against the
    # exact raw blob, allowing only the documented CRLF-to-LF transformation.
    variants = [source_bytes, source_bytes.replace(b"\r\n", b"\n"),
                source_bytes.replace(b"\r\n", b"\n").replace(b"\n", b"\r\n")]
    expected = "02a1cea8a4850352d4c1f924de3009db2ede7aef"
    matched = next((v for v in variants if git_blob(v) == expected), None)
    if matched is None:
        raise ValueError("model source differs from pinned repository blob")
    receipt["model_source"] = {"repository": "trottar/mc-single-arm",
        "commit": "2a27ddb7442ab2148e47fe4f3e998cf1be0a9fa8",
        "git_blob_sha1": expected, "local_sha256": hashlib.sha256(source_bytes).hexdigest(),
        "line_ending_only_difference": source_bytes != matched}
    with (ROOT / "data/resonance_table_I.csv").open(newline="") as f:
        rows = list(csv.DictReader(f))
    assert len(rows) == 4
    points = []
    for i, row in enumerate(rows, 1):
        points.append((i, row["E_GeV"], row["Q2_GeV2"], row["W_GeV"]))
    for i, row in enumerate(rows, 1):
        for j in range(201):
            points.append((1000*i+j, row["E_GeV"], row["Q2_GeV2"], f"{1.1+0.005*j:.6f}"))
    (ROOT / "build/points.dat").write_text("".join(" ".join(map(str,p))+"\n" for p in points))
    receipt["point_count"] = len(points)
    receipt["data_source"] = {"url": "https://arxiv.org/pdf/1304.7741v2",
        "location": "Table I, physics asymmetry results; manually transcribed and column checked",
        "sha256": hashlib.sha256((ROOT / "data/resonance_table_I.csv").read_bytes()).hexdigest()}
    (ROOT / "results").mkdir(exist_ok=True)
    (ROOT / "results/input_receipt.json").write_text(json.dumps(receipt, indent=2)+"\n")
    print(f"PASS: pinned model text, 3 exact GitHub table blobs, {len(points)} kinematic points")

if __name__ == "__main__":
    main()
