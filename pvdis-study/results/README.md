# Result and evidence files

- `resonance_benchmark.csv`: four published means, model values, measured errors, and diagnostic columns.
- `resonance_W_scan.csv`: 804 model-grid points; no additional measured points are implied.
- `benchmark_receipt.json`: physical constants, assumptions, calculation checks, limitations, and raw-output identity.
- `input_receipt.json`: exact model/table/data provenance.
- `generation_receipt.json`: binds current numerical source, inputs, CSVs, and figures to their hashes. The collector rejects stale files.
- `collection_negative_controls.json`: tests on temporary copies showing that stale or internally resealed incorrect evidence is rejected.
- `artifact_hashes.json`: **historical manifest of the first review bundle**, retained as provenance. It is not the current build manifest.

Use `runs/index.json` and each snapshot's `manifest.json` for the cumulative history. Use `output/pdf/pdf_build.json`, `pdf_qa.json`, and `bundle_manifest.json` for the current report. The original snapshot files are not edited when this working results directory is refreshed.
