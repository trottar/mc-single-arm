# Local study artifacts and ZIP sharing

User instruction, October 25, 2026: keep the growing PVDIS data and diagnostic artifacts local. Share a ZIP bundle with chat only when requested. Do not upload data, diagnostics, figures, reports, bundles, raw model output, execution logs, copied notebooks or evidence snapshots to GitHub. Existing project memory, roadmap and source code may remain in the repository; future Git publication still requires an explicit user request.

The one-time cleanup removes previously published study artifacts from the current `pvdis` branch tree and adds ignore rules. It does not rewrite Git history or revoke previously downloaded copies. Historical commit URLs can still expose the earlier published files; this cleanup must not be described as a historical purge.

All authoritative local data, original inputs and sealed studies remain preserved. The publication checkout is a separate copy. Its removal commit does not delete the local study evidence.

Use `python pvdis-study/tools/build_study.py --mode bundle` to verify existing study snapshots and package retained outputs without rerunning science or rebuilding the PDF. The ZIP includes the current report, raw numerical data, CSV/JSON diagnostics, exact input/source copies, manifests, execution logs and CHAT_ANALYSIS_PROMPT.md. It excludes compiler binaries, bytecode and other ZIPs. The ZIP CRC and every manifested file hash are checked.

Use the included CHAT_ANALYSIS_PROMPT.md when uploading the bundle to chat. It requests explicit plots with physics interpretations and preserves the distinction between the baseline, conditional component diagnostic and a full physical weak-resonance prediction.

The PDF is the latest science-report snapshot; this policy change does not alter its numerical results. The bundle's current CURRENT.md and ROADMAP.md carry the latest publication policy.
