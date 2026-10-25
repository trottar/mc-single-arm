"""Append reviewed-notebook evidence to the cumulative PDF."""
import json
from xml.sax.saxutils import escape
from reportlab.platypus import PageBreak,Spacer,Preformatted
from verify_notebook_review import verify_review

def append_notebook_reviews(story,root,heading,p,table,image,styles):
    index_path=root/"notebook_reviews/index.json"
    if not index_path.exists():return
    index=json.loads(index_path.read_text())
    for entry in index["reviews"]:
        folder=verify_review(entry)
        meta=json.loads((folder/"review.json").read_text())
        diag=json.loads((folder/"numeric_diagnostics.json").read_text())
        findings=json.loads((folder/"findings.json").read_text())["findings"]
        fit=diag["fresh_linear_fit"];probe=diag["probes"][1];corr=diag["correlation_fraction_probe"]
        story += [PageBreak(),heading("Fit notebook: received evidence and scope"),
            p("Review ID: "+escape(entry["review_id"]),"SmallText"),
            p("The supplied notebook has 39 cells. Its byte-exact source, exported cell text, saved-output inventory, inline data, numerical probes, and findings are preserved in a separate sealed notebook-review folder. Cell references below are zero-based."),
            table(["Evidence","Status"],[["Uploaded notebook","Received; source SHA-256: "+meta["source_sha256"]],
                ["Full notebook run","Not performed. Drive was not mounted, and external projection files were not loaded."],
                ["Selected numerical check","Only reviewed function definitions from cells 6, 14, 15 and the inline Hall A DIS dataset assignment from cell 5 were executed."],
                ["External data","34 files referenced: 28 Hall C setting tables, four HRS tables, solid11.dat, and pvdis_11gev_CJ15nlo.csv. Not supplied in the notebook."],
                ["Cached output","Historical evidence only. All execution_count fields are null; retained outputs do not establish a clean, ordered execution."],
                ["Clean-run blockers","Undefined setting groups and HallC1 result variables are additional blockers even if external files become available."]],[117,399]),Spacer(1,14),
            p("Hadronic templates are now explicit","SectionTitle"),
            Preformatted("A_i = (c1_i*X + c2_i*Y) * (1 + beta_HT/((1-x_i)^3*Q2_i)\n                                   + beta_CSV*x_i^2)\nX = 2*C1u - C1d;  Y = 2*C2u - C2d\nbeta_HT: GeV^2; beta_CSV: dimensionless",styles["CodeBlock"]),
            p("None means a free nuisance parameter; numeric values fix it. Configurations fit neither term, HT only, CSV only, or both. Both templates multiply the entire asymmetry and are applied where the second coupling coefficient is nonzero. No nuisance prior or parameter bound is supplied in the function."),
            p("This supplies an implemented ansatz, not its scientific validation. The covariance uses the X-Y submatrix of the full inverse local information matrix, which includes local nuisance marginalization when the fit is identifiable. It is not a global nonlinear profile scan."),
            PageBreak(),heading("Notebook audit: numerical findings"),
            table(["Check","Executed result"],[["Joint confidence convention",f'Delta chi2=1 encloses {100*probe["probability"]:.3f}% for df=2. Correct 68.3% joint threshold: {probe["joint_683_delta"]:.6f}; multiply both axes by {probe["joint_683_axis_multiplier"]:.6f}.'],
                ["Confidence argument ignored","Requests for 68.3% and 95% return identical axes because delta is overwritten with 1.0."],
                ["Hall C time scaling","At a 20% low-angle allocation, notebook factor is 0.632456; inverse-square-root-time factor relative to the same-total 50:50 reference is 1.581139. All setting files currently receive the low-angle factor."],
                ["SoLID correlation split",f'Notebook f={corr["notebook_f"]:.6f}; variance fraction for the listed sources is {corr["variance_fraction_for_listed_sources"]:.6f}. Off-diagonal variance ratio: {corr["ratio"]:.6f}.'],
                ["Parser rejection probe","An audit-owned two-row input with one malformed numeric row returns one row with no rejection diagnostic."],
                ["Recentered parameter records","Display centers overwrite fitted params; a four-element fitted-parameter vector can become a two-element display vector while fitted/cov_full/chi2 retain the original fit."]],[139,377]),Spacer(1,14),
            p("Inline Hall A DIS implementation check","SectionTitle"),
            p(f'Fresh unshifted fit: X={fit["params"][0]:.9f}, Y={fit["params"][1]:.9f}; sigma_X={fit["Xerr"]:.9f}, sigma_Y={fit["Yerr"]:.9f}; correlation={fit["param_corr"]:.6f}.'),
            p("Independent GLS reproduces this center and covariance. The nonlinear function with both betas fixed at zero reproduces the same result. The fresh center agrees with the saved printed notebook center to its displayed rounding. These are checks of the inline implementation, not an independently validated world fit. These two DIS points are separate from the four resonance means earlier in this report."),
            p("Thirty-day budget","SectionTitle"),
            p("The time-factor comparison above is conditional on a same-total 50:50 reference. It is not itself a new 30-day projection. If inherited errors are for 30 days per setting, each must instead be scaled by sqrt(30 days / actual allocated time). An explicit schedule must partition the user's 30 days and account for simultaneous spectrometers."),
            PageBreak(),heading("Confidence-convention diagnostic"),image(folder/"ellipse_convention.png"),Spacer(1,12),
            p("Only the two inline Hall A DIS measurements and the notebook's within-dataset covariance were used here. Both contours share the fresh unshifted fit center. The larger contour uses the correct joint 68.3% threshold; the smaller reproduces the notebook's Delta chi2=1 convention. No Hall C or SoLID projection was rerun."),
            p("A covariance ellipse with Delta chi2=1 can be reported as that convention, with its projections equal to one-standard-deviation parameter errors. It must not be labeled a joint 68.3% confidence region in two dimensions."),
            PageBreak(),heading("Notebook findings and prioritized corrections"),
            table(["ID / cells","Finding and correction"],[[f['id']+" / "+",".join(map(str,f["cells"])),f["title"]+". "+f["recommendation"]] for f in findings],[113,403],compact=True),Spacer(1,12),
            p("Next inputs and reproduction","SectionTitle"),
            p("Obtain the 34 named files in external_inputs.csv with their source hashes and exposure definitions. Resolve confidence labels, inverse time scaling, systematic-source covariance and missing setting definitions before accepting Hall C sensitivity rankings. Verify external-prior provenance and possible shared information before claiming a combined fit."),
            Preformatted("python tools/verify_notebook_review.py --recompute\npython tools/build_study.py --mode report",styles["CodeBlock"]),
            p("The notebook-review manifest pins original source, harness, numerical outputs and plots separately from the resonance-run registry. No original notebook or previous resonance snapshot was modified.","SmallText")]
