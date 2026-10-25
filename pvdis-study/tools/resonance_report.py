"""Lead the cumulative report with the active inelastic/resonance result."""
import json
import math
from reportlab.platypus import PageBreak,Spacer,Preformatted
from collect_run import read_csv
from verify_resonance_study import verify_study

def append_resonance_studies(story,root,heading,p,table,image,styles):
    index_path=root/"resonance_studies/index.json"
    if not index_path.exists():return
    index=json.loads(index_path.read_text())
    for entry in index["studies"]:
        folder=verify_study(entry)
        meta=json.loads((folder/"study.json").read_text());config=json.loads((folder/"physics_config.json").read_text())
        comparison=read_csv(folder/"benchmark_response_comparison.csv")
        representatives=read_csv(folder/"representative_kinematics.csv")
        isospin=read_csv(folder/"isospin_component_limits.csv")
        selected=[r for r in representatives if any(math.isclose(float(r["theta_deg"]),a,abs_tol=1e-8) for a in [12.5,20.,30.])]
        fmt=lambda r,k,n=3:f'{float(r[k]):.{n}f}'
        story += [PageBreak(),heading("PRIMARY RESULT: inelastic/resonance A_PV"),
            p("Fresh study ID: "+entry["study_id"],"SmallText"),
            p("The immediate goal is inelastic/resonance A_PV for preliminary MOLLER workshop results Monday afternoon October 26, with more complete work Tuesday October 27. Coupling-ellipse optimization is secondary. The earlier point benchmark and notebook audit do not constitute completion of this goal."),
            table(["New calculation","Executed scope"],[["11 GeV coverage",str(meta["scan_points"])+" fresh model points at eight angles from 12.5 to 30 degrees, W=1.10-2.20 GeV."],
                ["Response probes",str(meta["response_scenarios"])+" combinations of axial and longitudinal response assumptions, including the four published means."],
                ["Reference kinematics",str(meta["representative_kinematics"])+" angle/W landmarks, including Delta-region and higher-mass points; not acceptance-integrated or optimized settings."],
                ["Reproducibility","Pinned electromagnetic source/tables; fresh compilation; repeated model execution; direct-response algebra checks; per-point JSON and CSV."],
                ["Statistical reach","Unavailable until LD2 rates and exposure normalization are verified. Total planning budget remains 30 PAC days."],
                ["Scientific status","Conditional light-valence reference curves. Weak resonance amplitudes, resonance/background isospin, nuclear weak response and radiative/acceptance folding remain unresolved."]],[117,399]),Spacer(1,14),
            p("What these calculations resolve","SectionTitle"),
            p("We now separate the vector and axial asymmetry terms, vary the unknown F3 and weak longitudinal responses explicitly, and calculate the 11 GeV kinematics. No data were fitted to choose the response multipliers. Scenario differences are not a probability band or calibrated theory uncertainty."),
            p("What they do not resolve","SectionTitle"),
            p("Electromagnetic F1 and F2 do not uniquely fix the weak resonance response. In particular, the smooth isoscalar light-valence transverse ratio differs from the pure isospin-3/2 transition limit. A complete physical prediction requires resonance/background information beyond the current construction."),
            PageBreak(),heading("Fresh 11 GeV curves and axial-response dependence"),image(folder/"Apv_11GeV_reference.png"),Spacer(1,9),
            p("The plotted examples are fixed scattering angles with Q2 changing along W. The 896-row grid includes all eight angles. Solid curves use the nominal light-valence common-response closure; dotted curves remove F3; dashed curves multiply F3 by 1.5. RgammaZ=Rgamma in this figure."),
            p("These curves describe the physical asymmetry before beam-polarization dilution. They are not measured asymmetries, accepted yields, or validated hardware settings. No statistical error is assigned from unverified rates."),
            PageBreak(),heading("Representative 11 GeV inelastic/resonance kinematics"),
            p("A selected subset of the 40 landmarks is shown below. All values use the declared nucleon kinematics and model assumptions. Momentum availability, magnet configuration, mechanical constraints and target acceptance have not been validated."),
            table(["Angle (deg)","W (GeV)","Eprime (GeV)","Q2 (GeV2)","Nominal A (ppm)","Axial share (%)"],
                [[fmt(r,"theta_deg",1),fmt(r,"W_GeV"),fmt(r,"Eprime_GeV"),fmt(r,"Q2_GeV2"),fmt(r,"A_closure_ppm",2),f'{100*float(r["axial_fraction_abs"]):.2f}'] for r in selected],
                [63,62,91,89,113,98],compact=True),Spacer(1,12),
            p("At W=1.232 GeV, the nominal physical asymmetry is approximately -351.50 ppm at 12.5 degrees, -691.97 ppm at 20 degrees, and -1049.54 ppm at 30 degrees. The assumed axial term contributes about 4.34%, 7.45%, and 9.98% of the magnitude, respectively."),
            p("Files: representative_kinematics.csv, nominal_predictions.csv, response_scenarios.csv, diagnostics.json and data_dictionary.json in this sealed resonance-study directory. Larger asymmetry alone does not establish better statistical reach; rate falls and acceptance must be considered separately.","SmallText"),
            PageBreak(),heading("Published resonance data: response-assumption tests"),image(folder/"benchmark_response_scenarios.png"),Spacer(1,10),
            table(["Setting","Nominal A","F3=0 A","Axial half shift","Rgz 0.8R shift","Rgz 1.2R shift"],
                [[r["setting"],fmt(r,"nominal_ppm"),fmt(r,"vector_only_ppm"),f'{float(r["axial_half_ppm"])-float(r["nominal_ppm"]):+.3f}',f'{float(r["RgammaZ_08R_ppm"])-float(r["nominal_ppm"]):+.3f}',f'{float(r["RgammaZ_12R_ppm"])-float(r["nominal_ppm"]):+.3f}'] for r in comparison],
                [45,85,85,101,100,100],compact=True),Spacer(1,10),
            p("All table values are ppm. Removing the assumed axial contribution changes the four means by roughly 2.79-6.42 ppm. A 20% longitudinal-ratio change shifts them by roughly 2.08-5.02 ppm. Such changes can be comparable to experimental precision, and the old agreement does not isolate which weak-response assumptions are correct."),
            PageBreak(),heading("Resonance isospin: a missing physical constraint"),
            p("In the light-valence isoscalar construction, F1gammaZ/F1gamma=(6/5)(2gVu-gVd)=0.86 at sin2(thetaW)=0.235. A pure isospin-3/2 nucleon transition instead has vector-response ratio 2-4*sin2(thetaW)=1.06. The latter follows from a purely isovector transition; see Gorchtein et al., arXiv:1102.3910, Appendix A after Eq. (A9)."),
            p("The following values are component-limit diagnostics only. Applying the pure-transition ratio to the entire inclusive electromagnetic strength is not a resonance/background mixture calculation and is not adopted as the nominal inclusive prediction."),
            table(["W (GeV)","Angle (deg)","Q2 (GeV2)","Valence vector only","Pure I=3/2 vector only"],
                [[fmt(r,"W_GeV"),fmt(r,"theta_deg",2),fmt(r,"Q2_GeV2"),fmt(r,"valence_vector_only_ppm",2),fmt(r,"pure_I32_vector_only_ppm",2)] for r in isospin if r["point_id"]=="1" or any(math.isclose(float(r["theta_deg"]),a,abs_tol=1e-8) for a in [12.5,20.,30.])],
                [64,82,90,135,145],compact=True),Spacer(1,12),
            p("This is a 23.26% change in the vector response ratio before specifying any axial amplitude. At the first published mean it moves the vector-only diagnostic from -73.43 to -90.51 ppm. This establishes a model-construction limitation, not a calibrated prediction error for inclusive deuterium."),
            p("Next scientific gate","SectionTitle"),
            p("Separate resonant and nonresonant electromagnetic strengths and supply their weak rotations, including relevant proton/neutron transition-amplitude information. Establish the axial response and weak longitudinal ratio. Carlson and Rislow, arXiv:1201.3323v2, describe why resonance helicity-amplitude content matters. Until those ingredients are supplied, retain the 11 GeV curves as conditional references."),
            p("Reproduction","SectionTitle"),
            Preformatted("python tools/run_resonance_study.py\npython tools/verify_resonance_study.py --recompute\npython tools/build_study.py --mode report",styles["CodeBlock"]),
            p("The first command performs a fresh local model build/run and creates a new sealed study. The second reconstructs 900 nominal points and 10,800 scenarios from its retained source and raw values. Existing resonance and notebook snapshots are preserved.","SmallText")]
