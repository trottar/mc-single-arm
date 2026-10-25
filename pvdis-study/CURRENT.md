# Current handoff - October 25, 2026

**Primary objective: inelastic/resonance A_PV for the MOLLER workshop.** Preliminary results are needed Monday afternoon October 26; more complete work follows Tuesday October 27. The broader PAC study remains the long-term project, with 30 PAC days total. Coupling ellipses and setting optimization follow the urgent resonance work.

## Established evidence

- Executed and preserved a light-valence F1F221-based comparison to four published deuteron resonance means plus 804 W-grid points.
- Independently verified model/table identities, repeatability, algebra and CSV consistency. These are local/model checks, not physical acceptance.
- Built a cumulative PDF with source, JSON, CSV, figures and logs.
- Reviewed the fit notebook and reproduced selected inline Hall A DIS code. Identified contour, time-scaling, covariance, hidden-state and recentering issues. This audit does not complete the resonance deliverable.

## Fresh physics result executed in this continuation

Study `resonance_studies/20261025T193128315988Z` contains 896 fresh 11 GeV scan points over eight angles and W=1.10-2.20 GeV, four re-evaluated published means, 10,800 controlled weak-response scenarios, and 40 representative angle/W landmarks. A second process reproduced the raw model output byte-for-byte. Reconstruction of all 900 nominal points and 10,800 scenarios from sealed source/raw data passed with zero nominal numerical difference.

The report leads with these reference curves and separates vector, axial and longitudinal-response effects. A pure-isospin-3/2 component-limit diagnostic exposes an unresolved vector-rotation issue: 1.06 versus the light-valence isoscalar ratio 0.86 at the chosen sin2(thetaW). This is not an inclusive uncertainty band or a replacement physical model.

## Active next scientific work

**Vector response first, with the current axial prescription fixed.** Obtain signed p/n helicity and longitudinal amplitudes, a defensible background/isospin prescription and coherent interference treatment over a documented W/Q2 region. Start at the Hall A benchmark kinematics before extending to high-Q2 11 GeV points. Axial resonance modeling is a separate subsequent calculation. Acceptance consistency follows where practical. Historical PDF-fit-table availability does not block this resonance work.

## Additive component investigation - October 25

Study `component_studies/20261025T210414159583Z` freshly evaluates 900 matched points and 7200 selected empirical fit terms. Inclusive output exactly reproduces the original baseline; eight terms sum to the EM response within 8.71e-16 relative accuracy. Source and earlier sealed snapshots are unchanged. Reconstructed all new CSV predictions/component rows from retained raw data successfully.

The new figure shows the reference, an explicitly conditional Delta-fit-term-only rotation, their difference and the EM strength decomposition. Only this term changes from 0.86 to 1.06; I=1/2-labeled terms, the high-mass term and background retain the baseline assumption. The axial term is exactly fixed. This is **not a complete improved weak prediction**: the fit has no p/n helicity interference products or coherent background phases, and even its Delta p/n EM strengths do not exactly obey a pure-isospin relation. See `WEAK_RESPONSE_ASSESSMENT.md` for the derivation and required theoretical inputs.

At Hall A mean I the diagnostic moves -76.857 to -88.590 ppm (data -68.97 +/- 9.09); the other shifts are -0.689, -0.178 and -0.131 ppm. These are mean-kinematics comparisons only, with no acceptance validation, combined significance or independent duality claim.

The new EM Born cross-section diagnostic exposes suppression at 11 GeV high-x resonance points. At W=1.232 GeV, 12.5 to 30 degrees changes Q2 from 4.35 to 12.22 GeV2 and the model cross section from 2.7821 to 0.00590385 nb/(sr GeV). These are not rates. Current-revision empirical support, actual HMS/SHMS acceptance and useful-rate validation remain separate unresolved gates. Numerical model domain checks do not establish reliability.

Monday's material remains preliminary and conditional; Tuesday's first priority is an input-supported vector construction, followed by acceptance consistency and a separate axial study. A missing amplitude set must remain an explicit limitation rather than a fabricated full alternative.

## Still required before a full physics claim

Weak axial resonance and longitudinal response inputs; treatment of resonance/background isospin; acceptance/radiative consistency; verified LD2 rates before statistical projections; concise workshop-ready scientific presentation. No arbitrary variation may be described as a calibrated model uncertainty.

## Preservation and Git

Do not modify original uploaded or synced references, earlier evidence snapshots, or the model source. Branch pvdis was created from 2a27ddb7442ab2148e47fe4f3e998cf1be0a9fa8. **Latest user instruction, October 25: no further data/diagnostic publication to GitHub.** Keep growing artifacts local and produce ZIP bundles for chat only when requested. Memory/roadmap and existing code may remain in Git. The earlier one-time bundle publication is superseded by this instruction; the current cleanup removes those artifacts from the branch tree and adds ignore rules. No historical purge is claimed. No future push, farm execution, merge to master or hosted CI spending is authorized. See SHARING_POLICY.md.
