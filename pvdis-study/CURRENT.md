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

Resolve resonance/background weak rotations and axial response inputs, then acceptance/radiative consistency and verified rate normalization. Fresh 11 GeV conditional reference curves and controlled response tests are now present; they do not complete weak-model validation or statistical reach. Historical fit-table availability is not a blocker for this resonance work.

## Still required before a full physics claim

Weak axial resonance and longitudinal response inputs; treatment of resonance/background isospin; acceptance/radiative consistency; verified LD2 rates before statistical projections; concise workshop-ready scientific presentation. No arbitrary variation may be described as a calibrated model uncertainty.

## Preservation and Git

Do not modify original uploaded or synced references, earlier evidence snapshots, or the model source. Branch pvdis was created from 2a27ddb7442ab2148e47fe4f3e998cf1be0a9fa8. On October 25 the user explicitly authorized committing and pushing this study with its reproducibility bundle to pvdis for further chat analysis. Publication is separate from the numerical build workflow; no farm execution, merge to master, or hosted CI spending is authorized.
