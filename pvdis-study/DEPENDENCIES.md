# Dependency audit and remaining gates

Audit date: 2026-10-25. Public mc-single-arm baseline: `2a27ddb7442ab2148e47fe4f3e998cf1be0a9fa8` (branch `pvdis` initially points here).

| Input/capability | Current evidence | Status / implication |
|---|---|---|
| F1F2IN21 proton and neutron electromagnetic fits | Byte-exact `src/F1F2IN21_v1.0.f` at pinned repository commit, blob `02a1cea8a4850352d4c1f924de3009db2ede7aef` | Executed unchanged with a new diagnostic driver. |
| Deuteron electromagnetic response | F1F2IN21 uses `wfn=2`, on-shell smearing, QE excluded; explicit GETFY initialization is required before the deuteron call | Executed only after initializing and verifying tables. |
| CD-Bonn tables | Three existing workspace files independently matched exact Git blob identities at `trottar/simc_gfortran` commit `641b8c07182fd045fd8a5c007cdc1b71a418f3ce` | Exact copies retained in `build/f1f2tables`; SHA-256 receipts retained. These tables are not present in the audited public mc-single-arm tree. |
| Compiler | Official Ubuntu packages extracted beneath `build/compiler`, using installed GCC-13 tools | GNU Fortran 13.3.0. Bounds checks and invalid/zero/overflow traps enabled. No system compiler installation. |
| Published resonance measurements | arXiv:1304.7741v2 Table I, physics asymmetries and uncertainties | Four mean-kinematics benchmarks transcribed and uncertainty columns checked. |
| Original `READPDF3` and MMHT table | Named in the supplied Fortran; user received only that script | Not recovered. No replacement PDF grid or presumed ratio definition has been used. |
| Original `r1998`, `cross_section`, `NextUn` | External routines called by the supplied script | Their original implementations remain unavailable in the supplied material. F1F2IN21 is an independently pinned electromagnetic model, not a claim that the original script has been reproduced. |
| Original HMS/SHMS rate inputs | Script references external `hallc_code/.../worksim/` paths | Not supplied. No statistical projections or rate normalization accepted. |
| Colab fit and its inputs | User supplied Copy of plot_c1c2.ipynb; byte-exact source and a separate sealed review are retained | 39 cells reviewed; selected inline-data functions reproduced. 34 external data files and several undefined setting groups prevent full projection reproduction. |
| Gamma-Z resonance axial response | Electromagnetic data cannot determine q-minus-antiquark or axial transition amplitudes | Absent. Current light-valence rotation/common-response construction is an explicit ansatz. |
| Gamma-Z longitudinal ratio | No independent weak longitudinal response supplied | Nominal equality to electromagnetic R; 0.8R/1.2R checks are illustrative scenarios, not an uncertainty estimate. |
| Weak radiative corrections and sea/PDF prescription | Uploaded script uses tree-level effective couplings at sin²θW=0.235 | Current benchmark retains that input; not the paper's corrected electroweak calculation. |
| Detector geometry, target realization, rate efficiency | Prior local software study exists but is separate from pinned public baseline | No hardware or LD2 target acceptance validation inferred here. |
| HT/CSV templates, external priors, systematic covariance | Notebook implements multiplicative beta_HT/((1-x)^3*Q2) and beta_CSV*x^2, no nuisance priors/bounds; scalar prior values and block covariance are explicit | Implementation now known. Scientific template validation, prior independence, correct time/correlation scaling, and projection data remain unresolved. No new Hall C fit or optimum is claimed. |
| PAC time | User clarified 30 days total | Sequential paired configurations partition 30 days; simultaneous spectrometers share elapsed exposure. No asymmetry benchmark depends on this time. |

## Next numerical work

1. **Primary goal: inelastic/resonance A_PV.** Quantify axial/longitudinal and resonance-isospin assumptions against the measured resonance means; supply resonance/background weak rotations and axial amplitudes where possible.
2. Fresh 11 GeV conditional curves and representative kinematics are now executed and sealed in resonance_studies/. They remain conditional references; this is not physical acceptance of the complete weak model.
3. Establish one verified LD2 rate chain, target geometry, current, polarization, usable exposure, and detector efficiency before assigning statistical errors under the 30-day total budget.
4. Acceptance/radiative consistency and a concise workshop-ready resonance result remain active gates.
5. Broader PAC work follows: select PDF/electroweak prescriptions, repair the reviewed fit implementation, recover external projection inputs and prior provenance, then optimize settings with nuisance marginalization. Missing historical fit tables do not block steps 1-4.
