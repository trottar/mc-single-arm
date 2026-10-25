"""Build the cumulative PDF exclusively from verified, indexed snapshots."""
from pathlib import Path
import csv
from datetime import datetime,timezone
import hashlib
import importlib.metadata
import json
import math
import sys
from xml.sax.saxutils import escape
from reportlab.lib import colors
from reportlab.lib.enums import TA_LEFT
from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import getSampleStyleSheet,ParagraphStyle
from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer,Table,TableStyle,PageBreak,Image,Preformatted
from collect_run import ROOT,verify_snapshot,sha,write_json

OUT=ROOT/"output/pdf"
WIDTH=516
NAVY=colors.HexColor("#173f60")
TEAL=colors.HexColor("#147e88")
STYLES=getSampleStyleSheet()
STYLES.add(ParagraphStyle(name="ReportTitle",fontName="Helvetica-Bold",fontSize=27,leading=31,textColor=NAVY,spaceAfter=16))
STYLES.add(ParagraphStyle(name="SectionTitle",fontName="Helvetica-Bold",fontSize=17,leading=21,textColor=NAVY,spaceAfter=13))
STYLES.add(ParagraphStyle(name="Text",fontName="Helvetica",fontSize=9.3,leading=13.1,spaceAfter=8))
STYLES.add(ParagraphStyle(name="SmallText",fontName="Helvetica",fontSize=8,leading=10.5,spaceAfter=6))
STYLES.add(ParagraphStyle(name="Cell",fontName="Helvetica",fontSize=7.7,leading=10,wordWrap="CJK"))
STYLES.add(ParagraphStyle(name="CellHead",fontName="Helvetica-Bold",fontSize=7.7,leading=10,textColor=colors.white))
STYLES.add(ParagraphStyle(name="CodeBlock",fontName="Courier",fontSize=7.7,leading=10,spaceAfter=9))

def p(text,style="Text"):return Paragraph(text,STYLES[style])
def heading(text):return p(escape(text),"SectionTitle")
def table(headers,rows,widths=None,compact=False):
    vals=[[p(escape(str(c)),"CellHead") for c in headers]]
    vals += [[p(escape(str(c)),"Cell") for c in row] for row in rows]
    t=Table(vals,colWidths=widths or [WIDTH/len(headers)]*len(headers),repeatRows=1,hAlign="LEFT")
    t.setStyle(TableStyle([("BACKGROUND",(0,0),(-1,0),NAVY),("VALIGN",(0,0),(-1,-1),"TOP"),
        ("LEFTPADDING",(0,0),(-1,-1),6),("RIGHTPADDING",(0,0),(-1,-1),6),
        ("TOPPADDING",(0,0),(-1,-1),3 if compact else 5),("BOTTOMPADDING",(0,0),(-1,-1),3 if compact else 5),
        ("ROWBACKGROUNDS",(0,1),(-1,-1),[colors.HexColor("#eef4f7"),colors.white]),
        ("LINEBELOW",(0,-1),(-1,-1),.4,colors.HexColor("#c3d1da"))]))
    return t
def csvrows(path):
    with path.open(newline="") as f:return list(csv.DictReader(f))
def num(row,key,n=3):return f"{float(row[key]):.{n}f}"
def image(path):
    from PIL import Image as PILImage
    with PILImage.open(path) as im:w,h=im.size
    return Image(str(path),width=WIDTH,height=WIDTH*h/w)
def page(canvas,doc):
    canvas.saveState();w,h=letter
    canvas.setStrokeColor(colors.HexColor("#ccd8df"));canvas.line(48,h-35,w-48,h-35)
    canvas.setFillColor(NAVY);canvas.setFont("Helvetica",8)
    canvas.drawString(48,h-26,"PVDIS PAC STUDY | CUMULATIVE EVIDENCE RECORD")
    canvas.setFillColor(colors.HexColor("#526574"));canvas.setFont("Helvetica",7.5)
    canvas.drawString(48,25,"Preliminary model comparisons. Total planned exposure: 30 PAC days.")
    canvas.drawRightString(w-48,25,f"Page {doc.page}");canvas.restoreState()

def render():
    config=json.loads((ROOT/"data/report_config.json").read_text())
    index=json.loads((ROOT/"runs/index.json").read_text())
    snapshots=[(entry,verify_snapshot(entry)) for entry in index["runs"]]
    if not snapshots:raise ValueError("no validated runs")
    resonance_index_path=ROOT/"resonance_studies/index.json"
    latest_resonance_id=json.loads(resonance_index_path.read_text())["studies"][-1]["study_id"] if resonance_index_path.exists() else "not yet collected"
    story=[p(config["title"],"ReportTitle"),p(config["subtitle"],"SectionTitle"),
        p("PRIMARY GOAL: "+escape(config.get("primary_goal","Inelastic/resonance A_PV"))),
        p("Study opened October 25, 2026. This document grows by adding validated run snapshots; earlier run files are retained with their own hashes."),
        table(["Scope","Current status"],[
            ["Executed physics","Fresh 11 GeV scan: 896 points, 10,800 controlled response cases, 40 landmarks; four published resonance means. Earlier 804-point W scans are retained as historical evidence."],
            ["Construction","Light-valence, charge-symmetric rotation of electromagnetic F1F2IN21; explicitly assumed weak/nuclear response."],
            ["Validated locally","Pinned input identities, repeated raw output, algebra and numerical CSV consistency; see per-run checks."],
            ["Not established","Detector acceptance, radiative migration, full weak resonance model, validated rates, parameter fit, or optimal PAC settings."],
            ["Exposure","30 PAC days total; sequential configurations partition this total; simultaneous spectrometers share exposure. Not used in this benchmark."],
        ],[113,403]),Spacer(1,15),p("Run register","SectionTitle"),
        table(["Capture UTC","Run / mode","Gates"],[[e["capture_utc"],e["label"]+"; "+e["mode"],str(e["diagnostic_gate_count"])+" PASS"] for e,_ in snapshots],[147,306,63]),Spacer(1,12),
        p("Current physics study: resonance_studies/"+escape(latest_resonance_id)+". Its separate manifest and diagnostics appear with the new result; the register above records the historical benchmark captures.","SmallText"),
        p("Reading order: fresh 11 GeV inelastic/resonance result and response constraints; historical benchmarks; received fit-notebook audit; machine-data dictionary; reproduction and references.")]
    from resonance_report import append_resonance_studies
    append_resonance_studies(story,ROOT,heading,p,table,image,STYLES)
    for entry,root in snapshots:
        run=json.loads((root/"run.json").read_text());receipt=json.loads((root/"results/benchmark_receipt.json").read_text())
        diag=json.loads((root/"diagnostics.json").read_text());inp=json.loads((root/"results/input_receipt.json").read_text())
        b=csvrows(root/"results/resonance_benchmark.csv");s=csvrows(root/"results/resonance_W_scan.csv")
        story += [PageBreak(),heading("Run results: "+entry["label"]),p("Run ID: "+escape(run["run_id"]),"SmallText"),
            p("Collection mode: "+escape(run["collection_mode"])+". This collection "+("did" if run["capture_did_rerun_fortran"] else "did not")+" rerun Fortran. The original raw values are retained and independently checked against the saved analysis."),
            table(["Setting","E (GeV)","Q2 (GeV2)","W (GeV)","Data A (ppm)","Stat.","Syst.","Total"],
                [[r["setting"],num(r,"E_GeV"),num(r,"Q2_GeV2"),num(r,"W_GeV"),num(r,"A_phys_ppm",2),num(r,"stat_ppm",2),num(r,"syst_ppm",2),num(r,"total_reported_ppm",2)] for r in b],
                [39,61,65,61,83,69,69,69]),Spacer(1,14),
            table(["Setting","Free p+n (ppm)","Closure (ppm)","Closure - data (ppm)","Residual / total"],
                [[r["setting"],num(r,"A_free_ppm"),num(r,"A_closure_ppm"),num(r,"residual_closure_ppm"),num(r,"residual_over_total")] for r in b],[45,111,110,136,114]),Spacer(1,14),
            p("All four free-p+n differences are within 0.96 of one quoted total experimental error. This is a preliminary compatibility check at mean kinematics. It is not an acceptance-averaged comparison, p-value, or proof of model accuracy."),
            p("Measured values and their published total errors are transcribed from Wang et al., arXiv:1304.7741v2, Table I [1]. Statistical and systematic columns remain separate in the CSV. No fitted parameter or additional theory error has been used to force agreement."),
            p("Files: results/resonance_benchmark.csv contains all four means and diagnostics; results/resonance_W_scan.csv contains the 804 scan rows. The original measured input table is in data/resonance_table_I.csv.","SmallText"),
            PageBreak(),heading("Measured and calculated asymmetries"),image(root/"results/resonance_benchmark.png"),Spacer(1,12),
            p("Black points: physics asymmetries, with statistical and total experimental errors. Blue diamonds: this study's free-p+n rotation. Orange crosses: this study's smeared common-response closure. Gray squares: the paper's published CJ results, included as references and not recomputed."),
            p("The residual panel divides model minus data by the quoted total experimental error. It is a diagnostic display, not a combined fit statistic. The missing covariance, acceptance average, and theory uncertainty prevent a goodness-of-fit claim."),
            PageBreak(),heading("Kinematics and electromagnetic diagnostics"),
            p("These are calculations at each published mean Q2 and W, with the quoted beam energy. The implied angle is calculated from these means; it should not be confused with the central spectrometer setting."),
            table(["Setting","x","y","Angle (deg)","Y3 free","Y3 closure"],
                [[r["setting"],num(r,"x",5),num(r,"y",5),num(r,"theta_deg",4),num(r,"Y3_free",5),num(r,"Y3_closure",5)] for r in b]),Spacer(1,15),
            table(["Setting","F1 p+n","F1 d","F1 change (%)","R free","R d"],
                [[r["setting"],num(r,"F1_free_sum",6),num(r,"F1d_smeared",6),f'{100*(float(r["F1d_smeared"])/float(r["F1_free_sum"])-1):.3f}',num(r,"R_free",5),num(r,"R_smeared",5)] for r in b]),Spacer(1,15),
            table(["Setting","Coefficient C1 (ppm)","Coefficient C2 (ppm)","Closure - free (ppm)"],
                [[r["setting"],num(r,"coeff_C1_closure_ppm",6),num(r,"coeff_C2_closure_ppm",6),f'{float(r["A_closure_ppm"])-float(r["A_free_ppm"]):.6f}'] for r in b],[48,158,158,152]),Spacer(1,14),
            p("Structure functions are dimensionless. F1 d and F2 d are per-deuteron electromagnetic quantities; the free reference is the sum of proton and neutron. The deuteron call uses CD-Bonn on-shell smearing, explicit GETFY initialization, and no quasielastic contribution."),
            p("The small closure-minus-free asymmetry is largely imposed by the shared response assumption. It is not a measured or independently derived bound on weak nuclear effects."),
            PageBreak(),heading("Physics prescription and conventions"),
            table(["Quantity","Value / definition"],[["Signed A","(sigma_R - sigma_L)/(sigma_R + sigma_L)"],
                ["GF; alpha; M",f'{receipt["GF_GeV_minus2"]} GeV^-2; {receipt["alpha"]}; {receipt["M_GeV"]} GeV'],
                ["sin2(thetaW)",str(receipt["sin2thetaW"])+"; retained effective tree-level baseline"],
                ["C1; C2","2C1u-C1d = -0.7166667; 2C2u-C2d = -0.09"],
                ["Couplings","gAe=-1/2; gVe=-1/2+2s2; gVu=1/2-4s2/3; gVd=-1/2+2s2/3; Ciq=2*g_e*g_q"],
                ["Flavor assumption","Light valence only; charge symmetry; qminus=qplus; no independently supplied sea or axial resonance amplitudes."]],[114,402]),Spacer(1,12),
            p("The effective flavor inversion is U=6(4F1p-F1n)/5 and D=6(4F1n-F1p)/5. In the isoscalar sum, transverse resonance shapes cancel from the ratios F1gammaZ/F1gamma=(6/5)(2gVu-gVd) and F3gammaZ/F1gamma=18/5. These ratios are assumptions of the effective light-valence construction, not extracted literal resonance PDFs."),
            Preformatted("r2 = 1 + 4*M*M*x*x/Q2\nR  = r2*F2/(2*x*F1) - 1\nD  = 1 + (1-y)^2 - y*y*(1-r2/(1+R)) - x*y*M/E\nY3 = r2/(1+R) * y*(2-y)/D\nA_ppm = 3*GF*Q2/(10*sqrt(2)*pi*alpha)*1e6*(C1*Y1+C2*Y3)",STYLES["CodeBlock"]),
            p("Nominally RgammaZ=Rgamma and Y1=1. The finite-mass term follows Eq. (13) of arXiv:0801.4791v4 [2], and is checked by a direct response contraction. This study uses the half-isospin coupling convention explicitly listed above, rather than importing differently normalized couplings from another formula."),
            p("The smeared alternative imposes the same transverse response on F1gamma, F1gammaZ and F3gammaZ. This is a common-response closure, not a calculated weak axial smearing model. No unsmeared weak numerator is divided by a smeared electromagnetic denominator."),
            PageBreak(),heading("Exploratory dependence on W"),image(root/"results/resonance_W_scan.png"),Spacer(1,10),
            table(["Setting","Grid W range (GeV)","Rows","Free A range (ppm)"],
                [[r["setting"],"1.100 to 2.100","201",f'{min(float(p["A_free_ppm"]) for p in s if int(p["point_id"])//1000==i):.3f} to {max(float(p["A_free_ppm"]) for p in s if int(p["point_id"])//1000==i):.3f}'] for i,r in enumerate(b,1)],[63,163,63,227]),Spacer(1,10),
            p("Each curve fixes E and Q2 and varies W. Each data point is one published Table I mean. This does not reconstruct the overlapping group-trigger bins of the paper's Fig. 1. Smooth behavior partly follows from the imposed isoscalar flavor cancellation and must not be presented as independent evidence that weak resonance structure is absent."),
            PageBreak(),heading("Input variations and interpretive limits"),
            p("These are controlled changes of assumptions, not a model uncertainty estimate or an uncertainty budget. The nominal value is the common-response closure."),
            table(["Setting","Nominal A","s2=0.232 shift","Rgz=0.8R shift","Rgz=1.2R shift"],
                [[r["setting"],num(r,"A_closure_ppm"),f'{float(r["A_closure_SW2_0232_ppm"])-float(r["A_closure_ppm"]):+.3f}',f'{float(r["A_Rgz_08R_ppm"])-float(r["A_closure_ppm"]):+.3f}',f'{float(r["A_Rgz_12R_ppm"])-float(r["A_closure_ppm"]):+.3f}'] for r in b],[48,108,120,120,120]),Spacer(1,15),
            p("All entries are in ppm. In the RgammaZ variations, the transverse flavor prescription and electromagnetic response are held fixed; only the weak longitudinal ratio changes. The sin2(thetaW) variation changes the declared effective couplings. Neither variation is tuned to the data."),
            table(["Interpretation","Boundary"],[["Supported","The specified model construction can be executed reproducibly and compared numerically with four measured means."],
                ["Not supported","A precision weak/nuclear model, pointwise duality, acceptance or radiative validation, a rate forecast, or a fitted coupling ellipse."],
                ["Physical limitations","Sea quarks and axial resonance amplitudes are absent. Common weak/nuclear response is imposed. Published corrected CJ theory is not reproduced."],
                ["Scheduling limitation","A 30-day total budget is recorded. No beam-time-dependent observable was computed in this benchmark."]],[113,403]),Spacer(1,14),
            p("Current missing inputs","SectionTitle"),table(["Input","Status / consequence"],[[x["input"],x["status"]+": "+x["impact"]] for x in config["pending"]],[150,366]),
            PageBreak(),heading("Diagnostics and source provenance"),
            p(f'{len(diag["gates"])} collection gates passed. Maximum absolute discrepancy when reconstructing saved numeric CSV columns from raw values: {diag["maximum_csv_reconstruction_abs_error"]:.3e}. This checks pipeline consistency, not independent physical correctness.'),
            table(["Gate","Status"],[[g["name"].replace("_"," "),g["status"]] for g in diag["gates"]],[446,70]),Spacer(1,12),
            p("Model commit: "+inp["model_source"]["commit"],"SmallText"),
            p("Model Git blob: "+inp["model_source"]["git_blob_sha1"],"SmallText"),
            p("Table repository commit: "+inp["table_source"]["commit"],"SmallText"),
            table(["Table","SHA-256"],[[name,meta["sha256"]] for name,meta in inp["tables"].items()],[163,353]),Spacer(1,9),
            p("Snapshot manifest SHA-256: "+entry["manifest_sha256"],"SmallText"),
            p("The complete manifest also pins code, raw output, published data, diagnostic JSON, logs, and figures. Before each report build, every indexed snapshot is checked against its manifest. Earlier files are never updated by the collection workflow.","SmallText")]
        negative=root/"results/collection_negative_controls.json"
        if negative.exists():
            controls=json.loads(negative.read_text())
            story += [Spacer(1,8),p("Negative controls on temporary copies","SmallText"),
                table(["Deliberately invalid evidence","Outcome"],[[x["test"].replace("_"," "),x["status"]+": rejected"] for x in controls["tests"]],[396,120],compact=True),
                p("These controls verify rejection of stale CSVs, changed source, and a wrong numerical value even after its output hash is resealed. Original snapshots are checked again afterward.","SmallText")]
    latest_entry,latest=snapshots[-1]
    from notebook_report import append_notebook_reviews
    append_notebook_reviews(story,ROOT,heading,p,table,image,STYLES)
    dictionary=json.loads((latest/"data_dictionary.json").read_text())
    # Explanatory definitions may improve without rewriting an old evidence
    # snapshot. The current definitions file is separately pinned in PDF metadata.
    for name,description in json.loads((ROOT/"data/column_descriptions.json").read_text())["descriptions"].items():
        if name in dictionary["columns"]:dictionary["columns"][name]["description"]=description
    columns=list(dictionary["columns"].items())
    story += [PageBreak(),heading("Machine-readable evidence and column dictionary"),
        p("Each runs/&lt;run_id&gt;/ directory is a portable snapshot. diagnostics.json includes derived responses at all 808 points; results/*.csv retain the data tables; run.json records scope and capture mode; manifest.json records SHA-256 and byte count for each artifact. JSON uses finite numeric values and schema-version labels."),
        p("Measured inputs, identifiers, and kinematics","SectionTitle"),
        table(["Column","Unit","Definition"],[[k,v["unit"],v["description"]] for k,v in columns[:13]],[162,75,279],compact=True),Spacer(1,10),
        p("For the scan, setting is encoded in point_id (1000*s+j, s=1..4, j=0..200); 1-4 identify the measured means. A_residual fields are present only for those means. The machine dictionary is retained in each snapshot.","SmallText"),
        PageBreak(),heading("Column dictionary: model responses and variations"),
        table(["Column","Unit","Definition"],[[k,v["unit"],v["description"]] for k,v in columns[13:]],[162,75,279],compact=True),Spacer(1,10),
        p("These explanatory definitions use data/column_descriptions.json, whose hash is recorded with this PDF build. They do not modify any numeric value or historical snapshot. Coefficients are in ppm per unit dimensionless coupling combination.","SmallText"),
        PageBreak(),heading("Reproduction and debugging"),
        p("The default build collects verified retained evidence and rebuilds the cumulative report; it does not rerun the Fortran model. Use refresh to re-execute analysis and plots from retained raw output, or full to compile and rerun the model. Each successful collection adds a new snapshot."),
        Preformatted("# From the study directory (Python with analysis dependencies):\npython tools/build_study.py --mode collect\npython tools/build_study.py --mode refresh\npython tools/build_study.py --mode full\n\n# Verify snapshots and regenerate PDF without adding a run:\npython tools/build_study.py --mode report\n\n# If reportlab is in a different Python installation:\npython tools/build_study.py --mode report --report-python <python>\n\n# Verify a snapshot from its own directory:\npython tools/verify_snapshot.py <snapshot-directory>",STYLES["CodeBlock"]),
        p("Install the pinned analysis and reporting requirements into appropriate environments. A single environment can be used if compatible package versions are resolved and recorded; the captured run used separate analysis and PDF interpreters. Full mode requires WSL on Windows, GCC-13 tools and the workspace compiler, or PVDIS_FC/PVDIS_RUNTIME_LIB overrides on Linux."),
        table(["Failure","First check"],[["stale generation","A source/input/output hash changed. Refresh only for analysis changes; use full for model/driver/build changes."],
            ["model/table mismatch","Check exact Git blob identity and table commit. No replacement table is accepted silently."],
            ["CSV reconstruction mismatch","Inspect point_id, raw sf_values.dat row, data_dictionary.json, and diagnostics for the affected row."],
            ["model execution fails","Read build.log, command logs, compiler version, and rejected-input logs. Model FP traps and bounds checks are enabled."],
            ["snapshot changed","Restore exact original snapshot bytes from the review bundle; do not rehash a changed snapshot as if it were the original."]],[125,391]),Spacer(1,12),
        p("Build metadata and the report's content hash are recorded beside the PDF. The review ZIP contains the PDF, all indexed snapshots, current tools, raw data, tables, and manifests; compiler package binaries are excluded to keep it portable. No commit, push, or hosted CI is part of this workflow.","SmallText"),
        PageBreak(),heading("Sources and continuing study record")]
    for i,r in enumerate(config["references"],1):
        story += [p(f'[{i}] '+escape(r["label"])),p('<link href="'+escape(r["url"])+ '">'+escape(r["url"])+"</link>","SmallText"),p(escape(r["use"]),"SmallText"),Spacer(1,8)]
    story += [p("Future additions","SectionTitle"),
        p("Add new measurements or model runs as explicitly labeled datasets with their own provenance, assumptions, units, and validation gates. Preserve measured, modeled, projected, and fitted quantities as distinct fields. Update the run registry and regenerate this PDF from the verified snapshots."),
        p("The fit notebook has now been reviewed separately, and its implemented HT/CSV templates are recorded. The external projection tables, rate-chain validation, template physics validation, detector folding, and external-prior independence remain unresolved. New evidence must not retrospectively change the assumptions of an earlier run.")]
    OUT.mkdir(parents=True,exist_ok=True)
    path=OUT/"pvdis_pac_study.pdf"
    doc=SimpleDocTemplate(str(path),pagesize=letter,leftMargin=48,rightMargin=48,topMargin=49,bottomMargin=44,
        title=config["title"]+" - Cumulative Evidence Record",author="PVDIS PAC study")
    doc.build(story,onFirstPage=page,onLaterPages=page)
    meta={"schema_version":"pvdis_pdf_build/v1","built_utc":datetime.now(timezone.utc).isoformat(),
        "pdf":"output/pdf/pvdis_pac_study.pdf","sha256":sha(path),"indexed_runs":index["runs"],
        "run_index_sha256":sha(ROOT/"runs/index.json"),"report_config_sha256":sha(ROOT/"data/report_config.json"),
        "column_descriptions_sha256":sha(ROOT/"data/column_descriptions.json"),
        "notebook_index_sha256":sha(ROOT/"notebook_reviews/index.json") if (ROOT/"notebook_reviews/index.json").exists() else None,
        "notebook_report_template_sha256":sha(ROOT/"tools/notebook_report.py"),
        "resonance_study_index_sha256":sha(ROOT/"resonance_studies/index.json") if (ROOT/"resonance_studies/index.json").exists() else None,
        "resonance_report_template_sha256":sha(ROOT/"tools/resonance_report.py"),
        "renderer_sha256":sha(Path(__file__)),"python":sys.executable,"python_version":sys.version,
        "packages":{n:importlib.metadata.version(n) for n in ["reportlab","pypdf","Pillow"]}}
    write_json(OUT/"pdf_build.json",meta)
    print("Built",path)
    return path

if __name__=="__main__":render()
