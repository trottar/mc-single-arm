"""PDF text/bounds checks plus full-page renders for visual inspection."""
from pathlib import Path
import json
import math
import pdfplumber
from pypdf import PdfReader
import pypdfium2 as pdfium
from PIL import Image,ImageDraw
from collect_run import ROOT,sha,write_json

def main():
    folder=ROOT/"output/pdf";path=folder/"pvdis_pac_study.pdf"
    reader=PdfReader(path);texts=[p.extract_text() or "" for p in reader.pages]
    index=json.loads((ROOT/"runs/index.json").read_text())
    if len(texts)<8 or any(len(t)<30 for t in texts):raise ValueError("missing or blank report pages")
    alltext="\n".join(texts)
    for entry in index["runs"]:
        if entry["run_id"] not in alltext:raise ValueError("run missing from PDF")
    for marker in ["30 PAC days","diagnostics.json","Reproduction and debugging","Column"]:
        if marker not in alltext:raise ValueError("required report content missing: "+marker)
    outside=[]
    with pdfplumber.open(path) as pdf:
        for i,p in enumerate(pdf.pages,1):
            for word in p.extract_words():
                if word["x0"]<15 or word["x1"]>p.width-15 or word["top"]<10 or word["bottom"]>p.height-10:
                    outside.append({"page":i,"word":word["text"]})
    if outside:raise ValueError("PDF text outside page bounds: "+repr(outside[:5]))
    renders=folder/"rendered";renders.mkdir(exist_ok=True)
    pdf=pdfium.PdfDocument(path);thumbs=[]
    for i,p in enumerate(pdf,1):
        im=p.render(scale=1.5).to_pil().convert("RGB")
        im.save(renders/f"page_{i:03d}.png")
        im.thumbnail((306,396));tile=Image.new("RGB",(326,426),"#dddddd")
        tile.paste(im,((326-im.width)//2,8));ImageDraw.Draw(tile).text((12,407),f"Page {i}",fill="black")
        thumbs.append(tile)
    cols=3;rows=math.ceil(len(thumbs)/cols)
    sheet=Image.new("RGB",(cols*326,rows*426),"#c4c4c4")
    for i,t in enumerate(thumbs):sheet.paste(t,((i%cols)*326,(i//cols)*426))
    sheet.save(renders/"contact_sheet.png")
    write_json(folder/"pdf_qa.json",{"schema_version":"pvdis_pdf_qa/v1","pdf_sha256":sha(path),
        "pages":len(texts),"indexed_runs_present":len(index["runs"]),"nonblank_text_pages":True,
        "text_outside_page_bounds":outside,"rendered_pages":len(thumbs),
        "visual_review":"Page images generated; automated checks do not establish visual correctness."})
    print("PASS:",len(texts),"nonblank PDF pages; run presence; text bounds; rendered all pages")

if __name__=="__main__":main()
