"""Static scientific plots from saved CSV outputs; no recomputation or fitting."""
from pathlib import Path
import csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"results"
plt.rcParams.update({"font.size":11,"axes.spines.top":False,"axes.spines.right":False,
    "savefig.dpi":180,"axes.labelsize":11,"figure.facecolor":"white"})

def load(name):
    with (OUT/name).open(newline="") as f:
        return list(csv.DictReader(f))

def col(rows,name):return np.array([float(r[name]) for r in rows])

def main():
    b=load("resonance_benchmark.csv")
    s=load("resonance_W_scan.csv")
    x=np.arange(4)
    fig,(ax,res)=plt.subplots(2,1,figsize=(10,7.4),sharex=True,
        gridspec_kw={"height_ratios":[2.2,1]},layout="constrained")
    ax.errorbar(x,col(b,"A_phys_ppm"),yerr=col(b,"total_reported_ppm"),fmt="none",
        color="#666666",elinewidth=1.2,capsize=5,label="Data: quoted total error")
    ax.errorbar(x,col(b,"A_phys_ppm"),yerr=col(b,"stat_ppm"),fmt="o",color="#161616",
        elinewidth=3,capsize=0,label="Data: statistical error")
    ax.scatter(x-.07,col(b,"A_free_ppm"),marker="D",s=55,facecolor="none",edgecolor="#1761a0",
        linewidth=1.7,label="Free p+n, light-valence rotation (this work)")
    ax.scatter(x+.07,col(b,"A_closure_ppm"),marker="x",s=65,color="#b34b00",linewidth=1.8,
        label="Smeared EM + common-response closure (this work)")
    ax.scatter(x+.14,col(b,"CJ_published_ppm"),marker="s",s=35,color="#777777",
        label="Published CJ calculation (Table I; not rerun)")
    ax.set_ylabel(r"$A_{\rm PV}$ (ppm)")
    ax.set_title("First resonance benchmark: deuterium, Table I mean kinematics",loc="left",fontweight="bold")
    ax.legend(loc="lower left",fontsize=9,frameon=False)
    ax.grid(axis="y",alpha=.18)
    for values,offset,color,marker in [(col(b,"A_free_ppm"),-.05,"#1761a0","D"),
            (col(b,"A_closure_ppm"),.05,"#b34b00","x")]:
        res.scatter(x+offset,(values-col(b,"A_phys_ppm"))/col(b,"total_reported_ppm"),
            color=color,marker=marker,s=45)
    res.axhspan(-1,1,color="#eeeeee",zorder=0)
    res.axhline(0,color="black",linewidth=.8)
    res.set_ylim(-1.5,1.5)
    res.set_ylabel("(model - data)\n/ quoted total error",fontsize=10)
    res.set_xticks(x,[f'{r["setting"]}\nQ²={float(r["Q2_GeV2"]):.3f} GeV²; W={float(r["W_GeV"]):.3f} GeV' for r in b],fontsize=9)
    fig.supxlabel("Mean-kinematics comparison; no acceptance average, weak radiative corrections, or model-error band.\n"
        "Tree-level effective inputs: sin²θW=0.235; light valence; RγZ=Rγ. Residuals are not a fit statistic.",fontsize=9)
    fig.savefig(OUT/"resonance_benchmark.png")
    fig.savefig(OUT/"resonance_benchmark.svg")
    plt.close(fig)
    fig,axes=plt.subplots(2,2,figsize=(11,7.6),sharex=True,sharey=True,layout="constrained")
    for i,(ax,r) in enumerate(zip(axes.flat,b),1):
        group=[p for p in s if i*1000<=int(p["point_id"])<i*1000+201]
        w=col(group,"W_GeV"); q=float(r["Q2_GeV2"])
        ax.plot(w,col(group,"A_free_ppm")/q,color="#1761a0",label="Free p+n")
        ax.plot(w,col(group,"A_closure_ppm")/q,color="#b34b00",linestyle="--",label="Smeared-response closure")
        ax.errorbar(float(r["W_GeV"]),float(r["A_phys_ppm"])/q,
            yerr=float(r["total_reported_ppm"])/q,fmt="o",color="#161616",capsize=4,
            label="Published mean (total error)")
        ax.set_title(f'Setting {r["setting"]}: E={float(r["E_GeV"]):.3f} GeV, Q²={q:.3f} GeV²',fontsize=11)
        ax.grid(alpha=.15)
        ax.set_xlim(1.1,2.1);ax.set_ylim(-115,-45)
    axes[0,0].legend(frameon=False,fontsize=9,loc="lower left")
    fig.suptitle("Exploratory W dependence under the light-valence ansatz",fontweight="bold")
    fig.supxlabel("W (GeV)\nCurves scan fixed E and Q²; each marker is one Table I mean, not the overlapping Fig. 1 group data.",fontsize=10)
    fig.supylabel(r"$A_{\rm PV}/Q^2$ (ppm / GeV²)")
    fig.savefig(OUT/"resonance_W_scan.png")
    fig.savefig(OUT/"resonance_W_scan.svg")
    plt.close(fig)

if __name__=="__main__":main()
