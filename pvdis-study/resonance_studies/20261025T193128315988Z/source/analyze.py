"""Point-kinematics resonance benchmark; no PDFs, rates, or parameter fit."""
from pathlib import Path
import csv
import hashlib
import json
import math
import sys

ROOT = Path(__file__).resolve().parents[1]
M = 0.938272
GF = 1.166e-5
ALPHA = 7.3e-3
SW2 = 0.235

def kinematics(e, q2, w):
    if not all(math.isfinite(v) for v in (e,q2,w)) or min(e,q2)<=0 or w<=M:
        raise ValueError("invalid kinematics")
    nu = (w*w-M*M+q2)/(2*M)
    ep = e-nu
    arg = q2/(4*e*ep) if ep>0 else math.inf
    if ep<=0 or not 0<arg<=1:
        raise ValueError("unphysical scattering point")
    x=q2/(2*M*nu)
    y=nu/e
    r2=1+4*M*M*x*x/q2
    return {"x":x,"y":y,"r2":r2,"Eprime_GeV":ep,
            "theta_deg":2*math.asin(math.sqrt(arg))*180/math.pi}

def y_factors(e,q2,w,r_gamma,r_gz):
    k=kinematics(e,q2,w)
    x,y,r2=k["x"],k["y"],k["r2"]
    def denom(r):
        return 1+(1-y)**2-y*y*(1-r2/(1+r))-x*y*M/e
    d=denom(r_gamma)
    if min(1+r_gamma,1+r_gz,d)<=0:
        raise ValueError("invalid longitudinal response")
    return (1+r_gz)/(1+r_gamma)*denom(r_gz)/d, r2/(1+r_gamma)*y*(2-y)/d

def asymmetry(e,q2,w,f1,f2,sw2=SW2):
    k=kinematics(e,q2,w)
    x,y,r2=k["x"],k["y"],k["r2"]
    r=r2*f2/(2*x*f1)-1
    if min(f1,f2)<=0 or r < -1e-10:
        raise ValueError("nonphysical electromagnetic response")
    r=max(0.,r)
    y1,y3=y_factors(e,q2,w,r,r)
    c1=-1.5+(10/3)*sw2
    c2=-1.5+6*sw2
    coeff1=3*GF*q2/(10*math.sqrt(2)*math.pi*ALPHA)*1e6*y1
    coeff2=3*GF*q2/(10*math.sqrt(2)*math.pi*ALPHA)*1e6*y3
    a=coeff1*c1+coeff2*c2
    # Independent response contraction, in the half-isospin coupling convention
    # of arXiv:1304.7741 Eq. (1). This also checks the finite-M term in Y3.
    gae=-.5; gve=-.5+2*sw2
    gvu=.5-4/3*sw2; gvd=-.5+2/3*sw2
    f1gz=f1*6/5*(2*gvu-gvd)
    f2gz=f2*6/5*(2*gvu-gvd)
    f3gz=f1*18/5  # valence + identical transverse nuclear response, if smeared
    l=1-y-M*x*y/(2*e)
    b=y*y*x*f1+l*f2
    bgz=y*y*x*f1gz+l*f2gz
    direct=GF*q2/(4*math.sqrt(2)*math.pi*ALPHA)*1e6*(2*gae*bgz/b+gve*x*y*(2-y)*f3gz/b)
    if not math.isclose(a,direct,rel_tol=1e-13,abs_tol=1e-11):
        raise AssertionError("Y-factor and direct response disagree")
    return dict(k,R_gamma=r,Y1=y1,Y3=y3,coeff_C1_ppm=coeff1,coeff_C2_ppm=coeff2,
                A_ppm=a,direct_response_ppm=direct)

def write_csv(path,rows):
    with path.open("w",newline="",encoding="utf-8") as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]))
        writer.writeheader();writer.writerows(rows)

def main():
    results=ROOT/"results"
    with (ROOT/"data/resonance_table_I.csv").open(newline="") as f:
        data=list(csv.DictReader(f))
    for d in data:
        for key in list(d):
            if key!="setting":d[key]=float(d[key])
        if abs(math.hypot(d["stat_ppm"],d["syst_ppm"])-d["total_reported_ppm"])>.011:
            raise AssertionError("transcribed uncertainty columns inconsistent")
    raw=(ROOT/"build/sf_values.dat").read_bytes()
    assert raw==(ROOT/"build/sf_values_repeat.dat").read_bytes()
    points=[]
    max_rotation_residual=0.
    for line in raw.decode().splitlines():
        fields=line.split();pid=int(fields[0]);v=list(map(float,fields[1:]))
        e,q,w,fp,f2p,fn,f2n,fd,f2d=v
        assert len(v)==9 and all(math.isfinite(a) for a in v)
        u=6/5*(4*fp-fn); d=6/5*(4*fn-fp)
        if min(u,d)<0:raise ValueError("negative effective valence density")
        # Reconstruct both electromagnetic nucleons before forming the isoscalar.
        rp=(4*u+d)/18; rn=(u+4*d)/18
        max_rotation_residual=max(max_rotation_residual,abs(rp-fp),abs(rn-fn))
        free=asymmetry(e,q,w,fp+fn,f2p+f2n)
        closure=asymmetry(e,q,w,fd,f2d)
        change_r_low=y_factors(e,q,w,closure["R_gamma"],.8*closure["R_gamma"])[0]
        change_r_high=y_factors(e,q,w,closure["R_gamma"],1.2*closure["R_gamma"])[0]
        c1=-1.5+10/3*SW2;c2=-1.5+6*SW2
        scale=3*GF*q/(10*math.sqrt(2)*math.pi*ALPHA)*1e6
        row={"point_id":pid,"E_GeV":e,"Q2_GeV2":q,"W_GeV":w,
             "x":free["x"],"y":free["y"],"theta_deg":free["theta_deg"],
             "F1p":fp,"F1n":fn,"F1_free_sum":fp+fn,"F1d_smeared":fd,
             "F2_free_sum":f2p+f2n,"F2d_smeared":f2d,
             "R_free":free["R_gamma"],"R_smeared":closure["R_gamma"],
             "Y3_free":free["Y3"],"Y3_closure":closure["Y3"],
             "A_free_ppm":free["A_ppm"],"A_closure_ppm":closure["A_ppm"],
             "A_closure_SW2_0232_ppm":asymmetry(e,q,w,fd,f2d,sw2=.232)["A_ppm"],
             "A_Rgz_08R_ppm":scale*(change_r_low*c1+closure["Y3"]*c2),
             "A_Rgz_12R_ppm":scale*(change_r_high*c1+closure["Y3"]*c2),
             "coeff_C1_closure_ppm":closure["coeff_C1_ppm"],
             "coeff_C2_closure_ppm":closure["coeff_C2_ppm"]}
        points.append(row)
    assert len(points)==808
    # Independent analytic checks: flavor cancellation, physical Bjorken limit,
    # and guards against unphysical input. No parameter was fitted to the data.
    assert max_rotation_residual<1e-12
    y1,y3=y_factors(1e8,1e6,math.sqrt(M*M+2*M*1e8*.4-1e6),0,0)
    fy=.4*(2-.4)/(1+(1-.4)**2)
    assert abs(y1-1)<1e-14 and abs(y3-fy)<1e-6
    rejected=0
    for args in [(4.867,-1,1.263),(1,.95,3.),(math.nan,.95,1.263)]:
        try:kinematics(*args)
        except ValueError:rejected+=1
    assert rejected==3
    benchmark=[]
    for i,d in enumerate(data,1):
        p=next(p for p in points if p["point_id"]==i)
        benchmark.append(dict(d,**{k:v for k,v in p.items() if k not in d},
            residual_closure_ppm=p["A_closure_ppm"]-d["A_phys_ppm"],
            residual_over_total=(p["A_closure_ppm"]-d["A_phys_ppm"])/d["total_reported_ppm"]))
    write_csv(results/"resonance_benchmark.csv",benchmark)
    write_csv(results/"resonance_W_scan.csv",[p for p in points if p["point_id"]>=1000])
    config={"GF_GeV_minus2":GF,"alpha":ALPHA,"sin2thetaW":SW2,"M_GeV":M,
        "total_PAC_days":30,"exposure_used_in_this_benchmark":False,
        "asymmetry_sign":"(sigma_R-sigma_L)/(sigma_R+sigma_L)",
        "couplings":"gAe=-1/2; gVe=-1/2+2s2; gVu=1/2-4s2/3; gVd=-1/2+2s2/3; gAu=1/2; gAd=-1/2; Ciq=2*g_e*g_q",
        "Y_reference":"arXiv:0801.4791v4 Eq.13: denominator contains -x*y*M/E",
        "model_assumptions":["charge symmetry", "light valence only; qminus=qplus", "RgammaZ=Rgamma nominally",
            "free proton+neutron rotation is reference model", "smeared alternative imposes the same transverse response on F1gamma, F1gammaZ, F3gammaZ; not derived nuclear weak smearing"],
        "limitations":["mean kinematics, no acceptance average", "tree-level effective couplings, not published radiatively corrected CJ calculation",
            "no weak axial resonance amplitudes", "no radiative migration", "no rates, statistical forecast, parameter fit, or model uncertainty band"],
        "checks":{"points":len(points),"rotation_max_abs_residual":max_rotation_residual,
            "direct_response_agreement":"passed at every point for free and closure responses",
            "fresh_process_repeat":"byte identical", "invalid_python_inputs_rejected":rejected,
            "bounds_and_FP_traps":"enabled for Fortran run"},
        "raw_output_sha256":hashlib.sha256(raw).hexdigest()}
    (results/"benchmark_receipt.json").write_text(json.dumps(config,indent=2)+"\n")
    print("setting,data_ppm,free_ppm,closure_ppm,residual_over_total")
    for r in benchmark:
        print(f'{r["setting"]},{r["A_phys_ppm"]:.2f},{r["A_free_ppm"]:.3f},{r["A_closure_ppm"]:.3f},{r["residual_over_total"]:.3f}')
    print("PASS: algebra, direct response, limits, source data, and kinematic checks")
    return benchmark,points

if __name__=="__main__":main()
