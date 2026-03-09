#!/usr/bin/env python3
import csv, math, pathlib, re, subprocess
ROOT = pathlib.Path(__file__).resolve().parents[1]
R = ROOT / 'Reports'

def run(cmd):
    p = subprocess.run(cmd, shell=True, cwd=ROOT, capture_output=True, text=True)
    if p.returncode != 0:
        raise RuntimeError(f'Failed:\n{cmd}\nSTDOUT:\n{p.stdout}\nSTDERR:\n{p.stderr}')

def read_eff(path):
    s=path.read_text()
    m=re.search(r'tagging_efficiency,([0-9.eE+-]+),tagging_efficiency_err,([0-9.eE+-]+)', s)
    return float(m.group(1)), float(m.group(2))

def run_step23(prefix, seed_csv, signal='auto', bkg='baseline', mass=(0.99,1.06), fit=(1.00,1.05), bins=280):
    mc=prefix+'_mc'; data=prefix+'_data'
    cmd_mc=("root -l -b -q 'Code/step2_mc_random_pairs.C("+
            f"\"Samples/merged_mc_v2.2.root\",\"{seed_csv}\",\"Reports\",\"{mc}\",\"MC\",\"{signal}\",\"{bkg}\",{mass[0]:.3f},{mass[1]:.3f},{fit[0]:.3f},{fit[1]:.3f},{bins})'")
    cmd_da=("root -l -b -q 'Code/step2_mc_random_pairs.C("+
            f"\"Samples/merged_data_v2.2.root\",\"{seed_csv}\",\"Reports\",\"{data}\",\"Data\",\"{signal}\",\"{bkg}\",{mass[0]:.3f},{mass[1]:.3f},{fit[0]:.3f},{fit[1]:.3f},{bins})'")
    run(cmd_mc); run(cmd_da)
    emc,semc=read_eff(R/f'{mc}_fit_summary.txt')
    ed,sed=read_eff(R/f'{data}_fit_summary.txt')
    sf=ed/emc
    sferr=sf*math.sqrt((sed/ed)**2+(semc/emc)**2)
    return emc,semc,ed,sed,sf,sferr

def summarize(stem, rows, base_sf, base_sferr):
    out_csv=R/f'{stem}.csv'
    with out_csv.open('w',newline='') as f:
        w=csv.writer(f)
        w.writerow(['source','variation','eff_mc','eff_mc_err','eff_data','eff_data_err','sf','sf_err','delta_sf_vs_baseline'])
        w.writerows(rows)
    by={}
    for r in rows:
        if r[0]=='baseline': continue
        by[r[0]]=max(by.get(r[0],0.0),abs(r[-1]))
    tot=math.sqrt(sum(v*v for v in by.values()))
    out_txt=R/f'{stem}_summary.txt'
    with out_txt.open('w') as f:
        f.write(f'Step2/3 SF Systematic Uncertainty Summary ({stem})\n')
        f.write(f'Baseline SF = {base_sf:.6f} +/- {base_sferr:.6f}\n\n')
        for k in ['signal_model','background_model','fit_window_binning','matching_angle']:
            f.write(f'{k}: max |delta SF| = {by.get(k,0.0):.6f}\n')
        f.write(f'\nTotal systematic (quadrature) = {tot:.6f}\n')
        f.write(f'Relative systematic = {100*tot/base_sf:.3f}%\n')
    return tot

def main():
    # matching angle alt seed
    run("root -l -b -q 'Code/step1_make_report.C(\"Samples/merged_mc_v2.2.root\",\"Reports\",0.025,\"step1_a0025_syst\")'")

    rows=[]
    base=run_step23('kk_syst_baseline','Reports/step1_best_fit_params.csv','auto','baseline',(0.99,1.06),(1.00,1.05),280)
    base_sf,base_sferr=base[4],base[5]
    rows.append(('baseline','baseline',*base,0.0))

    for v,sig in [('sig_gauss','Gauss'),('sig_doublegauss','DoubleGauss')]:
        vals=run_step23('kk_'+v,'Reports/step1_best_fit_params.csv',sig,'baseline',(0.99,1.06),(1.00,1.05),280)
        rows.append(('signal_model',v,*vals,vals[4]-base_sf))

    for v,bkg in [('bkg_allpol2','all_pol2'),('bkg_allpol3','all_pol3')]:
        vals=run_step23('kk_'+v,'Reports/step1_best_fit_params.csv','auto',bkg,(0.99,1.06),(1.00,1.05),280)
        rows.append(('background_model',v,*vals,vals[4]-base_sf))

    for v,fit,bins in [('fwbin_narrow240',(1.005,1.045),240),('fwbin_alt320',(1.00,1.05),320),('fwbin_wide200',(0.995,1.055),200)]:
        vals=run_step23('kk_'+v,'Reports/step1_best_fit_params.csv','auto','baseline',(0.99,1.06),fit,bins)
        rows.append(('fit_window_binning',v,*vals,vals[4]-base_sf))

    vals=run_step23('kk_match_angle0025','Reports/step1_a0025_syst_best_fit_params.csv','auto','baseline',(0.99,1.06),(1.00,1.05),280)
    rows.append(('matching_angle','angle_0p025_vs_0p01',*vals,vals[4]-base_sf))

    # recommended excludes wide200 stress point
    rows_rec=[r for r in rows if r[1] != 'fwbin_wide200']
    summarize('step23_sf_systematics_recommended', rows_rec, base_sf, base_sferr)
    summarize('step23_sf_systematics_stress', rows, base_sf, base_sferr)
    print('done')

if __name__=='__main__':
    main()
