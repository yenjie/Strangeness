#!/usr/bin/env python3
import csv, math, pathlib, re, subprocess

ROOT = pathlib.Path(__file__).resolve().parents[1]
REPORTS = ROOT / 'Reports'

CHI2_MAX = 1.6
PULL_RMS_MAX = 1.8
PULL_MAXABS_MAX = 5.5


def run_root(cmd: str):
    p = subprocess.run(cmd, shell=True, cwd=ROOT, capture_output=True, text=True)
    if p.returncode != 0:
        raise RuntimeError(f'Command failed:\n{cmd}\nSTDOUT:\n{p.stdout}\nSTDERR:\n{p.stderr}')


def parse_summary(path: pathlib.Path):
    s = path.read_text()
    m_eff = re.search(r'tagging_efficiency,([0-9.eE+-]+),tagging_efficiency_err,([0-9.eE+-]+)', s)
    if not m_eff:
        raise RuntimeError(f'Cannot parse efficiency from {path}')
    eff, eff_err = float(m_eff.group(1)), float(m_eff.group(2))

    cats = {}
    for c in ['0tag', '1tag', '2tag']:
        m = re.search(rf'{c},model,([^,]+),chi2ndf,([0-9.eE+-]+),mean,([0-9.eE+-]+),meanErr,([0-9.eE+-]+),yield,([0-9.eE+-]+),yieldErr,([0-9.eE+-]+),pullRms,([0-9.eE+-]+),pullMeanAbs,([0-9.eE+-]+),pullMaxAbs,([0-9.eE+-]+)', s)
        if not m:
            raise RuntimeError(f'Cannot parse {c} metrics from {path}')
        cats[c] = {
            'model': m.group(1),
            'chi2': float(m.group(2)),
            'pull_rms': float(m.group(7)),
            'pull_max': float(m.group(9)),
        }
    return eff, eff_err, cats


def pass_quality(cats):
    bad = []
    for c in ['1tag', '2tag']:
        if cats[c]['chi2'] > CHI2_MAX:
            bad.append(f'{c}:chi2={cats[c]["chi2"]:.3f}')
        if cats[c]['pull_rms'] > PULL_RMS_MAX:
            bad.append(f'{c}:pullRMS={cats[c]["pull_rms"]:.3f}')
        if cats[c]['pull_max'] > PULL_MAXABS_MAX:
            bad.append(f'{c}:maxPull={cats[c]["pull_max"]:.3f}')
    return len(bad) == 0, '; '.join(bad)


def run_step56(prefix: str, seed_csv: str, signal_model: str, bkg_mode: str,
               mass_min: float, mass_max: float, fit_min: float, fit_max: float, n_bins: int):
    mc_prefix = f'{prefix}_mc'
    da_prefix = f'{prefix}_data'

    cmd_mc = (
        "root -l -b -q 'Code/step5_mc_random_pairs_kshort.C("
        f"\"Samples/merged_mc_v2.2.root\",\"{seed_csv}\",\"Reports\",\"{mc_prefix}\",\"MC\","
        f"\"{signal_model}\",\"{bkg_mode}\",{mass_min:.3f},{mass_max:.3f},{fit_min:.3f},{fit_max:.3f},{n_bins})'"
    )
    cmd_da = (
        "root -l -b -q 'Code/step5_mc_random_pairs_kshort.C("
        f"\"Samples/merged_data_v2.2.root\",\"{seed_csv}\",\"Reports\",\"{da_prefix}\",\"Data\","
        f"\"{signal_model}\",\"{bkg_mode}\",{mass_min:.3f},{mass_max:.3f},{fit_min:.3f},{fit_max:.3f},{n_bins})'"
    )
    run_root(cmd_mc)
    run_root(cmd_da)

    emc, semc, cmc = parse_summary(REPORTS / f'{mc_prefix}_fit_summary.txt')
    eda, seda, cda = parse_summary(REPORTS / f'{da_prefix}_fit_summary.txt')

    sf = eda / emc
    sf_err = sf * math.sqrt((seda / eda) ** 2 + (semc / emc) ** 2)

    q_mc_ok, q_mc_msg = pass_quality(cmc)
    q_da_ok, q_da_msg = pass_quality(cda)
    q_ok = q_mc_ok and q_da_ok
    q_msg = '; '.join(x for x in [q_mc_msg, q_da_msg] if x)

    return emc, semc, eda, seda, sf, sf_err, q_ok, q_msg


def main():
    baseline_seed = 'Reports/step4_a0025_best_fit_params.csv'

    variations = [
        ('baseline', 'baseline', 'auto', 'baseline', (0.40, 0.70, 0.40, 0.70, 280)),

        ('signal_model', 'sig_doublegauss', 'DoubleGauss', 'baseline', (0.40, 0.70, 0.40, 0.70, 280)),
        ('signal_model', 'sig_voigt', 'Voigt', 'baseline', (0.40, 0.70, 0.40, 0.70, 280)),

        ('background_model', 'bkg_allpol3', 'auto', 'all_pol3', (0.40, 0.70, 0.40, 0.70, 280)),
        ('background_model', 'bkg_allpol3plain', 'auto', 'all_pol3_plain', (0.40, 0.70, 0.40, 0.70, 280)),
        ('background_model', 'bkg_allpol4', 'auto', 'all_pol4', (0.40, 0.70, 0.40, 0.70, 280)),

        ('fit_window_binning', 'fwbin_narrow240', 'auto', 'baseline', (0.40, 0.70, 0.42, 0.68, 240)),
        ('fit_window_binning', 'fwbin_alt320', 'auto', 'baseline', (0.40, 0.70, 0.45, 0.58, 320)),
    ]

    # matching angle seed (0.01)
    run_root("root -l -b -q 'Code/step4_kshort_make_report.C(\"Samples/merged_mc_v2.2.root\",\"Reports\",0.01,\"step4_a0010_qc\",0.40,0.60)'"
    )
    variations.append(('matching_angle', 'match_angle0010', 'auto', 'baseline', (0.40, 0.70, 0.40, 0.70, 280)))

    rows = []
    base_sf = None
    base_sferr = None

    for source, name, sig, bkg, pars in variations:
        seed = baseline_seed if source != 'matching_angle' else 'Reports/step4_a0010_qc_best_fit_params.csv'
        emc, semc, eda, seda, sf, sferr, q_ok, q_msg = run_step56(name, seed, sig, bkg, *pars)
        if source == 'baseline':
            base_sf, base_sferr = sf, sferr
            delta = 0.0
        else:
            delta = sf - base_sf
        rows.append([source, name, sig, bkg, *pars, emc, semc, eda, seda, sf, sferr, delta, int(q_ok), q_msg])

    # keep only quality-pass variations for systematics, plus baseline
    selected = [r for r in rows if r[0] == 'baseline' or r[-2] == 1]

    # write full QC table
    out_all = REPORTS / 'step56_sf_systematics_qc_all.csv'
    with out_all.open('w', newline='') as f:
        w = csv.writer(f)
        w.writerow([
            'source','variation','signal_model','bkg_mode','mass_min','mass_max','fit_min','fit_max','n_bins',
            'eff_mc','eff_mc_err','eff_data','eff_data_err','sf','sf_err','delta_sf_vs_baseline','quality_pass','quality_notes'
        ])
        w.writerows(rows)

    # write selected table and summary
    out_sel = REPORTS / 'step56_sf_systematics_qc_selected.csv'
    with out_sel.open('w', newline='') as f:
        w = csv.writer(f)
        w.writerow([
            'source','variation','signal_model','bkg_mode','mass_min','mass_max','fit_min','fit_max','n_bins',
            'eff_mc','eff_mc_err','eff_data','eff_data_err','sf','sf_err','delta_sf_vs_baseline','quality_pass','quality_notes'
        ])
        w.writerows(selected)

    by = {}
    for r in selected:
        src = r[0]
        if src == 'baseline':
            continue
        by[src] = max(by.get(src, 0.0), abs(float(r[15])))

    total = math.sqrt(sum(v*v for v in by.values()))

    out_txt = REPORTS / 'step56_sf_systematics_qc_selected_summary.txt'
    with out_txt.open('w') as f:
        f.write('Step5/6 SF Systematics with Fit-Quality Selection\n')
        f.write(f'Quality thresholds: chi2/ndf<{CHI2_MAX}, pullRMS<{PULL_RMS_MAX}, max|pull|<{PULL_MAXABS_MAX} for 1tag/2tag in both MC and Data\n\n')
        f.write(f'Baseline SF = {base_sf:.6f} +/- {base_sferr:.6f}\n\n')
        for src in ['signal_model','background_model','fit_window_binning','matching_angle']:
            f.write(f'{src}: max |delta SF| = {by.get(src,0.0):.6f}\n')
        f.write(f'\nTotal systematic (quadrature) = {total:.6f}\n')
        f.write(f'Relative systematic = {100.0*total/base_sf:.3f}%\n')

    print(f'Wrote {out_all}')
    print(f'Wrote {out_sel}')
    print(f'Wrote {out_txt}')


if __name__ == '__main__':
    main()
