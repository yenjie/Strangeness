#!/usr/bin/env python3
import csv
import math
import pathlib
import re
import subprocess

ROOT = pathlib.Path(__file__).resolve().parents[1]
REPORTS = ROOT / "Reports"


def run_root(cmd: str):
    p = subprocess.run(cmd, shell=True, cwd=ROOT, capture_output=True, text=True)
    if p.returncode != 0:
        raise RuntimeError(f"Command failed:\n{cmd}\nSTDOUT:\n{p.stdout}\nSTDERR:\n{p.stderr}")


def read_eff(summary_path: pathlib.Path):
    s = summary_path.read_text()
    m = re.search(r"tagging_efficiency,([0-9.eE+-]+),tagging_efficiency_err,([0-9.eE+-]+)", s)
    if not m:
        raise RuntimeError(f"Could not parse efficiency in {summary_path}")
    return float(m.group(1)), float(m.group(2))


def run_step56(prefix: str, seed_csv: str, signal_model: str, bkg_mode: str,
               mass_min: float, mass_max: float, fit_min: float, fit_max: float, n_bins: int):
    mc_prefix = f"{prefix}_mc"
    data_prefix = f"{prefix}_data"

    cmd_mc = (
        "root -l -b -q 'Code/step5_mc_random_pairs_kshort.C("
        f"\"Samples/merged_mc_v2.2.root\",\"{seed_csv}\",\"Reports\",\"{mc_prefix}\",\"MC\","
        f"\"{signal_model}\",\"{bkg_mode}\",{mass_min:.3f},{mass_max:.3f},{fit_min:.3f},{fit_max:.3f},{n_bins})'"
    )
    cmd_data = (
        "root -l -b -q 'Code/step5_mc_random_pairs_kshort.C("
        f"\"Samples/merged_data_v2.2.root\",\"{seed_csv}\",\"Reports\",\"{data_prefix}\",\"Data\","
        f"\"{signal_model}\",\"{bkg_mode}\",{mass_min:.3f},{mass_max:.3f},{fit_min:.3f},{fit_max:.3f},{n_bins})'"
    )
    run_root(cmd_mc)
    run_root(cmd_data)

    emc, semc = read_eff(REPORTS / f"{mc_prefix}_fit_summary.txt")
    ed, sed = read_eff(REPORTS / f"{data_prefix}_fit_summary.txt")
    sf = ed / emc
    sf_err = sf * math.sqrt((sed / ed) ** 2 + (semc / emc) ** 2)
    return emc, semc, ed, sed, sf, sf_err


def main():
    rows = []

    # Baseline seed from step4 angle<0.025 run (fit 0.4-0.6)
    baseline_seed = "Reports/step4_a0025_best_fit_params.csv"

    # Baseline step5/6 configuration
    base = run_step56(
        prefix="syst_baseline",
        seed_csv=baseline_seed,
        signal_model="auto",
        bkg_mode="baseline",
        mass_min=0.40,
        mass_max=0.70,
        fit_min=0.40,
        fit_max=0.70,
        n_bins=280,
    )
    base_sf = base[4]
    rows.append(["baseline", "baseline", *base, 0.0])

    # 1) Signal model choice
    for name, sig in [("sig_gauss", "Gauss"), ("sig_doublegauss", "DoubleGauss")]:
        vals = run_step56(name, baseline_seed, sig, "baseline", 0.40, 0.70, 0.40, 0.70, 280)
        rows.append(["signal_model", name, *vals, vals[4] - base_sf])

    # 2) Background model choice
    for name, bkg in [("bkg_allpol2", "all_pol2"), ("bkg_allpol3", "all_pol3")]:
        vals = run_step56(name, baseline_seed, "auto", bkg, 0.40, 0.70, 0.40, 0.70, 280)
        rows.append(["background_model", name, *vals, vals[4] - base_sf])

    # 3) Fit window and binning
    for name, pars in [
        ("fwbin_narrow240", (0.40, 0.70, 0.42, 0.68, 240)),
        ("fwbin_wide200", (0.40, 0.70, 0.40, 0.70, 200)),
    ]:
        vals = run_step56(name, baseline_seed, "auto", "baseline", *pars)
        rows.append(["fit_window_binning", name, *vals, vals[4] - base_sf])

    # 4) Matching angle variation: regenerate step4 seed with angle<0.01
    run_root(
        "root -l -b -q 'Code/step4_kshort_make_report.C("
        "\"Samples/merged_mc_v2.2.root\",\"Reports\",0.01,\"step4_a0010_syst\",0.40,0.60)'"
    )
    vals = run_step56(
        prefix="match_angle0010",
        seed_csv="Reports/step4_a0010_syst_best_fit_params.csv",
        signal_model="auto",
        bkg_mode="baseline",
        mass_min=0.40,
        mass_max=0.70,
        fit_min=0.40,
        fit_max=0.70,
        n_bins=280,
    )
    rows.append(["matching_angle", "angle_0p01_vs_0p025", *vals, vals[4] - base_sf])

    out_csv = REPORTS / "step56_sf_systematics.csv"
    with out_csv.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow([
            "source", "variation", "eff_mc", "eff_mc_err", "eff_data", "eff_data_err",
            "sf", "sf_err", "delta_sf_vs_baseline"
        ])
        w.writerows(rows)

    # source-wise max absolute delta and total
    by_source = {}
    for r in rows:
        src = r[0]
        if src == "baseline":
            continue
        d = abs(float(r[-1]))
        by_source[src] = max(by_source.get(src, 0.0), d)

    total = math.sqrt(sum(v * v for v in by_source.values()))

    out_txt = REPORTS / "step56_sf_systematics_summary.txt"
    with out_txt.open("w") as f:
        f.write("Step5/6 SF Systematic Uncertainty Summary\n")
        f.write("Baseline SF = {:.6f} +/- {:.6f}\n\n".format(base[4], base[5]))
        for k in ["signal_model", "background_model", "fit_window_binning", "matching_angle"]:
            f.write(f"{k}: max |delta SF| = {by_source.get(k,0.0):.6f}\n")
        f.write("\n")
        f.write(f"Total systematic (quadrature) = {total:.6f}\n")
        f.write(f"Relative systematic = {100.0*total/base[4]:.3f}%\n")

    print(f"Wrote {out_csv}")
    print(f"Wrote {out_txt}")


if __name__ == "__main__":
    main()
