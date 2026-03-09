#!/usr/bin/env python3
import csv
import math
import pathlib
import re
import subprocess

ROOT = pathlib.Path(__file__).resolve().parents[1]
REPORTS = ROOT / "Reports"

CHI2_MAX = 1.8


def run_root(cmd: str, timeout_sec: int = 180):
    p = subprocess.run(cmd, shell=True, cwd=ROOT, capture_output=True, text=True, timeout=timeout_sec)
    if p.returncode != 0:
        raise RuntimeError(f"Command failed:\n{cmd}\nSTDOUT:\n{p.stdout}\nSTDERR:\n{p.stderr}")


def parse_summary(path: pathlib.Path):
    s = path.read_text()
    m_eff = re.search(r"tagging_efficiency,([0-9.eE+-]+),tagging_efficiency_err,([0-9.eE+-]+)", s)
    if not m_eff:
        raise RuntimeError(f"Cannot parse efficiency from {path}")
    eff, eff_err = float(m_eff.group(1)), float(m_eff.group(2))

    cats = {}
    for c in ["0tag", "1tag", "2tag"]:
        m = re.search(
            rf"{c},model,([^,]+),chi2ndf,([0-9.eE+-]+),mean,([0-9.eE+-]+),meanErr,([0-9.eE+-]+),yield,([0-9.eE+-]+),yieldErr,([0-9.eE+-]+)",
            s,
        )
        if not m:
            raise RuntimeError(f"Cannot parse {c} metrics from {path}")
        cats[c] = {
            "model": m.group(1),
            "chi2": float(m.group(2)),
            "yield": float(m.group(5)),
            "yield_err": float(m.group(6)),
        }
    return eff, eff_err, cats


def pass_quality(cats):
    bad = []
    for c in ["1tag", "2tag"]:
        if cats[c]["chi2"] > CHI2_MAX:
            bad.append(f"{c}:chi2={cats[c]['chi2']:.3f}")
    return len(bad) == 0, "; ".join(bad)


def run_step23_mombin(prefix: str, seed_csv: str, signal_model: str, bkg_mode: str,
                      mass_min: float, mass_max: float, fit_min: float, fit_max: float, n_bins: int,
                      p_min: float, p_max: float):
    mc_prefix = f"{prefix}_mc"
    da_prefix = f"{prefix}_data"
    mc_summary = REPORTS / f"{mc_prefix}_fit_summary.txt"
    da_summary = REPORTS / f"{da_prefix}_fit_summary.txt"

    if mc_summary.exists() and da_summary.exists():
        emc, semc, cmc = parse_summary(mc_summary)
        eda, seda, cda = parse_summary(da_summary)
        sf = eda / emc
        sf_err = sf * math.sqrt((seda / eda) ** 2 + (semc / emc) ** 2)
        q_mc_ok, q_mc_msg = pass_quality(cmc)
        q_da_ok, q_da_msg = pass_quality(cda)
        q_ok = q_mc_ok and q_da_ok
        q_msg = "; ".join(x for x in [q_mc_msg, q_da_msg] if x)
        return emc, semc, eda, seda, sf, sf_err, q_ok, q_msg

    cmd_mc = (
        "root -l -b -q 'Code/step2_mc_random_pairs.C("
        f"\"Samples/merged_mc_v2.2.root\",\"{seed_csv}\",\"Reports\",\"{mc_prefix}\",\"MC\","
        f"\"{signal_model}\",\"{bkg_mode}\",{mass_min:.3f},{mass_max:.3f},{fit_min:.3f},{fit_max:.3f},{n_bins},"
        f"{p_min:.3f},{p_max:.3f})'"
    )
    cmd_da = (
        "root -l -b -q 'Code/step2_mc_random_pairs.C("
        f"\"Samples/merged_data_v2.2.root\",\"{seed_csv}\",\"Reports\",\"{da_prefix}\",\"Data\","
        f"\"{signal_model}\",\"{bkg_mode}\",{mass_min:.3f},{mass_max:.3f},{fit_min:.3f},{fit_max:.3f},{n_bins},"
        f"{p_min:.3f},{p_max:.3f})'"
    )
    run_root(cmd_mc)
    run_root(cmd_da)

    emc, semc, cmc = parse_summary(mc_summary)
    eda, seda, cda = parse_summary(da_summary)

    sf = eda / emc
    sf_err = sf * math.sqrt((seda / eda) ** 2 + (semc / emc) ** 2)

    q_mc_ok, q_mc_msg = pass_quality(cmc)
    q_da_ok, q_da_msg = pass_quality(cda)
    q_ok = q_mc_ok and q_da_ok
    q_msg = "; ".join(x for x in [q_mc_msg, q_da_msg] if x)

    return emc, semc, eda, seda, sf, sf_err, q_ok, q_msg


def write_summary(rows, out_txt: pathlib.Path, baseline_by_bin):
    grouped = {}
    for r in rows:
        grouped.setdefault(r["bin_label"], []).append(r)

    with out_txt.open("w") as f:
        f.write("Phi->KK momentum-binned SF systematics\n")
        f.write(f"Quality thresholds: chi2/ndf<{CHI2_MAX} for 1tag/2tag in both MC and Data\n\n")
        for bin_label in sorted(grouped.keys()):
            pmin = grouped[bin_label][0]["p_min"]
            pmax = grouped[bin_label][0]["p_max"]
            base_sf, base_sf_err = baseline_by_bin[bin_label]
            f.write(f"{bin_label} [{pmin:.2f},{pmax:.2f}) GeV\n")
            f.write(f"  Baseline SF = {base_sf:.6f} +/- {base_sf_err:.6f}\n")
            by_source = {}
            for row in grouped[bin_label]:
                if row["source"] == "baseline":
                    continue
                if not math.isfinite(row["delta_sf_vs_baseline"]):
                    continue
                src = row["source"]
                by_source[src] = max(by_source.get(src, 0.0), abs(row["delta_sf_vs_baseline"]))
            total = math.sqrt(sum(v * v for v in by_source.values()))
            for src in ["signal_model", "background_model", "fit_window_binning", "matching_angle"]:
                f.write(f"  {src}: max |delta SF| = {by_source.get(src, 0.0):.6f}\n")
            rel = 100.0 * total / base_sf if base_sf > 0 else 0.0
            f.write(f"  Total systematic (quadrature) = {total:.6f}\n")
            f.write(f"  Relative systematic = {rel:.3f}%\n\n")


def main():
    bins = [
        (0.0, 1.2, "p00_12"),
        (1.2, 2.2, "p12_22"),
        (2.2, 5.0, "p22_50"),
    ]

    variations = [
        ("baseline", "baseline", "auto", "baseline", (0.99, 1.06, 1.00, 1.05, 280), "Reports/step1_best_fit_params.csv"),

        ("signal_model", "sig_gauss", "Gauss", "baseline", (0.99, 1.06, 1.00, 1.05, 280), "Reports/step1_best_fit_params.csv"),
        ("signal_model", "sig_doublegauss", "DoubleGauss", "baseline", (0.99, 1.06, 1.00, 1.05, 280), "Reports/step1_best_fit_params.csv"),

        ("background_model", "bkg_allpol2", "auto", "all_pol2", (0.99, 1.06, 1.00, 1.05, 280), "Reports/step1_best_fit_params.csv"),
        ("background_model", "bkg_allpol3", "auto", "all_pol3", (0.99, 1.06, 1.00, 1.05, 280), "Reports/step1_best_fit_params.csv"),

        ("fit_window_binning", "fwbin_narrow240", "auto", "baseline", (0.99, 1.06, 1.005, 1.045, 240), "Reports/step1_best_fit_params.csv"),
        ("fit_window_binning", "fwbin_alt320", "auto", "baseline", (0.99, 1.06, 1.00, 1.05, 320), "Reports/step1_best_fit_params.csv"),
        ("fit_window_binning", "fwbin_wide200", "auto", "baseline", (0.99, 1.06, 0.995, 1.055, 200), "Reports/step1_best_fit_params.csv"),

        ("matching_angle", "match_angle0025", "auto", "baseline", (0.99, 1.06, 1.00, 1.05, 280), "Reports/step1_a0025_syst_best_fit_params.csv"),
    ]

    run_root(
        "root -l -b -q 'Code/step1_make_report.C("
        "\"Samples/merged_mc_v2.2.root\",\"Reports\",0.025,\"step1_a0025_syst\")'"
    )

    rows = []
    baseline_by_bin = {}

    for pmin, pmax, blabel in bins:
        for source, name, sig, bkg, pars, seed in variations:
            prefix = f"step23_pbin_{blabel}_{name}"
            print(f"[phi] {blabel} :: {name}")
            try:
                emc, semc, eda, seda, sf, sf_err, q_ok, q_msg = run_step23_mombin(
                    prefix, seed, sig, bkg, *pars, pmin, pmax
                )
            except Exception as e:
                if source == "baseline":
                    raise RuntimeError(f"Baseline failed in bin {blabel}: {e}") from e
                emc = semc = eda = seda = sf = sf_err = float("nan")
                q_ok = False
                q_msg = f"fit_failed: {str(e).splitlines()[0][:160]}"

            if source == "baseline":
                baseline_by_bin[blabel] = (sf, sf_err)
                delta = 0.0
            else:
                delta = sf - baseline_by_bin[blabel][0] if math.isfinite(sf) else float("nan")
            rows.append({
                "bin_label": blabel,
                "p_min": pmin,
                "p_max": pmax,
                "source": source,
                "variation": name,
                "signal_model": sig,
                "bkg_mode": bkg,
                "mass_min": pars[0],
                "mass_max": pars[1],
                "fit_min": pars[2],
                "fit_max": pars[3],
                "n_bins": pars[4],
                "eff_mc": emc,
                "eff_mc_err": semc,
                "eff_data": eda,
                "eff_data_err": seda,
                "sf": sf,
                "sf_err": sf_err,
                "delta_sf_vs_baseline": delta,
                "quality_pass": int(q_ok),
                "quality_notes": q_msg,
            })

    out_all = REPORTS / "step23_phi_momentum_binned_systematics_qc_all.csv"
    with out_all.open("w", newline="") as f:
        w = csv.DictWriter(
            f,
            fieldnames=[
                "bin_label", "p_min", "p_max",
                "source", "variation", "signal_model", "bkg_mode",
                "mass_min", "mass_max", "fit_min", "fit_max", "n_bins",
                "eff_mc", "eff_mc_err", "eff_data", "eff_data_err",
                "sf", "sf_err", "delta_sf_vs_baseline",
                "quality_pass", "quality_notes",
            ],
        )
        w.writeheader()
        w.writerows(rows)

    curated = [r for r in rows if r["source"] == "baseline" or r["quality_pass"] == 1]
    out_cur = REPORTS / "step23_phi_momentum_binned_systematics_qc_curated.csv"
    with out_cur.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(curated)

    write_summary(
        curated,
        REPORTS / "step23_phi_momentum_binned_systematics_qc_curated_summary.txt",
        baseline_by_bin,
    )

    print(f"Wrote {out_all}")
    print(f"Wrote {out_cur}")
    print(f"Wrote {REPORTS / 'step23_phi_momentum_binned_systematics_qc_curated_summary.txt'}")


if __name__ == "__main__":
    main()
