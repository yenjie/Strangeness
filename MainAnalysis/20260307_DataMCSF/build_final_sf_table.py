#!/usr/bin/env python3
import csv
import math
import pathlib
import re


ROOT = pathlib.Path(__file__).resolve().parents[1]
REPORTS = ROOT / "Reports"


def read_curated_csv(path: pathlib.Path):
    with path.open() as f:
        rows = list(csv.DictReader(f))
    baseline = next((r for r in rows if r.get("variation") == "baseline"), None)
    if baseline is None:
        raise RuntimeError(f"No baseline row in {path}")
    return {
        "eff_mc": float(baseline["eff_mc"]),
        "eff_mc_err": float(baseline["eff_mc_err"]),
        "eff_data": float(baseline["eff_data"]),
        "eff_data_err": float(baseline["eff_data_err"]),
        "sf": float(baseline["sf"]),
        "sf_stat_err": float(baseline["sf_err"]),
    }


def read_syst(summary_path: pathlib.Path):
    text = summary_path.read_text()
    m_tot = re.search(r"Total systematic \(quadrature\)\s*=\s*([0-9.eE+-]+)", text)
    m_rel = re.search(r"Relative systematic\s*=\s*([0-9.eE+-]+)%", text)
    if not m_tot or not m_rel:
        raise RuntimeError(f"Could not parse systematics from {summary_path}")
    return float(m_tot.group(1)), float(m_rel.group(1))


def main():
    channels = [
        (
            "phi_to_KK",
            REPORTS / "step23_sf_systematics_qc_curated_v3.csv",
            REPORTS / "step23_sf_systematics_qc_curated_v3_summary.txt",
        ),
        (
            "kshort_to_pipi",
            REPORTS / "step56_sf_systematics_qc_curated_v2.csv",
            REPORTS / "step56_sf_systematics_qc_curated_v2_summary.txt",
        ),
    ]

    out_rows = []
    for channel, csv_path, summary_path in channels:
        base = read_curated_csv(csv_path)
        sf_syst_err, rel_syst_percent = read_syst(summary_path)
        sf_total_err = math.sqrt(base["sf_stat_err"] ** 2 + sf_syst_err ** 2)
        out_rows.append({
            "channel": channel,
            "eff_mc": base["eff_mc"],
            "eff_mc_err": base["eff_mc_err"],
            "eff_data": base["eff_data"],
            "eff_data_err": base["eff_data_err"],
            "sf": base["sf"],
            "sf_stat_err": base["sf_stat_err"],
            "sf_syst_err": sf_syst_err,
            "sf_total_err": sf_total_err,
            "rel_syst_percent": rel_syst_percent,
            "curated_csv": csv_path.name,
            "curated_summary": summary_path.name,
        })

    out_csv = REPORTS / "final_scale_factors.csv"
    with out_csv.open("w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "channel",
                "eff_mc",
                "eff_mc_err",
                "eff_data",
                "eff_data_err",
                "sf",
                "sf_stat_err",
                "sf_syst_err",
                "sf_total_err",
                "rel_syst_percent",
                "curated_csv",
                "curated_summary",
            ],
        )
        writer.writeheader()
        writer.writerows(out_rows)

    out_txt = REPORTS / "final_scale_factors_summary.txt"
    with out_txt.open("w") as f:
        f.write("Final Tagging-Efficiency Scale Factors\n")
        for row in out_rows:
            f.write(
                "\n"
                f"{row['channel']}: SF={row['sf']:.6f} +/- {row['sf_stat_err']:.6f} (stat)"
                f" +/- {row['sf_syst_err']:.6f} (syst)"
                f" => +/- {row['sf_total_err']:.6f} (total)\n"
                f"  rel_syst={row['rel_syst_percent']:.3f}%\n"
                f"  source_csv={row['curated_csv']}\n"
                f"  source_summary={row['curated_summary']}\n"
            )

    print(f"Wrote {out_csv}")
    print(f"Wrote {out_txt}")


if __name__ == "__main__":
    main()
