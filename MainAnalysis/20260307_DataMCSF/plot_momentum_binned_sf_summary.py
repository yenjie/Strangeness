#!/usr/bin/env python3
import csv
import math
import pathlib

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

ROOT = pathlib.Path(__file__).resolve().parents[1]
REPORTS = ROOT / "Reports"


def read_sf_rows(csv_path: pathlib.Path):
    rows = []
    with csv_path.open() as f:
        r = csv.DictReader(f)
        for row in r:
            rows.append(
                {
                    "bin_label": row["bin_label"],
                    "p_min": float(row["p_min"]),
                    "p_max": float(row["p_max"]),
                    "sf": float(row["sf"]),
                    "sf_err": float(row["sf_err"]),
                }
            )
    return rows


def read_systematic_by_bin(csv_path: pathlib.Path):
    rows = list(csv.DictReader(csv_path.open()))
    by_bin = {}
    for r in rows:
        b = r["bin_label"]
        by_bin.setdefault(b, []).append(r)

    out = {}
    for b, rr in by_bin.items():
        by_source = {}
        for r in rr:
            if r["source"] == "baseline":
                continue
            d = abs(float(r["delta_sf_vs_baseline"]))
            src = r["source"]
            by_source[src] = max(by_source.get(src, 0.0), d)
        out[b] = math.sqrt(sum(v * v for v in by_source.values()))
    return out


def make_plot(sf_rows, syst_by_bin, title, x_title, out_pdf):
    fig, ax = plt.subplots(figsize=(8.2, 5.2))

    yvals = []
    for row in sf_rows:
        x0, x1 = row["p_min"], row["p_max"]
        y = row["sf"]
        ystat = row["sf_err"]
        ysyst = syst_by_bin.get(row["bin_label"], 0.0)
        yvals.extend([y - ysyst, y + ysyst, y - ystat, y + ystat])

        rect = Rectangle(
            (x0, y - ysyst),
            x1 - x0,
            2.0 * ysyst,
            facecolor="#FFD84D",
            edgecolor="none",
            alpha=0.45,
            zorder=1,
        )
        ax.add_patch(rect)

        xc = 0.5 * (x0 + x1)
        ax.errorbar(
            [xc],
            [y],
            yerr=[ystat],
            fmt="o",
            color="black",
            markersize=4.5,
            capsize=3,
            lw=1.2,
            zorder=3,
        )
        ax.hlines(y, x0, x1, colors="black", lw=1.0, zorder=2)

    ax.set_title(title)
    ax.set_xlabel(x_title)
    ax.set_ylabel("Scale Factor")
    xmin = min(r["p_min"] for r in sf_rows)
    xmax = max(r["p_max"] for r in sf_rows)
    ax.set_xlim(xmin, xmax)
    if yvals:
        ymin = min(yvals)
        ymax = max(yvals)
        pad = 0.08 * (ymax - ymin if ymax > ymin else 1.0)
        ax.set_ylim(ymin - pad, ymax + pad)

    legend_handles = [
        Rectangle((0, 0), 1, 1, facecolor="#FFD84D", edgecolor="none", alpha=0.45, label="Systematic uncertainty"),
    ]
    ax.legend(handles=legend_handles, loc="best", frameon=True)
    ax.grid(alpha=0.25, linestyle="--", linewidth=0.6)
    fig.tight_layout()
    fig.savefig(out_pdf)
    plt.close(fig)


def main():
    phi_sf = read_sf_rows(REPORTS / "step23_phi_momentum_binned_sf.csv")
    phi_syst = read_systematic_by_bin(REPORTS / "step23_phi_momentum_binned_systematics_qc_curated.csv")
    make_plot(
        phi_sf,
        phi_syst,
        r"$\phi\to K^+K^-$ Momentum-Binned SF",
        r"Kaon momentum $p_{K}$ [GeV]",
        REPORTS / "step23_phi_momentum_binned_sf_summary_plot.pdf",
    )

    ks_sf = read_sf_rows(REPORTS / "step56_kshort_momentum_binned_sf_quadgauss.csv")
    ks_syst = read_systematic_by_bin(REPORTS / "step56_kshort_momentum_binned_systematics_quadgauss_qc_curated.csv")
    make_plot(
        ks_sf,
        ks_syst,
        r"$K^0_S\to\pi^+\pi^-$ Momentum-Binned SF",
        r"Pion momentum $p_{\pi}$ [GeV]",
        REPORTS / "step56_kshort_momentum_binned_sf_summary_plot.pdf",
    )

    print("Wrote", REPORTS / "step23_phi_momentum_binned_sf_summary_plot.pdf")
    print("Wrote", REPORTS / "step56_kshort_momentum_binned_sf_summary_plot.pdf")


if __name__ == "__main__":
    main()
