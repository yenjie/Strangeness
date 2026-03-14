#!/usr/bin/env python3
import csv
import math
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import ROOT

ROOT.gROOT.SetBatch(True)

BASE = os.path.abspath(os.path.dirname(__file__))
DATE = "20260314"
OUTDIR = os.path.join(BASE, f"result/{DATE}/dndy_toy_coverage")
ROOT_PATH = os.path.join(BASE, "output/systematics_20260314_dndy/nominal_unfold_dndy.root")
NTOYS = int(os.environ.get("KPI_DNDY_NTOYS", "400"))

os.makedirs(OUTDIR, exist_ok=True)


def clone_hist(path, name):
    f = ROOT.TFile.Open(path, "READ")
    if not f or f.IsZombie():
        raise RuntimeError(f"Cannot open {path}")
    h = f.Get(name)
    if not h:
        raise RuntimeError(f"Missing {name} in {path}")
    out = h.Clone(f"{name}_{os.path.basename(path).replace('.', '_')}")
    out.SetDirectory(0)
    f.Close()
    return out


def hist_arrays(h):
    n = h.GetNbinsX()
    centers = np.array([h.GetXaxis().GetBinCenter(i) for i in range(1, n + 1)], dtype=float)
    vals = np.array([h.GetBinContent(i) for i in range(1, n + 1)], dtype=float)
    errs = np.array([h.GetBinError(i) for i in range(1, n + 1)], dtype=float)
    return centers, vals, errs


def matrix_array(h2):
    nx = h2.GetNbinsX()
    ny = h2.GetNbinsY()
    arr = np.zeros((nx, ny), dtype=float)
    for ix in range(1, nx + 1):
        for iy in range(1, ny + 1):
            arr[ix - 1, iy - 1] = h2.GetBinContent(ix, iy)
    return arr


def build_ratio(num, num_err, den, den_err):
    ratio = np.zeros_like(num, dtype=float)
    err = np.zeros_like(num, dtype=float)
    mask = den > 0
    ratio[mask] = num[mask] / den[mask]
    mask2 = mask & (num > 0)
    err[mask2] = ratio[mask2] * np.sqrt(
        (num_err[mask2] / num[mask2]) ** 2 + (den_err[mask2] / den[mask2]) ** 2
    )
    return ratio, err


def iterative_bayes_unfold(meas, response_true_reco, prior, niter):
    ntrue, nreco = response_true_reco.shape
    prior_vec = np.clip(np.asarray(prior, dtype=float), 0.0, None)
    if prior_vec.sum() <= 0:
        prior_vec = np.full(ntrue, 1.0 / ntrue, dtype=float)
    else:
        prior_vec = prior_vec / prior_vec.sum()

    p = np.clip(np.asarray(response_true_reco, dtype=float), 0.0, None)
    unfolded = np.zeros(ntrue, dtype=float)
    for _ in range(niter):
        unfolded[:] = 0.0
        for r in range(nreco):
            mr = max(0.0, meas[r])
            if mr == 0.0:
                continue
            norm = np.dot(p[:, r], prior_vec)
            if norm <= 0.0:
                continue
            weights = (p[:, r] * prior_vec) / norm
            unfolded += weights * mr
        total = np.clip(unfolded, 0.0, None).sum()
        if total <= 0.0:
            break
        prior_vec = np.clip(unfolded, 0.0, None) / total
    return np.clip(unfolded, 0.0, None), np.sqrt(np.clip(unfolded, 0.0, None))


def write_csv(path, rows, fieldnames):
    with open(path, "w", newline="", encoding="ascii") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main():
    rng = np.random.default_rng(12345)

    resp_k = matrix_array(clone_hist(ROOT_PATH, "hDNdYResponseKRebinned"))
    resp_pi = matrix_array(clone_hist(ROOT_PATH, "hDNdYResponsePiRebinned"))
    centers, k_truth, _ = hist_arrays(clone_hist(ROOT_PATH, "hKPrior_dNdY"))
    _, pi_truth, _ = hist_arrays(clone_hist(ROOT_PATH, "hPiPrior_dNdY"))
    _, k_reco, _ = hist_arrays(clone_hist(ROOT_PATH, "hKMcReco"))
    _, pi_reco, _ = hist_arrays(clone_hist(ROOT_PATH, "hPiMcReco"))

    truth_ratio, _ = build_ratio(
        k_truth,
        np.sqrt(np.clip(k_truth, 0.0, None)),
        pi_truth,
        np.sqrt(np.clip(pi_truth, 0.0, None)),
    )

    nbins = len(truth_ratio)
    pulls = np.full((NTOYS, nbins), np.nan, dtype=float)
    ratios = np.full((NTOYS, nbins), np.nan, dtype=float)

    for itoy in range(NTOYS):
        k_meas = rng.poisson(np.clip(k_reco, 0.0, None))
        pi_meas = rng.poisson(np.clip(pi_reco, 0.0, None))
        k_unf, k_err = iterative_bayes_unfold(k_meas, resp_k, k_truth, 1)
        pi_unf, pi_err = iterative_bayes_unfold(pi_meas, resp_pi, pi_truth, 1)
        ratio, ratio_err = build_ratio(k_unf, k_err, pi_unf, pi_err)
        valid = ratio_err > 0
        ratios[itoy, valid] = ratio[valid]
        pulls[itoy, valid] = (ratio[valid] - truth_ratio[valid]) / ratio_err[valid]

    mean_pull = np.nanmean(pulls, axis=0)
    width_pull = np.nanstd(pulls, axis=0)
    cov1 = np.nanmean(np.abs(pulls) < 1.0, axis=0)
    cov2 = np.nanmean(np.abs(pulls) < 2.0, axis=0)
    ratio_mean = np.nanmean(ratios, axis=0)
    ratio_rms_abs = np.nanstd(ratios, axis=0)
    ratio_bias_abs = ratio_mean - truth_ratio
    ratio_rmse_abs = np.sqrt(np.nanmean((ratios - truth_ratio[np.newaxis, :]) ** 2, axis=0))

    summary_rows = [{
        "axis": "dNdY",
        "mean_abs_pull_bias": float(np.nanmean(np.abs(mean_pull))),
        "mean_pull_width": float(np.nanmean(width_pull)),
        "min_cov1": float(np.nanmin(cov1)),
        "min_cov2": float(np.nanmin(cov2)),
        "max_pull_width": float(np.nanmax(width_pull)),
        "n_toys": NTOYS,
    }]
    bin_rows = []
    for ib in range(nbins):
        n_valid = int(np.sum(np.isfinite(ratios[:, ib])))
        bin_rows.append({
            "axis": "dNdY",
            "bin": ib + 1,
            "center": float(centers[ib]),
            "truth_ratio": float(truth_ratio[ib]),
            "toy_mean": float(ratio_mean[ib]) if math.isfinite(ratio_mean[ib]) else 0.0,
            "toy_rms_abs": float(ratio_rms_abs[ib]) if math.isfinite(ratio_rms_abs[ib]) else 0.0,
            "toy_bias_abs": float(ratio_bias_abs[ib]) if math.isfinite(ratio_bias_abs[ib]) else 0.0,
            "toy_rmse_abs": float(ratio_rmse_abs[ib]) if math.isfinite(ratio_rmse_abs[ib]) else 0.0,
            "pull_mean": float(mean_pull[ib]) if math.isfinite(mean_pull[ib]) else 0.0,
            "pull_width": float(width_pull[ib]) if math.isfinite(width_pull[ib]) else 0.0,
            "cov1": float(cov1[ib]) if math.isfinite(cov1[ib]) else 0.0,
            "cov2": float(cov2[ib]) if math.isfinite(cov2[ib]) else 0.0,
            "n_valid": n_valid,
        })

    write_csv(os.path.join(OUTDIR, "toy_coverage_summary.csv"), summary_rows, list(summary_rows[0].keys()))
    write_csv(os.path.join(OUTDIR, "toy_coverage_bins.csv"), bin_rows, list(bin_rows[0].keys()))

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.6), constrained_layout=True)
    axes[0].axhline(0.0, color="black", linestyle="--", linewidth=1.5)
    axes[0].plot(centers, mean_pull, marker="o")
    axes[0].set_title(r"$dN_{ch}/dy$ toy pull mean")
    axes[0].set_xlabel(r"True $dN_{ch}/dy$ ($|y_T|<0.5$)")
    axes[0].set_ylabel("Mean pull")
    axes[0].grid(alpha=0.25)

    axes[1].axhline(1.0, color="black", linestyle="--", linewidth=1.5)
    axes[1].plot(centers, width_pull, marker="o", label="Pull width")
    axes[1].plot(centers, cov1, marker="s", label="Coverage |pull|<1")
    axes[1].plot(centers, cov2, marker="^", label="Coverage |pull|<2")
    axes[1].set_title(r"$dN_{ch}/dy$ toy pull width / coverage")
    axes[1].set_xlabel(r"True $dN_{ch}/dy$ ($|y_T|<0.5$)")
    axes[1].set_ylabel("Value")
    axes[1].grid(alpha=0.25)
    axes[1].legend(frameon=False, fontsize=9)

    fig.savefig(os.path.join(OUTDIR, "DNdYToyCoverage_Summary.pdf"))
    fig.savefig(os.path.join(OUTDIR, "DNdYToyCoverage_Summary.png"), dpi=160)
    plt.close(fig)

    print("Wrote", os.path.join(OUTDIR, "toy_coverage_bins.csv"))
    print("Wrote", os.path.join(OUTDIR, "DNdYToyCoverage_Summary.pdf"))


if __name__ == "__main__":
    main()
