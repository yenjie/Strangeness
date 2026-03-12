#!/usr/bin/env python3
import math
import os

import ROOT
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetEndErrorSize(0)

BASE = os.path.abspath(os.path.dirname(__file__))
OUTDIR = os.path.join(BASE, "result/20260310/unfolding_validation")
os.makedirs(OUTDIR, exist_ok=True)


def get_hist(path, name):
    f = ROOT.TFile.Open(path, "READ")
    if not f or f.IsZombie():
        raise RuntimeError(f"Cannot open {path}")
    h = f.Get(name)
    if not h:
        raise RuntimeError(f"Missing {name} in {path}")
    c = h.Clone(f"{name}_{os.path.basename(path).replace('.', '_')}")
    c.SetDirectory(0)
    f.Close()
    return c


def style_hist(h, color, marker):
    h.SetLineColor(color)
    h.SetMarkerColor(color)
    h.SetMarkerStyle(marker)
    h.SetMarkerSize(1.0)
    h.SetLineWidth(2)


def style_axis_top(h, ytitle):
    h.GetYaxis().SetTitle(ytitle)
    h.GetXaxis().SetLabelSize(0.0)
    h.GetXaxis().SetTitleSize(0.0)
    h.GetYaxis().SetTitleSize(0.07)
    h.GetYaxis().SetLabelSize(0.06)
    h.GetYaxis().SetTitleOffset(0.95)


def style_axis_ratio(h, xtitle):
    h.GetYaxis().SetTitle("Ratio")
    h.GetXaxis().SetTitle(xtitle)
    h.GetYaxis().SetTitleSize(0.10)
    h.GetYaxis().SetLabelSize(0.09)
    h.GetYaxis().SetTitleOffset(0.50)
    h.GetYaxis().SetNdivisions(505)
    h.GetXaxis().SetTitleSize(0.10)
    h.GetXaxis().SetLabelSize(0.09)
    h.GetXaxis().SetTitleOffset(0.95)


def build_ratio_hist(num, den, name):
    h = num.Clone(name)
    h.SetDirectory(0)
    h.Divide(den)
    return h


def draw_species_refold_column(
    canvas,
    x1,
    x2,
    measured_mc,
    refolded_mc,
    measured_data,
    refolded_data,
    title,
    xtitle,
    ratio_min=0.90,
    ratio_max=1.10,
):
    top = ROOT.TPad(f"top_{title}_{x1}", "", x1, 0.30, x2, 1.0)
    bot = ROOT.TPad(f"bot_{title}_{x1}", "", x1, 0.00, x2, 0.30)
    top.SetBottomMargin(0.02)
    top.SetTopMargin(0.08)
    top.SetLeftMargin(0.13 if x1 < 0.1 else 0.10)
    top.SetRightMargin(0.05 if x2 > 0.9 else 0.02)
    bot.SetTopMargin(0.02)
    bot.SetBottomMargin(0.34)
    bot.SetLeftMargin(top.GetLeftMargin())
    bot.SetRightMargin(top.GetRightMargin())
    top.Draw()
    bot.Draw()

    top.cd()
    top.SetTicks(1, 1)
    mmc = measured_mc.Clone(f"{measured_mc.GetName()}_{title}_top")
    rmc = refolded_mc.Clone(f"{refolded_mc.GetName()}_{title}_top")
    mdata = measured_data.Clone(f"{measured_data.GetName()}_{title}_top")
    rdata = refolded_data.Clone(f"{refolded_data.GetName()}_{title}_top")
    style_hist(mmc, ROOT.kBlack, 20)
    style_hist(rmc, ROOT.kAzure + 1, 24)
    style_hist(mdata, ROOT.kRed + 1, 21)
    style_hist(rdata, ROOT.kGreen + 2, 25)
    max_y = 1.25 * max(mmc.GetMaximum(), rmc.GetMaximum(), mdata.GetMaximum(), rdata.GetMaximum())
    min_y = 0.85 * min_positive(mmc, rmc, mdata, rdata)
    mmc.SetMinimum(min_y)
    mmc.SetMaximum(max_y)
    mmc.SetTitle(title)
    style_axis_top(mmc, "Species yield in reco space")
    mmc.Draw("E1")
    rmc.Draw("E1 SAME")
    mdata.Draw("E1 SAME")
    rdata.Draw("E1 SAME")
    lab = ROOT.TLatex()
    lab.SetNDC()
    lab.SetTextSize(0.055)
    lab.DrawLatex(0.15, 0.92, title)
    leg = ROOT.TLegend(0.18, 0.68, 0.90, 0.89)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.042)
    leg.SetNColumns(2)
    leg.AddEntry(mmc, "MC measured reco", "lep")
    leg.AddEntry(rmc, "MC refolded", "lep")
    leg.AddEntry(mdata, "Data measured reco", "lep")
    leg.AddEntry(rdata, "Data refolded", "lep")
    leg.Draw()

    bot.cd()
    bot.SetTicks(1, 1)
    rmc_ratio = build_ratio_hist(refolded_mc, measured_mc, f"{refolded_mc.GetName()}_{title}_ratio_mc")
    rdata_ratio = build_ratio_hist(refolded_data, measured_data, f"{refolded_data.GetName()}_{title}_ratio_data")
    style_hist(rmc_ratio, ROOT.kAzure + 1, 20)
    style_hist(rdata_ratio, ROOT.kRed + 1, 21)
    rmc_ratio.SetMinimum(ratio_min)
    rmc_ratio.SetMaximum(ratio_max)
    rmc_ratio.SetTitle("")
    style_axis_ratio(rmc_ratio, xtitle)
    rmc_ratio.Draw("E1")
    line = ROOT.TLine(rmc_ratio.GetXaxis().GetXmin(), 1.0, rmc_ratio.GetXaxis().GetXmax(), 1.0)
    line.SetLineStyle(2)
    line.SetLineWidth(2)
    line.Draw("SAME")
    rmc_ratio.Draw("E1 SAME")
    rdata_ratio.Draw("E1 SAME")
    canvas.cd()
    return {
        "mc_rms": rms_distance_from_unity(rmc_ratio),
        "data_rms": rms_distance_from_unity(rdata_ratio),
    }, [top, bot, mmc, rmc, mdata, rdata, lab, leg, rmc_ratio, rdata_ratio, line]


def make_refolding_figure(root_path, xtitle, out_name, names, ratio_min=0.90, ratio_max=1.10):
    c = ROOT.TCanvas(f"c_refold_{out_name}", "", 1400, 700)
    keep = []
    metrics_k, keep_k = draw_species_refold_column(
        c, 0.00, 0.50,
        get_hist(root_path, names["k_measured_mc"]),
        get_hist(root_path, names["k_refold_mc"]),
        get_hist(root_path, names["k_measured_data"]),
        get_hist(root_path, names["k_refold_data"]),
        "Kaon refolding", xtitle,
        ratio_min=ratio_min,
        ratio_max=ratio_max,
    )
    keep.extend(keep_k)
    metrics_pi, keep_pi = draw_species_refold_column(
        c, 0.50, 1.00,
        get_hist(root_path, names["pi_measured_mc"]),
        get_hist(root_path, names["pi_refold_mc"]),
        get_hist(root_path, names["pi_measured_data"]),
        get_hist(root_path, names["pi_refold_data"]),
        "Pion refolding", xtitle,
        ratio_min=ratio_min,
        ratio_max=ratio_max,
    )
    keep.extend(keep_pi)
    c._keepalive = keep
    c.Update()
    c.SaveAs(os.path.join(OUTDIR, out_name + ".pdf"))
    c.SaveAs(os.path.join(OUTDIR, out_name + ".png"))

    return {
        "K": metrics_k,
        "Pi": metrics_pi,
    }


def make_stress_figure(root_path_ntag, root_path_dndeta):
    c = ROOT.TCanvas("c_stress", "", 1400, 700)
    keep = []
    keep.extend(draw_stress_column(
        c, 0.00, 0.50, root_path_ntag,
        "hRatioStressTruth", "hRatioStressUnfold", "hStressClosure",
        "Injected shape test: N_{ch}^{tag}", "True N_{ch}^{tag}",
    ))
    keep.extend(draw_stress_column(
        c, 0.50, 1.00, root_path_dndeta,
        "hRatioStressTruth_dNdEta", "hRatioStressUnfold_dNdEta", "hStressClosure_dNdEta",
        "Injected shape test: dN_{ch}/d#eta", "True dN_{ch}/d#eta (|#eta|<0.5)",
    ))
    c._keepalive = keep
    c.Update()
    c.SaveAs(os.path.join(OUTDIR, "UnfoldingStressTest_Combined.pdf"))
    c.SaveAs(os.path.join(OUTDIR, "UnfoldingStressTest_Combined.png"))


def draw_stress_column(canvas, x1, x2, root_path, truth_name, unfold_name, ratio_name, title, xtitle):
    truth = get_hist(root_path, truth_name)
    unfold = get_hist(root_path, unfold_name)
    ratio = get_hist(root_path, ratio_name)

    top = ROOT.TPad(f"top_stress_{title}_{x1}", "", x1, 0.30, x2, 1.0)
    bot = ROOT.TPad(f"bot_stress_{title}_{x1}", "", x1, 0.00, x2, 0.30)
    top.SetBottomMargin(0.02)
    top.SetTopMargin(0.08)
    top.SetLeftMargin(0.13 if x1 < 0.1 else 0.10)
    top.SetRightMargin(0.05 if x2 > 0.9 else 0.02)
    bot.SetTopMargin(0.02)
    bot.SetBottomMargin(0.34)
    bot.SetLeftMargin(top.GetLeftMargin())
    bot.SetRightMargin(top.GetRightMargin())
    top.Draw()
    bot.Draw()

    top.cd()
    style_hist(truth, ROOT.kBlack, 20)
    style_hist(unfold, ROOT.kBlue + 1, 24)
    truth.SetMinimum(0.85 * min_positive(truth, unfold))
    truth.SetMaximum(1.25 * max(truth.GetMaximum(), unfold.GetMaximum()))
    truth.SetTitle(title)
    style_axis_top(truth, "K/#pi on true axis")
    truth.Draw("E1")
    unfold.Draw("E1 SAME")
    lab = ROOT.TLatex()
    lab.SetNDC()
    lab.SetTextSize(0.055)
    lab.DrawLatex(0.15, 0.92, title)
    leg = ROOT.TLegend(0.18, 0.74, 0.88, 0.89)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.050)
    leg.AddEntry(truth, "Injected truth", "lep")
    leg.AddEntry(unfold, "Unfolded result", "lep")
    leg.Draw()

    bot.cd()
    style_hist(ratio, ROOT.kBlue + 1, 20)
    ratio.SetMinimum(0.85)
    ratio.SetMaximum(1.15)
    ratio.SetTitle("")
    style_axis_ratio(ratio, xtitle)
    ratio.Draw("E1")
    line = ROOT.TLine(ratio.GetXaxis().GetXmin(), 1.0, ratio.GetXaxis().GetXmax(), 1.0)
    line.SetLineStyle(2)
    line.SetLineWidth(2)
    line.Draw("SAME")
    ratio.Draw("E1 SAME")
    canvas.cd()
    return [top, bot, truth, unfold, ratio, lab, leg, line]


def min_positive(*hists):
    values = []
    for h in hists:
        for ib in range(1, h.GetNbinsX() + 1):
            val = h.GetBinContent(ib)
            if val > 0:
                values.append(val)
    return min(values) if values else 0.0


def rms_distance_from_unity(h):
    vals = []
    for ib in range(1, h.GetNbinsX() + 1):
        y = h.GetBinContent(ib)
        ey = h.GetBinError(ib)
        if y == 0.0 and ey == 0.0:
            continue
        vals.append((y - 1.0) ** 2)
    return math.sqrt(sum(vals) / len(vals)) if vals else 0.0


def build_migration_metrics(response_path, hist_k, hist_pi, axis_tag):
    hk = get_hist(response_path, hist_k)
    hpi = get_hist(response_path, hist_pi)
    rows = []
    for label, h in [("K", hk), ("Pi", hpi)]:
        nb = min(h.GetNbinsX(), h.GetNbinsY())
        for ib in range(1, nb + 1):
            diag = h.GetBinContent(ib, ib)
            sum_true = 0.0
            sum_reco = 0.0
            for iy in range(1, h.GetNbinsY() + 1):
                sum_true += h.GetBinContent(ib, iy)
            for ix in range(1, h.GetNbinsX() + 1):
                sum_reco += h.GetBinContent(ix, ib)
            stability = diag / sum_true if sum_true > 0 else 0.0
            purity = diag / sum_reco if sum_reco > 0 else 0.0
            rows.append({
                "axis": axis_tag,
                "species": label,
                "bin": ib,
                "center": h.GetXaxis().GetBinCenter(ib),
                "purity": purity,
                "stability": stability,
            })
    return rows


def make_migration_metrics_plot(rows):
    grouped = {}
    for row in rows:
        grouped.setdefault(row["axis"], {}).setdefault(row["species"], []).append(row)

    fig, axes = plt.subplots(1, 2, figsize=(14, 6), constrained_layout=True)
    style_map = {
        ("K", "purity"): dict(color="#1f77b4", marker="o", linestyle="-", label="K purity"),
        ("K", "stability"): dict(color="#1f77b4", marker="s", linestyle="--", label="K stability"),
        ("Pi", "purity"): dict(color="#d62728", marker="o", linestyle="-", label="Pi purity"),
        ("Pi", "stability"): dict(color="#d62728", marker="s", linestyle="--", label="Pi stability"),
    }
    for ax, axis in zip(axes, ["Ntag", "dNdEta"]):
        axis_rows = grouped[axis]
        for species in ["K", "Pi"]:
            xs = [row["center"] for row in axis_rows[species]]
            for metric in ["purity", "stability"]:
                ys = [row[metric] for row in axis_rows[species]]
                ax.plot(xs, ys, **style_map[(species, metric)])
        ax.set_ylim(0.0, 1.05)
        ax.set_ylabel("Fraction")
        ax.set_xlabel("True $N_{ch}^{tag}$" if axis == "Ntag" else "True $dN_{ch}/d\\eta$ ($|\\eta|<0.5$)")
        ax.set_title(f"{axis} migration metrics")
        ax.grid(alpha=0.25)
        ax.legend(frameon=False, ncols=2, fontsize=10)
    fig.savefig(os.path.join(OUTDIR, "MigrationMetrics_PurityStability.pdf"))
    fig.savefig(os.path.join(OUTDIR, "MigrationMetrics_PurityStability.png"), dpi=150)
    plt.close(fig)


def write_status(rows, refold_ntag, refold_dndeta, stress_ntag_path, stress_dndeta_path):
    stress_ntag = get_hist(stress_ntag_path, "hStressClosure")
    stress_dndeta = get_hist(stress_dndeta_path, "hStressClosure_dNdEta")
    lines = [
        "# 2026-03-10 iteration 4 unfolding-validation status",
        "",
        "Implemented in this pass:",
        "- refolding validation for `Ntag` unfolding",
        "- refolding validation for `dNch/deta` unfolding",
        "- mismatched-shape closure stress test for `Ntag` and `dNch/deta`",
        "- migration purity/stability summary for K and pi responses",
        "",
        "Quantitative summary:",
        f"- `Ntag` refolding RMS(refold/reco - 1): K MC = {refold_ntag['K']['mc_rms']:.4f}, K data = {refold_ntag['K']['data_rms']:.4f}, Pi MC = {refold_ntag['Pi']['mc_rms']:.4f}, Pi data = {refold_ntag['Pi']['data_rms']:.4f}",
        f"- `dNch/deta` refolding RMS(refold/reco - 1): K MC = {refold_dndeta['K']['mc_rms']:.4f}, K data = {refold_dndeta['K']['data_rms']:.4f}, Pi MC = {refold_dndeta['Pi']['mc_rms']:.4f}, Pi data = {refold_dndeta['Pi']['data_rms']:.4f}",
        f"- `Ntag` mismatched-shape RMS(unfolded/injected - 1) = {rms_distance_from_unity(stress_ntag):.4f}",
        f"- `dNch/deta` mismatched-shape RMS(unfolded/injected - 1) = {rms_distance_from_unity(stress_dndeta):.4f}",
        "",
        "Pending reviewer items:",
        "- pseudo-experiment / bootstrap coverage study",
        "- explicit SVD `kReg` stability scan and singular-value summary",
        "",
        "Generated outputs:",
        "- `result/20260310/unfolding_validation/NtagUnfolding_RefoldingValidation.pdf`",
        "- `result/20260310/unfolding_validation/DNdEtaUnfolding_RefoldingValidation.pdf`",
        "- `result/20260310/unfolding_validation/UnfoldingStressTest_Combined.pdf`",
        "- `result/20260310/unfolding_validation/MigrationMetrics_PurityStability.pdf`",
    ]
    out_md = os.path.join(BASE, "report/20260310_iteration4_unfolding_validation_status.md")
    with open(out_md, "w", encoding="ascii") as f:
        f.write("\n".join(lines) + "\n")

    out_csv = os.path.join(OUTDIR, "migration_metrics.csv")
    with open(out_csv, "w", encoding="ascii") as f:
        f.write("axis,species,bin,center,purity,stability\n")
        for row in rows:
            f.write(
                f"{row['axis']},{row['species']},{row['bin']},{row['center']:.6f},"
                f"{row['purity']:.6f},{row['stability']:.6f}\n"
            )


def main():
    ntag_root = os.path.join(BASE, "output/NtagUnfolding_BayesSVD.root")
    dndeta_root = os.path.join(BASE, "output/systematics_20260306_dndeta/nominal_unfold_dndeta.root")

    refold_ntag = make_refolding_figure(
        ntag_root,
        "N_{ch}^{tag,reco}",
        "NtagUnfolding_RefoldingValidation",
        {
            "k_measured_mc": "hKMcReco",
            "k_refold_mc": "hKMcBayesRefold",
            "k_measured_data": "hKDataReco",
            "k_refold_data": "hKDataBayesRefold",
            "pi_measured_mc": "hPiMcReco",
            "pi_refold_mc": "hPiMcBayesRefold",
            "pi_measured_data": "hPiDataReco",
            "pi_refold_data": "hPiDataBayesRefold",
        },
        ratio_min=0.70,
        ratio_max=1.30,
    )
    refold_dndeta = make_refolding_figure(
        dndeta_root,
        "dN_{ch}^{reco}/d#eta (|#eta|<0.5)",
        "DNdEtaUnfolding_RefoldingValidation",
        {
            "k_measured_mc": "hKMcReco",
            "k_refold_mc": "hKMcBayesRefold_dNdEta",
            "k_measured_data": "hKDataReco",
            "k_refold_data": "hKDataBayesRefold_dNdEta",
            "pi_measured_mc": "hPiMcReco",
            "pi_refold_mc": "hPiMcBayesRefold_dNdEta",
            "pi_measured_data": "hPiDataReco",
            "pi_refold_data": "hPiDataBayesRefold_dNdEta",
        },
        ratio_min=0.00,
        ratio_max=1.20,
    )
    make_stress_figure(ntag_root, dndeta_root)

    rows = []
    rows.extend(build_migration_metrics(ntag_root, "hNtagResponseK", "hNtagResponsePi", "Ntag"))
    rows.extend(build_migration_metrics(dndeta_root, "hDNdEtaResponseKRebinned", "hDNdEtaResponsePiRebinned", "dNdEta"))
    make_migration_metrics_plot(rows)
    write_status(rows, refold_ntag, refold_dndeta, ntag_root, dndeta_root)
    print("Wrote unfolding validation package to", OUTDIR)


if __name__ == "__main__":
    main()
