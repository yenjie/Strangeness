#!/usr/bin/env python3
import csv
import math
import os

import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetEndErrorSize(0)

SF_K = 0.9101
SF_PI = 0.9571
SF_K_ERR = 0.052118
SF_PI_ERR = 0.035532
SF_RATIO = SF_K / SF_PI
SF_RATIO_REL_ERR = math.sqrt((SF_K_ERR / SF_K) ** 2 + (SF_PI_ERR / SF_PI) ** 2)

TOY_COVERAGE_CSV = os.environ.get(
    "KPI_TOY_COVERAGE_CSV",
    "result/20260312/reviewer_followup_validation/toy_coverage/toy_coverage_bins.csv",
)
DNDY_DIR = os.environ.get("KPI_DNDY_DIR", "output/systematics_20260314_dndy")
OUT_DIR = os.environ.get("KPI_DNDY_TOP_PLOTS_DIR", "result/20260314/top_plots_dndy")
LOGO_PATH = os.environ.get("EEA_LOGO_PATH", "assets/eea_logo.png")


def clone_hist(path, name):
    f = ROOT.TFile.Open(path, "READ")
    if not f or f.IsZombie():
        raise RuntimeError(f"Cannot open {path}")
    h = f.Get(name)
    if h is None:
        raise RuntimeError(f"Missing {name} in {path}")
    hc = h.Clone(f"{name}_{os.path.basename(path).replace('.', '_')}")
    hc.SetDirectory(0)
    f.Close()
    return hc


def load_toy_stat_errors(axis_name):
    if not os.path.exists(TOY_COVERAGE_CSV):
        return {}
    out = {}
    with open(TOY_COVERAGE_CSV, newline="", encoding="ascii") as f:
        for row in csv.DictReader(f):
            if row["axis"] != axis_name:
                continue
            value = row.get("toy_rmse_abs", row.get("toy_rms_abs", "0"))
            out[int(row["bin"])] = float(value)
    return out


def apply_stat_override(h, stat_map):
    if not stat_map:
        return False
    for ib in range(1, h.GetNbinsX() + 1):
        if ib in stat_map:
            h.SetBinError(ib, stat_map[ib])
    return True


def apply_residual_correction(h, h_closure):
    hc = h.Clone(h.GetName() + "_residcorr")
    hc.SetDirectory(0)
    for ib in range(1, hc.GetNbinsX() + 1):
        cval = h_closure.GetBinContent(ib)
        if abs(cval) > 1e-12:
            hc.SetBinContent(ib, hc.GetBinContent(ib) / cval)
            hc.SetBinError(ib, hc.GetBinError(ib) / abs(cval))
    return hc


def envelope_sys(nominal, var_hists):
    out = [0.0] * nominal.GetNbinsX()
    for ib in range(1, nominal.GetNbinsX() + 1):
        nom = nominal.GetBinContent(ib)
        out[ib - 1] = max(abs(var.GetBinContent(ib) - nom) for var in var_hists) if var_hists else 0.0
    return out


def build_graph_with_sys(h_nom, total_sys, color, marker):
    n = h_nom.GetNbinsX()
    g_stat = ROOT.TGraphErrors(n)
    g_sys = ROOT.TGraphErrors(n)
    idx = 0
    for ib in range(1, n + 1):
        y = h_nom.GetBinContent(ib)
        ey = h_nom.GetBinError(ib)
        if y == 0.0 and ey == 0.0:
            continue
        x = h_nom.GetXaxis().GetBinCenter(ib)
        ex = 0.45 * h_nom.GetXaxis().GetBinWidth(ib)
        es = total_sys[ib - 1]
        g_stat.SetPoint(idx, x, y)
        g_stat.SetPointError(idx, ex, ey)
        g_sys.SetPoint(idx, x, y)
        g_sys.SetPointError(idx, ex, es)
        idx += 1
    g_stat.Set(idx)
    g_sys.Set(idx)
    g_stat.SetMarkerStyle(marker)
    g_stat.SetMarkerSize(1.2)
    g_stat.SetMarkerColor(color)
    g_stat.SetLineColor(color)
    g_stat.SetLineWidth(2)
    g_sys.SetFillColorAlpha(color, 0.20)
    g_sys.SetLineColor(color)
    g_sys.SetLineWidth(1)
    g_sys.SetMarkerStyle(1)
    return g_stat, g_sys


def hist_to_curve_graph(h, color, width=3, style=1):
    g = ROOT.TGraph()
    j = 0
    for ib in range(1, h.GetNbinsX() + 1):
        y = h.GetBinContent(ib)
        if y == 0.0:
            continue
        x = h.GetXaxis().GetBinCenter(ib)
        g.SetPoint(j, x, y)
        j += 1
    g.SetLineColor(color)
    g.SetLineWidth(width)
    g.SetLineStyle(style)
    return g


def style_frame(h):
    h.GetXaxis().SetTitleSize(0.052)
    h.GetYaxis().SetTitleSize(0.052)
    h.GetXaxis().SetLabelSize(0.044)
    h.GetYaxis().SetLabelSize(0.044)
    h.GetXaxis().SetTitleOffset(1.05)
    h.GetYaxis().SetTitleOffset(1.15)


def style_legend(leg):
    leg.SetBorderSize(0)
    leg.SetFillColor(0)
    leg.SetFillStyle(0)
    leg.SetLineColor(0)


def draw_logo(canvas, x1=0.14, y1=0.79, x2=0.29, y2=0.93):
    if not os.path.exists(LOGO_PATH):
        return
    canvas.cd()
    pad = ROOT.TPad(f"{canvas.GetName()}_logo", "", x1, y1, x2, y2)
    pad.SetFillStyle(4000)
    pad.SetFrameFillStyle(4000)
    pad.SetMargin(0.0, 0.0, 0.0, 0.0)
    pad.Draw()
    pad.cd()
    image = ROOT.TImage.Open(LOGO_PATH)
    if image:
        image.Draw()
        if not hasattr(canvas, "_logo_refs"):
            canvas._logo_refs = []
        canvas._logo_refs.extend([pad, image])
    canvas.cd()


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    nominal_path = os.path.join(DNDY_DIR, "nominal_unfold_dndy.root")

    h_data_raw = clone_hist(nominal_path, "hRatioDataBayes_dNdY")
    h_mc_true = clone_hist(nominal_path, "hRatioMcTrue_dNdY")
    h_mc_bayes_raw = clone_hist(nominal_path, "hRatioMcBayes_dNdY")
    h_double = clone_hist(nominal_path, "hDataOverMcBayes_dNdY")
    h_closure = clone_hist(nominal_path, "hClosureBayes_dNdY")
    h_prior_diff = clone_hist(nominal_path, "hBayesPriorVariationDiff_dNdY")
    h_iter_diff = clone_hist(nominal_path, "hBayesIterVariationDiff_dNdY")

    toy_stat_used = False
    toy_stat = load_toy_stat_errors("dNdY")
    if toy_stat:
        toy_stat_used = apply_stat_override(h_data_raw, toy_stat)
        apply_stat_override(h_mc_bayes_raw, toy_stat)
        apply_stat_override(h_double, toy_stat)

    h_data = apply_residual_correction(h_data_raw, h_closure)
    h_mc_bayes = apply_residual_correction(h_mc_bayes_raw, h_closure)
    h_double = apply_residual_correction(h_double, h_closure)
    h_data.Scale(SF_RATIO)
    h_double.Scale(SF_RATIO)

    var_data = []
    var_mc = []
    for name in ("npt10", "npt14"):
        path = os.path.join(DNDY_DIR, f"{name}_unfold_dndy.root")
        var_data_hist = apply_residual_correction(clone_hist(path, "hRatioDataBayes_dNdY"), h_closure)
        var_data_hist.Scale(SF_RATIO)
        var_mc_hist = apply_residual_correction(clone_hist(path, "hRatioMcBayes_dNdY"), h_closure)
        var_data.append(var_data_hist)
        var_mc.append(var_mc_hist)

    comp_data = envelope_sys(h_data, var_data)
    comp_mc = envelope_sys(h_mc_bayes, var_mc)
    total_sys_data = [0.0] * h_data.GetNbinsX()
    total_sys_mc = [0.0] * h_data.GetNbinsX()

    for ib in range(1, h_data.GetNbinsX() + 1):
        d_nom = h_data.GetBinContent(ib)
        m_nom = h_mc_bayes.GetBinContent(ib)
        r_nom = h_double.GetBinContent(ib)

        cval = h_closure.GetBinContent(ib)
        rel_residual = 0.5 * abs((1.0 / cval) - 1.0) if abs(cval) > 1e-12 else 0.0
        abs_residual_data = rel_residual * abs(d_nom)
        abs_residual_mc = rel_residual * abs(m_nom)

        rel_unf = 0.0
        if abs(r_nom) > 1e-12:
            scale = (1.0 / abs(cval)) if abs(cval) > 1e-12 else 0.0
            rel_unf = max(abs(h_prior_diff.GetBinContent(ib)) * scale, abs(h_iter_diff.GetBinContent(ib)) * scale) / abs(r_nom)
        abs_unf_data = rel_unf * abs(d_nom)
        abs_unf_mc = rel_unf * abs(m_nom)

        pid_sf_data = abs(d_nom) * SF_RATIO_REL_ERR
        total_sys_data[ib - 1] = math.sqrt(comp_data[ib - 1] ** 2 + abs_residual_data ** 2 + abs_unf_data ** 2 + pid_sf_data ** 2)
        total_sys_mc[ib - 1] = math.sqrt(comp_mc[ib - 1] ** 2 + abs_residual_mc ** 2 + abs_unf_mc ** 2)

    g_data_stat, g_data_sys = build_graph_with_sys(h_data, total_sys_data, ROOT.kBlack, 20)
    g_mc_stat, g_mc_sys = build_graph_with_sys(h_mc_bayes, total_sys_mc, ROOT.kAzure + 2, 26)
    g_mc_true_curve = hist_to_curve_graph(h_mc_true, ROOT.kRed + 1)

    ymax = 0.0
    for ib in range(1, h_data.GetNbinsX() + 1):
        ymax = max(
            ymax,
            h_data.GetBinContent(ib) + h_data.GetBinError(ib) + total_sys_data[ib - 1],
            h_mc_bayes.GetBinContent(ib) + h_mc_bayes.GetBinError(ib) + total_sys_mc[ib - 1],
            h_mc_true.GetBinContent(ib),
        )

    c = ROOT.TCanvas("cDataMcDNdY", "", 900, 700)
    frame = h_data.Clone("frame_dndy_data_mc")
    frame.Reset()
    frame.SetTitle("")
    frame.SetMinimum(0.08)
    frame.SetMaximum(max(0.18, ymax * 1.30))
    frame.GetXaxis().SetLimits(0.0, 31.0)
    frame.GetXaxis().SetTitle("dN_{ch}/dy (|y_{T}|<0.5)")
    frame.GetYaxis().SetTitle("K/#pi")
    style_frame(frame)
    frame.Draw("AXIS")
    g_mc_sys.Draw("E2 SAME")
    g_data_sys.Draw("E2 SAME")
    g_mc_true_curve.Draw("L SAME")
    g_mc_stat.Draw("P SAME")
    g_data_stat.Draw("P SAME")

    leg = ROOT.TLegend(0.22, 0.70, 0.59, 0.89)
    style_legend(leg)
    leg.SetTextSize(0.028)
    leg.AddEntry(g_data_stat, "Data unfolded (stat)", "p")
    leg.AddEntry(g_data_sys, "Data systematic", "f")
    leg.AddEntry(g_mc_true_curve, "PYTHIA8 v2.5 truth", "l")
    leg.AddEntry(g_mc_stat, "PYTHIA8 v2.5 unfolded", "p")
    leg.AddEntry(g_mc_sys, "PYTHIA8 v2.5 systematic", "f")
    leg.Draw()

    header = ROOT.TLatex()
    header.SetNDC()
    header.SetTextAlign(13)
    header.SetTextSize(0.040)
    header.DrawLatex(0.14, 0.94, "K/#pi vs thrust-axis dN_{ch}/dy")
    header.SetTextSize(0.028)
    stat_label = "Stat: toy-calibrated" if toy_stat_used else "Stat: Bayes diagonal propagation"
    header.DrawLatex(0.14, 0.895, stat_label)
    draw_logo(c)

    out_pdf = os.path.join(OUT_DIR, "KtoPi_vs_dNdY_DataMC_with_Systematics.pdf")
    out_png = os.path.join(OUT_DIR, "KtoPi_vs_dNdY_DataMC_with_Systematics.png")
    c.SaveAs(out_pdf)
    c.SaveAs(out_png)

    fout = ROOT.TFile.Open(os.path.join(OUT_DIR, "KtoPi_vs_dNdY_Comparison.root"), "RECREATE")
    for obj in [
        h_data,
        h_mc_true,
        h_mc_bayes,
        h_double,
        h_closure,
        h_prior_diff,
        h_iter_diff,
        g_data_stat,
        g_data_sys,
        g_mc_stat,
        g_mc_sys,
        g_mc_true_curve,
    ]:
        obj.Write()
    fout.Close()
    print("Wrote", out_pdf)


if __name__ == "__main__":
    main()
