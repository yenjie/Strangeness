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
    "KPI_DNDY_TOY_COVERAGE_CSV",
    os.environ.get(
        "KPI_TOY_COVERAGE_CSV",
        "result/20260314/dndy_toy_coverage/toy_coverage_bins.csv",
    ),
)
DNDY_DIR = os.environ.get("KPI_DNDY_DIR", "output/systematics_20260314_dndy")
OUT_DIR = os.environ.get("KPI_DNDY_TOP_PLOTS_DIR", "result/20260314/top_plots_dndy")
LOGO_PATH = os.environ.get("EEA_LOGO_PATH", "assets/eea_logo.png")

PYTHIA8_TRUTH_ROOT = os.environ.get(
    "KPI_DNDY_PYTHIA8_TRUTH_ROOT",
    "result/20260314/pythia8_truth/pythia8_zpole_truth_400k_dndy.root",
)
PYTHIA8_ROPE_TRUTH_ROOT = os.environ.get(
    "KPI_DNDY_PYTHIA8_ROPE_TRUTH_ROOT",
    "result/20260314/pythia8_rope_truth/pythia8_zpole_truth_rope_400k_dndy.root",
)
DIRE_TRUTH_ROOT = os.environ.get(
    "KPI_DNDY_DIRE_TRUTH_ROOT",
    "result/20260314/pythia8_dire_truth/pythia8_zpole_truth_dire_400k_dndy.root",
)
HERWIG_TRUTH_ROOT = os.environ.get(
    "KPI_DNDY_HERWIG_TRUTH_ROOT",
    "result/20260314/herwig_truth/herwig_ee_zpole_truth_400k_dndy.root",
)
SHERPA_TRUTH_ROOT = os.environ.get(
    "KPI_DNDY_SHERPA_TRUTH_ROOT",
    "result/20260314/sherpa_truth/sherpa_ee_zpole_truth_400k_dndy.root",
)

PYTHIA8_LABEL = "PYTHIA v8.317, Tune:ee=7"
PYTHIA8_ROPE_LABEL = "PYTHIA v8.317, DIPSY flavour ropes"
DIRE_LABEL = "PYTHIA v8.315, Dire shower"
HERWIG_LABEL = "HERWIG v7.3.0 cluster"
SHERPA_LABEL = "SHERPA v3.0.3, CSS+Ahadic"


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
    out = {}
    if not os.path.exists(TOY_COVERAGE_CSV):
        return out
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
    g.SetMarkerStyle(0)
    return g


def normalize_hist(h):
    hn = h.Clone(h.GetName() + "_norm")
    hn.SetDirectory(0)
    if hn.Integral() > 0:
        hn.Scale(1.0 / hn.Integral())
    return hn


def style_hist(h, color, style=1, width=3):
    h.SetLineColor(color)
    h.SetLineStyle(style)
    h.SetLineWidth(width)
    h.SetMarkerSize(0)


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


def load_curve(root_path, color, style, width=3):
    if not os.path.exists(root_path):
        return None
    h = clone_hist(root_path, "hKPiTruth_dNdY")
    return hist_to_curve_graph(h, color, width, style)


def load_shape(root_path, color, style, width=3):
    if not os.path.exists(root_path):
        return None
    h = normalize_hist(clone_hist(root_path, "hNChY05Truth"))
    style_hist(h, color, style, width)
    return h


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    nominal_path = os.path.join(DNDY_DIR, "nominal_unfold_dndy.root")

    h_data_raw = clone_hist(nominal_path, "hRatioDataBayes_dNdY")
    h_mc_true = clone_hist(nominal_path, "hRatioMcTrue_dNdY")
    h_closure = clone_hist(nominal_path, "hClosureBayes_dNdY")
    h_prior_diff = clone_hist(nominal_path, "hBayesPriorVariationDiff_dNdY")
    h_iter_diff = clone_hist(nominal_path, "hBayesIterVariationDiff_dNdY")

    toy_stat = load_toy_stat_errors("dNdY")
    toy_used = apply_stat_override(h_data_raw, toy_stat)

    h_data = apply_residual_correction(h_data_raw, h_closure)
    h_data.Scale(SF_RATIO)

    var_data = []
    for name in ("npt10", "npt14"):
        path = os.path.join(DNDY_DIR, f"{name}_unfold_dndy.root")
        var_hist = apply_residual_correction(clone_hist(path, "hRatioDataBayes_dNdY"), clone_hist(path, "hClosureBayes_dNdY"))
        var_hist.Scale(SF_RATIO)
        var_data.append(var_hist)

    comp_data = envelope_sys(h_data, var_data)
    total_sys_data = [0.0] * h_data.GetNbinsX()
    for ib in range(1, h_data.GetNbinsX() + 1):
        d_nom = h_data.GetBinContent(ib)
        cval = h_closure.GetBinContent(ib)
        rel_residual = 0.5 * abs((1.0 / cval) - 1.0) if abs(cval) > 1e-12 else 0.0
        abs_residual = rel_residual * abs(d_nom)
        rel_unf = 0.0
        mc_ratio = h_mc_true.GetBinContent(ib)
        if abs(mc_ratio) > 1e-12:
            scale = (1.0 / abs(cval)) if abs(cval) > 1e-12 else 0.0
            rel_unf = max(abs(h_prior_diff.GetBinContent(ib)) * scale, abs(h_iter_diff.GetBinContent(ib)) * scale) / abs(mc_ratio)
        abs_unf = rel_unf * abs(d_nom)
        pid_sf = abs(d_nom) * SF_RATIO_REL_ERR
        total_sys_data[ib - 1] = math.sqrt(comp_data[ib - 1] ** 2 + abs_residual ** 2 + abs_unf ** 2 + pid_sf ** 2)

    g_data_stat, g_data_sys = build_graph_with_sys(h_data, total_sys_data, ROOT.kBlack, 20)
    g_mc_true_curve = hist_to_curve_graph(h_mc_true, ROOT.kRed + 1)

    g_py8 = load_curve(PYTHIA8_TRUTH_ROOT, ROOT.kBlue + 2, 7)
    g_py8_rope = load_curve(PYTHIA8_ROPE_TRUTH_ROOT, ROOT.kCyan + 2, 3)
    g_dire = load_curve(DIRE_TRUTH_ROOT, ROOT.kMagenta + 1, 9)
    g_herwig = load_curve(HERWIG_TRUTH_ROOT, ROOT.kGreen + 2, 4)
    g_sherpa = load_curve(SHERPA_TRUTH_ROOT, ROOT.kOrange + 7, 5)

    c1 = ROOT.TCanvas("cDNdYGenerators", "", 900, 700)
    frame = h_data.Clone("frame_dndy_generators")
    frame.Reset()
    frame.SetTitle("")
    frame.SetMinimum(0.08)
    frame.SetMaximum(0.26)
    frame.GetXaxis().SetLimits(0.0, 31.0)
    frame.GetXaxis().SetTitle("dN_{ch}/dy (|y_{T}|<0.5)")
    frame.GetYaxis().SetTitle("K/#pi")
    style_frame(frame)
    frame.Draw("AXIS")

    g_data_sys.Draw("E2 SAME")
    g_data_stat.Draw("P SAME")
    g_mc_true_curve.Draw("L SAME")
    for g in [g_py8, g_py8_rope, g_dire, g_herwig, g_sherpa]:
        if g:
            g.Draw("L SAME")

    leg = ROOT.TLegend(0.20, 0.64, 0.80, 0.90)
    style_legend(leg)
    leg.SetNColumns(2)
    leg.SetColumnSeparation(0.030)
    leg.SetMargin(0.20)
    leg.SetTextSize(0.025)
    leg.AddEntry(g_data_stat, "DELPHI data (unfolded, stat)", "p")
    leg.AddEntry(g_data_sys, "DELPHI data systematic", "f")
    leg.AddEntry(g_mc_true_curve, "DELPHI nominal PYTHIA8 gen-level", "l")
    if g_py8:
        leg.AddEntry(g_py8, PYTHIA8_LABEL, "l")
    if g_py8_rope:
        leg.AddEntry(g_py8_rope, PYTHIA8_ROPE_LABEL, "l")
    if g_dire:
        leg.AddEntry(g_dire, DIRE_LABEL, "l")
    if g_herwig:
        leg.AddEntry(g_herwig, HERWIG_LABEL, "l")
    if g_sherpa:
        leg.AddEntry(g_sherpa, SHERPA_LABEL, "l")
    leg.Draw()

    header = ROOT.TLatex()
    header.SetNDC()
    header.SetTextAlign(13)
    header.SetTextSize(0.040)
    header.DrawLatex(0.14, 0.94, "K/#pi vs thrust-axis dN_{ch}/dy")
    header.SetTextSize(0.028)
    header.DrawLatex(0.14, 0.895, "Stat: toy-calibrated" if toy_used else "Stat: Bayes diagonal propagation")
    header.DrawLatex(0.14, 0.860, "PYTHIA6 and X-SCAPE omitted here: current local outputs lack particle-level four-vectors")
    draw_logo(c1)
    c1.SaveAs(os.path.join(OUT_DIR, "KtoPi_vs_dNdY_DELPHI_vs_Generators.pdf"))
    c1.SaveAs(os.path.join(OUT_DIR, "KtoPi_vs_dNdY_DELPHI_vs_Generators.png"))

    h_ref = load_shape(PYTHIA8_TRUTH_ROOT, ROOT.kBlue + 2, 7)
    h_rope = load_shape(PYTHIA8_ROPE_TRUTH_ROOT, ROOT.kCyan + 2, 3)
    h_dire = load_shape(DIRE_TRUTH_ROOT, ROOT.kMagenta + 1, 9)
    h_herwig = load_shape(HERWIG_TRUTH_ROOT, ROOT.kGreen + 2, 4)
    h_sherpa = load_shape(SHERPA_TRUTH_ROOT, ROOT.kOrange + 7, 5)
    if h_ref is None:
        raise RuntimeError("Missing standalone PYTHIA8 dN/dy truth shape")
    h_ref.SetLineColor(ROOT.kBlack)
    h_ref.SetLineStyle(1)
    h_ref.SetLineWidth(4)

    c2 = ROOT.TCanvas("cGenDNdY", "", 860, 820)
    pad1 = ROOT.TPad("pad1_dndy", "", 0.0, 0.30, 1.0, 1.0)
    pad2 = ROOT.TPad("pad2_dndy", "", 0.0, 0.0, 1.0, 0.30)
    pad1.SetBottomMargin(0.015)
    pad1.SetLeftMargin(0.13)
    pad1.SetRightMargin(0.04)
    pad2.SetTopMargin(0.02)
    pad2.SetBottomMargin(0.30)
    pad2.SetLeftMargin(0.13)
    pad2.SetRightMargin(0.04)
    pad1.Draw()
    pad2.Draw()

    pad1.cd()
    frame2 = h_ref.Clone("frameGenDNdY")
    frame2.Reset()
    frame2.SetTitle(";dN_{ch}/dy (|y_{T}|<0.5);Normalized events")
    frame2.SetMinimum(0.0)
    ymax = max(h.GetMaximum() for h in [h_ref, h_rope, h_dire, h_herwig, h_sherpa] if h is not None)
    frame2.SetMaximum(ymax * 1.35)
    frame2.GetXaxis().SetLimits(0.0, 31.0)
    frame2.GetXaxis().SetLabelSize(0)
    frame2.GetYaxis().SetTitleSize(0.058)
    frame2.GetYaxis().SetLabelSize(0.048)
    frame2.GetYaxis().SetTitleOffset(1.05)
    frame2.Draw("HIST")
    h_ref.Draw("HIST SAME")
    for h in [h_rope, h_dire, h_herwig, h_sherpa]:
        if h:
            h.Draw("HIST SAME")

    leg2 = ROOT.TLegend(0.44, 0.61, 0.88, 0.88)
    style_legend(leg2)
    leg2.SetTextSize(0.029)
    leg2.AddEntry(h_ref, PYTHIA8_LABEL + " (reference)", "l")
    if h_rope:
        leg2.AddEntry(h_rope, PYTHIA8_ROPE_LABEL, "l")
    if h_dire:
        leg2.AddEntry(h_dire, DIRE_LABEL, "l")
    if h_herwig:
        leg2.AddEntry(h_herwig, HERWIG_LABEL, "l")
    if h_sherpa:
        leg2.AddEntry(h_sherpa, SHERPA_LABEL, "l")
    leg2.Draw()
    label = ROOT.TLatex()
    label.SetNDC()
    label.SetTextSize(0.032)
    label.DrawLatex(0.15, 0.92, "Truth-level thrust-axis activity shape")

    pad2.cd()
    ratio_frame = frame2.Clone("ratioFrameGenDNdY")
    ratio_frame.Reset()
    ratio_frame.SetTitle("")
    ratio_frame.SetMinimum(0.7)
    ratio_frame.SetMaximum(1.3)
    ratio_frame.GetXaxis().SetTitle("dN_{ch}/dy (|y_{T}|<0.5)")
    ratio_frame.GetYaxis().SetTitle("Gen/Ref")
    ratio_frame.GetXaxis().SetTitleSize(0.11)
    ratio_frame.GetXaxis().SetLabelSize(0.10)
    ratio_frame.GetYaxis().SetTitleSize(0.10)
    ratio_frame.GetYaxis().SetLabelSize(0.09)
    ratio_frame.GetYaxis().SetTitleOffset(0.55)
    ratio_frame.GetYaxis().SetNdivisions(505)
    ratio_frame.GetXaxis().SetLimits(0.0, 31.0)
    ratio_frame.Draw("HIST")

    line = ROOT.TLine(0.0, 1.0, 31.0, 1.0)
    line.SetLineStyle(2)
    line.SetLineWidth(2)
    line.Draw()

    keep = [h_ref, h_rope, h_dire, h_herwig, h_sherpa]
    for src in [h_rope, h_dire, h_herwig, h_sherpa]:
        if src is None:
            continue
        ratio = src.Clone(src.GetName() + "_ratio")
        ratio.SetDirectory(0)
        ratio.Divide(h_ref)
        ratio.Draw("HIST SAME")
        keep.append(ratio)

    c2.SaveAs(os.path.join(OUT_DIR, "Generator_dNdY_Comparison.pdf"))
    c2.SaveAs(os.path.join(OUT_DIR, "Generator_dNdY_Comparison.png"))

    fout = ROOT.TFile.Open(os.path.join(OUT_DIR, "KtoPi_vs_dNdY_GeneratorComparison.root"), "RECREATE")
    for obj in [h_data, h_mc_true, g_data_stat, g_data_sys, g_mc_true_curve, h_ref] + keep:
        if obj:
            obj.Write()
    for obj in [g_py8, g_py8_rope, g_dire, g_herwig, g_sherpa]:
        if obj:
            obj.Write()
    fout.Close()
    print("Wrote", os.path.join(OUT_DIR, "KtoPi_vs_dNdY_DELPHI_vs_Generators.pdf"))
    print("Wrote", os.path.join(OUT_DIR, "Generator_dNdY_Comparison.pdf"))


if __name__ == "__main__":
    main()
