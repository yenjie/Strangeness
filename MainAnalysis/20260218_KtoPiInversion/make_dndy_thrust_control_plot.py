#!/usr/bin/env python3
import os

import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetEndErrorSize(0)

INPUT_DIR = os.environ.get("KPI_DNDY_INPUT_DIR", "output/systematics_20260314_dndy_inputs")
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


def normalize(h):
    total = h.Integral(1, h.GetNbinsX())
    if total > 0.0:
        h.Scale(1.0 / total)
    return h


def ratio_hist(num, den, name):
    out = num.Clone(name)
    out.SetDirectory(0)
    for ib in range(1, out.GetNbinsX() + 1):
        n = num.GetBinContent(ib)
        en = num.GetBinError(ib)
        d = den.GetBinContent(ib)
        ed = den.GetBinError(ib)
        if d <= 0.0:
            out.SetBinContent(ib, 0.0)
            out.SetBinError(ib, 0.0)
            continue
        val = n / d
        rel2 = 0.0
        if n > 0.0:
            rel2 += (en / n) ** 2
        if d > 0.0:
            rel2 += (ed / d) ** 2
        out.SetBinContent(ib, val)
        out.SetBinError(ib, abs(val) * rel2 ** 0.5 if rel2 > 0.0 else 0.0)
    return out


def style_axis(h, upper):
    if upper:
        h.GetXaxis().SetLabelSize(0.0)
        h.GetXaxis().SetTitleSize(0.0)
        h.GetYaxis().SetTitle("Unit-normalized events")
        h.GetYaxis().SetTitleSize(0.060)
        h.GetYaxis().SetLabelSize(0.048)
        h.GetYaxis().SetTitleOffset(1.00)
    else:
        h.GetXaxis().SetTitle("dN_{ch}/dy (|y_{T}|<0.5)")
        h.GetYaxis().SetTitle("Ratio to MC reco")
        h.GetXaxis().SetTitleSize(0.115)
        h.GetYaxis().SetTitleSize(0.100)
        h.GetXaxis().SetLabelSize(0.100)
        h.GetYaxis().SetLabelSize(0.085)
        h.GetXaxis().SetTitleOffset(0.98)
        h.GetYaxis().SetTitleOffset(0.52)
        h.GetYaxis().SetNdivisions(505)


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
    data_path = os.path.join(INPUT_DIR, "nominal_data.root")
    mc_path = os.path.join(INPUT_DIR, "nominal_mc.root")

    h_data = clone_hist(data_path, "hDNdYReco")
    h_mc_reco = clone_hist(mc_path, "hDNdYReco")
    h_mc_gen = clone_hist(mc_path, "hDNdYTrue")
    normalize(h_data)
    normalize(h_mc_reco)
    normalize(h_mc_gen)

    h_data.SetMarkerStyle(20)
    h_data.SetMarkerSize(1.15)
    h_data.SetMarkerColor(ROOT.kBlack)
    h_data.SetLineColor(ROOT.kBlack)
    h_data.SetLineWidth(2)

    h_mc_reco.SetLineColor(ROOT.kRed + 1)
    h_mc_reco.SetLineWidth(3)
    h_mc_reco.SetFillColorAlpha(ROOT.kRed + 1, 0.16)

    h_mc_gen.SetLineColor(ROOT.kBlue + 2)
    h_mc_gen.SetLineWidth(3)
    h_mc_gen.SetLineStyle(2)

    ymax = max(h_data.GetMaximum(), h_mc_reco.GetMaximum(), h_mc_gen.GetMaximum())

    c = ROOT.TCanvas("cDNdYControl", "dN/dy wrt thrust axis", 880, 860)
    pad1 = ROOT.TPad("pad1", "", 0.0, 0.31, 1.0, 1.0)
    pad2 = ROOT.TPad("pad2", "", 0.0, 0.0, 1.0, 0.31)
    pad1.SetBottomMargin(0.02)
    pad1.SetLeftMargin(0.13)
    pad1.SetRightMargin(0.05)
    pad1.SetTopMargin(0.08)
    pad2.SetTopMargin(0.03)
    pad2.SetBottomMargin(0.34)
    pad2.SetLeftMargin(0.13)
    pad2.SetRightMargin(0.05)
    pad1.Draw()
    pad2.Draw()

    pad1.cd()
    frame = h_mc_reco.Clone("frame_dndy_control")
    frame.Reset()
    frame.SetTitle("")
    frame.SetMinimum(0.0)
    frame.SetMaximum(ymax * 1.35)
    style_axis(frame, True)
    frame.Draw("AXIS")
    h_mc_reco.Draw("HIST SAME")
    h_mc_gen.Draw("HIST SAME")
    h_data.Draw("E1 SAME")

    leg = ROOT.TLegend(0.50, 0.69, 0.89, 0.89)
    leg.SetBorderSize(0)
    leg.SetFillColor(0)
    leg.SetFillStyle(0)
    leg.SetLineColor(0)
    leg.SetTextSize(0.040)
    leg.AddEntry(h_data, "Data reco", "lep")
    leg.AddEntry(h_mc_reco, "PYTHIA8 v2.5 MC reco", "lf")
    leg.AddEntry(h_mc_gen, "PYTHIA8 v2.5 MC gen", "l")
    leg.Draw()

    header = ROOT.TLatex()
    header.SetNDC()
    header.SetTextAlign(13)
    header.SetTextSize(0.040)
    header.DrawLatex(0.14, 0.94, "Charged-particle dN_{ch}/dy wrt thrust axis")
    header.SetTextSize(0.028)
    header.DrawLatex(0.14, 0.895, "Selection: PassAll, charged tracks, |y_{T}|<0.5")
    draw_logo(c)

    pad2.cd()
    h_ratio_data = ratio_hist(h_data, h_mc_reco, "hRatioDataToMcReco_dNdY")
    h_ratio_gen = ratio_hist(h_mc_gen, h_mc_reco, "hRatioMcGenToReco_dNdY")
    h_ratio_data.SetMarkerStyle(20)
    h_ratio_data.SetMarkerSize(1.0)
    h_ratio_data.SetMarkerColor(ROOT.kBlack)
    h_ratio_data.SetLineColor(ROOT.kBlack)
    h_ratio_data.SetLineWidth(2)
    h_ratio_gen.SetLineColor(ROOT.kBlue + 2)
    h_ratio_gen.SetLineWidth(3)
    h_ratio_gen.SetLineStyle(2)

    rframe = h_ratio_data.Clone("rframe_dndy_control")
    rframe.Reset()
    rframe.SetMinimum(0.6)
    rframe.SetMaximum(1.4)
    style_axis(rframe, False)
    rframe.Draw("AXIS")
    line = ROOT.TLine(rframe.GetXaxis().GetXmin(), 1.0, rframe.GetXaxis().GetXmax(), 1.0)
    line.SetLineStyle(2)
    line.SetLineWidth(2)
    line.Draw("SAME")
    h_ratio_gen.Draw("HIST SAME")
    h_ratio_data.Draw("E1 SAME")

    out_pdf = os.path.join(OUT_DIR, "ControlPlot_DNdYThrust_DataMC.pdf")
    out_png = os.path.join(OUT_DIR, "ControlPlot_DNdYThrust_DataMC.png")
    c.SaveAs(out_pdf)
    c.SaveAs(out_png)

    fout = ROOT.TFile.Open(os.path.join(OUT_DIR, "ControlPlot_DNdYThrust_DataMC.root"), "RECREATE")
    for obj in [h_data, h_mc_reco, h_mc_gen, h_ratio_data, h_ratio_gen]:
        obj.Write()
    fout.Close()
    print("Wrote", out_pdf)


if __name__ == "__main__":
    main()
