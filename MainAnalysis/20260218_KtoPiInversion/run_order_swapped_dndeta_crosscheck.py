#!/usr/bin/env python3
import math
import os
import sys

import ROOT
import numpy as np

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetEndErrorSize(0)

SF_K = 0.9101
SF_PI = 0.9571
SF_K_ERR = 0.052118
SF_PI_ERR = 0.035532
SF_RATIO = SF_K / SF_PI
KEEP_BINS = 8
N_ITER_BAYES = 1


def clone_empty_like(h, name):
    out = h.Clone(name)
    out.SetDirectory(0)
    out.Reset()
    out.Sumw2()
    return out


def collapse_tail_1d(src, keep_bins, name):
    keep_bins = max(1, min(keep_bins, src.GetNbinsX()))
    out = ROOT.TH1D(name, src.GetTitle(), keep_bins, src.GetXaxis().GetXmin(), src.GetXaxis().GetBinUpEdge(keep_bins))
    out.SetDirectory(0)
    out.Sumw2()
    for ib in range(1, src.GetNbinsX() + 1):
        b = min(ib, keep_bins)
        c = src.GetBinContent(ib)
        e = src.GetBinError(ib)
        out.SetBinContent(b, out.GetBinContent(b) + c)
        out.SetBinError(b, math.sqrt(out.GetBinError(b) ** 2 + e ** 2))
    return out


def collapse_tail_2d(src, keep_x, keep_y, name):
    keep_x = max(1, min(keep_x, src.GetNbinsX()))
    keep_y = max(1, min(keep_y, src.GetNbinsY()))
    out = ROOT.TH2D(name, src.GetTitle(), keep_x, src.GetXaxis().GetXmin(), src.GetXaxis().GetBinUpEdge(keep_x), keep_y, src.GetYaxis().GetXmin(), src.GetYaxis().GetBinUpEdge(keep_y))
    out.SetDirectory(0)
    out.Sumw2()
    for ix in range(1, src.GetNbinsX() + 1):
        bx = min(ix, keep_x)
        for iy in range(1, src.GetNbinsY() + 1):
            by = min(iy, keep_y)
            c = src.GetBinContent(ix, iy)
            e = src.GetBinError(ix, iy)
            out.SetBinContent(bx, by, out.GetBinContent(bx, by) + c)
            out.SetBinError(bx, by, math.sqrt(out.GetBinError(bx, by) ** 2 + e ** 2))
    return out


def project_pt_bin_2d(h2, pt_bin, name):
    h = h2.ProjectionX(name, pt_bin, pt_bin, "e")
    h.SetDirectory(0)
    return h


def project_pt_bin_3d(h3, pt_bin, name):
    zaxis = h3.GetZaxis()
    old_first = zaxis.GetFirst()
    old_last = zaxis.GetLast()
    zaxis.SetRange(pt_bin, pt_bin)
    h = h3.Project3D("yx")
    h.SetName(name)
    h.SetDirectory(0)
    zaxis.SetRange(old_first, old_last)
    return h


def iterative_bayes(meas, resp_true_reco, prior_hist, n_iter, name):
    n_true = resp_true_reco.GetNbinsX()
    n_reco = resp_true_reco.GetNbinsY()
    prior = np.array([max(0.0, prior_hist.GetBinContent(i + 1)) for i in range(n_true)], dtype=float)
    s = prior.sum()
    if s <= 0.0:
        prior[:] = 1.0 / max(1, n_true)
    else:
        prior /= s
    P = np.array([[max(0.0, resp_true_reco.GetBinContent(t + 1, r + 1)) for r in range(n_reco)] for t in range(n_true)], dtype=float)
    unfolded = np.zeros(n_true, dtype=float)
    for _ in range(n_iter):
        unfolded[:] = 0.0
        for r in range(n_reco):
            mr = max(0.0, meas.GetBinContent(r + 1))
            if mr == 0.0:
                continue
            norm = np.dot(P[:, r], prior)
            if norm <= 0.0:
                continue
            unfolded += (P[:, r] * prior / norm) * mr
        s = np.clip(unfolded, 0.0, None).sum()
        if s <= 0.0:
            break
        prior = np.clip(unfolded, 0.0, None) / s
    h = ROOT.TH1D(name, name, n_true, prior_hist.GetXaxis().GetXmin(), prior_hist.GetXaxis().GetXmax())
    h.SetDirectory(0)
    h.Sumw2()
    for i in range(n_true):
        h.SetBinContent(i + 1, max(0.0, float(unfolded[i])))
        h.SetBinError(i + 1, math.sqrt(max(0.0, float(unfolded[i]))))
    return h


def build_ratio(num, den, name, title):
    h = num.Clone(name)
    h.SetDirectory(0)
    h.SetTitle(title)
    h.Divide(den)
    return h


def load_obj(f, name):
    obj = f.Get(name)
    if not obj:
        raise RuntimeError(f"Missing {name} in {f.GetName()}")
    out = obj.Clone(f"{name}_{os.path.basename(f.GetName()).replace('.', '_')}")
    out.SetDirectory(0)
    return out


def invert_matrix(m):
    try:
        cond = np.linalg.cond(m)
        if not np.isfinite(cond) or cond > 1e8:
            return np.linalg.pinv(m), False
        return np.linalg.inv(m), True
    except np.linalg.LinAlgError:
        return np.linalg.pinv(m), False


def summarize_shift(nominal, cross):
    max_abs = 0.0
    max_rel = 0.0
    for ib in range(1, nominal.GetNbinsX() + 1):
        nom = nominal.GetBinContent(ib)
        alt = cross.GetBinContent(ib)
        if nom == 0.0 and alt == 0.0:
            continue
        max_abs = max(max_abs, abs(alt - nom))
        if nom != 0.0:
            max_rel = max(max_rel, abs(alt / nom - 1.0))
    return max_abs, max_rel


def style_frame(h):
    h.GetXaxis().SetTitleSize(0.052)
    h.GetYaxis().SetTitleSize(0.052)
    h.GetXaxis().SetLabelSize(0.044)
    h.GetYaxis().SetLabelSize(0.044)
    h.GetXaxis().SetTitleOffset(1.05)
    h.GetYaxis().SetTitleOffset(1.20)


def style_legend(leg):
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)


def main():
    if len(sys.argv) != 6:
        raise SystemExit("Usage: run_order_swapped_dndeta_crosscheck.py <mc_inputs.root> <data_inputs.root> <nominal_summary.root> <nominal_unfold.root> <out_dir>")

    mc_inputs, data_inputs, nominal_summary_path, nominal_unfold_path, out_dir = sys.argv[1:]
    os.makedirs(out_dir, exist_ok=True)

    f_mc = ROOT.TFile.Open(mc_inputs, "READ")
    f_data = ROOT.TFile.Open(data_inputs, "READ")
    f_nom = ROOT.TFile.Open(nominal_summary_path, "READ")
    f_nom_unfold = ROOT.TFile.Open(nominal_unfold_path, "READ")
    if any((not f) or f.IsZombie() for f in [f_mc, f_data, f_nom, f_nom_unfold]):
        raise RuntimeError("Failed to open one or more inputs")

    h_tag_mc = {
        "K": load_obj(f_mc, "hTagKFakeCorrRecoDNdEta"),
        "Pi": load_obj(f_mc, "hTagPiFakeCorrRecoDNdEta"),
        "P": load_obj(f_mc, "hTagPFakeCorrRecoDNdEta"),
    }
    h_tag_data = {
        "K": load_obj(f_data, "hTagKFakeCorrRecoDNdEta"),
        "Pi": load_obj(f_data, "hTagPiFakeCorrRecoDNdEta"),
        "P": load_obj(f_data, "hTagPFakeCorrRecoDNdEta"),
    }
    h_resp = {
        "K": load_obj(f_mc, "hRespTagKFakeCorrDNdEta"),
        "Pi": load_obj(f_mc, "hRespTagPiFakeCorrDNdEta"),
        "P": load_obj(f_mc, "hRespTagPFakeCorrDNdEta"),
    }
    h_true_gen = {
        "K": load_obj(f_mc, "hGenAllKDNdEtaPt"),
        "Pi": load_obj(f_mc, "hGenAllPiDNdEtaPt"),
        "P": load_obj(f_mc, "hGenAllPDNdEtaPt"),
    }
    h_matched_den = {
        "K": load_obj(f_mc, "hMatchedDenKDNdEtaPt"),
        "Pi": load_obj(f_mc, "hMatchedDenPiDNdEtaPt"),
        "P": load_obj(f_mc, "hMatchedDenPDNdEtaPt"),
    }
    h_num = {
        ("K", "K"): load_obj(f_mc, "hTrueKTagAsKDNdEtaPt"),
        ("K", "Pi"): load_obj(f_mc, "hTrueKTagAsPiDNdEtaPt"),
        ("K", "P"): load_obj(f_mc, "hTrueKTagAsPDNdEtaPt"),
        ("Pi", "K"): load_obj(f_mc, "hTruePiTagAsKDNdEtaPt"),
        ("Pi", "Pi"): load_obj(f_mc, "hTruePiTagAsPiDNdEtaPt"),
        ("Pi", "P"): load_obj(f_mc, "hTruePiTagAsPDNdEtaPt"),
        ("P", "K"): load_obj(f_mc, "hTruePTagAsKDNdEtaPt"),
        ("P", "Pi"): load_obj(f_mc, "hTruePTagAsPiDNdEtaPt"),
        ("P", "P"): load_obj(f_mc, "hTruePTagAsPDNdEtaPt"),
    }

    n_pt = h_true_gen["K"].GetNbinsY()
    n_true_fine = h_true_gen["K"].GetNbinsX()

    unfolded_mc = {k: clone_empty_like(collapse_tail_2d(h_true_gen["K"], KEEP_BINS, n_pt, f"hEmptyMc_{k}"), f"hUnfoldedTagMc_{k}") for k in ["K", "Pi", "P"]}
    unfolded_data = {k: clone_empty_like(collapse_tail_2d(h_true_gen["K"], KEEP_BINS, n_pt, f"hEmptyData_{k}"), f"hUnfoldedTagData_{k}") for k in ["K", "Pi", "P"]}

    for obs in ["K", "Pi", "P"]:
        for ipt in range(1, n_pt + 1):
            meas_mc = collapse_tail_1d(project_pt_bin_2d(h_tag_mc[obs], ipt, f"hTagMc_{obs}_pt{ipt}"), KEEP_BINS, f"hTagMc_{obs}_pt{ipt}_collapsed")
            meas_data = collapse_tail_1d(project_pt_bin_2d(h_tag_data[obs], ipt, f"hTagData_{obs}_pt{ipt}"), KEEP_BINS, f"hTagData_{obs}_pt{ipt}_collapsed")
            resp = collapse_tail_2d(project_pt_bin_3d(h_resp[obs], ipt, f"hResp_{obs}_pt{ipt}"), KEEP_BINS, KEEP_BINS, f"hResp_{obs}_pt{ipt}_collapsed")
            prior = resp.ProjectionX(f"hPrior_{obs}_pt{ipt}")
            prior.SetDirectory(0)
            u_mc = iterative_bayes(meas_mc, resp, prior, N_ITER_BAYES, f"hUnfoldMc_{obs}_pt{ipt}")
            u_data = iterative_bayes(meas_data, resp, prior, N_ITER_BAYES, f"hUnfoldData_{obs}_pt{ipt}")
            for ib in range(1, KEEP_BINS + 1):
                unfolded_mc[obs].SetBinContent(ib, ipt, u_mc.GetBinContent(ib))
                unfolded_mc[obs].SetBinError(ib, ipt, u_mc.GetBinError(ib))
                unfolded_data[obs].SetBinContent(ib, ipt, u_data.GetBinContent(ib))
                unfolded_data[obs].SetBinError(ib, ipt, u_data.GetBinError(ib))

    h_true_gen_collapsed = {k: collapse_tail_2d(h_true_gen[k], KEEP_BINS, n_pt, f"{h_true_gen[k].GetName()}_collapsed") for k in ["K", "Pi", "P"]}
    h_matched_den_collapsed = {k: collapse_tail_2d(h_matched_den[k], KEEP_BINS, n_pt, f"{h_matched_den[k].GetName()}_collapsed") for k in ["K", "Pi", "P"]}
    h_num_collapsed = {key: collapse_tail_2d(h_num[key], KEEP_BINS, n_pt, f"{h_num[key].GetName()}_collapsed") for key in h_num}

    corrected_mc = {k: clone_empty_like(h_true_gen_collapsed["K"], f"hOrderSwapCorrectedMc_{k}") for k in ["K", "Pi", "P"]}
    corrected_data = {k: clone_empty_like(h_true_gen_collapsed["K"], f"hOrderSwapCorrectedData_{k}") for k in ["K", "Pi", "P"]}
    singular_cells = 0
    pseudo_inverse_cells = 0

    species = ["K", "Pi", "P"]
    obs_order = ["K", "Pi", "P"]
    for ib in range(1, KEEP_BINS + 1):
        for ipt in range(1, n_pt + 1):
            m = np.zeros((3, 3), dtype=float)
            ge = np.zeros(3, dtype=float)
            for ti, true_species in enumerate(species):
                den_matched = h_matched_den_collapsed[true_species].GetBinContent(ib, ipt)
                den_gen = h_true_gen_collapsed[true_species].GetBinContent(ib, ipt)
                ge[ti] = (den_matched / den_gen) if den_gen > 0.0 else 0.0
                if den_matched > 0.0:
                    for oj, obs_species in enumerate(obs_order):
                        m[oj, ti] = h_num_collapsed[(true_species, obs_species)].GetBinContent(ib, ipt) / den_matched
            if np.allclose(m, 0.0):
                singular_cells += 1
                continue
            minv, exact = invert_matrix(m)
            if not exact:
                pseudo_inverse_cells += 1
            y_mc = np.array([unfolded_mc[o].GetBinContent(ib, ipt) for o in obs_order], dtype=float)
            y_data = np.array([unfolded_data[o].GetBinContent(ib, ipt) for o in obs_order], dtype=float)
            x_matched_mc = minv.dot(y_mc)
            x_matched_data = minv.dot(y_data)
            for ti, true_species in enumerate(species):
                val_mc = max(0.0, float(x_matched_mc[ti]))
                val_data = max(0.0, float(x_matched_data[ti]))
                eff = ge[ti]
                if eff > 1e-12:
                    val_mc /= eff
                    val_data /= eff
                else:
                    val_mc = 0.0
                    val_data = 0.0
                corrected_mc[true_species].SetBinContent(ib, ipt, val_mc)
                corrected_data[true_species].SetBinContent(ib, ipt, val_data)
                corrected_mc[true_species].SetBinError(ib, ipt, math.sqrt(max(0.0, val_mc)))
                corrected_data[true_species].SetBinError(ib, ipt, math.sqrt(max(0.0, val_data)))

    def integrate_pt(h2, name):
        h = h2.ProjectionX(name, 1, h2.GetNbinsY(), "e")
        h.SetDirectory(0)
        return h

    h_mc_ratio_true = build_ratio(integrate_pt(h_true_gen_collapsed["K"], "hKTrue1D"), integrate_pt(h_true_gen_collapsed["Pi"], "hPiTrue1D"), "hRatioMcTruthOrderSwap", "MC truth K/#pi;dN_{ch}/d#eta (|#eta|<0.5);K/#pi")
    h_mc_ratio_cross = build_ratio(integrate_pt(corrected_mc["K"], "hKMcCross1D"), integrate_pt(corrected_mc["Pi"], "hPiMcCross1D"), "hRatioMcOrderSwap", "MC reordered K/#pi;dN_{ch}/d#eta (|#eta|<0.5);K/#pi")
    h_data_ratio_cross = build_ratio(integrate_pt(corrected_data["K"], "hKDataCross1D"), integrate_pt(corrected_data["Pi"], "hPiDataCross1D"), "hRatioDataOrderSwap", "Data reordered K/#pi;dN_{ch}/d#eta (|#eta|<0.5);K/#pi")
    h_closure = build_ratio(h_mc_ratio_cross, h_mc_ratio_true, "hClosureOrderSwap", "Reordered closure;dN_{ch}/d#eta (|#eta|<0.5);Unfolded / MC truth")
    h_double_ratio = build_ratio(h_data_ratio_cross, h_mc_ratio_cross, "hDoubleRatioOrderSwapPreResidual", "Reordered double ratio;dN_{ch}/d#eta (|#eta|<0.5);(K/#pi)_{Data}/(K/#pi)_{MC}")
    h_double_ratio_resid = h_double_ratio.Clone("hDoubleRatioOrderSwap")
    h_double_ratio_resid.SetDirectory(0)
    for ib in range(1, h_double_ratio_resid.GetNbinsX() + 1):
        c = h_closure.GetBinContent(ib)
        if abs(c) > 1e-12:
            h_double_ratio_resid.SetBinContent(ib, h_double_ratio_resid.GetBinContent(ib) / c)
            h_double_ratio_resid.SetBinError(ib, h_double_ratio_resid.GetBinError(ib) / abs(c))
    h_double_ratio_resid.Scale(1.0 / SF_RATIO)

    h_nominal = load_obj(f_nom, "hDoubleRatioNominal_dNdEta")
    h_nominal_sys = load_obj(f_nom, "hSysTotal_dNdEta")
    h_compare_ratio = build_ratio(h_double_ratio_resid, h_nominal, "hOrderSwapOverNominal", ";dN_{ch}/d#eta (|#eta|<0.5);Order-swapped / nominal")
    max_abs, max_rel = summarize_shift(h_nominal, h_double_ratio_resid)

    fout = ROOT.TFile.Open(os.path.join(out_dir, "order_swapped_dndeta_crosscheck.root"), "RECREATE")
    ROOT.TNamed("crosscheck_definition", "Observed tagged yields are fake-corrected and unfolded in dN/deta first, then untangled with a truth-matched 3x3 matrix and corrected with a truth-matched generator efficiency.").Write()
    ROOT.TNamed("crosscheck_scope", "Approximate reordered dN/deta cross-check. Activity is unfolded before the 3x3 species inversion within the existing factorized pT-binned framework; this is not the full Yi Chen 3D procedure.").Write()
    ROOT.TNamed("keepBins", str(KEEP_BINS)).Write()
    ROOT.TNamed("nIterBayes", str(N_ITER_BAYES)).Write()
    ROOT.TNamed("singularCells", str(singular_cells)).Write()
    ROOT.TNamed("pseudoInverseCells", str(pseudo_inverse_cells)).Write()
    ROOT.TNamed("maxAbsShift", f"{max_abs:.6f}").Write()
    ROOT.TNamed("maxRelShift", f"{max_rel:.6f}").Write()
    for obj in [h_mc_ratio_true, h_mc_ratio_cross, h_data_ratio_cross, h_closure, h_double_ratio, h_double_ratio_resid, h_nominal, h_compare_ratio]:
        obj.Write()
    h_nominal_sys.Write("hNominalSysTotal_dNdEta")
    for k in ["K", "Pi", "P"]:
        unfolded_mc[k].Write()
        unfolded_data[k].Write()
        corrected_mc[k].Write()
        corrected_data[k].Write()
    fout.Close()

    c = ROOT.TCanvas("c_order_swap_compare", "c_order_swap_compare", 900, 900)
    c.Divide(1, 2)
    pad1 = c.cd(1)
    pad1.SetPad(0.0, 0.33, 1.0, 1.0)
    pad1.SetBottomMargin(0.02)
    frame = h_nominal.Clone("hFrameOrderSwap")
    frame.Reset()
    frame.SetTitle("Order-swapped dN_{ch}/d#eta cross-check;dN_{ch}/d#eta (|#eta|<0.5);(K/#pi)_{Data}/(K/#pi)_{MC}")
    frame.SetMinimum(0.7)
    frame.SetMaximum(1.35)
    style_frame(frame)
    frame.Draw()
    band = ROOT.TGraphErrors(h_nominal.GetNbinsX())
    for ib in range(1, h_nominal.GetNbinsX() + 1):
        x = h_nominal.GetXaxis().GetBinCenter(ib)
        ex = 0.5 * h_nominal.GetXaxis().GetBinWidth(ib)
        band.SetPoint(ib - 1, x, h_nominal.GetBinContent(ib))
        band.SetPointError(ib - 1, ex, h_nominal_sys.GetBinContent(ib))
    band.SetFillColorAlpha(ROOT.kAzure - 9, 0.35)
    band.SetLineColor(ROOT.kAzure - 9)
    band.Draw("2 same")
    h_nominal.SetMarkerStyle(20)
    h_nominal.SetLineColor(ROOT.kBlack)
    h_nominal.SetMarkerColor(ROOT.kBlack)
    h_nominal.Draw("E1 same")
    h_double_ratio_resid.SetMarkerStyle(24)
    h_double_ratio_resid.SetMarkerSize(1.1)
    h_double_ratio_resid.SetMarkerColor(ROOT.kRed + 1)
    h_double_ratio_resid.SetLineColor(ROOT.kRed + 1)
    h_double_ratio_resid.Draw("E1 same")
    leg = ROOT.TLegend(0.52, 0.72, 0.88, 0.88)
    style_legend(leg)
    leg.AddEntry(h_nominal, "Nominal result", "lep")
    leg.AddEntry(band, "Nominal syst.", "f")
    leg.AddEntry(h_double_ratio_resid, "Order-swapped cross-check", "lep")
    leg.Draw()

    pad2 = c.cd(2)
    pad2.SetPad(0.0, 0.0, 1.0, 0.33)
    pad2.SetTopMargin(0.03)
    pad2.SetBottomMargin(0.28)
    ratio_frame = h_compare_ratio.Clone("hFrameOrderSwapRatio")
    ratio_frame.Reset()
    ratio_frame.SetTitle(";dN_{ch}/d#eta (|#eta|<0.5);Cross-check / nominal")
    ratio_frame.SetMinimum(0.85)
    ratio_frame.SetMaximum(1.15)
    style_frame(ratio_frame)
    ratio_frame.GetYaxis().SetNdivisions(505)
    ratio_frame.Draw()
    line = ROOT.TLine(ratio_frame.GetXaxis().GetXmin(), 1.0, ratio_frame.GetXaxis().GetXmax(), 1.0)
    line.SetLineStyle(2)
    line.Draw("same")
    h_compare_ratio.SetMarkerStyle(24)
    h_compare_ratio.SetMarkerColor(ROOT.kRed + 1)
    h_compare_ratio.SetLineColor(ROOT.kRed + 1)
    h_compare_ratio.Draw("E1 same")
    c.SaveAs(os.path.join(out_dir, "OrderSwapped_DNdEta_Comparison.pdf"))
    c.SaveAs(os.path.join(out_dir, "OrderSwapped_DNdEta_Comparison.png"))

    with open(os.path.join(out_dir, "order_swapped_dndeta_crosscheck.txt"), "w", encoding="ascii") as out:
        out.write("definition = fake-correct observed tags -> unfold activity -> truth-matched 3x3 inversion -> truth-matched gen-efficiency correction\n")
        out.write("scope = approximate reordered dN/deta cross-check within existing factorized framework; not full Yi 3D procedure\n")
        out.write(f"keepBins = {KEEP_BINS}\n")
        out.write(f"nIterBayes = {N_ITER_BAYES}\n")
        out.write(f"singularCells = {singular_cells}\n")
        out.write(f"pseudoInverseCells = {pseudo_inverse_cells}\n")
        out.write(f"maxAbsShift = {max_abs:.6f}\n")
        out.write(f"maxRelShift = {max_rel:.6f}\n")
        out.write("bin center nominal order_swapped ratio closure\n")
        for ib in range(1, h_nominal.GetNbinsX() + 1):
            out.write(
                f"{ib} {h_nominal.GetXaxis().GetBinCenter(ib):.3f} {h_nominal.GetBinContent(ib):.6f} {h_double_ratio_resid.GetBinContent(ib):.6f} {h_compare_ratio.GetBinContent(ib):.6f} {h_closure.GetBinContent(ib):.6f}\n"
            )

    print(f"Wrote {os.path.join(out_dir, 'order_swapped_dndeta_crosscheck.root')}")
    print(f"Wrote {os.path.join(out_dir, 'OrderSwapped_DNdEta_Comparison.pdf')}")
    print(f"Max absolute shift: {max_abs:.6f}")
    print(f"Max relative shift: {max_rel:.6f}")


if __name__ == "__main__":
    main()
