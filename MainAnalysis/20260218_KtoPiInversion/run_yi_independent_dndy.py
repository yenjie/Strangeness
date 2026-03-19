#!/usr/bin/env python3
import math
import os
import sys

import numpy as np
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetEndErrorSize(0)

SF_K = 0.9101
SF_PI = 0.9571
SF_RATIO = SF_K / SF_PI
N_ITER_BAYES = 1
P_GROUPS = ((0, 1), (2, 3), (4, 5), (6, 7))
COS_GROUPS = ((0, 1), (2, 3))


def load_obj(root_file, name):
    obj = root_file.Get(name)
    if not obj:
        raise RuntimeError(f"Missing {name} in {root_file.GetName()}")
    out = obj.Clone(f"{name}_{os.path.basename(root_file.GetName()).replace('.', '_')}")
    out.SetDirectory(0)
    return out


def extract_edges_1d(h):
    return [h.GetXaxis().GetBinLowEdge(i) for i in range(1, h.GetNbinsX() + 1)] + [h.GetXaxis().GetBinUpEdge(h.GetNbinsX())]


def collapse_tail_to_template(src, template, name):
    edges = np.array(extract_edges_1d(template), dtype="d")
    out = ROOT.TH1D(name, src.GetTitle(), template.GetNbinsX(), edges)
    out.SetDirectory(0)
    out.Sumw2()
    for ib in range(1, src.GetNbinsX() + 1):
        x_low = src.GetXaxis().GetBinLowEdge(ib)
        dest = out.FindBin(x_low + 1e-6)
        if dest < 1:
            dest = 1
        if dest > out.GetNbinsX():
            dest = out.GetNbinsX()
        c = src.GetBinContent(ib)
        e = src.GetBinError(ib)
        out.SetBinContent(dest, out.GetBinContent(dest) + c)
        out.SetBinError(dest, math.sqrt(out.GetBinError(dest) ** 2 + e ** 2))
    return out


def build_ratio(num, den, name, title):
    out = num.Clone(name)
    out.SetDirectory(0)
    out.SetTitle(title)
    out.Divide(den)
    return out


def invert_matrix(m):
    try:
        cond = np.linalg.cond(m)
        if not np.isfinite(cond) or cond > 1e8:
            return np.linalg.pinv(m), False
        return np.linalg.inv(m), True
    except np.linalg.LinAlgError:
        return np.linalg.pinv(m), False


def iterative_bayes(meas, resp_true_reco, prior_hist, n_iter, name):
    n_true = resp_true_reco.GetNbinsX()
    n_reco = resp_true_reco.GetNbinsY()
    prior = np.array([max(0.0, prior_hist.GetBinContent(i + 1)) for i in range(n_true)], dtype=float)
    s = prior.sum()
    if s <= 0.0:
        prior[:] = 1.0 / max(1, n_true)
    else:
        prior /= s

    p = np.array(
        [[max(0.0, resp_true_reco.GetBinContent(t + 1, r + 1)) for r in range(n_reco)] for t in range(n_true)],
        dtype=float,
    )

    unfolded = np.zeros(n_true, dtype=float)
    for _ in range(n_iter):
        unfolded[:] = 0.0
        for r in range(n_reco):
            mr = max(0.0, meas.GetBinContent(r + 1))
            if mr == 0.0:
                continue
            norm = np.dot(p[:, r], prior)
            if norm <= 0.0:
                continue
            unfolded += (p[:, r] * prior / norm) * mr
        s = np.clip(unfolded, 0.0, None).sum()
        if s <= 0.0:
            break
        prior = np.clip(unfolded, 0.0, None) / s

    out = ROOT.TH1D(name, name, n_true, prior_hist.GetXaxis().GetXmin(), prior_hist.GetXaxis().GetXmax())
    out.SetDirectory(0)
    out.Sumw2()
    for i in range(n_true):
        value = max(0.0, float(unfolded[i]))
        out.SetBinContent(i + 1, value)
        out.SetBinError(i + 1, math.sqrt(value))
    return out


def rms_distance_from_unity(h):
    vals = []
    for ib in range(1, h.GetNbinsX() + 1):
        y = h.GetBinContent(ib)
        ey = h.GetBinError(ib)
        if y == 0.0 and ey == 0.0:
            continue
        vals.append((y - 1.0) ** 2)
    return math.sqrt(sum(vals) / len(vals)) if vals else 0.0


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


def build_coarse_mapping(n_act, n_p_fine, n_cos_fine, p_groups, cos_groups):
    p_map = {p: i for i, group in enumerate(p_groups) for p in group}
    cos_map = {c: i for i, group in enumerate(cos_groups) for c in group}
    n_pc_fine = n_p_fine * n_cos_fine
    n_pc_coarse = len(p_groups) * len(cos_groups)
    fine_to_coarse = []
    for iflat in range(n_act * n_pc_fine):
        act = iflat // n_pc_fine
        rem = iflat % n_pc_fine
        p_bin = rem // n_cos_fine
        cos_bin = rem % n_cos_fine
        coarse = (act * n_pc_coarse) + (p_map[p_bin] * len(cos_groups)) + cos_map[cos_bin]
        fine_to_coarse.append(coarse)
    return fine_to_coarse, len(p_groups), len(cos_groups), n_act * n_pc_coarse


def aggregate_flat_1d(src, fine_to_coarse, n_coarse, name):
    out = ROOT.TH1D(name, src.GetTitle(), n_coarse, -0.5, n_coarse - 0.5)
    out.SetDirectory(0)
    out.Sumw2()
    for iflat, coarse in enumerate(fine_to_coarse):
        c = src.GetBinContent(iflat + 1)
        e = src.GetBinError(iflat + 1)
        out.SetBinContent(coarse + 1, out.GetBinContent(coarse + 1) + c)
        out.SetBinError(coarse + 1, math.sqrt(out.GetBinError(coarse + 1) ** 2 + e ** 2))
    return out


def aggregate_flat_2d(src, fine_to_coarse, n_coarse, name):
    out = ROOT.TH2D(name, src.GetTitle(), n_coarse, -0.5, n_coarse - 0.5, n_coarse, -0.5, n_coarse - 0.5)
    out.SetDirectory(0)
    out.Sumw2()
    n_fine = len(fine_to_coarse)
    for ix in range(n_fine):
        coarse_x = fine_to_coarse[ix]
        for iy in range(n_fine):
            c = src.GetBinContent(ix + 1, iy + 1)
            if c == 0.0:
                continue
            coarse_y = fine_to_coarse[iy]
            e = src.GetBinError(ix + 1, iy + 1)
            out.SetBinContent(coarse_x + 1, coarse_y + 1, out.GetBinContent(coarse_x + 1, coarse_y + 1) + c)
            out.SetBinError(
                coarse_x + 1,
                coarse_y + 1,
                math.sqrt(out.GetBinError(coarse_x + 1, coarse_y + 1) ** 2 + e ** 2),
            )
    return out


def make_activity_projection(flat_hist, activity_template, n_act, n_p, n_cos, name, title):
    edges = np.array(extract_edges_1d(activity_template), dtype="d")
    out = ROOT.TH1D(name, title, n_act, edges)
    out.SetDirectory(0)
    out.Sumw2()
    block = n_p * n_cos
    for iflat in range(n_act * block):
        act = iflat // block
        c = flat_hist.GetBinContent(iflat + 1)
        e = flat_hist.GetBinError(iflat + 1)
        out.SetBinContent(act + 1, out.GetBinContent(act + 1) + c)
        out.SetBinError(act + 1, math.sqrt(out.GetBinError(act + 1) ** 2 + e ** 2))
    return out


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


def draw_comparison(h_nominal, h_nom_sys, h_alt, out_path, title):
    c = ROOT.TCanvas("c_yi_compare", "", 900, 900)
    c.Divide(1, 2)
    pad1 = c.cd(1)
    pad1.SetPad(0.0, 0.33, 1.0, 1.0)
    pad1.SetBottomMargin(0.02)

    frame = h_nominal.Clone("hFrameYiCompare")
    frame.Reset()
    frame.SetTitle(title)
    frame.SetMinimum(0.70)
    frame.SetMaximum(1.35)
    style_frame(frame)
    frame.Draw()

    band = ROOT.TGraphErrors(h_nominal.GetNbinsX())
    for ib in range(1, h_nominal.GetNbinsX() + 1):
        x = h_nominal.GetXaxis().GetBinCenter(ib)
        ex = 0.5 * h_nominal.GetXaxis().GetBinWidth(ib)
        band.SetPoint(ib - 1, x, h_nominal.GetBinContent(ib))
        band.SetPointError(ib - 1, ex, h_nom_sys.GetBinContent(ib))
    band.SetFillColorAlpha(ROOT.kAzure - 9, 0.35)
    band.SetLineColor(ROOT.kAzure - 9)
    band.Draw("2 same")

    h_nominal.SetLineColor(ROOT.kBlack)
    h_nominal.SetMarkerColor(ROOT.kBlack)
    h_nominal.SetMarkerStyle(20)
    h_nominal.Draw("E1 same")

    h_alt.SetLineColor(ROOT.kRed + 1)
    h_alt.SetMarkerColor(ROOT.kRed + 1)
    h_alt.SetMarkerStyle(24)
    h_alt.SetMarkerSize(1.1)
    h_alt.Draw("E1 same")

    leg = ROOT.TLegend(0.50, 0.72, 0.88, 0.88)
    style_legend(leg)
    leg.AddEntry(h_nominal, "Nominal result", "lep")
    leg.AddEntry(band, "Nominal syst.", "f")
    leg.AddEntry(h_alt, "Yi-style independent", "lep")
    leg.Draw()

    pad2 = c.cd(2)
    pad2.SetPad(0.0, 0.0, 1.0, 0.33)
    pad2.SetTopMargin(0.03)
    pad2.SetBottomMargin(0.28)
    ratio = build_ratio(h_alt, h_nominal, "hYiOverNominal", ";dN_{ch}/dy (thrust axis, |y_{T}|<0.5);Yi / nominal")
    rframe = ratio.Clone("hFrameYiRatio")
    rframe.Reset()
    rframe.SetMinimum(0.85)
    rframe.SetMaximum(1.15)
    style_frame(rframe)
    rframe.GetYaxis().SetNdivisions(505)
    rframe.Draw()
    line = ROOT.TLine(rframe.GetXaxis().GetXmin(), 1.0, rframe.GetXaxis().GetXmax(), 1.0)
    line.SetLineStyle(2)
    line.Draw("same")
    ratio.SetLineColor(ROOT.kRed + 1)
    ratio.SetMarkerColor(ROOT.kRed + 1)
    ratio.SetMarkerStyle(24)
    ratio.Draw("E1 same")

    c.SaveAs(out_path)
    c.SaveAs(out_path.replace(".pdf", ".png"))
    return ratio


def draw_closure_comparison(h_nominal_closure, h_alt_closure, out_path):
    c = ROOT.TCanvas("c_yi_closure", "", 900, 900)
    c.Divide(1, 2)
    pad1 = c.cd(1)
    pad1.SetPad(0.0, 0.33, 1.0, 1.0)
    pad1.SetBottomMargin(0.02)

    frame = h_nominal_closure.Clone("hFrameYiClosure")
    frame.Reset()
    frame.SetTitle(";dN_{ch}/dy (thrust axis, |y_{T}|<0.5);Unfolded / MC truth")
    frame.SetMinimum(0.90)
    frame.SetMaximum(1.15)
    style_frame(frame)
    frame.Draw()
    line = ROOT.TLine(frame.GetXaxis().GetXmin(), 1.0, frame.GetXaxis().GetXmax(), 1.0)
    line.SetLineStyle(2)
    line.SetLineWidth(2)
    line.Draw("same")

    h_nominal_closure.SetLineColor(ROOT.kBlack)
    h_nominal_closure.SetMarkerColor(ROOT.kBlack)
    h_nominal_closure.SetMarkerStyle(20)
    h_nominal_closure.SetLineWidth(2)
    h_nominal_closure.Draw("E1 same")

    h_alt_closure.SetLineColor(ROOT.kRed + 1)
    h_alt_closure.SetMarkerColor(ROOT.kRed + 1)
    h_alt_closure.SetMarkerStyle(24)
    h_alt_closure.SetLineWidth(2)
    h_alt_closure.Draw("E1 same")

    leg = ROOT.TLegend(0.18, 0.74, 0.58, 0.88)
    style_legend(leg)
    leg.AddEntry(h_nominal_closure, "Nominal closure", "lep")
    leg.AddEntry(h_alt_closure, "Yi-style closure", "lep")
    leg.Draw()

    lab = ROOT.TLatex()
    lab.SetNDC()
    lab.SetTextSize(0.035)
    lab.DrawLatex(0.60, 0.84, f"Nominal RMS = {rms_distance_from_unity(h_nominal_closure):.3f}")
    lab.DrawLatex(0.60, 0.79, f"Yi-style RMS = {rms_distance_from_unity(h_alt_closure):.3f}")

    pad2 = c.cd(2)
    pad2.SetPad(0.0, 0.0, 1.0, 0.33)
    pad2.SetTopMargin(0.03)
    pad2.SetBottomMargin(0.28)
    ratio = build_ratio(
        h_alt_closure,
        h_nominal_closure,
        "hYiClosureOverNominal",
        ";dN_{ch}/dy (thrust axis, |y_{T}|<0.5);Yi / nominal",
    )
    rframe = ratio.Clone("hFrameYiClosureRatio")
    rframe.Reset()
    rframe.SetMinimum(0.90)
    rframe.SetMaximum(1.10)
    style_frame(rframe)
    rframe.GetYaxis().SetNdivisions(505)
    rframe.Draw()
    line2 = ROOT.TLine(rframe.GetXaxis().GetXmin(), 1.0, rframe.GetXaxis().GetXmax(), 1.0)
    line2.SetLineStyle(2)
    line2.Draw("same")
    ratio.SetLineColor(ROOT.kRed + 1)
    ratio.SetMarkerColor(ROOT.kRed + 1)
    ratio.SetMarkerStyle(24)
    ratio.Draw("E1 same")
    c.SaveAs(out_path)
    c.SaveAs(out_path.replace(".pdf", ".png"))
    return ratio


def main():
    if len(sys.argv) != 6:
        raise SystemExit(
            "Usage: run_yi_independent_dndy.py "
            "<mc_inputs.root> <data_inputs.root> <nominal_summary.root> <nominal_unfold.root> <out_dir>"
        )

    mc_inputs, data_inputs, nominal_summary_path, nominal_unfold_path, out_dir = sys.argv[1:]
    os.makedirs(out_dir, exist_ok=True)

    f_mc = ROOT.TFile.Open(mc_inputs, "READ")
    f_data = ROOT.TFile.Open(data_inputs, "READ")
    f_nom = ROOT.TFile.Open(nominal_summary_path, "READ")
    f_nom_unfold = ROOT.TFile.Open(nominal_unfold_path, "READ")
    if any((not f) or f.IsZombie() for f in [f_mc, f_data, f_nom, f_nom_unfold]):
        raise RuntimeError("Failed to open one or more inputs")

    h_tag_mc = {
        "K": load_obj(f_mc, "hRawTagKRecoFlat"),
        "Pi": load_obj(f_mc, "hRawTagPiRecoFlat"),
        "P": load_obj(f_mc, "hRawTagPRecoFlat"),
    }
    h_tag_data = {
        "K": load_obj(f_data, "hRawTagKRecoFlat"),
        "Pi": load_obj(f_data, "hRawTagPiRecoFlat"),
        "P": load_obj(f_data, "hRawTagPRecoFlat"),
    }
    h_resp = {
        "K": load_obj(f_mc, "hRespTagKFlat"),
        "Pi": load_obj(f_mc, "hRespTagPiFlat"),
        "P": load_obj(f_mc, "hRespTagPFlat"),
    }
    h_true_gen = {
        "K": load_obj(f_mc, "hGenAllKFlat"),
        "Pi": load_obj(f_mc, "hGenAllPiFlat"),
        "P": load_obj(f_mc, "hGenAllPFlat"),
    }
    h_matched_den = {
        "K": load_obj(f_mc, "hMatchedDenKFlat"),
        "Pi": load_obj(f_mc, "hMatchedDenPiFlat"),
        "P": load_obj(f_mc, "hMatchedDenPFlat"),
    }
    h_num = {
        ("K", "K"): load_obj(f_mc, "hTrueKTagAsKFlat"),
        ("K", "Pi"): load_obj(f_mc, "hTrueKTagAsPiFlat"),
        ("K", "P"): load_obj(f_mc, "hTrueKTagAsPFlat"),
        ("Pi", "K"): load_obj(f_mc, "hTruePiTagAsKFlat"),
        ("Pi", "Pi"): load_obj(f_mc, "hTruePiTagAsPiFlat"),
        ("Pi", "P"): load_obj(f_mc, "hTruePiTagAsPFlat"),
        ("P", "K"): load_obj(f_mc, "hTruePTagAsKFlat"),
        ("P", "Pi"): load_obj(f_mc, "hTruePTagAsPiFlat"),
        ("P", "P"): load_obj(f_mc, "hTruePTagAsPFlat"),
    }
    h_reco_counts = load_obj(f_mc, "hRecoCountsDNdY")
    p_edges = load_obj(f_mc, "hPEdges")
    cos_edges = load_obj(f_mc, "hAbsCosEdges")

    n_act = h_reco_counts.GetNbinsX()
    n_p_fine = p_edges.GetNbinsX()
    n_cos_fine = cos_edges.GetNbinsX()
    n_flat_fine = h_true_gen["K"].GetNbinsX()
    if n_flat_fine != n_act * n_p_fine * n_cos_fine:
        raise RuntimeError(
            f"Inconsistent flat size: nFlat={n_flat_fine}, nAct={n_act}, nP={n_p_fine}, nCos={n_cos_fine}"
        )
    fine_to_coarse, n_p, n_cos, n_flat = build_coarse_mapping(
        n_act, n_p_fine, n_cos_fine, P_GROUPS, COS_GROUPS
    )

    h_tag_mc = {
        k: aggregate_flat_1d(v, fine_to_coarse, n_flat, f"{v.GetName()}_coarse") for k, v in h_tag_mc.items()
    }
    h_tag_data = {
        k: aggregate_flat_1d(v, fine_to_coarse, n_flat, f"{v.GetName()}_coarse") for k, v in h_tag_data.items()
    }
    h_resp = {
        k: aggregate_flat_2d(v, fine_to_coarse, n_flat, f"{v.GetName()}_coarse") for k, v in h_resp.items()
    }
    h_true_gen = {
        k: aggregate_flat_1d(v, fine_to_coarse, n_flat, f"{v.GetName()}_coarse") for k, v in h_true_gen.items()
    }
    h_matched_den = {
        k: aggregate_flat_1d(v, fine_to_coarse, n_flat, f"{v.GetName()}_coarse") for k, v in h_matched_den.items()
    }
    h_num = {
        key: aggregate_flat_1d(v, fine_to_coarse, n_flat, f"{v.GetName()}_coarse") for key, v in h_num.items()
    }

    unfolded_mc = {}
    unfolded_data = {}
    for obs in ["K", "Pi", "P"]:
        prior = h_resp[obs].ProjectionX(f"hPrior_{obs}")
        prior.SetDirectory(0)
        unfolded_mc[obs] = iterative_bayes(
            h_tag_mc[obs], h_resp[obs], prior, N_ITER_BAYES, f"hUnfoldedTagMc_{obs}"
        )
        unfolded_data[obs] = iterative_bayes(
            h_tag_data[obs], h_resp[obs], prior, N_ITER_BAYES, f"hUnfoldedTagData_{obs}"
        )

    corrected_mc = {k: h_true_gen["K"].Clone(f"hYiCorrectedMc_{k}") for k in ["K", "Pi", "P"]}
    corrected_data = {k: h_true_gen["K"].Clone(f"hYiCorrectedData_{k}") for k in ["K", "Pi", "P"]}
    for h in list(corrected_mc.values()) + list(corrected_data.values()):
        h.SetDirectory(0)
        h.Reset()
        h.Sumw2()

    species = ["K", "Pi", "P"]
    obs_order = ["K", "Pi", "P"]
    singular_cells = 0
    pseudo_inverse_cells = 0
    zero_gen_eff_cells = 0

    for ib in range(1, n_flat + 1):
        m = np.zeros((3, 3), dtype=float)
        ge = np.zeros(3, dtype=float)
        for ti, true_species in enumerate(species):
            den_matched = h_matched_den[true_species].GetBinContent(ib)
            den_gen = h_true_gen[true_species].GetBinContent(ib)
            ge[ti] = (den_matched / den_gen) if den_gen > 0.0 else 0.0
            if den_gen <= 0.0:
                zero_gen_eff_cells += 1
            if den_matched > 0.0:
                for oj, obs_species in enumerate(obs_order):
                    m[oj, ti] = h_num[(true_species, obs_species)].GetBinContent(ib) / den_matched
        if np.allclose(m, 0.0):
            singular_cells += 1
            continue

        minv, exact = invert_matrix(m)
        if not exact:
            pseudo_inverse_cells += 1

        y_mc = np.array([unfolded_mc[o].GetBinContent(ib) for o in obs_order], dtype=float)
        y_data = np.array([unfolded_data[o].GetBinContent(ib) for o in obs_order], dtype=float)
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
            corrected_mc[true_species].SetBinContent(ib, val_mc)
            corrected_mc[true_species].SetBinError(ib, math.sqrt(val_mc))
            corrected_data[true_species].SetBinContent(ib, val_data)
            corrected_data[true_species].SetBinError(ib, math.sqrt(val_data))

    h_k_true_fine = make_activity_projection(
        h_true_gen["K"],
        h_reco_counts,
        n_act,
        n_p,
        n_cos,
        "hYiKTruthFine",
        "MC truth K;dN_{ch}/dy (thrust axis, |y_{T}|<0.5);Counts",
    )
    h_pi_true_fine = make_activity_projection(
        h_true_gen["Pi"],
        h_reco_counts,
        n_act,
        n_p,
        n_cos,
        "hYiPiTruthFine",
        "MC truth #pi;dN_{ch}/dy (thrust axis, |y_{T}|<0.5);Counts",
    )
    h_k_mc_fine = make_activity_projection(
        corrected_mc["K"],
        h_reco_counts,
        n_act,
        n_p,
        n_cos,
        "hYiKMcFine",
        "Yi-style MC K;dN_{ch}/dy (thrust axis, |y_{T}|<0.5);Counts",
    )
    h_pi_mc_fine = make_activity_projection(
        corrected_mc["Pi"],
        h_reco_counts,
        n_act,
        n_p,
        n_cos,
        "hYiPiMcFine",
        "Yi-style MC #pi;dN_{ch}/dy (thrust axis, |y_{T}|<0.5);Counts",
    )
    h_k_data_fine = make_activity_projection(
        corrected_data["K"],
        h_reco_counts,
        n_act,
        n_p,
        n_cos,
        "hYiKDataFine",
        "Yi-style data K;dN_{ch}/dy (thrust axis, |y_{T}|<0.5);Counts",
    )
    h_pi_data_fine = make_activity_projection(
        corrected_data["Pi"],
        h_reco_counts,
        n_act,
        n_p,
        n_cos,
        "hYiPiDataFine",
        "Yi-style data #pi;dN_{ch}/dy (thrust axis, |y_{T}|<0.5);Counts",
    )

    h_nominal = load_obj(f_nom, "hDoubleRatioNominal_dNdY")
    h_nominal_sys = load_obj(f_nom, "hSysTotal_dNdY")
    h_nominal_closure = load_obj(f_nom_unfold, "hClosureBayes_dNdY")
    h_k_true = collapse_tail_to_template(h_k_true_fine, h_nominal, "hYiKTruth")
    h_pi_true = collapse_tail_to_template(h_pi_true_fine, h_nominal, "hYiPiTruth")
    h_k_mc = collapse_tail_to_template(h_k_mc_fine, h_nominal, "hYiKMc")
    h_pi_mc = collapse_tail_to_template(h_pi_mc_fine, h_nominal, "hYiPiMc")
    h_k_data = collapse_tail_to_template(h_k_data_fine, h_nominal, "hYiKData")
    h_pi_data = collapse_tail_to_template(h_pi_data_fine, h_nominal, "hYiPiData")

    h_ratio_mc_truth = build_ratio(h_k_true, h_pi_true, "hYiRatioMcTruth", "MC truth K/#pi")
    h_ratio_mc = build_ratio(h_k_mc, h_pi_mc, "hYiRatioMc", "Yi-style MC K/#pi")
    h_ratio_data = build_ratio(h_k_data, h_pi_data, "hYiRatioData", "Yi-style data K/#pi")
    h_closure = build_ratio(h_ratio_mc, h_ratio_mc_truth, "hYiClosure", "Yi-style closure")

    h_double_ratio = build_ratio(
        h_ratio_data,
        h_ratio_mc,
        "hYiDoubleRatioPreResidual",
        "Yi-style double ratio;dN_{ch}/dy (|#eta|<0.5);(K/#pi)_{Data}/(K/#pi)_{MC}",
    )
    h_double_ratio_resid = h_double_ratio.Clone("hYiDoubleRatio")
    h_double_ratio_resid.SetDirectory(0)
    for ib in range(1, h_double_ratio_resid.GetNbinsX() + 1):
        c = h_closure.GetBinContent(ib)
        if abs(c) > 1e-12:
            h_double_ratio_resid.SetBinContent(ib, h_double_ratio_resid.GetBinContent(ib) / c)
            h_double_ratio_resid.SetBinError(ib, h_double_ratio_resid.GetBinError(ib) / abs(c))
    h_double_ratio_resid.Scale(1.0 / SF_RATIO)

    h_compare_ratio = draw_comparison(
        h_nominal,
        h_nominal_sys,
        h_double_ratio_resid,
        os.path.join(out_dir, "YiIndependent_DNdY_Comparison.pdf"),
        "Yi-style independent dN_{ch}/dy cross-check;dN_{ch}/dy (thrust axis, |y_{T}|<0.5);(K/#pi)_{Data}/(K/#pi)_{MC}",
    )
    h_closure_compare_ratio = draw_closure_comparison(
        h_nominal_closure,
        h_closure,
        os.path.join(out_dir, "YiIndependent_DNdY_ClosureComparison.pdf"),
    )

    max_abs, max_rel = summarize_shift(h_nominal, h_double_ratio_resid)
    max_shift_over_sys = 0.0
    max_shift_over_sys_bin = 0
    h_shift = h_nominal.Clone("hYiMethodShift")
    h_shift.SetDirectory(0)
    h_shift.Reset()
    h_shift_over_sys = h_nominal.Clone("hYiShiftOverNominalSys")
    h_shift_over_sys.SetDirectory(0)
    h_shift_over_sys.Reset()
    for ib in range(1, h_nominal.GetNbinsX() + 1):
        shift = abs(h_double_ratio_resid.GetBinContent(ib) - h_nominal.GetBinContent(ib))
        nominal_sys = h_nominal_sys.GetBinContent(ib)
        h_shift.SetBinContent(ib, shift)
        ratio = shift / nominal_sys if nominal_sys > 0.0 else 0.0
        h_shift_over_sys.SetBinContent(ib, ratio)
        if ratio > max_shift_over_sys:
            max_shift_over_sys = ratio
            max_shift_over_sys_bin = ib

    closure_rms_nom = rms_distance_from_unity(h_nominal_closure)
    closure_rms_yi = rms_distance_from_unity(h_closure)

    fout = ROOT.TFile.Open(os.path.join(out_dir, "yi_independent_dndy.root"), "RECREATE")
    ROOT.TNamed(
        "definition",
        "Independent Yi-style dN/dy chain: fake-correct observed tags in full flat mu=(dN/dy,p,|cos(theta)|), unfold observed tags in flat mu, apply truth-cell 3x3 inversion, then divide by truth-cell matched/gen efficiency.",
    ).Write()
    ROOT.TNamed("nIterBayes", str(N_ITER_BAYES)).Write()
    ROOT.TNamed("nActivityBinsFine", str(n_act)).Write()
    ROOT.TNamed("nPBins", str(n_p)).Write()
    ROOT.TNamed("nAbsCosBins", str(n_cos)).Write()
    ROOT.TNamed("pGroups", str(P_GROUPS)).Write()
    ROOT.TNamed("absCosGroups", str(COS_GROUPS)).Write()
    ROOT.TNamed("singularCells", str(singular_cells)).Write()
    ROOT.TNamed("pseudoInverseCells", str(pseudo_inverse_cells)).Write()
    ROOT.TNamed("zeroGenEffCells", str(zero_gen_eff_cells)).Write()
    ROOT.TNamed("maxAbsShift", f"{max_abs:.6f}").Write()
    ROOT.TNamed("maxRelShift", f"{max_rel:.6f}").Write()
    for obj in [
        h_ratio_mc_truth,
        h_ratio_mc,
        h_ratio_data,
        h_closure,
        h_double_ratio,
        h_double_ratio_resid,
        h_nominal,
        h_compare_ratio,
        h_nominal_closure,
        h_closure_compare_ratio,
        h_nominal_sys,
        h_shift,
        h_shift_over_sys,
        h_k_true_fine,
        h_pi_true_fine,
        h_k_mc_fine,
        h_pi_mc_fine,
        h_k_data_fine,
        h_pi_data_fine,
        h_k_true,
        h_pi_true,
        h_k_mc,
        h_pi_mc,
        h_k_data,
        h_pi_data,
    ]:
        obj.Write()
    for obs in ["K", "Pi", "P"]:
        unfolded_mc[obs].Write()
        unfolded_data[obs].Write()
    for sp in ["K", "Pi", "P"]:
        corrected_mc[sp].Write()
        corrected_data[sp].Write()
    fout.Close()

    with open(os.path.join(out_dir, "yi_independent_dndy.txt"), "w", encoding="ascii") as out:
        out.write("definition = fake-correct observed tags in flat mu -> unfold observed tags in flat mu -> truth-cell 3x3 inversion -> truth-cell gen-eff correction\n")
        out.write(f"nIterBayes = {N_ITER_BAYES}\n")
        out.write(f"nActivityBinsFine = {n_act}\n")
        out.write(f"nPBins = {n_p}\n")
        out.write(f"nAbsCosBins = {n_cos}\n")
        out.write(f"singularCells = {singular_cells}\n")
        out.write(f"pseudoInverseCells = {pseudo_inverse_cells}\n")
        out.write(f"zeroGenEffCells = {zero_gen_eff_cells}\n")
        out.write(f"nominalClosureRMS = {closure_rms_nom:.6f}\n")
        out.write(f"yiClosureRMS = {closure_rms_yi:.6f}\n")
        out.write(f"maxAbsShift = {max_abs:.6f}\n")
        out.write(f"maxRelShift = {max_rel:.6f}\n")
        out.write(f"maxShiftOverNominalSys = {max_shift_over_sys:.6f}\n")
        out.write(f"maxShiftOverNominalSysBin = {max_shift_over_sys_bin}\n")
        out.write("bin center nominal yi ratio closure shiftOverSys\n")
        for ib in range(1, h_nominal.GetNbinsX() + 1):
            out.write(
                f"{ib} {h_nominal.GetXaxis().GetBinCenter(ib):.3f} "
                f"{h_nominal.GetBinContent(ib):.6f} {h_double_ratio_resid.GetBinContent(ib):.6f} "
                f"{h_compare_ratio.GetBinContent(ib):.6f} {h_closure.GetBinContent(ib):.6f} "
                f"{h_shift_over_sys.GetBinContent(ib):.6f}\n"
            )

    print(f"Wrote {os.path.join(out_dir, 'yi_independent_dndy.root')}")
    print(f"Wrote {os.path.join(out_dir, 'YiIndependent_DNdY_Comparison.pdf')}")
    print(f"Wrote {os.path.join(out_dir, 'YiIndependent_DNdY_ClosureComparison.pdf')}")
    print(f"Yi closure RMS: {closure_rms_yi:.6f}")
    print(f"Max absolute shift: {max_abs:.6f}")
    print(f"Max relative shift: {max_rel:.6f}")


if __name__ == "__main__":
    main()
