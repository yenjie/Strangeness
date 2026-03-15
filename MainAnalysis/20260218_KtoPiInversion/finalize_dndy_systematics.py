#!/usr/bin/env python3
import math
import os
import ROOT
import csv

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
ACTIVE_KEEPBINS = int(os.environ.get("KEEPBINS_OVERRIDE", "-1"))
OUT_DIR = os.environ.get("KPI_DNDY_DIR", "output/systematics_20260314_dndy")
LOGO_PATH = os.environ.get("EEA_LOGO_PATH", "assets/eea_logo.png")


def style_frame(h):
    h.GetXaxis().SetTitleSize(0.052)
    h.GetYaxis().SetTitleSize(0.052)
    h.GetXaxis().SetLabelSize(0.044)
    h.GetYaxis().SetLabelSize(0.044)
    h.GetXaxis().SetTitleOffset(1.05)
    h.GetYaxis().SetTitleOffset(1.15)


def style_legend(leg):
    leg.SetBorderSize(0)
    leg.SetFillColor(ROOT.kWhite)
    leg.SetFillStyle(1001)
    leg.SetLineColor(ROOT.kWhite)


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


def load_toy_stat_errors(axis_name):
    out = {}
    if not os.path.exists(TOY_COVERAGE_CSV):
        return out
    with open(TOY_COVERAGE_CSV, newline="", encoding="ascii") as f:
        reader = csv.DictReader(f)
        for row in reader:
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

out_dir = OUT_DIR
variations = [
    ('nominal', 'nominal'),
    ('npt10', 'binning'),
    ('npt14', 'binning'),
]


def read_hist(path, name):
    f = ROOT.TFile.Open(path, 'READ')
    if not f or f.IsZombie():
        raise RuntimeError(f'Cannot open {path}')
    h = f.Get(name)
    if not h:
        raise RuntimeError(f'Missing {name} in {path}')
    hc = h.Clone(f'{name}_{os.path.basename(path).replace(".","_")}')
    hc.SetDirectory(0)
    f.Close()
    return hc


def read_keepbins_used(path):
    f = ROOT.TFile.Open(path, "READ")
    if not f or f.IsZombie():
        return ACTIVE_KEEPBINS
    obj = f.Get("keepBinsUsed_dNdY")
    value = ACTIVE_KEEPBINS
    if obj and hasattr(obj, "GetVal"):
        value = int(obj.GetVal())
    f.Close()
    return value


def apply_residual_correction(h, h_closure):
    hc = h.Clone(h.GetName() + "_residcorr")
    hc.SetDirectory(0)
    for ib in range(1, hc.GetNbinsX() + 1):
        cval = h_closure.GetBinContent(ib)
        if abs(cval) > 1e-12:
            hc.SetBinContent(ib, hc.GetBinContent(ib) / cval)
            hc.SetBinError(ib, hc.GetBinError(ib) / abs(cval))
    return hc

NOMINAL_PATH = os.path.join(out_dir, 'nominal_unfold_dndy.root')
KEEPBINS_USED = read_keepbins_used(NOMINAL_PATH)

h_double = {}
h_closure_by_var = {}
for v, _ in variations:
    p = os.path.join(out_dir, f'{v}_unfold_dndy.root')
    h_double[v] = read_hist(p, 'hDataOverMcBayes_dNdY')
    h_closure_by_var[v] = read_hist(p, 'hClosureBayes_dNdY')

h_closure = h_closure_by_var['nominal']
for v, _ in variations:
    h_double[v] = apply_residual_correction(h_double[v], h_closure_by_var[v])
    h_double[v].Scale(SF_RATIO)

nominal = h_double['nominal']
toy_override_used = apply_stat_override(nominal, load_toy_stat_errors("dNdY"))
nbin = nominal.GetNbinsX()
last_active = 0
for ib in range(1, nbin + 1):
    if nominal.GetBinContent(ib) != 0.0 or nominal.GetBinError(ib) != 0.0:
        last_active = ib

source = {
    'binning': [0.0] * nbin,
    'residual': [0.0] * nbin,
    'unfolding_prior': [0.0] * nbin,
    'unfolding_iter': [0.0] * nbin,
    'unfolding': [0.0] * nbin,
    'pid_sf': [0.0] * nbin,
}

groups = {'binning': []}
for name, g in variations:
    if g in groups:
        groups[g].append(name)

for g, names in groups.items():
    for ib in range(1, nbin + 1):
        nom = nominal.GetBinContent(ib)
        mx = 0.0
        for n in names:
            mx = max(mx, abs(h_double[n].GetBinContent(ib) - nom))
        source[g][ib - 1] = mx

h_prior = read_hist(os.path.join(out_dir, 'nominal_unfold_dndy.root'), 'hBayesPriorVariationDiff_dNdY')
h_iter = read_hist(os.path.join(out_dir, 'nominal_unfold_dndy.root'), 'hBayesIterVariationDiff_dNdY')
for ib in range(1, nbin + 1):
    if ib > last_active:
        continue
    cval = h_closure.GetBinContent(ib)
    inv_scale = (1.0 / abs(cval)) if abs(cval) > 1e-12 else 0.0
    rel_residual = 0.5 * abs((1.0 / cval) - 1.0) if abs(cval) > 1e-12 else 0.0
    source['residual'][ib - 1] = rel_residual * abs(nominal.GetBinContent(ib))
    source['unfolding_prior'][ib - 1] = abs(h_prior.GetBinContent(ib)) * inv_scale * SF_RATIO
    source['unfolding_iter'][ib - 1] = abs(h_iter.GetBinContent(ib)) * inv_scale * SF_RATIO
    source['unfolding'][ib - 1] = max(source['unfolding_prior'][ib - 1], source['unfolding_iter'][ib - 1])
    source['pid_sf'][ib - 1] = abs(nominal.GetBinContent(ib)) * SF_RATIO_REL_ERR

source['total_sys'] = [0.0] * nbin
for i in range(nbin):
    source['total_sys'][i] = math.sqrt(
        source['binning'][i] ** 2
        + source['residual'][i] ** 2
        + source['unfolding'][i] ** 2
        + source['pid_sf'][i] ** 2
    )

os.makedirs(out_dir, exist_ok=True)

# table
with open(os.path.join(out_dir, 'systematics_dndy_table.txt'), 'w', encoding='ascii') as f:
    f.write(f'# dN/dy status branch: keepBinsOverride={ACTIVE_KEEPBINS}, keepBinsUsed={KEEPBINS_USED}\n')
    f.write('# authoritative builder: finalize_dndy_systematics.py\n')
    if toy_override_used:
        f.write('# stat treatment: toy-calibrated per-bin RMSE override\n')
    else:
        f.write('# stat treatment: nominal Bayes propagated bin errors (no dedicated dN/dy toy coverage yet)\n')
    f.write('bin dndy_center nominal stat binning residual unfolding pid_sf total_sys total\n')
    for ib in range(1, last_active + 1):
        x = nominal.GetXaxis().GetBinCenter(ib)
        nom = nominal.GetBinContent(ib)
        stat = nominal.GetBinError(ib)
        bng = source['binning'][ib - 1]
        res = source['residual'][ib - 1]
        unf = source['unfolding'][ib - 1]
        pid = source['pid_sf'][ib - 1]
        ts = source['total_sys'][ib - 1]
        tt = math.sqrt(stat * stat + ts * ts)
        f.write(f'{ib} {x:.1f} {nom:.6f} {stat:.6f} {bng:.6f} {res:.6f} {unf:.6f} {pid:.6f} {ts:.6f} {tt:.6f}\n')

# root summary
fout = ROOT.TFile.Open(os.path.join(out_dir, 'systematics_dndy_summary.root'), 'RECREATE')
ROOT.TNamed('dNdYStatusBranch', f'keepBinsOverride={ACTIVE_KEEPBINS}; keepBinsUsed={KEEPBINS_USED}; builder=finalize_dndy_systematics.py').Write()
ROOT.TNamed('dNdYKeepBinsUsed', str(KEEPBINS_USED)).Write()
ROOT.TNamed(
    'dNdYStatTreatment',
    'toy_override' if toy_override_used else 'bayes_diagonal_errors_no_dedicated_toy_coverage',
).Write()
nominal.Write('hDoubleRatioNominal_dNdY')
for k, h in h_double.items():
    h.Write(f'hDoubleRatio_{k}_dNdY')
h_sys = nominal.Clone('hSysTotal_dNdY')
h_sys.Reset()
for ib in range(1, nbin + 1):
    h_sys.SetBinContent(ib, source['total_sys'][ib - 1])
    h_sys.SetBinError(ib, 0.0)
h_sys.Write()
fout.Close()

# component plot
c1 = ROOT.TCanvas('c1_dndy', 'c1', 900, 650)
frame = nominal.Clone('hFrame_dndy')
frame.Reset()
frame.SetTitle('Systematic components;dN_{ch}/dy (|y_{T}|<0.5);Absolute uncertainty')
frame.SetMinimum(0.0)
plot_keys = ['binning', 'residual', 'unfolding_prior', 'unfolding_iter', 'unfolding', 'pid_sf', 'total_sys']
ymax = max(max(source[k]) for k in plot_keys)
frame.SetMaximum(0.2)
style_frame(frame)
frame.Draw()
styles = {
    'binning': (ROOT.kGreen + 2, 20),
    'residual': (ROOT.kMagenta + 1, 21),
    'unfolding_prior': (ROOT.kCyan + 2, 26),
    'unfolding_iter': (ROOT.kViolet + 2, 32),
    'unfolding': (ROOT.kOrange + 7, 33),
    'pid_sf': (ROOT.kRed + 1, 25),
    'total_sys': (ROOT.kBlack, 24),
}
leg = ROOT.TLegend(0.36, 0.57, 0.70, 0.89)
style_legend(leg)
leg.SetTextSize(0.034)
keep = []
for key in plot_keys:
    h = nominal.Clone(f'h_{key}_dndy_plot')
    h.Reset()
    col, ms = styles[key]
    for ib in range(1, nbin + 1):
        h.SetBinContent(ib, source[key][ib - 1])
    h.SetLineColor(col)
    h.SetMarkerColor(col)
    h.SetMarkerStyle(ms)
    h.SetLineWidth(2)
    h.Draw('PL SAME')
    label_map = {
        'binning': 'Binning',
        'residual': 'Residual correction (50%)',
        'unfolding_prior': 'Unfolding prior variation',
        'unfolding_iter': 'Unfolding iteration variation',
        'unfolding': 'Bayesian unfolding (envelope)',
        'pid_sf': 'PID SF propagation',
        'total_sys': 'Total syst (v3)',
    }
    leg.AddEntry(h, label_map.get(key, key.replace('_', ' ')), 'lp')
    keep.append(h)
leg.Draw()
c1.SaveAs(os.path.join(out_dir, 'Systematics_Components_vs_dNdY.pdf'))
c1.SaveAs(os.path.join(out_dir, 'Systematics_Components_vs_dNdY.png'))

# final plot
c2 = ROOT.TCanvas('c2_dndy', 'c2', 900, 650)
nominal.SetTitle('Unfolded double ratio vs dN_{ch}/dy;dN_{ch}/dy (|y_{T}|<0.5);(K/#pi)_{Data}/(K/#pi)_{MC}')
nominal.SetMarkerStyle(20)
nominal.SetLineWidth(2)
nominal.SetMinimum(0.6)
nominal.SetMaximum(1.4)
style_frame(nominal)
h_tot = nominal.Clone('hSystematicBand_dndy')
for ib in range(1, nbin + 1):
    sy = source['total_sys'][ib - 1]
    h_tot.SetBinError(ib, sy)
nominal.Draw('E1')
hline = ROOT.TLine(nominal.GetXaxis().GetXmin(), 1.0, nominal.GetXaxis().GetXmax(), 1.0)
hline.SetLineStyle(2)
hline.SetLineWidth(2)
hline.Draw('SAME')
h_tot.SetFillColorAlpha(ROOT.kAzure + 1, 0.30)
h_tot.SetLineColor(ROOT.kAzure + 1)
h_tot.Draw('E2 SAME')
nominal.Draw('E1 SAME')
leg2 = ROOT.TLegend(0.52, 0.74, 0.89, 0.89)
style_legend(leg2)
leg2.AddEntry(nominal, 'Nominal (Bayes, stat)', 'lep')
leg2.AddEntry(h_tot, 'Systematic uncertainty', 'f')
leg2.Draw()
draw_logo(c2)
c2.SaveAs(os.path.join(out_dir, 'KtoPi_DoubleRatio_vs_dNdY_with_Systematics.pdf'))
c2.SaveAs(os.path.join(out_dir, 'KtoPi_DoubleRatio_vs_dNdY_with_Systematics.png'))

print('Wrote dN/dy summary to', out_dir)
