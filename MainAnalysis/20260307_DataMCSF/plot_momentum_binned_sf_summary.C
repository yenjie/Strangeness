#include <TBox.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH1F.h>
#include <TLegend.h>
#include <TLine.h>
#include <TROOT.h>
#include <TStyle.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace {
struct SFBin {
  std::string label;
  double pmin = 0.0;
  double pmax = 0.0;
  double sf = 0.0;
  double stat = 0.0;
};

std::vector<std::string> splitCsv(const std::string& s) {
  std::vector<std::string> out;
  std::stringstream ss(s);
  std::string item;
  while (std::getline(ss, item, ',')) out.push_back(item);
  return out;
}

std::vector<SFBin> readSFBins(const std::string& path) {
  std::vector<SFBin> out;
  std::ifstream in(path);
  if (!in) return out;
  std::string line;
  if (!std::getline(in, line)) return out;
  while (std::getline(in, line)) {
    if (line.empty()) continue;
    auto c = splitCsv(line);
    if (c.size() < 15) continue;
    SFBin b;
    b.label = c[0];
    b.pmin = std::atof(c[1].c_str());
    b.pmax = std::atof(c[2].c_str());
    b.sf = std::atof(c[13].c_str());
    b.stat = std::atof(c[14].c_str());
    out.push_back(b);
  }
  return out;
}

std::map<std::string, double> readSystematicsByBin(const std::string& path) {
  std::ifstream in(path);
  std::map<std::string, std::map<std::string, double>> maxByBinSource;
  if (!in) return {};
  std::string line;
  if (!std::getline(in, line)) return {};
  while (std::getline(in, line)) {
    if (line.empty()) continue;
    auto c = splitCsv(line);
    if (c.size() < 20) continue;
    const std::string bin = c[0];
    const std::string source = c[3];
    if (source == "baseline") continue;
    const double d = std::fabs(std::atof(c[18].c_str()));
    auto& ref = maxByBinSource[bin][source];
    ref = std::max(ref, d);
  }

  std::map<std::string, double> sysByBin;
  for (const auto& kv : maxByBinSource) {
    double q2 = 0.0;
    for (const auto& s : kv.second) q2 += s.second * s.second;
    sysByBin[kv.first] = std::sqrt(q2);
  }
  return sysByBin;
}

void drawSummary(const std::vector<SFBin>& bins,
                 const std::map<std::string, double>& sysByBin,
                 const char* xtitle,
                 const char* title,
                 double sfInclusive,
                 double sfInclusiveTotErr,
                 const char* outPdf) {
  if (bins.empty()) return;

  gStyle->SetOptStat(0);
  TCanvas c("c", "c", 900, 600);
  c.SetLeftMargin(0.12);
  c.SetRightMargin(0.04);
  c.SetBottomMargin(0.12);
  c.SetTopMargin(0.08);

  double xmin = 1e9, xmax = -1e9;
  double ymin = 1e9, ymax = -1e9;
  for (const auto& b : bins) {
    xmin = std::min(xmin, b.pmin);
    xmax = std::max(xmax, b.pmax);
    const double sys = sysByBin.count(b.label) ? sysByBin.at(b.label) : 0.0;
    ymin = std::min(ymin, b.sf - std::max(b.stat, sys));
    ymax = std::max(ymax, b.sf + std::max(b.stat, sys));
  }
  const double pad = 0.08 * std::max(1e-6, ymax - ymin);
  ymin -= pad;
  ymax += pad;

  TH1F* frame = c.DrawFrame(xmin, ymin, xmax, ymax);
  frame->SetTitle(title);
  frame->GetXaxis()->SetTitle(xtitle);
  frame->GetYaxis()->SetTitle("Scale Factor");

  // Momentum-independent SF band: central value +/- total uncertainty.
  TBox* incBand = new TBox(xmin, sfInclusive - sfInclusiveTotErr, xmax, sfInclusive + sfInclusiveTotErr);
  incBand->SetFillColor(kGreen + 1);
  incBand->SetFillStyle(3002);
  incBand->SetLineColor(kGreen + 3);
  incBand->Draw("same");

  TLegend leg(0.58, 0.78, 0.92, 0.92);
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
  leg.AddEntry(incBand, "Inclusive SF #pm total unc.", "f");

  TBox sampleBox;
  sampleBox.SetFillColor(kYellow);
  sampleBox.SetFillStyle(3001);
  sampleBox.SetLineColor(kYellow + 2);
  bool addedBox = false;

  std::vector<double> xc, yc, ex, ey;
  for (const auto& b : bins) {
    const double sys = sysByBin.count(b.label) ? sysByBin.at(b.label) : 0.0;
    TBox* box = new TBox(b.pmin, b.sf - sys, b.pmax, b.sf + sys);
    box->SetFillColor(kYellow);
    box->SetFillStyle(3001);
    box->SetLineColor(kYellow + 2);
    box->Draw("same");
    if (!addedBox) {
      leg.AddEntry(box, "Systematic uncertainty", "f");
      addedBox = true;
    }

    TLine* h = new TLine(b.pmin, b.sf, b.pmax, b.sf);
    h->SetLineColor(kBlack);
    h->SetLineWidth(1);
    h->Draw("same");

    xc.push_back(0.5 * (b.pmin + b.pmax));
    yc.push_back(b.sf);
    ex.push_back(0.0);
    ey.push_back(b.stat);
  }

  TGraphErrors gr((int)xc.size());
  for (int i = 0; i < (int)xc.size(); ++i) {
    gr.SetPoint(i, xc[i], yc[i]);
    gr.SetPointError(i, ex[i], ey[i]);
  }
  gr.SetMarkerStyle(20);
  gr.SetMarkerSize(1.0);
  gr.SetLineColor(kBlack);
  gr.SetLineWidth(2);
  gr.Draw("P same");
  leg.AddEntry(&gr, "Statistical uncertainty", "lep");

  leg.Draw();
  c.SaveAs(outPdf);
}
}  // namespace

void plot_momentum_binned_sf_summary() {
  gROOT->SetBatch(kTRUE);

  auto phiBins = readSFBins("Reports/step23_phi_momentum_binned_sf.csv");
  auto phiSys = readSystematicsByBin("Reports/step23_phi_momentum_binned_systematics_qc_curated.csv");
  const double phiSfInclusive = 0.893777;
  const double phiTot = std::sqrt(0.008899 * 0.008899 + 0.042283 * 0.042283);
  drawSummary(phiBins, phiSys,
              "Kaon momentum p_{K} [GeV]",
              "#phi#rightarrowK^{+}K^{-} momentum-binned SF",
              phiSfInclusive,
              phiTot,
              "Reports/step23_phi_momentum_binned_sf_summary_plot.pdf");

  auto ksBins = readSFBins("Reports/step56_kshort_momentum_binned_sf_quadgauss.csv");
  auto ksSys = readSystematicsByBin("Reports/step56_kshort_momentum_binned_systematics_quadgauss_qc_curated.csv");
  const double ksSfInclusive = 1.019310;
  const double ksTot = std::sqrt(0.007145 * 0.007145 + 0.093109 * 0.093109);
  drawSummary(ksBins, ksSys,
              "Pion momentum p_{#pi} [GeV]",
              "K^{0}_{S}#rightarrow#pi^{+}#pi^{-} momentum-binned SF",
              ksSfInclusive,
              ksTot,
              "Reports/step56_kshort_momentum_binned_sf_summary_plot.pdf");
}
