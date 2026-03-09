#include <TCanvas.h>
#include <TFile.h>
#include <TF1.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <TPad.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TTree.h>

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

namespace {
constexpr double kPionMass = 0.13957039;
constexpr double kAbsCosMin = 0.15;
constexpr double kAbsCosMax = 0.675;
constexpr int kMaxReco = 256;
constexpr int kMaxKShort = 128;
constexpr int kBins = 280;

double buildMass(double px1, double py1, double pz1,
                 double px2, double py2, double pz2) {
  const double p1sq = px1 * px1 + py1 * py1 + pz1 * pz1;
  const double p2sq = px2 * px2 + py2 * py2 + pz2 * pz2;
  const double e1 = std::sqrt(p1sq + kPionMass * kPionMass);
  const double e2 = std::sqrt(p2sq + kPionMass * kPionMass);
  const double e = e1 + e2;
  const double px = px1 + px2;
  const double py = py1 + py2;
  const double pz = pz1 + pz2;
  const double m2 = e * e - (px * px + py * py + pz * pz);
  return (m2 > 0.0 ? std::sqrt(m2) : 0.0);
}

bool passTrackAcceptance(double px, double py, double pz) {
  const double p = std::sqrt(px * px + py * py + pz * pz);
  if (p <= 0.0) return false;
  const double abscos = std::fabs(pz / p);
  return (abscos >= kAbsCosMin && abscos <= kAbsCosMax);
}

struct FitPage {
  std::string model;
  double chi2ndf = 1e9;
  TF1* fn = nullptr;
};

TF1* buildModel(const std::string& name, const std::string& fname,
                double fitMin, double fitMax, double maxY, double b0) {
  TF1* f = nullptr;
  if (name == "Gauss") {
    f = new TF1(fname.c_str(), "gaus(0)", fitMin, fitMax);
    f->SetParameters(maxY, 0.4976, 0.006);
    f->SetParLimits(1, 0.485, 0.510);
    f->SetParLimits(2, 0.001, 0.04);
  } else if (name == "Voigt") {
    f = new TF1(fname.c_str(), "[0]*TMath::Voigt(x-[1],[2],[3],4)", fitMin, fitMax);
    f->SetParameters(maxY, 0.4976, 0.002, 0.006);
    f->SetParLimits(1, 0.485, 0.510);
    f->SetParLimits(2, 0.0005, 0.03);
    f->SetParLimits(3, 0.001, 0.03);
  } else if (name == "DoubleGauss") {
    f = new TF1(fname.c_str(),
                "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]*exp(-0.5*((x-[1])/[4])^2)",
                fitMin, fitMax);
    f->SetParameters(0.65 * maxY, 0.4976, 0.004, 0.35 * maxY, 0.012);
    f->SetParLimits(1, 0.485, 0.510);
    f->SetParLimits(2, 0.001, 0.04);
    f->SetParLimits(4, 0.001, 0.08);
  } else if (name == "TripleGauss") {
    f = new TF1(fname.c_str(),
                "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]*exp(-0.5*((x-[1])/[4])^2)+"
                "[5]*exp(-0.5*((x-[1])/[6])^2)",
                fitMin, fitMax);
    f->SetParameters(0.5 * maxY, 0.4976, 0.003, 0.3 * maxY, 0.009, 0.2 * maxY, 0.03);
    f->SetParLimits(1, 0.485, 0.510);
    f->SetParLimits(2, 0.001, 0.04);
    f->SetParLimits(4, 0.001, 0.08);
    f->SetParLimits(6, 0.001, 0.12);
  } else if (name == "QuadGauss") {
    f = new TF1(fname.c_str(),
                "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]*exp(-0.5*((x-[1])/[4])^2)+"
                "[5]*exp(-0.5*((x-[1])/[6])^2)+[7]*exp(-0.5*((x-[1])/[8])^2)",
                fitMin, fitMax);
    f->SetParameters(0.42 * maxY, 0.4976, 0.003, 0.28 * maxY, 0.008, 0.20 * maxY, 0.02, 0.10 * maxY, 0.05);
    f->SetParLimits(1, 0.485, 0.510);
    f->SetParLimits(2, 0.001, 0.04);
    f->SetParLimits(4, 0.001, 0.08);
    f->SetParLimits(6, 0.001, 0.12);
    f->SetParLimits(8, 0.001, 0.20);
  } else if (name == "BifurGauss") {
    f = new TF1(fname.c_str(),
                "[0]*(x<[1]?exp(-0.5*((x-[1])/[2])^2):exp(-0.5*((x-[1])/[3])^2))",
                fitMin, fitMax);
    f->SetParameters(maxY, 0.4976, 0.004, 0.012);
    f->SetParLimits(1, 0.485, 0.510);
    f->SetParLimits(2, 0.001, 0.05);
    f->SetParLimits(3, 0.001, 0.08);
  } else if (name == "VoigtPlusGauss") {
    f = new TF1(fname.c_str(),
                "[0]*TMath::Voigt(x-[1],[2],[3],4)+gaus(4)",
                fitMin, fitMax);
    f->SetParameters(0.7 * maxY, 0.4976, 0.002, 0.006, 0.3 * maxY, 0.4976, 0.012);
    f->SetParLimits(1, 0.485, 0.510);
    f->SetParLimits(2, 0.0005, 0.03);
    f->SetParLimits(3, 0.001, 0.03);
    f->SetParLimits(5, 0.485, 0.510);
    f->SetParLimits(6, 0.001, 0.08);
  }
  if (f) {
    f->SetLineColor(kRed + 1);
    f->SetLineWidth(2);
  }
  return f;
}

std::vector<FitPage> fitAllModels(TH1D* h, const std::string& cat,
                                  const std::vector<std::string>& models,
                                  double fitMin, double fitMax) {
  std::vector<FitPage> out;
  const double maxY = std::max(100.0, h->GetMaximum());
  const double b0 = 0.5 * (h->GetBinContent(1) + h->GetBinContent(h->GetNbinsX()));
  for (const auto& m : models) {
    FitPage p;
    p.model = m;
    p.fn = buildModel(m, "f_" + cat + "_" + m, fitMin, fitMax, maxY, b0);
    if (!p.fn) continue;
    TH1D* hc = (TH1D*)h->Clone(("hc_" + cat + "_" + m).c_str());
    hc->Fit(p.fn, "RQ0");
    if (p.fn->GetNDF() > 0) p.chi2ndf = p.fn->GetChisquare() / p.fn->GetNDF();
    delete hc;
    out.push_back(p);
  }
  return out;
}
}  // namespace

void step4_kshort_model_scan_report(const char* input = "Samples/merged_mc_v2.2.root",
                                    const char* outPdf = "Reports/step4_modelscan_fr030_100_peaky_report.pdf",
                                    double angleCut = 0.025,
                                    double massMin = 0.30,
                                    double massMax = 1.00,
                                    double fitMin = 0.30,
                                    double fitMax = 1.00) {
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  TFile f(input, "READ");
  if (f.IsZombie()) {
    std::cerr << "Failed to open input: " << input << std::endl;
    return;
  }
  TTree* tree = nullptr;
  f.GetObject("Tree", tree);
  if (!tree) {
    std::cerr << "Cannot find Tree in " << input << std::endl;
    return;
  }

  Long64_t nReco = 0, nKShort = 0;
  double recoPx[kMaxReco], recoPy[kMaxReco], recoPz[kMaxReco], recoCharge[kMaxReco];
  Long64_t recoGoodTrack[kMaxReco], recoPIDPion[kMaxReco];
  Long64_t ksReco1ID[kMaxKShort], ksReco2ID[kMaxKShort];
  double ksReco1Angle[kMaxKShort], ksReco2Angle[kMaxKShort];

  tree->SetBranchAddress("NReco", &nReco);
  tree->SetBranchAddress("NKShort", &nKShort);
  tree->SetBranchAddress("RecoPx", recoPx);
  tree->SetBranchAddress("RecoPy", recoPy);
  tree->SetBranchAddress("RecoPz", recoPz);
  tree->SetBranchAddress("RecoCharge", recoCharge);
  tree->SetBranchAddress("RecoGoodTrack", recoGoodTrack);
  tree->SetBranchAddress("RecoPIDPion", recoPIDPion);
  tree->SetBranchAddress("KShortReco1ID[NKShort]", ksReco1ID);
  tree->SetBranchAddress("KShortReco2ID[NKShort]", ksReco2ID);
  tree->SetBranchAddress("KShortReco1Angle[NKShort]", ksReco1Angle);
  tree->SetBranchAddress("KShortReco2Angle[NKShort]", ksReco2Angle);

  TH1D h0("hScan0", "0-tag", kBins, massMin, massMax);
  TH1D h1("hScan1", "1-tag", kBins, massMin, massMax);
  TH1D h2("hScan2", "2-tag", kBins, massMin, massMax);

  const Long64_t nEntries = tree->GetEntries();
  for (Long64_t ie = 0; ie < nEntries; ++ie) {
    tree->GetEntry(ie);
    if (nReco > kMaxReco || nKShort > kMaxKShort) continue;
    for (Long64_t i = 0; i < nKShort; ++i) {
      const Long64_t i1 = ksReco1ID[i];
      const Long64_t i2 = ksReco2ID[i];
      if (i1 < 0 || i2 < 0 || i1 >= nReco || i2 >= nReco) continue;
      if (!(ksReco1Angle[i] < angleCut && ksReco2Angle[i] < angleCut)) continue;
      if (recoGoodTrack[i1] == 0 || recoGoodTrack[i2] == 0) continue;
      if (!passTrackAcceptance(recoPx[i1], recoPy[i1], recoPz[i1])) continue;
      if (!passTrackAcceptance(recoPx[i2], recoPy[i2], recoPz[i2])) continue;
      if (recoCharge[i1] * recoCharge[i2] >= 0) continue;
      const double m = buildMass(recoPx[i1], recoPy[i1], recoPz[i1], recoPx[i2], recoPy[i2], recoPz[i2]);
      if (m < massMin || m > massMax) continue;
      const int nTag = (recoPIDPion[i1] >= 2 ? 1 : 0) + (recoPIDPion[i2] >= 2 ? 1 : 0);
      if (nTag == 0) h0.Fill(m);
      if (nTag == 1) h1.Fill(m);
      if (nTag == 2) h2.Fill(m);
    }
  }

  std::vector<std::string> models = {
      "Gauss", "Voigt", "DoubleGauss", "TripleGauss", "QuadGauss", "BifurGauss", "VoigtPlusGauss"};

  auto p0 = fitAllModels(&h0, "0tag", models, fitMin, fitMax);
  auto p1 = fitAllModels(&h1, "1tag", models, fitMin, fitMax);
  auto p2 = fitAllModels(&h2, "2tag", models, fitMin, fitMax);

  TCanvas c("c_modelscan", "c_modelscan", 1100, 850);
  c.Print((std::string(outPdf) + "[").c_str());

  auto drawModelWithPull = [&](TH1D* h, const FitPage& page, const std::string& cat) {
    c.Clear();
    TPad pTop("pTop", "pTop", 0.0, 0.32, 1.0, 1.0);
    TPad pBot("pBot", "pBot", 0.0, 0.0, 1.0, 0.32);
    pTop.SetLeftMargin(0.12);
    pTop.SetRightMargin(0.03);
    pTop.SetBottomMargin(0.02);
    pBot.SetLeftMargin(0.12);
    pBot.SetRightMargin(0.03);
    pBot.SetTopMargin(0.04);
    pBot.SetBottomMargin(0.33);
    pTop.Draw();
    pBot.Draw();

    pTop.cd();
    h->SetStats(0);
    h->SetLineColor(kBlack);
    h->SetMarkerStyle(20);
    h->SetMarkerSize(0.6);
    h->GetYaxis()->SetTitle("Candidates / bin");
    h->SetTitle((cat + " : " + page.model).c_str());
    h->Draw("E");
    if (page.fn) page.fn->Draw("same");
    TLatex t;
    t.SetNDC();
    t.SetTextSize(0.045);
    t.DrawLatex(0.14, 0.86, Form("#chi^{2}/ndf = %.3f", page.chi2ndf));
    TLegend leg(0.60, 0.76, 0.96, 0.92);
    leg.SetBorderSize(0);
    leg.AddEntry(h, "MC points", "lep");
    leg.AddEntry(page.fn, page.model.c_str(), "l");
    leg.Draw();

    pBot.cd();
    TH1D hPull(("hPull_" + cat + "_" + page.model).c_str(), "", h->GetNbinsX(),
               h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());
    hPull.SetStats(0);
    hPull.SetMarkerStyle(20);
    hPull.SetMarkerSize(0.5);
    hPull.SetLineColor(kBlack);
    for (int i = 1; i <= h->GetNbinsX(); ++i) {
      const double x = h->GetBinCenter(i);
      if (x < fitMin || x > fitMax) continue;
      const double y = h->GetBinContent(i);
      const double ey = h->GetBinError(i);
      if (ey <= 0.0 || !page.fn) continue;
      const double yfit = page.fn->Eval(x);
      hPull.SetBinContent(i, (y - yfit) / ey);
      hPull.SetBinError(i, 1.0);
    }
    hPull.GetYaxis()->SetTitle("Pull ((Obs-Fit)/#sigma)");
    hPull.GetXaxis()->SetTitle("m(#pi^{+}#pi^{-}) [GeV]");
    hPull.GetYaxis()->SetNdivisions(505);
    hPull.GetYaxis()->SetTitleSize(0.11);
    hPull.GetYaxis()->SetLabelSize(0.085);
    hPull.GetYaxis()->SetTitleOffset(0.50);
    hPull.GetXaxis()->SetTitleSize(0.11);
    hPull.GetXaxis()->SetLabelSize(0.085);
    hPull.GetXaxis()->SetTitleOffset(1.0);
    hPull.SetMinimum(-5.0);
    hPull.SetMaximum(5.0);
    hPull.Draw("E1P");
    TLine l0(massMin, 0.0, massMax, 0.0);
    l0.SetLineStyle(2);
    l0.Draw();

    c.Print(outPdf);
  };

  for (const auto& p : p0) drawModelWithPull(&h0, p, "0-tag");
  for (const auto& p : p1) drawModelWithPull(&h1, p, "1-tag");
  for (const auto& p : p2) drawModelWithPull(&h2, p, "2-tag");

  c.Print((std::string(outPdf) + "]").c_str());

  for (auto& p : p0) delete p.fn;
  for (auto& p : p1) delete p.fn;
  for (auto& p : p2) delete p.fn;

  std::cout << "Wrote " << outPdf << std::endl;
}
