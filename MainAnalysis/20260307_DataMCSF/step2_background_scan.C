#include <TCanvas.h>
#include <TFile.h>
#include <TF1.h>
#include <TFitResultPtr.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TLine.h>
#include <TPaveText.h>
#include <TRandom3.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace {
constexpr double kFitMin = 1.00;
constexpr double kFitMax = 1.05;
constexpr double kSideLo1 = 1.000;
constexpr double kSideHi1 = 1.010;
constexpr double kSideLo2 = 1.028;
constexpr double kSideHi2 = 1.050;

struct Seed {
  std::string model = "TripleGauss";
  double mean = 1.0195;
  double s1 = 0.0025;
  double s2 = 0.0060;
  double s3 = 0.0120;
  double frac1 = 0.5;
  double frac2 = 0.3;
  double frac3 = 0.2;
};

struct FitMetrics {
  std::string category;
  std::string model;
  int nFree = 0;
  int nBinsUsed = 0;
  int fitStatus = -1;
  double chi2 = 1e99;
  double ndf = 1.0;
  double chi2ndf = 1e99;
  double aic = 1e99;
  double bic = 1e99;
  double pullMean = 0.0;
  double pullRms = 0.0;
  double yield = 0.0;
  double yieldErr = 0.0;
  std::vector<double> pars;
};

std::vector<std::string> splitCsv(const std::string& s) {
  std::vector<std::string> out;
  std::stringstream ss(s);
  std::string item;
  while (std::getline(ss, item, ',')) out.push_back(item);
  return out;
}

std::map<std::string, Seed> readSeeds(const std::string& csvPath) {
  std::map<std::string, Seed> out;
  std::ifstream in(csvPath);
  if (!in) return out;
  std::string line;
  std::getline(in, line);
  while (std::getline(in, line)) {
    if (line.empty()) continue;
    auto c = splitCsv(line);
    if (c.size() < 16) continue;
    Seed s;
    s.model = c[1];
    s.mean = std::atof(c[3].c_str());
    s.s1 = std::atof(c[5].c_str());
    s.s2 = std::atof(c[7].c_str());
    s.s3 = std::atof(c[9].c_str());
    if (c.size() > 18) {
      s.frac1 = std::atof(c[16].c_str());
      s.frac2 = std::atof(c[17].c_str());
      s.frac3 = std::atof(c[18].c_str());
    }
    const double fs = s.frac1 + s.frac2 + s.frac3;
    if (fs > 0) {
      s.frac1 /= fs;
      s.frac2 /= fs;
      s.frac3 /= fs;
    }
    out[c[0]] = s;
  }
  return out;
}

std::map<std::string, double> readBaselineYields(const std::string& csvPath) {
  std::map<std::string, double> out;
  std::ifstream in(csvPath);
  if (!in) return out;
  std::string line;
  std::getline(in, line);
  while (std::getline(in, line)) {
    if (line.empty()) continue;
    auto c = splitCsv(line);
    if (c.size() < 6) continue;
    if (c[0] == "0tag" || c[0] == "1tag" || c[0] == "2tag") {
      out[c[0]] = std::atof(c[5].c_str());
    }
  }
  return out;
}

double signalEval(double x, double N, double k, const Seed& s) {
  const double z1 = (x - s.mean) / (k * s.s1);
  const double z2 = (x - s.mean) / (k * s.s2);
  const double z3 = (x - s.mean) / (k * s.s3);
  return N * (s.frac1 * std::exp(-0.5 * z1 * z1) +
              s.frac2 * std::exp(-0.5 * z2 * z2) +
              s.frac3 * std::exp(-0.5 * z3 * z3));
}

int nFreeForModel(const std::string& m) {
  if (m == "pol2_g1") return 8;   // N,k + 3 poly + 3 gaus
  if (m == "pol3_g1") return 9;   // +4 poly
  if (m == "pol4_g1") return 10;  // +5 poly
  if (m == "exp2_g1") return 8;   // +3 exp +3 gaus
  if (m == "pol3_g1_gs") return 12;  // +4 poly + 3 gaus +3 shoulder
  return 999;
}

std::vector<std::string> modelsForCategory(const std::string& cat) {
  if (cat == "0tag") return {"pol3_g1", "pol4_g1", "exp2_g1"};
  if (cat == "1tag") return {"pol3_g1", "pol4_g1"};
  return {"pol2_g1", "pol3_g1"};
}

void fillSidebandGraph(TH1D* h, TGraphErrors& g) {
  int p = 0;
  for (int i = 1; i <= h->GetNbinsX(); ++i) {
    const double x = h->GetBinCenter(i);
    if (!((x >= kSideLo1 && x <= kSideHi1) || (x >= kSideLo2 && x <= kSideHi2))) continue;
    const double y = h->GetBinContent(i);
    const double ey = h->GetBinError(i);
    if (ey <= 0) continue;
    g.SetPoint(p, x, y);
    g.SetPointError(p, 0.0, ey);
    ++p;
  }
}

TF1* makeBkgOnly(const std::string& name, const std::string& tag) {
  if (name == "pol2_g1") return new TF1((tag + "_bkgSB").c_str(), "pol2(0)+gaus(3)", kFitMin, kFitMax);
  if (name == "pol3_g1") return new TF1((tag + "_bkgSB").c_str(), "pol3(0)+gaus(4)", kFitMin, kFitMax);
  if (name == "pol4_g1") return new TF1((tag + "_bkgSB").c_str(), "pol4(0)+gaus(5)", kFitMin, kFitMax);
  if (name == "exp2_g1") return new TF1((tag + "_bkgSB").c_str(), "exp([0]+[1]*x+[2]*x*x)+gaus(3)", kFitMin, kFitMax);
  if (name == "pol3_g1_gs") return new TF1((tag + "_bkgSB").c_str(), "pol3(0)+gaus(4)+gaus(7)", kFitMin, kFitMax);
  return nullptr;
}

TF1* makeTotal(const std::string& name, const Seed& s, const std::string& tag) {
  const std::string sig = Form("[0]*(%.15g*exp(-0.5*((x-%.15g)/(%.15g*[1]))^2)+"
                               "%.15g*exp(-0.5*((x-%.15g)/(%.15g*[1]))^2)+"
                               "%.15g*exp(-0.5*((x-%.15g)/(%.15g*[1]))^2))",
                               s.frac1, s.mean, s.s1,
                               s.frac2, s.mean, s.s2,
                               s.frac3, s.mean, s.s3);
  std::string expr;
  if (name == "pol2_g1") expr = sig + "+pol2(2)+gaus(5)";
  if (name == "pol3_g1") expr = sig + "+pol3(2)+gaus(6)";
  if (name == "pol4_g1") expr = sig + "+pol4(2)+gaus(7)";
  if (name == "exp2_g1") expr = sig + "+exp([2]+[3]*x+[4]*x*x)+gaus(5)";
  if (name == "pol3_g1_gs") expr = sig + "+pol3(2)+gaus(6)+gaus(9)";
  return new TF1((tag + "_tot").c_str(), expr.c_str(), kFitMin, kFitMax);
}

void configureParLimits(const std::string& name, TF1* f, double maxY) {
  f->SetParLimits(0, 0.0, 1e9);     // N
  f->SetParLimits(1, 0.5, 2.0);     // common width scale

  int gidx = -1;
  if (name == "pol2_g1") gidx = 5;
  if (name == "pol3_g1") gidx = 6;
  if (name == "pol4_g1") gidx = 7;
  if (name == "exp2_g1") gidx = 5;
  if (name == "pol3_g1_gs") gidx = 6;

  if (gidx >= 0) {
    f->SetParLimits(gidx + 0, 0.0, 1e9);
    f->SetParLimits(gidx + 1, 0.996, 1.004);
    f->SetParLimits(gidx + 2, 0.0008, 0.01);
  }

  if (name == "pol3_g1_gs") {
    f->SetParLimits(9, 0.0, 1e9);
    f->SetParLimits(10, 1.026, 1.036);
    f->SetParLimits(11, 0.0010, 0.010);
  }
}

void seedFromSideband(const std::string& name, TF1* f, TF1* fb, double maxY) {
  f->SetParameter(0, 0.2 * maxY);
  f->SetParameter(1, 1.0);

  if (name == "pol2_g1") {
    for (int i=0;i<6;++i) f->SetParameter(2+i, fb->GetParameter(i));
  } else if (name == "pol3_g1") {
    for (int i=0;i<7;++i) f->SetParameter(2+i, fb->GetParameter(i));
  } else if (name == "pol4_g1") {
    for (int i=0;i<8;++i) f->SetParameter(2+i, fb->GetParameter(i));
  } else if (name == "exp2_g1") {
    for (int i=0;i<6;++i) f->SetParameter(2+i, fb->GetParameter(i));
  } else if (name == "pol3_g1_gs") {
    for (int i=0;i<10;++i) f->SetParameter(2+i, fb->GetParameter(i));
  }
}

FitMetrics fitCandidate(TH1D* h, const std::string& category, const std::string& name, const Seed& seed, TRandom3& rng) {
  FitMetrics best;
  best.category = category;
  best.model = name;
  best.nFree = nFreeForModel(name);

  TGraphErrors gsb;
  fillSidebandGraph(h, gsb);
  TF1* fb = makeBkgOnly(name, "sb_" + category + "_" + name);

  const double maxY = std::max(10.0, h->GetMaximum());
  const double b0 = 0.5 * (h->GetBinContent(1) + h->GetBinContent(h->GetNbinsX()));

  if (name == "pol2_g1") fb->SetParameters(b0, 0, 0, 0.05*maxY, 1.0, 0.002);
  if (name == "pol3_g1") fb->SetParameters(b0, 0, 0, 0, 0.05*maxY, 1.0, 0.002);
  if (name == "pol4_g1") fb->SetParameters(b0, 0, 0, 0, 0, 0.05*maxY, 1.0, 0.002);
  if (name == "exp2_g1") fb->SetParameters(std::log(std::max(1.0,b0)), 0, 0, 0.05*maxY, 1.0, 0.002);
  if (name == "pol3_g1_gs") fb->SetParameters(b0, 0, 0, 0, 0.05*maxY, 1.0, 0.002, 0.02*maxY, 1.03, 0.004);

  gsb.Fit(fb, "Q0");

  std::vector<double> bestPars;
  const int nTry = (category == "1tag" && name == "pol4_g1") ? 20 : 5;
  for (int itry = 0; itry < nTry; ++itry) {
    TF1* f = makeTotal(name, seed, "fit_" + category + "_" + name + "_" + std::to_string(itry));
    seedFromSideband(name, f, fb, maxY);
    configureParLimits(name, f, maxY);

    // Randomized starts for robustness.
    f->SetParameter(0, std::max(1.0, f->GetParameter(0) * (1.0 + rng.Gaus(0.0, 0.35))));
    f->SetParameter(1, std::max(0.55, std::min(1.8, 1.0 + rng.Gaus(0.0, 0.15))));
    if (category == "1tag" && name == "pol4_g1") {
      // Stabilize the dedicated 1.00 GeV background peak for high-order polynomial fits.
      f->SetParameter(7, std::max(1.0, 0.08 * maxY * (1.0 + rng.Gaus(0.0, 0.25))));
      f->SetParameter(8, 1.000 + rng.Gaus(0.0, 0.0006));
      f->SetParameter(9, std::max(0.0010, std::min(0.0045, 0.0020 + rng.Gaus(0.0, 0.0004))));
    }

    TFitResultPtr fr = h->Fit(f, "RSQ0");
    const int status = int(fr);
    int statusFinal = status;
    if (statusFinal != 0) {
      // Retry with a more robust minimization path when Migrad struggles.
      fr = h->Fit(f, "RMSQ0");
      statusFinal = int(fr);
    }
    if (statusFinal != 0 || f->GetNDF() <= 0) { delete f; continue; }

    const double chi2ndf = f->GetChisquare() / f->GetNDF();
    if (chi2ndf < best.chi2ndf) {
      best.fitStatus = statusFinal;
      best.chi2 = f->GetChisquare();
      best.ndf = f->GetNDF();
      best.chi2ndf = chi2ndf;
      bestPars.resize(f->GetNpar());
      for (int ip = 0; ip < f->GetNpar(); ++ip) bestPars[ip] = f->GetParameter(ip);
    }
    delete f;
  }

  delete fb;

  if (bestPars.empty()) return best;

  TF1* fbest = makeTotal(name, seed, "best_" + category + "_" + name);
  for (int ip = 0; ip < (int)bestPars.size(); ++ip) fbest->SetParameter(ip, bestPars[ip]);

  // Metrics
  int n = 0;
  double pSum = 0.0, p2Sum = 0.0;
  double ySig = 0.0;
  for (int i = 1; i <= h->GetNbinsX(); ++i) {
    const double x = h->GetBinCenter(i);
    if (x < kFitMin || x > kFitMax) continue;
    const double ey = h->GetBinError(i);
    if (ey <= 0) continue;
    const double pull = (h->GetBinContent(i) - fbest->Eval(x)) / ey;
    pSum += pull;
    p2Sum += pull * pull;
    ++n;

    ySig += std::max(0.0, signalEval(x, bestPars[0], bestPars[1], seed));
  }
  best.nBinsUsed = n;
  if (n > 0) {
    best.pullMean = pSum / n;
    const double m2 = p2Sum / n - best.pullMean * best.pullMean;
    best.pullRms = (m2 > 0.0 ? std::sqrt(m2) : 0.0);
  }
  best.aic = best.chi2 + 2.0 * best.nFree;
  best.bic = best.chi2 + best.nFree * std::log(std::max(2, n));
  best.yield = ySig;
  best.yieldErr = std::sqrt(std::max(0.0, ySig));
  best.pars = bestPars;

  delete fbest;
  return best;
}

TF1* makeBestTF1(const FitMetrics& m, const Seed& s) {
  TF1* f = makeTotal(m.model, s, "plot_" + m.category + "_" + m.model);
  for (int ip = 0; ip < (int)m.pars.size(); ++ip) f->SetParameter(ip, m.pars[ip]);
  return f;
}

TF1* makeBestBkg(const FitMetrics& m) {
  TF1* b = nullptr;
  if (m.model == "pol2_g1") {
    b = new TF1(("b_"+m.category).c_str(), "pol2(0)+gaus(3)", kFitMin, kFitMax);
    b->SetParameters(m.pars[2],m.pars[3],m.pars[4],m.pars[5],m.pars[6],m.pars[7]);
  } else if (m.model == "pol3_g1") {
    b = new TF1(("b_"+m.category).c_str(), "pol3(0)+gaus(4)", kFitMin, kFitMax);
    b->SetParameters(m.pars[2],m.pars[3],m.pars[4],m.pars[5],m.pars[6],m.pars[7],m.pars[8]);
  } else if (m.model == "pol4_g1") {
    b = new TF1(("b_"+m.category).c_str(), "pol4(0)+gaus(5)", kFitMin, kFitMax);
    b->SetParameters(m.pars[2],m.pars[3],m.pars[4],m.pars[5],m.pars[6],m.pars[7],m.pars[8],m.pars[9]);
  } else if (m.model == "exp2_g1") {
    b = new TF1(("b_"+m.category).c_str(), "exp([0]+[1]*x+[2]*x*x)+gaus(3)", kFitMin, kFitMax);
    b->SetParameters(m.pars[2],m.pars[3],m.pars[4],m.pars[5],m.pars[6],m.pars[7]);
  } else {
    b = new TF1(("b_"+m.category).c_str(), "pol3(0)+gaus(4)+gaus(7)", kFitMin, kFitMax);
    b->SetParameters(m.pars[2],m.pars[3],m.pars[4],m.pars[5],m.pars[6],m.pars[7],m.pars[8],m.pars[9],m.pars[10],m.pars[11]);
  }
  return b;
}

void drawBestPage(TCanvas& c, TH1D* h, const FitMetrics& m, const Seed& s, double baselineYield, const std::string& pdf) {
  TF1* f = makeBestTF1(m, s);
  TF1* b = makeBestBkg(m);

  c.Clear();
  c.Divide(1,2);
  c.cd(1);
  gPad->SetPad(0.0,0.32,1.0,1.0);
  gPad->SetLeftMargin(0.14);
  gPad->SetRightMargin(0.04);
  gPad->SetBottomMargin(0.03);
  h->SetStats(0);
  h->SetLineWidth(2);
  h->GetXaxis()->SetRangeUser(kFitMin, kFitMax);
  if (m.category == "1tag") {
    double localMax = 0.0;
    for (int i = 1; i <= h->GetNbinsX(); ++i) {
      const double x = h->GetBinCenter(i);
      if (x < 1.012 || x > 1.040) continue;
      localMax = std::max(localMax, h->GetBinContent(i));
    }
    if (localMax > 0.0) {
      h->SetMinimum(0.0);
      h->SetMaximum(1.35 * localMax);
    }
  }
  h->Draw("E");
  f->SetLineColor(kBlue+1);
  f->SetLineWidth(2);
  f->Draw("same");

  TF1 sig(("sig_"+m.category).c_str(), "[0]", kFitMin, kFitMax);
  sig.SetLineColor(kRed+1);
  sig.SetLineWidth(2);
  sig.SetNpx(1000);

  TGraph* gSig = new TGraph();
  int p=0;
  for (int i=0;i<400;++i){
    const double x = kFitMin + (kFitMax-kFitMin)*i/399.0;
    gSig->SetPoint(p++, x, signalEval(x,m.pars[0],m.pars[1],s));
  }
  gSig->SetLineColor(kRed+1);
  gSig->SetLineWidth(2);
  gSig->Draw("L same");

  b->SetLineColor(kGreen+2);
  b->SetLineStyle(2);
  b->SetLineWidth(2);
  b->Draw("same");

  TLegend leg(0.58,0.62,0.89,0.88);
  leg.SetBorderSize(0);
  leg.AddEntry(h, (m.category+" data").c_str(), "lep");
  leg.AddEntry(f, "Total", "l");
  leg.AddEntry(gSig, "Signal", "l");
  leg.AddEntry(b, "Background", "l");
  leg.Draw();

  TPaveText txt(0.12,0.58,0.55,0.88,"NDC");
  txt.SetBorderSize(0);
  txt.SetFillStyle(0);
  txt.AddText(("Best model: "+m.model).c_str());
  txt.AddText(Form("chi2/ndf=%.3f  AIC=%.1f  BIC=%.1f",m.chi2ndf,m.aic,m.bic));
  txt.AddText(Form("Pull mean=%.3f  RMS=%.3f",m.pullMean,m.pullRms));
  txt.AddText(Form("Yield=%.0f",m.yield));
  if (baselineYield > 0) txt.AddText(Form("Shift vs baseline: %+0.2f%%",100.0*(m.yield-baselineYield)/baselineYield));
  txt.Draw();

  c.cd(2);
  gPad->SetPad(0.0,0.0,1.0,0.32);
  gPad->SetLeftMargin(0.14);
  gPad->SetRightMargin(0.04);
  gPad->SetTopMargin(0.03);
  gPad->SetBottomMargin(0.34);
  TH1D hp(("hp_"+m.category).c_str(),"",h->GetNbinsX(),h->GetXaxis()->GetXmin(),h->GetXaxis()->GetXmax());
  hp.SetStats(0);
  hp.SetMarkerStyle(20);
  hp.SetMarkerSize(0.6);
  for (int i=1;i<=h->GetNbinsX();++i){
    const double x = h->GetBinCenter(i);
    if (x<kFitMin||x>kFitMax) continue;
    const double ey = h->GetBinError(i);
    if (ey<=0) continue;
    hp.SetBinContent(i,(h->GetBinContent(i)-f->Eval(x))/ey);
    hp.SetBinError(i,1.0);
  }
  hp.SetMinimum(-5); hp.SetMaximum(5);
  hp.GetYaxis()->SetTitle("Pull ((Obs-Fit)/#sigma)");
  hp.GetXaxis()->SetTitle("m(K^{+}K^{-}) [GeV]");
  hp.GetYaxis()->SetTitleSize(0.11); hp.GetYaxis()->SetLabelSize(0.085);
  hp.GetYaxis()->SetTitleOffset(0.55);
  hp.GetXaxis()->SetTitleSize(0.11); hp.GetXaxis()->SetLabelSize(0.085);
  hp.GetXaxis()->SetTitleOffset(1.05);
  hp.Draw("E1P");
  TLine l0(kFitMin,0,kFitMax,0); l0.SetLineStyle(2); l0.Draw();

  c.Print(pdf.c_str());

  delete gSig;
  delete f;
  delete b;
}

}  // namespace

void step2_background_scan(const char* histsFile = "Reports/step2/step2_hists.root",
                           const char* seedCsv = "Reports/step1_best_fit_params.csv",
                           const char* baselineCsv = "Reports/step2_fit_results.csv",
                           const char* outDir = "Reports/step2_bgscan") {
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);
  gSystem->mkdir(outDir, kTRUE);

  TFile f(histsFile, "READ");
  if (f.IsZombie()) {
    std::cerr << "Cannot open: " << histsFile << std::endl;
    return;
  }

  TH1D* h0 = dynamic_cast<TH1D*>(f.Get("hStep2_0tag"));
  TH1D* h1 = dynamic_cast<TH1D*>(f.Get("hStep2_1tag"));
  TH1D* h2 = dynamic_cast<TH1D*>(f.Get("hStep2_2tag"));
  if (!h0 || !h1 || !h2) {
    std::cerr << "Missing hStep2_* histograms in " << histsFile << std::endl;
    return;
  }

  auto seeds = readSeeds(seedCsv);
  auto baseline = readBaselineYields(baselineCsv);

  std::map<std::string, TH1D*> hmap{{"0tag",h0},{"1tag",h1},{"2tag",h2}};
  std::map<std::string, FitMetrics> bestByCat;
  std::vector<FitMetrics> all;

  TRandom3 rng(12345);

  for (const auto& cat : {std::string("0tag"), std::string("1tag"), std::string("2tag")}) {
    Seed s = seeds.count(cat) ? seeds[cat] : Seed{};
    auto models = modelsForCategory(cat);
    FitMetrics best;
    best.bic = 1e99;

    for (const auto& m : models) {
      FitMetrics r = fitCandidate(hmap[cat], cat, m, s, rng);
      all.push_back(r);
      if (r.fitStatus == 0 && std::isfinite(r.bic) && r.bic < best.bic) best = r;
    }

    if (!std::isfinite(best.bic)) {
      std::cerr << "No converged fit for category " << cat << std::endl;
      continue;
    }
    bestByCat[cat] = best;
  }

  const std::string outDirStr(outDir);
  const std::string csvPath = outDirStr + "/step2_background_model_scan.csv";
  const std::string txtPath = outDirStr + "/step2_background_model_scan.txt";
  const std::string tunedPath = outDirStr + "/step2_tuned_fit_summary.txt";
  const std::string pdfPath = outDirStr + "/step2_tuned_report.pdf";

  std::ofstream csv(csvPath);
  csv << "category,model,fit_status,n_free,n_bins,chi2,ndf,chi2ndf,aic,bic,pull_mean,pull_rms,yield,yield_err\n";
  for (const auto& r : all) {
    csv << r.category << "," << r.model << "," << r.fitStatus << "," << r.nFree << "," << r.nBinsUsed
        << "," << r.chi2 << "," << r.ndf << "," << r.chi2ndf << "," << r.aic << "," << r.bic
        << "," << r.pullMean << "," << r.pullRms << "," << r.yield << "," << r.yieldErr << "\n";
  }
  csv.close();

  std::ofstream txt(txtPath);
  txt << std::fixed << std::setprecision(6);
  txt << "Step2 background scan results\n";
  txt << "Input hists: " << histsFile << "\n";
  txt << "Seed CSV: " << seedCsv << "\n\n";
  for (const auto& cat : {std::string("0tag"), std::string("1tag"), std::string("2tag")}) {
    txt << "Category " << cat << "\n";
    std::vector<FitMetrics> v;
    for (const auto& r : all) if (r.category == cat) v.push_back(r);
    std::sort(v.begin(), v.end(), [](const FitMetrics& a, const FitMetrics& b){ return a.bic < b.bic; });
    for (const auto& r : v) {
      txt << "  " << r.model << "  status=" << r.fitStatus << "  chi2/ndf=" << r.chi2ndf
          << "  AIC=" << r.aic << "  BIC=" << r.bic << "  pullRMS=" << r.pullRms
          << "  yield=" << r.yield << "\n";
    }
    if (bestByCat.count(cat)) {
      const auto& b = bestByCat[cat];
      txt << "  -> selected: " << b.model << " (BIC=" << b.bic << ")\n\n";
    }
  }
  txt.close();

  std::ofstream tuned(tunedPath);
  tuned << std::fixed << std::setprecision(6);
  tuned << "category,selected_model,chi2ndf,aic,bic,pull_mean,pull_rms,yield,yield_err,baseline_yield,yield_shift_percent\n";
  for (const auto& cat : {std::string("0tag"), std::string("1tag"), std::string("2tag")}) {
    if (!bestByCat.count(cat)) continue;
    const auto& b = bestByCat[cat];
    const double base = baseline.count(cat) ? baseline[cat] : 0.0;
    const double shift = (base > 0.0) ? 100.0 * (b.yield - base) / base : 0.0;
    tuned << cat << "," << b.model << "," << b.chi2ndf << "," << b.aic << "," << b.bic << ","
          << b.pullMean << "," << b.pullRms << "," << b.yield << "," << b.yieldErr << ","
          << base << "," << shift << "\n";
  }
  tuned.close();

  TCanvas c("c_scan", "c_scan", 1100, 800);
  c.Print((pdfPath + "[").c_str());

  c.Clear();
  TPaveText cover(0.06,0.08,0.94,0.92,"NDC");
  cover.SetBorderSize(0);
  cover.SetFillStyle(0);
  cover.SetTextAlign(12);
  cover.AddText("Step 2 Background Shape Optimization");
  cover.AddText(Form("Input histograms: %s", histsFile));
  cover.AddText(Form("Signal seed (fixed shape): %s", seedCsv));
  cover.AddText("Selection objective: lowest BIC with stable convergence");
  cover.AddText("Robustness: 5 randomized refits per candidate");
  cover.Draw();
  c.Print(pdfPath.c_str());

  for (const auto& cat : {std::string("0tag"), std::string("1tag"), std::string("2tag")}) {
    c.Clear();
    TPaveText p(0.05,0.05,0.95,0.95,"NDC");
    p.SetBorderSize(0);
    p.SetFillStyle(0);
    p.SetTextAlign(12);
    p.AddText(("Category " + cat + " model ranking").c_str());
    p.AddText("model | status | chi2/ndf | AIC | BIC | pullMean | pullRMS | yield");

    std::vector<FitMetrics> v;
    for (const auto& r : all) if (r.category == cat) v.push_back(r);
    std::sort(v.begin(), v.end(), [](const FitMetrics& a, const FitMetrics& b){ return a.bic < b.bic; });

    for (const auto& r : v) {
      p.AddText(Form("%s | %d | %.3f | %.1f | %.1f | %.3f | %.3f | %.0f",
                     r.model.c_str(), r.fitStatus, r.chi2ndf, r.aic, r.bic, r.pullMean, r.pullRms, r.yield));
    }
    if (bestByCat.count(cat)) p.AddText(("Selected: " + bestByCat[cat].model).c_str());
    p.Draw();
    c.Print(pdfPath.c_str());
  }

  for (const auto& cat : {std::string("0tag"), std::string("1tag"), std::string("2tag")}) {
    if (!bestByCat.count(cat)) continue;
    Seed s = seeds.count(cat) ? seeds[cat] : Seed{};
    const double base = baseline.count(cat) ? baseline[cat] : 0.0;
    drawBestPage(c, hmap[cat], bestByCat[cat], s, base, pdfPath);
  }

  c.Print((pdfPath + "]").c_str());

  std::cout << "Wrote scan table: " << csvPath << std::endl;
  std::cout << "Wrote scan text: " << txtPath << std::endl;
  std::cout << "Wrote tuned summary: " << tunedPath << std::endl;
  std::cout << "Wrote tuned report: " << pdfPath << std::endl;
}
