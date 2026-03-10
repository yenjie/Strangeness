#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPad.h"
#include "TPaveText.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"

namespace {
double gFitMin = 0.70;
double gFitMax = 1.10;
constexpr int kDisplayRebin = 4;
constexpr double kXRef = 0.82;

struct FitResult {
  std::string category;
  std::string model;
  double chi2 = 0.0;
  double ndf = 0.0;
  double chi2ndf = 0.0;
  bool converged = false;
};

std::string getArgument(int argc, char* argv[], const std::string& option,
                        const std::string& defaultValue) {
  for (int i = 1; i + 1 < argc; ++i)
    if (argv[i] == option) return argv[i + 1];
  return defaultValue;
}

double getDoubleArgument(int argc, char* argv[], const std::string& option, double defaultValue) {
  const std::string value = getArgument(argc, argv, option, "");
  return value.empty() ? defaultValue : std::stod(value);
}

double logisticRise(double x, double x0, double w) {
  return 1.0 / (1.0 + std::exp(-(x - x0) / w));
}

double logisticFall(double x, double x0, double w) {
  return 1.0 / (1.0 + std::exp((x - x0) / w));
}

double DoubleLogisticExp(double* x, double* p) {
  const double xx = x[0];
  const double env = p[0] * logisticRise(xx, p[1], p[2]) * logisticFall(xx, p[3], p[4]);
  return env * std::exp(p[5] * (xx - kXRef));
}

double DoubleLogisticPoly1(double* x, double* p) {
  const double xx = x[0];
  const double dx = xx - kXRef;
  const double env = p[0] * logisticRise(xx, p[1], p[2]) * logisticFall(xx, p[3], p[4]);
  const double poly = 1.0 + p[5] * dx;
  return env * std::max(0.0, poly);
}

double DoubleLogisticPoly2(double* x, double* p) {
  const double xx = x[0];
  const double dx = xx - kXRef;
  const double env = p[0] * logisticRise(xx, p[1], p[2]) * logisticFall(xx, p[3], p[4]);
  const double poly = 1.0 + p[5] * dx + p[6] * dx * dx;
  return env * std::max(0.0, poly);
}

double DoubleLogisticExpPoly1(double* x, double* p) {
  const double xx = x[0];
  const double dx = xx - kXRef;
  const double env = p[0] * logisticRise(xx, p[1], p[2]) * logisticFall(xx, p[3], p[4]);
  const double poly = 1.0 + p[6] * dx;
  return env * std::exp(p[5] * dx) * std::max(0.0, poly);
}

double SharpRiseLogisticExp(double* x, double* p) {
  const double xx = x[0];
  const double rise = std::pow(logisticRise(xx, p[1], p[2]), p[5]);
  const double fall = logisticFall(xx, p[3], p[4]);
  return p[0] * rise * fall * std::exp(p[6] * (xx - kXRef));
}

double SharpRiseLogisticPoly1(double* x, double* p) {
  const double xx = x[0];
  const double dx = xx - kXRef;
  const double rise = std::pow(logisticRise(xx, p[1], p[2]), p[5]);
  const double fall = logisticFall(xx, p[3], p[4]);
  const double poly = 1.0 + p[6] * dx;
  return p[0] * rise * fall * std::max(0.0, poly);
}

double DoubleLogisticExpPlusGauss(double* x, double* p) {
  const double xx = x[0];
  const double dx = xx - kXRef;
  const double env = p[0] * logisticRise(xx, p[1], p[2]) * logisticFall(xx, p[3], p[4]);
  const double base = env * std::exp(p[5] * dx);
  const double gauss = p[6] * std::exp(-0.5 * std::pow((xx - p[7]) / p[8], 2));
  return base + gauss;
}

double TwoDoubleLogisticExp(double* x, double* p) {
  const double xx = x[0];
  const double dx = xx - kXRef;
  const double c1 = p[0] * logisticRise(xx, p[1], p[2]) * logisticFall(xx, p[3], p[4]) * std::exp(p[5] * dx);
  const double c2 = p[6] * logisticRise(xx, p[7], p[8]) * logisticFall(xx, p[9], p[10]) * std::exp(p[11] * dx);
  return c1 + c2;
}

double Poly1Exp(double* x, double* p) {
  const double dx = x[0] - kXRef;
  return p[0] * (1.0 + p[1] * dx) * std::exp(p[2] * dx);
}

double Poly2Exp(double* x, double* p) {
  const double dx = x[0] - kXRef;
  return p[0] * (1.0 + p[1] * dx + p[2] * dx * dx) * std::exp(p[3] * dx);
}

TF1 makeModel(const std::string& name, const std::string& model, TH1D* h) {
  const double maxY = std::max(1.0, h->GetMaximum());

  if (model == "DoubleLogisticExp") {
    TF1 f(name.c_str(), DoubleLogisticExp, gFitMin, gFitMax, 6);
    f.SetParameters(maxY * 1.2, 0.73, 0.02, 0.92, 0.08, -4.0);
    f.SetParLimits(0, 0.0, 1e9);
    f.SetParLimits(1, 0.68, 0.82);
    f.SetParLimits(2, 0.002, 0.10);
    f.SetParLimits(3, 0.76, 1.08);
    f.SetParLimits(4, 0.005, 0.25);
    f.SetParLimits(5, -40.0, 20.0);
    return f;
  }
  if (model == "DoubleLogisticPoly1") {
    TF1 f(name.c_str(), DoubleLogisticPoly1, gFitMin, gFitMax, 6);
    f.SetParameters(maxY * 1.2, 0.73, 0.02, 0.92, 0.08, -3.0);
    f.SetParLimits(0, 0.0, 1e9);
    f.SetParLimits(1, 0.68, 0.82);
    f.SetParLimits(2, 0.002, 0.10);
    f.SetParLimits(3, 0.76, 1.08);
    f.SetParLimits(4, 0.005, 0.25);
    f.SetParLimits(5, -80.0, 80.0);
    return f;
  }
  if (model == "DoubleLogisticPoly2") {
    TF1 f(name.c_str(), DoubleLogisticPoly2, gFitMin, gFitMax, 7);
    f.SetParameters(maxY * 1.2, 0.73, 0.02, 0.92, 0.08, -3.0, 15.0);
    f.SetParLimits(0, 0.0, 1e9);
    f.SetParLimits(1, 0.68, 0.82);
    f.SetParLimits(2, 0.002, 0.10);
    f.SetParLimits(3, 0.76, 1.08);
    f.SetParLimits(4, 0.005, 0.25);
    f.SetParLimits(5, -150.0, 150.0);
    f.SetParLimits(6, -600.0, 600.0);
    return f;
  }
  if (model == "Poly1Exp") {
    TF1 f(name.c_str(), Poly1Exp, gFitMin, gFitMax, 3);
    f.SetParameters(maxY, -3.0, -4.0);
    f.SetParLimits(0, 0.0, 1e9);
    f.SetParLimits(1, -80.0, 80.0);
    f.SetParLimits(2, -40.0, 20.0);
    return f;
  }
  if (model == "Poly2Exp") {
    TF1 f(name.c_str(), Poly2Exp, gFitMin, gFitMax, 4);
    f.SetParameters(maxY, -3.0, 10.0, -4.0);
    f.SetParLimits(0, 0.0, 1e9);
    f.SetParLimits(1, -120.0, 120.0);
    f.SetParLimits(2, -600.0, 600.0);
    f.SetParLimits(3, -40.0, 20.0);
    return f;
  }
  if (model == "SharpRiseLogisticExp") {
    TF1 f(name.c_str(), SharpRiseLogisticExp, gFitMin, gFitMax, 7);
    f.SetParameters(maxY * 1.4, 0.735, 0.02, 0.92, 0.08, 3.0, -4.0);
    f.SetParLimits(0, 0.0, 1e9);
    f.SetParLimits(1, 0.70, 0.80);
    f.SetParLimits(2, 0.002, 0.06);
    f.SetParLimits(3, 0.78, 1.08);
    f.SetParLimits(4, 0.005, 0.25);
    f.SetParLimits(5, 1.0, 20.0);
    f.SetParLimits(6, -40.0, 20.0);
    return f;
  }
  if (model == "SharpRiseLogisticPoly1") {
    TF1 f(name.c_str(), SharpRiseLogisticPoly1, gFitMin, gFitMax, 7);
    f.SetParameters(maxY * 1.4, 0.735, 0.02, 0.92, 0.08, 3.0, -3.0);
    f.SetParLimits(0, 0.0, 1e9);
    f.SetParLimits(1, 0.70, 0.80);
    f.SetParLimits(2, 0.002, 0.06);
    f.SetParLimits(3, 0.78, 1.08);
    f.SetParLimits(4, 0.005, 0.25);
    f.SetParLimits(5, 1.0, 20.0);
    f.SetParLimits(6, -80.0, 80.0);
    return f;
  }
  if (model == "DoubleLogisticExpPlusGauss") {
    TF1 f(name.c_str(), DoubleLogisticExpPlusGauss, gFitMin, gFitMax, 9);
    f.SetParameters(maxY * 1.2, 0.73, 0.02, 0.92, 0.08, -4.0,
                    0.15 * maxY, 0.75, 0.02);
    f.SetParLimits(0, 0.0, 1e9);
    f.SetParLimits(1, 0.68, 0.82);
    f.SetParLimits(2, 0.002, 0.10);
    f.SetParLimits(3, 0.76, 1.08);
    f.SetParLimits(4, 0.005, 0.25);
    f.SetParLimits(5, -40.0, 20.0);
    f.SetParLimits(6, 0.0, 1e9);
    f.SetParLimits(7, 0.70, 0.82);
    f.SetParLimits(8, 0.003, 0.08);
    return f;
  }
  if (model == "TwoDoubleLogisticExp") {
    TF1 f(name.c_str(), TwoDoubleLogisticExp, gFitMin, gFitMax, 12);
    double params[12] = {maxY, 0.73, 0.02, 0.92, 0.08, -4.0,
                         0.20 * maxY, 0.745, 0.01, 0.82, 0.05, -2.0};
    f.SetParameters(params);
    f.SetParLimits(0, 0.0, 1e9);
    f.SetParLimits(1, 0.68, 0.82);
    f.SetParLimits(2, 0.002, 0.10);
    f.SetParLimits(3, 0.76, 1.08);
    f.SetParLimits(4, 0.005, 0.25);
    f.SetParLimits(5, -40.0, 20.0);
    f.SetParLimits(6, 0.0, 1e9);
    f.SetParLimits(7, 0.70, 0.80);
    f.SetParLimits(8, 0.002, 0.05);
    f.SetParLimits(9, 0.75, 0.95);
    f.SetParLimits(10, 0.005, 0.20);
    f.SetParLimits(11, -40.0, 20.0);
    return f;
  }

  TF1 f(name.c_str(), DoubleLogisticExpPoly1, gFitMin, gFitMax, 7);
  f.SetParameters(maxY * 1.2, 0.73, 0.02, 0.92, 0.08, -3.0, -3.0);
  f.SetParLimits(0, 0.0, 1e9);
  f.SetParLimits(1, 0.68, 0.82);
  f.SetParLimits(2, 0.002, 0.10);
  f.SetParLimits(3, 0.76, 1.08);
  f.SetParLimits(4, 0.005, 0.25);
  f.SetParLimits(5, -40.0, 20.0);
  f.SetParLimits(6, -80.0, 80.0);
  return f;
}

FitResult runFit(TH1D* h, const std::string& category, const std::string& model,
                 std::vector<TF1>& drawnFunctions, int color, TH1D* hDisp, TPaveText& text) {
  TF1 f = makeModel("f_" + category + "_" + model, model, h);
  const int status = h->Fit(&f, "RQ0");

  TF1 draw(f);
  draw.SetName(("draw_" + category + "_" + model).c_str());
  draw.SetLineColor(color);
  draw.SetLineWidth(3);
  const double scale = hDisp->GetXaxis()->GetBinWidth(1) / h->GetXaxis()->GetBinWidth(1);
  draw.SetParameter(0, draw.GetParameter(0) * scale);
  drawnFunctions.push_back(draw);

  const double chi2ndf = (f.GetNDF() > 0) ? f.GetChisquare() / f.GetNDF() : 0.0;
  text.AddText(Form("%s: %.3f", model.c_str(), chi2ndf));

  FitResult result;
  result.category = category;
  result.model = model;
  result.chi2 = f.GetChisquare();
  result.ndf = f.GetNDF();
  result.chi2ndf = chi2ndf;
  result.converged = (status == 0);
  return result;
}

void drawCategory(TH1D* h, const std::string& category, const std::string& outputDir,
                  std::vector<FitResult>& results) {
  TH1D* hDisp = static_cast<TH1D*>(h->Clone((std::string(h->GetName()) + "_disp").c_str()));
  if (kDisplayRebin > 1) hDisp->Rebin(kDisplayRebin);

  TCanvas c(("c_" + category).c_str(), ("c_" + category).c_str(), 900, 800);
  c.Divide(1, 2);

  c.cd(1);
  gPad->SetPad(0.0, 0.32, 1.0, 1.0);
  gPad->SetLeftMargin(0.14);
  gPad->SetRightMargin(0.04);
  gPad->SetBottomMargin(0.03);
  hDisp->SetStats(0);
  hDisp->SetLineWidth(2);
  hDisp->SetTitle((category + " KShort treated as K#pi fit scan").c_str());
  hDisp->GetXaxis()->SetTitle("m(K#pi) [GeV]");
  hDisp->GetYaxis()->SetTitle("Candidates / bin");
  hDisp->Draw("E");

  TLegend leg(0.56, 0.61, 0.89, 0.88);
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
  leg.AddEntry(hDisp, "Histogram", "lep");

  TPaveText text(0.16, 0.13, 0.50, 0.34, "NDC");
  text.SetBorderSize(0);
  text.SetFillStyle(0);

  const std::vector<std::string> models = {
      "DoubleLogisticExp",
      "DoubleLogisticPoly1",
      "DoubleLogisticPoly2",
      "DoubleLogisticExpPoly1",
      "Poly1Exp",
      "Poly2Exp",
      "SharpRiseLogisticExp",
      "SharpRiseLogisticPoly1",
      "DoubleLogisticExpPlusGauss",
      "TwoDoubleLogisticExp"};
  const std::vector<int> colors = {kRed + 1, kBlue + 1, kGreen + 2, kMagenta + 1, kOrange + 7, kCyan + 2, kViolet + 1, kAzure + 2, kPink + 1, kSpring + 5};

  std::vector<TF1> drawnFunctions;
  drawnFunctions.reserve(models.size());

  for (size_t i = 0; i < models.size(); ++i) {
    FitResult result = runFit(h, category, models[i], drawnFunctions, colors[i], hDisp, text);
    results.push_back(result);
  }

  for (size_t i = 0; i < models.size(); ++i) {
    drawnFunctions[i].Draw("same");
    leg.AddEntry(&drawnFunctions[i], models[i].c_str(), "l");
  }
  leg.Draw();
  text.Draw();

  auto begin = results.end() - static_cast<long long>(models.size());
  auto bestIt = std::min_element(begin, results.end(),
                                 [](const FitResult& a, const FitResult& b) {
                                   return a.chi2ndf < b.chi2ndf;
                                 });

  c.cd(2);
  gPad->SetPad(0.0, 0.0, 1.0, 0.32);
  gPad->SetLeftMargin(0.14);
  gPad->SetRightMargin(0.04);
  gPad->SetTopMargin(0.03);
  gPad->SetBottomMargin(0.34);

  TF1 bestFit = makeModel("best_" + category, bestIt->model, h);
  h->Fit(&bestFit, "RQ0");

  TH1D hPull(("hPull_" + category).c_str(), "", hDisp->GetNbinsX(),
             hDisp->GetXaxis()->GetXmin(), hDisp->GetXaxis()->GetXmax());
  hPull.SetStats(0);
  hPull.SetMarkerStyle(20);
  hPull.SetMarkerSize(0.65);
  hPull.GetYaxis()->SetTitle("Pull");
  hPull.GetXaxis()->SetTitle("m(K#pi) [GeV]");
  hPull.GetYaxis()->SetNdivisions(505);
  hPull.GetYaxis()->SetTitleSize(0.11);
  hPull.GetYaxis()->SetLabelSize(0.085);
  hPull.GetYaxis()->SetTitleOffset(0.55);
  hPull.GetXaxis()->SetTitleSize(0.11);
  hPull.GetXaxis()->SetLabelSize(0.085);
  hPull.GetXaxis()->SetTitleOffset(1.05);
  hPull.SetMinimum(-5.0);
  hPull.SetMaximum(5.0);

  for (int i = 1; i <= hDisp->GetNbinsX(); ++i) {
    const double xMin = hDisp->GetBinLowEdge(i);
    const double xMax = xMin + hDisp->GetBinWidth(i);
    const double y = hDisp->GetBinContent(i);
    const double ey = hDisp->GetBinError(i);
    if (ey <= 0.0) continue;
    const double yfit = bestFit.Integral(xMin, xMax) / h->GetXaxis()->GetBinWidth(1);
    hPull.SetBinContent(i, (y - yfit) / ey);
    hPull.SetBinError(i, 1.0);
  }
  hPull.Draw("E");
  TLine line0(gFitMin, 0.0, gFitMax, 0.0);
  line0.SetLineStyle(2);
  line0.Draw("same");

  c.SaveAs((outputDir + "/kshort_wrong_as_kstar_" + category + "_fit_scan.pdf").c_str());
  delete hDisp;
}
}  // namespace

int main(int argc, char* argv[]) {
  const std::string inputFileName =
      getArgument(argc, argv, "--input", "KShortWrongAsKStarHistograms.root");
  const std::string outputDir =
      getArgument(argc, argv, "--output-dir", "KShortWrongAsKStarFitResults");
  gFitMin = getDoubleArgument(argc, argv, "--fit-min", gFitMin);
  gFitMax = getDoubleArgument(argc, argv, "--fit-max", gFitMax);

  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);
  gSystem->mkdir(outputDir.c_str(), true);

  TFile inputFile(inputFileName.c_str(), "READ");
  if (inputFile.IsZombie()) {
    std::cerr << "Cannot open " << inputFileName << std::endl;
    return 1;
  }

  TH1D* hNoPionTag = nullptr;
  TH1D* hWithPionTag = nullptr;
  inputFile.GetObject("hKShortWrongAsKStarMassKaonTag", hNoPionTag);
  inputFile.GetObject("hKShortWrongAsKStarMassKaonPionTag", hWithPionTag);
  if (hNoPionTag == nullptr || hWithPionTag == nullptr) {
    std::cerr << "Missing histograms in " << inputFileName << std::endl;
    return 1;
  }

  std::vector<FitResult> results;
  drawCategory(hNoPionTag, "no_pion_tag", outputDir, results);
  drawCategory(hWithPionTag, "with_pion_tag", outputDir, results);

  std::ofstream out(outputDir + "/kshort_wrong_as_kstar_fit_summary.csv");
  out << "category,model,chi2,ndf,chi2ndf,converged\n";
  for (const FitResult& r : results)
    out << r.category << "," << r.model << "," << r.chi2 << "," << r.ndf << ","
        << r.chi2ndf << "," << (r.converged ? "yes" : "no") << "\n";
  out.close();

  for (const std::string category : {"no_pion_tag", "with_pion_tag"}) {
    auto begin = std::find_if(results.begin(), results.end(),
                              [&](const FitResult& r) { return r.category == category; });
    auto end = std::find_if(results.rbegin(), results.rend(),
                            [&](const FitResult& r) { return r.category == category; }).base();
    auto bestIt = std::min_element(begin, end,
                                   [](const FitResult& a, const FitResult& b) {
                                     return a.chi2ndf < b.chi2ndf;
                                   });
    if (bestIt != end)
      std::cout << category << " best " << bestIt->model << " chi2/ndf=" << bestIt->chi2ndf << std::endl;
  }

  return 0;
}
