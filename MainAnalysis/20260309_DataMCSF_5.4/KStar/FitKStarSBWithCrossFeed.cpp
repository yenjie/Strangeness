#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>

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
constexpr double kKaonMass = 0.493677;
constexpr double kPionMass = 0.13957039;
constexpr double kThreshold = kKaonMass + kPionMass;
constexpr int kDisplayRebin = 4;
double gFitMin = 0.70;
double gFitMax = 1.10;

struct SignalShape {
  double mean = 0.8955;
  double sigma = 0.020;
  double alphaL = 1.5;
  double nL = 5.0;
  double alphaR = 1.5;
  double nR = 5.0;
  double gaussFrac = 0.2;
  double gaussSigma = 0.05;
  double gaussFrac2 = 0.0;
  double gaussSigma2 = 0.0;
  double gaussFrac3 = 0.0;
  double gaussSigma3 = 0.0;
};

struct PhiShape {
  double x0L = 0.735;
  double wL = 0.02;
  double x0R = 0.84;
  double wR = 0.05;
  double c1 = -3.0;
  double c2 = 10.0;
};

struct KShortShape {
  double x0L = 0.73;
  double wL = 0.02;
  double x0R = 0.92;
  double wR = 0.08;
  double slope = -4.0;
};

struct FitSummary {
  std::string category;
  double chi2 = 0.0;
  double ndf = 0.0;
  double chi2ndf = 0.0;
  double signalAmp = 0.0;
  double phiAmp = 0.0;
  double kshortAmp = 0.0;
  double thresholdNorm = 0.0;
  double thresholdPower = 0.0;
  double thresholdSlope = 0.0;
};

SignalShape gSignalShape;
PhiShape gPhiShape;
KShortShape gKShortShape;
std::string gSignalModel = "DoubleSidedCBPlusGauss";

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

double DoubleSidedCrystalBallUnit(double x, double mean, double sigma,
                                  double alphaL, double nL,
                                  double alphaR, double nR) {
  if (sigma <= 0.0 || alphaL <= 0.0 || nL <= 1.0 || alphaR <= 0.0 || nR <= 1.0) return 0.0;
  const double t = (x - mean) / sigma;
  if (t > -alphaL && t < alphaR) return std::exp(-0.5 * t * t);
  if (t <= -alphaL) {
    const double A = std::pow(nL / alphaL, nL) * std::exp(-0.5 * alphaL * alphaL);
    const double B = nL / alphaL - alphaL;
    return A / std::pow(B - t, nL);
  }
  const double A = std::pow(nR / alphaR, nR) * std::exp(-0.5 * alphaR * alphaR);
  const double B = nR / alphaR - alphaR;
  return A / std::pow(B + t, nR);
}

double SignalUnit(double x) {
  if (gSignalModel == "QuadGaussian") {
    const double g1 = std::exp(-0.5 * std::pow((x - gSignalShape.mean) / gSignalShape.sigma, 2));
    const double g2 = std::exp(-0.5 * std::pow((x - gSignalShape.mean) / gSignalShape.gaussSigma, 2));
    const double g3 = std::exp(-0.5 * std::pow((x - gSignalShape.mean) / gSignalShape.gaussSigma2, 2));
    const double g4 = std::exp(-0.5 * std::pow((x - gSignalShape.mean) / gSignalShape.gaussSigma3, 2));
    return g1 + gSignalShape.gaussFrac * g2
        + gSignalShape.gaussFrac2 * g3 + gSignalShape.gaussFrac3 * g4;
  }
  const double dscb =
      DoubleSidedCrystalBallUnit(x, gSignalShape.mean, gSignalShape.sigma,
                                 gSignalShape.alphaL, gSignalShape.nL,
                                 gSignalShape.alphaR, gSignalShape.nR);
  if (gSignalModel == "DoubleSidedCB") return dscb;
  const double gauss = std::exp(-0.5 * std::pow((x - gSignalShape.mean) / gSignalShape.gaussSigma, 2));
  return dscb + gSignalShape.gaussFrac * gauss;
}

double PhiUnit(double x) {
  const double dx = x - 0.80;
  const double env = logisticRise(x, gPhiShape.x0L, gPhiShape.wL)
      * logisticFall(x, gPhiShape.x0R, gPhiShape.wR);
  const double poly = std::max(0.0, 1.0 + gPhiShape.c1 * dx + gPhiShape.c2 * dx * dx);
  return env * poly;
}

double KShortUnit(double x) {
  const double dx = x - 0.82;
  const double env = logisticRise(x, gKShortShape.x0L, gKShortShape.wL)
      * logisticFall(x, gKShortShape.x0R, gKShortShape.wR);
  return env * std::exp(gKShortShape.slope * dx);
}

double ThresholdExp1(double x, double norm, double power, double slope) {
  if (x <= kThreshold) return 0.0;
  return norm * std::pow(x - kThreshold, power) * std::exp(slope * x);
}

double TotalKaonTag(double* x, double* p) {
  return p[0] * SignalUnit(x[0]) + p[1] * PhiUnit(x[0]) + p[2] * KShortUnit(x[0])
      + ThresholdExp1(x[0], p[3], p[4], p[5]);
}

double TotalKaonPionTag(double* x, double* p) {
  return p[0] * SignalUnit(x[0]) + p[1] * KShortUnit(x[0])
      + ThresholdExp1(x[0], p[2], p[3], p[4]);
}

SignalShape deriveSignalShape(TH1D* h) {
  if (gSignalModel == "QuadGaussian") {
    TF1 f("fSignalShapeKStarQuadG",
          "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]*exp(-0.5*((x-[1])/[4])^2)+"
          "[5]*exp(-0.5*((x-[1])/[6])^2)+[7]*exp(-0.5*((x-[1])/[8])^2)",
          gFitMin, gFitMax);
    f.SetParameters(std::max(100.0, 0.40 * h->GetMaximum()), 0.8955, 0.012,
                    std::max(50.0, 0.28 * h->GetMaximum()), 0.025,
                    std::max(30.0, 0.20 * h->GetMaximum()), 0.045,
                    std::max(10.0, 0.12 * h->GetMaximum()), 0.080);
    f.SetParLimits(1, 0.86, 0.93);
    f.SetParLimits(2, 0.002, 0.08);
    f.SetParLimits(4, 0.002, 0.12);
    f.SetParLimits(6, 0.002, 0.18);
    f.SetParLimits(8, 0.002, 0.24);
    h->Fit(&f, "RQ0");

    SignalShape shape;
    shape.mean = f.GetParameter(1);
    shape.sigma = f.GetParameter(2);
    shape.gaussSigma = f.GetParameter(4);
    shape.gaussSigma2 = f.GetParameter(6);
    shape.gaussSigma3 = f.GetParameter(8);
    shape.gaussFrac = (f.GetParameter(0) != 0.0) ? f.GetParameter(3) / f.GetParameter(0) : 0.0;
    shape.gaussFrac2 = (f.GetParameter(0) != 0.0) ? f.GetParameter(5) / f.GetParameter(0) : 0.0;
    shape.gaussFrac3 = (f.GetParameter(0) != 0.0) ? f.GetParameter(7) / f.GetParameter(0) : 0.0;
    shape.alphaL = shape.alphaR = 1.0;
    shape.nL = shape.nR = 2.0;
    return shape;
  }

  if (gSignalModel == "DoubleSidedCB") {
    TF1 f("fSignalShapeKStarDSCB", [](double* x, double* p) {
      return p[0] * DoubleSidedCrystalBallUnit(x[0], p[1], p[2], p[3], p[4], p[5], p[6]);
    }, gFitMin, gFitMax, 7);
    f.SetParameters(std::max(100.0, 0.8 * h->GetMaximum()), 0.8955, 0.020, 1.5, 5.0, 1.5, 5.0);
    f.SetParLimits(1, 0.86, 0.93);
    f.SetParLimits(2, 0.002, 0.10);
    f.SetParLimits(3, 0.2, 8.0);
    f.SetParLimits(4, 1.2, 80.0);
    f.SetParLimits(5, 0.2, 8.0);
    f.SetParLimits(6, 1.2, 80.0);
    h->Fit(&f, "RQ0");

    SignalShape shape;
    shape.mean = f.GetParameter(1);
    shape.sigma = f.GetParameter(2);
    shape.alphaL = f.GetParameter(3);
    shape.nL = f.GetParameter(4);
    shape.alphaR = f.GetParameter(5);
    shape.nR = f.GetParameter(6);
    shape.gaussFrac = 0.0;
    shape.gaussSigma = shape.sigma;
    return shape;
  }

  TF1 f("fSignalShapeKStar", [](double* x, double* p) {
    const double xx = x[0];
    const double dscb = p[0] * DoubleSidedCrystalBallUnit(xx, p[1], p[2], p[3], p[4], p[5], p[6]);
    const double gauss = p[7] * std::exp(-0.5 * std::pow((xx - p[1]) / p[8], 2));
    return dscb + gauss;
  }, gFitMin, gFitMax, 9);
  f.SetParameters(std::max(100.0, 0.8 * h->GetMaximum()), 0.8955, 0.020, 1.5, 5.0, 1.5, 5.0,
                  std::max(30.0, 0.2 * h->GetMaximum()), 0.05);
  f.SetParLimits(1, 0.86, 0.93);
  f.SetParLimits(2, 0.002, 0.10);
  f.SetParLimits(3, 0.2, 8.0);
  f.SetParLimits(4, 1.2, 80.0);
  f.SetParLimits(5, 0.2, 8.0);
  f.SetParLimits(6, 1.2, 80.0);
  f.SetParLimits(8, 0.002, 0.20);
  h->Fit(&f, "RQ0");

  SignalShape shape;
  shape.mean = f.GetParameter(1);
  shape.sigma = f.GetParameter(2);
  shape.alphaL = f.GetParameter(3);
  shape.nL = f.GetParameter(4);
  shape.alphaR = f.GetParameter(5);
  shape.nR = f.GetParameter(6);
  shape.gaussFrac = (f.GetParameter(0) != 0.0) ? f.GetParameter(7) / f.GetParameter(0) : 0.0;
  shape.gaussSigma = f.GetParameter(8);
  return shape;
}

PhiShape derivePhiShape(TH1D* h) {
  TF1 f("fPhiShape", [](double* x, double* p) {
    const double xx = x[0];
    const double dx = xx - 0.80;
    const double env = p[0] * logisticRise(xx, p[1], p[2]) * logisticFall(xx, p[3], p[4]);
    const double poly = std::max(0.0, 1.0 + p[5] * dx + p[6] * dx * dx);
    return env * poly;
  }, gFitMin, gFitMax, 7);
  f.SetParameters(std::max(50.0, 1.2 * h->GetMaximum()), 0.73, 0.02, 0.84, 0.05, -3.0, 10.0);
  f.SetParLimits(0, 0.0, 1e9);
  f.SetParLimits(1, 0.68, 0.80);
  f.SetParLimits(2, 0.002, 0.08);
  f.SetParLimits(3, 0.76, 1.05);
  f.SetParLimits(4, 0.005, 0.20);
  f.SetParLimits(5, -150.0, 150.0);
  f.SetParLimits(6, -600.0, 600.0);
  h->Fit(&f, "RQ0");

  PhiShape shape;
  shape.x0L = f.GetParameter(1);
  shape.wL = f.GetParameter(2);
  shape.x0R = f.GetParameter(3);
  shape.wR = f.GetParameter(4);
  shape.c1 = f.GetParameter(5);
  shape.c2 = f.GetParameter(6);
  return shape;
}

KShortShape deriveKShortShape(TH1D* h) {
  TF1 f("fKShortShape", [](double* x, double* p) {
    const double xx = x[0];
    const double dx = xx - 0.82;
    const double env = p[0] * logisticRise(xx, p[1], p[2]) * logisticFall(xx, p[3], p[4]);
    return env * std::exp(p[5] * dx);
  }, gFitMin, gFitMax, 6);
  f.SetParameters(std::max(50.0, 1.2 * h->GetMaximum()), 0.73, 0.02, 0.92, 0.08, -4.0);
  f.SetParLimits(0, 0.0, 1e9);
  f.SetParLimits(1, 0.68, 0.82);
  f.SetParLimits(2, 0.002, 0.10);
  f.SetParLimits(3, 0.76, 1.08);
  f.SetParLimits(4, 0.005, 0.25);
  f.SetParLimits(5, -40.0, 20.0);
  h->Fit(&f, "RQ0");

  KShortShape shape;
  shape.x0L = f.GetParameter(1);
  shape.wL = f.GetParameter(2);
  shape.x0R = f.GetParameter(3);
  shape.wR = f.GetParameter(4);
  shape.slope = f.GetParameter(5);
  return shape;
}

FitSummary fitKaonTag(TH1D* hSignal, TH1D* hSB, TH1D* hPhi, TH1D* hKShort,
                      const std::string& outputDir) {
  gSignalShape = deriveSignalShape(hSignal);
  gPhiShape = derivePhiShape(hPhi);
  gKShortShape = deriveKShortShape(hKShort);

  TF1 total("fTotalKaonTagCrossFeed", TotalKaonTag, gFitMin, gFitMax, 6);
  total.SetParNames("S", "Phi", "KShort", "N", "p", "b1");
  const double maxY = std::max(1000.0, hSB->GetMaximum());
  total.SetParameters(0.06 * maxY, 0.05 * maxY, 0.03 * maxY, 0.3 * maxY, 0.8, -4.0);
  total.SetParLimits(0, 0.0, 1e9);
  total.SetParLimits(1, 0.0, 1e9);
  total.SetParLimits(2, 0.0, 1e9);
  total.SetParLimits(3, 0.0, 1e12);
  total.SetParLimits(4, 0.0, 10.0);
  total.SetParLimits(5, -30.0, 10.0);
  hSB->Fit(&total, "RQ0");

  TH1D* hDisp = static_cast<TH1D*>(hSB->Clone("hKaonTagCrossFeedDisp"));
  if (kDisplayRebin > 1) hDisp->Rebin(kDisplayRebin);
  const double scale = hDisp->GetXaxis()->GetBinWidth(1) / hSB->GetXaxis()->GetBinWidth(1);

  TF1 signalDraw("fSignalDrawKaonTag", [](double* x, double* p) { return p[0] * SignalUnit(x[0]); }, gFitMin, gFitMax, 1);
  signalDraw.SetParameter(0, total.GetParameter(0) * scale);
  signalDraw.SetLineColor(kBlue + 1);
  signalDraw.SetLineWidth(3);

  TF1 phiDraw("fPhiDrawKaonTag", [](double* x, double* p) { return p[0] * PhiUnit(x[0]); }, gFitMin, gFitMax, 1);
  phiDraw.SetParameter(0, total.GetParameter(1) * scale);
  phiDraw.SetLineColor(kMagenta + 1);
  phiDraw.SetLineWidth(3);
  phiDraw.SetLineStyle(2);

  TF1 kshortDraw("fKShortDrawKaonTag", [](double* x, double* p) { return p[0] * KShortUnit(x[0]); }, gFitMin, gFitMax, 1);
  kshortDraw.SetParameter(0, total.GetParameter(2) * scale);
  kshortDraw.SetLineColor(kOrange + 7);
  kshortDraw.SetLineWidth(3);
  kshortDraw.SetLineStyle(3);

  TF1 thresholdDraw("fThresholdDrawKaonTag", [](double* x, double* p) {
    return ThresholdExp1(x[0], p[0], p[1], p[2]);
  }, gFitMin, gFitMax, 3);
  thresholdDraw.SetParameters(total.GetParameter(3) * scale, total.GetParameter(4), total.GetParameter(5));
  thresholdDraw.SetLineColor(kGreen + 2);
  thresholdDraw.SetLineWidth(3);
  thresholdDraw.SetLineStyle(4);

  TF1 backgroundDraw("fBackgroundDrawKaonTag", [](double* x, double* p) {
    return p[0] * PhiUnit(x[0]) + p[1] * KShortUnit(x[0])
        + ThresholdExp1(x[0], p[2], p[3], p[4]);
  }, gFitMin, gFitMax, 5);
  backgroundDraw.SetParameters(total.GetParameter(1) * scale, total.GetParameter(2) * scale,
                               total.GetParameter(3) * scale, total.GetParameter(4),
                               total.GetParameter(5));
  backgroundDraw.SetLineColor(kGray + 2);
  backgroundDraw.SetLineWidth(3);
  backgroundDraw.SetLineStyle(7);

  TF1 totalDraw(total);
  totalDraw.SetParameter(0, total.GetParameter(0) * scale);
  totalDraw.SetParameter(1, total.GetParameter(1) * scale);
  totalDraw.SetParameter(2, total.GetParameter(2) * scale);
  totalDraw.SetParameter(3, total.GetParameter(3) * scale);
  totalDraw.SetLineColor(kRed + 1);
  totalDraw.SetLineWidth(4);

  TCanvas c("cKaonTagCrossFeed", "cKaonTagCrossFeed", 900, 800);
  c.Divide(1, 2);
  c.cd(1);
  gPad->SetPad(0.0, 0.32, 1.0, 1.0);
  gPad->SetLeftMargin(0.14);
  gPad->SetRightMargin(0.04);
  gPad->SetBottomMargin(0.03);
  hDisp->SetStats(0);
  hDisp->SetLineWidth(2);
  hDisp->SetTitle("kaon-tag reco-only S+B fit with #phi/K^{0}_{S} cross-feed");
  hDisp->GetXaxis()->SetTitle("m(K#pi) [GeV]");
  hDisp->GetYaxis()->SetTitle("Assignments / bin");
  hDisp->Draw("E");
  totalDraw.Draw("same");
  backgroundDraw.Draw("same");
  signalDraw.Draw("same");
  phiDraw.Draw("same");
  kshortDraw.Draw("same");
  thresholdDraw.Draw("same");

  TLegend leg(0.45, 0.54, 0.89, 0.88);
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
  leg.AddEntry(hDisp, "Reco-only same-event assignments", "lep");
  leg.AddEntry(&totalDraw, "Total fit", "l");
  leg.AddEntry(&backgroundDraw, "Total background", "l");
  leg.AddEntry(&signalDraw, "K^{*} signal", "l");
  leg.AddEntry(&phiDraw, "#phi wrong-treatment shape", "l");
  leg.AddEntry(&kshortDraw, "K^{0}_{S} wrong-treatment shape", "l");
  leg.AddEntry(&thresholdDraw, "Threshold background", "l");
  leg.Draw();

  TPaveText txt(0.14, 0.12, 0.44, 0.30, "NDC");
  txt.SetBorderSize(0);
  txt.SetFillStyle(0);
  txt.AddText(Form("#chi^{2}/ndf = %.3f", total.GetChisquare() / total.GetNDF()));
  txt.AddText(Form("S = %.1f, #phi = %.1f", total.GetParameter(0), total.GetParameter(1)));
  txt.AddText(Form("K^{0}_{S} = %.1f", total.GetParameter(2)));
  txt.Draw();

  c.cd(2);
  gPad->SetPad(0.0, 0.0, 1.0, 0.32);
  gPad->SetLeftMargin(0.14);
  gPad->SetRightMargin(0.04);
  gPad->SetTopMargin(0.03);
  gPad->SetBottomMargin(0.34);
  TH1D hPull("hPullKaonTagCrossFeed", "", hDisp->GetNbinsX(), hDisp->GetXaxis()->GetXmin(), hDisp->GetXaxis()->GetXmax());
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
    const double yfit = total.Integral(xMin, xMax) / hSB->GetXaxis()->GetBinWidth(1);
    hPull.SetBinContent(i, (y - yfit) / ey);
    hPull.SetBinError(i, 1.0);
  }
  hPull.Draw("E");
  TLine line0(gFitMin, 0.0, gFitMax, 0.0);
  line0.SetLineStyle(2);
  line0.Draw("same");

  c.SaveAs((outputDir + "/kstar_sb_kaon_tag_crossfeed_fit.pdf").c_str());
  delete hDisp;

  FitSummary s;
  s.category = "kaon_tag";
  s.chi2 = total.GetChisquare();
  s.ndf = total.GetNDF();
  s.chi2ndf = total.GetChisquare() / total.GetNDF();
  s.signalAmp = total.GetParameter(0);
  s.phiAmp = total.GetParameter(1);
  s.kshortAmp = total.GetParameter(2);
  s.thresholdNorm = total.GetParameter(3);
  s.thresholdPower = total.GetParameter(4);
  s.thresholdSlope = total.GetParameter(5);
  return s;
}

FitSummary fitKaonPionTag(TH1D* hSignal, TH1D* hSB, TH1D* hKShort,
                          const std::string& outputDir) {
  gSignalShape = deriveSignalShape(hSignal);
  gKShortShape = deriveKShortShape(hKShort);

  TF1 total("fTotalKaonPionTagCrossFeed", TotalKaonPionTag, gFitMin, gFitMax, 5);
  total.SetParNames("S", "KShort", "N", "p", "b1");
  const double maxY = std::max(1000.0, hSB->GetMaximum());
  total.SetParameters(0.08 * maxY, 0.03 * maxY, 0.3 * maxY, 0.8, -4.0);
  total.SetParLimits(0, 0.0, 1e9);
  total.SetParLimits(1, 0.0, 1e9);
  total.SetParLimits(2, 0.0, 1e12);
  total.SetParLimits(3, 0.0, 10.0);
  total.SetParLimits(4, -30.0, 10.0);
  hSB->Fit(&total, "RQ0");

  TH1D* hDisp = static_cast<TH1D*>(hSB->Clone("hKaonPionTagCrossFeedDisp"));
  if (kDisplayRebin > 1) hDisp->Rebin(kDisplayRebin);
  const double scale = hDisp->GetXaxis()->GetBinWidth(1) / hSB->GetXaxis()->GetBinWidth(1);

  TF1 signalDraw("fSignalDrawKaonPionTag", [](double* x, double* p) { return p[0] * SignalUnit(x[0]); }, gFitMin, gFitMax, 1);
  signalDraw.SetParameter(0, total.GetParameter(0) * scale);
  signalDraw.SetLineColor(kBlue + 1);
  signalDraw.SetLineWidth(3);

  TF1 kshortDraw("fKShortDrawKaonPionTag", [](double* x, double* p) { return p[0] * KShortUnit(x[0]); }, gFitMin, gFitMax, 1);
  kshortDraw.SetParameter(0, total.GetParameter(1) * scale);
  kshortDraw.SetLineColor(kOrange + 7);
  kshortDraw.SetLineWidth(3);
  kshortDraw.SetLineStyle(3);

  TF1 thresholdDraw("fThresholdDrawKaonPionTag", [](double* x, double* p) {
    return ThresholdExp1(x[0], p[0], p[1], p[2]);
  }, gFitMin, gFitMax, 3);
  thresholdDraw.SetParameters(total.GetParameter(2) * scale, total.GetParameter(3), total.GetParameter(4));
  thresholdDraw.SetLineColor(kGreen + 2);
  thresholdDraw.SetLineWidth(3);
  thresholdDraw.SetLineStyle(4);

  TF1 backgroundDraw("fBackgroundDrawKaonPionTag", [](double* x, double* p) {
    return p[0] * KShortUnit(x[0]) + ThresholdExp1(x[0], p[1], p[2], p[3]);
  }, gFitMin, gFitMax, 4);
  backgroundDraw.SetParameters(total.GetParameter(1) * scale, total.GetParameter(2) * scale,
                               total.GetParameter(3), total.GetParameter(4));
  backgroundDraw.SetLineColor(kGray + 2);
  backgroundDraw.SetLineWidth(3);
  backgroundDraw.SetLineStyle(7);

  TF1 totalDraw(total);
  totalDraw.SetParameter(0, total.GetParameter(0) * scale);
  totalDraw.SetParameter(1, total.GetParameter(1) * scale);
  totalDraw.SetParameter(2, total.GetParameter(2) * scale);
  totalDraw.SetLineColor(kRed + 1);
  totalDraw.SetLineWidth(4);

  TCanvas c("cKaonPionTagCrossFeed", "cKaonPionTagCrossFeed", 900, 800);
  c.Divide(1, 2);
  c.cd(1);
  gPad->SetPad(0.0, 0.32, 1.0, 1.0);
  gPad->SetLeftMargin(0.14);
  gPad->SetRightMargin(0.04);
  gPad->SetBottomMargin(0.03);
  hDisp->SetStats(0);
  hDisp->SetLineWidth(2);
  hDisp->SetTitle("kaon+pion-tag reco-only S+B fit with K^{0}_{S} cross-feed");
  hDisp->GetXaxis()->SetTitle("m(K#pi) [GeV]");
  hDisp->GetYaxis()->SetTitle("Assignments / bin");
  hDisp->Draw("E");
  totalDraw.Draw("same");
  backgroundDraw.Draw("same");
  signalDraw.Draw("same");
  kshortDraw.Draw("same");
  thresholdDraw.Draw("same");

  TLegend leg(0.45, 0.57, 0.89, 0.88);
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
  leg.AddEntry(hDisp, "Reco-only same-event assignments", "lep");
  leg.AddEntry(&totalDraw, "Total fit", "l");
  leg.AddEntry(&backgroundDraw, "Total background", "l");
  leg.AddEntry(&signalDraw, "K^{*} signal", "l");
  leg.AddEntry(&kshortDraw, "K^{0}_{S} wrong-treatment shape", "l");
  leg.AddEntry(&thresholdDraw, "Threshold background", "l");
  leg.Draw();

  TPaveText txt(0.14, 0.12, 0.44, 0.28, "NDC");
  txt.SetBorderSize(0);
  txt.SetFillStyle(0);
  txt.AddText(Form("#chi^{2}/ndf = %.3f", total.GetChisquare() / total.GetNDF()));
  txt.AddText(Form("S = %.1f, K^{0}_{S} = %.1f", total.GetParameter(0), total.GetParameter(1)));
  txt.Draw();

  c.cd(2);
  gPad->SetPad(0.0, 0.0, 1.0, 0.32);
  gPad->SetLeftMargin(0.14);
  gPad->SetRightMargin(0.04);
  gPad->SetTopMargin(0.03);
  gPad->SetBottomMargin(0.34);
  TH1D hPull("hPullKaonPionTagCrossFeed", "", hDisp->GetNbinsX(), hDisp->GetXaxis()->GetXmin(), hDisp->GetXaxis()->GetXmax());
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
    const double yfit = total.Integral(xMin, xMax) / hSB->GetXaxis()->GetBinWidth(1);
    hPull.SetBinContent(i, (y - yfit) / ey);
    hPull.SetBinError(i, 1.0);
  }
  hPull.Draw("E");
  TLine line0(gFitMin, 0.0, gFitMax, 0.0);
  line0.SetLineStyle(2);
  line0.Draw("same");

  c.SaveAs((outputDir + "/kstar_sb_kaon_pion_tag_crossfeed_fit.pdf").c_str());
  delete hDisp;

  FitSummary s;
  s.category = "kaon_pion_tag";
  s.chi2 = total.GetChisquare();
  s.ndf = total.GetNDF();
  s.chi2ndf = total.GetChisquare() / total.GetNDF();
  s.signalAmp = total.GetParameter(0);
  s.kshortAmp = total.GetParameter(1);
  s.thresholdNorm = total.GetParameter(2);
  s.thresholdPower = total.GetParameter(3);
  s.thresholdSlope = total.GetParameter(4);
  return s;
}
}  // namespace

int main(int argc, char* argv[]) {
  const std::string signalInput =
      getArgument(argc, argv, "--signal-input", "KStarCombinedAssignmentHistograms.root");
  const std::string sbInput =
      getArgument(argc, argv, "--sb-input", "KStarSBHistograms.root");
  const std::string phiInput =
      getArgument(argc, argv, "--phi-input", "PhiWrongAsKStarHistograms.root");
  const std::string kshortInput =
      getArgument(argc, argv, "--kshort-input", "KShortWrongAsKStarHistograms.root");
  const std::string outputDir =
      getArgument(argc, argv, "--output-dir", "SBFitResultsCrossFeed");
  gSignalModel = getArgument(argc, argv, "--signal-model", gSignalModel);
  gFitMin = getDoubleArgument(argc, argv, "--fit-min", gFitMin);
  gFitMax = getDoubleArgument(argc, argv, "--fit-max", gFitMax);

  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);
  gSystem->mkdir(outputDir.c_str(), true);

  TFile signalFile(signalInput.c_str(), "READ");
  TFile sbFile(sbInput.c_str(), "READ");
  TFile phiFile(phiInput.c_str(), "READ");
  TFile kshortFile(kshortInput.c_str(), "READ");
  if (signalFile.IsZombie() || sbFile.IsZombie() || phiFile.IsZombie() || kshortFile.IsZombie()) {
    std::cerr << "Error opening input files" << std::endl;
    return 1;
  }

  TH1D *hSignal1 = nullptr, *hSignal2 = nullptr, *hSB1 = nullptr, *hSB2 = nullptr;
  TH1D *hPhi1 = nullptr, *hPhi2 = nullptr, *hKShort1 = nullptr, *hKShort2 = nullptr;
  signalFile.GetObject("hKStarMassKaonTag", hSignal1);
  signalFile.GetObject("hKStarMassKaonPionTag", hSignal2);
  sbFile.GetObject("hKStarSBMassKaonTag", hSB1);
  sbFile.GetObject("hKStarSBMassKaonPionTag", hSB2);
  phiFile.GetObject("hPhiWrongAsKStarMassKaonTag", hPhi1);
  phiFile.GetObject("hPhiWrongAsKStarMassKaonPionTag", hPhi2);
  kshortFile.GetObject("hKShortWrongAsKStarMassKaonTag", hKShort1);
  kshortFile.GetObject("hKShortWrongAsKStarMassKaonPionTag", hKShort2);

  if (hSignal1 == nullptr || hSignal2 == nullptr || hSB1 == nullptr || hSB2 == nullptr ||
      hPhi1 == nullptr || hPhi2 == nullptr || hKShort1 == nullptr || hKShort2 == nullptr) {
    std::cerr << "Missing required histograms" << std::endl;
    return 1;
  }

  FitSummary kaonTag = fitKaonTag(hSignal1, hSB1, hPhi1, hKShort1, outputDir);
  FitSummary kaonPionTag = fitKaonPionTag(hSignal2, hSB2, hKShort2, outputDir);

  std::ofstream out(outputDir + "/kstar_sb_crossfeed_summary.csv");
  out << "category,chi2,ndf,chi2ndf,signalAmp,phiAmp,kshortAmp,thresholdNorm,thresholdPower,thresholdSlope\n";
  for (const FitSummary& s : {kaonTag, kaonPionTag}) {
    out << s.category << "," << s.chi2 << "," << s.ndf << "," << s.chi2ndf << ","
        << s.signalAmp << "," << s.phiAmp << "," << s.kshortAmp << ","
        << s.thresholdNorm << "," << s.thresholdPower << "," << s.thresholdSlope << "\n";
  }
  out.close();

  std::cout << "kaon_tag cross-feed chi2/ndf = " << kaonTag.chi2ndf << std::endl;
  std::cout << "kaon_pion_tag cross-feed chi2/ndf = " << kaonPionTag.chi2ndf << std::endl;
  return 0;
}
