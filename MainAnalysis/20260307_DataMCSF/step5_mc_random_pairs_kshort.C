#include <TCanvas.h>
#include <TFile.h>
#include <TF1.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TLine.h>
#include <TPaveText.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace {
constexpr double kPionMass = 0.13957039;
constexpr double kAbsCosMin = 0.15;
constexpr double kAbsCosMax = 0.675;
constexpr int kMaxReco = 256;
constexpr int kDisplayRebin = 4;
double gMassMin = 0.40;
double gMassMax = 0.70;
double gFitMin = 0.40;
double gFitMax = 0.70;
int gBins = 280;

struct Step4Seed {
  std::string model = "TripleGauss";
  double mean = 0.4976;
  double sigma1 = 0.006;
  double sigma2 = 0.015;
  double sigma3 = 0.030;
  double gamma = 0.004;
  double frac1 = 0.5;
  double frac2 = 0.3;
  double frac3 = 0.2;
};

struct FitResult {
  std::string model;
  double chi2ndf = 1e9;
  double mean = 0.0;
  double meanErr = 0.0;
  double yield = 0.0;
  double yieldErr = 0.0;
  double pullRms = 0.0;
  double pullMeanAbs = 0.0;
  double pullMaxAbs = 0.0;
  TF1* total = nullptr;
  TF1* signal = nullptr;
  TF1* background = nullptr;
};

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
  return m2 > 0.0 ? std::sqrt(m2) : 0.0;
}

bool passTrackAcceptance(double px, double py, double pz) {
  const double p = std::sqrt(px * px + py * py + pz * pz);
  if (p <= 0.0) return false;
  const double abscos = std::fabs(pz / p);
  return (abscos >= kAbsCosMin && abscos <= kAbsCosMax);
}

std::vector<std::string> splitCsv(const std::string& s) {
  std::vector<std::string> out;
  std::stringstream ss(s);
  std::string item;
  while (std::getline(ss, item, ',')) out.push_back(item);
  return out;
}

std::map<std::string, Step4Seed> loadStep4Seeds(const std::string& csvPath) {
  std::map<std::string, Step4Seed> seeds;
  std::ifstream in(csvPath);
  if (!in) {
    std::cerr << "Warning: cannot open " << csvPath
              << ", falling back to defaults for all categories." << std::endl;
    return seeds;
  }

  std::string line;
  if (!std::getline(in, line)) return seeds;
  while (std::getline(in, line)) {
    if (line.empty()) continue;
    auto cols = splitCsv(line);
    if (cols.size() < 16) continue;

    Step4Seed s;
    const std::string category = cols[0];
    s.model = cols[1];
    s.mean = std::atof(cols[3].c_str());
    s.sigma1 = std::atof(cols[5].c_str());
    s.sigma2 = std::atof(cols[7].c_str());
    s.sigma3 = std::atof(cols[9].c_str());
    s.gamma = std::atof(cols[11].c_str());
    if (cols.size() > 18) {
      s.frac1 = std::atof(cols[16].c_str());
      s.frac2 = std::atof(cols[17].c_str());
      s.frac3 = std::atof(cols[18].c_str());
    }
    const double fsum = s.frac1 + s.frac2 + s.frac3;
    if (fsum > 0.0) {
      s.frac1 /= fsum;
      s.frac2 /= fsum;
      s.frac3 /= fsum;
    } else {
      s.frac1 = 0.5;
      s.frac2 = 0.3;
      s.frac3 = 0.2;
    }
    seeds[category] = s;
  }
  return seeds;
}

FitResult fitCategory(TH1D* h, const Step4Seed& seed, const std::string& cat, const std::string& modelPrefix,
                      const std::string& signalModelOverride, const std::string& bkgMode) {
  FitResult r;
  if (!signalModelOverride.empty() && signalModelOverride != "auto") {
    r.model = signalModelOverride;
  } else {
    r.model = (seed.model.empty() ? std::string("TripleGauss") : seed.model);
  }
  bool usePol3Bkg = (cat == "2tag");
  const bool usePol4Bkg = (bkgMode == "all_pol4");
  const bool usePol3PlainBkg = (bkgMode == "all_pol3_plain");
  const bool useExpoPol2Bkg = (bkgMode == "all_expopol2");
  if (bkgMode == "all_pol2") usePol3Bkg = false;
  if (bkgMode == "all_pol3") usePol3Bkg = true;
  if (bkgMode == "baseline") usePol3Bkg = (cat == "2tag");

  const double maxY = std::max(10.0, h->GetMaximum());
  const double b0 = 0.5 * (h->GetBinContent(1) + h->GetBinContent(h->GetNbinsX()));
  const double meanLo = std::max(0.485, seed.mean - 0.008);
  const double meanHi = std::min(0.510, seed.mean + 0.008);

  if (r.model == "Gauss") {
    if (usePol3Bkg) {
      r.total = new TF1((modelPrefix + "_tot").c_str(), "gaus(0)+pol3(3)+gaus(7)", gFitMin, gFitMax);
      r.total->SetParameters(maxY, seed.mean, std::max(0.001, seed.sigma1), b0, 0.0, 0.0, 0.0, 0.10 * maxY, 0.400, 0.010);
      r.total->SetParLimits(7, 0.0, 1e9);
      r.total->SetParLimits(8, 0.390, 0.410);
      r.total->SetParLimits(9, 0.001, 0.05);
    } else {
      r.total = new TF1((modelPrefix + "_tot").c_str(), "gaus(0)+pol2(3)", gFitMin, gFitMax);
      r.total->SetParameters(maxY, seed.mean, std::max(0.001, seed.sigma1), b0, 0.0, 0.0);
    }
    r.total->SetParLimits(0, 0.0, 1e9);
    r.total->SetParLimits(1, meanLo, meanHi);
    r.total->SetParLimits(2, 0.001, 0.05);

    h->Fit(r.total, "RQ0");

    r.signal = new TF1((modelPrefix + "_sig").c_str(), "gaus(0)", gFitMin, gFitMax);
    r.signal->SetParameters(r.total->GetParameter(0), r.total->GetParameter(1), r.total->GetParameter(2));

    if (usePol3Bkg) {
      r.background = new TF1((modelPrefix + "_bkg").c_str(), "pol3(0)+gaus(4)", gFitMin, gFitMax);
      r.background->SetParameters(r.total->GetParameter(3), r.total->GetParameter(4), r.total->GetParameter(5), r.total->GetParameter(6),
                                  r.total->GetParameter(7), r.total->GetParameter(8), r.total->GetParameter(9));
    } else {
      r.background = new TF1((modelPrefix + "_bkg").c_str(), "pol2(0)", gFitMin, gFitMax);
      r.background->SetParameters(r.total->GetParameter(3), r.total->GetParameter(4), r.total->GetParameter(5));
    }
  } else if (r.model == "Voigt") {
    if (usePol3Bkg) {
      r.total = new TF1((modelPrefix + "_tot").c_str(), "[0]*TMath::Voigt(x-[1],[2],[3],4)+pol3(4)+gaus(8)", gFitMin, gFitMax);
      r.total->SetParameters(maxY, seed.mean, std::max(0.0005, seed.sigma1), std::max(0.001, seed.gamma), b0, 0.0, 0.0, 0.0,
                             0.10 * maxY, 0.400, 0.010);
      r.total->SetParLimits(8, 0.0, 1e9);
      r.total->SetParLimits(9, 0.390, 0.410);
      r.total->SetParLimits(10, 0.001, 0.05);
    } else {
      r.total = new TF1((modelPrefix + "_tot").c_str(), "[0]*TMath::Voigt(x-[1],[2],[3],4)+pol2(4)", gFitMin, gFitMax);
      r.total->SetParameters(maxY, seed.mean, std::max(0.0005, seed.sigma1), std::max(0.001, seed.gamma), b0, 0.0, 0.0);
    }
    r.total->SetParLimits(0, 0.0, 1e9);
    r.total->SetParLimits(1, meanLo, meanHi);
    r.total->SetParLimits(2, 0.0003, 0.03);
    r.total->SetParLimits(3, 0.001, 0.03);

    h->Fit(r.total, "RQ0");

    r.signal = new TF1((modelPrefix + "_sig").c_str(), "[0]*TMath::Voigt(x-[1],[2],[3],4)", gFitMin, gFitMax);
    r.signal->SetParameters(r.total->GetParameter(0), r.total->GetParameter(1), r.total->GetParameter(2), r.total->GetParameter(3));

    if (usePol3Bkg) {
      r.background = new TF1((modelPrefix + "_bkg").c_str(), "pol3(0)+gaus(4)", gFitMin, gFitMax);
      r.background->SetParameters(r.total->GetParameter(4), r.total->GetParameter(5), r.total->GetParameter(6), r.total->GetParameter(7),
                                  r.total->GetParameter(8), r.total->GetParameter(9), r.total->GetParameter(10));
    } else {
      r.background = new TF1((modelPrefix + "_bkg").c_str(), "pol2(0)", gFitMin, gFitMax);
      r.background->SetParameters(r.total->GetParameter(4), r.total->GetParameter(5), r.total->GetParameter(6));
    }
  } else if (r.model == "DoubleGauss") {
    if (usePol3Bkg) {
      r.total = new TF1((modelPrefix + "_tot").c_str(),
                        "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]*exp(-0.5*((x-[1])/[4])^2)+pol3(5)+gaus(9)",
                        gFitMin, gFitMax);
      r.total->SetParameters(0.6 * maxY, seed.mean, std::max(0.001, seed.sigma1),
                             0.4 * maxY, std::max(0.001, seed.sigma2 > 0 ? seed.sigma2 : 0.01),
                             b0, 0.0, 0.0, 0.0, 0.10 * maxY, 0.400);
      r.total->SetParameter(11, 0.010);
      r.total->SetParLimits(9, 0.0, 1e9);
      r.total->SetParLimits(10, 0.390, 0.410);
      r.total->SetParLimits(11, 0.001, 0.05);
    } else {
      r.total = new TF1((modelPrefix + "_tot").c_str(),
                        "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]*exp(-0.5*((x-[1])/[4])^2)+pol2(5)",
                        gFitMin, gFitMax);
      r.total->SetParameters(0.6 * maxY, seed.mean, std::max(0.001, seed.sigma1),
                             0.4 * maxY, std::max(0.001, seed.sigma2 > 0 ? seed.sigma2 : 0.01),
                             b0, 0.0, 0.0);
    }
    r.total->SetParLimits(0, 0.0, 1e9);
    r.total->SetParLimits(3, 0.0, 1e9);
    r.total->SetParLimits(1, meanLo, meanHi);
    r.total->SetParLimits(2, 0.001, 0.05);
    r.total->SetParLimits(4, 0.001, 0.08);

    h->Fit(r.total, "RQ0");

    r.signal = new TF1((modelPrefix + "_sig").c_str(),
                       "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]*exp(-0.5*((x-[1])/[4])^2)",
                       gFitMin, gFitMax);
    r.signal->SetParameters(r.total->GetParameter(0), r.total->GetParameter(1), r.total->GetParameter(2),
                            r.total->GetParameter(3), r.total->GetParameter(4));

    if (usePol3Bkg) {
      r.background = new TF1((modelPrefix + "_bkg").c_str(), "pol3(0)+gaus(4)", gFitMin, gFitMax);
      r.background->SetParameters(r.total->GetParameter(5), r.total->GetParameter(6), r.total->GetParameter(7), r.total->GetParameter(8),
                                  r.total->GetParameter(9), r.total->GetParameter(10), r.total->GetParameter(11));
    } else {
      r.background = new TF1((modelPrefix + "_bkg").c_str(), "pol2(0)", gFitMin, gFitMax);
      r.background->SetParameters(r.total->GetParameter(5), r.total->GetParameter(6), r.total->GetParameter(7));
    }
  } else if (r.model == "QuadGauss") {
    const double s1 = std::max(0.001, seed.sigma1);
    const double s2 = std::max(0.001, seed.sigma2 > 0 ? seed.sigma2 : 0.01);
    const double s3 = std::max(0.001, seed.sigma3 > 0 ? seed.sigma3 : 0.02);
    const double s4 = std::max(0.001, std::max(0.04, 2.0 * s3));

    // Build fixed relative component weights from step4 seed fractions plus one broad tail component.
    const double frac4_seed = 0.08;
    const double f1n = seed.frac1 * (1.0 - frac4_seed);
    const double f2n = seed.frac2 * (1.0 - frac4_seed);
    const double f3n = seed.frac3 * (1.0 - frac4_seed);
    const double f4n = frac4_seed;
    const double c1raw = f1n / s1;
    const double c2raw = f2n / s2;
    const double c3raw = f3n / s3;
    const double c4raw = f4n / s4;
    const double csum = c1raw + c2raw + c3raw + c4raw;
    const double c1 = (csum > 0.0 ? c1raw / csum : 0.45);
    const double c2 = (csum > 0.0 ? c2raw / csum : 0.27);
    const double c3 = (csum > 0.0 ? c3raw / csum : 0.20);
    const double c4 = (csum > 0.0 ? c4raw / csum : 0.08);

    if (usePol4Bkg) {
      r.total = new TF1((modelPrefix + "_tot").c_str(),
                        "[0]*([1]*exp(-0.5*((x-[2])/([3]*[10]))^2)+[4]*exp(-0.5*((x-[2])/([5]*[10]))^2)+"
                        "[6]*exp(-0.5*((x-[2])/([7]*[10]))^2)+[8]*exp(-0.5*((x-[2])/([9]*[10]))^2))+pol4(11)",
                        gFitMin, gFitMax);
      r.total->SetParameters(0.5 * maxY, c1, seed.mean, s1, c2, s2, c3, s3, c4, s4, 1.0);
      r.total->SetParameter(11, b0);
      r.total->SetParameter(12, 0.0);
      r.total->SetParameter(13, 0.0);
      r.total->SetParameter(14, 0.0);
    } else if (useExpoPol2Bkg) {
      r.total = new TF1((modelPrefix + "_tot").c_str(),
                        "[0]*([1]*exp(-0.5*((x-[2])/([3]*[10]))^2)+[4]*exp(-0.5*((x-[2])/([5]*[10]))^2)+"
                        "[6]*exp(-0.5*((x-[2])/([7]*[10]))^2)+[8]*exp(-0.5*((x-[2])/([9]*[10]))^2))+[11]*exp([12]*x)+pol2(13)",
                        gFitMin, gFitMax);
      r.total->SetParameters(0.5 * maxY, c1, seed.mean, s1, c2, s2, c3, s3, c4, s4, 1.0);
      r.total->SetParameter(11, b0);
      r.total->SetParameter(12, -2.0);
      r.total->SetParameter(13, b0);
      r.total->SetParameter(14, 0.0);
      r.total->SetParameter(15, 0.0);
      r.total->SetParLimits(11, -1e9, 1e9);
      r.total->SetParLimits(12, -50.0, 10.0);
    } else if (usePol3PlainBkg) {
      r.total = new TF1((modelPrefix + "_tot").c_str(),
                        "[0]*([1]*exp(-0.5*((x-[2])/([3]*[10]))^2)+[4]*exp(-0.5*((x-[2])/([5]*[10]))^2)+"
                        "[6]*exp(-0.5*((x-[2])/([7]*[10]))^2)+[8]*exp(-0.5*((x-[2])/([9]*[10]))^2))+pol3(11)",
                        gFitMin, gFitMax);
      r.total->SetParameters(0.5 * maxY, c1, seed.mean, s1, c2, s2, c3, s3, c4, s4, 1.0);
      r.total->SetParameter(11, b0);
      r.total->SetParameter(12, 0.0);
      r.total->SetParameter(13, 0.0);
      r.total->SetParameter(14, 0.0);
    } else if (usePol3Bkg) {
      r.total = new TF1((modelPrefix + "_tot").c_str(),
                        "[0]*([1]*exp(-0.5*((x-[2])/([3]*[10]))^2)+[4]*exp(-0.5*((x-[2])/([5]*[10]))^2)+"
                        "[6]*exp(-0.5*((x-[2])/([7]*[10]))^2)+[8]*exp(-0.5*((x-[2])/([9]*[10]))^2))+pol3(11)+gaus(15)",
                        gFitMin, gFitMax);
      r.total->SetParameters(0.5 * maxY, c1, seed.mean, s1, c2, s2, c3, s3, c4, s4, 1.0);
      r.total->SetParameter(11, b0);
      r.total->SetParameter(12, 0.0);
      r.total->SetParameter(13, 0.0);
      r.total->SetParameter(14, 0.0);
      r.total->SetParameter(15, 0.10 * maxY);
      r.total->SetParameter(16, 0.400);
      r.total->SetParameter(17, 0.010);
      r.total->SetParLimits(15, 0.0, 1e9);
      r.total->SetParLimits(16, 0.390, 0.410);
      r.total->SetParLimits(17, 0.001, 0.05);
    } else {
      r.total = new TF1((modelPrefix + "_tot").c_str(),
                        "[0]*([1]*exp(-0.5*((x-[2])/([3]*[10]))^2)+[4]*exp(-0.5*((x-[2])/([5]*[10]))^2)+"
                        "[6]*exp(-0.5*((x-[2])/([7]*[10]))^2)+[8]*exp(-0.5*((x-[2])/([9]*[10]))^2))+pol2(11)",
                        gFitMin, gFitMax);
      r.total->SetParameters(0.5 * maxY, c1, seed.mean, s1, c2, s2, c3, s3, c4, s4, 1.0);
      r.total->SetParameter(11, b0);
      r.total->SetParameter(12, 0.0);
      r.total->SetParameter(13, 0.0);
    }
    r.total->FixParameter(1, c1);
    r.total->FixParameter(2, seed.mean);
    r.total->FixParameter(3, s1);
    r.total->FixParameter(4, c2);
    r.total->FixParameter(5, s2);
    r.total->FixParameter(6, c3);
    r.total->FixParameter(7, s3);
    r.total->FixParameter(8, c4);
    r.total->FixParameter(9, s4);
    r.total->SetParLimits(0, 0.0, 1e9);
    r.total->SetParLimits(10, 0.5, 2.5);

    h->Fit(r.total, "RQ0");

    r.signal = new TF1((modelPrefix + "_sig").c_str(),
                       "[0]*([1]*exp(-0.5*((x-[2])/([3]*[10]))^2)+[4]*exp(-0.5*((x-[2])/([5]*[10]))^2)+"
                       "[6]*exp(-0.5*((x-[2])/([7]*[10]))^2)+[8]*exp(-0.5*((x-[2])/([9]*[10]))^2))",
                       gFitMin, gFitMax);
    r.signal->SetParameters(r.total->GetParameter(0), r.total->GetParameter(1), r.total->GetParameter(2),
                            r.total->GetParameter(3), r.total->GetParameter(4), r.total->GetParameter(5),
                            r.total->GetParameter(6), r.total->GetParameter(7), r.total->GetParameter(8),
                            r.total->GetParameter(9), r.total->GetParameter(10));

    if (usePol4Bkg) {
      r.background = new TF1((modelPrefix + "_bkg").c_str(), "pol4(0)", gFitMin, gFitMax);
      r.background->SetParameters(r.total->GetParameter(11), r.total->GetParameter(12), r.total->GetParameter(13),
                                  r.total->GetParameter(14), r.total->GetParameter(15));
    } else if (useExpoPol2Bkg) {
      r.background = new TF1((modelPrefix + "_bkg").c_str(), "[0]*exp([1]*x)+pol2(2)", gFitMin, gFitMax);
      r.background->SetParameters(r.total->GetParameter(11), r.total->GetParameter(12),
                                  r.total->GetParameter(13), r.total->GetParameter(14), r.total->GetParameter(15));
    } else if (usePol3PlainBkg) {
      r.background = new TF1((modelPrefix + "_bkg").c_str(), "pol3(0)", gFitMin, gFitMax);
      r.background->SetParameters(r.total->GetParameter(11), r.total->GetParameter(12),
                                  r.total->GetParameter(13), r.total->GetParameter(14));
    } else if (usePol3Bkg) {
      r.background = new TF1((modelPrefix + "_bkg").c_str(), "pol3(0)+gaus(4)", gFitMin, gFitMax);
      r.background->SetParameters(r.total->GetParameter(11), r.total->GetParameter(12), r.total->GetParameter(13), r.total->GetParameter(14),
                                  r.total->GetParameter(15), r.total->GetParameter(16), r.total->GetParameter(17));
    } else {
      r.background = new TF1((modelPrefix + "_bkg").c_str(), "pol2(0)", gFitMin, gFitMax);
      r.background->SetParameters(r.total->GetParameter(11), r.total->GetParameter(12), r.total->GetParameter(13));
    }
  } else {
    r.model = "TripleGauss";

    const double s1 = std::max(0.001, seed.sigma1);
    const double s2 = std::max(0.001, seed.sigma2 > 0 ? seed.sigma2 : 0.01);
    const double s3 = std::max(0.001, seed.sigma3 > 0 ? seed.sigma3 : 0.02);

    const double c1raw = seed.frac1 / s1;
    const double c2raw = seed.frac2 / s2;
    const double c3raw = seed.frac3 / s3;
    const double csum = c1raw + c2raw + c3raw;
    const double c1 = (csum > 0.0 ? c1raw / csum : 0.5);
    const double c2 = (csum > 0.0 ? c2raw / csum : 0.3);
    const double c3 = (csum > 0.0 ? c3raw / csum : 0.2);

    if (usePol4Bkg) {
      r.total = new TF1((modelPrefix + "_tot").c_str(),
                        "[0]*([1]*exp(-0.5*((x-[2])/([3]*[8]))^2)+[4]*exp(-0.5*((x-[2])/([5]*[8]))^2)+"
                        "[6]*exp(-0.5*((x-[2])/([7]*[8]))^2))+pol4(9)",
                        gFitMin, gFitMax);
      r.total->SetParameters(0.5 * maxY, c1, seed.mean, s1, c2, s2, c3, s3, 1.0, b0, 0.0);
      r.total->SetParameter(11, 0.0);
      r.total->SetParameter(12, 0.0);
      r.total->SetParameter(13, 0.0);
    } else if (useExpoPol2Bkg) {
      r.total = new TF1((modelPrefix + "_tot").c_str(),
                        "[0]*([1]*exp(-0.5*((x-[2])/([3]*[8]))^2)+[4]*exp(-0.5*((x-[2])/([5]*[8]))^2)+"
                        "[6]*exp(-0.5*((x-[2])/([7]*[8]))^2))+[9]*exp([10]*x)+pol2(11)",
                        gFitMin, gFitMax);
      r.total->SetParameters(0.5 * maxY, c1, seed.mean, s1, c2, s2, c3, s3, 1.0, b0, -2.0);
      r.total->SetParameter(11, b0);
      r.total->SetParameter(12, 0.0);
      r.total->SetParameter(13, 0.0);
      r.total->SetParLimits(9, -1e9, 1e9);
      r.total->SetParLimits(10, -50.0, 10.0);
    } else if (usePol3PlainBkg) {
      r.total = new TF1((modelPrefix + "_tot").c_str(),
                        "[0]*([1]*exp(-0.5*((x-[2])/([3]*[8]))^2)+[4]*exp(-0.5*((x-[2])/([5]*[8]))^2)+"
                        "[6]*exp(-0.5*((x-[2])/([7]*[8]))^2))+pol3(9)",
                        gFitMin, gFitMax);
      r.total->SetParameters(0.5 * maxY, c1, seed.mean, s1, c2, s2, c3, s3, 1.0, b0, 0.0);
      r.total->SetParameter(11, 0.0);
      r.total->SetParameter(12, 0.0);
    } else if (usePol3Bkg) {
      r.total = new TF1((modelPrefix + "_tot").c_str(),
                        "[0]*([1]*exp(-0.5*((x-[2])/([3]*[8]))^2)+[4]*exp(-0.5*((x-[2])/([5]*[8]))^2)+"
                        "[6]*exp(-0.5*((x-[2])/([7]*[8]))^2))+pol3(9)+gaus(13)",
                        gFitMin, gFitMax);
      r.total->SetParameters(0.5 * maxY, c1, seed.mean, s1, c2, s2, c3, s3, 1.0, b0, 0.0);
      r.total->SetParameter(11, 0.0);
      r.total->SetParameter(12, 0.0);
      r.total->SetParameter(13, 0.10 * maxY);
      r.total->SetParameter(14, 0.400);
      r.total->SetParameter(15, 0.010);
      r.total->SetParLimits(13, 0.0, 1e9);
      r.total->SetParLimits(14, 0.390, 0.410);
      r.total->SetParLimits(15, 0.001, 0.05);
    } else {
      r.total = new TF1((modelPrefix + "_tot").c_str(),
                        "[0]*([1]*exp(-0.5*((x-[2])/([3]*[8]))^2)+[4]*exp(-0.5*((x-[2])/([5]*[8]))^2)+"
                        "[6]*exp(-0.5*((x-[2])/([7]*[8]))^2))+pol2(9)",
                        gFitMin, gFitMax);
      r.total->SetParameters(0.5 * maxY, c1, seed.mean, s1, c2, s2, c3, s3, 1.0, b0, 0.0);
      r.total->SetParameter(11, 0.0);
    }
    r.total->FixParameter(1, c1);
    r.total->FixParameter(2, seed.mean);
    r.total->FixParameter(3, s1);
    r.total->FixParameter(4, c2);
    r.total->FixParameter(5, s2);
    r.total->FixParameter(6, c3);
    r.total->FixParameter(7, s3);
    r.total->SetParLimits(0, 0.0, 1e9);
    r.total->SetParLimits(8, 0.5, 2.5);

    h->Fit(r.total, "RQ0");

    r.signal = new TF1((modelPrefix + "_sig").c_str(),
                       "[0]*([1]*exp(-0.5*((x-[2])/([3]*[8]))^2)+[4]*exp(-0.5*((x-[2])/([5]*[8]))^2)+"
                       "[6]*exp(-0.5*((x-[2])/([7]*[8]))^2))",
                       gFitMin, gFitMax);
    r.signal->SetParameters(r.total->GetParameter(0), r.total->GetParameter(1), r.total->GetParameter(2),
                            r.total->GetParameter(3), r.total->GetParameter(4), r.total->GetParameter(5),
                            r.total->GetParameter(6), r.total->GetParameter(7), r.total->GetParameter(8));

    if (usePol4Bkg) {
      r.background = new TF1((modelPrefix + "_bkg").c_str(), "pol4(0)", gFitMin, gFitMax);
      r.background->SetParameters(r.total->GetParameter(9), r.total->GetParameter(10), r.total->GetParameter(11),
                                  r.total->GetParameter(12), r.total->GetParameter(13));
    } else if (useExpoPol2Bkg) {
      r.background = new TF1((modelPrefix + "_bkg").c_str(), "[0]*exp([1]*x)+pol2(2)", gFitMin, gFitMax);
      r.background->SetParameters(r.total->GetParameter(9), r.total->GetParameter(10),
                                  r.total->GetParameter(11), r.total->GetParameter(12), r.total->GetParameter(13));
    } else if (usePol3PlainBkg) {
      r.background = new TF1((modelPrefix + "_bkg").c_str(), "pol3(0)", gFitMin, gFitMax);
      r.background->SetParameters(r.total->GetParameter(9), r.total->GetParameter(10),
                                  r.total->GetParameter(11), r.total->GetParameter(12));
    } else if (usePol3Bkg) {
      r.background = new TF1((modelPrefix + "_bkg").c_str(), "pol3(0)+gaus(4)", gFitMin, gFitMax);
      r.background->SetParameters(r.total->GetParameter(9), r.total->GetParameter(10), r.total->GetParameter(11), r.total->GetParameter(12),
                                  r.total->GetParameter(13), r.total->GetParameter(14), r.total->GetParameter(15));
    } else {
      r.background = new TF1((modelPrefix + "_bkg").c_str(), "pol2(0)", gFitMin, gFitMax);
      r.background->SetParameters(r.total->GetParameter(9), r.total->GetParameter(10), r.total->GetParameter(11));
    }
  }

  if (r.total->GetNDF() > 0) r.chi2ndf = r.total->GetChisquare() / r.total->GetNDF();
  if (r.model == "TripleGauss" || r.model == "QuadGauss") {
    r.mean = r.total->GetParameter(2);
    r.meanErr = r.total->GetParError(2);
  } else {
    r.mean = r.total->GetParameter(1);
    r.meanErr = r.total->GetParError(1);
  }

  double y = 0.0;
  double pull2 = 0.0;
  double pullAbs = 0.0;
  int npull = 0;
  for (int i = 1; i <= h->GetNbinsX(); ++i) {
    const double x = h->GetBinCenter(i);
    if (x < gFitMin || x > gFitMax) continue;
    y += std::max(0.0, r.signal->Eval(x));
    const double obs = h->GetBinContent(i);
    const double err = h->GetBinError(i);
    if (err > 0.0) {
      const double p = (obs - r.total->Eval(x)) / err;
      pull2 += p * p;
      pullAbs += std::fabs(p);
      r.pullMaxAbs = std::max(r.pullMaxAbs, std::fabs(p));
      ++npull;
    }
  }
  r.yield = y;
  r.yieldErr = std::sqrt(std::max(0.0, y));
  if (npull > 0) {
    r.pullRms = std::sqrt(pull2 / npull);
    r.pullMeanAbs = pullAbs / npull;
  }

  return r;
}

void drawFitPage(TCanvas& c, TH1D* h, const FitResult& r, const std::string& title, const std::string& pdfPath) {
  TH1D* hDisp = (TH1D*)h->Clone((std::string(h->GetName()) + "_display").c_str());
  if (kDisplayRebin > 1 && hDisp->GetNbinsX() >= kDisplayRebin) hDisp->Rebin(kDisplayRebin);
  c.Clear();
  c.Divide(1, 2);

  c.cd(1);
  gPad->SetPad(0.0, 0.32, 1.0, 1.0);
  gPad->SetLeftMargin(0.14);
  gPad->SetRightMargin(0.04);
  gPad->SetBottomMargin(0.03);
  hDisp->SetStats(0);
  hDisp->SetLineWidth(2);
  hDisp->SetTitle(title.c_str());
  hDisp->GetXaxis()->SetTitle("m(#pi^{+}#pi^{-}) [GeV]");
  hDisp->GetYaxis()->SetTitle("Opposite-charge pairs / rebinned bin");
  hDisp->Draw("E");

  r.total->SetLineColor(kBlue + 1);
  r.total->SetLineWidth(2);
  r.total->Draw("same");

  r.signal->SetLineColor(kRed + 1);
  r.signal->SetLineWidth(2);
  r.signal->SetLineStyle(1);
  r.signal->Draw("same");

  r.background->SetLineColor(kGreen + 2);
  r.background->SetLineWidth(2);
  r.background->SetLineStyle(2);
  r.background->Draw("same");

  TLegend leg(0.58, 0.62, 0.89, 0.88);
  leg.SetBorderSize(0);
  leg.AddEntry(hDisp, "MC pairs", "lep");
  leg.AddEntry(r.total, "Total fit", "l");
  leg.AddEntry(r.signal, "Signal component", "l");
  leg.AddEntry(r.background, "Background component", "l");
  leg.Draw();

  TPaveText txt(0.12, 0.62, 0.54, 0.88, "NDC");
  txt.SetBorderSize(0);
  txt.SetFillStyle(0);
  txt.AddText(("Signal model: " + r.model).c_str());
  txt.AddText(Form("#chi^{2}/ndf = %.3f", r.chi2ndf));
  txt.AddText(Form("pull RMS = %.3f, max|pull| = %.3f", r.pullRms, r.pullMaxAbs));
  txt.AddText(Form("mean = %.6f #pm %.6f GeV", r.mean, r.meanErr));
  txt.AddText(Form("signal yield (fit window) = %.0f", r.yield));
  txt.Draw();

  c.cd(2);
  gPad->SetPad(0.0, 0.0, 1.0, 0.32);
  gPad->SetLeftMargin(0.14);
  gPad->SetRightMargin(0.04);
  gPad->SetTopMargin(0.03);
  gPad->SetBottomMargin(0.34);
  TH1D hPull((std::string("hPull_") + h->GetName()).c_str(), "", hDisp->GetNbinsX(), hDisp->GetXaxis()->GetXmin(), hDisp->GetXaxis()->GetXmax());
  hPull.SetStats(0);
  hPull.SetMarkerStyle(20);
  hPull.SetMarkerSize(0.65);
  hPull.SetMarkerColor(kBlack);
  for (int i = 1; i <= hDisp->GetNbinsX(); ++i) {
    const double x = hDisp->GetBinCenter(i);
    if (x < gFitMin || x > gFitMax) continue;
    const double y = hDisp->GetBinContent(i);
    const double ey = hDisp->GetBinError(i);
    if (ey <= 0.0) continue;
    hPull.SetBinContent(i, (y - r.total->Eval(x)) / ey);
    hPull.SetBinError(i, 1.0);
  }
  hPull.SetMinimum(-5.0);
  hPull.SetMaximum(5.0);
  hPull.GetYaxis()->SetTitle("Pull ((Obs-Fit)/#sigma)");
  hPull.GetXaxis()->SetTitle("m(#pi^{+}#pi^{-}) [GeV]");
  hPull.GetYaxis()->SetNdivisions(505);
  hPull.GetYaxis()->SetTitleSize(0.11);
  hPull.GetYaxis()->SetLabelSize(0.085);
  hPull.GetYaxis()->SetTitleOffset(0.55);
  hPull.GetXaxis()->SetTitleSize(0.11);
  hPull.GetXaxis()->SetLabelSize(0.085);
  hPull.GetXaxis()->SetTitleOffset(1.05);
  hPull.Draw("E1P");
  TLine l0(gMassMin, 0.0, gMassMax, 0.0);
  l0.SetLineStyle(2);
  l0.Draw();

  c.Print(pdfPath.c_str());
  delete hDisp;
}

}  // namespace

void step5_mc_random_pairs_kshort(const char* input = "Samples/merged_mc_v2.2.root",
                                  const char* step4SeedCsv = "Reports/step4_a0025_best_fit_params.csv",
                                  const char* outDir = "Reports",
                                  const char* reportPrefix = "step5",
                                  const char* datasetLabel = "MC",
                                  const char* signalModelOverride = "auto",
                                  const char* bkgMode = "baseline",
                                  double massMin = 0.40,
                                  double massMax = 0.70,
                                  double fitMin = 0.40,
                                  double fitMax = 0.70,
                                  int nBins = 280,
                                  double piPMin = -1.0,
                                  double piPMax = -1.0) {
  gMassMin = massMin;
  gMassMax = massMax;
  gFitMin = fitMin;
  gFitMax = fitMax;
  gBins = nBins;
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gSystem->mkdir(outDir, kTRUE);
  const std::string outDirStr(outDir);
  const std::string prefix(reportPrefix);
  gSystem->mkdir((outDirStr + "/" + prefix).c_str(), kTRUE);

  const auto seeds = loadStep4Seeds(step4SeedCsv);

  auto getSeed = [&](const std::string& cat) {
    auto it = seeds.find(cat);
    if (it != seeds.end()) return it->second;
    Step4Seed s;
    return s;
  };

  TFile f(input, "READ");
  if (f.IsZombie()) {
    std::cerr << "Failed to open input file: " << input << std::endl;
    return;
  }

  TTree* tree = nullptr;
  f.GetObject("Tree", tree);
  if (!tree) {
    std::cerr << "Cannot find TTree 'Tree' in " << input << std::endl;
    return;
  }

  Long64_t nReco = 0;
  double recoPx[kMaxReco], recoPy[kMaxReco], recoPz[kMaxReco], recoCharge[kMaxReco];
  Long64_t recoGoodTrack[kMaxReco], recoPIDPion[kMaxReco];

  tree->SetBranchAddress("NReco", &nReco);
  tree->SetBranchAddress("RecoPx", recoPx);
  tree->SetBranchAddress("RecoPy", recoPy);
  tree->SetBranchAddress("RecoPz", recoPz);
  tree->SetBranchAddress("RecoCharge", recoCharge);
  tree->SetBranchAddress("RecoGoodTrack", recoGoodTrack);
  tree->SetBranchAddress("RecoPIDPion", recoPIDPion);

  TH1D h0("hStep5_0tag", "Step5 0-tag; m(#pi^{+}#pi^{-}) [GeV]; Pairs / bin", gBins, gMassMin, gMassMax);
  TH1D h1("hStep5_1tag", "Step5 1-tag; m(#pi^{+}#pi^{-}) [GeV]; Pairs / bin", gBins, gMassMin, gMassMax);
  TH1D h2("hStep5_2tag", "Step5 2-tag; m(#pi^{+}#pi^{-}) [GeV]; Pairs / bin", gBins, gMassMin, gMassMax);

  Long64_t pairTotal = 0;
  Long64_t pairInMass = 0;

  const Long64_t nEntries = tree->GetEntries();
  for (Long64_t ie = 0; ie < nEntries; ++ie) {
    tree->GetEntry(ie);
    if (nReco > kMaxReco) continue;

    std::vector<int> pos;
    std::vector<int> neg;
    pos.reserve(32);
    neg.reserve(32);

    for (int i = 0; i < nReco; ++i) {
      if (recoGoodTrack[i] == 0) continue;
      if (!passTrackAcceptance(recoPx[i], recoPy[i], recoPz[i])) continue;
      if (recoCharge[i] > 0) pos.push_back(i);
      if (recoCharge[i] < 0) neg.push_back(i);
    }

    for (int ip : pos) {
      for (int in : neg) {
        const double pPlus = std::sqrt(recoPx[ip] * recoPx[ip] + recoPy[ip] * recoPy[ip] + recoPz[ip] * recoPz[ip]);
        const double pMinus = std::sqrt(recoPx[in] * recoPx[in] + recoPy[in] * recoPy[in] + recoPz[in] * recoPz[in]);
        if (piPMin >= 0.0 && piPMax > piPMin) {
          // Strict momentum-binned mode: both daughters must be in the same requested pion-momentum bin.
          if (pPlus < piPMin || pPlus >= piPMax) continue;
          if (pMinus < piPMin || pMinus >= piPMax) continue;
        }
        ++pairTotal;
        const double m = buildMass(recoPx[ip], recoPy[ip], recoPz[ip], recoPx[in], recoPy[in], recoPz[in]);
        if (m < gMassMin || m > gMassMax) continue;
        ++pairInMass;

        const int nTag = (recoPIDPion[ip] >= 2 ? 1 : 0) + (recoPIDPion[in] >= 2 ? 1 : 0);
        if (nTag == 0) h0.Fill(m);
        if (nTag == 1) h1.Fill(m);
        if (nTag == 2) h2.Fill(m);
      }
    }
  }

  FitResult r0 = fitCategory(&h0, getSeed("0tag"), "0tag", "step5_0tag", signalModelOverride, bkgMode);
  FitResult r1 = fitCategory(&h1, getSeed("1tag"), "1tag", "step5_1tag", signalModelOverride, bkgMode);
  FitResult r2 = fitCategory(&h2, getSeed("2tag"), "2tag", "step5_2tag", signalModelOverride, bkgMode);

  const double n1 = r1.yield;
  const double n2 = r2.yield;
  const double d = n1 + 2.0 * n2;
  const double eff = (d > 0.0 ? (2.0 * n2 / d) : 0.0);

  const double s1 = std::max(1.0, r1.yieldErr);
  const double s2 = std::max(1.0, r2.yieldErr);
  const double de_dn1 = (d > 0.0 ? (-2.0 * n2 / (d * d)) : 0.0);
  const double de_dn2 = (d > 0.0 ? (2.0 * n1 / (d * d)) : 0.0);
  const double effErr = std::sqrt(de_dn1 * de_dn1 * s1 * s1 + de_dn2 * de_dn2 * s2 * s2);

  const std::string pdfPath = outDirStr + "/" + prefix + "_kshort_random_pairs_report.pdf";

  TCanvas c("cStep5", "cStep5", 1100, 800);
  c.Print((pdfPath + "[").c_str());

  c.Clear();
  TPaveText cover(0.06, 0.08, 0.94, 0.92, "NDC");
  cover.SetBorderSize(0);
  cover.SetFillStyle(0);
  cover.AddText(Form("%s (%s): all opposite-charge good-track pairs with pion mass hypothesis", reportPrefix, datasetLabel));
  cover.AddText(Form("Input: %s", input));
  cover.AddText("Selection: RecoGoodTrack, 0.15<=|cos#theta|<=0.675, opposite charge; no KShort branches used");
  if (piPMin >= 0.0 && piPMax > piPMin) {
    cover.AddText(Form("Momentum bin (strict): both daughters %.3f <= p_{#pi} < %.3f GeV", piPMin, piPMax));
  } else {
    cover.AddText("Momentum bin: inclusive pion momentum");
  }
  cover.AddText("Tag WP: RecoPIDPion>=2 on each daughter");
  cover.AddText(Form("Signal override: %s; Background mode: %s", signalModelOverride, bkgMode));
  cover.AddText("Fit model: step4-seeded signal + category background");
  cover.AddText(Form("Step4 seed file: %s", step4SeedCsv));
  cover.Draw();
  c.Print(pdfPath.c_str());

  drawFitPage(c, &h0, r0, "Step5 fit: 0-tag", pdfPath);
  drawFitPage(c, &h1, r1, "Step5 fit: 1-tag", pdfPath);
  drawFitPage(c, &h2, r2, "Step5 fit: 2-tag", pdfPath);

  c.Clear();
  TPaveText summ(0.06, 0.08, 0.94, 0.92, "NDC");
  summ.SetBorderSize(0);
  summ.SetFillStyle(0);
  summ.SetTextAlign(12);
  summ.AddText("Pair counts");
  summ.AddText(Form("  Total opposite-charge good-track pairs: %lld", pairTotal));
  summ.AddText(Form("  Pairs in mass histogram range [%.2f, %.2f] GeV: %lld", gMassMin, gMassMax, pairInMass));
  summ.AddText(" ");
  summ.AddText(Form("Fitted signal yields in fit window [%.2f, %.2f] GeV", gFitMin, gFitMax));
  summ.AddText(Form("  0-tag: N = %.0f", r0.yield));
  summ.AddText(Form("  1-tag: N = %.0f", r1.yield));
  summ.AddText(Form("  2-tag: N = %.0f", r2.yield));
  summ.AddText(" ");
  summ.AddText(Form("Tagging efficiency: eps = 2*N2/(N1+2*N2) = %.5f +/- %.5f", eff, effErr));
  summ.Draw();
  c.Print(pdfPath.c_str());

  c.Print((pdfPath + "]").c_str());

  std::ofstream txt(outDirStr + "/" + prefix + "_fit_summary.txt");
  txt << std::fixed << std::setprecision(6);
  txt << "Input," << input << "\n";
  txt << "Step4Seed," << step4SeedCsv << "\n";
  txt << "Selection,RecoGoodTrack&&0.15<=|cos(theta)|<=0.675&&opposite_charge\n";
  if (piPMin >= 0.0 && piPMax > piPMin) {
    txt << "MomentumBin,strict_both_daughters," << piPMin << "," << piPMax << "\n";
  } else {
    txt << "MomentumBin,inclusive\n";
  }
  txt << "TagWP,RecoPIDPion>=2\n";
  txt << "SignalOverride," << signalModelOverride << "\n";
  txt << "BackgroundMode," << bkgMode << "\n";
  txt << "MassRange," << gMassMin << "," << gMassMax << "\n";
  txt << "FitRange," << gFitMin << "," << gFitMax << "\n";
  txt << "NBins," << gBins << "\n";
  txt << "PairTotal," << pairTotal << "\n";
  txt << "PairInMassRange," << pairInMass << "\n\n";

  auto dump = [&](const std::string& cat, const FitResult& r) {
    txt << cat << ",model," << r.model << ",chi2ndf," << r.chi2ndf << ",mean," << r.mean
        << ",meanErr," << r.meanErr << ",yield," << r.yield << ",yieldErr," << r.yieldErr
        << ",pullRms," << r.pullRms << ",pullMeanAbs," << r.pullMeanAbs << ",pullMaxAbs," << r.pullMaxAbs << "\n";
  };
  dump("0tag", r0);
  dump("1tag", r1);
  dump("2tag", r2);
  txt << "\ntagging_efficiency," << eff << ",tagging_efficiency_err," << effErr << "\n";
  txt.close();

  std::ofstream csv(outDirStr + "/" + prefix + "_fit_results.csv");
  csv << "category,signal_model,chi2ndf,mean,mean_err,signal_yield,signal_yield_err,pull_rms,pull_mean_abs,pull_max_abs\n";
  csv << "0tag," << r0.model << "," << r0.chi2ndf << "," << r0.mean << "," << r0.meanErr << "," << r0.yield << "," << r0.yieldErr
      << "," << r0.pullRms << "," << r0.pullMeanAbs << "," << r0.pullMaxAbs << "\n";
  csv << "1tag," << r1.model << "," << r1.chi2ndf << "," << r1.mean << "," << r1.meanErr << "," << r1.yield << "," << r1.yieldErr
      << "," << r1.pullRms << "," << r1.pullMeanAbs << "," << r1.pullMaxAbs << "\n";
  csv << "2tag," << r2.model << "," << r2.chi2ndf << "," << r2.mean << "," << r2.meanErr << "," << r2.yield << "," << r2.yieldErr
      << "," << r2.pullRms << "," << r2.pullMeanAbs << "," << r2.pullMaxAbs << "\n";
  csv << "tagging_efficiency,,," << eff << "," << effErr << ",,\n";
  csv.close();

  TFile hout((outDirStr + "/" + prefix + "/" + prefix + "_hists.root").c_str(), "RECREATE");
  h0.Write();
  h1.Write();
  h2.Write();
  hout.Close();

  std::cout << "Step 5 complete. Report: " << pdfPath << std::endl;

  delete r0.total;
  delete r0.signal;
  delete r0.background;
  delete r1.total;
  delete r1.signal;
  delete r1.background;
  delete r2.total;
  delete r2.signal;
  delete r2.background;
}
