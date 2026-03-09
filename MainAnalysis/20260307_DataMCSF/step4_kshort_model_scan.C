#include <TFile.h>
#include <TF1.h>
#include <TH1D.h>
#include <TROOT.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <fstream>
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

struct FitRes {
  std::string category;
  std::string model;
  double chi2ndf = 1e9;
  int ndf = -1;
  double mean = 0.0;
  double sigma1 = 0.0;
  double sigma2 = 0.0;
  double sigma3 = 0.0;
};

FitRes runModel(TH1D* h, const std::string& cat, const std::string& model,
                double fitMin, double fitMax) {
  FitRes r;
  r.category = cat;
  r.model = model;

  const double maxY = std::max(100.0, h->GetMaximum());
  const double b0 = 0.5 * (h->GetBinContent(1) + h->GetBinContent(h->GetNbinsX()));
  const std::string fname = "f_" + cat + "_" + model;
  TF1 f;

  if (model == "Gauss") {
    f = TF1(fname.c_str(), "gaus(0)", fitMin, fitMax);
    f.SetParameters(maxY, 0.4976, 0.006);
    f.SetParLimits(1, 0.485, 0.510);
    f.SetParLimits(2, 0.001, 0.04);
  } else if (model == "Voigt") {
    f = TF1(fname.c_str(), "[0]*TMath::Voigt(x-[1],[2],[3],4)", fitMin, fitMax);
    f.SetParameters(maxY, 0.4976, 0.002, 0.006);
    f.SetParLimits(1, 0.485, 0.510);
    f.SetParLimits(2, 0.0005, 0.03);
    f.SetParLimits(3, 0.001, 0.03);
  } else if (model == "DoubleGauss") {
    f = TF1(fname.c_str(),
            "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]*exp(-0.5*((x-[1])/[4])^2)",
            fitMin, fitMax);
    f.SetParameters(0.65 * maxY, 0.4976, 0.004, 0.35 * maxY, 0.012);
    f.SetParLimits(1, 0.485, 0.510);
    f.SetParLimits(2, 0.001, 0.04);
    f.SetParLimits(4, 0.001, 0.08);
  } else if (model == "TripleGauss") {
    f = TF1(fname.c_str(),
            "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]*exp(-0.5*((x-[1])/[4])^2)+"
            "[5]*exp(-0.5*((x-[1])/[6])^2)",
            fitMin, fitMax);
    f.SetParameters(0.5 * maxY, 0.4976, 0.003, 0.3 * maxY, 0.009, 0.2 * maxY, 0.03);
    f.SetParLimits(1, 0.485, 0.510);
    f.SetParLimits(2, 0.001, 0.04);
    f.SetParLimits(4, 0.001, 0.08);
    f.SetParLimits(6, 0.001, 0.12);
  } else if (model == "QuadGauss") {
    f = TF1(fname.c_str(),
            "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]*exp(-0.5*((x-[1])/[4])^2)+"
            "[5]*exp(-0.5*((x-[1])/[6])^2)+[7]*exp(-0.5*((x-[1])/[8])^2)",
            fitMin, fitMax);
    f.SetParameters(0.42 * maxY, 0.4976, 0.003, 0.28 * maxY, 0.008, 0.20 * maxY, 0.02, 0.10 * maxY, 0.05);
    f.SetParLimits(1, 0.485, 0.510);
    f.SetParLimits(2, 0.001, 0.04);
    f.SetParLimits(4, 0.001, 0.08);
    f.SetParLimits(6, 0.001, 0.12);
    f.SetParLimits(8, 0.001, 0.20);
  } else if (model == "BifurGauss") {
    f = TF1(fname.c_str(),
            "[0]*(x<[1]?exp(-0.5*((x-[1])/[2])^2):exp(-0.5*((x-[1])/[3])^2))",
            fitMin, fitMax);
    f.SetParameters(maxY, 0.4976, 0.004, 0.012);
    f.SetParLimits(1, 0.485, 0.510);
    f.SetParLimits(2, 0.001, 0.05);
    f.SetParLimits(3, 0.001, 0.08);
  } else if (model == "VoigtPlusGauss") {
    f = TF1(fname.c_str(),
            "[0]*TMath::Voigt(x-[1],[2],[3],4)+gaus(4)",
            fitMin, fitMax);
    f.SetParameters(0.7 * maxY, 0.4976, 0.002, 0.006, 0.3 * maxY, 0.4976, 0.012);
    f.SetParLimits(1, 0.485, 0.510);
    f.SetParLimits(2, 0.0005, 0.03);
    f.SetParLimits(3, 0.001, 0.03);
    f.SetParLimits(5, 0.485, 0.510);
    f.SetParLimits(6, 0.001, 0.08);
  } else if (model == "GaussPol2") {
    f = TF1(fname.c_str(), "gaus(0)+pol2(3)", fitMin, fitMax);
    f.SetParameters(maxY, 0.4976, 0.006, b0, 0.0, 0.0);
    f.SetParLimits(1, 0.485, 0.510);
    f.SetParLimits(2, 0.001, 0.04);
  } else if (model == "VoigtPol2") {
    f = TF1(fname.c_str(), "[0]*TMath::Voigt(x-[1],[2],[3],4)+pol2(4)", fitMin, fitMax);
    f.SetParameters(maxY, 0.4976, 0.002, 0.006, b0, 0.0, 0.0);
    f.SetParLimits(1, 0.485, 0.510);
    f.SetParLimits(2, 0.0005, 0.03);
    f.SetParLimits(3, 0.001, 0.03);
  } else if (model == "DoubleGaussPol2") {
    f = TF1(fname.c_str(),
            "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]*exp(-0.5*((x-[1])/[4])^2)+pol2(5)",
            fitMin, fitMax);
    f.SetParameters(0.65 * maxY, 0.4976, 0.004, 0.35 * maxY, 0.012, b0, 0.0, 0.0);
    f.SetParLimits(1, 0.485, 0.510);
    f.SetParLimits(2, 0.001, 0.04);
    f.SetParLimits(4, 0.001, 0.08);
  } else if (model == "TripleGaussPol2") {
    f = TF1(fname.c_str(),
            "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]*exp(-0.5*((x-[1])/[4])^2)+"
            "[5]*exp(-0.5*((x-[1])/[6])^2)+pol2(7)",
            fitMin, fitMax);
    f.SetParameters(0.5 * maxY, 0.4976, 0.003, 0.3 * maxY, 0.009, 0.2 * maxY, 0.03, b0, 0.0, 0.0);
    f.SetParLimits(1, 0.485, 0.510);
    f.SetParLimits(2, 0.001, 0.04);
    f.SetParLimits(4, 0.001, 0.08);
    f.SetParLimits(6, 0.001, 0.12);
  } else if (model == "TripleGaussPol3") {
    f = TF1(fname.c_str(),
            "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]*exp(-0.5*((x-[1])/[4])^2)+"
            "[5]*exp(-0.5*((x-[1])/[6])^2)+pol3(7)",
            fitMin, fitMax);
    f.SetParameters(0.5 * maxY, 0.4976, 0.003, 0.3 * maxY, 0.009, 0.2 * maxY, 0.03, b0, 0.0, 0.0, 0.0);
    f.SetParLimits(1, 0.485, 0.510);
    f.SetParLimits(2, 0.001, 0.04);
    f.SetParLimits(4, 0.001, 0.08);
    f.SetParLimits(6, 0.001, 0.12);
  } else {
    return r;
  }

  h->Fit(&f, "RQ0");
  if (f.GetNDF() > 0) {
    r.ndf = f.GetNDF();
    r.chi2ndf = f.GetChisquare() / f.GetNDF();
    if (f.GetNpar() > 1) r.mean = f.GetParameter(1);
    if (f.GetNpar() > 2) r.sigma1 = std::fabs(f.GetParameter(2));
    if (model.find("DoubleGauss") != std::string::npos || model.find("TripleGauss") != std::string::npos) {
      if (f.GetNpar() > 4) r.sigma2 = std::fabs(f.GetParameter(4));
    }
    if (model.find("TripleGauss") != std::string::npos && f.GetNpar() > 6) {
      r.sigma3 = std::fabs(f.GetParameter(6));
    }
  }
  return r;
}
}  // namespace

void step4_kshort_model_scan(const char* input = "Samples/merged_mc_v2.2.root",
                             const char* outDir = "Reports",
                             double angleCut = 0.025,
                             const char* prefix = "step4_modelscan_fr030_100",
                             double massMin = 0.30,
                             double massMax = 1.00,
                             double fitMin = 0.30,
                             double fitMax = 1.00) {
  gROOT->SetBatch(kTRUE);

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

  TH1D h0("hScan0", "0tag", kBins, massMin, massMax);
  TH1D h1("hScan1", "1tag", kBins, massMin, massMax);
  TH1D h2("hScan2", "2tag", kBins, massMin, massMax);

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
      const double m = buildMass(recoPx[i1], recoPy[i1], recoPz[i1],
                                 recoPx[i2], recoPy[i2], recoPz[i2]);
      if (m < massMin || m > massMax) continue;
      const int nTag = (recoPIDPion[i1] >= 2 ? 1 : 0) + (recoPIDPion[i2] >= 2 ? 1 : 0);
      if (nTag == 0) h0.Fill(m);
      if (nTag == 1) h1.Fill(m);
      if (nTag == 2) h2.Fill(m);
    }
  }

  std::vector<std::string> models = {
      "Gauss", "Voigt", "DoubleGauss", "TripleGauss", "QuadGauss",
      "BifurGauss", "VoigtPlusGauss"};

  std::vector<FitRes> all;
  for (const auto& m : models) all.push_back(runModel(&h0, "0tag", m, fitMin, fitMax));
  for (const auto& m : models) all.push_back(runModel(&h1, "1tag", m, fitMin, fitMax));
  for (const auto& m : models) all.push_back(runModel(&h2, "2tag", m, fitMin, fitMax));

  std::string outCsv = std::string(outDir) + "/" + prefix + "_results.csv";
  std::ofstream csv(outCsv);
  csv << "category,model,chi2ndf,ndf,mean,sigma1,sigma2,sigma3\n";
  for (const auto& r : all) {
    csv << r.category << "," << r.model << "," << r.chi2ndf << "," << r.ndf << ","
        << r.mean << "," << r.sigma1 << "," << r.sigma2 << "," << r.sigma3 << "\n";
  }
  csv.close();

  std::string outTxt = std::string(outDir) + "/" + prefix + "_summary.txt";
  std::ofstream txt(outTxt);
  txt << "Step4 Kshort model scan\n";
  txt << "Input," << input << "\n";
  txt << "AngleCut," << angleCut << "\n";
  txt << "MassRange," << massMin << "," << massMax << "\n";
  txt << "FitRange," << fitMin << "," << fitMax << "\n\n";
  for (const char* cat : {"0tag", "1tag", "2tag"}) {
    std::vector<FitRes> v;
    for (const auto& r : all) if (r.category == cat && r.ndf > 0) v.push_back(r);
    std::sort(v.begin(), v.end(), [](const FitRes& a, const FitRes& b) { return a.chi2ndf < b.chi2ndf; });
    txt << cat << " best models:\n";
    for (size_t i = 0; i < v.size() && i < 5; ++i) {
      txt << "  " << (i + 1) << ". " << v[i].model << " chi2/ndf=" << v[i].chi2ndf << "\n";
    }
    txt << "\n";
  }
  txt.close();

  std::cout << "Wrote " << outCsv << std::endl;
  std::cout << "Wrote " << outTxt << std::endl;
}
