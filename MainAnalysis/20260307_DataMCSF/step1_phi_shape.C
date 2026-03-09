#include <TCanvas.h>
#include <TFile.h>
#include <TF1.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMath.h>
#include <TPaveText.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

namespace {
constexpr double kKaonMass = 0.493677;
constexpr double kAbsCosMin = 0.15;
constexpr double kAbsCosMax = 0.675;
constexpr double kMassMin = 0.99;
constexpr double kMassMax = 1.06;
constexpr int kBins = 280;
constexpr double kFitMin = 1.00;
constexpr double kFitMax = 1.05;
constexpr int kMaxReco = 256;
constexpr int kMaxPhi = 16;

struct FitOutcome {
  std::string model;
  double chi2ndf = 1e9;
  double mean = 0.0;
  double meanErr = 0.0;
  double width = 0.0;
  double widthErr = 0.0;
};

double buildMassFromRecoTracks(double px1, double py1, double pz1,
                               double px2, double py2, double pz2) {
  const double p1sq = px1 * px1 + py1 * py1 + pz1 * pz1;
  const double p2sq = px2 * px2 + py2 * py2 + pz2 * pz2;
  const double e1 = std::sqrt(p1sq + kKaonMass * kKaonMass);
  const double e2 = std::sqrt(p2sq + kKaonMass * kKaonMass);

  const double e = e1 + e2;
  const double px = px1 + px2;
  const double py = py1 + py2;
  const double pz = pz1 + pz2;
  const double m2 = e * e - (px * px + py * py + pz * pz);
  return (m2 > 0.0) ? std::sqrt(m2) : 0.0;
}

bool passTrackAcceptance(double px, double py, double pz) {
  const double p = std::sqrt(px * px + py * py + pz * pz);
  if (p <= 0.0) return false;
  const double abscos = std::fabs(pz / p);
  return (abscos >= kAbsCosMin && abscos <= kAbsCosMax);
}

FitOutcome fitAndPickBest(TH1D* h, const std::string& tag, const std::string& outDir) {
  FitOutcome best;

  TF1 fGauss("fGauss", "gaus(0)", kFitMin, kFitMax);
  fGauss.SetParNames("A", "mean", "sigma");
  fGauss.SetParameters(std::max(100.0, h->GetMaximum()), 1.0195, 0.0035);
  fGauss.SetParLimits(1, 1.015, 1.024);
  fGauss.SetParLimits(2, 0.001, 0.02);

  TH1D* hClone1 = (TH1D*)h->Clone((tag + "_g").c_str());
  hClone1->Fit(&fGauss, "RQ0");
  if (fGauss.GetNDF() > 0) {
    const double c = fGauss.GetChisquare() / fGauss.GetNDF();
    if (c < best.chi2ndf) {
      best.model = "Gauss";
      best.chi2ndf = c;
      best.mean = fGauss.GetParameter(1);
      best.meanErr = fGauss.GetParError(1);
      best.width = fGauss.GetParameter(2);
      best.widthErr = fGauss.GetParError(2);
    }
  }

  TF1 fVoigt("fVoigt", "[0]*TMath::Voigt(x-[1],[2],[3],4)", kFitMin, kFitMax);
  fVoigt.SetParNames("A", "mean", "sigma", "gamma");
  fVoigt.SetParameters(std::max(100.0, h->GetMaximum()), 1.0195, 0.0018, 0.0043);
  fVoigt.SetParLimits(1, 1.015, 1.024);
  fVoigt.SetParLimits(2, 0.0003, 0.01);
  fVoigt.SetParLimits(3, 0.001, 0.01);

  TH1D* hClone2 = (TH1D*)h->Clone((tag + "_v").c_str());
  hClone2->Fit(&fVoigt, "RQ0");
  if (fVoigt.GetNDF() > 0) {
    const double c = fVoigt.GetChisquare() / fVoigt.GetNDF();
    if (c < best.chi2ndf) {
      best.model = "Voigt";
      best.chi2ndf = c;
      best.mean = fVoigt.GetParameter(1);
      best.meanErr = fVoigt.GetParError(1);
      best.width = fVoigt.GetParameter(2);
      best.widthErr = fVoigt.GetParError(2);
    }
  }

  TF1 fDouble("fDouble",
              "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]*exp(-0.5*((x-[1])/[4])^2)",
              kFitMin, kFitMax);
  fDouble.SetParNames("A1", "mean", "sigma1", "A2", "sigma2");
  fDouble.SetParameters(std::max(60.0, 0.6 * h->GetMaximum()), 1.0195, 0.0025,
                        std::max(40.0, 0.4 * h->GetMaximum()), 0.0055);
  fDouble.SetParLimits(1, 1.015, 1.024);
  fDouble.SetParLimits(2, 0.0008, 0.02);
  fDouble.SetParLimits(4, 0.0008, 0.03);

  TH1D* hClone3 = (TH1D*)h->Clone((tag + "_d").c_str());
  hClone3->Fit(&fDouble, "RQ0");
  if (fDouble.GetNDF() > 0) {
    const double c = fDouble.GetChisquare() / fDouble.GetNDF();
    if (c < best.chi2ndf) {
      best.model = "DoubleGauss";
      best.chi2ndf = c;
      best.mean = fDouble.GetParameter(1);
      best.meanErr = fDouble.GetParError(1);
      best.width = fDouble.GetParameter(2);
      best.widthErr = fDouble.GetParError(2);
    }
  }

  TF1 fTriple("fTriple",
              "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]*exp(-0.5*((x-[1])/[4])^2)+"
              "[5]*exp(-0.5*((x-[1])/[6])^2)",
              kFitMin, kFitMax);
  fTriple.SetParNames("A1", "mean", "sigma1", "A2", "sigma2", "A3", "sigma3");
  fTriple.SetParameters(std::max(45.0, 0.5 * h->GetMaximum()), 1.0195, 0.0022,
                        std::max(35.0, 0.3 * h->GetMaximum()), 0.0045,
                        std::max(20.0, 0.2 * h->GetMaximum()), 0.0075);
  fTriple.SetParLimits(1, 1.015, 1.024);
  fTriple.SetParLimits(2, 0.0008, 0.02);
  fTriple.SetParLimits(4, 0.0008, 0.03);
  fTriple.SetParLimits(6, 0.0010, 0.05);

  TH1D* hClone4 = (TH1D*)h->Clone((tag + "_t").c_str());
  hClone4->Fit(&fTriple, "RQ0");
  if (fTriple.GetNDF() > 0) {
    const double c = fTriple.GetChisquare() / fTriple.GetNDF();
    if (c < best.chi2ndf) {
      best.model = "TripleGauss";
      best.chi2ndf = c;
      best.mean = fTriple.GetParameter(1);
      best.meanErr = fTriple.GetParError(1);
      best.width = fTriple.GetParameter(2);
      best.widthErr = fTriple.GetParError(2);
    }
  }

  TCanvas c(("c_" + tag).c_str(), ("c_" + tag).c_str(), 900, 800);
  c.Divide(1, 2);
  c.cd(1);
  gPad->SetPad(0.0, 0.32, 1.0, 1.0);
  gPad->SetLeftMargin(0.14);
  gPad->SetRightMargin(0.04);
  gPad->SetBottomMargin(0.03);
  h->SetLineWidth(2);
  h->SetStats(0);
  h->GetXaxis()->SetTitle("m(K^{+}K^{-}) [GeV]");
  h->GetYaxis()->SetTitle("Candidates / bin");
  h->Draw("E");

  fGauss.SetLineColor(kBlue + 1);
  fGauss.SetLineWidth(2);
  fGauss.Draw("same");

  fVoigt.SetLineColor(kRed + 1);
  fVoigt.SetLineWidth(2);
  fVoigt.Draw("same");

  fDouble.SetLineColor(kGreen + 2);
  fDouble.SetLineWidth(2);
  fDouble.Draw("same");

  fTriple.SetLineColor(kMagenta + 1);
  fTriple.SetLineWidth(2);
  fTriple.Draw("same");

  TLegend leg(0.60, 0.68, 0.89, 0.88);
  leg.SetBorderSize(0);
  leg.AddEntry(h, tag.c_str(), "lep");
  leg.AddEntry(&fGauss, "Gauss", "l");
  leg.AddEntry(&fVoigt, "Voigt", "l");
  leg.AddEntry(&fDouble, "DoubleGauss", "l");
  leg.AddEntry(&fTriple, "TripleGauss", "l");
  leg.Draw();

  TPaveText txt(0.14, 0.74, 0.54, 0.88, "NDC");
  txt.SetBorderSize(0);
  txt.SetFillStyle(0);
  txt.AddText(("Best: " + best.model).c_str());
  txt.AddText(Form("#chi^{2}/ndf = %.3f", best.chi2ndf));
  txt.AddText(Form("mean = %.6f #pm %.6f GeV", best.mean, best.meanErr));
  txt.AddText(Form("width = %.6f #pm %.6f GeV", best.width, best.widthErr));
  txt.Draw();

  c.cd(2);
  gPad->SetPad(0.0, 0.0, 1.0, 0.32);
  gPad->SetLeftMargin(0.14);
  gPad->SetRightMargin(0.04);
  gPad->SetTopMargin(0.03);
  gPad->SetBottomMargin(0.34);
  TH1D hPull(("hPull_" + tag).c_str(), "", h->GetNbinsX(), h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());
  hPull.SetStats(0);
  hPull.SetMarkerStyle(20);
  hPull.SetMarkerSize(0.65);
  hPull.SetLineColor(kBlack);
  hPull.SetMarkerColor(kBlack);
  for (int i = 1; i <= h->GetNbinsX(); ++i) {
    const double x = h->GetBinCenter(i);
    if (x < kFitMin || x > kFitMax) continue;
    const double y = h->GetBinContent(i);
    const double ey = h->GetBinError(i);
    if (ey <= 0.0) continue;
    double yfit = fGauss.Eval(x);
    if (best.model == "Voigt") yfit = fVoigt.Eval(x);
    if (best.model == "DoubleGauss") yfit = fDouble.Eval(x);
    if (best.model == "TripleGauss") yfit = fTriple.Eval(x);
    hPull.SetBinContent(i, (y - yfit) / ey);
    hPull.SetBinError(i, 1.0);
  }
  hPull.GetYaxis()->SetTitle("Pull ((Obs-Fit)/#sigma)");
  hPull.GetXaxis()->SetTitle("m(K^{+}K^{-}) [GeV]");
  hPull.GetYaxis()->SetNdivisions(505);
  hPull.GetYaxis()->SetTitleSize(0.11);
  hPull.GetYaxis()->SetLabelSize(0.085);
  hPull.GetYaxis()->SetTitleOffset(0.55);
  hPull.GetXaxis()->SetTitleSize(0.11);
  hPull.GetXaxis()->SetLabelSize(0.085);
  hPull.GetXaxis()->SetTitleOffset(1.05);
  hPull.SetMinimum(-5.0);
  hPull.SetMaximum(5.0);
  hPull.Draw("E1P");
  TLine l0(kMassMin, 0.0, kMassMax, 0.0);
  l0.SetLineStyle(2);
  l0.Draw();

  c.SaveAs((outDir + "/step1_" + tag + ".pdf").c_str());
  c.SaveAs((outDir + "/step1_" + tag + ".png").c_str());

  delete hClone1;
  delete hClone2;
  delete hClone3;
  delete hClone4;
  return best;
}

}  // namespace

void step1_phi_shape(const char* input = "Samples/merged_mc_v2.2.root",
                     const char* outDir = "Reports/step1") {
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  gSystem->mkdir(outDir, kTRUE);

  TFile f(input, "READ");
  if (f.IsZombie()) {
    std::cerr << "Failed to open input file: " << input << std::endl;
    return;
  }

  TTree* tree = nullptr;
  f.GetObject("Tree", tree);
  if (!tree) {
    std::cerr << "Cannot find TTree named 'Tree' in " << input << std::endl;
    return;
  }

  Long64_t nReco = 0;
  Long64_t nPhi = 0;
  double recoPx[kMaxReco], recoPy[kMaxReco], recoPz[kMaxReco], recoCharge[kMaxReco];
  Long64_t recoGoodTrack[kMaxReco], recoPIDKaon[kMaxReco];
  Long64_t phiReco1ID[kMaxPhi], phiReco2ID[kMaxPhi];
  double phiReco1Angle[kMaxPhi], phiReco2Angle[kMaxPhi];

  tree->SetBranchAddress("NReco", &nReco);
  tree->SetBranchAddress("NPhi", &nPhi);
  tree->SetBranchAddress("RecoPx", recoPx);
  tree->SetBranchAddress("RecoPy", recoPy);
  tree->SetBranchAddress("RecoPz", recoPz);
  tree->SetBranchAddress("RecoCharge", recoCharge);
  tree->SetBranchAddress("RecoGoodTrack", recoGoodTrack);
  tree->SetBranchAddress("RecoPIDKaon", recoPIDKaon);
  tree->SetBranchAddress("PhiReco1ID[NPhi]", phiReco1ID);
  tree->SetBranchAddress("PhiReco2ID[NPhi]", phiReco2ID);
  tree->SetBranchAddress("PhiReco1Angle[NPhi]", phiReco1Angle);
  tree->SetBranchAddress("PhiReco2Angle[NPhi]", phiReco2Angle);

  TH1D h0("h0", "0-tag: neither daughter passes K PID (RecoPIDKaon>=2)", kBins, kMassMin, kMassMax);
  TH1D h1("h1", "1-tag: exactly one daughter passes K PID (RecoPIDKaon>=2)", kBins, kMassMin, kMassMax);
  TH1D h2("h2", "2-tag: both daughters pass K PID (RecoPIDKaon>=2)", kBins, kMassMin, kMassMax);

  Long64_t nPhiTotal = 0;
  Long64_t nPhiAccepted = 0;
  const Long64_t nEntries = tree->GetEntries();

  for (Long64_t ie = 0; ie < nEntries; ++ie) {
    tree->GetEntry(ie);

    if (nReco > kMaxReco || nPhi > kMaxPhi) {
      std::cerr << "Array bounds exceeded at event " << ie << ": NReco=" << nReco
                << " NPhi=" << nPhi << std::endl;
      continue;
    }

    for (Long64_t i = 0; i < nPhi; ++i) {
      ++nPhiTotal;

      const Long64_t i1 = phiReco1ID[i];
      const Long64_t i2 = phiReco2ID[i];
      if (i1 < 0 || i2 < 0) continue;
      if (i1 >= nReco || i2 >= nReco) continue;

      if (phiReco1Angle[i] >= 0.01 || phiReco2Angle[i] >= 0.01) continue;
      if (recoGoodTrack[i1] == 0 || recoGoodTrack[i2] == 0) continue;
      if (!passTrackAcceptance(recoPx[i1], recoPy[i1], recoPz[i1])) continue;
      if (!passTrackAcceptance(recoPx[i2], recoPy[i2], recoPz[i2])) continue;
      if (recoCharge[i1] * recoCharge[i2] >= 0) continue;

      const double m = buildMassFromRecoTracks(recoPx[i1], recoPy[i1], recoPz[i1],
                                               recoPx[i2], recoPy[i2], recoPz[i2]);
      if (m <= 0.0) continue;

      ++nPhiAccepted;
      const int nTagged = (recoPIDKaon[i1] >= 2 ? 1 : 0) + (recoPIDKaon[i2] >= 2 ? 1 : 0);
      if (nTagged == 0) h0.Fill(m);
      if (nTagged == 1) h1.Fill(m);
      if (nTagged == 2) h2.Fill(m);
    }
  }

  const std::string outDirStr = outDir;
  FitOutcome r0 = fitAndPickBest(&h0, "pid0", outDirStr);
  FitOutcome r1 = fitAndPickBest(&h1, "pid1", outDirStr);
  FitOutcome r2 = fitAndPickBest(&h2, "pid2", outDirStr);

  std::ofstream summary(outDirStr + "/step1_fit_summary.txt");
  summary << std::fixed << std::setprecision(6);
  summary << "Input: " << input << "\n";
  summary << "Selection: valid PhiReco IDs, PhiReco1Angle<0.01, PhiReco2Angle<0.01, both RecoGoodTrack, "
             "0.15<=|cos(theta)|<=0.675 on both daughters, opposite charge, tag if RecoPIDKaon>=2\n";
  summary << "Fit mode: signal-only models (Gauss, Voigt, DoubleGauss, TripleGauss), fit range 1.00-1.05 GeV\n";
  summary << "Total phi candidates scanned: " << nPhiTotal << "\n";
  summary << "Accepted gen-matched phi candidates: " << nPhiAccepted << "\n\n";

  auto dump = [&](const char* label, TH1D& h, const FitOutcome& r) {
    summary << "Category " << label << "\n";
    summary << "  Entries: " << h.GetEntries() << "\n";
    summary << "  Best model: " << r.model << "\n";
    summary << "  Chi2/NDF: " << r.chi2ndf << "\n";
    summary << "  Mean [GeV]: " << r.mean << " +/- " << r.meanErr << "\n";
    summary << "  Width [GeV]: " << r.width << " +/- " << r.widthErr << "\n\n";
  };

  dump("0-tag", h0, r0);
  dump("1-tag", h1, r1);
  dump("2-tag", h2, r2);
  summary.close();

  TFile out((outDirStr + "/step1_hists.root").c_str(), "RECREATE");
  h0.Write();
  h1.Write();
  h2.Write();
  out.Close();

  std::cout << "Step 1 finished. Outputs in: " << outDir << std::endl;
}
