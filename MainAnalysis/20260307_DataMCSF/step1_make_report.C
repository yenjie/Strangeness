#include <TCanvas.h>
#include <TFile.h>
#include <TF1.h>
#include <TH1D.h>
#include <TH1I.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMath.h>
#include <TPad.h>
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
constexpr int kDisplayRebin = 4;
constexpr int kMaxReco = 256;
constexpr int kMaxPhi = 16;

struct ModelFit {
  std::string name;
  double chi2ndf = 1e9;
  double p0 = 0.0;
  double p1 = 0.0;
  double p2 = 0.0;
  double p3 = 0.0;
  double p4 = 0.0;
  double p5 = 0.0;
  double p6 = 0.0;
  double e1 = 0.0;
  double e2 = 0.0;
  double e3 = 0.0;
  double e4 = 0.0;
};

struct FitSet {
  ModelFit gauss;
  ModelFit voigt;
  ModelFit dgauss;
  ModelFit tgauss;
  ModelFit best;
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

FitSet fitMass(TH1D* h) {
  FitSet out;
  out.gauss.name = "Gauss";
  out.voigt.name = "Voigt";
  out.dgauss.name = "DoubleGauss";
  out.tgauss.name = "TripleGauss";

  TF1 fG("fG", "gaus(0)", kFitMin, kFitMax);
  fG.SetParameters(std::max(100.0, h->GetMaximum()), 1.0195, 0.0035);
  fG.SetParLimits(1, 1.015, 1.024);
  fG.SetParLimits(2, 0.001, 0.02);

  TH1D* hc1 = (TH1D*)h->Clone("hc1_tmp");
  hc1->Fit(&fG, "RQ0");
  if (fG.GetNDF() > 0) {
    out.gauss.chi2ndf = fG.GetChisquare() / fG.GetNDF();
    out.gauss.p0 = fG.GetParameter(0);
    out.gauss.p1 = fG.GetParameter(1);
    out.gauss.p2 = fG.GetParameter(2);
    out.gauss.e1 = fG.GetParError(1);
    out.gauss.e2 = fG.GetParError(2);
  }

  TF1 fV("fV", "[0]*TMath::Voigt(x-[1],[2],[3],4)", kFitMin, kFitMax);
  fV.SetParameters(std::max(100.0, h->GetMaximum()), 1.0195, 0.0018, 0.0043);
  fV.SetParLimits(1, 1.015, 1.024);
  fV.SetParLimits(2, 0.0003, 0.01);
  fV.SetParLimits(3, 0.001, 0.01);

  TH1D* hc2 = (TH1D*)h->Clone("hc2_tmp");
  hc2->Fit(&fV, "RQ0");
  if (fV.GetNDF() > 0) {
    out.voigt.chi2ndf = fV.GetChisquare() / fV.GetNDF();
    out.voigt.p0 = fV.GetParameter(0);
    out.voigt.p1 = fV.GetParameter(1);
    out.voigt.p2 = fV.GetParameter(2);
    out.voigt.p3 = fV.GetParameter(3);
    out.voigt.e1 = fV.GetParError(1);
    out.voigt.e2 = fV.GetParError(2);
    out.voigt.e3 = fV.GetParError(3);
  }

  TF1 fD("fD",
         "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]*exp(-0.5*((x-[1])/[4])^2)",
         kFitMin, kFitMax);
  fD.SetParameters(std::max(60.0, 0.6 * h->GetMaximum()), 1.0195, 0.0025,
                   std::max(40.0, 0.4 * h->GetMaximum()), 0.0055);
  fD.SetParLimits(1, 1.015, 1.024);
  fD.SetParLimits(2, 0.0008, 0.02);
  fD.SetParLimits(4, 0.0008, 0.03);

  TH1D* hc3 = (TH1D*)h->Clone("hc3_tmp");
  hc3->Fit(&fD, "RQ0");
  if (fD.GetNDF() > 0) {
    out.dgauss.chi2ndf = fD.GetChisquare() / fD.GetNDF();
    out.dgauss.p0 = fD.GetParameter(0);
    out.dgauss.p1 = fD.GetParameter(1);
    out.dgauss.p2 = fD.GetParameter(2);
    out.dgauss.p3 = fD.GetParameter(3);
    out.dgauss.p4 = fD.GetParameter(4);
    out.dgauss.e1 = fD.GetParError(1);
    out.dgauss.e2 = fD.GetParError(2);
    out.dgauss.e3 = fD.GetParError(4);
  }

  TF1 fT("fT",
         "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]*exp(-0.5*((x-[1])/[4])^2)+"
         "[5]*exp(-0.5*((x-[1])/[6])^2)",
         kFitMin, kFitMax);
  fT.SetParameters(std::max(45.0, 0.5 * h->GetMaximum()), 1.0195, 0.0022,
                   std::max(35.0, 0.3 * h->GetMaximum()), 0.0045,
                   std::max(20.0, 0.2 * h->GetMaximum()), 0.0075);
  fT.SetParLimits(1, 1.015, 1.024);
  fT.SetParLimits(2, 0.0008, 0.02);
  fT.SetParLimits(4, 0.0008, 0.03);
  fT.SetParLimits(6, 0.0010, 0.05);

  TH1D* hc4 = (TH1D*)h->Clone("hc4_tmp");
  hc4->Fit(&fT, "RQ0");
  if (fT.GetNDF() > 0) {
    out.tgauss.chi2ndf = fT.GetChisquare() / fT.GetNDF();
    out.tgauss.p0 = fT.GetParameter(0);
    out.tgauss.p1 = fT.GetParameter(1);
    out.tgauss.p2 = fT.GetParameter(2);
    out.tgauss.p3 = fT.GetParameter(3);
    out.tgauss.p4 = fT.GetParameter(4);
    out.tgauss.p5 = fT.GetParameter(5);
    out.tgauss.p6 = fT.GetParameter(6);
    out.tgauss.e1 = fT.GetParError(1);
    out.tgauss.e2 = fT.GetParError(2);
    out.tgauss.e3 = fT.GetParError(4);
    out.tgauss.e4 = fT.GetParError(6);
  }

  out.best = out.gauss;
  if (out.voigt.chi2ndf < out.best.chi2ndf) out.best = out.voigt;
  if (out.dgauss.chi2ndf < out.best.chi2ndf) out.best = out.dgauss;
  if (out.tgauss.chi2ndf < out.best.chi2ndf) out.best = out.tgauss;

  delete hc1;
  delete hc2;
  delete hc3;
  delete hc4;
  return out;
}

void drawMassPage(TCanvas& c, TH1D* h, const FitSet& fit, const std::string& title,
                  const std::string& pdfPath) {
  TH1D* hDisp = static_cast<TH1D*>(h->Clone((std::string(h->GetName()) + "_disp").c_str()));
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
  hDisp->GetXaxis()->SetTitle("m(K^{+}K^{-}) [GeV]");
  hDisp->GetYaxis()->SetTitle("Candidates / bin");
  hDisp->Draw("E");

  TF1 fG("fG_draw", "gaus(0)", kFitMin, kFitMax);
  fG.SetParameters(fit.gauss.p0, fit.gauss.p1, fit.gauss.p2);
  fG.SetLineColor(kBlue + 1);
  fG.SetLineWidth(2);
  fG.Draw("same");

  TF1 fV("fV_draw", "[0]*TMath::Voigt(x-[1],[2],[3],4)", kFitMin, kFitMax);
  fV.SetParameters(fit.voigt.p0, fit.voigt.p1, fit.voigt.p2, fit.voigt.p3);
  fV.SetLineColor(kRed + 1);
  fV.SetLineWidth(2);
  fV.Draw("same");

  TF1 fD("fD_draw",
         "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]*exp(-0.5*((x-[1])/[4])^2)",
         kFitMin, kFitMax);
  fD.SetParameters(fit.dgauss.p0, fit.dgauss.p1, fit.dgauss.p2, fit.dgauss.p3, fit.dgauss.p4);
  fD.SetLineColor(kGreen + 2);
  fD.SetLineWidth(2);
  fD.Draw("same");

  TF1 fT("fT_draw",
         "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]*exp(-0.5*((x-[1])/[4])^2)+"
         "[5]*exp(-0.5*((x-[1])/[6])^2)",
         kFitMin, kFitMax);
  fT.SetParameters(fit.tgauss.p0, fit.tgauss.p1, fit.tgauss.p2, fit.tgauss.p3, fit.tgauss.p4,
                   fit.tgauss.p5, fit.tgauss.p6);
  fT.SetLineColor(kMagenta + 1);
  fT.SetLineWidth(2);
  fT.Draw("same");

  TLegend leg(0.58, 0.57, 0.89, 0.88);
  leg.SetBorderSize(0);
  leg.AddEntry(hDisp, "MC points", "lep");
  leg.AddEntry(&fG, "Gauss", "l");
  leg.AddEntry(&fV, "Voigt", "l");
  leg.AddEntry(&fD, "DoubleGauss", "l");
  leg.AddEntry(&fT, "TripleGauss", "l");
  leg.Draw();

  TPaveText txt(0.12, 0.50, 0.58, 0.88, "NDC");
  txt.SetBorderSize(0);
  txt.SetFillStyle(0);
  txt.AddText(Form("Entries: %.0f", hDisp->GetEntries()));
  txt.AddText(Form("Gauss #chi^{2}/ndf = %.3f", fit.gauss.chi2ndf));
  txt.AddText(Form("Voigt #chi^{2}/ndf = %.3f", fit.voigt.chi2ndf));
  txt.AddText(Form("DoubleGauss #chi^{2}/ndf = %.3f", fit.dgauss.chi2ndf));
  txt.AddText(Form("TripleGauss #chi^{2}/ndf = %.3f", fit.tgauss.chi2ndf));
  txt.AddText(("Best: " + fit.best.name).c_str());
  txt.Draw();

  c.cd(2);
  gPad->SetPad(0.0, 0.0, 1.0, 0.32);
  gPad->SetLeftMargin(0.14);
  gPad->SetRightMargin(0.04);
  gPad->SetTopMargin(0.03);
  gPad->SetBottomMargin(0.34);
  TH1D hPull("hPull", "", hDisp->GetNbinsX(), hDisp->GetXaxis()->GetXmin(), hDisp->GetXaxis()->GetXmax());
  hPull.SetStats(0);
  hPull.SetMarkerStyle(20);
  hPull.SetMarkerSize(0.65);
  hPull.SetMarkerColor(kBlack);
  for (int i = 1; i <= hDisp->GetNbinsX(); ++i) {
    const double x = hDisp->GetBinCenter(i);
    if (x < kFitMin || x > kFitMax) continue;
    const double y = hDisp->GetBinContent(i);
    const double ey = hDisp->GetBinError(i);
    if (ey <= 0) continue;
    double yfit = 0.0;
    if (fit.best.name == "Voigt") {
      yfit = fV.Eval(x);
    } else if (fit.best.name == "DoubleGauss") {
      yfit = fD.Eval(x);
    } else if (fit.best.name == "TripleGauss") {
      yfit = fT.Eval(x);
    } else {
      yfit = fG.Eval(x);
    }
    hPull.SetBinContent(i, (y - yfit) / ey);
    hPull.SetBinError(i, 1.0);
  }
  hPull.SetLineColor(kBlack);
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

  c.Print(pdfPath.c_str());
  delete hDisp;
}

}  // namespace

void step1_make_report(const char* input = "Samples/merged_mc_v2.2.root",
                       const char* outDir = "Reports",
                       double angleCut = 0.01,
                       const char* prefix = "step1") {
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gSystem->mkdir(outDir, kTRUE);
  const std::string outDirStr(outDir);
  const std::string prefixStr(prefix);
  gSystem->mkdir((outDirStr + "/" + prefixStr).c_str(), kTRUE);

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

  TH1D hMass0("hMass0", "0-tag (RecoPIDKaon>=2)", kBins, kMassMin, kMassMax);
  TH1D hMass1("hMass1", "1-tag (RecoPIDKaon>=2)", kBins, kMassMin, kMassMax);
  TH1D hMass2("hMass2", "2-tag (RecoPIDKaon>=2)", kBins, kMassMin, kMassMax);
  TH1D hMassAll("hMassAll", "Accepted #phi candidates", kBins, kMassMin, kMassMax);
  TH1I hNTag("hNTag", "Kaon PID tags per #phi candidate;N tagged daughters;Candidates", 3, -0.5, 2.5);

  TH1D hPhiP("hPhiP", "#phi momentum;p_{#phi} [GeV];Candidates", 120, 0, 6);
  TH1D hPhiPt("hPhiPt", "#phi p_{T};p_{T,#phi} [GeV];Candidates", 120, 0, 4);
  TH1D hPhiCos("hPhiCos", "#phi cos#theta;cos#theta_{#phi};Candidates", 100, -1.0, 1.0);

  TH1D hKP("hKP", "Daughter K momentum;p_{K} [GeV];Tracks", 120, 0, 5);
  TH1D hKPt("hKPt", "Daughter K p_{T};p_{T,K} [GeV];Tracks", 120, 0, 3);
  TH1D hKCos("hKCos", "Daughter K cos#theta;cos#theta_{K};Tracks", 100, -1.0, 1.0);

  TH1D hAng1All("hAng1All", "PhiReco1Angle;PhiReco1Angle;Candidates", 100, 0.0, 0.05);
  TH1D hAng2All("hAng2All", "PhiReco2Angle;PhiReco2Angle;Candidates", 100, 0.0, 0.05);
  TH1D hAng1Pass("hAng1Pass", "PhiReco1Angle (pass);PhiReco1Angle;Candidates", 100, 0.0, 0.02);
  TH1D hAng2Pass("hAng2Pass", "PhiReco2Angle (pass);PhiReco2Angle;Candidates", 100, 0.0, 0.02);

  Long64_t nPhiTotal = 0;
  Long64_t nIdValid = 0;
  Long64_t nAnglePass = 0;
  Long64_t nGoodTrackPass = 0;
  Long64_t nAcceptancePass = 0;
  Long64_t nOppChargePass = 0;
  Long64_t nAccepted = 0;

  const Long64_t nEntries = tree->GetEntries();
  for (Long64_t ie = 0; ie < nEntries; ++ie) {
    tree->GetEntry(ie);
    if (nReco > kMaxReco || nPhi > kMaxPhi) continue;

    for (Long64_t i = 0; i < nPhi; ++i) {
      ++nPhiTotal;
      const Long64_t i1 = phiReco1ID[i];
      const Long64_t i2 = phiReco2ID[i];
      if (i1 < 0 || i2 < 0 || i1 >= nReco || i2 >= nReco) continue;
      ++nIdValid;

      hAng1All.Fill(phiReco1Angle[i]);
      hAng2All.Fill(phiReco2Angle[i]);
      const bool passAngle = (phiReco1Angle[i] < angleCut && phiReco2Angle[i] < angleCut);
      if (!passAngle) continue;
      ++nAnglePass;
      hAng1Pass.Fill(phiReco1Angle[i]);
      hAng2Pass.Fill(phiReco2Angle[i]);

      if (recoGoodTrack[i1] == 0 || recoGoodTrack[i2] == 0) continue;
      ++nGoodTrackPass;
      if (!passTrackAcceptance(recoPx[i1], recoPy[i1], recoPz[i1])) continue;
      if (!passTrackAcceptance(recoPx[i2], recoPy[i2], recoPz[i2])) continue;
      ++nAcceptancePass;
      if (recoCharge[i1] * recoCharge[i2] >= 0) continue;
      ++nOppChargePass;

      const double m = buildMassFromRecoTracks(recoPx[i1], recoPy[i1], recoPz[i1],
                                               recoPx[i2], recoPy[i2], recoPz[i2]);
      if (m <= 0.0) continue;
      ++nAccepted;

      const int nTag = (recoPIDKaon[i1] >= 2 ? 1 : 0) + (recoPIDKaon[i2] >= 2 ? 1 : 0);
      hNTag.Fill(nTag);
      hMassAll.Fill(m);
      if (nTag == 0) hMass0.Fill(m);
      if (nTag == 1) hMass1.Fill(m);
      if (nTag == 2) hMass2.Fill(m);

      const double px1 = recoPx[i1], py1 = recoPy[i1], pz1 = recoPz[i1];
      const double px2 = recoPx[i2], py2 = recoPy[i2], pz2 = recoPz[i2];
      const double p1 = std::sqrt(px1 * px1 + py1 * py1 + pz1 * pz1);
      const double p2 = std::sqrt(px2 * px2 + py2 * py2 + pz2 * pz2);
      const double pt1 = std::sqrt(px1 * px1 + py1 * py1);
      const double pt2 = std::sqrt(px2 * px2 + py2 * py2);
      if (p1 > 0) hKCos.Fill(pz1 / p1);
      if (p2 > 0) hKCos.Fill(pz2 / p2);
      hKP.Fill(p1);
      hKP.Fill(p2);
      hKPt.Fill(pt1);
      hKPt.Fill(pt2);

      const double px = px1 + px2;
      const double py = py1 + py2;
      const double pz = pz1 + pz2;
      const double p = std::sqrt(px * px + py * py + pz * pz);
      const double pt = std::sqrt(px * px + py * py);
      hPhiP.Fill(p);
      hPhiPt.Fill(pt);
      if (p > 0) hPhiCos.Fill(pz / p);
    }
  }

  FitSet f0 = fitMass(&hMass0);
  FitSet f1 = fitMass(&hMass1);
  FitSet f2 = fitMass(&hMass2);

  const std::string pdfPath = outDirStr + "/" + prefixStr + "_phi_tagging_report.pdf";
  TCanvas c("cReport", "cReport", 1100, 800);
  c.Print((pdfPath + "[").c_str());

  c.Clear();
  TPaveText cover(0.07, 0.08, 0.93, 0.92, "NDC");
  cover.SetBorderSize(0);
  cover.SetFillStyle(0);
  cover.AddText("Step 1 Report: #phi #rightarrow K^{+}K^{-} signal shape from gen-matched candidates (MC)");
  cover.AddText(Form("Input file: %s", input));
  cover.AddText("Selection:");
  cover.AddText(Form("  valid PhiReco1ID/PhiReco2ID, PhiReco1Angle<%.3f, PhiReco2Angle<%.3f,", angleCut, angleCut));
  cover.AddText("  both RecoGoodTrack, 0.15<=|cos#theta|<=0.675 for both daughters,");
  cover.AddText("  opposite charge, kaon-mass hypothesis for m(KK)");
  cover.AddText("PID split by RecoPIDKaon>=2 on two daughters: 0-tag / 1-tag / 2-tag");
  cover.AddText("Fit range: 1.00-1.05 GeV");
  cover.AddText("Signal-only models compared: Gauss, Voigt, DoubleGauss, TripleGauss");
  cover.Draw();
  c.Print(pdfPath.c_str());

  drawMassPage(c, &hMass0, f0, "Mass fit: 0-tag category", pdfPath);
  drawMassPage(c, &hMass1, f1, "Mass fit: 1-tag category", pdfPath);
  drawMassPage(c, &hMass2, f2, "Mass fit: 2-tag category", pdfPath);

  c.Clear();
  TPaveText summary(0.05, 0.05, 0.95, 0.95, "NDC");
  summary.SetBorderSize(0);
  summary.SetFillStyle(0);
  summary.SetTextAlign(12);
  summary.AddText("Event-flow counts");
  summary.AddText(Form("  Scanned #phi candidates: %lld", nPhiTotal));
  summary.AddText(Form("  + valid daughter IDs: %lld", nIdValid));
  summary.AddText(Form("  + angle cut (<%.3f on both): %lld", angleCut, nAnglePass));
  summary.AddText(Form("  + both RecoGoodTrack: %lld", nGoodTrackPass));
  summary.AddText(Form("  + acceptance 0.15<=|cos#theta|<=0.675 on both: %lld", nAcceptancePass));
  summary.AddText(Form("  + opposite charge: %lld", nOppChargePass));
  summary.AddText(Form("  Accepted for fits: %lld", nAccepted));
  summary.AddText(" ");
  summary.AddText("Best-fit parameters by PID category");
  summary.AddText(Form("  0-tag: %s, chi2/ndf=%.3f, mean=%.6f+-%.6f GeV, width=%.6f+-%.6f GeV",
                       f0.best.name.c_str(), f0.best.chi2ndf, f0.best.p1, f0.best.e1, f0.best.p2, f0.best.e2));
  summary.AddText(Form("  1-tag: %s, chi2/ndf=%.3f, mean=%.6f+-%.6f GeV, width=%.6f+-%.6f GeV",
                       f1.best.name.c_str(), f1.best.chi2ndf, f1.best.p1, f1.best.e1, f1.best.p2, f1.best.e2));
  summary.AddText(Form("  2-tag: %s, chi2/ndf=%.3f, mean=%.6f+-%.6f GeV, width=%.6f+-%.6f GeV",
                       f2.best.name.c_str(), f2.best.chi2ndf, f2.best.p1, f2.best.e1, f2.best.p2, f2.best.e2));
  summary.Draw();
  c.Print(pdfPath.c_str());

  c.Clear();
  c.Divide(2, 1);
  c.cd(1);
  hMassAll.SetLineWidth(2);
  hMassAll.Draw("hist");
  c.cd(2);
  hNTag.GetXaxis()->SetBinLabel(1, "0");
  hNTag.GetXaxis()->SetBinLabel(2, "1");
  hNTag.GetXaxis()->SetBinLabel(3, "2");
  hNTag.SetFillColor(kAzure - 9);
  hNTag.Draw("hist");
  c.Print(pdfPath.c_str());

  c.Clear();
  c.Divide(2, 2);
  c.cd(1);
  hPhiP.SetLineColor(kBlue + 1);
  hPhiP.SetLineWidth(2);
  hPhiP.Draw("hist");
  c.cd(2);
  hPhiPt.SetLineColor(kBlue + 1);
  hPhiPt.SetLineWidth(2);
  hPhiPt.Draw("hist");
  c.cd(3);
  hPhiCos.SetLineColor(kBlue + 1);
  hPhiCos.SetLineWidth(2);
  hPhiCos.Draw("hist");
  c.Print(pdfPath.c_str());

  c.Clear();
  c.Divide(2, 2);
  c.cd(1);
  hKP.SetLineColor(kRed + 1);
  hKP.SetLineWidth(2);
  hKP.Draw("hist");
  c.cd(2);
  hKPt.SetLineColor(kRed + 1);
  hKPt.SetLineWidth(2);
  hKPt.Draw("hist");
  c.cd(3);
  hKCos.SetLineColor(kRed + 1);
  hKCos.SetLineWidth(2);
  hKCos.Draw("hist");
  c.Print(pdfPath.c_str());

  c.Clear();
  c.Divide(2, 2);
  c.cd(1);
  gPad->SetLogy();
  hAng1All.SetLineWidth(2);
  hAng1All.Draw("hist");
  TLine lcut1(angleCut, 1e-1, angleCut, std::max(1.0, hAng1All.GetMaximum() * 1.2));
  lcut1.SetLineColor(kRed + 1);
  lcut1.SetLineStyle(2);
  lcut1.Draw();
  c.cd(2);
  gPad->SetLogy();
  hAng2All.SetLineWidth(2);
  hAng2All.Draw("hist");
  TLine lcut2(angleCut, 1e-1, angleCut, std::max(1.0, hAng2All.GetMaximum() * 1.2));
  lcut2.SetLineColor(kRed + 1);
  lcut2.SetLineStyle(2);
  lcut2.Draw();
  c.cd(3);
  gPad->SetLogy(0);
  hAng1Pass.SetLineColor(kGreen + 2);
  hAng1Pass.SetLineWidth(2);
  hAng1Pass.Draw("hist");
  c.cd(4);
  hAng2Pass.SetLineColor(kGreen + 2);
  hAng2Pass.SetLineWidth(2);
  hAng2Pass.Draw("hist");
  c.Print(pdfPath.c_str());

  c.Print((pdfPath + "]").c_str());

  std::ofstream txt(outDirStr + "/" + prefixStr + "_fit_parameters.txt");
  txt << std::fixed << std::setprecision(6);
  txt << "Input: " << input << "\n";
  txt << "Selection: valid IDs, PhiReco1Angle<" << angleCut << ", PhiReco2Angle<" << angleCut << ", both RecoGoodTrack, "
         "0.15<=|cos(theta)|<=0.675 on both daughters, opposite charge, tag if RecoPIDKaon>=2\n\n";
  txt << "Fit mode: signal-only models (Gauss, Voigt, DoubleGauss, TripleGauss), fit range 1.00-1.05 GeV\n\n";
  txt << "Flow counts\n";
  txt << "scanned_phi," << nPhiTotal << "\n";
  txt << "valid_id," << nIdValid << "\n";
  txt << "pass_angle," << nAnglePass << "\n";
  txt << "pass_goodtrack," << nGoodTrackPass << "\n";
  txt << "pass_acceptance_cosTheta," << nAcceptancePass << "\n";
  txt << "pass_opposite_charge," << nOppChargePass << "\n";
  txt << "accepted," << nAccepted << "\n\n";

  auto writeModel = [&](const std::string& cat, const ModelFit& m) {
    txt << cat << "," << m.name << ",chi2ndf," << m.chi2ndf << ",mean," << m.p1 << ",meanErr," << m.e1
        << ",width," << m.p2 << ",widthErr," << m.e2;
    if (m.name == "Voigt") txt << ",gamma," << m.p3 << ",gammaErr," << m.e3;
    if (m.name == "DoubleGauss") txt << ",sigma2," << m.p4 << ",sigma2Err," << m.e3;
    if (m.name == "TripleGauss") {
      txt << ",sigma2," << m.p4 << ",sigma2Err," << m.e3 << ",sigma3," << m.p6
          << ",sigma3Err," << m.e4;
    }
    txt << "\n";
  };

  txt << "Fits\n";
  writeModel("0tag-gauss", f0.gauss);
  writeModel("0tag-voigt", f0.voigt);
  writeModel("0tag-dgauss", f0.dgauss);
  writeModel("0tag-tgauss", f0.tgauss);
  writeModel("0tag-best", f0.best);
  writeModel("1tag-gauss", f1.gauss);
  writeModel("1tag-voigt", f1.voigt);
  writeModel("1tag-dgauss", f1.dgauss);
  writeModel("1tag-tgauss", f1.tgauss);
  writeModel("1tag-best", f1.best);
  writeModel("2tag-gauss", f2.gauss);
  writeModel("2tag-voigt", f2.voigt);
  writeModel("2tag-dgauss", f2.dgauss);
  writeModel("2tag-tgauss", f2.tgauss);
  writeModel("2tag-best", f2.best);
  txt.close();

  std::ofstream bestCsv(outDirStr + "/" + prefixStr + "_best_fit_params.csv");
  bestCsv << "category,best_model,chi2ndf,mean,mean_err,sigma1,sigma1_err,sigma2,sigma2_err,sigma3,sigma3_err,gamma,gamma_err,fit_min,fit_max,tag_wp,frac1,frac2,frac3\n";
  auto writeBest = [&](const std::string& cat, const ModelFit& m) {
    double sigma2 = 0.0, sigma2err = 0.0, sigma3 = 0.0, sigma3err = 0.0, gamma = 0.0, gammaErr = 0.0;
    double frac1 = 0.0, frac2 = 0.0, frac3 = 0.0;
    if (m.name == "Voigt") {
      gamma = m.p3;
      gammaErr = m.e3;
    }
    if (m.name == "DoubleGauss") {
      sigma2 = m.p4;
      sigma2err = m.e3;
    }
    if (m.name == "TripleGauss") {
      sigma2 = m.p4;
      sigma2err = m.e3;
      sigma3 = m.p6;
      sigma3err = m.e4;
      const double i1 = m.p0 * m.p2;
      const double i2 = m.p3 * m.p4;
      const double i3 = m.p5 * m.p6;
      const double it = i1 + i2 + i3;
      if (it > 0.0) {
        frac1 = i1 / it;
        frac2 = i2 / it;
        frac3 = i3 / it;
      }
    }
    bestCsv << cat << "," << m.name << "," << m.chi2ndf << "," << m.p1 << "," << m.e1 << "," << m.p2 << ","
            << m.e2 << "," << sigma2 << "," << sigma2err << "," << sigma3 << "," << sigma3err << ","
            << gamma << "," << gammaErr << "," << kFitMin << "," << kFitMax << ",RecoPIDKaon>=2,"
            << frac1 << "," << frac2 << "," << frac3 << "\n";
  };
  writeBest("0tag", f0.best);
  writeBest("1tag", f1.best);
  writeBest("2tag", f2.best);
  bestCsv.close();

  TFile hout((outDirStr + "/" + prefixStr + "/" + prefixStr + "_qc_hists.root").c_str(), "RECREATE");
  hMass0.Write();
  hMass1.Write();
  hMass2.Write();
  hMassAll.Write();
  hNTag.Write();
  hPhiP.Write();
  hPhiPt.Write();
  hPhiCos.Write();
  hKP.Write();
  hKPt.Write();
  hKCos.Write();
  hAng1All.Write();
  hAng2All.Write();
  hAng1Pass.Write();
  hAng2Pass.Write();
  hout.Close();

  std::cout << "Created report PDF: " << pdfPath << std::endl;
}
