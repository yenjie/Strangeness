#include <iostream>
#include <string>

#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"

void StyleHistogram(TH1D* h, int color) {
  h->SetStats(0);
  h->SetLineColor(color);
  h->SetLineWidth(2);
}

void DrawSingle(TH1D* h, const std::string& title, int color, const std::string& outputName) {
  TCanvas c("c_single", "c_single", 800, 600);
  c.SetLeftMargin(0.13);
  c.SetRightMargin(0.04);
  c.SetBottomMargin(0.12);

  TH1D* hc = static_cast<TH1D*>(h->Clone((std::string(h->GetName()) + "_single").c_str()));
  StyleHistogram(hc, color);
  hc->SetTitle(title.c_str());
  hc->GetXaxis()->SetTitle("m(K#pi) [GeV]");
  hc->GetYaxis()->SetTitle("Candidates / bin");
  hc->SetMinimum(0.0);
  hc->Draw("hist");
  c.SaveAs(outputName.c_str());
  delete hc;
}

int PlotPhiWrongAsKStar(std::string input = "PhiWrongAsKStarHistograms.root",
                        std::string outputDir = "PlotsPhiWrongAsKStar") {
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);
  gSystem->mkdir(outputDir.c_str(), true);

  TFile file(input.c_str(), "READ");
  if (file.IsZombie()) {
    std::cerr << "Cannot open " << input << std::endl;
    return 1;
  }

  TH1D* hAccepted = nullptr;
  TH1D* hKaonTag = nullptr;
  TH1D* hKaonPionTag = nullptr;
  file.GetObject("hPhiWrongAsKStarMassAccepted", hAccepted);
  file.GetObject("hPhiWrongAsKStarMassKaonTag", hKaonTag);
  file.GetObject("hPhiWrongAsKStarMassKaonPionTag", hKaonPionTag);
  if (hAccepted == nullptr || hKaonTag == nullptr || hKaonPionTag == nullptr) {
    std::cerr << "Missing histograms in " << input << std::endl;
    return 1;
  }

  TH1D* hK = static_cast<TH1D*>(hKaonTag->Clone("hK_draw"));
  TH1D* hKP = static_cast<TH1D*>(hKaonPionTag->Clone("hKP_draw"));
  StyleHistogram(hK, kBlue + 1);
  StyleHistogram(hKP, kRed + 1);
  TH1D* hAcc = static_cast<TH1D*>(hAccepted->Clone("hAcc_draw"));
  StyleHistogram(hAcc, kBlack);

  DrawSingle(hAcc,
             "#phi #rightarrow K^{+}K^{-} treated as K#pi: accepted in K^{*} window",
             kBlack,
             outputDir + "/kstar_phi_wrong_treatment_accepted.pdf");
  DrawSingle(hK,
             "#phi #rightarrow K^{+}K^{-} treated as K#pi: kaon tag, no pion tag",
             kBlue + 1,
             outputDir + "/kstar_phi_wrong_treatment_no_pion_tag.pdf");
  DrawSingle(hKP,
             "#phi #rightarrow K^{+}K^{-} treated as K#pi: kaon+pion tag",
             kRed + 1,
             outputDir + "/kstar_phi_wrong_treatment_with_pion_tag.pdf");

  TCanvas c("c", "c", 800, 600);
  c.SetLeftMargin(0.13);
  c.SetRightMargin(0.04);
  c.SetBottomMargin(0.12);
  hK->GetXaxis()->SetTitle("m(K#pi) [GeV]");
  hK->GetYaxis()->SetTitle("Candidates / bin");
  hK->SetMinimum(0.0);
  hK->SetMaximum(std::max(hK->GetMaximum(), hKP->GetMaximum()) * 1.18);
  hK->SetTitle("#phi #rightarrow K^{+}K^{-} treated as K#pi in the K^{*} workflow");
  hK->Draw("hist");
  hKP->Draw("hist same");

  TLegend leg2(0.55, 0.62, 0.88, 0.87);
  leg2.SetBorderSize(0);
  leg2.SetFillStyle(0);
  leg2.AddEntry(hK, "RecoPIDKaon #geq 2 on assumed kaon", "l");
  leg2.AddEntry(hKP, "Also require RecoPIDPion #geq 2 on assumed pion", "l");
  leg2.Draw();

  c.SaveAs((outputDir + "/phi_wrong_as_kstar.pdf").c_str());
  delete hK;
  delete hKP;
  delete hAcc;
  return 0;
}
