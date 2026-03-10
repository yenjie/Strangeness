#include <iostream>
#include <string>

#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"

void StyleHistogram(TH1D* h, int color) {
  h->SetStats(0);
  h->SetLineColor(color);
  h->SetLineWidth(2);
}

void DrawSingleKShort(TH1D* h, const std::string& title, int color, const std::string& outputName) {
  TCanvas c("c_single_ks", "c_single_ks", 800, 600);
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

int PlotKShortWrongAsKStar(std::string input = "KShortWrongAsKStarHistograms.root",
                           std::string outputDir = "PlotsKShortWrongAsKStar") {
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
  file.GetObject("hKShortWrongAsKStarMassAccepted", hAccepted);
  file.GetObject("hKShortWrongAsKStarMassKaonTag", hKaonTag);
  file.GetObject("hKShortWrongAsKStarMassKaonPionTag", hKaonPionTag);
  if (hAccepted == nullptr || hKaonTag == nullptr || hKaonPionTag == nullptr) {
    std::cerr << "Missing histograms in " << input << std::endl;
    return 1;
  }

  DrawSingleKShort(hAccepted,
                   "K^{0}_{S} #rightarrow #pi^{+}#pi^{-} treated as K#pi: accepted in K^{*} window",
                   kBlack,
                   outputDir + "/kstar_kshort_wrong_treatment_accepted.pdf");
  DrawSingleKShort(hKaonTag,
                   "K^{0}_{S} #rightarrow #pi^{+}#pi^{-} treated as K#pi: kaon tag, no pion tag",
                   kBlue + 1,
                   outputDir + "/kstar_kshort_wrong_treatment_no_pion_tag.pdf");
  DrawSingleKShort(hKaonPionTag,
                   "K^{0}_{S} #rightarrow #pi^{+}#pi^{-} treated as K#pi: kaon+pion tag",
                   kRed + 1,
                   outputDir + "/kstar_kshort_wrong_treatment_with_pion_tag.pdf");

  return 0;
}
