void StyleHistogram(TH2D *H, const char *ZTitle)
{
   H->SetStats(0);
   H->SetMinimum(0.01);
   H->SetMaximum(1.0);
   H->SetContour(255);
   H->GetXaxis()->SetTitle("cos#theta");
   H->GetYaxis()->SetTitle("p (GeV)");
   H->GetZaxis()->SetTitle(ZTitle);
   H->GetYaxis()->SetRangeUser(0.15, 100.0);
   H->GetXaxis()->CenterTitle();
   H->GetYaxis()->CenterTitle();
   H->GetZaxis()->CenterTitle();
}

void DrawOne(TFile &File, const char *DenominatorName, const char *NumeratorName,
   const char *HistogramName, const char *Title, const char *OutputBase)
{
   TH2D *Denominator = (TH2D *)File.Get(DenominatorName);
   TH2D *Numerator = (TH2D *)File.Get(NumeratorName);

   TH2D *Efficiency = (TH2D *)Numerator->Clone(HistogramName);
   Efficiency->SetDirectory(nullptr);
   Efficiency->Divide(Numerator, Denominator, 1.0, 1.0, "B");
   Efficiency->SetTitle(Title);
   StyleHistogram(Efficiency, "Matching efficiency");

   TCanvas Canvas("Canvas", "", 900, 700);
   Canvas.SetLogy();
   Canvas.SetLogz();
   Canvas.SetRightMargin(0.16);
   Canvas.SetLeftMargin(0.11);
   Canvas.SetBottomMargin(0.12);

   Efficiency->Draw("colz");
   Canvas.SaveAs(Form("%s.pdf", OutputBase));
   Canvas.SaveAs(Form("%s.png", OutputBase));

   delete Efficiency;
}

void make_matching_efficiency_plots()
{
   gROOT->SetBatch(kTRUE);
   gStyle->SetOptStat(0);
   gStyle->SetPalette(kBird);
   TGaxis::SetMaxDigits(3);

   TFile File("20260304_MC_Merged_Matched_EfficiencyVariableBinning_FakeRate.root");

   DrawOne(File, "HGenPion", "HGenPionMatched",
      "HEffPionMatch",
      "Generated #pi #rightarrow matched reco track",
      "matching_efficiency_pion");
   DrawOne(File, "HGenKaon", "HGenKaonMatched",
      "HEffKaonMatch",
      "Generated K #rightarrow matched reco track",
      "matching_efficiency_kaon");
   DrawOne(File, "HGenProton", "HGenProtonMatched",
      "HEffProtonMatch",
      "Generated p #rightarrow matched reco track",
      "matching_efficiency_proton");
}
