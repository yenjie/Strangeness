void StylePIDHistogram(TH2D *H, const char *ZTitle)
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

TGraph *MakePTThresholdGraph(double PTMin)
{
   const int N = 400;
   TGraph *G = new TGraph;
   for(int i = 0; i < N; i++)
   {
      double CosTheta = -0.999 + (1.998 / (N - 1)) * i;
      double SinTheta = sqrt(1.0 - CosTheta * CosTheta);
      if(SinTheta <= 0)
         continue;
      double P = PTMin / SinTheta;
      G->SetPoint(G->GetN(), CosTheta, P);
   }

   G->SetLineColor(kRed + 1);
   G->SetLineWidth(2);
   G->SetLineStyle(2);
   return G;
}

void DrawPIDOne(TFile &File, const char *DenominatorName, const char *NumeratorName,
   const char *HistogramName, const char *Title, const char *OutputBase)
{
   TH2D *Denominator = (TH2D *)File.Get(DenominatorName);
   TH2D *Numerator = (TH2D *)File.Get(NumeratorName);

   TH2D *Efficiency = (TH2D *)Numerator->Clone(HistogramName);
   Efficiency->SetDirectory(nullptr);
   Efficiency->Divide(Numerator, Denominator, 1.0, 1.0, "B");
   Efficiency->SetTitle(Title);
   StylePIDHistogram(Efficiency, "Tag probability");

   TCanvas Canvas("Canvas", "", 900, 700);
   Canvas.SetLogy();
   Canvas.SetLogz();
   Canvas.SetRightMargin(0.16);
   Canvas.SetLeftMargin(0.11);
   Canvas.SetBottomMargin(0.12);

   Efficiency->Draw("colz");
   TGraph *Threshold = MakePTThresholdGraph(0.4);
   Threshold->Draw("l same");
   Canvas.SaveAs(Form("%s.pdf", OutputBase));
   Canvas.SaveAs(Form("%s.png", OutputBase));

   delete Threshold;
   delete Efficiency;
}

void make_pid_response_plots()
{
   gROOT->SetBatch(kTRUE);
   gStyle->SetOptStat(0);
   gStyle->SetPalette(kBird);
   TGaxis::SetMaxDigits(3);

   TFile File("20260304_MC_Merged_Matched_EfficiencyVariableBinning_FakeRate.root");

   DrawPIDOne(File, "HGenPionMatched", "HGenPionMatchedPionTagged",
      "HPionAsPion", "Matched generated #pi tagged as #pi", "pid_response_pion_as_pion");
   DrawPIDOne(File, "HGenPionMatched", "HGenPionMatchedKaonTagged",
      "HPionAsKaon", "Matched generated #pi tagged as K", "pid_response_pion_as_kaon");
   DrawPIDOne(File, "HGenPionMatched", "HGenPionMatchedProtonTagged",
      "HPionAsProton", "Matched generated #pi tagged as p", "pid_response_pion_as_proton");

   DrawPIDOne(File, "HGenKaonMatched", "HGenKaonMatchedPionTagged",
      "HKaonAsPion", "Matched generated K tagged as #pi", "pid_response_kaon_as_pion");
   DrawPIDOne(File, "HGenKaonMatched", "HGenKaonMatchedKaonTagged",
      "HKaonAsKaon", "Matched generated K tagged as K", "pid_response_kaon_as_kaon");
   DrawPIDOne(File, "HGenKaonMatched", "HGenKaonMatchedProtonTagged",
      "HKaonAsProton", "Matched generated K tagged as p", "pid_response_kaon_as_proton");

   DrawPIDOne(File, "HGenProtonMatched", "HGenProtonMatchedPionTagged",
      "HProtonAsPion", "Matched generated p tagged as #pi", "pid_response_proton_as_pion");
   DrawPIDOne(File, "HGenProtonMatched", "HGenProtonMatchedKaonTagged",
      "HProtonAsKaon", "Matched generated p tagged as K", "pid_response_proton_as_kaon");
   DrawPIDOne(File, "HGenProtonMatched", "HGenProtonMatchedProtonTagged",
      "HProtonAsProton", "Matched generated p tagged as p", "pid_response_proton_as_proton");
}
