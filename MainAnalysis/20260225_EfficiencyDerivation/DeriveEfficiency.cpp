#include <cmath>
#include <iostream>
using namespace std;

#include "TH2D.h"
#include "TFile.h"

#include "StrangenessMessenger.h"
#include "CommandLine.h"
#include "ProgressBar.h"

int main(int argc, char *argv[])
{
   CommandLine CL(argc, argv);

   string InputFileName = CL.Get("Input", "Trees/merged_mc_v2.root");
   string OutputFileName = CL.Get("Output", "Efficiency.root");
   double Fraction      = CL.GetDouble("Fraction", 1.00);

   int NBinsX = 51;
   int NBinsY = 31;
   double BinsX[] = {-1, -0.94, -0.91, -0.82, -0.70, -0.675, -0.65, -0.625, -0.575, -0.55, -0.525, -0.5, -0.475, -0.45, -0.4, -0.375, -0.35, -0.325, -0.3, -0.275, -0.25, -0.225, -0.2, -0.175, -0.15, -0.05, 0.05, 0.15, 0.175, 0.2, 0.225, 0.25, 0.275, 0.3, 0.325, 0.35, 0.375, 0.4, 0.425, 0.45, 0.475, 0.5, 0.525, 0.55, 0.575, 0.6, 0.625, 0.65, 0.675, 0.70, 0.82, 0.91, 0.94, 1.00};
   double BinsY[] = {0, 0.15, 0.25, 0.35, 0.4, 0.5, 0.6, 0.718, 0.8, 0.9, 1.00, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.8, 2.0, 2.3, 2.5, 2.6, 2.8, 2.9, 3.0, 3.5, 4.0, 5.0, 6.0, 7.0, 8.0, 100};

   TFile InputFile(InputFileName.c_str());
   TFile OutputFile(OutputFileName.c_str(), "RECREATE");
   
   TH2D HGenPion("HGenPion", ";;", NBinsX, BinsX, NBinsY, BinsY);
   TH2D HGenPionMatched("HGenPionMatched", ";;", NBinsX, BinsX, NBinsY, BinsY);
   TH2D HGenPionMatchedPionTagged("HGenPionMatchedPionTagged", ";;", NBinsX, BinsX, NBinsY, BinsY);
   TH2D HGenPionMatchedKaonTagged("HGenPionMatchedKaonTagged", ";;", NBinsX, BinsX, NBinsY, BinsY);
   TH2D HGenPionMatchedProtonTagged("HGenPionMatchedProtonTagged", ";;", NBinsX, BinsX, NBinsY, BinsY);
   TH2D HGenKaon("HGenKaon", ";;", NBinsX, BinsX, NBinsY, BinsY);
   TH2D HGenKaonMatched("HGenKaonMatched", ";;", NBinsX, BinsX, NBinsY, BinsY);
   TH2D HGenKaonMatchedPionTagged("HGenKaonMatchedPionTagged", ";;", NBinsX, BinsX, NBinsY, BinsY);
   TH2D HGenKaonMatchedKaonTagged("HGenKaonMatchedKaonTagged", ";;", NBinsX, BinsX, NBinsY, BinsY);
   TH2D HGenKaonMatchedProtonTagged("HGenKaonMatchedProtonTagged", ";;", NBinsX, BinsX, NBinsY, BinsY);
   TH2D HGenProton("HGenProton", ";;", NBinsX, BinsX, NBinsY, BinsY);
   TH2D HGenProtonMatched("HGenProtonMatched", ";;", NBinsX, BinsX, NBinsY, BinsY);
   TH2D HGenProtonMatchedPionTagged("HGenProtonMatchedPionTagged", ";;", NBinsX, BinsX, NBinsY, BinsY);
   TH2D HGenProtonMatchedKaonTagged("HGenProtonMatchedKaonTagged", ";;", NBinsX, BinsX, NBinsY, BinsY);
   TH2D HGenProtonMatchedProtonTagged("HGenProtonMatchedProtonTagged", ";;", NBinsX, BinsX, NBinsY, BinsY);

   TH2D HRecoPion("HRecoPion", ";;", NBinsX, BinsX, NBinsY, BinsY);
   TH2D HRecoPionMatched("HRecoPionMatched", ";;", NBinsX, BinsX, NBinsY, BinsY);
   TH2D HRecoKaon("HRecoKaon", ";;", NBinsX, BinsX, NBinsY, BinsY);
   TH2D HRecoKaonMatched("HRecoKaonMatched", ";;", NBinsX, BinsX, NBinsY, BinsY);
   TH2D HRecoProton("HRecoProton", ";;", NBinsX, BinsX, NBinsY, BinsY);
   TH2D HRecoProtonMatched("HRecoProtonMatched", ";;", NBinsX, BinsX, NBinsY, BinsY);

   StrangenessTreeMessenger M(InputFile);

   int EntryCount = M.GetEntries() * Fraction;
   for(int iE = 0; iE < EntryCount; iE++)
   {
      M.GetEntry(iE);

      if(M.PassAll == false)
         continue;

      vector<bool> RecoMatched(M.NReco, false);

      for(int iG = 0; iG < M.NGen; iG++)
      {
         if(M.GenStatus[iG] != 1)
            continue;
         if(M.GenID[iG] != 211 && M.GenID[iG] != -211
            && M.GenID[iG] != 321 && M.GenID[iG] != -321
            && M.GenID[iG] != 2212 && M.GenID[iG] != -2212)
            continue;

         double P = sqrt(M.GenPx[iG] * M.GenPx[iG] + M.GenPy[iG] * M.GenPy[iG] + M.GenPz[iG] * M.GenPz[iG]);
         double CosTheta = M.GenPz[iG] / P;

         bool Matched = M.GenMatchAngle[iG] < 0.01;
         if(Matched && M.RecoGoodTrack[M.GenMatchIndex[iG]] == false)
            Matched = false;

         if(Matched == true)
            RecoMatched[M.GenMatchIndex[iG]] = true;

         if(M.GenID[iG] == 211 || M.GenID[iG] == -211)
         {
            HGenPion.Fill(CosTheta, P);
            if(Matched == true)
               HGenPionMatched.Fill(CosTheta, P);
            if(Matched == true && M.RecoPIDPion[M.GenMatchIndex[iG]] >= 2)
               HGenPionMatchedPionTagged.Fill(CosTheta, P);
            if(Matched == true && M.RecoPIDKaon[M.GenMatchIndex[iG]] >= 2)
               HGenPionMatchedKaonTagged.Fill(CosTheta, P);
            if(Matched == true && M.RecoPIDProton[M.GenMatchIndex[iG]] >= 2)
               HGenPionMatchedProtonTagged.Fill(CosTheta, P);
         }
         if(M.GenID[iG] == 321 || M.GenID[iG] == -321)
         {
            HGenKaon.Fill(CosTheta, P);
            if(Matched == true)
               HGenKaonMatched.Fill(CosTheta, P);
            if(Matched == true && M.RecoPIDPion[M.GenMatchIndex[iG]] >= 2)
               HGenKaonMatchedPionTagged.Fill(CosTheta, P);
            if(Matched == true && M.RecoPIDKaon[M.GenMatchIndex[iG]] >= 2)
               HGenKaonMatchedKaonTagged.Fill(CosTheta, P);
            if(Matched == true && M.RecoPIDProton[M.GenMatchIndex[iG]] >= 2)
               HGenKaonMatchedProtonTagged.Fill(CosTheta, P);
         }
         if(M.GenID[iG] == 2212 || M.GenID[iG] == -2212)
         {
            HGenProton.Fill(CosTheta, P);
            if(Matched == true)
               HGenProtonMatched.Fill(CosTheta, P);
            if(Matched == true && M.RecoPIDPion[M.GenMatchIndex[iG]] >= 2)
               HGenProtonMatchedPionTagged.Fill(CosTheta, P);
            if(Matched == true && M.RecoPIDKaon[M.GenMatchIndex[iG]] >= 2)
               HGenProtonMatchedKaonTagged.Fill(CosTheta, P);
            if(Matched == true && M.RecoPIDProton[M.GenMatchIndex[iG]] >= 2)
               HGenProtonMatchedProtonTagged.Fill(CosTheta, P);
         }
      }

      for(int iR = 0; iR < M.NReco; iR++)
      {
         if(M.RecoGoodTrack[iR] != 1)
            continue;

         double P = sqrt(M.RecoPx[iR] * M.RecoPx[iR] + M.RecoPy[iR] * M.RecoPy[iR] + M.RecoPz[iR] * M.RecoPz[iR]);
         double CosTheta = M.RecoPz[iR] / P;

         if(M.RecoPIDPion[iR] >= 2)
         {
            HRecoPion.Fill(CosTheta, P);
            if(RecoMatched[iR] == true)
               HRecoPionMatched.Fill(CosTheta, P);
         }
         if(M.RecoPIDKaon[iR] >= 2)
         {
            HRecoKaon.Fill(CosTheta, P);
            if(RecoMatched[iR] == true)
               HRecoKaonMatched.Fill(CosTheta, P);
         }
         if(M.RecoPIDProton[iR] >= 2)
         {
            HRecoProton.Fill(CosTheta, P);
            if(RecoMatched[iR] == true)
               HRecoProtonMatched.Fill(CosTheta, P);
         }
      }
   }

   HGenPion.Write();
   HGenPionMatched.Write();
   HGenPionMatchedPionTagged.Write();
   HGenPionMatchedKaonTagged.Write();
   HGenPionMatchedProtonTagged.Write();

   HGenKaon.Write();
   HGenKaonMatched.Write();
   HGenKaonMatchedPionTagged.Write();
   HGenKaonMatchedKaonTagged.Write();
   HGenKaonMatchedProtonTagged.Write();
   
   HGenProton.Write();
   HGenProtonMatched.Write();
   HGenProtonMatchedPionTagged.Write();
   HGenProtonMatchedKaonTagged.Write();
   HGenProtonMatchedProtonTagged.Write();

   HRecoPion.Write();
   HRecoPionMatched.Write();
   HRecoKaon.Write();
   HRecoKaonMatched.Write();
   HRecoProton.Write();
   HRecoProtonMatched.Write();

   TH2D *HGenPionEfficiency = (TH2D *)HGenPionMatched.Clone("HGenPionEfficiency");
   TH2D *HGenPionPionTaggedEfficiency = (TH2D *)HGenPionMatchedPionTagged.Clone("HGenPionEfficiencyPionTagged");
   TH2D *HGenPionKaonTaggedEfficiency = (TH2D *)HGenPionMatchedKaonTagged.Clone("HGenPionEfficiencyKaonTagged");
   TH2D *HGenPionProtonTaggedEfficiency = (TH2D *)HGenPionMatchedProtonTagged.Clone("HGenPionEfficiencyProtonTagged");
   TH2D *HGenKaonEfficiency = (TH2D *)HGenKaonMatched.Clone("HGenKaonEfficiency");
   TH2D *HGenKaonPionTaggedEfficiency = (TH2D *)HGenKaonMatchedPionTagged.Clone("HGenKaonEfficiencyPionTagged");
   TH2D *HGenKaonKaonTaggedEfficiency = (TH2D *)HGenKaonMatchedKaonTagged.Clone("HGenKaonEfficiencyKaonTagged");
   TH2D *HGenKaonProtonTaggedEfficiency = (TH2D *)HGenKaonMatchedProtonTagged.Clone("HGenKaonEfficiencyProtonTagged");
   TH2D *HGenProtonEfficiency = (TH2D *)HGenProtonMatched.Clone("HGenProtonEfficiency");
   TH2D *HGenProtonPionTaggedEfficiency = (TH2D *)HGenProtonMatchedPionTagged.Clone("HGenProtonEfficiencyPionTagged");
   TH2D *HGenProtonKaonTaggedEfficiency = (TH2D *)HGenProtonMatchedKaonTagged.Clone("HGenProtonEfficiencyKaonTagged");
   TH2D *HGenProtonProtonTaggedEfficiency = (TH2D *)HGenProtonMatchedProtonTagged.Clone("HGenProtonEfficiencyProtonTagged");

   TH2D *HRecoPionEfficiency = (TH2D *)HRecoPionMatched.Clone("HRecoPionEfficiency");
   TH2D *HRecoKaonEfficiency = (TH2D *)HRecoKaonMatched.Clone("HRecoKaonEfficiency");
   TH2D *HRecoProtonEfficiency = (TH2D *)HRecoProtonMatched.Clone("HRecoProtonEfficiency");
   
   HGenPionEfficiency->Divide(&HGenPion);
   HGenPionPionTaggedEfficiency->Divide(&HGenPionMatched);
   HGenPionKaonTaggedEfficiency->Divide(&HGenPionMatched);
   HGenPionProtonTaggedEfficiency->Divide(&HGenPionMatched);
   HGenKaonEfficiency->Divide(&HGenKaon);
   HGenKaonPionTaggedEfficiency->Divide(&HGenKaonMatched);
   HGenKaonKaonTaggedEfficiency->Divide(&HGenKaonMatched);
   HGenKaonProtonTaggedEfficiency->Divide(&HGenKaonMatched);
   HGenProtonEfficiency->Divide(&HGenProton);
   HGenProtonPionTaggedEfficiency->Divide(&HGenProtonMatched);
   HGenProtonKaonTaggedEfficiency->Divide(&HGenProtonMatched);
   HGenProtonProtonTaggedEfficiency->Divide(&HGenProtonMatched);

   HRecoPionEfficiency->Divide(&HRecoPion);
   HRecoKaonEfficiency->Divide(&HRecoKaon);
   HRecoProtonEfficiency->Divide(&HRecoProton);

   HGenPionEfficiency->Write();
   HGenPionPionTaggedEfficiency->Write();
   HGenPionKaonTaggedEfficiency->Write();
   HGenPionProtonTaggedEfficiency->Write();
   HGenKaonEfficiency->Write();
   HGenKaonPionTaggedEfficiency->Write();
   HGenKaonKaonTaggedEfficiency->Write();
   HGenKaonProtonTaggedEfficiency->Write();
   HGenProtonEfficiency->Write();
   HGenProtonPionTaggedEfficiency->Write();
   HGenProtonKaonTaggedEfficiency->Write();
   HGenProtonProtonTaggedEfficiency->Write();

   HRecoPionEfficiency->Write();
   HRecoKaonEfficiency->Write();
   HRecoProtonEfficiency->Write();

   OutputFile.Close();
   InputFile.Close();

   return 0;
}


