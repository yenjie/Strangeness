#include <cmath>
#include <iostream>
using namespace std;

#include "TH1D.h"
#include "TFile.h"

#include "StrangenessMessenger.h"
#include "CommandLine.h"
#include "ProgressBar.h"

int main(int argc, char *argv[])
{
   CommandLine CL(argc, argv);

   string InputFileName = CL.Get("Input", "Trees/merged_mc_v2.2_partial.root");
   string OutputFileName = CL.Get("Output", "EfficiencyClosure.root");
   double Fraction      = CL.GetDouble("Fraction", 1.00);

   TFile InputFile(InputFileName.c_str());
   TFile OutputFile(OutputFileName.c_str(), "RECREATE");
   
   TH1D HGenPion("HGenPion", ";;", 50, 0, 5);
   TH1D HGenKaon("HGenKaon", ";;", 50, 0, 5);
   TH1D HGenProton("HGenProton", ";;", 50, 0, 5);

   TH1D HRecoPion("HRecoPion", ";;", 50, 0, 5);
   TH1D HRecoKaon("HRecoKaon", ";;", 50, 0, 5);
   TH1D HRecoProton("HRecoProton", ";;", 50, 0, 5);
   
   TH1D HRecoPionCorrected("HRecoPionCorrected", ";;", 50, 0, 5);
   TH1D HRecoKaonCorrected("HRecoKaonCorrected", ";;", 50, 0, 5);
   TH1D HRecoProtonCorrected("HRecoProtonCorrected", ";;", 50, 0, 5);

   TH1D HRecoPionMistagAsKaon("HRecoPionMistagAsKaon", ";;", 50, 0, 5);
   
   StrangenessTreeMessenger M(InputFile);

   int EntryCount = M.GetEntries() * Fraction;
   for(int iE = 0; iE < EntryCount; iE++)
   {
      M.GetEntry(iE);

      if(M.PassAll == false)
         continue;

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

         if(abs(CosTheta) < 0.2 || abs(CosTheta) > 0.6)
            continue;

         if(M.GenID[iG] == 211 || M.GenID[iG] == -211)
            HGenPion.Fill(P);
         if(M.GenID[iG] == 321 || M.GenID[iG] == -321)
            HGenKaon.Fill(P);
         if(M.GenID[iG] == 2212 || M.GenID[iG] == -2212)
            HGenProton.Fill(P);
      }

      for(int iR = 0; iR < M.NReco; iR++)
      {
         if(M.RecoGoodTrack[iR] != 1)
            continue;

         double P = sqrt(M.RecoPx[iR] * M.RecoPx[iR] + M.RecoPy[iR] * M.RecoPy[iR] + M.RecoPz[iR] * M.RecoPz[iR]);
         double CosTheta = M.RecoPz[iR] / P;

         if(abs(CosTheta) < 0.2 || abs(CosTheta) > 0.6)
            continue;

         if(M.RecoPIDPion[iR] >= 2)
         {
            HRecoPion.Fill(P);
            HRecoPionCorrected.Fill(P, M.RecoEfficiencyPi[iR] / M.RecoEfficiencyPiAsPi[iR]);
            HRecoPionMistagAsKaon.Fill(P, M.RecoEfficiencyPi[iR] / M.RecoEfficiencyPiAsPi[iR] * M.RecoEfficiencyPiAsK[iR]);
         }
         if(M.RecoPIDKaon[iR] >= 2)
         {
            HRecoKaon.Fill(P);
            HRecoKaonCorrected.Fill(P, M.RecoEfficiencyK[iR] / M.RecoEfficiencyKAsK[iR]);
         }
         if(M.RecoPIDProton[iR] >= 2)
         {
            HRecoProton.Fill(P);
            HRecoProtonCorrected.Fill(P, M.RecoEfficiencyP[iR] / M.RecoEfficiencyPAsP[iR]);
         }
      }  
   }

   HGenPion.Write();
   HGenKaon.Write();
   HGenProton.Write();
   
   HRecoPion.Write();
   HRecoKaon.Write();
   HRecoProton.Write();
   
   HRecoPionCorrected.Write();
   HRecoKaonCorrected.Write();
   HRecoProtonCorrected.Write();

   HRecoPionMistagAsKaon.Write();

   OutputFile.Close();
   InputFile.Close();

   return 0;
}








