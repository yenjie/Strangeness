#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TNamed.h"
#include "TParameter.h"

#include "StrangenessMessenger.h"
#include "TruthCountingPolicy.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

namespace
{
constexpr int kMaxNchTag = 60;
constexpr int kNBinsNch = kMaxNchTag / 4 + 1;
constexpr double kNchMin = -0.5;
constexpr double kNchMax = kMaxNchTag + 0.5;
constexpr double kPtMin = 0.4;
constexpr double kPtMax = 5.0;
constexpr int kNPtBins = 12;
constexpr double kNtagPtMin = 0.2;
constexpr double kGenMatchAngleMax = 0.01;

std::vector<double> BuildPtEdges()
{
   std::vector<double> edges(kNPtBins + 1, 0.0);
   const double dpt = (kPtMax - kPtMin) / kNPtBins;
   for (int i = 0; i <= kNPtBins; ++i)
      edges[i] = kPtMin + dpt * i;
   return edges;
}

bool ComputeEta(double px, double py, double pz, double &eta)
{
   const double pt = std::sqrt(px * px + py * py);
   if (pt <= 0.0)
      return false;
   eta = std::asinh(pz / pt);
   return std::isfinite(eta);
}

int ClassifyExclusiveTag(const StrangenessTreeMessenger &M, int i)
{
   const int kTag = static_cast<int>(M.RecoPIDKaon[i]);
   const int piTag = static_cast<int>(M.RecoPIDPion[i]);
   const int pTag = static_cast<int>(M.RecoPIDProton[i]);
   const bool passKaonTag = (kTag >= 2);
   const bool passPionTag = (piTag >= 2);
   const bool passProtonTag = (pTag >= 2);
   const bool passTag = (passKaonTag || passPionTag || passProtonTag);
   if (!passTag)
      return 3;

   const int best = std::max(kTag, std::max(piTag, pTag));
   int obsCat = 0; // K > pi > p legacy tie rule
   if (piTag > kTag && piTag >= pTag)
      obsCat = 1;
   if (pTag > kTag && pTag > piTag)
      obsCat = 2;
   if (best < 2)
      obsCat = 3;
   return obsCat;
}

bool PassRecoAcceptedTrack(const StrangenessTreeMessenger &M, int i)
{
   if (M.RecoGoodTrack[i] != 1)
      return false;
   if (M.RecoCharge[i] == 0.0)
      return false;
   const double px = M.RecoPx[i];
   const double py = M.RecoPy[i];
   const double pz = M.RecoPz[i];
   if (!TruthCountingPolicy::PassPIDFiducialFromMom(px, py, pz, true, 0.15, 0.675))
      return false;
   const double pt = std::sqrt(px * px + py * py);
   if (pt < kPtMin || pt >= kPtMax)
      return false;
   return true;
}

bool PassTruthSelectedParticle(const StrangenessTreeMessenger &M, int i)
{
   if (M.GenStatus[i] != 1)
      return false;
   const long long apdg = std::llabs(M.GenID[i]);
   if (apdg != 211 && apdg != 321 && apdg != 2212)
      return false;
   const double px = M.GenPx[i];
   const double py = M.GenPy[i];
   const double pz = M.GenPz[i];
   const double pt = std::sqrt(px * px + py * py);
   if (pt < kPtMin || pt >= kPtMax)
      return false;
   if (!TruthCountingPolicy::PassPIDFiducialFromMom(px, py, pz, true, 0.15, 0.675))
      return false;
   return true;
}

void AddWeight(TH2D *h, double x, double y, double w)
{
   h->Fill(x, y, w);
}

void AddWeight(TH3D *h, double x, double y, double z, double w)
{
   h->Fill(x, y, z, w);
}
}

int main(int argc, char *argv[])
{
   if (argc < 4)
   {
      std::cerr << "Usage: " << argv[0]
                << " <input.root> <output.root> <mode:data|mc>\n";
      return 1;
   }

   const std::string inputPath = argv[1];
   const std::string outputPath = argv[2];
   const std::string mode = argv[3];
   const bool isMC = (mode == "mc" || mode == "MC");

   TFile inputFile(inputPath.c_str(), "READ");
   if (inputFile.IsZombie())
   {
      std::cerr << "Cannot open input file: " << inputPath << "\n";
      return 1;
   }

   StrangenessTreeMessenger M(inputFile, "Tree");
   if (M.Tree == nullptr)
   {
      std::cerr << "Missing Tree in " << inputPath << "\n";
      return 1;
   }

   const bool hasRecoEff = (M.Tree->GetBranch("RecoEfficiencyK") != nullptr &&
                            M.Tree->GetBranch("RecoEfficiencyPi") != nullptr &&
                            M.Tree->GetBranch("RecoEfficiencyP") != nullptr);

   const std::vector<double> ptEdges = BuildPtEdges();
   const double *ptEdgeArray = &(ptEdges[0]);

   TFile outputFile(outputPath.c_str(), "RECREATE");
   if (outputFile.IsZombie())
   {
      std::cerr << "Cannot create output file: " << outputPath << "\n";
      return 1;
   }

   TH2D hTagKFakeCorrReco("hTagKFakeCorrRecoDNdEta", "Fake-corrected K-tag reco yield;dN_{ch}/d#eta (reco, |#eta|<0.5);p_{T} (GeV/c)",
                          kNBinsNch, kNchMin, kNchMax, kNPtBins, ptEdgeArray);
   TH2D hTagPiFakeCorrReco("hTagPiFakeCorrRecoDNdEta", "Fake-corrected #pi-tag reco yield;dN_{ch}/d#eta (reco, |#eta|<0.5);p_{T} (GeV/c)",
                           kNBinsNch, kNchMin, kNchMax, kNPtBins, ptEdgeArray);
   TH2D hTagPFakeCorrReco("hTagPFakeCorrRecoDNdEta", "Fake-corrected p-tag reco yield;dN_{ch}/d#eta (reco, |#eta|<0.5);p_{T} (GeV/c)",
                          kNBinsNch, kNchMin, kNchMax, kNPtBins, ptEdgeArray);
   TH1D hRecoCounts("hRecoCountsDNdEta", "Reco dN_{ch}/d#eta counts;dN_{ch}/d#eta (reco, |#eta|<0.5);Events", kNBinsNch, kNchMin, kNchMax);
   hTagKFakeCorrReco.Sumw2();
   hTagPiFakeCorrReco.Sumw2();
   hTagPFakeCorrReco.Sumw2();
   hRecoCounts.Sumw2();

   TH3D hRespTagKFakeCorr("hRespTagKFakeCorrDNdEta", "Fake-corrected K-tag activity response;dN_{ch}/d#eta (true, |#eta|<0.5);dN_{ch}/d#eta (reco, |#eta|<0.5);p_{T} (GeV/c)",
                          kNBinsNch, kNchMin, kNchMax, kNBinsNch, kNchMin, kNchMax, kNPtBins, kPtMin, kPtMax);
   TH3D hRespTagPiFakeCorr("hRespTagPiFakeCorrDNdEta", "Fake-corrected #pi-tag activity response;dN_{ch}/d#eta (true, |#eta|<0.5);dN_{ch}/d#eta (reco, |#eta|<0.5);p_{T} (GeV/c)",
                           kNBinsNch, kNchMin, kNchMax, kNBinsNch, kNchMin, kNchMax, kNPtBins, kPtMin, kPtMax);
   TH3D hRespTagPFakeCorr("hRespTagPFakeCorrDNdEta", "Fake-corrected p-tag activity response;dN_{ch}/d#eta (true, |#eta|<0.5);dN_{ch}/d#eta (reco, |#eta|<0.5);p_{T} (GeV/c)",
                          kNBinsNch, kNchMin, kNchMax, kNBinsNch, kNchMin, kNchMax, kNPtBins, kPtMin, kPtMax);
   TH2D hGenAllK("hGenAllKDNdEtaPt", "All truth K;dN_{ch}/d#eta (true, |#eta|<0.5);p_{T}^{truth} (GeV/c)",
                 kNBinsNch, kNchMin, kNchMax, kNPtBins, ptEdgeArray);
   TH2D hGenAllPi("hGenAllPiDNdEtaPt", "All truth #pi;dN_{ch}/d#eta (true, |#eta|<0.5);p_{T}^{truth} (GeV/c)",
                  kNBinsNch, kNchMin, kNchMax, kNPtBins, ptEdgeArray);
   TH2D hGenAllP("hGenAllPDNdEtaPt", "All truth p;dN_{ch}/d#eta (true, |#eta|<0.5);p_{T}^{truth} (GeV/c)",
                 kNBinsNch, kNchMin, kNchMax, kNPtBins, ptEdgeArray);
   TH2D hMatchedDenK("hMatchedDenKDNdEtaPt", "Matched accepted truth K;dN_{ch}/d#eta (true, |#eta|<0.5);p_{T}^{truth} (GeV/c)",
                     kNBinsNch, kNchMin, kNchMax, kNPtBins, ptEdgeArray);
   TH2D hMatchedDenPi("hMatchedDenPiDNdEtaPt", "Matched accepted truth #pi;dN_{ch}/d#eta (true, |#eta|<0.5);p_{T}^{truth} (GeV/c)",
                      kNBinsNch, kNchMin, kNchMax, kNPtBins, ptEdgeArray);
   TH2D hMatchedDenP("hMatchedDenPDNdEtaPt", "Matched accepted truth p;dN_{ch}/d#eta (true, |#eta|<0.5);p_{T}^{truth} (GeV/c)",
                     kNBinsNch, kNchMin, kNchMax, kNPtBins, ptEdgeArray);
   TH2D hTrueKTagAsK("hTrueKTagAsKDNdEtaPt", "True K tagged as K;dN_{ch}/d#eta (true, |#eta|<0.5);p_{T}^{truth} (GeV/c)",
                     kNBinsNch, kNchMin, kNchMax, kNPtBins, ptEdgeArray);
   TH2D hTrueKTagAsPi("hTrueKTagAsPiDNdEtaPt", "True K tagged as #pi;dN_{ch}/d#eta (true, |#eta|<0.5);p_{T}^{truth} (GeV/c)",
                      kNBinsNch, kNchMin, kNchMax, kNPtBins, ptEdgeArray);
   TH2D hTrueKTagAsP("hTrueKTagAsPDNdEtaPt", "True K tagged as p;dN_{ch}/d#eta (true, |#eta|<0.5);p_{T}^{truth} (GeV/c)",
                     kNBinsNch, kNchMin, kNchMax, kNPtBins, ptEdgeArray);
   TH2D hTruePiTagAsK("hTruePiTagAsKDNdEtaPt", "True #pi tagged as K;dN_{ch}/d#eta (true, |#eta|<0.5);p_{T}^{truth} (GeV/c)",
                      kNBinsNch, kNchMin, kNchMax, kNPtBins, ptEdgeArray);
   TH2D hTruePiTagAsPi("hTruePiTagAsPiDNdEtaPt", "True #pi tagged as #pi;dN_{ch}/d#eta (true, |#eta|<0.5);p_{T}^{truth} (GeV/c)",
                       kNBinsNch, kNchMin, kNchMax, kNPtBins, ptEdgeArray);
   TH2D hTruePiTagAsP("hTruePiTagAsPDNdEtaPt", "True #pi tagged as p;dN_{ch}/d#eta (true, |#eta|<0.5);p_{T}^{truth} (GeV/c)",
                      kNBinsNch, kNchMin, kNchMax, kNPtBins, ptEdgeArray);
   TH2D hTruePTagAsK("hTruePTagAsKDNdEtaPt", "True p tagged as K;dN_{ch}/d#eta (true, |#eta|<0.5);p_{T}^{truth} (GeV/c)",
                     kNBinsNch, kNchMin, kNchMax, kNPtBins, ptEdgeArray);
   TH2D hTruePTagAsPi("hTruePTagAsPiDNdEtaPt", "True p tagged as #pi;dN_{ch}/d#eta (true, |#eta|<0.5);p_{T}^{truth} (GeV/c)",
                      kNBinsNch, kNchMin, kNchMax, kNPtBins, ptEdgeArray);
   TH2D hTruePTagAsP("hTruePTagAsPDNdEtaPt", "True p tagged as p;dN_{ch}/d#eta (true, |#eta|<0.5);p_{T}^{truth} (GeV/c)",
                     kNBinsNch, kNchMin, kNchMax, kNPtBins, ptEdgeArray);
   hRespTagKFakeCorr.Sumw2();
   hRespTagPiFakeCorr.Sumw2();
   hRespTagPFakeCorr.Sumw2();
   hGenAllK.Sumw2();
   hGenAllPi.Sumw2();
   hGenAllP.Sumw2();
   hMatchedDenK.Sumw2();
   hMatchedDenPi.Sumw2();
   hMatchedDenP.Sumw2();
   hTrueKTagAsK.Sumw2();
   hTrueKTagAsPi.Sumw2();
   hTrueKTagAsP.Sumw2();
   hTruePiTagAsK.Sumw2();
   hTruePiTagAsPi.Sumw2();
   hTruePiTagAsP.Sumw2();
   hTruePTagAsK.Sumw2();
   hTruePTagAsPi.Sumw2();
   hTruePTagAsP.Sumw2();

   long long nProcessed = 0;
   long long nPassAll = 0;
   const long long nEntries = M.GetEntries();

   for (long long ievt = 0; ievt < nEntries; ++ievt)
   {
      M.GetEntry(ievt);
      ++nProcessed;
      if (M.PassAll != 1)
         continue;
      ++nPassAll;

      int nreco = static_cast<int>(std::min<long long>(M.NReco, STRANGE_MAX_RECO));
      int ngen = static_cast<int>(std::min<long long>(M.NGen, STRANGE_MAX_GEN));

      int nChEta05Reco = 0;
      for (int i = 0; i < nreco; ++i)
      {
         if (M.RecoGoodTrack[i] != 1)
            continue;
         if (M.RecoCharge[i] == 0.0)
            continue;
         double eta = 0.0;
         if (ComputeEta(M.RecoPx[i], M.RecoPy[i], M.RecoPz[i], eta) && std::abs(eta) < 0.5)
            ++nChEta05Reco;
      }
      if (nChEta05Reco > kMaxNchTag)
         nChEta05Reco = kMaxNchTag;

      int nChEta05True = 0;
      if (isMC)
      {
         for (int i = 0; i < ngen; ++i)
         {
            if (M.GenStatus[i] != 1)
               continue;
            if (!TruthCountingPolicy::IsCountedChargedForActivity(M.GenID[i]))
               continue;
            double eta = 0.0;
            if (ComputeEta(M.GenPx[i], M.GenPy[i], M.GenPz[i], eta) && std::abs(eta) < 0.5)
               ++nChEta05True;
         }
         if (nChEta05True > kMaxNchTag)
            nChEta05True = kMaxNchTag;
      }

      hRecoCounts.Fill(nChEta05Reco);

      for (int i = 0; i < nreco; ++i)
      {
         if (!PassRecoAcceptedTrack(M, i))
            continue;

         const int obsCat = ClassifyExclusiveTag(M, i);
         if (obsCat > 2)
            continue;

         const double pt = std::sqrt(M.RecoPx[i] * M.RecoPx[i] + M.RecoPy[i] * M.RecoPy[i]);
         double fakeWeight = 1.0;
         if (hasRecoEff)
         {
            if (obsCat == 0)
               fakeWeight = M.RecoEfficiencyK[i];
            else if (obsCat == 1)
               fakeWeight = M.RecoEfficiencyPi[i];
            else
               fakeWeight = M.RecoEfficiencyP[i];
         }

         if (obsCat == 0)
         {
            AddWeight(&hTagKFakeCorrReco, nChEta05Reco, pt, fakeWeight);
            if (isMC)
               AddWeight(&hRespTagKFakeCorr, nChEta05True, nChEta05Reco, pt, fakeWeight);
         }
         else if (obsCat == 1)
         {
            AddWeight(&hTagPiFakeCorrReco, nChEta05Reco, pt, fakeWeight);
            if (isMC)
               AddWeight(&hRespTagPiFakeCorr, nChEta05True, nChEta05Reco, pt, fakeWeight);
         }
         else if (obsCat == 2)
         {
            AddWeight(&hTagPFakeCorrReco, nChEta05Reco, pt, fakeWeight);
            if (isMC)
               AddWeight(&hRespTagPFakeCorr, nChEta05True, nChEta05Reco, pt, fakeWeight);
         }
      }

      if (!isMC)
         continue;

      for (int i = 0; i < ngen; ++i)
      {
         if (!PassTruthSelectedParticle(M, i))
            continue;

         const long long apdg = std::llabs(M.GenID[i]);
         const double genPt = std::sqrt(M.GenPx[i] * M.GenPx[i] + M.GenPy[i] * M.GenPy[i]);
         if (apdg == 321)
            hGenAllK.Fill(nChEta05True, genPt);
         else if (apdg == 211)
            hGenAllPi.Fill(nChEta05True, genPt);
         else if (apdg == 2212)
            hGenAllP.Fill(nChEta05True, genPt);

         const int matchIndex = static_cast<int>(M.GenMatchIndex[i]);
         if (matchIndex < 0 || matchIndex >= nreco)
            continue;
         if (M.GenMatchAngle[i] >= kGenMatchAngleMax)
            continue;
         if (!PassRecoAcceptedTrack(M, matchIndex))
            continue;

         if (apdg == 321)
            hMatchedDenK.Fill(nChEta05True, genPt);
         else if (apdg == 211)
            hMatchedDenPi.Fill(nChEta05True, genPt);
         else if (apdg == 2212)
            hMatchedDenP.Fill(nChEta05True, genPt);

         const int obsCat = ClassifyExclusiveTag(M, matchIndex);
         if (obsCat == 0)
         {
            if (apdg == 321) hTrueKTagAsK.Fill(nChEta05True, genPt);
            if (apdg == 211) hTruePiTagAsK.Fill(nChEta05True, genPt);
            if (apdg == 2212) hTruePTagAsK.Fill(nChEta05True, genPt);
         }
         else if (obsCat == 1)
         {
            if (apdg == 321) hTrueKTagAsPi.Fill(nChEta05True, genPt);
            if (apdg == 211) hTruePiTagAsPi.Fill(nChEta05True, genPt);
            if (apdg == 2212) hTruePTagAsPi.Fill(nChEta05True, genPt);
         }
         else if (obsCat == 2)
         {
            if (apdg == 321) hTrueKTagAsP.Fill(nChEta05True, genPt);
            if (apdg == 211) hTruePiTagAsP.Fill(nChEta05True, genPt);
            if (apdg == 2212) hTruePTagAsP.Fill(nChEta05True, genPt);
         }
      }
   }

   outputFile.cd();
   hTagKFakeCorrReco.Write();
   hTagPiFakeCorrReco.Write();
   hTagPFakeCorrReco.Write();
   hRecoCounts.Write();
   if (isMC)
   {
      hRespTagKFakeCorr.Write();
      hRespTagPiFakeCorr.Write();
      hRespTagPFakeCorr.Write();
      hGenAllK.Write();
      hGenAllPi.Write();
      hGenAllP.Write();
      hMatchedDenK.Write();
      hMatchedDenPi.Write();
      hMatchedDenP.Write();
      hTrueKTagAsK.Write();
      hTrueKTagAsPi.Write();
      hTrueKTagAsP.Write();
      hTruePiTagAsK.Write();
      hTruePiTagAsPi.Write();
      hTruePiTagAsP.Write();
      hTruePTagAsK.Write();
      hTruePTagAsPi.Write();
      hTruePTagAsP.Write();
   }
   TNamed("crosscheck_method", "Order-swapped dN/deta study: fake-correct observed tags -> activity unfolding -> truth-matched 3x3 inversion -> truth-matched gen-efficiency correction").Write();
   TNamed("crosscheck_scope", "Approximate reordered cross-check within the existing factorized dN/deta framework; not the full Yi Chen 3D (p,|cos(theta)|,activity) procedure").Write();
   TParameter<long long>("nProcessed", nProcessed).Write();
   TParameter<long long>("nPassAll", nPassAll).Write();
   TParameter<int>("hasRecoEfficiencyBranches", hasRecoEff ? 1 : 0).Write();
   outputFile.Close();

   std::cout << "Wrote " << outputPath << "\n";
   std::cout << "Processed entries: " << nProcessed << ", PassAll: " << nPassAll << "\n";
   return 0;
}
