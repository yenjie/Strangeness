#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TNamed.h"

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
constexpr double kGenMatchAngleMax = 0.01;

std::vector<double> BuildPEdges()
{
   return {0.4, 0.8, 1.2, 1.6, 2.0, 3.0, 5.0, 10.0, 50.0};
}

std::vector<double> BuildAbsCosEdges()
{
   return {0.15, 0.25, 0.35, 0.50, 0.675};
}

bool ComputeEta(double px, double py, double pz, double &eta)
{
   const double pt = std::sqrt(px * px + py * py);
   if (pt <= 0.0)
      return false;
   eta = std::asinh(pz / pt);
   return std::isfinite(eta);
}

double ComputeP(double px, double py, double pz)
{
   return std::sqrt(px * px + py * py + pz * pz);
}

bool ComputeAbsCosTheta(double px, double py, double pz, double &absCos)
{
   const double p = ComputeP(px, py, pz);
   if (p <= 0.0)
      return false;
   absCos = std::abs(pz) / p;
   return std::isfinite(absCos);
}

int FindBin(const std::vector<double> &edges, double value)
{
   if (edges.size() < 2)
      return -1;
   if (value < edges.front() || value >= edges.back())
      return -1;
   for (int i = 0; i + 1 < static_cast<int>(edges.size()); ++i)
      if (value >= edges[i] && value < edges[i + 1])
         return i;
   return -1;
}

int FindActivityBin(double nCh)
{
   if (nCh < kNchMin || nCh >= kNchMax)
      return -1;
   const double width = (kNchMax - kNchMin) / kNBinsNch;
   const int bin = static_cast<int>((nCh - kNchMin) / width);
   if (bin < 0 || bin >= kNBinsNch)
      return -1;
   return bin;
}

int FlatIndex(int actBin, int pBin, int cosBin, int nP, int nCos)
{
   return (actBin * nP + pBin) * nCos + cosBin;
}

int ClassifyExclusiveTag(const StrangenessTreeMessenger &M, int i)
{
   const int kTag = static_cast<int>(M.RecoPIDKaon[i]);
   const int piTag = static_cast<int>(M.RecoPIDPion[i]);
   const int pTag = static_cast<int>(M.RecoPIDProton[i]);
   const bool passKaonTag = (kTag >= 2);
   const bool passPionTag = (piTag >= 2);
   const bool passProtonTag = (pTag >= 2);
   if (!(passKaonTag || passPionTag || passProtonTag))
      return 3;

   const int best = std::max(kTag, std::max(piTag, pTag));
   int obsCat = 0; // K > pi > p
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

double GetFakeWeight(const StrangenessTreeMessenger &M, int recoIndex, int obsCat, bool hasRecoEff)
{
   if (!hasRecoEff || recoIndex < 0)
      return 1.0;
   if (obsCat == 0)
      return M.RecoEfficiencyK[recoIndex];
   if (obsCat == 1)
      return M.RecoEfficiencyPi[recoIndex];
   if (obsCat == 2)
      return M.RecoEfficiencyP[recoIndex];
   return 1.0;
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

   const std::vector<double> pEdges = BuildPEdges();
   const std::vector<double> cosEdges = BuildAbsCosEdges();
   const int nP = static_cast<int>(pEdges.size()) - 1;
   const int nCos = static_cast<int>(cosEdges.size()) - 1;
   const int nFlat = kNBinsNch * nP * nCos;

   TFile outputFile(outputPath.c_str(), "RECREATE");
   if (outputFile.IsZombie())
   {
      std::cerr << "Cannot create output file: " << outputPath << "\n";
      return 1;
   }

   TH1D hRawTagK("hRawTagKRecoFlat", "Fake-corrected K-tag reco counts;Flat reco #mu bin;Weighted counts", nFlat, -0.5, nFlat - 0.5);
   TH1D hRawTagPi("hRawTagPiRecoFlat", "Fake-corrected #pi-tag reco counts;Flat reco #mu bin;Weighted counts", nFlat, -0.5, nFlat - 0.5);
   TH1D hRawTagP("hRawTagPRecoFlat", "Fake-corrected p-tag reco counts;Flat reco #mu bin;Weighted counts", nFlat, -0.5, nFlat - 0.5);
   TH2D hRespTagK("hRespTagKFlat", "K-tag 3D response;Flat true #mu bin;Flat reco #mu bin", nFlat, -0.5, nFlat - 0.5, nFlat, -0.5, nFlat - 0.5);
   TH2D hRespTagPi("hRespTagPiFlat", "#pi-tag 3D response;Flat true #mu bin;Flat reco #mu bin", nFlat, -0.5, nFlat - 0.5, nFlat, -0.5, nFlat - 0.5);
   TH2D hRespTagP("hRespTagPFlat", "p-tag 3D response;Flat true #mu bin;Flat reco #mu bin", nFlat, -0.5, nFlat - 0.5, nFlat, -0.5, nFlat - 0.5);
   TH1D hGenAllK("hGenAllKFlat", "All truth K;Flat true #mu bin;Counts", nFlat, -0.5, nFlat - 0.5);
   TH1D hGenAllPi("hGenAllPiFlat", "All truth #pi;Flat true #mu bin;Counts", nFlat, -0.5, nFlat - 0.5);
   TH1D hGenAllP("hGenAllPFlat", "All truth p;Flat true #mu bin;Counts", nFlat, -0.5, nFlat - 0.5);
   TH1D hMatchedDenK("hMatchedDenKFlat", "Matched truth K;Flat true #mu bin;Counts", nFlat, -0.5, nFlat - 0.5);
   TH1D hMatchedDenPi("hMatchedDenPiFlat", "Matched truth #pi;Flat true #mu bin;Counts", nFlat, -0.5, nFlat - 0.5);
   TH1D hMatchedDenP("hMatchedDenPFlat", "Matched truth p;Flat true #mu bin;Counts", nFlat, -0.5, nFlat - 0.5);
   TH1D hTrueKTagAsK("hTrueKTagAsKFlat", "Weighted true K tagged as K;Flat true #mu bin;Weighted counts", nFlat, -0.5, nFlat - 0.5);
   TH1D hTrueKTagAsPi("hTrueKTagAsPiFlat", "Weighted true K tagged as #pi;Flat true #mu bin;Weighted counts", nFlat, -0.5, nFlat - 0.5);
   TH1D hTrueKTagAsP("hTrueKTagAsPFlat", "Weighted true K tagged as p;Flat true #mu bin;Weighted counts", nFlat, -0.5, nFlat - 0.5);
   TH1D hTruePiTagAsK("hTruePiTagAsKFlat", "Weighted true #pi tagged as K;Flat true #mu bin;Weighted counts", nFlat, -0.5, nFlat - 0.5);
   TH1D hTruePiTagAsPi("hTruePiTagAsPiFlat", "Weighted true #pi tagged as #pi;Flat true #mu bin;Weighted counts", nFlat, -0.5, nFlat - 0.5);
   TH1D hTruePiTagAsP("hTruePiTagAsPFlat", "Weighted true #pi tagged as p;Flat true #mu bin;Weighted counts", nFlat, -0.5, nFlat - 0.5);
   TH1D hTruePTagAsK("hTruePTagAsKFlat", "Weighted true p tagged as K;Flat true #mu bin;Weighted counts", nFlat, -0.5, nFlat - 0.5);
   TH1D hTruePTagAsPi("hTruePTagAsPiFlat", "Weighted true p tagged as #pi;Flat true #mu bin;Weighted counts", nFlat, -0.5, nFlat - 0.5);
   TH1D hTruePTagAsP("hTruePTagAsPFlat", "Weighted true p tagged as p;Flat true #mu bin;Weighted counts", nFlat, -0.5, nFlat - 0.5);
   TH1D hRecoCounts("hRecoCountsDNdEta", "Reco dN_{ch}/d#eta counts;dN_{ch}/d#eta (reco, |#eta|<0.5);Events", kNBinsNch, kNchMin, kNchMax);

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

      const int nreco = static_cast<int>(std::min<long long>(M.NReco, STRANGE_MAX_RECO));
      const int ngen = static_cast<int>(std::min<long long>(M.NGen, STRANGE_MAX_GEN));

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
      const int actRecoBin = FindActivityBin(static_cast<double>(nChEta05Reco));
      if (actRecoBin < 0)
         continue;
      hRecoCounts.Fill(nChEta05Reco);

      int nChEta05True = 0;
      std::vector<int> recoToGen(nreco, -1);
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
         for (int i = 0; i < ngen; ++i)
         {
            const int recoIndex = static_cast<int>(M.GenMatchIndex[i]);
            if (recoIndex < 0 || recoIndex >= nreco)
               continue;
            if (M.GenMatchAngle[i] >= kGenMatchAngleMax)
               continue;
            recoToGen[recoIndex] = i;
         }
      }
      const int actTrueBin = isMC ? FindActivityBin(static_cast<double>(nChEta05True)) : -1;

      for (int i = 0; i < nreco; ++i)
      {
         if (!PassRecoAcceptedTrack(M, i))
            continue;
         const int obsCat = ClassifyExclusiveTag(M, i);
         if (obsCat > 2)
            continue;

         double recoAbsCos = 0.0;
         if (!ComputeAbsCosTheta(M.RecoPx[i], M.RecoPy[i], M.RecoPz[i], recoAbsCos))
            continue;
         const double recoP = ComputeP(M.RecoPx[i], M.RecoPy[i], M.RecoPz[i]);
         const int pRecoBin = FindBin(pEdges, recoP);
         const int cosRecoBin = FindBin(cosEdges, recoAbsCos);
         if (pRecoBin < 0 || cosRecoBin < 0)
            continue;
         const int recoFlat = FlatIndex(actRecoBin, pRecoBin, cosRecoBin, nP, nCos);
         const double fakeWeight = GetFakeWeight(M, i, obsCat, hasRecoEff);

         if (obsCat == 0)
            hRawTagK.Fill(recoFlat, fakeWeight);
         else if (obsCat == 1)
            hRawTagPi.Fill(recoFlat, fakeWeight);
         else
            hRawTagP.Fill(recoFlat, fakeWeight);

         if (!isMC || actTrueBin < 0)
            continue;
         const int genIndex = recoToGen[i];
         if (genIndex < 0)
            continue;
         if (!PassTruthSelectedParticle(M, genIndex))
            continue;

         double genAbsCos = 0.0;
         if (!ComputeAbsCosTheta(M.GenPx[genIndex], M.GenPy[genIndex], M.GenPz[genIndex], genAbsCos))
            continue;
         const double genP = ComputeP(M.GenPx[genIndex], M.GenPy[genIndex], M.GenPz[genIndex]);
         const int pTrueBin = FindBin(pEdges, genP);
         const int cosTrueBin = FindBin(cosEdges, genAbsCos);
         if (pTrueBin < 0 || cosTrueBin < 0)
            continue;
         const int trueFlat = FlatIndex(actTrueBin, pTrueBin, cosTrueBin, nP, nCos);

         if (obsCat == 0)
            hRespTagK.Fill(trueFlat, recoFlat, 1.0);
         else if (obsCat == 1)
            hRespTagPi.Fill(trueFlat, recoFlat, 1.0);
         else
            hRespTagP.Fill(trueFlat, recoFlat, 1.0);
      }

      if (!isMC || actTrueBin < 0)
         continue;

      for (int i = 0; i < ngen; ++i)
      {
         if (!PassTruthSelectedParticle(M, i))
            continue;
         double genAbsCos = 0.0;
         if (!ComputeAbsCosTheta(M.GenPx[i], M.GenPy[i], M.GenPz[i], genAbsCos))
            continue;
         const double genP = ComputeP(M.GenPx[i], M.GenPy[i], M.GenPz[i]);
         const int pTrueBin = FindBin(pEdges, genP);
         const int cosTrueBin = FindBin(cosEdges, genAbsCos);
         if (pTrueBin < 0 || cosTrueBin < 0)
            continue;
         const int trueFlat = FlatIndex(actTrueBin, pTrueBin, cosTrueBin, nP, nCos);

         const long long apdg = std::llabs(M.GenID[i]);
         if (apdg == 321)
            hGenAllK.Fill(trueFlat);
         else if (apdg == 211)
            hGenAllPi.Fill(trueFlat);
         else if (apdg == 2212)
            hGenAllP.Fill(trueFlat);

         const int matchIndex = static_cast<int>(M.GenMatchIndex[i]);
         if (matchIndex < 0 || matchIndex >= nreco)
            continue;
         if (M.GenMatchAngle[i] >= kGenMatchAngleMax)
            continue;
         if (!PassRecoAcceptedTrack(M, matchIndex))
            continue;

         if (apdg == 321)
            hMatchedDenK.Fill(trueFlat);
         else if (apdg == 211)
            hMatchedDenPi.Fill(trueFlat);
         else if (apdg == 2212)
            hMatchedDenP.Fill(trueFlat);

         const int obsCat = ClassifyExclusiveTag(M, matchIndex);
         if (obsCat == 0)
         {
            if (apdg == 321) hTrueKTagAsK.Fill(trueFlat, 1.0);
            if (apdg == 211) hTruePiTagAsK.Fill(trueFlat, 1.0);
            if (apdg == 2212) hTruePTagAsK.Fill(trueFlat, 1.0);
         }
         else if (obsCat == 1)
         {
            if (apdg == 321) hTrueKTagAsPi.Fill(trueFlat, 1.0);
            if (apdg == 211) hTruePiTagAsPi.Fill(trueFlat, 1.0);
            if (apdg == 2212) hTruePTagAsPi.Fill(trueFlat, 1.0);
         }
         else if (obsCat == 2)
         {
            if (apdg == 321) hTrueKTagAsP.Fill(trueFlat, 1.0);
            if (apdg == 211) hTruePiTagAsP.Fill(trueFlat, 1.0);
            if (apdg == 2212) hTruePTagAsP.Fill(trueFlat, 1.0);
         }
      }
   }

   outputFile.cd();
   hRawTagK.Write();
   hRawTagPi.Write();
   hRawTagP.Write();
   hRecoCounts.Write();
   if (isMC)
   {
      hRespTagK.Write();
      hRespTagPi.Write();
      hRespTagP.Write();
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
   TNamed("yi_definition", "Yi-style independent dN/deta inputs: fake-correct observed tag counts in reco (p,|cos(theta)|,dN/deta), tag-specific 3D responses, weighted truth tag numerators, and unweighted matched/gen denominators.").Write();
   TNamed("flat_axes", "mu=(dN/deta,p,|cos(theta)|); activity bins=16 uniform, p bins variable, |cos(theta)| bins variable").Write();
   TH1D hPEdges("hPEdges", "p bin edges;edge index;p (GeV/c)", static_cast<int>(pEdges.size()) - 1, 0.5, static_cast<double>(pEdges.size()) - 0.5);
   for (int i = 0; i + 1 < static_cast<int>(pEdges.size()); ++i)
      hPEdges.SetBinContent(i + 1, pEdges[i + 1]);
   hPEdges.Write();
   TH1D hAbsCosEdges("hAbsCosEdges", "|cos(theta)| bin edges;edge index;|cos(theta)|", static_cast<int>(cosEdges.size()) - 1, 0.5, static_cast<double>(cosEdges.size()) - 0.5);
   for (int i = 0; i + 1 < static_cast<int>(cosEdges.size()); ++i)
      hAbsCosEdges.SetBinContent(i + 1, cosEdges[i + 1]);
   hAbsCosEdges.Write();
   outputFile.Close();

   std::cout << "Wrote " << outputPath << "\n";
   std::cout << "Processed entries: " << nProcessed << ", PassAll: " << nPassAll << "\n";
   return 0;
}
