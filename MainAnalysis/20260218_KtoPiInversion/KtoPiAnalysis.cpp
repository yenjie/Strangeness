//============================================================
// KtoPiAnalysis.cpp  (updated for 3x3, pT-dependent PID matrix)
//
// Computes kaon, pion and proton yields vs Nch_tag and
// K/pi and p/pi ratios using the Strangeness tree messenger.
//
// New in this version
// -------------------
//  * Proton included in the PID unfolding (3×3 matrix).
//  * All PID corrections are done differential in pT.
//  * Flexible pT binning controlled from the command line:
//      - NPtBins / PtMin / PtMax  -> uniform binning
//      - PtBinEdges=a,b,c,...     -> arbitrary bin edges
//  * For each Nch_tag bin we build
//      - raw reconstructed pT spectra (K, pi, p)
//      - PID–corrected pT spectra via 3×3 matrix inversion
//      - pT–integrated yields and K/pi, p/pi vs Nch_tag
//  * Optional PID observation mode:
//      - exclusive: one observed tag per track (legacy K > pi > p tie rule)
//      - inclusive: every species with PID score >= 2 is filled, so duplicate
//        candidate counts are allowed in the raw spectra
//
// The existing generator–level mode (IsGen=true) is kept:
// we still count K/π/p at truth level as a cross–check.
//============================================================

#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <sstream>
#include <cstdlib>
#include <algorithm>
#include <unordered_map>

// ROOT
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TNtuple.h"

// Project common code
#include "utilities.h"      // smartWrite, etc.
#include "helpMessage.h"    // printHelpMessage()
#include "CommandLine.h"    // CommandLine parser
#include "ProgressBar.h"    // nice progress bar

// Strangeness tree messenger
#include "StrangenessMessenger.h"

using namespace std;

//------------------------------------------------------------
// Small helper: parse comma-separated list of doubles
// Example: "0.2,0.4,0.6,1.0,2.0,5.0"
//------------------------------------------------------------
static std::vector<double> ParseDoubleList(const std::string &input)
{
   std::vector<double> values;
   std::stringstream ss(input);
   std::string token;

   while (std::getline(ss, token, ','))
   {
      if (token.empty())
         continue;
      values.push_back(std::atof(token.c_str()));
   }

   return values;
}

//============================================================
// Simple parameter container for this analysis
//============================================================
struct KtoPiParameters
{
   std::string input;
   std::string output;

   int    MaxNchTag;   // max Nch_tag, overflow goes into last bin
   int    MaxEvents;   // max events to process (-1 = all)
   double EcmRef;      // reference energy in GeV (91.2)
   int    MinNch;      // Nch >= MinNch
   double MinTheta;    // in radians
   double MaxTheta;    // in radians

   bool   IsGen;       // if true, count K/pi/p at generator level
   bool   UseMCTruthMatrix;  // if true, build PID matrix from MC truth matching (closure mode)
   bool   UsePassAllSelection; // if true, use archived event-selection bit instead of recomputing cuts
   bool   UseCentralEtaNtag; // if true, define reco Ntag using |eta|<0.5 tracks
   bool   UsePIDFiducial;    // if true, apply track-level PID fiducial in |cos(theta_track)|
   double PIDTrackAbsCosMin; // lower edge for |cos(theta_track)|
   double PIDTrackAbsCosMax; // upper edge for |cos(theta_track)|
   int    PIDTieMode;        // 0=legacy priority K>pi>p, 1=ties become untagged
   bool   UseInclusivePIDObservation; // if true, fill all K/pi/p candidates with PID>=2

   // pT binning
   int    NPtBins;           // number of pT bins (uniform mode)
   double PtMin;             // min pT (GeV/c)
   double PtMax;             // max pT (GeV/c)
   double NtagPtMin;         // min pT for Ntag counting
   std::vector<double> PtBinEdges;  // if non-empty, overrides NPtBins/PtMin/PtMax

   KtoPiParameters()
      : input("sample/Strangeness/merged_pythia_v2.5.root")
      , output("output/KtoPi.root")
      , MaxNchTag(60)
      , MaxEvents(-1)
      , EcmRef(91.2)
      , MinNch(7)
      , MinTheta(30.0 * TMath::Pi() / 180.0)
      , MaxTheta(150.0 * TMath::Pi() / 180.0)
      , IsGen(false)
      , UseMCTruthMatrix(false)
      , UsePassAllSelection(true)
      , UseCentralEtaNtag(false)
      , UsePIDFiducial(true)
      , PIDTrackAbsCosMin(0.15)
      , PIDTrackAbsCosMax(0.675)
      , PIDTieMode(0)
      , UseInclusivePIDObservation(false)
      , NPtBins(12)
      , PtMin(0.4)
      , PtMax(5.0)
      , NtagPtMin(0.2)
   {
   }
};

static bool IsChargedPDG(long long pdg)
{
   const long long apdg = (pdg >= 0 ? pdg : -pdg);
   static std::unordered_map<long long, bool> cache;
   auto it = cache.find(apdg);
   if (it != cache.end())
      return it->second;

   bool charged = false;
   // Common stable charged particles and long-lived charged hadrons.
   if (apdg == 11 || apdg == 13 || apdg == 15 ||          // e, mu, tau
       apdg == 211 || apdg == 321 || apdg == 2212 ||      // pi, K, p
       apdg == 3112 || apdg == 3222 || apdg == 3312 ||    // Sigma-, Sigma+, Xi-
       apdg == 3334 || apdg == 411 || apdg == 431 ||      // Omega-, D+, Ds+
       apdg == 521 || apdg == 541 || apdg == 24)          // B+, Bc+, W+
      charged = true;

   cache[apdg] = charged;
   return charged;
}

static bool ComputeAxisRapidity(double px, double py, double pz, double e,
                                double ax, double ay, double az, double &rapidity)
{
   const double norm = std::sqrt(ax * ax + ay * ay + az * az);
   if (norm <= 0.0 || e <= 0.0)
      return false;

   const double pLong = (px * ax + py * ay + pz * az) / norm;
   const double plus = e + pLong;
   const double minus = e - pLong;
   if (plus <= 0.0 || minus <= 0.0)
      return false;

   rapidity = 0.5 * std::log(plus / minus);
   return std::isfinite(rapidity);
}

//============================================================
// Analyzer class
//============================================================
class KtoPiAnalyzer
{
public:
   TFile *inf;
   TFile *outf;
   StrangenessTreeMessenger *M;
   bool HasRecoMatchingBranches;
   bool HasGenMatchingBranches;
   double RecoEfficiencyKExtra[STRANGE_MAX_RECO];
   double RecoEfficiencyPiExtra[STRANGE_MAX_RECO];
   double RecoEfficiencyPExtra[STRANGE_MAX_RECO];
   double RecoGenEfficiencyKExtra[STRANGE_MAX_RECO];
   double RecoGenEfficiencyPiExtra[STRANGE_MAX_RECO];
   double RecoGenEfficiencyPExtra[STRANGE_MAX_RECO];

   // 1D raw (uncorrected) yields vs Nch_tag
   TH1D *hK;
   TH1D *hPi;
   TH1D *hP;
   TH1D *hKoverPi;
   TH1D *hPoverPi;

   // 1D PID–corrected yields vs Nch_tag (pT–integrated)
   TH1D *hKCorrected;
   TH1D *hPiCorrected;
   TH1D *hPCorrected;
   TH1D *hKoverPiCorrected;
   TH1D *hPoverPiCorrected;

   // 1D raw / corrected yields vs reco dNch/deta(|eta|<0.5)
   TH1D *hKDNdEta;
   TH1D *hPiDNdEta;
   TH1D *hPDNdEta;
   TH1D *hKoverPiDNdEta;
   TH1D *hPoverPiDNdEta;
   TH1D *hKCorrectedDNdEta;
   TH1D *hPiCorrectedDNdEta;
   TH1D *hPCorrectedDNdEta;
   TH1D *hKoverPiCorrectedDNdEta;
   TH1D *hPoverPiCorrectedDNdEta;
   TH1D *hKDNdY;
   TH1D *hPiDNdY;
   TH1D *hPDNdY;
   TH1D *hKoverPiDNdY;
   TH1D *hPoverPiDNdY;
   TH1D *hKCorrectedDNdY;
   TH1D *hPiCorrectedDNdY;
   TH1D *hPCorrectedDNdY;
   TH1D *hKoverPiCorrectedDNdY;
   TH1D *hPoverPiCorrectedDNdY;

   // 2D pT spectra: x = Nch_tag, y = pT
   TH2D *hKPt;
   TH2D *hPiPt;
   TH2D *hPPt;
   TH2D *hUPt;

   // 2D PID–corrected pT spectra
   TH2D *hKPtCorrected;
   TH2D *hPiPtCorrected;
   TH2D *hPPtCorrected;

   // 2D pT spectra: x = reco dNch/deta(|eta|<0.5), y = pT
   TH2D *hKPtDNdEta;
   TH2D *hPiPtDNdEta;
   TH2D *hPPtDNdEta;
   TH2D *hUPtDNdEta;
   TH2D *hKPtCorrectedDNdEta;
   TH2D *hPiPtCorrectedDNdEta;
   TH2D *hPPtCorrectedDNdEta;
   TH2D *hKPtDNdY;
   TH2D *hPiPtDNdY;
   TH2D *hPPtDNdY;
   TH2D *hUPtDNdY;
   TH2D *hKPtCorrectedDNdY;
   TH2D *hPiPtCorrectedDNdY;
   TH2D *hPPtCorrectedDNdY;

   // MC-only response and truth-yield histograms for Ntag unfolding
   TH2D *hNtagResponse;   // x=true Ntag, y=reco Ntag
   TH2D *hNtagResponseK;  // x=true Ntag, y=reco Ntag (K-yield weighted)
   TH2D *hNtagResponsePi; // x=true Ntag, y=reco Ntag (pi-yield weighted)
   TH2D *hNtagResponseP;  // x=true Ntag, y=reco Ntag (p-yield weighted)
   TH1D *hNtagTrue;
   TH1D *hNtagReco;
   TH1D *hKTrueNtag;
   TH1D *hPiTrueNtag;
   TH1D *hPTrueNtag;
   TH2D *hDNdEtaResponse;   // x=true dNch/deta(|eta|<0.5), y=reco dNch/deta(|eta|<0.5)
   TH2D *hDNdEtaResponseK;  // K-yield weighted
   TH2D *hDNdEtaResponsePi; // pi-yield weighted
   TH2D *hDNdEtaResponseP;  // p-yield weighted
   TH1D *hDNdEtaTrue;
   TH1D *hDNdEtaReco;
   TH1D *hKTruedNdEta;
   TH1D *hPiTruedNdEta;
   TH1D *hPTruedNdEta;
   TH2D *hDNdYResponse;
   TH2D *hDNdYResponseK;
   TH2D *hDNdYResponsePi;
   TH2D *hDNdYResponseP;
   TH1D *hDNdYTrue;
   TH1D *hDNdYReco;
   TH1D *hKTruedNdY;
   TH1D *hPiTruedNdY;
   TH1D *hPTruedNdY;

   // Parameters
   KtoPiParameters par;

   // pT binning (copied from par and finalized in ctor)
   std::vector<double> PtBinEdges; // size = NPtBins+1
   int NPtBins;
   int NNchBins;                   // number of Nch_tag bins in histograms

   // Per-(Nch_tag, pT) averages of PID efficiencies (3×3 matrix)
   //
   // Indexing convention:
   //   flatIndex = (iNchBin-1)*NPtBins + (iPtBin-1)
   //
   // Stored entries correspond to:
   //   KAsK  etc = Prob(tag==K | true==K) etc.
   //
   std::vector<double> SumKAsK;
   std::vector<double> SumKAsPi;
   std::vector<double> SumKAsP;

   std::vector<double> SumPiAsK;
   std::vector<double> SumPiAsPi;
   std::vector<double> SumPiAsP;

   std::vector<double> SumPAsK;
   std::vector<double> SumPAsPi;
   std::vector<double> SumPAsP;
   std::vector<double> SumRecoEffK;
   std::vector<double> SumRecoEffPi;
   std::vector<double> SumRecoEffP;
   std::vector<double> SumGenEffK;
   std::vector<double> SumGenEffPi;
   std::vector<double> SumGenEffP;

   std::vector<double> SumKAsKDNdEta;
   std::vector<double> SumKAsPiDNdEta;
   std::vector<double> SumKAsPDNdEta;

   std::vector<double> SumPiAsKDNdEta;
   std::vector<double> SumPiAsPiDNdEta;
   std::vector<double> SumPiAsPDNdEta;

   std::vector<double> SumPAsKDNdEta;
   std::vector<double> SumPAsPiDNdEta;
   std::vector<double> SumPAsPDNdEta;
   std::vector<double> SumRecoEffKDNdEta;
   std::vector<double> SumRecoEffPiDNdEta;
   std::vector<double> SumRecoEffPDNdEta;
   std::vector<double> SumGenEffKDNdEta;
   std::vector<double> SumGenEffPiDNdEta;
   std::vector<double> SumGenEffPDNdEta;
   std::vector<double> SumKAsKDNdY;
   std::vector<double> SumKAsPiDNdY;
   std::vector<double> SumKAsPDNdY;
   std::vector<double> SumPiAsKDNdY;
   std::vector<double> SumPiAsPiDNdY;
   std::vector<double> SumPiAsPDNdY;
   std::vector<double> SumPAsKDNdY;
   std::vector<double> SumPAsPiDNdY;
   std::vector<double> SumPAsPDNdY;
   std::vector<double> SumRecoEffKDNdY;
   std::vector<double> SumRecoEffPiDNdY;
   std::vector<double> SumRecoEffPDNdY;
   std::vector<double> SumGenEffKDNdY;
   std::vector<double> SumGenEffPiDNdY;
   std::vector<double> SumGenEffPDNdY;

   std::vector<long long> CountEffTracks;
   std::vector<long long> CountTrueK;
   std::vector<long long> CountTruePi;
   std::vector<long long> CountTrueP;
   std::vector<long long> CountGenK;
   std::vector<long long> CountGenPi;
   std::vector<long long> CountGenP;
   std::vector<long long> CountEffTracksDNdEta;
   std::vector<long long> CountEffTracksDNdY;
   long long NPIDPassTagTracks;
   long long NPIDTieTracks;

public:
   KtoPiAnalyzer(const KtoPiParameters &apar)
      : inf(nullptr)
      , outf(nullptr)
      , M(nullptr)
      , HasRecoMatchingBranches(false)
      , HasGenMatchingBranches(false)
      , hK(nullptr)
      , hPi(nullptr)
      , hP(nullptr)
      , hKoverPi(nullptr)
      , hPoverPi(nullptr)
      , hKCorrected(nullptr)
      , hPiCorrected(nullptr)
      , hPCorrected(nullptr)
      , hKoverPiCorrected(nullptr)
      , hPoverPiCorrected(nullptr)
      , hKDNdEta(nullptr)
      , hPiDNdEta(nullptr)
      , hPDNdEta(nullptr)
      , hKoverPiDNdEta(nullptr)
      , hPoverPiDNdEta(nullptr)
      , hKCorrectedDNdEta(nullptr)
      , hPiCorrectedDNdEta(nullptr)
      , hPCorrectedDNdEta(nullptr)
      , hKoverPiCorrectedDNdEta(nullptr)
      , hPoverPiCorrectedDNdEta(nullptr)
      , hKDNdY(nullptr)
      , hPiDNdY(nullptr)
      , hPDNdY(nullptr)
      , hKoverPiDNdY(nullptr)
      , hPoverPiDNdY(nullptr)
      , hKCorrectedDNdY(nullptr)
      , hPiCorrectedDNdY(nullptr)
      , hPCorrectedDNdY(nullptr)
      , hKoverPiCorrectedDNdY(nullptr)
      , hPoverPiCorrectedDNdY(nullptr)
      , hKPt(nullptr)
      , hPiPt(nullptr)
      , hPPt(nullptr)
      , hUPt(nullptr)
      , hKPtCorrected(nullptr)
      , hPiPtCorrected(nullptr)
      , hPPtCorrected(nullptr)
      , hKPtDNdEta(nullptr)
      , hPiPtDNdEta(nullptr)
      , hPPtDNdEta(nullptr)
      , hUPtDNdEta(nullptr)
      , hKPtCorrectedDNdEta(nullptr)
      , hPiPtCorrectedDNdEta(nullptr)
      , hPPtCorrectedDNdEta(nullptr)
      , hKPtDNdY(nullptr)
      , hPiPtDNdY(nullptr)
      , hPPtDNdY(nullptr)
      , hUPtDNdY(nullptr)
      , hKPtCorrectedDNdY(nullptr)
      , hPiPtCorrectedDNdY(nullptr)
      , hPPtCorrectedDNdY(nullptr)
      , hNtagResponse(nullptr)
      , hNtagResponseK(nullptr)
      , hNtagResponsePi(nullptr)
      , hNtagResponseP(nullptr)
      , hNtagTrue(nullptr)
      , hNtagReco(nullptr)
      , hKTrueNtag(nullptr)
      , hPiTrueNtag(nullptr)
      , hPTrueNtag(nullptr)
      , hDNdEtaResponse(nullptr)
      , hDNdEtaResponseK(nullptr)
      , hDNdEtaResponsePi(nullptr)
      , hDNdEtaResponseP(nullptr)
      , hDNdEtaTrue(nullptr)
      , hDNdEtaReco(nullptr)
      , hKTruedNdEta(nullptr)
      , hPiTruedNdEta(nullptr)
      , hPTruedNdEta(nullptr)
      , hDNdYResponse(nullptr)
      , hDNdYResponseK(nullptr)
      , hDNdYResponsePi(nullptr)
      , hDNdYResponseP(nullptr)
      , hDNdYTrue(nullptr)
      , hDNdYReco(nullptr)
      , hKTruedNdY(nullptr)
      , hPiTruedNdY(nullptr)
      , hPTruedNdY(nullptr)
      , NPIDPassTagTracks(0)
      , NPIDTieTracks(0)
      , par(apar)
      , PtBinEdges()
      , NPtBins(0)
      , NNchBins(0)
   {
      // Open input
      inf = new TFile(par.input.c_str());
      if (inf == nullptr || inf->IsZombie())
      {
         cerr << "Error: cannot open input file '" << par.input << "'" << endl;
         return;
      }

      // Attach messenger to tree "Tree"
      M = new StrangenessTreeMessenger(*inf, std::string("Tree"));

      // Open output
      outf = new TFile(par.output.c_str(), "RECREATE");
      if (outf == nullptr || outf->IsZombie())
      {
         cerr << "Error: cannot create output file '" << par.output << "'" << endl;
         return;
      }
      outf->cd();

      // v2.2 trees provide species-dependent matching efficiencies.
      if (M != nullptr && M->Tree != nullptr)
      {
         bool hasGenK  = (M->Tree->GetBranch("RecoGenEfficiencyK")  != nullptr);
         bool hasGenPi = (M->Tree->GetBranch("RecoGenEfficiencyPi") != nullptr);
         bool hasGenP  = (M->Tree->GetBranch("RecoGenEfficiencyP")  != nullptr);
         bool hasRecK  = (M->Tree->GetBranch("RecoEfficiencyK")  != nullptr);
         bool hasRecPi = (M->Tree->GetBranch("RecoEfficiencyPi") != nullptr);
         bool hasRecP  = (M->Tree->GetBranch("RecoEfficiencyP")  != nullptr);

         if (hasRecK && hasRecPi && hasRecP)
         {
            HasRecoMatchingBranches = true;
            M->Tree->SetBranchAddress("RecoEfficiencyK",  RecoEfficiencyKExtra);
            M->Tree->SetBranchAddress("RecoEfficiencyPi", RecoEfficiencyPiExtra);
            M->Tree->SetBranchAddress("RecoEfficiencyP",  RecoEfficiencyPExtra);
         }
         if (hasGenK && hasGenPi && hasGenP)
         {
            HasGenMatchingBranches = true;
            M->Tree->SetBranchAddress("RecoGenEfficiencyK",  RecoGenEfficiencyKExtra);
            M->Tree->SetBranchAddress("RecoGenEfficiencyPi", RecoGenEfficiencyPiExtra);
            M->Tree->SetBranchAddress("RecoGenEfficiencyP",  RecoGenEfficiencyPExtra);
         }
      }

      //--------------------------------------------------
      // Finalize pT binning
      //--------------------------------------------------
      if (!par.PtBinEdges.empty())
      {
         PtBinEdges = par.PtBinEdges;
         std::sort(PtBinEdges.begin(), PtBinEdges.end());
         if (PtBinEdges.size() < 2)
         {
            cerr << "Warning: PtBinEdges has fewer than 2 entries. "
                 << "Falling back to uniform pT binning." << endl;
            PtBinEdges.clear();
         }
      }

      if (PtBinEdges.empty())
      {
         // Build uniform binning from PtMin / PtMax / NPtBins
         PtBinEdges.resize(par.NPtBins + 1);
         const double dPt = (par.PtMax - par.PtMin) / par.NPtBins;
         for (int i = 0; i <= par.NPtBins; ++i)
            PtBinEdges[i] = par.PtMin + i * dPt;
      }

      NPtBins = static_cast<int>(PtBinEdges.size()) - 1;

      //--------------------------------------------------
      // Book 1D histograms vs Nch_tag
      //--------------------------------------------------
      const int maxNchTag = par.MaxNchTag;
      const int nbinsNch  = maxNchTag / 4 + 1;   // same choice as original macro

      NNchBins = nbinsNch;

      // Kaon, Proton and Pion Spectra vs NchTag
      hK = new TH1D("hK",
                    "Kaon candidates vs N_{ch}^{tag};N_{ch}^{tag};Yield (sum over tracks)",
                    nbinsNch, -0.5, maxNchTag + 0.5);

      hPi = (TH1D *)hK->Clone("hPi");
      hPi->SetTitle("Pion candidates vs N_{ch}^{tag};N_{ch}^{tag};Yield (sum over tracks)");

      hP = (TH1D *)hK->Clone("hP");
      hP->SetTitle("Proton candidates vs N_{ch}^{tag};N_{ch}^{tag};Yield (sum over tracks)");

      hK->Sumw2();
      hPi->Sumw2();
      hP->Sumw2();

      hKoverPi = nullptr;
      hPoverPi = nullptr;

      // Corrected 1D histograms (will be filled after matrix inversion)
      hKCorrected = (TH1D *)hK->Clone("hKCorrected");
      hKCorrected->SetTitle("PID-corrected K yield vs N_{ch}^{tag};N_{ch}^{tag};Corrected N_{K}");
      hKCorrected->Reset();
      hKCorrected->Sumw2();

      hPiCorrected = (TH1D *)hK->Clone("hPiCorrected");
      hPiCorrected->SetTitle("PID-corrected #pi yield vs N_{ch}^{tag};N_{ch}^{tag};Corrected N_{#pi}");
      hPiCorrected->Reset();
      hPiCorrected->Sumw2();

      hPCorrected = (TH1D *)hK->Clone("hPCorrected");
      hPCorrected->SetTitle("PID-corrected p yield vs N_{ch}^{tag};N_{ch}^{tag};Corrected N_{p}");
      hPCorrected->Reset();
      hPCorrected->Sumw2();

      hKoverPiCorrected = nullptr;
      hPoverPiCorrected = nullptr;

      hKDNdEta = (TH1D *)hK->Clone("hKDNdEta");
      hKDNdEta->SetTitle("Kaon candidates vs reco dN_{ch}/d#eta(|#eta|<0.5);dN_{ch}/d#eta (reco, |#eta|<0.5);Yield (sum over tracks)");
      hKDNdEta->Reset();
      hKDNdEta->Sumw2();

      hPiDNdEta = (TH1D *)hKDNdEta->Clone("hPiDNdEta");
      hPiDNdEta->SetTitle("Pion candidates vs reco dN_{ch}/d#eta(|#eta|<0.5);dN_{ch}/d#eta (reco, |#eta|<0.5);Yield (sum over tracks)");

      hPDNdEta = (TH1D *)hKDNdEta->Clone("hPDNdEta");
      hPDNdEta->SetTitle("Proton candidates vs reco dN_{ch}/d#eta(|#eta|<0.5);dN_{ch}/d#eta (reco, |#eta|<0.5);Yield (sum over tracks)");

      hKoverPiDNdEta = nullptr;
      hPoverPiDNdEta = nullptr;

      hKCorrectedDNdEta = (TH1D *)hKDNdEta->Clone("hKCorrectedDNdEta");
      hKCorrectedDNdEta->SetTitle("PID-corrected K yield vs reco dN_{ch}/d#eta(|#eta|<0.5);dN_{ch}/d#eta (reco, |#eta|<0.5);Corrected N_{K}");
      hKCorrectedDNdEta->Reset();
      hKCorrectedDNdEta->Sumw2();

      hPiCorrectedDNdEta = (TH1D *)hKDNdEta->Clone("hPiCorrectedDNdEta");
      hPiCorrectedDNdEta->SetTitle("PID-corrected #pi yield vs reco dN_{ch}/d#eta(|#eta|<0.5);dN_{ch}/d#eta (reco, |#eta|<0.5);Corrected N_{#pi}");
      hPiCorrectedDNdEta->Reset();
      hPiCorrectedDNdEta->Sumw2();

      hPCorrectedDNdEta = (TH1D *)hKDNdEta->Clone("hPCorrectedDNdEta");
      hPCorrectedDNdEta->SetTitle("PID-corrected p yield vs reco dN_{ch}/d#eta(|#eta|<0.5);dN_{ch}/d#eta (reco, |#eta|<0.5);Corrected N_{p}");
      hPCorrectedDNdEta->Reset();
      hPCorrectedDNdEta->Sumw2();

      hKoverPiCorrectedDNdEta = nullptr;
      hPoverPiCorrectedDNdEta = nullptr;

      hKDNdY = (TH1D *)hK->Clone("hKDNdY");
      hKDNdY->SetTitle("Kaon candidates vs reco dN_{ch}/dy(|y_{T}|<0.5);dN_{ch}/dy (reco, thrust axis, |y_{T}|<0.5);Yield (sum over tracks)");
      hKDNdY->Reset();
      hKDNdY->Sumw2();

      hPiDNdY = (TH1D *)hKDNdY->Clone("hPiDNdY");
      hPiDNdY->SetTitle("Pion candidates vs reco dN_{ch}/dy(|y_{T}|<0.5);dN_{ch}/dy (reco, thrust axis, |y_{T}|<0.5);Yield (sum over tracks)");

      hPDNdY = (TH1D *)hKDNdY->Clone("hPDNdY");
      hPDNdY->SetTitle("Proton candidates vs reco dN_{ch}/dy(|y_{T}|<0.5);dN_{ch}/dy (reco, thrust axis, |y_{T}|<0.5);Yield (sum over tracks)");

      hKoverPiDNdY = nullptr;
      hPoverPiDNdY = nullptr;

      hKCorrectedDNdY = (TH1D *)hKDNdY->Clone("hKCorrectedDNdY");
      hKCorrectedDNdY->SetTitle("PID-corrected K yield vs reco dN_{ch}/dy(|y_{T}|<0.5);dN_{ch}/dy (reco, thrust axis, |y_{T}|<0.5);Corrected N_{K}");
      hKCorrectedDNdY->Reset();
      hKCorrectedDNdY->Sumw2();

      hPiCorrectedDNdY = (TH1D *)hKDNdY->Clone("hPiCorrectedDNdY");
      hPiCorrectedDNdY->SetTitle("PID-corrected #pi yield vs reco dN_{ch}/dy(|y_{T}|<0.5);dN_{ch}/dy (reco, thrust axis, |y_{T}|<0.5);Corrected N_{#pi}");
      hPiCorrectedDNdY->Reset();
      hPiCorrectedDNdY->Sumw2();

      hPCorrectedDNdY = (TH1D *)hKDNdY->Clone("hPCorrectedDNdY");
      hPCorrectedDNdY->SetTitle("PID-corrected p yield vs reco dN_{ch}/dy(|y_{T}|<0.5);dN_{ch}/dy (reco, thrust axis, |y_{T}|<0.5);Corrected N_{p}");
      hPCorrectedDNdY->Reset();
      hPCorrectedDNdY->Sumw2();

      hKoverPiCorrectedDNdY = nullptr;
      hPoverPiCorrectedDNdY = nullptr;

      //--------------------------------------------------
      // Book 2D pT spectra: x = Nch_tag, y = pT
      //--------------------------------------------------
      const double nchMin = hK->GetXaxis()->GetXmin();
      const double nchMax = hK->GetXaxis()->GetXmax();
      const double *ptEdgesArray = &(PtBinEdges[0]);

      hKPt = new TH2D("hKPt",
                      "Kaon candidates;N_{ch}^{tag};p_{T} (GeV/c)",
                      nbinsNch, nchMin, nchMax,
                      NPtBins, ptEdgesArray);

      hPiPt = new TH2D("hPiPt",
                       "Pion candidates;N_{ch}^{tag};p_{T} (GeV/c)",
                       nbinsNch, nchMin, nchMax,
                       NPtBins, ptEdgesArray);

      hPPt = new TH2D("hPPt",
                      "Proton candidates;N_{ch}^{tag};p_{T} (GeV/c)",
                      nbinsNch, nchMin, nchMax,
                      NPtBins, ptEdgesArray);
      hUPt = new TH2D("hUPt",
                      "Untagged charged tracks;N_{ch}^{tag};p_{T} (GeV/c)",
                      nbinsNch, nchMin, nchMax,
                      NPtBins, ptEdgesArray);

      hKPt->Sumw2();
      hPiPt->Sumw2();
      hPPt->Sumw2();
      hUPt->Sumw2();

      hKPtCorrected = (TH2D *)hKPt->Clone("hKPtCorrected");
      hKPtCorrected->SetTitle("PID-corrected K p_{T} spectrum;N_{ch}^{tag};p_{T} (GeV/c)");
      hKPtCorrected->Reset();
      hKPtCorrected->Sumw2();

      hPiPtCorrected = (TH2D *)hPiPt->Clone("hPiPtCorrected");
      hPiPtCorrected->SetTitle("PID-corrected #pi p_{T} spectrum;N_{ch}^{tag};p_{T} (GeV/c)");
      hPiPtCorrected->Reset();
      hPiPtCorrected->Sumw2();

      hPPtCorrected = (TH2D *)hPPt->Clone("hPPtCorrected");
      hPPtCorrected->SetTitle("PID-corrected p p_{T} spectrum;N_{ch}^{tag};p_{T} (GeV/c)");
      hPPtCorrected->Reset();
      hPPtCorrected->Sumw2();

      hKPtDNdEta = (TH2D *)hKPt->Clone("hKPtDNdEta");
      hKPtDNdEta->SetTitle("Kaon candidates;dN_{ch}/d#eta (reco, |#eta|<0.5);p_{T} (GeV/c)");
      hKPtDNdEta->Reset();
      hKPtDNdEta->Sumw2();

      hPiPtDNdEta = (TH2D *)hPiPt->Clone("hPiPtDNdEta");
      hPiPtDNdEta->SetTitle("Pion candidates;dN_{ch}/d#eta (reco, |#eta|<0.5);p_{T} (GeV/c)");
      hPiPtDNdEta->Reset();
      hPiPtDNdEta->Sumw2();

      hPPtDNdEta = (TH2D *)hPPt->Clone("hPPtDNdEta");
      hPPtDNdEta->SetTitle("Proton candidates;dN_{ch}/d#eta (reco, |#eta|<0.5);p_{T} (GeV/c)");
      hPPtDNdEta->Reset();
      hPPtDNdEta->Sumw2();

      hUPtDNdEta = (TH2D *)hUPt->Clone("hUPtDNdEta");
      hUPtDNdEta->SetTitle("Untagged charged tracks;dN_{ch}/d#eta (reco, |#eta|<0.5);p_{T} (GeV/c)");
      hUPtDNdEta->Reset();
      hUPtDNdEta->Sumw2();

      hKPtCorrectedDNdEta = (TH2D *)hKPtDNdEta->Clone("hKPtCorrectedDNdEta");
      hKPtCorrectedDNdEta->SetTitle("PID-corrected K p_{T} spectrum;dN_{ch}/d#eta (reco, |#eta|<0.5);p_{T} (GeV/c)");
      hKPtCorrectedDNdEta->Reset();
      hKPtCorrectedDNdEta->Sumw2();

      hPiPtCorrectedDNdEta = (TH2D *)hPiPtDNdEta->Clone("hPiPtCorrectedDNdEta");
      hPiPtCorrectedDNdEta->SetTitle("PID-corrected #pi p_{T} spectrum;dN_{ch}/d#eta (reco, |#eta|<0.5);p_{T} (GeV/c)");
      hPiPtCorrectedDNdEta->Reset();
      hPiPtCorrectedDNdEta->Sumw2();

      hPPtCorrectedDNdEta = (TH2D *)hPPtDNdEta->Clone("hPPtCorrectedDNdEta");
      hPPtCorrectedDNdEta->SetTitle("PID-corrected p p_{T} spectrum;dN_{ch}/d#eta (reco, |#eta|<0.5);p_{T} (GeV/c)");
      hPPtCorrectedDNdEta->Reset();
      hPPtCorrectedDNdEta->Sumw2();

      hKPtDNdY = (TH2D *)hKPt->Clone("hKPtDNdY");
      hKPtDNdY->SetTitle("Kaon candidates;dN_{ch}/dy (reco, thrust axis, |y_{T}|<0.5);p_{T} (GeV/c)");
      hKPtDNdY->Reset();
      hKPtDNdY->Sumw2();

      hPiPtDNdY = (TH2D *)hPiPt->Clone("hPiPtDNdY");
      hPiPtDNdY->SetTitle("Pion candidates;dN_{ch}/dy (reco, thrust axis, |y_{T}|<0.5);p_{T} (GeV/c)");
      hPiPtDNdY->Reset();
      hPiPtDNdY->Sumw2();

      hPPtDNdY = (TH2D *)hPPt->Clone("hPPtDNdY");
      hPPtDNdY->SetTitle("Proton candidates;dN_{ch}/dy (reco, thrust axis, |y_{T}|<0.5);p_{T} (GeV/c)");
      hPPtDNdY->Reset();
      hPPtDNdY->Sumw2();

      hUPtDNdY = (TH2D *)hUPt->Clone("hUPtDNdY");
      hUPtDNdY->SetTitle("Untagged charged tracks;dN_{ch}/dy (reco, thrust axis, |y_{T}|<0.5);p_{T} (GeV/c)");
      hUPtDNdY->Reset();
      hUPtDNdY->Sumw2();

      hKPtCorrectedDNdY = (TH2D *)hKPtDNdY->Clone("hKPtCorrectedDNdY");
      hKPtCorrectedDNdY->SetTitle("PID-corrected K p_{T} spectrum;dN_{ch}/dy (reco, thrust axis, |y_{T}|<0.5);p_{T} (GeV/c)");
      hKPtCorrectedDNdY->Reset();
      hKPtCorrectedDNdY->Sumw2();

      hPiPtCorrectedDNdY = (TH2D *)hPiPtDNdY->Clone("hPiPtCorrectedDNdY");
      hPiPtCorrectedDNdY->SetTitle("PID-corrected #pi p_{T} spectrum;dN_{ch}/dy (reco, thrust axis, |y_{T}|<0.5);p_{T} (GeV/c)");
      hPiPtCorrectedDNdY->Reset();
      hPiPtCorrectedDNdY->Sumw2();

      hPPtCorrectedDNdY = (TH2D *)hPPtDNdY->Clone("hPPtCorrectedDNdY");
      hPPtCorrectedDNdY->SetTitle("PID-corrected p p_{T} spectrum;dN_{ch}/dy (reco, thrust axis, |y_{T}|<0.5);p_{T} (GeV/c)");
      hPPtCorrectedDNdY->Reset();
      hPPtCorrectedDNdY->Sumw2();

      // MC-only histograms for multiplicity response/unfolding support
      hNtagResponse = new TH2D("hNtagResponse",
                               "N_{tag}^{ch} response;N_{tag,true}^{ch};N_{tag,reco}^{ch}",
                               nbinsNch, nchMin, nchMax,
                               nbinsNch, nchMin, nchMax);
      hNtagResponse->Sumw2();

      hNtagResponseK = (TH2D *)hNtagResponse->Clone("hNtagResponseK");
      hNtagResponseK->SetTitle("K-weighted N_{tag}^{ch} response;N_{tag,true}^{ch};N_{tag,reco}^{ch}");
      hNtagResponseK->Reset();
      hNtagResponseK->Sumw2();

      hNtagResponsePi = (TH2D *)hNtagResponse->Clone("hNtagResponsePi");
      hNtagResponsePi->SetTitle("#pi-weighted N_{tag}^{ch} response;N_{tag,true}^{ch};N_{tag,reco}^{ch}");
      hNtagResponsePi->Reset();
      hNtagResponsePi->Sumw2();

      hNtagResponseP = (TH2D *)hNtagResponse->Clone("hNtagResponseP");
      hNtagResponseP->SetTitle("p-weighted N_{tag}^{ch} response;N_{tag,true}^{ch};N_{tag,reco}^{ch}");
      hNtagResponseP->Reset();
      hNtagResponseP->Sumw2();

      hNtagTrue = (TH1D *)hK->Clone("hNtagTrue");
      hNtagTrue->SetTitle("True N_{tag}^{ch} distribution;N_{tag,true}^{ch};Events");
      hNtagTrue->Reset();
      hNtagTrue->Sumw2();

      hNtagReco = (TH1D *)hK->Clone("hNtagReco");
      hNtagReco->SetTitle("Reco N_{tag}^{ch} distribution;N_{tag,reco}^{ch};Events");
      hNtagReco->Reset();
      hNtagReco->Sumw2();

      hKTrueNtag = (TH1D *)hK->Clone("hKTrueNtag");
      hKTrueNtag->SetTitle("Generator-level K yield vs true N_{tag}^{ch};N_{tag,true}^{ch};N_{K}^{gen}");
      hKTrueNtag->Reset();
      hKTrueNtag->Sumw2();

      hPiTrueNtag = (TH1D *)hK->Clone("hPiTrueNtag");
      hPiTrueNtag->SetTitle("Generator-level #pi yield vs true N_{tag}^{ch};N_{tag,true}^{ch};N_{#pi}^{gen}");
      hPiTrueNtag->Reset();
      hPiTrueNtag->Sumw2();

      hPTrueNtag = (TH1D *)hK->Clone("hPTrueNtag");
      hPTrueNtag->SetTitle("Generator-level p yield vs true N_{tag}^{ch};N_{tag,true}^{ch};N_{p}^{gen}");
      hPTrueNtag->Reset();
      hPTrueNtag->Sumw2();

      hDNdEtaResponse = new TH2D("hDNdEtaResponse",
                                 "dN_{ch}/d#eta response;dN_{ch}/d#eta (true, |#eta|<0.5);dN_{ch}/d#eta (reco, |#eta|<0.5)",
                                 nbinsNch, nchMin, nchMax,
                                 nbinsNch, nchMin, nchMax);
      hDNdEtaResponse->Sumw2();

      hDNdEtaResponseK = (TH2D *)hDNdEtaResponse->Clone("hDNdEtaResponseK");
      hDNdEtaResponseK->SetTitle("K-weighted dN_{ch}/d#eta response;dN_{ch}/d#eta (true, |#eta|<0.5);dN_{ch}/d#eta (reco, |#eta|<0.5)");
      hDNdEtaResponseK->Reset();
      hDNdEtaResponseK->Sumw2();

      hDNdEtaResponsePi = (TH2D *)hDNdEtaResponse->Clone("hDNdEtaResponsePi");
      hDNdEtaResponsePi->SetTitle("#pi-weighted dN_{ch}/d#eta response;dN_{ch}/d#eta (true, |#eta|<0.5);dN_{ch}/d#eta (reco, |#eta|<0.5)");
      hDNdEtaResponsePi->Reset();
      hDNdEtaResponsePi->Sumw2();

      hDNdEtaResponseP = (TH2D *)hDNdEtaResponse->Clone("hDNdEtaResponseP");
      hDNdEtaResponseP->SetTitle("p-weighted dN_{ch}/d#eta response;dN_{ch}/d#eta (true, |#eta|<0.5);dN_{ch}/d#eta (reco, |#eta|<0.5)");
      hDNdEtaResponseP->Reset();
      hDNdEtaResponseP->Sumw2();

      hDNdEtaTrue = (TH1D *)hK->Clone("hDNdEtaTrue");
      hDNdEtaTrue->SetTitle("True dN_{ch}/d#eta distribution (|#eta|<0.5);dN_{ch}/d#eta (true, |#eta|<0.5);Events");
      hDNdEtaTrue->Reset();
      hDNdEtaTrue->Sumw2();

      hDNdEtaReco = (TH1D *)hKDNdEta->Clone("hDNdEtaReco");
      hDNdEtaReco->SetTitle("Reco dN_{ch}/d#eta distribution (|#eta|<0.5);dN_{ch}/d#eta (reco, |#eta|<0.5);Events");
      hDNdEtaReco->Reset();
      hDNdEtaReco->Sumw2();

      hKTruedNdEta = (TH1D *)hK->Clone("hKTruedNdEta");
      hKTruedNdEta->SetTitle("Generator-level K yield vs true dN_{ch}/d#eta;dN_{ch}/d#eta (true, |#eta|<0.5);N_{K}^{gen}");
      hKTruedNdEta->Reset();
      hKTruedNdEta->Sumw2();

      hPiTruedNdEta = (TH1D *)hK->Clone("hPiTruedNdEta");
      hPiTruedNdEta->SetTitle("Generator-level #pi yield vs true dN_{ch}/d#eta;dN_{ch}/d#eta (true, |#eta|<0.5);N_{#pi}^{gen}");
      hPiTruedNdEta->Reset();
      hPiTruedNdEta->Sumw2();

      hPTruedNdEta = (TH1D *)hK->Clone("hPTruedNdEta");
      hPTruedNdEta->SetTitle("Generator-level p yield vs true dN_{ch}/d#eta;dN_{ch}/d#eta (true, |#eta|<0.5);N_{p}^{gen}");
      hPTruedNdEta->Reset();
      hPTruedNdEta->Sumw2();

      hDNdYResponse = new TH2D("hDNdYResponse",
                               "dN_{ch}/dy response wrt thrust axis;dN_{ch}/dy (true, thrust axis, |y_{T}|<0.5);dN_{ch}/dy (reco, thrust axis, |y_{T}|<0.5)",
                               nbinsNch, nchMin, nchMax,
                               nbinsNch, nchMin, nchMax);
      hDNdYResponse->Sumw2();

      hDNdYResponseK = (TH2D *)hDNdYResponse->Clone("hDNdYResponseK");
      hDNdYResponseK->SetTitle("K-weighted dN_{ch}/dy response;dN_{ch}/dy (true, thrust axis, |y_{T}|<0.5);dN_{ch}/dy (reco, thrust axis, |y_{T}|<0.5)");
      hDNdYResponseK->Reset();
      hDNdYResponseK->Sumw2();

      hDNdYResponsePi = (TH2D *)hDNdYResponse->Clone("hDNdYResponsePi");
      hDNdYResponsePi->SetTitle("#pi-weighted dN_{ch}/dy response;dN_{ch}/dy (true, thrust axis, |y_{T}|<0.5);dN_{ch}/dy (reco, thrust axis, |y_{T}|<0.5)");
      hDNdYResponsePi->Reset();
      hDNdYResponsePi->Sumw2();

      hDNdYResponseP = (TH2D *)hDNdYResponse->Clone("hDNdYResponseP");
      hDNdYResponseP->SetTitle("p-weighted dN_{ch}/dy response;dN_{ch}/dy (true, thrust axis, |y_{T}|<0.5);dN_{ch}/dy (reco, thrust axis, |y_{T}|<0.5)");
      hDNdYResponseP->Reset();
      hDNdYResponseP->Sumw2();

      hDNdYTrue = (TH1D *)hK->Clone("hDNdYTrue");
      hDNdYTrue->SetTitle("True dN_{ch}/dy distribution wrt thrust axis (|y_{T}|<0.5);dN_{ch}/dy (true, thrust axis, |y_{T}|<0.5);Events");
      hDNdYTrue->Reset();
      hDNdYTrue->Sumw2();

      hDNdYReco = (TH1D *)hKDNdY->Clone("hDNdYReco");
      hDNdYReco->SetTitle("Reco dN_{ch}/dy distribution wrt thrust axis (|y_{T}|<0.5);dN_{ch}/dy (reco, thrust axis, |y_{T}|<0.5);Events");
      hDNdYReco->Reset();
      hDNdYReco->Sumw2();

      hKTruedNdY = (TH1D *)hK->Clone("hKTruedNdY");
      hKTruedNdY->SetTitle("Generator-level K yield vs true dN_{ch}/dy;dN_{ch}/dy (true, thrust axis, |y_{T}|<0.5);N_{K}^{gen}");
      hKTruedNdY->Reset();
      hKTruedNdY->Sumw2();

      hPiTruedNdY = (TH1D *)hK->Clone("hPiTruedNdY");
      hPiTruedNdY->SetTitle("Generator-level #pi yield vs true dN_{ch}/dy;dN_{ch}/dy (true, thrust axis, |y_{T}|<0.5);N_{#pi}^{gen}");
      hPiTruedNdY->Reset();
      hPiTruedNdY->Sumw2();

      hPTruedNdY = (TH1D *)hK->Clone("hPTruedNdY");
      hPTruedNdY->SetTitle("Generator-level p yield vs true dN_{ch}/dy;dN_{ch}/dy (true, thrust axis, |y_{T}|<0.5);N_{p}^{gen}");
      hPTruedNdY->Reset();
      hPTruedNdY->Sumw2();

      //--------------------------------------------------
      // Allocate per-(Nch_tag, pT) efficiency accumulators
      //--------------------------------------------------
      const int nCells = NNchBins * NPtBins;

      SumKAsK.assign(nCells + 1, 0.0);
      SumKAsPi.assign(nCells + 1, 0.0);
      SumKAsP.assign(nCells + 1, 0.0);

      SumPiAsK.assign(nCells + 1, 0.0);
      SumPiAsPi.assign(nCells + 1, 0.0);
      SumPiAsP.assign(nCells + 1, 0.0);

      SumPAsK.assign(nCells + 1, 0.0);
      SumPAsPi.assign(nCells + 1, 0.0);
      SumPAsP.assign(nCells + 1, 0.0);
      SumRecoEffK.assign(nCells + 1, 0.0);
      SumRecoEffPi.assign(nCells + 1, 0.0);
      SumRecoEffP.assign(nCells + 1, 0.0);
      SumGenEffK.assign(nCells + 1, 0.0);
      SumGenEffPi.assign(nCells + 1, 0.0);
      SumGenEffP.assign(nCells + 1, 0.0);
      SumKAsKDNdEta.assign(nCells + 1, 0.0);
      SumKAsPiDNdEta.assign(nCells + 1, 0.0);
      SumKAsPDNdEta.assign(nCells + 1, 0.0);
      SumPiAsKDNdEta.assign(nCells + 1, 0.0);
      SumPiAsPiDNdEta.assign(nCells + 1, 0.0);
      SumPiAsPDNdEta.assign(nCells + 1, 0.0);
      SumPAsKDNdEta.assign(nCells + 1, 0.0);
      SumPAsPiDNdEta.assign(nCells + 1, 0.0);
      SumPAsPDNdEta.assign(nCells + 1, 0.0);
      SumRecoEffKDNdEta.assign(nCells + 1, 0.0);
      SumRecoEffPiDNdEta.assign(nCells + 1, 0.0);
      SumRecoEffPDNdEta.assign(nCells + 1, 0.0);
      SumGenEffKDNdEta.assign(nCells + 1, 0.0);
      SumGenEffPiDNdEta.assign(nCells + 1, 0.0);
      SumGenEffPDNdEta.assign(nCells + 1, 0.0);
      SumKAsKDNdY.assign(nCells + 1, 0.0);
      SumKAsPiDNdY.assign(nCells + 1, 0.0);
      SumKAsPDNdY.assign(nCells + 1, 0.0);
      SumPiAsKDNdY.assign(nCells + 1, 0.0);
      SumPiAsPiDNdY.assign(nCells + 1, 0.0);
      SumPiAsPDNdY.assign(nCells + 1, 0.0);
      SumPAsKDNdY.assign(nCells + 1, 0.0);
      SumPAsPiDNdY.assign(nCells + 1, 0.0);
      SumPAsPDNdY.assign(nCells + 1, 0.0);
      SumRecoEffKDNdY.assign(nCells + 1, 0.0);
      SumRecoEffPiDNdY.assign(nCells + 1, 0.0);
      SumRecoEffPDNdY.assign(nCells + 1, 0.0);
      SumGenEffKDNdY.assign(nCells + 1, 0.0);
      SumGenEffPiDNdY.assign(nCells + 1, 0.0);
      SumGenEffPDNdY.assign(nCells + 1, 0.0);

      CountEffTracks.assign(nCells + 1, 0);
      CountEffTracksDNdEta.assign(nCells + 1, 0);
      CountEffTracksDNdY.assign(nCells + 1, 0);
      CountTrueK.assign(nCells + 1, 0);
      CountTruePi.assign(nCells + 1, 0);
      CountTrueP.assign(nCells + 1, 0);
      CountGenK.assign(nCells + 1, 0);
      CountGenPi.assign(nCells + 1, 0);
      CountGenP.assign(nCells + 1, 0);
   }

   ~KtoPiAnalyzer()
   {
      delete hK;
      delete hPi;
      delete hP;
      delete hKoverPi;
      delete hPoverPi;
      delete hKDNdEta;
      delete hPiDNdEta;
      delete hPDNdEta;
      delete hKoverPiDNdEta;
      delete hPoverPiDNdEta;

      delete hKCorrected;
      delete hPiCorrected;
      delete hPCorrected;
      delete hKoverPiCorrected;
      delete hPoverPiCorrected;
      delete hKCorrectedDNdEta;
      delete hPiCorrectedDNdEta;
      delete hPCorrectedDNdEta;
      delete hKoverPiCorrectedDNdEta;
      delete hPoverPiCorrectedDNdEta;
      delete hKDNdY;
      delete hPiDNdY;
      delete hPDNdY;
      delete hKoverPiDNdY;
      delete hPoverPiDNdY;
      delete hKCorrectedDNdY;
      delete hPiCorrectedDNdY;
      delete hPCorrectedDNdY;
      delete hKoverPiCorrectedDNdY;
      delete hPoverPiCorrectedDNdY;

      delete hKPt;
      delete hPiPt;
      delete hPPt;
      delete hUPt;
      delete hKPtCorrected;
      delete hPiPtCorrected;
      delete hPPtCorrected;
      delete hKPtDNdEta;
      delete hPiPtDNdEta;
      delete hPPtDNdEta;
      delete hUPtDNdEta;
      delete hKPtCorrectedDNdEta;
      delete hPiPtCorrectedDNdEta;
      delete hPPtCorrectedDNdEta;
      delete hKPtDNdY;
      delete hPiPtDNdY;
      delete hPPtDNdY;
      delete hUPtDNdY;
      delete hKPtCorrectedDNdY;
      delete hPiPtCorrectedDNdY;
      delete hPPtCorrectedDNdY;
      delete hNtagResponse;
      delete hNtagResponseK;
      delete hNtagResponsePi;
      delete hNtagResponseP;
      delete hNtagTrue;
      delete hNtagReco;
      delete hKTrueNtag;
      delete hPiTrueNtag;
      delete hPTrueNtag;
      delete hDNdEtaResponse;
      delete hDNdEtaResponseK;
      delete hDNdEtaResponsePi;
      delete hDNdEtaResponseP;
      delete hDNdEtaTrue;
      delete hDNdEtaReco;
      delete hKTruedNdEta;
      delete hPiTruedNdEta;
      delete hPTruedNdEta;
      delete hDNdYResponse;
      delete hDNdYResponseK;
      delete hDNdYResponsePi;
      delete hDNdYResponseP;
      delete hDNdYTrue;
      delete hDNdYReco;
      delete hKTruedNdY;
      delete hPiTruedNdY;
      delete hPTruedNdY;

      if (inf)
      {
         inf->Close();
         delete inf;
      }
      if (outf)
      {
         outf->Close();
         delete outf;
      }
      delete M;
   }

   // Map (Nch_bin, pT_bin) -> flat index into vectors
   inline int flatIndex(int iNchBin, int iPtBin) const
   {
      // Safety clamp (bins are 1..NNchBins and 1..NPtBins for histograms)
      if (iNchBin < 1)        iNchBin = 1;
      if (iNchBin > NNchBins) iNchBin = NNchBins;
      if (iPtBin < 1)         iPtBin = 1;
      if (iPtBin > NPtBins)   iPtBin = NPtBins;

      return (iNchBin - 1) * NPtBins + (iPtBin - 1);
   }

   // Invert a 3×3 matrix using cofactors.
   // M[row][col] with row,col = 0..2.
   // Returns false if determinant is (numerically) zero.
   bool invert3x3(const double Mmat[3][3], double Minv[3][3], double &det) const
   {
      const double a = Mmat[0][0];
      const double b = Mmat[0][1];
      const double c = Mmat[0][2];
      const double d = Mmat[1][0];
      const double e = Mmat[1][1];
      const double f = Mmat[1][2];
      const double g = Mmat[2][0];
      const double h = Mmat[2][1];
      const double i = Mmat[2][2];

      det = a*(e*i - f*h) - b*(d*i - f*g) + c*(d*h - e*g);

      if (std::fabs(det) < 1.0e-10)
         return false;

      // Cofactor matrix
      const double C00 =  (e*i - f*h);
      const double C01 = -(d*i - f*g);
      const double C02 =  (d*h - e*g);

      const double C10 = -(b*i - c*h);
      const double C11 =  (a*i - c*g);
      const double C12 = -(a*h - b*g);

      const double C20 =  (b*f - c*e);
      const double C21 = -(a*f - c*d);
      const double C22 =  (a*e - b*d);

      // Inverse = (1/det) * Cofactor^T
      const double invDet = 1.0 / det;

      Minv[0][0] = C00 * invDet;
      Minv[0][1] = C10 * invDet;
      Minv[0][2] = C20 * invDet;

      Minv[1][0] = C01 * invDet;
      Minv[1][1] = C11 * invDet;
      Minv[1][2] = C21 * invDet;

      Minv[2][0] = C02 * invDet;
      Minv[2][1] = C12 * invDet;
      Minv[2][2] = C22 * invDet;

      return true;
   }

   void analyze()
   {
      if (M == nullptr || inf == nullptr || outf == nullptr)
         return;

      //---------------------------------------------------
      // Event loop
      //---------------------------------------------------
      long long nEntries = M->GetEntries();
      if (par.MaxEvents > 0 && par.MaxEvents < nEntries)
         nEntries = par.MaxEvents;

      cout << "Total entries to process: " << nEntries << endl;
      cout << "Using 3-step correction (reco-match -> 3x3 tagging -> gen-match)." << endl;
      cout << "  Reco matching branches: " << (HasRecoMatchingBranches ? "RecoEfficiency*" : "fallback=1") << endl;
      cout << "  Gen matching branches : " << (HasGenMatchingBranches ? "RecoGenEfficiency*" : "fallback=1") << endl;

      auto passPIDFiducialFromMom = [&](double px, double py, double pz) -> bool
      {
         if (!par.UsePIDFiducial)
            return true;
         const double p2 = px * px + py * py + pz * pz;
         if (p2 <= 0.0)
            return false;
         const double absCos = std::abs(pz / std::sqrt(p2));
         return (absCos > par.PIDTrackAbsCosMin && absCos < par.PIDTrackAbsCosMax);
      };

      ProgressBar Bar(cout, nEntries);
      Bar.SetStyle(1);
      long long deltaI = nEntries / 100 + 1;

      const double EcmRef   = par.EcmRef;
      const int    MinNch   = par.MinNch;
      const double MinTheta = par.MinTheta;
      const double MaxTheta = par.MaxTheta;

      for (long long ievt = 0; ievt < nEntries; ++ievt)
      {
         M->GetEntry(ievt);

         if (ievt % deltaI == 0)
         {
            Bar.Update(ievt);
            Bar.Print();
         }

         // Safety: cap NReco to messenger array size
         if (M->NReco > STRANGE_MAX_RECO)
         {
            cerr << "Warning: NReco = " << M->NReco
                 << " > STRANGE_MAX_RECO = " << STRANGE_MAX_RECO
                 << " at entry " << ievt << ".  Clipping to STRANGE_MAX_RECO."
                 << endl;
         }
         int nreco = (M->NReco < STRANGE_MAX_RECO)
                        ? static_cast<int>(M->NReco)
                        : STRANGE_MAX_RECO;

         // Prepare NGen when available (MC), needed for:
         // - IsGen mode
         // - closure matrix mode
         // - Ntag response / truth-prior histograms in reco mode
         int ngen = 0;
         if (M->NGen > 0)
         {
            if (M->NGen > STRANGE_MAX_GEN)
            {
               cerr << "Warning: NGen = " << M->NGen
                    << " > STRANGE_MAX_GEN = " << STRANGE_MAX_GEN
                    << " at entry " << ievt << ".  Clipping to STRANGE_MAX_GEN."
                    << endl;
            }
            ngen = (M->NGen < STRANGE_MAX_GEN)
                      ? static_cast<int>(M->NGen)
                      : STRANGE_MAX_GEN;
         }

         //-------------------------
         // Event selection
         //-------------------------

         if (par.UsePassAllSelection)
         {
            // Nominal v6+ event preselection: use the archived event-selection bit
            // written into the open-data trees instead of recomputing the component cuts.
            if (M->PassAll != 1)
               continue;
         }
         else
         {
            // Legacy fallback: recompute the historical selection from the stored
            // reco observables for debugging and compatibility studies.
            double sumRecoE = 0.0;
            for (int i = 0; i < nreco; ++i)
               sumRecoE += M->RecoE[i];

            if (sumRecoE / EcmRef <= 0.5)
               continue;
            if (M->Nch < MinNch)
               continue;

            double theta = std::acos(M->ThrustZ);
            if (theta <= MinTheta)
               continue;
            if (theta >= MaxTheta)
               continue;
         }

         const double thrustNorm = std::sqrt(M->ThrustX * M->ThrustX + M->ThrustY * M->ThrustY + M->ThrustZ * M->ThrustZ);
         const bool hasThrustAxis = (thrustNorm > 0.0);
         const double thrustX = hasThrustAxis ? (M->ThrustX / thrustNorm) : 0.0;
         const double thrustY = hasThrustAxis ? (M->ThrustY / thrustNorm) : 0.0;
         const double thrustZ = hasThrustAxis ? (M->ThrustZ / thrustNorm) : 1.0;

         //-------------------------
         // Compute Nch_tag for this event
         //-------------------------
         int NchTag = 0;
         int NchEta05Reco = 0;
         int NchY05Reco = 0;
         for (int i = 0; i < nreco; ++i)
         {
            if (M->RecoGoodTrack[i] != 1)
               continue;
            if (M->RecoCharge[i] == 0.0)
               continue;

            const double px = M->RecoPx[i];
            const double py = M->RecoPy[i];
            const double pz = M->RecoPz[i];
            const double pt = std::sqrt(px * px + py * py);

            // Same-variable reco activity estimator for the dNch/deta unfolding.
            // This count uses all charged good tracks in |eta|<0.5 with no PID or pT threshold.
            if (pt > 0.0)
            {
               const double eta = std::asinh(pz / pt);
               if (std::abs(eta) < 0.5)
                  ++NchEta05Reco;
            }
            if (hasThrustAxis)
            {
               double yThrust = 0.0;
               if (ComputeAxisRapidity(px, py, pz, M->RecoE[i], thrustX, thrustY, thrustZ, yThrust) &&
                   std::abs(yThrust) < 0.5)
               {
                  ++NchY05Reco;
               }
            }

            if (pt < par.NtagPtMin)
               continue;

            if (par.UseCentralEtaNtag)
            {
               if (pt <= 0.0)
                  continue;
               const double eta = std::asinh(pz / pt);
               if (std::abs(eta) >= 0.5)
                  continue;
            }
            ++NchTag;
         }

         // Put overflow NchTag into the last visible bin
         if (NchTag > par.MaxNchTag)
            NchTag = par.MaxNchTag;
         if (NchEta05Reco > par.MaxNchTag)
            NchEta05Reco = par.MaxNchTag;
         if (NchY05Reco > par.MaxNchTag)
            NchY05Reco = par.MaxNchTag;

         // Bin index along Nch_tag axis (1..NNchBins)
         int nchBin = hK->GetXaxis()->FindBin(static_cast<double>(NchTag));
         if (nchBin < 1)
            nchBin = 1;
         if (nchBin > NNchBins)
            nchBin = NNchBins;

         int dndetaBin = hKDNdEta->GetXaxis()->FindBin(static_cast<double>(NchEta05Reco));
         if (dndetaBin < 1)
            dndetaBin = 1;
         if (dndetaBin > NNchBins)
            dndetaBin = NNchBins;

         int dndyBin = hKDNdY->GetXaxis()->FindBin(static_cast<double>(NchY05Reco));
         if (dndyBin < 1)
            dndyBin = 1;
         if (dndyBin > NNchBins)
            dndyBin = NNchBins;

         // Build true multiplicity and truth yields (MC only) for response/unfolding support.
         // The truth-side identified yields must use the same fiducial definition as the
         // standalone generator reference so that closure compares identical quantities.
         int NchTagTrue = 0;
         int nKgenEvt = 0;
         int nPigenEvt = 0;
         int nPgenEvt = 0;
         int nChEta05True = 0;
         int nChY05True = 0;
         if (ngen > 0)
         {
            for (int i = 0; i < ngen; ++i)
            {
               const long long pdg = M->GenID[i];
               const long long absPdg = (pdg >= 0 ? pdg : -pdg);
               const long long status = M->GenStatus[i];
               if (status != 1)
                  continue;
               if (!IsChargedPDG(pdg))
                  continue;

               const double genPx = M->GenPx[i];
               const double genPy = M->GenPy[i];
               const double genPz = M->GenPz[i];
               const double genPtAll = std::sqrt(genPx * genPx + genPy * genPy);
               if (genPtAll > 0.0)
               {
                  const double eta = std::asinh(genPz / genPtAll);
                  if (std::abs(eta) < 0.5)
                     ++nChEta05True;

                  if (par.UseCentralEtaNtag && std::abs(eta) >= 0.5)
                     continue;
               }
               if (hasThrustAxis)
               {
                  double yThrust = 0.0;
                  if (ComputeAxisRapidity(genPx, genPy, genPz, M->GenE[i], thrustX, thrustY, thrustZ, yThrust) &&
                      std::abs(yThrust) < 0.5)
                  {
                     ++nChY05True;
                  }
               }

               if (genPtAll < par.NtagPtMin)
                  continue;
               ++NchTagTrue;

               if (absPdg != 211 && absPdg != 321 && absPdg != 2212)
                  continue;
               const double genPt = genPtAll;
               if (genPt < PtBinEdges.front() || genPt >= PtBinEdges.back())
                  continue;
               if (!passPIDFiducialFromMom(genPx, genPy, genPz))
                  continue;
               if (absPdg == 321)  ++nKgenEvt;
               if (absPdg == 211)  ++nPigenEvt;
               if (absPdg == 2212) ++nPgenEvt;
            }

            if (NchTagTrue > par.MaxNchTag)
               NchTagTrue = par.MaxNchTag;
         }
         if (nChEta05True > par.MaxNchTag)
            nChEta05True = par.MaxNchTag;
         if (nChY05True > par.MaxNchTag)
            nChY05True = par.MaxNchTag;

         // Note: reco correction always uses efficiency branches from the tree.
         // MC generator truth is produced in a separate IsGen=true run.

         //-------------------------
         // Fill generator-level yields (IsGen mode)
         //-------------------------
         if (par.IsGen)
         {
            int nKgen  = 0;
            int nPigen = 0;
            int nPgen  = 0;
            const int trueNchAxis = NchTagTrue;

            for (int i = 0; i < ngen; ++i)
            {
               const long long pdg    = M->GenID[i];
               const long long absPdg = (pdg >= 0 ? pdg : -pdg);
               const long long status = M->GenStatus[i];

               // Use stable final-state truth for a meaningful closure target
               if (status != 1)
                  continue;

               const double genPt = std::sqrt(M->GenPx[i] * M->GenPx[i] + M->GenPy[i] * M->GenPy[i]);
               if (genPt < PtBinEdges.front() || genPt >= PtBinEdges.back())
                  continue;
               if (!passPIDFiducialFromMom(M->GenPx[i], M->GenPy[i], M->GenPz[i]))
                  continue;

               // Charged kaons: K+, K-
               if (absPdg == 321)
               {
                  ++nKgen;
                  hKPt->Fill(trueNchAxis, genPt);
               }

               // Charged pions: pi+, pi-
               if (absPdg == 211)
               {
                  ++nPigen;
                  hPiPt->Fill(trueNchAxis, genPt);
               }

               // Protons / anti-protons
               if (absPdg == 2212)
               {
                  ++nPgen;
                  hPPt->Fill(trueNchAxis, genPt);
               }
            }

            hK->Fill(trueNchAxis, nKgen);
            hPi->Fill(trueNchAxis, nPigen);
            hP->Fill(trueNchAxis, nPgen);

            continue; // nothing more to do for this event
         }

         //-------------------------
         // Reco-based PID counting and pT spectra (IsGen == false)
         //-------------------------
         int nK  = 0;
         int nPi = 0;
         int nP  = 0;


         for (int i = 0; i < nreco; ++i)
         {
            if (M->RecoGoodTrack[i] != 1)
               continue;

            const int kTag = static_cast<int>(M->RecoPIDKaon[i]);
            const int piTag = static_cast<int>(M->RecoPIDPion[i]);
            const int pTag = static_cast<int>(M->RecoPIDProton[i]);
            const bool passKaonTag = (kTag >= 2);
            const bool passPionTag = (piTag >= 2);
            const bool passProtonTag = (pTag >= 2);
            const bool passTag = (passKaonTag || passPionTag || passProtonTag);
            if (!passPIDFiducialFromMom(M->RecoPx[i], M->RecoPy[i], M->RecoPz[i]))
               continue;

            bool isKaonTag = false;
            bool isPionTag = false;
            bool isProtonTag = false;
            bool isUntagged = false;

            if (passTag)
            {
               ++NPIDPassTagTracks;
               const int best = std::max(kTag, std::max(piTag, pTag));
               const int nBest = (kTag == best) + (piTag == best) + (pTag == best);
               if (nBest > 1)
                  ++NPIDTieTracks;
            }

            if (par.UseInclusivePIDObservation)
            {
               // v5-style observed spectra: every species with PID score >= 2
               // contributes to its raw spectrum, so duplicate tag counts are
               // allowed across K/pi/p for the same track.
               isKaonTag = passKaonTag;
               isPionTag = passPionTag;
               isProtonTag = passProtonTag;
               isUntagged = !passTag;
            }
            else
            {
               // Exclusive observed PID category: K, pi, p, untagged
               int obsCat = 3; // untagged
               if (passTag)
               {
                  const int best = std::max(kTag, std::max(piTag, pTag));
                  const int nBest = (kTag == best) + (piTag == best) + (pTag == best);
                  if (best < 2)
                  {
                     obsCat = 3;
                  }
                  else if (nBest > 1 && par.PIDTieMode == 1)
                  {
                     obsCat = 3;
                  }
                  else
                  {
                     // Legacy deterministic tie handling (priority K > pi > p).
                     obsCat = 0;
                     if (piTag > kTag && piTag >= pTag)
                        obsCat = 1;
                     if (pTag > kTag && pTag > piTag)
                        obsCat = 2;
                  }
               }
               isKaonTag = (obsCat == 0);
               isPionTag = (obsCat == 1);
               isProtonTag = (obsCat == 2);
               isUntagged = (obsCat == 3);
            }

            // NOTE: you may need to adapt the pT definition below to your
            // StrangenessMessenger. If you do not have a direct RecoPT[],
            // replace the line with something like
            //
            //   double pt = std::sqrt(M->RecoPx[i]*M->RecoPx[i] +
            //                         M->RecoPy[i]*M->RecoPy[i]);
            //
            // according to your tree content.
            double pt = sqrt(M->RecoPx[i]*M->RecoPx[i]+M->RecoPy[i]*M->RecoPy[i]);

            // Restrict to the configured pT range
            if (pt < PtBinEdges.front() || pt >= PtBinEdges.back())
               continue;

            if (isKaonTag)   ++nK;
            if (isPionTag)   ++nPi;
            if (isProtonTag) ++nP;

            int ptBin = hKPt->GetYaxis()->FindBin(pt);
            if (ptBin < 1 || ptBin > NPtBins)
               continue;

            // Raw reconstructed pT spectra
            if (isKaonTag)
               hKPt->Fill(NchTag, pt);
            if (isPionTag)
               hPiPt->Fill(NchTag, pt);
            if (isProtonTag)
               hPPt->Fill(NchTag, pt);
            if (isUntagged)
               hUPt->Fill(NchTag, pt);

            if (isKaonTag)
               hKPtDNdEta->Fill(NchEta05Reco, pt);
            if (isPionTag)
               hPiPtDNdEta->Fill(NchEta05Reco, pt);
            if (isProtonTag)
               hPPtDNdEta->Fill(NchEta05Reco, pt);
            if (isUntagged)
               hUPtDNdEta->Fill(NchEta05Reco, pt);
            if (isKaonTag)
               hKPtDNdY->Fill(NchY05Reco, pt);
            if (isPionTag)
               hPiPtDNdY->Fill(NchY05Reco, pt);
            if (isProtonTag)
               hPPtDNdY->Fill(NchY05Reco, pt);
            if (isUntagged)
               hUPtDNdY->Fill(NchY05Reco, pt);

            // Accumulate PID efficiencies / fake rates
            if (M->RecoCharge[i] == 0.0)
               continue;  // only charged tracks are taggable

            const int idx = flatIndex(nchBin, ptBin);
            const int idxDNdEta = flatIndex(dndetaBin, ptBin);
            const int idxDNdY = flatIndex(dndyBin, ptBin);

            // *** IMPORTANT ***
            // The names below assume that the messenger provides the
            // full 3×3 set of calibration arrays:
            //
            //   RecoEfficiencyKAsK, RecoEfficiencyKAsPi, RecoEfficiencyKAsP
            //   RecoEfficiencyPiAsK, RecoEfficiencyPiAsPi, RecoEfficiencyPiAsP
            //   RecoEfficiencyPAsK, RecoEfficiencyPAsPi, RecoEfficiencyPAsP
            //
            // If your tree uses slightly different names, just edit the
            // corresponding lines here.

            // Unified mode: always use per-track efficiency parametrization
            SumKAsK[idx]   += M->RecoEfficiencyKAsK[i];
            SumKAsPi[idx]  += M->RecoEfficiencyKAsPi[i];
            SumKAsP[idx]   += M->RecoEfficiencyKAsP[i];
            SumKAsKDNdEta[idxDNdEta]  += M->RecoEfficiencyKAsK[i];
            SumKAsPiDNdEta[idxDNdEta] += M->RecoEfficiencyKAsPi[i];
            SumKAsPDNdEta[idxDNdEta]  += M->RecoEfficiencyKAsP[i];
            SumKAsKDNdY[idxDNdY]  += M->RecoEfficiencyKAsK[i];
            SumKAsPiDNdY[idxDNdY] += M->RecoEfficiencyKAsPi[i];
            SumKAsPDNdY[idxDNdY]  += M->RecoEfficiencyKAsP[i];

            SumPiAsK[idx]  += M->RecoEfficiencyPiAsK[i];
            SumPiAsPi[idx] += M->RecoEfficiencyPiAsPi[i];
            SumPiAsP[idx]  += M->RecoEfficiencyPiAsP[i];
            SumPiAsKDNdEta[idxDNdEta]  += M->RecoEfficiencyPiAsK[i];
            SumPiAsPiDNdEta[idxDNdEta] += M->RecoEfficiencyPiAsPi[i];
            SumPiAsPDNdEta[idxDNdEta]  += M->RecoEfficiencyPiAsP[i];
            SumPiAsKDNdY[idxDNdY]  += M->RecoEfficiencyPiAsK[i];
            SumPiAsPiDNdY[idxDNdY] += M->RecoEfficiencyPiAsPi[i];
            SumPiAsPDNdY[idxDNdY]  += M->RecoEfficiencyPiAsP[i];

            SumPAsK[idx]   += M->RecoEfficiencyPAsK[i];
            SumPAsPi[idx]  += M->RecoEfficiencyPAsPi[i];
            SumPAsP[idx]   += M->RecoEfficiencyPAsP[i];
            SumPAsKDNdEta[idxDNdEta]  += M->RecoEfficiencyPAsK[i];
            SumPAsPiDNdEta[idxDNdEta] += M->RecoEfficiencyPAsPi[i];
            SumPAsPDNdEta[idxDNdEta]  += M->RecoEfficiencyPAsP[i];
            SumPAsKDNdY[idxDNdY]  += M->RecoEfficiencyPAsK[i];
            SumPAsPiDNdY[idxDNdY] += M->RecoEfficiencyPAsPi[i];
            SumPAsPDNdY[idxDNdY]  += M->RecoEfficiencyPAsP[i];

            if (HasRecoMatchingBranches)
            {
               SumRecoEffK[idx]  += RecoEfficiencyKExtra[i];
               SumRecoEffPi[idx] += RecoEfficiencyPiExtra[i];
               SumRecoEffP[idx]  += RecoEfficiencyPExtra[i];
               SumRecoEffKDNdEta[idxDNdEta]  += RecoEfficiencyKExtra[i];
               SumRecoEffPiDNdEta[idxDNdEta] += RecoEfficiencyPiExtra[i];
               SumRecoEffPDNdEta[idxDNdEta]  += RecoEfficiencyPExtra[i];
               SumRecoEffKDNdY[idxDNdY]  += RecoEfficiencyKExtra[i];
               SumRecoEffPiDNdY[idxDNdY] += RecoEfficiencyPiExtra[i];
               SumRecoEffPDNdY[idxDNdY]  += RecoEfficiencyPExtra[i];
            }
            else
            {
               SumRecoEffK[idx]  += 1.0;
               SumRecoEffPi[idx] += 1.0;
               SumRecoEffP[idx]  += 1.0;
               SumRecoEffKDNdEta[idxDNdEta]  += 1.0;
               SumRecoEffPiDNdEta[idxDNdEta] += 1.0;
               SumRecoEffPDNdEta[idxDNdEta]  += 1.0;
               SumRecoEffKDNdY[idxDNdY]  += 1.0;
               SumRecoEffPiDNdY[idxDNdY] += 1.0;
               SumRecoEffPDNdY[idxDNdY]  += 1.0;
            }
            if (HasGenMatchingBranches)
            {
               SumGenEffK[idx]  += RecoGenEfficiencyKExtra[i];
               SumGenEffPi[idx] += RecoGenEfficiencyPiExtra[i];
               SumGenEffP[idx]  += RecoGenEfficiencyPExtra[i];
               SumGenEffKDNdEta[idxDNdEta]  += RecoGenEfficiencyKExtra[i];
               SumGenEffPiDNdEta[idxDNdEta] += RecoGenEfficiencyPiExtra[i];
               SumGenEffPDNdEta[idxDNdEta]  += RecoGenEfficiencyPExtra[i];
               SumGenEffKDNdY[idxDNdY]  += RecoGenEfficiencyKExtra[i];
               SumGenEffPiDNdY[idxDNdY] += RecoGenEfficiencyPiExtra[i];
               SumGenEffPDNdY[idxDNdY]  += RecoGenEfficiencyPExtra[i];
            }
            else
            {
               SumGenEffK[idx]  += 1.0;
               SumGenEffPi[idx] += 1.0;
               SumGenEffP[idx]  += 1.0;
               SumGenEffKDNdEta[idxDNdEta]  += 1.0;
               SumGenEffPiDNdEta[idxDNdEta] += 1.0;
               SumGenEffPDNdEta[idxDNdEta]  += 1.0;
               SumGenEffKDNdY[idxDNdY]  += 1.0;
               SumGenEffPiDNdY[idxDNdY] += 1.0;
               SumGenEffPDNdY[idxDNdY]  += 1.0;
            }

            CountEffTracks[idx]++;
            CountEffTracksDNdEta[idxDNdEta]++;
            CountEffTracksDNdY[idxDNdY]++;
         }

         // Event-wise raw yields integrated over pT (sanity check)
         hK->Fill(NchTag, nK);
         hPi->Fill(NchTag, nPi);
         hP->Fill(NchTag, nP);
         hKDNdEta->Fill(NchEta05Reco, nK);
         hPiDNdEta->Fill(NchEta05Reco, nPi);
         hPDNdEta->Fill(NchEta05Reco, nP);
         hKDNdY->Fill(NchY05Reco, nK);
         hPiDNdY->Fill(NchY05Reco, nPi);
         hPDNdY->Fill(NchY05Reco, nP);
         if (hNtagReco != nullptr)
            hNtagReco->Fill(NchTag);
         if (hDNdEtaReco != nullptr)
            hDNdEtaReco->Fill(NchEta05Reco);
         if (hDNdYReco != nullptr)
            hDNdYReco->Fill(NchY05Reco);

         // MC-only response bookkeeping in reco mode
         if (ngen > 0 && hNtagResponse != nullptr)
         {
            const double dNdEtaTrue = static_cast<double>(nChEta05True);
            const double dNdYTrue = static_cast<double>(nChY05True);
            hNtagResponse->Fill(NchTagTrue, NchTag);
            hNtagResponseK->Fill(NchTagTrue, NchTag, nKgenEvt);
            hNtagResponsePi->Fill(NchTagTrue, NchTag, nPigenEvt);
            hNtagResponseP->Fill(NchTagTrue, NchTag, nPgenEvt);
            hNtagTrue->Fill(NchTagTrue);
            hKTrueNtag->Fill(NchTagTrue, nKgenEvt);
            hPiTrueNtag->Fill(NchTagTrue, nPigenEvt);
            hPTrueNtag->Fill(NchTagTrue, nPgenEvt);
            hDNdEtaResponse->Fill(dNdEtaTrue, NchEta05Reco);
            hDNdEtaResponseK->Fill(dNdEtaTrue, NchEta05Reco, nKgenEvt);
            hDNdEtaResponsePi->Fill(dNdEtaTrue, NchEta05Reco, nPigenEvt);
            hDNdEtaResponseP->Fill(dNdEtaTrue, NchEta05Reco, nPgenEvt);
            hDNdEtaTrue->Fill(dNdEtaTrue);
            hKTruedNdEta->Fill(dNdEtaTrue, nKgenEvt);
            hPiTruedNdEta->Fill(dNdEtaTrue, nPigenEvt);
            hPTruedNdEta->Fill(dNdEtaTrue, nPgenEvt);
            hDNdYResponse->Fill(dNdYTrue, NchY05Reco);
            hDNdYResponseK->Fill(dNdYTrue, NchY05Reco, nKgenEvt);
            hDNdYResponsePi->Fill(dNdYTrue, NchY05Reco, nPigenEvt);
            hDNdYResponseP->Fill(dNdYTrue, NchY05Reco, nPgenEvt);
            hDNdYTrue->Fill(dNdYTrue);
            hKTruedNdY->Fill(dNdYTrue, nKgenEvt);
            hPiTruedNdY->Fill(dNdYTrue, nPigenEvt);
            hPTruedNdY->Fill(dNdYTrue, nPgenEvt);
         }
      }

      cout << endl << "Event loop finished." << endl;

      //-------------------------
      // Update titles depending on mode
      //-------------------------
      if (par.IsGen)
      {
         hK->SetTitle ("Generator-level kaons vs N_{ch}^{tag};N_{ch}^{tag};N_{K}^{gen}");
         hPi->SetTitle("Generator-level pions vs N_{ch}^{tag};N_{ch}^{tag};N_{#pi}^{gen}");
         hP->SetTitle ("Generator-level protons vs N_{ch}^{tag};N_{ch}^{tag};N_{p}^{gen}");
      }
      else
      {
         hK->SetTitle ("Kaon candidates vs N_{ch}^{tag};N_{ch}^{tag};Yield (sum over tracks)");
         hPi->SetTitle("Pion candidates vs N_{ch}^{tag};N_{ch}^{tag};Yield (sum over tracks)");
         hP->SetTitle ("Proton candidates vs N_{ch}^{tag};N_{ch}^{tag};Yield (sum over tracks)");
      }

      //-------------------------
      // Build ratios & apply PID correction
      //-------------------------

      // Generator-level case: build raw K/pi and p/pi and stop.
      if (par.IsGen)
      {
         hKoverPi = (TH1D *)hK->Clone("hKoverPi");
         hPoverPi = (TH1D *)hP->Clone("hPoverPi");

         hKoverPi->SetTitle("Generator-level K/#pi yield ratio vs N_{ch}^{tag};N_{ch}^{tag};K/#pi (gen)");
         hPoverPi->SetTitle("Generator-level p/#pi yield ratio vs N_{ch}^{tag};N_{ch}^{tag};p/#pi (gen)");

         hKoverPi->Divide(hPi);
         hPoverPi->Divide(hPi);

         return;  // no PID unfolding in generator mode
      }

      //-------------------------------------------------
      // 3-step PID correction, pT-dependent (reco mode)
      //-------------------------------------------------
      const int nNchBinsLocal = NNchBins;

      // Reset corrected histograms just in case
      hKPtCorrected->Reset();
      hPiPtCorrected->Reset();
      hPPtCorrected->Reset();
      hKCorrected->Reset();
      hPiCorrected->Reset();
      hPCorrected->Reset();
      hKPtCorrectedDNdEta->Reset();
      hPiPtCorrectedDNdEta->Reset();
      hPPtCorrectedDNdEta->Reset();
      hKCorrectedDNdEta->Reset();
      hPiCorrectedDNdEta->Reset();
      hPCorrectedDNdEta->Reset();
      hKPtCorrectedDNdY->Reset();
      hPiPtCorrectedDNdY->Reset();
      hPPtCorrectedDNdY->Reset();
      hKCorrectedDNdY->Reset();
      hPiCorrectedDNdY->Reset();
      hPCorrectedDNdY->Reset();

      auto correctAxis = [&](int axisMode)
      {
         TH2D *hRawK2D = (axisMode == 1) ? hKPtDNdEta : ((axisMode == 2) ? hKPtDNdY : hKPt);
         TH2D *hRawPi2D = (axisMode == 1) ? hPiPtDNdEta : ((axisMode == 2) ? hPiPtDNdY : hPiPt);
         TH2D *hRawP2D = (axisMode == 1) ? hPPtDNdEta : ((axisMode == 2) ? hPPtDNdY : hPPt);
         TH2D *hCorrK2D = (axisMode == 1) ? hKPtCorrectedDNdEta : ((axisMode == 2) ? hKPtCorrectedDNdY : hKPtCorrected);
         TH2D *hCorrPi2D = (axisMode == 1) ? hPiPtCorrectedDNdEta : ((axisMode == 2) ? hPiPtCorrectedDNdY : hPiPtCorrected);
         TH2D *hCorrP2D = (axisMode == 1) ? hPPtCorrectedDNdEta : ((axisMode == 2) ? hPPtCorrectedDNdY : hPPtCorrected);
         TH1D *hRawK1D = (axisMode == 1) ? hKDNdEta : ((axisMode == 2) ? hKDNdY : hK);
         TH1D *hRawPi1D = (axisMode == 1) ? hPiDNdEta : ((axisMode == 2) ? hPiDNdY : hPi);
         TH1D *hRawP1D = (axisMode == 1) ? hPDNdEta : ((axisMode == 2) ? hPDNdY : hP);
         TH1D *hCorrK1D = (axisMode == 1) ? hKCorrectedDNdEta : ((axisMode == 2) ? hKCorrectedDNdY : hKCorrected);
         TH1D *hCorrPi1D = (axisMode == 1) ? hPiCorrectedDNdEta : ((axisMode == 2) ? hPiCorrectedDNdY : hPiCorrected);
         TH1D *hCorrP1D = (axisMode == 1) ? hPCorrectedDNdEta : ((axisMode == 2) ? hPCorrectedDNdY : hPCorrected);
         const std::vector<double> *vKAsK = (axisMode == 1) ? &SumKAsKDNdEta : ((axisMode == 2) ? &SumKAsKDNdY : &SumKAsK);
         const std::vector<double> *vKAsPi = (axisMode == 1) ? &SumKAsPiDNdEta : ((axisMode == 2) ? &SumKAsPiDNdY : &SumKAsPi);
         const std::vector<double> *vKAsP = (axisMode == 1) ? &SumKAsPDNdEta : ((axisMode == 2) ? &SumKAsPDNdY : &SumKAsP);
         const std::vector<double> *vPiAsK = (axisMode == 1) ? &SumPiAsKDNdEta : ((axisMode == 2) ? &SumPiAsKDNdY : &SumPiAsK);
         const std::vector<double> *vPiAsPi = (axisMode == 1) ? &SumPiAsPiDNdEta : ((axisMode == 2) ? &SumPiAsPiDNdY : &SumPiAsPi);
         const std::vector<double> *vPiAsP = (axisMode == 1) ? &SumPiAsPDNdEta : ((axisMode == 2) ? &SumPiAsPDNdY : &SumPiAsP);
         const std::vector<double> *vPAsK = (axisMode == 1) ? &SumPAsKDNdEta : ((axisMode == 2) ? &SumPAsKDNdY : &SumPAsK);
         const std::vector<double> *vPAsPi = (axisMode == 1) ? &SumPAsPiDNdEta : ((axisMode == 2) ? &SumPAsPiDNdY : &SumPAsPi);
         const std::vector<double> *vPAsP = (axisMode == 1) ? &SumPAsPDNdEta : ((axisMode == 2) ? &SumPAsPDNdY : &SumPAsP);
         const std::vector<double> *vRecoK = (axisMode == 1) ? &SumRecoEffKDNdEta : ((axisMode == 2) ? &SumRecoEffKDNdY : &SumRecoEffK);
         const std::vector<double> *vRecoPi = (axisMode == 1) ? &SumRecoEffPiDNdEta : ((axisMode == 2) ? &SumRecoEffPiDNdY : &SumRecoEffPi);
         const std::vector<double> *vRecoP = (axisMode == 1) ? &SumRecoEffPDNdEta : ((axisMode == 2) ? &SumRecoEffPDNdY : &SumRecoEffP);
         const std::vector<double> *vGenK = (axisMode == 1) ? &SumGenEffKDNdEta : ((axisMode == 2) ? &SumGenEffKDNdY : &SumGenEffK);
         const std::vector<double> *vGenPi = (axisMode == 1) ? &SumGenEffPiDNdEta : ((axisMode == 2) ? &SumGenEffPiDNdY : &SumGenEffPi);
         const std::vector<double> *vGenP = (axisMode == 1) ? &SumGenEffPDNdEta : ((axisMode == 2) ? &SumGenEffPDNdY : &SumGenEffP);
         const std::vector<long long> *vCount = (axisMode == 1) ? &CountEffTracksDNdEta : ((axisMode == 2) ? &CountEffTracksDNdY : &CountEffTracks);
         const char *axisLabel = (axisMode == 1) ? "reco dNch/deta" : ((axisMode == 2) ? "reco dNch/dy" : "NchTag");

         for (int iNch = 1; iNch <= nNchBinsLocal; ++iNch)
         {
            for (int iPt = 1; iPt <= NPtBins; ++iPt)
            {
               const int idx = flatIndex(iNch, iPt);
               if ((*vCount)[idx] <= 0)
                  continue;

               const double den = static_cast<double>((*vCount)[idx]);
               const double eKAsK = (*vKAsK)[idx] / den;
               const double eKAsPi = (*vKAsPi)[idx] / den;
               const double eKAsP = (*vKAsP)[idx] / den;
               const double ePiAsK = (*vPiAsK)[idx] / den;
               const double ePiAsPi = (*vPiAsPi)[idx] / den;
               const double ePiAsP = (*vPiAsP)[idx] / den;
               const double ePAsK = (*vPAsK)[idx] / den;
               const double ePAsPi = (*vPAsPi)[idx] / den;
               const double ePAsP = (*vPAsP)[idx] / den;
               const double eRecoK = (*vRecoK)[idx] / den;
               const double eRecoPi = (*vRecoPi)[idx] / den;
               const double eRecoP = (*vRecoP)[idx] / den;
               const double eGenK = (*vGenK)[idx] / den;
               const double eGenPi = (*vGenPi)[idx] / den;
               const double eGenP = (*vGenP)[idx] / den;

               const double NKtag = hRawK2D->GetBinContent(iNch, iPt);
               const double NPiTag = hRawPi2D->GetBinContent(iNch, iPt);
               const double NPtag = hRawP2D->GetBinContent(iNch, iPt);

               double NtrueReco[3] = {0.0, 0.0, 0.0};
               double eNtrueReco[3] = {0.0, 0.0, 0.0};
               double Ntrue[3] = {0.0, 0.0, 0.0};
               double eNtrue[3] = {0.0, 0.0, 0.0};

               double recoMatch[3] = {
                  std::clamp(eRecoK, 0.0, 1.0),
                  std::clamp(eRecoPi, 0.0, 1.0),
                  std::clamp(eRecoP, 0.0, 1.0)
               };
               const double Y[3] = {NKtag * recoMatch[0], NPiTag * recoMatch[1], NPtag * recoMatch[2]};

               double Mmat[3][3];
               Mmat[0][0] = eKAsK;  Mmat[0][1] = ePiAsK;  Mmat[0][2] = ePAsK;
               Mmat[1][0] = eKAsPi; Mmat[1][1] = ePiAsPi; Mmat[1][2] = ePAsPi;
               Mmat[2][0] = eKAsP;  Mmat[2][1] = ePiAsP;  Mmat[2][2] = ePAsP;

               double Minv[3][3];
               double det = 0.0;
               if (!invert3x3(Mmat, Minv, det))
               {
                  cerr << "Warning: 3x3 tagging matrix near-singular in "
                       << axisLabel << " bin " << iNch << ", pT bin " << iPt
                       << " (det = " << det << "). Skipping correction for this bin."
                       << endl;
                  continue;
               }

               for (int s = 0; s < 3; ++s)
               {
                  NtrueReco[s] = Minv[s][0] * Y[0] + Minv[s][1] * Y[1] + Minv[s][2] * Y[2];
                  if (NtrueReco[s] < 0.0)
                     NtrueReco[s] = 0.0;
               }

               const double varY[3] = {
                  (Y[0] > 0.0 ? Y[0] : 1.0),
                  (Y[1] > 0.0 ? Y[1] : 1.0),
                  (Y[2] > 0.0 ? Y[2] : 1.0)
               };
               for (int s = 0; s < 3; ++s)
               {
                  double v = 0.0;
                  for (int i = 0; i < 3; ++i)
                     v += Minv[s][i] * Minv[s][i] * varY[i];
                  eNtrueReco[s] = (v > 0.0 ? std::sqrt(v) : 0.0);
               }

               double gMatch[3] = {
                  std::clamp(eGenK, 0.0, 1.0),
                  std::clamp(eGenPi, 0.0, 1.0),
                  std::clamp(eGenP, 0.0, 1.0)
               };
               for (int s = 0; s < 3; ++s)
               {
                  if (gMatch[s] > 1e-12)
                  {
                     Ntrue[s] = NtrueReco[s] / gMatch[s];
                     eNtrue[s] = eNtrueReco[s] / gMatch[s];
                  }
                  else
                  {
                     Ntrue[s] = 0.0;
                     eNtrue[s] = 0.0;
                  }
               }

               hCorrK2D->SetBinContent(iNch, iPt, Ntrue[0]);
               hCorrK2D->SetBinError(iNch, iPt, eNtrue[0]);
               hCorrPi2D->SetBinContent(iNch, iPt, Ntrue[1]);
               hCorrPi2D->SetBinError(iNch, iPt, eNtrue[1]);
               hCorrP2D->SetBinContent(iNch, iPt, Ntrue[2]);
               hCorrP2D->SetBinError(iNch, iPt, eNtrue[2]);
            }

            double sumKcorr = 0.0;
            double sumPicorr = 0.0;
            double sumPcorr = 0.0;
            double err2Kcorr = 0.0;
            double err2Picorr = 0.0;
            double err2Pcorr = 0.0;
            double sumKraw = 0.0;
            double sumPiraw = 0.0;
            double sumPraw = 0.0;
            double err2Kraw = 0.0;
            double err2Piraw = 0.0;
            double err2Praw = 0.0;

            for (int iPt = 1; iPt <= NPtBins; ++iPt)
            {
               const double vKc = hCorrK2D->GetBinContent(iNch, iPt);
               const double vPic = hCorrPi2D->GetBinContent(iNch, iPt);
               const double vPc = hCorrP2D->GetBinContent(iNch, iPt);
               const double eKc = hCorrK2D->GetBinError(iNch, iPt);
               const double ePic = hCorrPi2D->GetBinError(iNch, iPt);
               const double ePc = hCorrP2D->GetBinError(iNch, iPt);

               sumKcorr += vKc;
               sumPicorr += vPic;
               sumPcorr += vPc;
               err2Kcorr += eKc * eKc;
               err2Picorr += ePic * ePic;
               err2Pcorr += ePc * ePc;

               const double vKr = hRawK2D->GetBinContent(iNch, iPt);
               const double vPir = hRawPi2D->GetBinContent(iNch, iPt);
               const double vPr = hRawP2D->GetBinContent(iNch, iPt);
               const double eKr = hRawK2D->GetBinError(iNch, iPt);
               const double ePir = hRawPi2D->GetBinError(iNch, iPt);
               const double ePr = hRawP2D->GetBinError(iNch, iPt);

               sumKraw += vKr;
               sumPiraw += vPir;
               sumPraw += vPr;
               err2Kraw += eKr * eKr;
               err2Piraw += ePir * ePir;
               err2Praw += ePr * ePr;
            }

            hCorrK1D->SetBinContent(iNch, sumKcorr);
            hCorrK1D->SetBinError(iNch, std::sqrt(err2Kcorr));
            hCorrPi1D->SetBinContent(iNch, sumPicorr);
            hCorrPi1D->SetBinError(iNch, std::sqrt(err2Picorr));
            hCorrP1D->SetBinContent(iNch, sumPcorr);
            hCorrP1D->SetBinError(iNch, std::sqrt(err2Pcorr));
            hRawK1D->SetBinContent(iNch, sumKraw);
            hRawK1D->SetBinError(iNch, std::sqrt(err2Kraw));
            hRawPi1D->SetBinContent(iNch, sumPiraw);
            hRawPi1D->SetBinError(iNch, std::sqrt(err2Piraw));
            hRawP1D->SetBinContent(iNch, sumPraw);
            hRawP1D->SetBinError(iNch, std::sqrt(err2Praw));
         }
      };

      correctAxis(0);
      correctAxis(1);
      correctAxis(2);

      //-------------------------------------------------
      // Raw K/pi and p/pi vs Nch_tag from p_{T}-integrated spectra
      //-------------------------------------------------
      hKoverPi = (TH1D *)hK->Clone("hKoverPi");
      hPoverPi = (TH1D *)hP->Clone("hPoverPi");

      hKoverPi->SetTitle("K/#pi yield ratio vs N_{ch}^{tag};N_{ch}^{tag};K/#pi (reco, raw, p_{T}-integrated)");
      hPoverPi->SetTitle("p/#pi yield ratio vs N_{ch}^{tag};N_{ch}^{tag};p/#pi (reco, raw, p_{T}-integrated)");

      hKoverPi->Divide(hPi);
      hPoverPi->Divide(hPi);

      //-------------------------------------------------
      // Corrected K/pi and p/pi vs Nch_tag (from 1D yields)
      //-------------------------------------------------
      hKoverPiCorrected = (TH1D *)hKCorrected->Clone("hKoverPiCorrected");
      hKoverPiCorrected->SetTitle("K/#pi vs N_{ch}^{tag};N_{ch}^{tag};K/#pi (PID-corrected, p_{T}-integrated)");
      hKoverPiCorrected->Divide(hPiCorrected);

      hPoverPiCorrected = (TH1D *)hPCorrected->Clone("hPoverPiCorrected");
      hPoverPiCorrected->SetTitle("p/#pi vs N_{ch}^{tag};N_{ch}^{tag};p/#pi (PID-corrected, p_{T}-integrated)");
      hPoverPiCorrected->Divide(hPiCorrected);

      hKoverPiDNdEta = (TH1D *)hKDNdEta->Clone("hKoverPiDNdEta");
      hPoverPiDNdEta = (TH1D *)hPDNdEta->Clone("hPoverPiDNdEta");
      hKoverPiDNdEta->SetTitle("K/#pi yield ratio vs reco dN_{ch}/d#eta;dN_{ch}/d#eta (reco, |#eta|<0.5);K/#pi (reco, raw, p_{T}-integrated)");
      hPoverPiDNdEta->SetTitle("p/#pi yield ratio vs reco dN_{ch}/d#eta;dN_{ch}/d#eta (reco, |#eta|<0.5);p/#pi (reco, raw, p_{T}-integrated)");
      hKoverPiDNdEta->Divide(hPiDNdEta);
      hPoverPiDNdEta->Divide(hPiDNdEta);

      hKoverPiCorrectedDNdEta = (TH1D *)hKCorrectedDNdEta->Clone("hKoverPiCorrectedDNdEta");
      hPoverPiCorrectedDNdEta = (TH1D *)hPCorrectedDNdEta->Clone("hPoverPiCorrectedDNdEta");
      hKoverPiCorrectedDNdEta->SetTitle("K/#pi vs reco dN_{ch}/d#eta;dN_{ch}/d#eta (reco, |#eta|<0.5);K/#pi (PID-corrected, p_{T}-integrated)");
      hPoverPiCorrectedDNdEta->SetTitle("p/#pi vs reco dN_{ch}/d#eta;dN_{ch}/d#eta (reco, |#eta|<0.5);p/#pi (PID-corrected, p_{T}-integrated)");
      hKoverPiCorrectedDNdEta->Divide(hPiCorrectedDNdEta);
      hPoverPiCorrectedDNdEta->Divide(hPiCorrectedDNdEta);

      hKoverPiDNdY = (TH1D *)hKDNdY->Clone("hKoverPiDNdY");
      hPoverPiDNdY = (TH1D *)hPDNdY->Clone("hPoverPiDNdY");
      hKoverPiDNdY->SetTitle("K/#pi yield ratio vs reco dN_{ch}/dy;dN_{ch}/dy (reco, thrust axis, |y_{T}|<0.5);K/#pi (reco, raw, p_{T}-integrated)");
      hPoverPiDNdY->SetTitle("p/#pi yield ratio vs reco dN_{ch}/dy;dN_{ch}/dy (reco, thrust axis, |y_{T}|<0.5);p/#pi (reco, raw, p_{T}-integrated)");
      hKoverPiDNdY->Divide(hPiDNdY);
      hPoverPiDNdY->Divide(hPiDNdY);

      hKoverPiCorrectedDNdY = (TH1D *)hKCorrectedDNdY->Clone("hKoverPiCorrectedDNdY");
      hPoverPiCorrectedDNdY = (TH1D *)hPCorrectedDNdY->Clone("hPoverPiCorrectedDNdY");
      hKoverPiCorrectedDNdY->SetTitle("K/#pi vs reco dN_{ch}/dy;dN_{ch}/dy (reco, thrust axis, |y_{T}|<0.5);K/#pi (PID-corrected, p_{T}-integrated)");
      hPoverPiCorrectedDNdY->SetTitle("p/#pi vs reco dN_{ch}/dy;dN_{ch}/dy (reco, thrust axis, |y_{T}|<0.5);p/#pi (PID-corrected, p_{T}-integrated)");
      hKoverPiCorrectedDNdY->Divide(hPiCorrectedDNdY);
      hPoverPiCorrectedDNdY->Divide(hPiCorrectedDNdY);

      if (!par.IsGen && NPIDPassTagTracks > 0)
      {
         const double frac = static_cast<double>(NPIDTieTracks) / static_cast<double>(NPIDPassTagTracks);
         cout << "PID overlap diagnostics:" << endl;
         cout << "  pass-tag tracks = " << NPIDPassTagTracks << endl;
         cout << "  tie tracks      = " << NPIDTieTracks
              << " (" << 100.0 * frac << "%)" << endl;
      }
   }

   void writeHistograms()
   {
      if (outf == nullptr)
         return;

      outf->cd();

      // 1D histograms vs Nch_tag
      smartWrite(hK);
      smartWrite(hPi);
      smartWrite(hP);
      smartWrite(hKoverPi);
      smartWrite(hPoverPi);
      smartWrite(hKDNdEta);
      smartWrite(hPiDNdEta);
      smartWrite(hPDNdEta);
      smartWrite(hKoverPiDNdEta);
      smartWrite(hPoverPiDNdEta);
      smartWrite(hKDNdY);
      smartWrite(hPiDNdY);
      smartWrite(hPDNdY);
      smartWrite(hKoverPiDNdY);
      smartWrite(hPoverPiDNdY);
      smartWrite(hKPt);
      smartWrite(hPiPt);
      smartWrite(hPPt);
      smartWrite(hKPtDNdEta);
      smartWrite(hPiPtDNdEta);
      smartWrite(hPPtDNdEta);
      smartWrite(hKPtDNdY);
      smartWrite(hPiPtDNdY);
      smartWrite(hPPtDNdY);

      if (!par.IsGen)
      {
         smartWrite(hKCorrected);
         smartWrite(hPiCorrected);
         smartWrite(hPCorrected);
         smartWrite(hKoverPiCorrected);
         smartWrite(hPoverPiCorrected);
         smartWrite(hKCorrectedDNdEta);
         smartWrite(hPiCorrectedDNdEta);
         smartWrite(hPCorrectedDNdEta);
         smartWrite(hKoverPiCorrectedDNdEta);
         smartWrite(hPoverPiCorrectedDNdEta);
         smartWrite(hKCorrectedDNdY);
         smartWrite(hPiCorrectedDNdY);
         smartWrite(hPCorrectedDNdY);
         smartWrite(hKoverPiCorrectedDNdY);
         smartWrite(hPoverPiCorrectedDNdY);

         // 2D pT spectra
         smartWrite(hUPt);
         smartWrite(hKPtCorrected);
         smartWrite(hPiPtCorrected);
         smartWrite(hPPtCorrected);
         smartWrite(hUPtDNdEta);
         smartWrite(hKPtCorrectedDNdEta);
         smartWrite(hPiPtCorrectedDNdEta);
         smartWrite(hPPtCorrectedDNdEta);
         smartWrite(hUPtDNdY);
         smartWrite(hKPtCorrectedDNdY);
         smartWrite(hPiPtCorrectedDNdY);
         smartWrite(hPPtCorrectedDNdY);

         // MC-only histograms used by Ntag unfolding (empty for data)
         smartWrite(hNtagResponse);
         smartWrite(hNtagResponseK);
         smartWrite(hNtagResponsePi);
         smartWrite(hNtagResponseP);
         smartWrite(hNtagTrue);
         smartWrite(hNtagReco);
         smartWrite(hKTrueNtag);
         smartWrite(hPiTrueNtag);
         smartWrite(hPTrueNtag);
         smartWrite(hDNdEtaResponse);
         smartWrite(hDNdEtaResponseK);
         smartWrite(hDNdEtaResponsePi);
         smartWrite(hDNdEtaResponseP);
         smartWrite(hDNdEtaTrue);
         smartWrite(hDNdEtaReco);
         smartWrite(hKTruedNdEta);
         smartWrite(hPiTruedNdEta);
         smartWrite(hPTruedNdEta);
         smartWrite(hDNdYResponse);
         smartWrite(hDNdYResponseK);
         smartWrite(hDNdYResponsePi);
         smartWrite(hDNdYResponseP);
         smartWrite(hDNdYTrue);
         smartWrite(hDNdYReco);
         smartWrite(hKTruedNdY);
         smartWrite(hPiTruedNdY);
         smartWrite(hPTruedNdY);
      }

      // Raw K/π canvas
      TCanvas c1("c1", "K/pi vs NchTag (raw)", 800, 600);
      hKoverPi->SetMarkerStyle(20);
      hKoverPi->SetMarkerSize(1.0);
      hKoverPi->Draw("E1");
      c1.Write();

      // Raw p/π canvas
      TCanvas c1b("c1b", "p/pi vs NchTag (raw)", 800, 600);
      hPoverPi->SetMarkerStyle(20);
      hPoverPi->SetMarkerSize(1.0);
      hPoverPi->Draw("E1");
      c1b.Write();

      // PID-corrected K/π and p/π canvases (reco mode only)
      if (!par.IsGen && hKoverPiCorrected != nullptr && hPoverPiCorrected != nullptr)
      {
         TCanvas c2("c2", "K/pi vs NchTag (PID-corrected)", 800, 600);
         hKoverPiCorrected->SetMarkerStyle(21);
         hKoverPiCorrected->SetMarkerSize(1.0);
         hKoverPiCorrected->Draw("E1");
         c2.Write();

         TCanvas c3("c3", "p/pi vs NchTag (PID-corrected)", 800, 600);
         hPoverPiCorrected->SetMarkerStyle(21);
         hPoverPiCorrected->SetMarkerSize(1.0);
         hPoverPiCorrected->Draw("E1");
         c3.Write();
      }
   }
};

//============================================================
// Main analysis
//============================================================
int main(int argc, char *argv[])
{
   // if (printHelpMessage(argc, argv))
   //    return 0;

   CommandLine CL(argc, argv);
   KtoPiParameters par;

   // Basic I/O
   par.input  = CL.Get("Input",  par.input);
   par.output = CL.Get("Output", par.output);

   // Physics / binning parameters
   par.MaxNchTag = CL.GetInt   ("MaxNchTag", par.MaxNchTag);
   par.MaxEvents = CL.GetInt   ("MaxEvents", par.MaxEvents);
   par.EcmRef    = CL.GetDouble("EcmRef",    par.EcmRef);
   par.MinNch    = CL.GetInt   ("MinNch",    par.MinNch);

   double MinThetaDeg = CL.GetDouble("MinThetaDeg", 30.0);
   double MaxThetaDeg = CL.GetDouble("MaxThetaDeg", 150.0);
   par.MinTheta       = MinThetaDeg * TMath::Pi() / 180.0;
   par.MaxTheta       = MaxThetaDeg * TMath::Pi() / 180.0;

   // Generator-level flag: IsGen=true/false or 1/0 or yes/no
   std::string isGenStr = CL.Get("IsGen", std::string("false"));
   if (isGenStr == "1" || isGenStr == "true" || isGenStr == "True" ||
       isGenStr == "TRUE" || isGenStr == "yes" || isGenStr == "Yes" ||
       isGenStr == "YES")
      par.IsGen = true;
   else
      par.IsGen = false;

   std::string useTruthMatrixStr = CL.Get("UseMCTruthMatrix", std::string("false"));
   if (useTruthMatrixStr == "1" || useTruthMatrixStr == "true" || useTruthMatrixStr == "True" ||
       useTruthMatrixStr == "TRUE" || useTruthMatrixStr == "yes" || useTruthMatrixStr == "Yes" ||
       useTruthMatrixStr == "YES")
      par.UseMCTruthMatrix = true;
   else
      par.UseMCTruthMatrix = false;

   std::string useCentralEtaNtagStr = CL.Get("UseCentralEtaNtag", std::string("false"));
   if (useCentralEtaNtagStr == "1" || useCentralEtaNtagStr == "true" || useCentralEtaNtagStr == "True" ||
       useCentralEtaNtagStr == "TRUE" || useCentralEtaNtagStr == "yes" || useCentralEtaNtagStr == "Yes" ||
       useCentralEtaNtagStr == "YES")
      par.UseCentralEtaNtag = true;
   else
      par.UseCentralEtaNtag = false;

   std::string usePassAllSelectionStr = CL.Get("UsePassAllSelection", std::string("true"));
   if (usePassAllSelectionStr == "1" || usePassAllSelectionStr == "true" || usePassAllSelectionStr == "True" ||
       usePassAllSelectionStr == "TRUE" || usePassAllSelectionStr == "yes" || usePassAllSelectionStr == "Yes" ||
       usePassAllSelectionStr == "YES")
      par.UsePassAllSelection = true;
   else
      par.UsePassAllSelection = false;

   std::string usePIDFiducialStr = CL.Get("UsePIDFiducial", std::string("true"));
   if (usePIDFiducialStr == "1" || usePIDFiducialStr == "true" || usePIDFiducialStr == "True" ||
       usePIDFiducialStr == "TRUE" || usePIDFiducialStr == "yes" || usePIDFiducialStr == "Yes" ||
       usePIDFiducialStr == "YES")
      par.UsePIDFiducial = true;
   else
      par.UsePIDFiducial = false;
   par.PIDTrackAbsCosMin = CL.GetDouble("PIDTrackAbsCosMin", par.PIDTrackAbsCosMin);
   par.PIDTrackAbsCosMax = CL.GetDouble("PIDTrackAbsCosMax", par.PIDTrackAbsCosMax);
   std::string pidTieMode = CL.Get("PIDTieMode", std::string("legacy"));
   if (pidTieMode == "untag" || pidTieMode == "Untag" || pidTieMode == "UNTAG")
      par.PIDTieMode = 1;
   else
      par.PIDTieMode = 0;
   std::string pidObservationMode = CL.Get("PIDObservationMode", std::string("exclusive"));
   std::string allowDuplicatePID = CL.Get("AllowDuplicatePIDCandidates", std::string(""));
   if (!allowDuplicatePID.empty())
   {
      if (allowDuplicatePID == "1" || allowDuplicatePID == "true" || allowDuplicatePID == "True" ||
          allowDuplicatePID == "TRUE" || allowDuplicatePID == "yes" || allowDuplicatePID == "Yes" ||
          allowDuplicatePID == "YES")
         par.UseInclusivePIDObservation = true;
      else
         par.UseInclusivePIDObservation = false;
   }
   else if (pidObservationMode == "inclusive" || pidObservationMode == "Inclusive" ||
            pidObservationMode == "INCLUSIVE" || pidObservationMode == "duplicate" ||
            pidObservationMode == "Duplicate" || pidObservationMode == "DUPLICATE")
   {
      par.UseInclusivePIDObservation = true;
   }
   else
   {
      par.UseInclusivePIDObservation = false;
   }

   // pT binning:
   //   Option 1 (default): uniform bins
   //     NPtBins, PtMin, PtMax
   //   Option 2: explicit edges
   //     PtBinEdges=0.2,0.4,0.6,1.0,2.0,3.0,5.0
   par.NPtBins = CL.GetInt   ("NPtBins", par.NPtBins);
   par.PtMin   = CL.GetDouble("PtMin",   par.PtMin);
   par.PtMax   = CL.GetDouble("PtMax",   par.PtMax);
   par.NtagPtMin = CL.GetDouble("NtagPtMin", par.NtagPtMin);

   const std::string ptEdgesStr = CL.Get("PtBinEdges", std::string(""));
   if (!ptEdgesStr.empty())
   {
      par.PtBinEdges = ParseDoubleList(ptEdgesStr);
   }

   cout << "Running KtoPiAnalysis with parameters:" << endl;
   cout << "  Input       = " << par.input      << endl;
   cout << "  Output      = " << par.output     << endl;
   cout << "  MaxNchTag   = " << par.MaxNchTag  << endl;
   cout << "  MaxEvents   = " << par.MaxEvents  << endl;
   cout << "  EcmRef      = " << par.EcmRef     << endl;
   cout << "  MinNch      = " << par.MinNch     << endl;
   cout << "  MinThetaDeg = " << MinThetaDeg    << endl;
   cout << "  MaxThetaDeg = " << MaxThetaDeg    << endl;
   cout << "  IsGen       = " << (par.IsGen ? "true" : "false") << endl;
   cout << "  UseMCTruthMatrix = " << (par.UseMCTruthMatrix ? "true" : "false") << endl;
   cout << "  UsePassAllSelection = " << (par.UsePassAllSelection ? "true" : "false") << endl;
   cout << "  UseCentralEtaNtag = " << (par.UseCentralEtaNtag ? "true" : "false") << endl;
   cout << "  UsePIDFiducial = " << (par.UsePIDFiducial ? "true" : "false") << endl;
   cout << "  PIDTrackAbsCosMin = " << par.PIDTrackAbsCosMin << endl;
   cout << "  PIDTrackAbsCosMax = " << par.PIDTrackAbsCosMax << endl;
   cout << "  PIDObservationMode = " << (par.UseInclusivePIDObservation ? "inclusive" : "exclusive") << endl;
   cout << "  PIDTieMode = " << (par.PIDTieMode == 1 ? "untag" : "legacy") << endl;
   cout << "  NtagPtMin   = " << par.NtagPtMin << endl;

   if (!par.PtBinEdges.empty())
   {
      cout << "  pT binning  = custom edges: ";
      for (size_t i = 0; i < par.PtBinEdges.size(); ++i)
      {
         cout << par.PtBinEdges[i];
         if (i + 1 < par.PtBinEdges.size()) cout << ", ";
      }
      cout << endl;
   }
   else
   {
      cout << "  pT binning  = uniform: NPtBins=" << par.NPtBins
           << ", PtMin=" << par.PtMin
           << ", PtMax=" << par.PtMax << endl;
   }

   KtoPiAnalyzer analyzer(par);
   analyzer.analyze();
   analyzer.writeHistograms();

   cout << "Done. Output written to " << par.output << endl;
   return 0;
}
