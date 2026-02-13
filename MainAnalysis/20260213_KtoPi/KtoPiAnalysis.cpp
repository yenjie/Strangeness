//============================================================
// KtoPiAnalysis.cpp
//
// Computes kaon and pion yields vs Nch_tag and K/pi ratio
// using the Strangeness tree messenger.
//
// Intended usage (following PhysicsZHadronEEC style):
//
//   ./KtoPiAnalysis \
//     Input=sample/Strangeness/merged_mc_v2.root \
//     Output=output/KtoPi.root \
//     MaxNchTag=60 MaxEvents=-1
//
// New option:
//
//   IsGen=true  -> count K/π at generator level using GenID (PDG ID)
//   IsGen=false -> (default) use reconstructed PID info
//
// In reco mode (IsGen=false), we now also apply a simple PID
// matrix correction using the nine RecoEfficiency*As* arrays:
//   - KAsK, KAsPi, PiAsK, PiAsPi
// are used to correct the kaon and pion yields for
// tagging efficiency and fake rates (2×2 K/π sub-matrix).
//============================================================

#include <iostream>
#include <string>

// ROOT
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
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

   bool   IsGen;       // if true, count K/π at generator level

   KtoPiParameters()
      : input("sample/Strangeness/merged_mc_v2.root")
      , output("output/KtoPi.root")
      , MaxNchTag(60)
      , MaxEvents(-1)
      , EcmRef(91.2)
      , MinNch(7)
      , MinTheta(30.0 * TMath::Pi() / 180.0)
      , MaxTheta(150.0 * TMath::Pi() / 180.0)
      , IsGen(false)
   {
   }
};

//============================================================
// Analyzer class
//============================================================
class KtoPiAnalyzer
{
public:
   TFile *inf;
   TFile *outf;
   StrangenessTreeMessenger *M;

   // Raw (uncorrected) yields
   TH1D *hK;
   TH1D *hPi;
   TH1D *hKoverPi;

   // PID–corrected yields (efficiency + fake-rate correction)
   TH1D *hKCorrected;
   TH1D *hPiCorrected;
   TH1D *hKoverPiCorrected;

   // Parameters
   KtoPiParameters par;

   // Accumulators to build a global 2×2 PID matrix (K/π only)
   double sumKAsK;
   double sumKAsPi;
   double sumPiAsK;
   double sumPiAsPi;
   long long countEffTracks;

public:
   KtoPiAnalyzer(const KtoPiParameters &apar)
      : inf(nullptr)
      , outf(nullptr)
      , M(nullptr)
      , hK(nullptr)
      , hPi(nullptr)
      , hKoverPi(nullptr)
      , hKCorrected(nullptr)
      , hPiCorrected(nullptr)
      , hKoverPiCorrected(nullptr)
      , par(apar)
      , sumKAsK(0.0)
      , sumKAsPi(0.0)
      , sumPiAsK(0.0)
      , sumPiAsPi(0.0)
      , countEffTracks(0)
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

      // Book histograms
      const int maxNchTag = par.MaxNchTag;
      const int nbins     = maxNchTag / 4 + 1;   // same as your macro

      hK = new TH1D("hK",
                    "Kaon candidates vs N_{ch}^{tag};N_{ch}^{tag};Yield (sum over events)",
                    nbins, -0.5, maxNchTag + 0.5);

      hPi = (TH1D *)hK->Clone("hPi");
      hPi->SetTitle("Pion candidates vs N_{ch}^{tag};N_{ch}^{tag};Yield (sum over events)");

      hK->Sumw2();
      hPi->Sumw2();

      // Corrected versions – initially empty, filled after matrix inversion
      hKCorrected = (TH1D *)hK->Clone("hKCorrected");
      hKCorrected->SetTitle("PID-corrected K yield vs N_{ch}^{tag};N_{ch}^{tag};Corrected K yield");
      hKCorrected->Reset();
      hKCorrected->Sumw2();

      hPiCorrected = (TH1D *)hK->Clone("hPiCorrected");
      hPiCorrected->SetTitle("PID-corrected #pi yield vs N_{ch}^{tag};N_{ch}^{tag};Corrected #pi yield");
      hPiCorrected->Reset();
      hPiCorrected->Sumw2();

      hKoverPi          = nullptr;  // will be built after event loop
      hKoverPiCorrected = nullptr;  // will be built after matrix correction
   }

   ~KtoPiAnalyzer()
   {
      // Clean up histograms (the file should own them after writing, but be safe)
      delete hK;
      delete hPi;
      delete hKoverPi;
      delete hKCorrected;
      delete hPiCorrected;
      delete hKoverPiCorrected;

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

   void analyze()
   {
      if (M == nullptr || inf == nullptr || outf == nullptr)
         return;

      // Event loop
      long long nEntries = M->GetEntries();
      if (par.MaxEvents > 0 && par.MaxEvents < nEntries)
         nEntries = par.MaxEvents;

      cout << "Total entries to process: " << nEntries << endl;

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

         // Optionally prepare NGen if we are doing generator-level counting
         int ngen = 0;
         if (par.IsGen)
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

         // Sum(RecoE)
         double sumRecoE = 0.0;
         for (int i = 0; i < nreco; ++i)
            sumRecoE += M->RecoE[i];

         // Sum(RecoE)/EcmRef > 0.5
         if (sumRecoE / EcmRef <= 0.5)
            continue;

         // Nch >= MinNch
         if (M->Nch < MinNch)
            continue;

         // 30° < acos(ThrustZ) < 150°
         double theta = std::acos(M->ThrustZ);
         if (theta <= MinTheta)
            continue;
         if (theta >= MaxTheta)
            continue;

         //-------------------------
         // Build NchTag, nK, nPi
         //-------------------------
         int NchTag = 0;
         int nK     = 0;
         int nPi    = 0;

         if (par.IsGen)
         {
            // NchTag is still defined by reconstructed tagged tracks
            for (int i = 0; i < nreco; ++i)
            {
               bool isKaonTag   = (M->RecoPIDKaon[i] >= 2);
               bool isPionTag   = (M->RecoPIDPion[i] >= 2);
               bool isProtonTag = (M->RecoPIDProton[i] >= 2);

               // NchTag = Sum(RecoPIDKaon>=2 || RecoPIDPion>=2 || RecoPIDProton>=2)
               if (isKaonTag || isPionTag || isProtonTag)
                  ++NchTag;
            }

            // Count generator-level kaons and pions using PDG ID
            for (int i = 0; i < ngen; ++i)
            {
               const int pdg    = static_cast<int>(M->GenID[i]);
               const int absPdg = (pdg >= 0 ? pdg : -pdg);

               // Charged kaons: K+, K-
               if (absPdg == 321)
                  ++nK;

               // Charged pions: pi+, pi-
               if (absPdg == 211)
                  ++nPi;
            }
         }
         else
         {
            // Reco-based PID counting (to be corrected later)
            for (int i = 0; i < nreco; ++i)
            {
               bool isKaonTag   = (M->RecoPIDKaon[i] >= 2);
               bool isPionTag   = (M->RecoPIDPion[i] >= 2);
               bool isProtonTag = (M->RecoPIDProton[i] >= 2);

               // NchTag = Sum(RecoPIDKaon>=2 || RecoPIDPion>=2 || RecoPIDProton>=2)
               if (isKaonTag || isPionTag || isProtonTag)
                  ++NchTag;

               if (isKaonTag)
                  ++nK;
               if (isPionTag)
                  ++nPi;

               // Accumulate PID efficiency / fake-rate info for global 2×2 matrix
               //
               // The RecoEfficiencyXAsY arrays are defined as:
               //   e.g. PiAsPi: probability for a *true pion* in this (pT, theta)
               //   bin to be tagged as a pion.
               //
               // They come from MC calibration and are stored per track as a
               // function of its kinematics.  We don't know the true species
               // in data, but we can average these calibration functions over
               // all charged tracks to build an effective K/π matrix.
               //
               // To keep things simple and general for both MC and future data,
               // we do *not* use truth here; we just average over all charged
               // tracks that we consider "taggable".
               if (M->RecoCharge[i] != 0.0)
               {
                  sumKAsK   += M->RecoEfficiencyKAsK[i];
                  sumKAsPi  += M->RecoEfficiencyKAsPi[i];
                  sumPiAsK  += M->RecoEfficiencyPiAsK[i];
                  sumPiAsPi += M->RecoEfficiencyPiAsPi[i];
                  ++countEffTracks;
               }
            }
         }

         // Put overflow NchTag into the last bin
         if (NchTag > par.MaxNchTag)
            NchTag = par.MaxNchTag;

         // Fill event-wise yields (raw)
         hK->Fill(NchTag, nK);
         hPi->Fill(NchTag, nPi);
      }

      cout << endl << "Event loop finished." << endl;

      //-------------------------
      // Update titles depending on mode
      //-------------------------
      if (par.IsGen)
      {
         hK->SetTitle ("Generator-level kaons vs N_{ch}^{tag};N_{ch}^{tag};N_{K}^{gen}");
         hPi->SetTitle("Generator-level pions vs N_{ch}^{tag};N_{ch}^{tag};N_{#pi}^{gen}");
      }
      else
      {
         hK->SetTitle ("Kaon candidates vs N_{ch}^{tag};N_{ch}^{tag};Yield (sum over events)");
         hPi->SetTitle("Pion candidates vs N_{ch}^{tag};N_{ch}^{tag};Yield (sum over events)");
      }

      //-------------------------
      // Build raw K/pi ratio histogram
      //-------------------------
      hKoverPi = (TH1D *)hK->Clone("hKoverPi");
      if (par.IsGen)
         hKoverPi->SetTitle("Generator-level K/#pi yield ratio vs N_{ch}^{tag};N_{ch}^{tag};K/#pi (gen)");
      else
         hKoverPi->SetTitle("K/#pi yield ratio vs N_{ch}^{tag};N_{ch}^{tag};K/#pi (reco)");

      hKoverPi->Divide(hPi);

      //-------------------------
      // Build PID-corrected yields & ratio (reco mode only)
      //-------------------------
      if (!par.IsGen && countEffTracks > 0)
      {
         const double eKAsK_avg   = sumKAsK   / countEffTracks;
         const double eKAsPi_avg  = sumKAsPi  / countEffTracks;
         const double ePiAsK_avg  = sumPiAsK  / countEffTracks;
         const double ePiAsPi_avg = sumPiAsPi / countEffTracks;

         cout << "Average K/π PID matrix (rows = tag K,tag π; cols = true K,true π)" << endl;
         cout << "  [tagK]  KAsK=" << eKAsK_avg  << "   PiAsK=" << ePiAsK_avg  << endl;
         cout << "  [tagPi] KAsPi=" << eKAsPi_avg << "   PiAsPi=" << ePiAsPi_avg << endl;

         // 2×2 matrix:
         //   [ N(tag K) ]   [ eKAsK   ePiAsK  ] [ N_true(K) ]
         //   [ N(tag π) ] = [ eKAsPi  ePiAsPi ] [ N_true(π) ]
         //
         // We invert this to estimate N_true(K), N_true(π) in each Nch_tag bin.
         const double det = eKAsK_avg * ePiAsPi_avg - ePiAsK_avg * eKAsPi_avg;

         if (TMath::Abs(det) < 1.0e-8)
         {
            cerr << "Warning: PID 2x2 K/π matrix determinant is tiny (" << det
                 << "). Skipping efficiency/fake-rate correction." << endl;
         }
         else
         {
            const int nbins = hK->GetNbinsX();

            for (int ib = 1; ib <= nbins; ++ib)
            {
               const double NKtag  = hK->GetBinContent(ib);
               const double NPiTag = hPi->GetBinContent(ib);

               // Matrix inversion:
               //
               // N_true(K) = ( ePiAsPi * NKtag - ePiAsK * NPiTag ) / det
               // N_true(π) = ( -eKAsPi * NKtag + eKAsK * NPiTag ) / det
               double NKtrue  = ( ePiAsPi_avg * NKtag - ePiAsK_avg * NPiTag ) / det;
               double NPItrue = ( -eKAsPi_avg * NKtag + eKAsK_avg * NPiTag ) / det;

               // Guard against small negative values from fluctuations / approximations
               if (NKtrue  < 0.0) NKtrue  = 0.0;
               if (NPItrue < 0.0) NPItrue = 0.0;

               hKCorrected->SetBinContent(ib, NKtrue);
               hPiCorrected->SetBinContent(ib, NPItrue);

               // Simple error propagation: same linear transformation applied
               const double eNKtag  = hK->GetBinError(ib);
               const double eNPitag = hPi->GetBinError(ib);

               const double eNKtrue = TMath::Sqrt(
                  TMath::Power(ePiAsPi_avg * eNKtag / det, 2) +
                  TMath::Power(ePiAsK_avg  * eNPitag / det, 2));

               const double eNPItrue = TMath::Sqrt(
                  TMath::Power(eKAsPi_avg * eNKtag / det, 2) +
                  TMath::Power(eKAsK_avg  * eNPitag / det, 2));

               hKCorrected->SetBinError(ib, eNKtrue);
               hPiCorrected->SetBinError(ib, eNPItrue);
            }

            hKoverPiCorrected = (TH1D *)hKCorrected->Clone("hKoverPiCorrected");
            hKoverPiCorrected->SetTitle("K/#pi vs N_{ch}^{tag};N_{ch}^{tag};K/#pi (PID-corrected)");
            hKoverPiCorrected->Divide(hPiCorrected);
         }
      }
      else if (!par.IsGen)
      {
         cerr << "Warning: no tracks accumulated for efficiency calibration; "
              << "PID-corrected histograms will remain empty." << endl;
      }
   }

   void writeHistograms()
   {
      if (outf == nullptr)
         return;

      outf->cd();

      smartWrite(hK);
      smartWrite(hPi);
      smartWrite(hKoverPi);

      if (!par.IsGen)
      {
         smartWrite(hKCorrected);
         smartWrite(hPiCorrected);
         smartWrite(hKoverPiCorrected);
      }

      // Raw K/π canvas
      TCanvas c1("c1", "K/pi vs NchTag (raw)", 800, 600);
      hKoverPi->SetMarkerStyle(20);
      hKoverPi->SetMarkerSize(1.0);
      hKoverPi->Draw("E1");
      c1.Write();  // store in ROOT file

      // PID-corrected K/π canvas (reco mode only)
      if (!par.IsGen && hKoverPiCorrected != nullptr)
      {
         TCanvas c2("c2", "K/pi vs NchTag (PID-corrected)", 800, 600);
         hKoverPiCorrected->SetMarkerStyle(21);
         hKoverPiCorrected->SetMarkerSize(1.0);
         hKoverPiCorrected->Draw("E1");
         c2.Write();
      }
   }
};

//============================================================
// Main analysis
//============================================================
int main(int argc, char *argv[])
{
   // Use the same help-message pattern as CorrelationAnalysis.cpp
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

   // Generator-level flag: IsGen=true/false or 1/0 or yes/no (case-insensitive)
   std::string isGenStr = CL.Get("IsGen", std::string("false"));
   if (isGenStr == "1" || isGenStr == "true" || isGenStr == "True" ||
       isGenStr == "TRUE" || isGenStr == "yes" || isGenStr == "Yes" ||
       isGenStr == "YES")
   {
      par.IsGen = true;
   }
   else
   {
      par.IsGen = false;
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

   KtoPiAnalyzer analyzer(par);
   analyzer.analyze();
   analyzer.writeHistograms();

   cout << "Done. Output written to " << par.output << endl;
   return 0;
}
