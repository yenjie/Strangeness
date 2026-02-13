//============================================================//
// KtoPiAnalysis.cpp
//
// Computes kaon and pion yields vs Nch_tag and K/pi ratio
// using the Strangeness tree messenger.
//
// Intended usage (following PhysicsZHadronEEC style):
//
//   ./KtoPiAnalysis \
//      Input=sample/Strangeness/merged_mc_v2.root \
//      Output=output/KtoPi.root \
//      MaxNchTag=60 MaxEvents=-1
//
//============================================================//

#include <iostream>
#include <cmath>
#include <string>
// ROOT
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TNtuple.h"

// Project common code (you may need to copy these from PhysicsZHadronEEC)
#include "utilities.h"      // smartWrite, etc.
#include "helpMessage.h"    // printHelpMessage()
#include "CommandLine.h"    // CommandLine parser
#include "ProgressBar.h"    // nice progress bar

// Strangeness tree messenger
#include "StrangenessMessenger.h"

using namespace std;

//============================================================//
// Simple parameter container for this analysis
//============================================================//

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

   KtoPiParameters()
      : input("sample/Strangeness/merged_mc_v2.root")
      , output("output/KtoPi.root")
      , MaxNchTag(60)
      , MaxEvents(-1)
      , EcmRef(91.2)
      , MinNch(7)
      , MinTheta(30.0 * TMath::Pi() / 180.0)
      , MaxTheta(150.0 * TMath::Pi() / 180.0)
   {
   }
};

//============================================================//
// Analyzer class 
//============================================================//

class KtoPiAnalyzer
{
public:
   TFile *inf;
   TFile *outf;

   StrangenessTreeMessenger *M;

   TH1D *hK;
   TH1D *hPi;
   TH1D *hKoverPi;

   KtoPiParameters par;

public:
   KtoPiAnalyzer(const KtoPiParameters &apar)
      : inf(nullptr)
      , outf(nullptr)
      , M(nullptr)
      , hK(nullptr)
      , hPi(nullptr)
      , hKoverPi(nullptr)
      , par(apar)
   {
      // Open input
      inf = new TFile(par.input.c_str());
      if (inf == nullptr || inf->IsZombie()) {
         cerr << "Error: cannot open input file '" << par.input << "'" << endl;
         return;
      }

      // Attach messenger to tree "Tree"
      M = new StrangenessTreeMessenger(*inf, std::string("Tree"));

      // Open output
      outf = new TFile(par.output.c_str(), "RECREATE");
      if (outf == nullptr || outf->IsZombie()) {
         cerr << "Error: cannot create output file '" << par.output << "'" << endl;
         return;
      }
      outf->cd();

      // Book histograms
      const int maxNchTag = par.MaxNchTag;
      const int nbins     = maxNchTag / 4 + 1;  // same as your macro

      hK = new TH1D("hK",
         "Kaon and pion yields vs N_{ch}^{tag};N_{ch}^{tag};Yield (sum over events)",
         nbins, -0.5, maxNchTag + 0.5);

      hPi = (TH1D*)hK->Clone("hPi");
      hPi->SetTitle("Pion candidates vs N_{ch}^{tag}");

      hK->Sumw2();
      hPi->Sumw2();

      hKoverPi = nullptr;  // will be built after the event loop
   }

   ~KtoPiAnalyzer()
   {
      // Clean up histograms (the file should own them after writing, but be safe)
      delete hK;
      delete hPi;
      delete hKoverPi;

      if (inf) {
         inf->Close();
         delete inf;
      }

      if (outf) {
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

      const double EcmRef = par.EcmRef;
      const int    MinNch = par.MinNch;
      const double MinTheta = par.MinTheta;
      const double MaxTheta = par.MaxTheta;

      for (long long ievt = 0; ievt < nEntries; ++ievt)
      {
         M->GetEntry(ievt);

         if (ievt % deltaI == 0) {
            Bar.Update(ievt);
            Bar.Print();
         }

         // Safety: cap NReco to messenger array size
         if (M->NReco > STRANGE_MAX_RECO) {
            cerr << "Warning: NReco = " << M->NReco
                 << " > STRANGE_MAX_RECO = " << STRANGE_MAX_RECO
                 << " at entry " << ievt
                 << ". Clipping to STRANGE_MAX_RECO." << endl;
         }
         int nreco = (M->NReco < STRANGE_MAX_RECO)
                     ? static_cast<int>(M->NReco)
                     : STRANGE_MAX_RECO;

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
         int nK = 0;
         int nPi = 0;

         for (int i = 0; i < nreco; ++i) {
            bool isKaonTag   = (M->RecoPIDKaon[i]   >= 2);
            bool isPionTag   = (M->RecoPIDPion[i]   >= 2);
            bool isProtonTag = (M->RecoPIDProton[i] >= 2);

            // NchTag = Sum(RecoPIDKaon>=2 || RecoPIDPion>=2 || RecoPIDProton>=2)
            if (isKaonTag || isPionTag || isProtonTag)
               ++NchTag;

            if (isKaonTag)
               ++nK;
            if (isPionTag)
               ++nPi;
         }

         // Put overflow NchTag into the last bin
         if (NchTag > par.MaxNchTag)
            NchTag = par.MaxNchTag;

         // Fill event-wise yields
         hK->Fill(NchTag,  nK);
         hPi->Fill(NchTag, nPi);
      }

      cout << endl << "Event loop finished." << endl;

      //-------------------------
      // Build K/pi ratio histogram
      //-------------------------

      hKoverPi = (TH1D*)hK->Clone("hKoverPi");
      hKoverPi->SetTitle("K/#pi yield ratio vs N_{ch}^{tag};N_{ch}^{tag};K/#pi");
      hKoverPi->Divide(hPi);
   }

   void writeHistograms()
   {
      if (outf == nullptr)
         return;

      outf->cd();

      smartWrite(hK);
      smartWrite(hPi);
      smartWrite(hKoverPi);

      TCanvas c1("c1", "K/pi vs NchTag", 800, 600);
      hKoverPi->SetMarkerStyle(20);
      hKoverPi->SetMarkerSize(1.0);
      hKoverPi->Draw("E1");
      c1.Write();  // store in ROOT file
   }
};

//============================================================//
// Main analysis
//============================================================//

int main(int argc, char *argv[])
{
   // Use the same help-message pattern as CorrelationAnalysis.cpp
//   if (printHelpMessage(argc, argv))
//      return 0;

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

   par.MinTheta = MinThetaDeg * TMath::Pi() / 180.0;
   par.MaxTheta = MaxThetaDeg * TMath::Pi() / 180.0;

   cout << "Running KtoPiAnalysis with parameters:" << endl;
   cout << "  Input       = " << par.input << endl;
   cout << "  Output      = " << par.output << endl;
   cout << "  MaxNchTag   = " << par.MaxNchTag << endl;
   cout << "  MaxEvents   = " << par.MaxEvents << endl;
   cout << "  EcmRef      = " << par.EcmRef << endl;
   cout << "  MinNch      = " << par.MinNch << endl;
   cout << "  MinThetaDeg = " << MinThetaDeg << endl;
   cout << "  MaxThetaDeg = " << MaxThetaDeg << endl;

   KtoPiAnalyzer analyzer(par);
   analyzer.analyze();
   analyzer.writeHistograms();

   cout << "Done. Output written to " << par.output << endl;

   return 0;
}
