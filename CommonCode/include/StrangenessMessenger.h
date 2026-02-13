#ifndef STRANGENESS_MESSENGER_H
#define STRANGENESS_MESSENGER_H

#include <string>
#include "TTree.h"
#include "TFile.h"

// These are generous upper bounds; adjust if you learn better maxima from the file
#define STRANGE_MAX_GEN     10000
#define STRANGE_MAX_RECO    10000
#define STRANGE_MAX_SIM     10000
#define STRANGE_MAX_KSHORT  4096
#define STRANGE_MAX_PHI     4096

class StrangenessTreeMessenger
{
public:
   TTree *Tree;

   // Event-level scalars
   double     Ecm;
   long long  Nch;
   long long  Run;
   long long  Event;
   long long  Fill;
   long long  GoodNch;
   long long  GoodNneu;
   double     TotalEch;
   double     TotalEneu;
   long long  PassNch;
   long long  PassThrust;
   long long  PassTotalE;
   long long  PassAll;
   double     Thrust;
   double     ThrustX;
   double     ThrustY;
   double     ThrustZ;
   double     ThrustTheta;

   // Generator-level (truth) particles
   long long  NGen;
   double     GenPx[STRANGE_MAX_GEN];
   double     GenPy[STRANGE_MAX_GEN];
   double     GenPz[STRANGE_MAX_GEN];
   double     GenE[STRANGE_MAX_GEN];
   double     GenM[STRANGE_MAX_GEN];
   long long  GenID[STRANGE_MAX_GEN];
   long long  GenStatus[STRANGE_MAX_GEN];
   long long  GenParent[STRANGE_MAX_GEN];
   long long  GenMatchIndex[STRANGE_MAX_GEN];
   double     GenMatchAngle[STRANGE_MAX_GEN];

   // Reconstructed particles
   long long  NReco;
   double     RecoPx[STRANGE_MAX_RECO];
   double     RecoPy[STRANGE_MAX_RECO];
   double     RecoPz[STRANGE_MAX_RECO];
   double     RecoE[STRANGE_MAX_RECO];
   double     RecoCharge[STRANGE_MAX_RECO];
   long long  RecoID[STRANGE_MAX_RECO];
   double     RecoTrackLength[STRANGE_MAX_RECO];
   double     RecoTrackD0[STRANGE_MAX_RECO];
   double     RecoTrackZ0[STRANGE_MAX_RECO];
   long long  RecoPIDElectron[STRANGE_MAX_RECO];
   long long  RecoPIDProton[STRANGE_MAX_RECO];
   long long  RecoPIDKaon[STRANGE_MAX_RECO];
   long long  RecoPIDPion[STRANGE_MAX_RECO];
   long long  RecoPIDHeavy[STRANGE_MAX_RECO];
   double     RecoPIDQProton[STRANGE_MAX_RECO];
   double     RecoPIDQKaon[STRANGE_MAX_RECO];
   long long  RecoMuID[STRANGE_MAX_RECO];
   long long  RecoEleID[STRANGE_MAX_RECO];
   long long  RecoConversionID[STRANGE_MAX_RECO];
   long long  RecoGoodTrack[STRANGE_MAX_RECO];
   long long  RecoGoodNeutral[STRANGE_MAX_RECO];
   double     RecoEfficiencyKAsK[STRANGE_MAX_RECO];
   double     RecoEfficiencyKAsPi[STRANGE_MAX_RECO];
   double     RecoEfficiencyKAsP[STRANGE_MAX_RECO];
   double     RecoEfficiencyPiAsK[STRANGE_MAX_RECO];
   double     RecoEfficiencyPiAsPi[STRANGE_MAX_RECO];
   double     RecoEfficiencyPiAsP[STRANGE_MAX_RECO];
   double     RecoEfficiencyPAsK[STRANGE_MAX_RECO];
   double     RecoEfficiencyPAsPi[STRANGE_MAX_RECO];
   double     RecoEfficiencyPAsP[STRANGE_MAX_RECO];

   // Simulation-level particles
   long long  NSim;
   double     SimPx[STRANGE_MAX_SIM];
   double     SimPy[STRANGE_MAX_SIM];
   double     SimPz[STRANGE_MAX_SIM];
   double     SimE[STRANGE_MAX_SIM];
   long long  SimID[STRANGE_MAX_SIM];

   // K0S candidates
   long long  NKShort;
   double     KShortPx[STRANGE_MAX_KSHORT];
   double     KShortPy[STRANGE_MAX_KSHORT];
   double     KShortPz[STRANGE_MAX_KSHORT];
   double     KShortE[STRANGE_MAX_KSHORT];
   long long  KShortSim1ID[STRANGE_MAX_KSHORT];
   long long  KShortSim2ID[STRANGE_MAX_KSHORT];
   long long  KShortReco1ID[STRANGE_MAX_KSHORT];
   long long  KShortReco2ID[STRANGE_MAX_KSHORT];
   double     KShortReco1Angle[STRANGE_MAX_KSHORT];
   double     KShortReco2Angle[STRANGE_MAX_KSHORT];
   double     KShortRecoPx[STRANGE_MAX_KSHORT];
   double     KShortRecoPy[STRANGE_MAX_KSHORT];
   double     KShortRecoPz[STRANGE_MAX_KSHORT];
   double     KShortRecoE[STRANGE_MAX_KSHORT];

   // Phi meson candidates
   long long  NPhi;
   double     PhiPx[STRANGE_MAX_PHI];
   double     PhiPy[STRANGE_MAX_PHI];
   double     PhiPz[STRANGE_MAX_PHI];
   double     PhiE[STRANGE_MAX_PHI];
   long long  PhiGen1ID[STRANGE_MAX_PHI];
   long long  PhiGen2ID[STRANGE_MAX_PHI];
   long long  PhiReco1ID[STRANGE_MAX_PHI];
   long long  PhiReco2ID[STRANGE_MAX_PHI];
   double     PhiReco1Angle[STRANGE_MAX_PHI];
   double     PhiReco2Angle[STRANGE_MAX_PHI];
   double     PhiRecoPx[STRANGE_MAX_PHI];
   double     PhiRecoPy[STRANGE_MAX_PHI];
   double     PhiRecoPz[STRANGE_MAX_PHI];
   double     PhiRecoE[STRANGE_MAX_PHI];

public:
   StrangenessTreeMessenger();
   StrangenessTreeMessenger(TFile &file, const std::string &treeName = "Tree");
   StrangenessTreeMessenger(TFile *file, const std::string &treeName = "Tree");
   StrangenessTreeMessenger(TTree *tree);

   bool Initialize(TTree *tree);   // attach to given tree and set branch addresses
   bool Initialize();              // reuse existing Tree pointer

   bool       GetEntry(long long iEntry);
   long long  GetEntries() const;
};

#endif
