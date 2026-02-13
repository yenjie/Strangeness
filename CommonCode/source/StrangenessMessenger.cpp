#include "StrangenessMessenger.h"
#include <iostream>

StrangenessTreeMessenger::StrangenessTreeMessenger()
   : Tree(nullptr)
{
}

StrangenessTreeMessenger::StrangenessTreeMessenger(TFile &file, const std::string &treeName)
   : Tree(nullptr)
{
   TTree *t = nullptr;
   file.GetObject(treeName.c_str(), t);
   Initialize(t);
}

StrangenessTreeMessenger::StrangenessTreeMessenger(TFile *file, const std::string &treeName)
   : Tree(nullptr)
{
   if(file == nullptr)
      return;

   TTree *t = nullptr;
   file->GetObject(treeName.c_str(), t);
   Initialize(t);
}

StrangenessTreeMessenger::StrangenessTreeMessenger(TTree *tree)
   : Tree(nullptr)
{
   Initialize(tree);
}

bool StrangenessTreeMessenger::Initialize(TTree *tree)
{
   if(tree == nullptr)
      return false;

   Tree = tree;

   // Event-level
   Tree->SetBranchAddress("Ecm",          &Ecm);
   Tree->SetBranchAddress("Nch",          &Nch);
   Tree->SetBranchAddress("Run",          &Run);
   Tree->SetBranchAddress("Event",        &Event);
   Tree->SetBranchAddress("Fill",         &Fill);
   Tree->SetBranchAddress("GoodNch",      &GoodNch);
   Tree->SetBranchAddress("GoodNneu",     &GoodNneu);
   Tree->SetBranchAddress("TotalEch",     &TotalEch);
   Tree->SetBranchAddress("TotalEneu",    &TotalEneu);
   Tree->SetBranchAddress("PassNch",      &PassNch);
   Tree->SetBranchAddress("PassThrust",   &PassThrust);
   Tree->SetBranchAddress("PassTotalE",   &PassTotalE);
   Tree->SetBranchAddress("PassAll",      &PassAll);
   Tree->SetBranchAddress("Thrust",       &Thrust);
   Tree->SetBranchAddress("ThrustX",      &ThrustX);
   Tree->SetBranchAddress("ThrustY",      &ThrustY);
   Tree->SetBranchAddress("ThrustZ",      &ThrustZ);
   Tree->SetBranchAddress("ThrustTheta",  &ThrustTheta);

   // Generator-level
   Tree->SetBranchAddress("NGen",         &NGen);
   Tree->SetBranchAddress("GenPx",        GenPx);
   Tree->SetBranchAddress("GenPy",        GenPy);
   Tree->SetBranchAddress("GenPz",        GenPz);
   Tree->SetBranchAddress("GenE",         GenE);
   Tree->SetBranchAddress("GenM",         GenM);
   Tree->SetBranchAddress("GenID",        GenID);
   Tree->SetBranchAddress("GenStatus",    GenStatus);
   Tree->SetBranchAddress("GenParent",    GenParent);
   Tree->SetBranchAddress("GenMatchIndex",GenMatchIndex);
   Tree->SetBranchAddress("GenMatchAngle",GenMatchAngle);

   // Reco-level
   Tree->SetBranchAddress("NReco",                 &NReco);
   Tree->SetBranchAddress("RecoPx",                RecoPx);
   Tree->SetBranchAddress("RecoPy",                RecoPy);
   Tree->SetBranchAddress("RecoPz",                RecoPz);
   Tree->SetBranchAddress("RecoE",                 RecoE);
   Tree->SetBranchAddress("RecoCharge",            RecoCharge);
   Tree->SetBranchAddress("RecoID",                RecoID);
   Tree->SetBranchAddress("RecoTrackLength",       RecoTrackLength);
   Tree->SetBranchAddress("RecoTrackD0",           RecoTrackD0);
   Tree->SetBranchAddress("RecoTrackZ0",           RecoTrackZ0);
   Tree->SetBranchAddress("RecoPIDElectron",       RecoPIDElectron);
   Tree->SetBranchAddress("RecoPIDProton",         RecoPIDProton);
   Tree->SetBranchAddress("RecoPIDKaon",           RecoPIDKaon);
   Tree->SetBranchAddress("RecoPIDPion",           RecoPIDPion);
   Tree->SetBranchAddress("RecoPIDHeavy",          RecoPIDHeavy);
   Tree->SetBranchAddress("RecoPIDQProton",        RecoPIDQProton);
   Tree->SetBranchAddress("RecoPIDQKaon",          RecoPIDQKaon);
   Tree->SetBranchAddress("RecoMuID",              RecoMuID);
   Tree->SetBranchAddress("RecoEleID",             RecoEleID);
   Tree->SetBranchAddress("RecoConversionID",      RecoConversionID);
   Tree->SetBranchAddress("RecoGoodTrack",         RecoGoodTrack);
   Tree->SetBranchAddress("RecoGoodNeutral",       RecoGoodNeutral);
   Tree->SetBranchAddress("RecoEfficiencyKAsK",    RecoEfficiencyKAsK);
   Tree->SetBranchAddress("RecoEfficiencyKAsPi",   RecoEfficiencyKAsPi);
   Tree->SetBranchAddress("RecoEfficiencyKAsP",    RecoEfficiencyKAsP);
   Tree->SetBranchAddress("RecoEfficiencyPiAsK",   RecoEfficiencyPiAsK);
   Tree->SetBranchAddress("RecoEfficiencyPiAsPi",  RecoEfficiencyPiAsPi);
   Tree->SetBranchAddress("RecoEfficiencyPiAsP",   RecoEfficiencyPiAsP);
   Tree->SetBranchAddress("RecoEfficiencyPAsK",    RecoEfficiencyPAsK);
   Tree->SetBranchAddress("RecoEfficiencyPAsPi",   RecoEfficiencyPAsPi);
   Tree->SetBranchAddress("RecoEfficiencyPAsP",    RecoEfficiencyPAsP);

   // Sim-level
   Tree->SetBranchAddress("NSim",         &NSim);
   Tree->SetBranchAddress("SimPx",        SimPx);
   Tree->SetBranchAddress("SimPy",        SimPy);
   Tree->SetBranchAddress("SimPz",        SimPz);
   Tree->SetBranchAddress("SimE",         SimE);
   Tree->SetBranchAddress("SimID",        SimID);

   // KShort candidates
   Tree->SetBranchAddress("NKShort",          &NKShort);
   Tree->SetBranchAddress("KShortPx",         KShortPx);
   Tree->SetBranchAddress("KShortPy",         KShortPy);
   Tree->SetBranchAddress("KShortPz",         KShortPz);
   Tree->SetBranchAddress("KShortE",          KShortE);
   Tree->SetBranchAddress("KShortSim1ID",     KShortSim1ID);
   Tree->SetBranchAddress("KShortSim2ID",     KShortSim2ID);
   Tree->SetBranchAddress("KShortReco1ID",    KShortReco1ID);
   Tree->SetBranchAddress("KShortReco2ID",    KShortReco2ID);
   Tree->SetBranchAddress("KShortReco1Angle", KShortReco1Angle);
   Tree->SetBranchAddress("KShortReco2Angle", KShortReco2Angle);
   Tree->SetBranchAddress("KShortRecoPx",     KShortRecoPx);
   Tree->SetBranchAddress("KShortRecoPy",     KShortRecoPy);
   Tree->SetBranchAddress("KShortRecoPz",     KShortRecoPz);
   Tree->SetBranchAddress("KShortRecoE",      KShortRecoE);

   // Phi candidates
   Tree->SetBranchAddress("NPhi",          &NPhi);
   Tree->SetBranchAddress("PhiPx",         PhiPx);
   Tree->SetBranchAddress("PhiPy",         PhiPy);
   Tree->SetBranchAddress("PhiPz",         PhiPz);
   Tree->SetBranchAddress("PhiE",          PhiE);
   Tree->SetBranchAddress("PhiGen1ID",     PhiGen1ID);
   Tree->SetBranchAddress("PhiGen2ID",     PhiGen2ID);
   Tree->SetBranchAddress("PhiReco1ID",    PhiReco1ID);
   Tree->SetBranchAddress("PhiReco2ID",    PhiReco2ID);
   Tree->SetBranchAddress("PhiReco1Angle", PhiReco1Angle);
   Tree->SetBranchAddress("PhiReco2Angle", PhiReco2Angle);
   Tree->SetBranchAddress("PhiRecoPx",     PhiRecoPx);
   Tree->SetBranchAddress("PhiRecoPy",     PhiRecoPy);
   Tree->SetBranchAddress("PhiRecoPz",     PhiRecoPz);
   Tree->SetBranchAddress("PhiRecoE",      PhiRecoE);

   return true;
}

bool StrangenessTreeMessenger::Initialize()
{
   if(Tree == nullptr)
      return false;
   return Initialize(Tree);
}

bool StrangenessTreeMessenger::GetEntry(long long iEntry)
{
   if(Tree == nullptr)
      return false;
   if(iEntry < 0)
      return false;
   if(iEntry >= Tree->GetEntries())
      return false;

   return (Tree->GetEntry(iEntry) > 0);
}

long long StrangenessTreeMessenger::GetEntries() const
{
   if(Tree == nullptr)
      return 0;
   return Tree->GetEntries();
}
