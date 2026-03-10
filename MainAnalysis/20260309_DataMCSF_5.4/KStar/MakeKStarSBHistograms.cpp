#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "TFile.h"
#include "TH1D.h"
#include "TNamed.h"
#include "TParameter.h"
#include "TTree.h"

namespace {
constexpr double kKaonMass = 0.493677;
constexpr double kPionMass = 0.13957039;
constexpr double kMassWindowMin = 0.70;
constexpr double kMassWindowMax = 1.10;
constexpr int kMassBins = 320;
constexpr double kAbsCosMin = 0.15;
constexpr double kAbsCosMax = 0.675;
constexpr long long kKaonTagThreshold = 2;
constexpr long long kPionTagThreshold = 2;
constexpr int kMaxReco = 10000;

struct TrackKinematics {
  double px = 0.0;
  double py = 0.0;
  double pz = 0.0;
  double charge = 0.0;
  long long kaonTag = 0;
  long long pionTag = 0;
};

double buildMass(const TrackKinematics& kaon, const TrackKinematics& pion) {
  const double pksq = kaon.px * kaon.px + kaon.py * kaon.py + kaon.pz * kaon.pz;
  const double ppisq = pion.px * pion.px + pion.py * pion.py + pion.pz * pion.pz;
  const double eK = std::sqrt(pksq + kKaonMass * kKaonMass);
  const double ePi = std::sqrt(ppisq + kPionMass * kPionMass);

  const double px = kaon.px + pion.px;
  const double py = kaon.py + pion.py;
  const double pz = kaon.pz + pion.pz;
  const double e = eK + ePi;
  const double m2 = e * e - (px * px + py * py + pz * pz);
  return (m2 > 0.0) ? std::sqrt(m2) : 0.0;
}

bool passAcceptance(const TrackKinematics& t) {
  const double p = std::sqrt(t.px * t.px + t.py * t.py + t.pz * t.pz);
  if (p <= 0.0) return false;
  const double absCosTheta = std::fabs(t.pz / p);
  return (absCosTheta >= kAbsCosMin && absCosTheta <= kAbsCosMax);
}

std::string getArgument(int argc, char* argv[], const std::string& option,
                        const std::string& defaultValue) {
  for (int i = 1; i + 1 < argc; ++i)
    if (argv[i] == option) return argv[i + 1];
  return defaultValue;
}

double getDoubleArgument(int argc, char* argv[], const std::string& option, double defaultValue) {
  const std::string value = getArgument(argc, argv, option, "");
  return value.empty() ? defaultValue : std::stod(value);
}
}  // namespace

int main(int argc, char* argv[]) {
  const std::string inputFileName =
      getArgument(argc, argv, "--input", "../../../../Samples/merged_mc_v2.3.root");
  const std::string outputFileName =
      getArgument(argc, argv, "--output", "KStarSBHistograms.root");
  const std::string treeName = getArgument(argc, argv, "--tree", "Tree");
  const double massMin = getDoubleArgument(argc, argv, "--mass-min", kMassWindowMin);
  const double massMax = getDoubleArgument(argc, argv, "--mass-max", kMassWindowMax);

  TFile inputFile(inputFileName.c_str(), "READ");
  if (inputFile.IsZombie()) {
    std::cerr << "Error: cannot open input file " << inputFileName << std::endl;
    return 1;
  }

  TTree* tree = nullptr;
  inputFile.GetObject(treeName.c_str(), tree);
  if (tree == nullptr) {
    std::cerr << "Error: cannot find tree '" << treeName << "' in " << inputFileName << std::endl;
    return 1;
  }

  long long nReco = 0;
  double recoPx[kMaxReco] = {0.0};
  double recoPy[kMaxReco] = {0.0};
  double recoPz[kMaxReco] = {0.0};
  double recoCharge[kMaxReco] = {0.0};
  long long recoPIDKaon[kMaxReco] = {0};
  long long recoPIDPion[kMaxReco] = {0};
  long long recoGoodTrack[kMaxReco] = {0};

  tree->SetBranchAddress("NReco", &nReco);
  tree->SetBranchAddress("RecoPx", recoPx);
  tree->SetBranchAddress("RecoPy", recoPy);
  tree->SetBranchAddress("RecoPz", recoPz);
  tree->SetBranchAddress("RecoCharge", recoCharge);
  tree->SetBranchAddress("RecoPIDKaon", recoPIDKaon);
  tree->SetBranchAddress("RecoPIDPion", recoPIDPion);
  tree->SetBranchAddress("RecoGoodTrack", recoGoodTrack);

  TH1D hMassKaonTag(
      "hKStarSBMassKaonTag",
      "K^{*} same-event reco pairs, kaon-tag; m(K#pi) [GeV]; Assignments / bin",
      kMassBins, massMin, massMax);
  TH1D hMassKaonPionTag(
      "hKStarSBMassKaonPionTag",
      "K^{*} same-event reco pairs, kaon+pion-tag; m(K#pi) [GeV]; Assignments / bin",
      kMassBins, massMin, massMax);
  TH1D hMassDoubleKaonTag(
      "hKStarSBMassDoubleKaonTag",
      "K^{*} same-event reco pairs, two-kaon-tag; m(K#pi) [GeV]; Assignments / bin",
      kMassBins, massMin, massMax);
  TH1D hMassAccepted(
      "hKStarSBMassAccepted",
      "K^{*} same-event reco OS pairs, accepted; m(K#pi) [GeV]; Assignments / bin",
      kMassBins, massMin, massMax);

  long long acceptedTracks = 0;
  long long oppositeSignPairs = 0;
  long long positiveKaonAssignments = 0;
  long long negativeKaonAssignments = 0;
  long long countKaonTag = 0;
  long long countKaonPionTag = 0;
  long long countDoubleKaonTag = 0;

  const long long entryCount = tree->GetEntries();
  for (long long entry = 0; entry < entryCount; ++entry) {
    tree->GetEntry(entry);

    std::vector<TrackKinematics> tracks;
    tracks.reserve(nReco);
    for (long long i = 0; i < nReco; ++i) {
      if (recoGoodTrack[i] != 1) continue;
      if (recoCharge[i] == 0.0) continue;
      TrackKinematics t{recoPx[i], recoPy[i], recoPz[i], recoCharge[i], recoPIDKaon[i], recoPIDPion[i]};
      if (!passAcceptance(t)) continue;
      tracks.push_back(t);
      acceptedTracks++;
    }

    const int trackCount = static_cast<int>(tracks.size());
    for (int i = 0; i < trackCount; ++i) {
      for (int j = i + 1; j < trackCount; ++j) {
        if (tracks[i].charge * tracks[j].charge >= 0.0) continue;
        oppositeSignPairs++;

        const TrackKinematics* positive = (tracks[i].charge > 0.0) ? &tracks[i] : &tracks[j];
        const TrackKinematics* negative = (tracks[i].charge > 0.0) ? &tracks[j] : &tracks[i];

        const double positiveKaonMass = buildMass(*positive, *negative);
        hMassAccepted.Fill(positiveKaonMass);
        positiveKaonAssignments++;
        if (positive->kaonTag >= kKaonTagThreshold) {
          hMassKaonTag.Fill(positiveKaonMass);
          countKaonTag++;
          if (negative->kaonTag >= kKaonTagThreshold) {
            hMassDoubleKaonTag.Fill(positiveKaonMass);
            countDoubleKaonTag++;
          }
          if (negative->pionTag >= kPionTagThreshold) {
            hMassKaonPionTag.Fill(positiveKaonMass);
            countKaonPionTag++;
          }
        }

        const double negativeKaonMass = buildMass(*negative, *positive);
        hMassAccepted.Fill(negativeKaonMass);
        negativeKaonAssignments++;
        if (negative->kaonTag >= kKaonTagThreshold) {
          hMassKaonTag.Fill(negativeKaonMass);
          countKaonTag++;
          if (positive->kaonTag >= kKaonTagThreshold) {
            hMassDoubleKaonTag.Fill(negativeKaonMass);
            countDoubleKaonTag++;
          }
          if (positive->pionTag >= kPionTagThreshold) {
            hMassKaonPionTag.Fill(negativeKaonMass);
            countKaonPionTag++;
          }
        }
      }
    }
  }

  TFile outputFile(outputFileName.c_str(), "RECREATE");
  hMassKaonTag.Write();
  hMassKaonPionTag.Write();
  hMassDoubleKaonTag.Write();
  hMassAccepted.Write();

  TNamed selection(
      "SelectionSummary",
      Form("Reco-only KStar same-event OS pairs: RecoGoodTrack==1, nonzero charge, "
           "0.15<=|cos(theta)|<=0.675 on both tracks. Each OS pair is filled twice, once with the "
           "positive track treated as kaon and once with the negative track treated as kaon. "
           "The assumed kaon track must satisfy RecoPIDKaon>=2 for the kaon-tag histogram, and the "
           "assumed pion track must also satisfy RecoPIDPion>=2 for the kaon+pion-tag histogram. "
           "Histogram range %.3f-%.3f GeV.",
           massMin, massMax));
  selection.Write();
  TParameter<long long>("AcceptedTracks", acceptedTracks).Write();
  TParameter<long long>("OppositeSignPairs", oppositeSignPairs).Write();
  TParameter<long long>("PositiveKaonAssignments", positiveKaonAssignments).Write();
  TParameter<long long>("NegativeKaonAssignments", negativeKaonAssignments).Write();
  TParameter<long long>("CountKaonTag", countKaonTag).Write();
  TParameter<long long>("CountKaonPionTag", countKaonPionTag).Write();
  TParameter<long long>("CountDoubleKaonTag", countDoubleKaonTag).Write();
  TParameter<double>("MassMin", massMin).Write();
  TParameter<double>("MassMax", massMax).Write();
  outputFile.Close();

  std::cout << "Wrote " << outputFileName << std::endl;
  std::cout << "  Accepted tracks:         " << acceptedTracks << std::endl;
  std::cout << "  Opposite-sign pairs:     " << oppositeSignPairs << std::endl;
  std::cout << "  Positive-kaon fills:     " << positiveKaonAssignments << std::endl;
  std::cout << "  Negative-kaon fills:     " << negativeKaonAssignments << std::endl;
  std::cout << "  Kaon-tag assignments:    " << countKaonTag << std::endl;
  std::cout << "  Kaon+pion assignments:   " << countKaonPionTag << std::endl;
  std::cout << "  Two-kaon assignments:    " << countDoubleKaonTag << std::endl;
  return 0;
}
