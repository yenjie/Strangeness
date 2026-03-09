#include <TFile.h>
#include <TROOT.h>
#include <TTree.h>

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace {
constexpr double kKaonMass = 0.493677;
constexpr double kAbsCosMin = 0.15;
constexpr double kAbsCosMax = 0.675;
constexpr int kMaxReco = 256;

double buildMass(double px1, double py1, double pz1,
                 double px2, double py2, double pz2) {
  const double p1sq = px1 * px1 + py1 * py1 + pz1 * pz1;
  const double p2sq = px2 * px2 + py2 * py2 + pz2 * pz2;
  const double e1 = std::sqrt(p1sq + kKaonMass * kKaonMass);
  const double e2 = std::sqrt(p2sq + kKaonMass * kKaonMass);
  const double e = e1 + e2;
  const double px = px1 + px2;
  const double py = py1 + py2;
  const double pz = pz1 + pz2;
  const double m2 = e * e - (px * px + py * py + pz * pz);
  return (m2 > 0.0 ? std::sqrt(m2) : 0.0);
}

bool passTrackAcceptance(double px, double py, double pz) {
  const double p = std::sqrt(px * px + py * py + pz * pz);
  if (p <= 0.0) return false;
  const double abscos = std::fabs(pz / p);
  return (abscos >= kAbsCosMin && abscos <= kAbsCosMax);
}

double trackP(double px, double py, double pz) {
  return std::sqrt(px * px + py * py + pz * pz);
}

int findBin(double p, const std::vector<double>& edges) {
  if (edges.size() < 2) return -1;
  if (p < edges.front() || p >= edges.back()) return -1;
  for (size_t i = 0; i + 1 < edges.size(); ++i) {
    if (p >= edges[i] && p < edges[i + 1]) return static_cast<int>(i);
  }
  return -1;
}
}  // namespace

void phi_momentum_bin_stats(const char* input = "Samples/merged_mc_v2.2.root",
                            const char* outDir = "Reports",
                            const char* outPrefix = "phi_pbin_mc",
                            double massMin = 0.99,
                            double massMax = 1.06) {
  gROOT->SetBatch(kTRUE);

  // Candidate working binning for kaon momentum.
  const std::vector<double> pEdges = {0.0, 0.6, 0.9, 1.2, 1.6, 2.2, 5.0};
  const int nBins = static_cast<int>(pEdges.size()) - 1;

  TFile f(input, "READ");
  if (f.IsZombie()) {
    std::cerr << "Failed to open input file: " << input << std::endl;
    return;
  }
  TTree* tree = nullptr;
  f.GetObject("Tree", tree);
  if (!tree) {
    std::cerr << "Cannot find TTree 'Tree' in " << input << std::endl;
    return;
  }

  Long64_t nReco = 0;
  double recoPx[kMaxReco], recoPy[kMaxReco], recoPz[kMaxReco], recoCharge[kMaxReco];
  Long64_t recoGoodTrack[kMaxReco], recoPIDKaon[kMaxReco];

  tree->SetBranchAddress("NReco", &nReco);
  tree->SetBranchAddress("RecoPx", recoPx);
  tree->SetBranchAddress("RecoPy", recoPy);
  tree->SetBranchAddress("RecoPz", recoPz);
  tree->SetBranchAddress("RecoCharge", recoCharge);
  tree->SetBranchAddress("RecoGoodTrack", recoGoodTrack);
  tree->SetBranchAddress("RecoPIDKaon", recoPIDKaon);

  std::vector<Long64_t> pairsBothInBin(nBins, 0);
  std::vector<Long64_t> pairsAnyInBin(nBins, 0);
  std::vector<Long64_t> tracksInBin(nBins, 0);
  std::vector<Long64_t> n0tagBoth(nBins, 0);
  std::vector<Long64_t> n1tagBoth(nBins, 0);
  std::vector<Long64_t> n2tagBoth(nBins, 0);

  Long64_t pairTotal = 0;
  Long64_t pairInMass = 0;

  const Long64_t nEntries = tree->GetEntries();
  for (Long64_t ie = 0; ie < nEntries; ++ie) {
    tree->GetEntry(ie);
    if (nReco > kMaxReco) continue;

    std::vector<int> pos, neg;
    pos.reserve(32);
    neg.reserve(32);

    for (int i = 0; i < nReco; ++i) {
      if (recoGoodTrack[i] == 0) continue;
      if (!passTrackAcceptance(recoPx[i], recoPy[i], recoPz[i])) continue;
      if (recoCharge[i] > 0) pos.push_back(i);
      if (recoCharge[i] < 0) neg.push_back(i);
    }

    for (int ip : pos) {
      for (int in : neg) {
        ++pairTotal;
        const double m = buildMass(recoPx[ip], recoPy[ip], recoPz[ip], recoPx[in], recoPy[in], recoPz[in]);
        if (m < massMin || m > massMax) continue;
        ++pairInMass;

        const double p1 = trackP(recoPx[ip], recoPy[ip], recoPz[ip]);
        const double p2 = trackP(recoPx[in], recoPy[in], recoPz[in]);
        const int b1 = findBin(p1, pEdges);
        const int b2 = findBin(p2, pEdges);

        if (b1 >= 0) ++tracksInBin[b1];
        if (b2 >= 0) ++tracksInBin[b2];
        if (b1 >= 0) ++pairsAnyInBin[b1];
        if (b2 >= 0 && b2 != b1) ++pairsAnyInBin[b2];

        if (b1 >= 0 && b1 == b2) {
          ++pairsBothInBin[b1];
          const int nTag = (recoPIDKaon[ip] >= 2 ? 1 : 0) + (recoPIDKaon[in] >= 2 ? 1 : 0);
          if (nTag == 0) ++n0tagBoth[b1];
          if (nTag == 1) ++n1tagBoth[b1];
          if (nTag == 2) ++n2tagBoth[b1];
        }
      }
    }
  }

  const std::string outCsv = std::string(outDir) + "/" + outPrefix + "_stats.csv";
  std::ofstream csv(outCsv);
  csv << "p_min,p_max,pairs_both_in_bin,pairs_with_any_daughter_in_bin,tracks_in_bin,bothbin_0tag,bothbin_1tag,bothbin_2tag\n";
  for (int i = 0; i < nBins; ++i) {
    csv << std::fixed << std::setprecision(3) << pEdges[i] << "," << pEdges[i + 1] << ","
        << pairsBothInBin[i] << "," << pairsAnyInBin[i] << "," << tracksInBin[i] << ","
        << n0tagBoth[i] << "," << n1tagBoth[i] << "," << n2tagBoth[i] << "\n";
  }
  csv.close();

  const std::string outTxt = std::string(outDir) + "/" + outPrefix + "_summary.txt";
  std::ofstream txt(outTxt);
  txt << "Input," << input << "\n";
  txt << "Selection,RecoGoodTrack&&0.15<=|cos(theta)|<=0.675&&opposite_charge\n";
  txt << "TagWP,RecoPIDKaon>=2\n";
  txt << "MassRange," << std::fixed << std::setprecision(3) << massMin << "," << massMax << "\n";
  txt << "PairTotal," << pairTotal << "\n";
  txt << "PairInMassRange," << pairInMass << "\n";
  txt << "OutputCSV," << outCsv << "\n\n";
  txt << "bin,p_min,p_max,pairs_both_in_bin,pairs_with_any_daughter_in_bin,tracks_in_bin,bothbin_0tag,bothbin_1tag,bothbin_2tag\n";
  for (int i = 0; i < nBins; ++i) {
    txt << i << "," << std::fixed << std::setprecision(3) << pEdges[i] << "," << pEdges[i + 1] << ","
        << pairsBothInBin[i] << "," << pairsAnyInBin[i] << "," << tracksInBin[i] << ","
        << n0tagBoth[i] << "," << n1tagBoth[i] << "," << n2tagBoth[i] << "\n";
  }
  txt.close();

  std::cout << "Wrote " << outCsv << std::endl;
  std::cout << "Wrote " << outTxt << std::endl;
}
