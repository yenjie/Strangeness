#include "Pythia8/Analysis.h"
#include "Pythia8/Event.h"

#include "TFile.h"
#include "TH1D.h"
#include "TNamed.h"
#include "TParameter.h"
#include "TProfile.h"
#include "TTree.h"

#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

using namespace Pythia8;

namespace {

double computeAxisRapidity(double px, double py, double pz, double e,
                           double ax, double ay, double az, bool &ok) {
  const double pLong = px * ax + py * ay + pz * az;
  const double plus = e + pLong;
  const double minus = e - pLong;
  if (plus <= 0.0 || minus <= 0.0) {
    ok = false;
    return 0.0;
  }
  const double y = 0.5 * std::log(plus / minus);
  ok = std::isfinite(y);
  return y;
}

}

int main(int argc, char *argv[]) {
  if (argc < 4) {
    std::cerr << "Usage: " << argv[0]
              << " <input_truth.root> <template_dndy.root> <output_dndy.root>\n";
    return 1;
  }

  const std::string inputPath = argv[1];
  const std::string templatePath = argv[2];
  const std::string outputPath = argv[3];

  std::unique_ptr<TFile> inputFile(TFile::Open(inputPath.c_str(), "READ"));
  if (!inputFile || inputFile->IsZombie()) {
    std::cerr << "Cannot open input file: " << inputPath << "\n";
    return 1;
  }
  TTree *tree = dynamic_cast<TTree *>(inputFile->Get("Events"));
  if (!tree) {
    std::cerr << "Missing Events tree in " << inputPath << "\n";
    return 1;
  }

  std::unique_ptr<TFile> templateFile(TFile::Open(templatePath.c_str(), "READ"));
  if (!templateFile || templateFile->IsZombie()) {
    std::cerr << "Cannot open template file: " << templatePath << "\n";
    return 1;
  }
  TH1D *templateHist = dynamic_cast<TH1D *>(templateFile->Get("hRatioMcTrue_dNdY"));
  if (!templateHist) {
    templateHist = dynamic_cast<TH1D *>(templateFile->Get("hRatioDataBayes_dNdY"));
  }
  if (!templateHist) {
    std::cerr << "Missing dNdY template histogram in " << templatePath << "\n";
    return 1;
  }

  auto hKCounts = std::unique_ptr<TH1D>(static_cast<TH1D *>(templateHist->Clone("hKCountsTruth_dNdY")));
  auto hPiCounts = std::unique_ptr<TH1D>(static_cast<TH1D *>(templateHist->Clone("hPiCountsTruth_dNdY")));
  auto hRatio = std::unique_ptr<TH1D>(static_cast<TH1D *>(templateHist->Clone("hKPiTruth_dNdY")));
  auto hDNdY = std::unique_ptr<TH1D>(static_cast<TH1D *>(templateHist->Clone("hNChY05Truth")));
  std::unique_ptr<TProfile> pRatio;
  if (templateHist->GetXaxis()->GetXbins()->GetSize() > 0) {
    pRatio = std::unique_ptr<TProfile>(new TProfile("pKPiVsDNdY",
        "K/#pi vs thrust-axis dN_{ch}/dy;dN_{ch}/dy (|y_{T}|<0.5);K/#pi",
        templateHist->GetNbinsX(),
        templateHist->GetXaxis()->GetXbins()->GetArray()));
  } else {
    pRatio = std::unique_ptr<TProfile>(new TProfile("pKPiVsDNdY",
        "K/#pi vs thrust-axis dN_{ch}/dy;dN_{ch}/dy (|y_{T}|<0.5);K/#pi",
        templateHist->GetNbinsX(),
        templateHist->GetXaxis()->GetXmin(),
        templateHist->GetXaxis()->GetXmax()));
  }

  hKCounts->Reset();
  hPiCounts->Reset();
  hRatio->Reset();
  hDNdY->Reset();
  hKCounts->Sumw2();
  hPiCounts->Sumw2();
  hRatio->Sumw2();
  hDNdY->Sumw2();

  std::vector<int> *pdg = nullptr;
  std::vector<float> *charge = nullptr;
  std::vector<float> *px = nullptr;
  std::vector<float> *py = nullptr;
  std::vector<float> *pz = nullptr;
  std::vector<float> *energy = nullptr;
  std::vector<float> *mass = nullptr;
  std::vector<float> *pt = nullptr;

  tree->SetBranchAddress("pdg", &pdg);
  tree->SetBranchAddress("charge", &charge);
  tree->SetBranchAddress("px", &px);
  tree->SetBranchAddress("py", &py);
  tree->SetBranchAddress("pz", &pz);
  tree->SetBranchAddress("e", &energy);
  tree->SetBranchAddress("m", &mass);
  tree->SetBranchAddress("pt", &pt);

  Event event;
  Thrust thrust(1);

  long long nEntries = tree->GetEntries();
  long long nAccepted = 0;
  long long nThrustFail = 0;

  for (long long i = 0; i < nEntries; ++i) {
    tree->GetEntry(i);
    if (!pdg || !charge || !px || !py || !pz || !energy || !mass || !pt) {
      continue;
    }

    event.reset();
    const int n = static_cast<int>(pdg->size());
    for (int j = 0; j < n; ++j) {
      event.append((*pdg)[j], 1, 0, 0, (*px)[j], (*py)[j], (*pz)[j], (*energy)[j], (*mass)[j]);
    }

    if (!thrust.analyze(event)) {
      ++nThrustFail;
      continue;
    }

    Vec4 axis = thrust.eventAxis(1);
    const double norm = std::sqrt(axis.px() * axis.px() + axis.py() * axis.py() + axis.pz() * axis.pz());
    if (norm <= 0.0) {
      ++nThrustFail;
      continue;
    }
    const double ax = axis.px() / norm;
    const double ay = axis.py() / norm;
    const double az = axis.pz() / norm;

    int nChY05 = 0;
    int nPiPt0405 = 0;
    int nKPt0405 = 0;
    for (int j = 0; j < n; ++j) {
      bool ok = false;
      const double yT = computeAxisRapidity((*px)[j], (*py)[j], (*pz)[j], (*energy)[j], ax, ay, az, ok);
      if ((*charge)[j] != 0.0f && ok && std::abs(yT) < 0.5) {
        ++nChY05;
      }
      const int apdg = std::abs((*pdg)[j]);
      const double ptv = (*pt)[j];
      if (ptv >= 0.4 && ptv < 5.0 && (*charge)[j] != 0.0f) {
        if (apdg == 211) ++nPiPt0405;
        if (apdg == 321) ++nKPt0405;
      }
    }

    hDNdY->Fill(static_cast<double>(nChY05));
    hKCounts->Fill(static_cast<double>(nChY05), nKPt0405);
    hPiCounts->Fill(static_cast<double>(nChY05), nPiPt0405);
    if (nPiPt0405 > 0) {
      pRatio->Fill(static_cast<double>(nChY05), static_cast<double>(nKPt0405) / static_cast<double>(nPiPt0405));
    }
    ++nAccepted;
  }

  hRatio->Divide(hKCounts.get(), hPiCounts.get(), 1.0, 1.0, "");

  TFile outputFile(outputPath.c_str(), "RECREATE");
  hKCounts->Write();
  hPiCounts->Write();
  hRatio->Write();
  hDNdY->Write();
  pRatio->Write();
  TNamed("sourceTruthFile", inputPath.c_str()).Write();
  TNamed("templateDNdYFile", templatePath.c_str()).Write();
  TNamed("thrustDefinition", "Pythia8::Thrust(select=1) on stored final-state particles").Write();
  TParameter<long long>("nProcessed", nEntries).Write();
  TParameter<long long>("nAccepted", nAccepted).Write();
  TParameter<long long>("nThrustFail", nThrustFail).Write();
  outputFile.Close();

  std::cout << "Wrote " << outputPath << "\n";
  std::cout << "Processed events: " << nEntries << "\n";
  std::cout << "Accepted events: " << nAccepted << "\n";
  std::cout << "Thrust failures: " << nThrustFail << "\n";
  return 0;
}
