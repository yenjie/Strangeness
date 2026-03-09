#include <TFile.h>
#include <TF1.h>
#include <TH1D.h>
#include <TMath.h>

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace {
constexpr double kFitMin = 1.00;
constexpr double kFitMax = 1.05;

struct Seed {
  std::string model = "TripleGauss";
  double mean = 1.0195;
  double s1 = 0.0025;
  double s2 = 0.0060;
  double s3 = 0.0120;
  double gamma = 0.004;
};

std::vector<std::string> splitCsv(const std::string& s) {
  std::vector<std::string> out;
  std::stringstream ss(s);
  std::string item;
  while (std::getline(ss, item, ',')) out.push_back(item);
  return out;
}

std::map<std::string, Seed> readSeeds(const std::string& csvPath) {
  std::map<std::string, Seed> out;
  std::ifstream in(csvPath);
  std::string line;
  if (!in) return out;
  std::getline(in, line);
  while (std::getline(in, line)) {
    if (line.empty()) continue;
    auto c = splitCsv(line);
    if (c.size() < 12) continue;
    Seed s;
    s.model = c[1];
    s.mean = std::atof(c[3].c_str());
    s.s1 = std::atof(c[5].c_str());
    s.s2 = std::atof(c[7].c_str());
    s.s3 = std::atof(c[9].c_str());
    s.gamma = std::atof(c[11].c_str());
    out[c[0]] = s;
  }
  return out;
}

void fitStep1Triple(TH1D* h, const std::string& label, std::ofstream& csv, std::ofstream& txt) {
  TF1 f(("f_s1_" + label).c_str(),
        "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]*exp(-0.5*((x-[1])/[4])^2)+[5]*exp(-0.5*((x-[1])/[6])^2)",
        kFitMin, kFitMax);
  double m = h->GetMaximum();
  f.SetParameters(0.5*m, 1.0196, 0.0025, 0.3*m, 0.006, 0.2*m, 0.012);
  f.SetParLimits(0, 0, 1e12);
  f.SetParLimits(1, 1.015, 1.024);
  f.SetParLimits(2, 0.0005, 0.02);
  f.SetParLimits(3, 0, 1e12);
  f.SetParLimits(4, 0.0005, 0.04);
  f.SetParLimits(5, 0, 1e12);
  f.SetParLimits(6, 0.0005, 0.06);
  h->Fit(&f, "RQ0");

  double A1=f.GetParameter(0), mean=f.GetParameter(1), s1=f.GetParameter(2);
  double A2=f.GetParameter(3), s2=f.GetParameter(4);
  double A3=f.GetParameter(5), s3=f.GetParameter(6);
  double I1=A1*s1, I2=A2*s2, I3=A3*s3;
  double It=I1+I2+I3;
  double f1 = It>0?I1/It:0, f2=It>0?I2/It:0, f3=It>0?I3/It:0;

  csv << "step1," << label << ",TripleGauss," << f.GetChisquare()/f.GetNDF()
      << "," << mean << "," << f.GetParError(1)
      << "," << A1 << "," << f.GetParError(0)
      << "," << s1 << "," << f.GetParError(2)
      << "," << A2 << "," << f.GetParError(3)
      << "," << s2 << "," << f.GetParError(4)
      << "," << A3 << "," << f.GetParError(5)
      << "," << s3 << "," << f.GetParError(6)
      << "," << f1 << "," << f2 << "," << f3
      << ",,,,,,,\n";

  txt << "Step1 " << label << " TripleGauss: mean=" << mean
      << " s1=" << s1 << " s2=" << s2 << " s3=" << s3
      << " A1=" << A1 << " A2=" << A2 << " A3=" << A3
      << " frac=" << f1 << "/" << f2 << "/" << f3 << "\n";
}

void fitStep23(TH1D* h, const Seed& sd, const std::string& step, const std::string& label, std::ofstream& csv, std::ofstream& txt) {
  double maxY = std::max(10.0, h->GetMaximum());
  double b0 = 0.5*(h->GetBinContent(1)+h->GetBinContent(h->GetNbinsX()));
  double meanLo = std::max(1.017, sd.mean-0.0015);
  double meanHi = std::min(1.022, sd.mean+0.0015);

  TF1 f(("f_"+step+"_"+label).c_str(),
        "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]*exp(-0.5*((x-[1])/[4])^2)+[5]*exp(-0.5*((x-[1])/[6])^2)+pol2(7)+gaus(10)",
        kFitMin, kFitMax);
  f.SetParameters(0.5*maxY, sd.mean, std::max(0.0008,sd.s1),
                  0.3*maxY, std::max(0.001, sd.s2>0?sd.s2:0.004),
                  0.2*maxY, std::max(0.001, sd.s3>0?sd.s3:0.008),
                  b0, 0.0, 0.0);
  f.SetParameter(10, 0.1*maxY);
  f.SetParameter(11, 1.0);
  f.SetParameter(12, 0.002);
  f.SetParLimits(0,0,1e12); f.SetParLimits(3,0,1e12); f.SetParLimits(5,0,1e12);
  f.SetParLimits(1,meanLo,meanHi);
  f.SetParLimits(2,0.0008,0.02); f.SetParLimits(4,0.0008,0.04); f.SetParLimits(6,0.001,0.06);
  f.SetParLimits(10,0,1e12); f.SetParLimits(11,0.996,1.004); f.SetParLimits(12,0.0008,0.01);
  h->Fit(&f, "RQ0");

  double A1=f.GetParameter(0), mean=f.GetParameter(1), s1=f.GetParameter(2);
  double A2=f.GetParameter(3), s2=f.GetParameter(4);
  double A3=f.GetParameter(5), s3=f.GetParameter(6);
  double I1=A1*s1, I2=A2*s2, I3=A3*s3;
  double It=I1+I2+I3;
  double f1 = It>0?I1/It:0, f2=It>0?I2/It:0, f3=It>0?I3/It:0;

  double p0=f.GetParameter(7), p1=f.GetParameter(8), p2=f.GetParameter(9);
  double Ag=f.GetParameter(10), mug=f.GetParameter(11), sigg=f.GetParameter(12);

  csv << step << "," << label << ",TripleGauss+pol2+gaus1.00," << f.GetChisquare()/f.GetNDF()
      << "," << mean << "," << f.GetParError(1)
      << "," << A1 << "," << f.GetParError(0)
      << "," << s1 << "," << f.GetParError(2)
      << "," << A2 << "," << f.GetParError(3)
      << "," << s2 << "," << f.GetParError(4)
      << "," << A3 << "," << f.GetParError(5)
      << "," << s3 << "," << f.GetParError(6)
      << "," << f1 << "," << f2 << "," << f3
      << "," << p0 << "," << p1 << "," << p2
      << "," << Ag << "," << mug << "," << sigg << "\n";

  txt << step << " " << label << " total: mean=" << mean
      << " s1=" << s1 << " s2=" << s2 << " s3=" << s3
      << " A1=" << A1 << " A2=" << A2 << " A3=" << A3
      << " frac=" << f1 << "/" << f2 << "/" << f3
      << " bkg(pol2)=" << p0 << "," << p1 << "," << p2
      << " peak1.00(A,mu,sig)=" << Ag << "," << mug << "," << sigg << "\n";
}

}

void export_all_fit_params(const char* step1HistFile = "Reports/step1/step1_hists.root",
                           const char* step2HistFile = "Reports/step2/step2_hists.root",
                           const char* step3HistFile = "Reports/step3/step2/step3_hists.root",
                           const char* seedCsv = "Reports/step1_best_fit_params.csv",
                           const char* outCsv = "Reports/all_step_fit_function_parameters.csv",
                           const char* outTxt = "Reports/all_step_fit_function_parameters.txt") {
  auto seeds = readSeeds(seedCsv);

  TFile f1(step1HistFile, "READ");
  TFile f2(step2HistFile, "READ");
  TFile f3(step3HistFile, "READ");

  if (f1.IsZombie() || f2.IsZombie() || f3.IsZombie()) {
    std::cerr << "Failed opening one input ROOT file" << std::endl;
    return;
  }

  auto* h10 = dynamic_cast<TH1D*>(f1.Get("h0"));
  auto* h11 = dynamic_cast<TH1D*>(f1.Get("h1"));
  auto* h12 = dynamic_cast<TH1D*>(f1.Get("h2"));

  auto* h20 = dynamic_cast<TH1D*>(f2.Get("hStep2_0tag"));
  auto* h21 = dynamic_cast<TH1D*>(f2.Get("hStep2_1tag"));
  auto* h22 = dynamic_cast<TH1D*>(f2.Get("hStep2_2tag"));

  auto* h30 = dynamic_cast<TH1D*>(f3.Get("hStep2_0tag"));
  auto* h31 = dynamic_cast<TH1D*>(f3.Get("hStep2_1tag"));
  auto* h32 = dynamic_cast<TH1D*>(f3.Get("hStep2_2tag"));

  if (!h10 || !h11 || !h12 || !h20 || !h21 || !h22 || !h30 || !h31 || !h32) {
    std::cerr << "Missing histograms in input ROOT files" << std::endl;
    return;
  }

  std::ofstream csv(outCsv);
  std::ofstream txt(outTxt);
  csv << std::setprecision(10);
  txt << std::setprecision(10);

  csv << "step,category,model,chi2ndf,mean,mean_err,A1,A1_err,sigma1,sigma1_err,A2,A2_err,sigma2,sigma2_err,A3,A3_err,sigma3,sigma3_err,frac1,frac2,frac3,bkg_p0,bkg_p1,bkg_p2,peakA,peakMean,peakSigma\n";

  fitStep1Triple(h10, "0tag", csv, txt);
  fitStep1Triple(h11, "1tag", csv, txt);
  fitStep1Triple(h12, "2tag", csv, txt);

  fitStep23(h20, seeds["0tag"], "step2_mc", "0tag", csv, txt);
  fitStep23(h21, seeds["1tag"], "step2_mc", "1tag", csv, txt);
  fitStep23(h22, seeds["2tag"], "step2_mc", "2tag", csv, txt);

  fitStep23(h30, seeds["0tag"], "step3_data", "0tag", csv, txt);
  fitStep23(h31, seeds["1tag"], "step3_data", "1tag", csv, txt);
  fitStep23(h32, seeds["2tag"], "step3_data", "2tag", csv, txt);

  csv.close();
  txt.close();

  std::cout << "Wrote: " << outCsv << std::endl;
  std::cout << "Wrote: " << outTxt << std::endl;
}
