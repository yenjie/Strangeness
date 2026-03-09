#include <TCanvas.h>
#include <TPaveText.h>
#include <TROOT.h>
#include <TStyle.h>

#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace {

std::vector<std::string> splitCsv(const std::string& s) {
  std::vector<std::string> out;
  std::stringstream ss(s);
  std::string item;
  while (std::getline(ss, item, ',')) out.push_back(item);
  return out;
}

struct Step1Best {
  std::string category;
  std::string model;
  double chi2ndf = 0.0;
  double mean = 0.0;
  double meanErr = 0.0;
  double sigma1 = 0.0;
  double sigma1Err = 0.0;
  double sigma2 = 0.0;
  double sigma2Err = 0.0;
  double sigma3 = 0.0;
  double sigma3Err = 0.0;
};

struct Step23Row {
  std::string category;
  std::string model;
  double chi2ndf = 0.0;
  double mean = 0.0;
  double meanErr = 0.0;
  double yield = 0.0;
  double yieldErr = 0.0;
};

struct Step23Summary {
  std::string input;
  double efficiency = 0.0;
  double efficiencyErr = 0.0;
  std::map<std::string, Step23Row> rows;
};

std::map<std::string, Step1Best> readStep1Best(const std::string& path) {
  std::ifstream in(path);
  std::map<std::string, Step1Best> out;
  if (!in) {
    std::cerr << "Cannot open: " << path << std::endl;
    return out;
  }

  std::string line;
  while (std::getline(in, line)) {
    if (line.empty()) continue;
    auto cols = splitCsv(line);
    if (cols.size() < 12) continue;
    if (cols[0].find("-best") == std::string::npos) continue;

    Step1Best b;
    b.category = cols[0].substr(0, cols[0].find("-best"));
    b.model = cols[1];

    for (size_t i = 2; i + 1 < cols.size(); i += 2) {
      const std::string& k = cols[i];
      const double v = std::atof(cols[i + 1].c_str());
      if (k == "chi2ndf") b.chi2ndf = v;
      if (k == "mean") b.mean = v;
      if (k == "meanErr") b.meanErr = v;
      if (k == "width") b.sigma1 = v;
      if (k == "widthErr") b.sigma1Err = v;
      if (k == "sigma2") b.sigma2 = v;
      if (k == "sigma2Err") b.sigma2Err = v;
      if (k == "sigma3") b.sigma3 = v;
      if (k == "sigma3Err") b.sigma3Err = v;
    }
    out[b.category] = b;
  }

  return out;
}

Step23Summary readStep23(const std::string& path) {
  Step23Summary s;
  std::ifstream in(path);
  if (!in) {
    std::cerr << "Cannot open: " << path << std::endl;
    return s;
  }

  std::string line;
  while (std::getline(in, line)) {
    if (line.empty()) continue;
    auto cols = splitCsv(line);
    if (cols.empty()) continue;

    if (cols[0] == "Input" && cols.size() > 1) {
      s.input = cols[1];
      continue;
    }

    if (cols[0] == "tagging_efficiency" || cols[0] == "MC_efficiency") {
      if (cols.size() >= 4) {
        s.efficiency = std::atof(cols[1].c_str());
        s.efficiencyErr = std::atof(cols[3].c_str());
      }
      continue;
    }

    if ((cols[0] == "0tag" || cols[0] == "1tag" || cols[0] == "2tag") && cols.size() >= 14) {
      Step23Row r;
      r.category = cols[0];
      r.model = cols[2];
      r.chi2ndf = std::atof(cols[4].c_str());
      r.mean = std::atof(cols[6].c_str());
      r.meanErr = std::atof(cols[8].c_str());
      r.yield = std::atof(cols[10].c_str());
      r.yieldErr = std::atof(cols[12].c_str());
      s.rows[r.category] = r;
    }
  }

  return s;
}

void addLine(TPaveText& t, const std::string& s) { t.AddText(s.c_str()); }

}  // namespace

void report_fit_params_step123(const char* step1Path = "Reports/step1_fit_parameters.txt",
                               const char* step2Path = "Reports/step2_fit_summary.txt",
                               const char* step3Path = "Reports/step3/step3_fit_summary.txt",
                               const char* outPdf = "Reports/fit_parameter_summary_step1_step2_step3.pdf") {
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);

  auto step1 = readStep1Best(step1Path);
  auto step2 = readStep23(step2Path);
  auto step3 = readStep23(step3Path);

  TCanvas c("c_summary", "c_summary", 1100, 800);
  c.Print((std::string(outPdf) + "[").c_str());

  c.Clear();
  TPaveText cover(0.05, 0.06, 0.95, 0.94, "NDC");
  cover.SetBorderSize(0);
  cover.SetFillStyle(0);
  cover.SetTextAlign(12);
  addLine(cover, "Fit Parameter Summary: Step 1, Step 2 (MC), Step 3 (Data)");
  addLine(cover, " ");
  addLine(cover, std::string("Step1 source: ") + step1Path);
  addLine(cover, std::string("Step2 source: ") + step2Path);
  addLine(cover, std::string("Step3 source: ") + step3Path);
  addLine(cover, " ");
  addLine(cover, "Tag working point: RecoPIDKaon >= 2");
  addLine(cover, "Fit range: 1.00 - 1.05 GeV");
  cover.Draw();
  c.Print(outPdf);

  c.Clear();
  TPaveText p1(0.04, 0.04, 0.96, 0.96, "NDC");
  p1.SetBorderSize(0);
  p1.SetFillStyle(0);
  p1.SetTextAlign(12);
  addLine(p1, "Step 1: Best Signal Shape Parameters (MC gen-matched phi)");
  addLine(p1, " ");
  addLine(p1, "Category | Model | chi2/ndf | mean +- err [GeV] | sigma1 +- err | sigma2 +- err | sigma3 +- err");
  for (const auto& cat : {std::string("0tag"), std::string("1tag"), std::string("2tag")}) {
    auto it = step1.find(cat);
    if (it == step1.end()) continue;
    const auto& b = it->second;
    std::ostringstream os;
    os << std::fixed << std::setprecision(6)
       << b.category << " | " << b.model
       << " | " << std::setprecision(3) << b.chi2ndf
       << " | " << std::setprecision(6) << b.mean << " +- " << b.meanErr
       << " | " << b.sigma1 << " +- " << b.sigma1Err
       << " | " << b.sigma2 << " +- " << b.sigma2Err
       << " | " << b.sigma3 << " +- " << b.sigma3Err;
    addLine(p1, os.str());
  }
  p1.Draw();
  c.Print(outPdf);

  c.Clear();
  TPaveText p2(0.04, 0.04, 0.96, 0.96, "NDC");
  p2.SetBorderSize(0);
  p2.SetFillStyle(0);
  p2.SetTextAlign(12);
  addLine(p2, "Step 2: MC Random-Pair Fit Results");
  addLine(p2, std::string("Input: ") + step2.input);
  addLine(p2, " ");
  addLine(p2, "Category | Model | chi2/ndf | mean +- err [GeV] | signal yield +- err");
  for (const auto& cat : {std::string("0tag"), std::string("1tag"), std::string("2tag")}) {
    auto it = step2.rows.find(cat);
    if (it == step2.rows.end()) continue;
    const auto& r = it->second;
    std::ostringstream os;
    os << std::fixed << std::setprecision(6)
       << r.category << " | " << r.model
       << " | " << std::setprecision(3) << r.chi2ndf
       << " | " << std::setprecision(6) << r.mean << " +- " << r.meanErr
       << " | " << std::setprecision(1) << r.yield << " +- " << r.yieldErr;
    addLine(p2, os.str());
  }
  {
    std::ostringstream os;
    os << std::fixed << std::setprecision(6)
       << "Tagging efficiency (MC): " << step2.efficiency << " +- " << step2.efficiencyErr;
    addLine(p2, " ");
    addLine(p2, os.str());
  }
  p2.Draw();
  c.Print(outPdf);

  c.Clear();
  TPaveText p3(0.04, 0.04, 0.96, 0.96, "NDC");
  p3.SetBorderSize(0);
  p3.SetFillStyle(0);
  p3.SetTextAlign(12);
  addLine(p3, "Step 3: Data Random-Pair Fit Results");
  addLine(p3, std::string("Input: ") + step3.input);
  addLine(p3, " ");
  addLine(p3, "Category | Model | chi2/ndf | mean +- err [GeV] | signal yield +- err");
  for (const auto& cat : {std::string("0tag"), std::string("1tag"), std::string("2tag")}) {
    auto it = step3.rows.find(cat);
    if (it == step3.rows.end()) continue;
    const auto& r = it->second;
    std::ostringstream os;
    os << std::fixed << std::setprecision(6)
       << r.category << " | " << r.model
       << " | " << std::setprecision(3) << r.chi2ndf
       << " | " << std::setprecision(6) << r.mean << " +- " << r.meanErr
       << " | " << std::setprecision(1) << r.yield << " +- " << r.yieldErr;
    addLine(p3, os.str());
  }
  {
    std::ostringstream os;
    os << std::fixed << std::setprecision(6)
       << "Tagging efficiency (Data): " << step3.efficiency << " +- " << step3.efficiencyErr;
    addLine(p3, " ");
    addLine(p3, os.str());
  }
  p3.Draw();
  c.Print(outPdf);

  c.Clear();
  TPaveText p4(0.04, 0.04, 0.96, 0.96, "NDC");
  p4.SetBorderSize(0);
  p4.SetFillStyle(0);
  p4.SetTextAlign(12);
  addLine(p4, "Combined Summary");
  addLine(p4, " ");
  double sf = 0.0;
  double sfErr = 0.0;
  if (step2.efficiency > 0.0 && step3.efficiency > 0.0) {
    sf = step3.efficiency / step2.efficiency;
    const double rel2 = std::pow(step2.efficiencyErr / step2.efficiency, 2) +
                        std::pow(step3.efficiencyErr / step3.efficiency, 2);
    sfErr = sf * std::sqrt(rel2);
  }
  {
    std::ostringstream os1;
    os1 << std::fixed << std::setprecision(6)
        << "Efficiency (MC):   " << step2.efficiency << " +- " << step2.efficiencyErr;
    addLine(p4, os1.str());

    std::ostringstream os2;
    os2 << std::fixed << std::setprecision(6)
        << "Efficiency (Data): " << step3.efficiency << " +- " << step3.efficiencyErr;
    addLine(p4, os2.str());

    std::ostringstream os3;
    os3 << std::fixed << std::setprecision(6)
        << "Data/MC scale factor: " << sf << " +- " << sfErr;
    addLine(p4, os3.str());
  }
  p4.Draw();
  c.Print(outPdf);

  c.Print((std::string(outPdf) + "]").c_str());
  std::cout << "Created: " << outPdf << std::endl;
}
