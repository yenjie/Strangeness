#include <algorithm>
#include <cmath>
#include <vector>

#include "TCanvas.h"
#include "TDecompSVD.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLegend.h"
#include "TLine.h"
#include "TMatrixD.h"
#include "TParameter.h"
#include "TColor.h"
#include "TStyle.h"
#include "TString.h"
#include "TVectorD.h"

namespace
{
   TH1D *CollapseTail1D(const TH1D *src, int keepBins, const char *name)
   {
      keepBins = std::max(1, std::min(keepBins, src->GetNbinsX()));
      const double xmin = src->GetXaxis()->GetXmin();
      const double xmax = src->GetXaxis()->GetBinUpEdge(keepBins);
      TH1D *out = new TH1D(name, src->GetTitle(), keepBins, xmin, xmax);
      out->SetDirectory(nullptr);
      out->Sumw2();

      for (int ib = 1; ib <= src->GetNbinsX(); ++ib)
      {
         const int b = std::min(ib, keepBins);
         const double c = src->GetBinContent(ib);
         const double e = src->GetBinError(ib);
         out->SetBinContent(b, out->GetBinContent(b) + c);
         const double err2 = out->GetBinError(b) * out->GetBinError(b) + e * e;
         out->SetBinError(b, std::sqrt(std::max(0.0, err2)));
      }
      return out;
   }

   TH2D *CollapseTail2D(const TH2D *src, int keepBinsX, int keepBinsY, const char *name)
   {
      keepBinsX = std::max(1, std::min(keepBinsX, src->GetNbinsX()));
      keepBinsY = std::max(1, std::min(keepBinsY, src->GetNbinsY()));
      const double xmin = src->GetXaxis()->GetXmin();
      const double xmax = src->GetXaxis()->GetBinUpEdge(keepBinsX);
      const double ymin = src->GetYaxis()->GetXmin();
      const double ymax = src->GetYaxis()->GetBinUpEdge(keepBinsY);
      TH2D *out = new TH2D(name, src->GetTitle(), keepBinsX, xmin, xmax, keepBinsY, ymin, ymax);
      out->SetDirectory(nullptr);
      out->Sumw2();

      for (int ix = 1; ix <= src->GetNbinsX(); ++ix)
      {
         const int bx = std::min(ix, keepBinsX);
         for (int iy = 1; iy <= src->GetNbinsY(); ++iy)
         {
            const int by = std::min(iy, keepBinsY);
            const double c = src->GetBinContent(ix, iy);
            const double e = src->GetBinError(ix, iy);
            out->SetBinContent(bx, by, out->GetBinContent(bx, by) + c);
            const double err2 = out->GetBinError(bx, by) * out->GetBinError(bx, by) + e * e;
            out->SetBinError(bx, by, std::sqrt(std::max(0.0, err2)));
         }
      }
      return out;
   }

   int DetermineOverflowKeepBins(const TH1D *recoCounts, double minEvents)
   {
      const int nb = recoCounts->GetNbinsX();
      int lastRegular = 0;
      for (int ib = nb; ib >= 1; --ib)
      {
         if (recoCounts->GetBinContent(ib) >= minEvents)
         {
            lastRegular = ib;
            break;
         }
      }
      if (lastRegular <= 0)
         return nb;
      return std::min(nb, lastRegular + 1);
   }

   TH1D *FoldTruth1D(const TH1D *truth, const TH2D *respTrueReco, const char *name, const char *title)
   {
      const int nTrue = respTrueReco->GetNbinsX();
      const int nReco = respTrueReco->GetNbinsY();
      TH1D *h = new TH1D(name, title, nReco,
                         respTrueReco->GetYaxis()->GetXmin(),
                         respTrueReco->GetYaxis()->GetXmax());
      h->SetDirectory(nullptr);
      h->Sumw2();

      for (int t = 1; t <= nTrue; ++t)
      {
         const double truthVal = truth->GetBinContent(t);
         const double truthErr = truth->GetBinError(t);
         double colSum = 0.0;
         for (int r = 1; r <= nReco; ++r)
            colSum += respTrueReco->GetBinContent(t, r);
         if (colSum <= 0.0)
            continue;
         for (int r = 1; r <= nReco; ++r)
         {
            const double prob = respTrueReco->GetBinContent(t, r) / colSum;
            h->SetBinContent(r, h->GetBinContent(r) + prob * truthVal);
            const double err2 = h->GetBinError(r) * h->GetBinError(r) + prob * prob * truthErr * truthErr;
            h->SetBinError(r, std::sqrt(std::max(0.0, err2)));
         }
      }
      return h;
   }

   TH1D *BuildWeightedTruth(const TH1D *src, double slope, const char *name, const char *title)
   {
      TH1D *h = (TH1D *)src->Clone(name);
      h->SetDirectory(nullptr);
      h->SetTitle(title);
      const double xmin = h->GetXaxis()->GetXmin();
      const double xmax = h->GetXaxis()->GetXmax();
      const double center = 0.5 * (xmin + xmax);
      const double halfRange = std::max(1e-9, 0.5 * (xmax - xmin));
      const double origIntegral = h->Integral();

      for (int ib = 1; ib <= h->GetNbinsX(); ++ib)
      {
         const double x = h->GetXaxis()->GetBinCenter(ib);
         const double xnorm = (x - center) / halfRange;
         const double weight = std::max(0.25, 1.0 + slope * xnorm);
         h->SetBinContent(ib, h->GetBinContent(ib) * weight);
         h->SetBinError(ib, h->GetBinError(ib) * weight);
      }

      const double newIntegral = h->Integral();
      if (origIntegral > 0.0 && newIntegral > 0.0)
         h->Scale(origIntegral / newIntegral);
      return h;
   }

   TH2D *RebinResponseToMeasurementBinning(const TH2D *respFine, const TH1D *measTemplate, const char *name)
   {
      const int nCoarse = measTemplate->GetNbinsX();
      TH2D *out = new TH2D(name,
                           respFine->GetTitle(),
                           nCoarse,
                           measTemplate->GetXaxis()->GetXmin(),
                           measTemplate->GetXaxis()->GetXmax(),
                           nCoarse,
                           measTemplate->GetXaxis()->GetXmin(),
                           measTemplate->GetXaxis()->GetXmax());
      out->SetDirectory(nullptr);
      out->Sumw2();

      for (int ix = 1; ix <= respFine->GetNbinsX(); ++ix)
      {
         const double x = respFine->GetXaxis()->GetBinCenter(ix);
         const int bx = out->GetXaxis()->FindBin(x);
         if (bx < 1 || bx > nCoarse)
            continue;
         for (int iy = 1; iy <= respFine->GetNbinsY(); ++iy)
         {
            const double y = respFine->GetYaxis()->GetBinCenter(iy);
            const int by = out->GetYaxis()->FindBin(y);
            if (by < 1 || by > nCoarse)
               continue;
            const double w = respFine->GetBinContent(ix, iy);
            if (w == 0.0)
               continue;
            out->SetBinContent(bx, by, out->GetBinContent(bx, by) + w);
         }
      }
      return out;
   }

   TH1D *RebinPriorToMeasurementBinning(const TH1D *priorFine, const TH1D *measTemplate, const char *name)
   {
      const int nCoarse = measTemplate->GetNbinsX();
      TH1D *out = new TH1D(name,
                           priorFine->GetTitle(),
                           nCoarse,
                           measTemplate->GetXaxis()->GetXmin(),
                           measTemplate->GetXaxis()->GetXmax());
      out->SetDirectory(nullptr);
      out->Sumw2();

      for (int ib = 1; ib <= priorFine->GetNbinsX(); ++ib)
      {
         const double x = priorFine->GetXaxis()->GetBinCenter(ib);
         const int b = out->GetXaxis()->FindBin(x);
         if (b < 1 || b > nCoarse)
            continue;
         const double w = priorFine->GetBinContent(ib);
         if (w == 0.0)
            continue;
         out->SetBinContent(b, out->GetBinContent(b) + w);
      }
      return out;
   }

   TH1D *CloneEmptyLike(const TH1D *h, const char *name)
   {
      TH1D *c = (TH1D *)h->Clone(name);
      c->SetDirectory(nullptr);
      c->Reset();
      return c;
   }

   TH1D *BuildRatio(const TH1D *num, const TH1D *den, const char *name, const char *title)
   {
      TH1D *r = (TH1D *)num->Clone(name);
      r->SetDirectory(nullptr);
      r->SetTitle(title);
      r->Divide((TH1D *)den);
      return r;
   }

   TH1D *IterativeBayesUnfold1D(const TH1D *meas, const TH2D *respTrueReco, const TH1D *priorHist, int nIter, const char *name)
   {
      const int nTrue = respTrueReco->GetNbinsX();
      const int nReco = respTrueReco->GetNbinsY();
      if (meas->GetNbinsX() != nReco || priorHist->GetNbinsX() != nTrue)
      {
         Error("IterativeBayesUnfold1D", "Histogram binning mismatch");
         return CloneEmptyLike(priorHist, name);
      }

      std::vector<double> prior(nTrue, 0.0);
      double sumPrior = 0.0;
      for (int t = 1; t <= nTrue; ++t)
      {
         prior[t - 1] = std::max(0.0, priorHist->GetBinContent(t));
         sumPrior += prior[t - 1];
      }
      if (sumPrior <= 0.0)
      {
         for (int t = 0; t < nTrue; ++t)
            prior[t] = 1.0 / nTrue;
      }
      else
      {
         for (int t = 0; t < nTrue; ++t)
            prior[t] /= sumPrior;
      }

      std::vector<std::vector<double>> P(nTrue, std::vector<double>(nReco, 0.0));
      for (int t = 1; t <= nTrue; ++t)
      {
         for (int r = 1; r <= nReco; ++r)
            P[t - 1][r - 1] = std::max(0.0, respTrueReco->GetBinContent(t, r));
      }

      std::vector<double> unfolded(nTrue, 0.0), newPrior(nTrue, 0.0);
      for (int iter = 0; iter < nIter; ++iter)
      {
         std::fill(unfolded.begin(), unfolded.end(), 0.0);

         for (int r = 0; r < nReco; ++r)
         {
            const double mr = std::max(0.0, meas->GetBinContent(r + 1));
            if (mr == 0.0)
               continue;

            double norm = 0.0;
            for (int t = 0; t < nTrue; ++t)
               norm += P[t][r] * prior[t];

            if (norm <= 0.0)
               continue;

            for (int t = 0; t < nTrue; ++t)
            {
               const double w = (P[t][r] * prior[t]) / norm;
               unfolded[t] += w * mr;
            }
         }

         double s = 0.0;
         for (double v : unfolded)
            s += std::max(0.0, v);
         if (s <= 0.0)
            break;
         for (int t = 0; t < nTrue; ++t)
            newPrior[t] = std::max(0.0, unfolded[t]) / s;
         prior.swap(newPrior);
      }

      TH1D *h = CloneEmptyLike(priorHist, name);
      for (int t = 1; t <= nTrue; ++t)
      {
         h->SetBinContent(t, std::max(0.0, unfolded[t - 1]));
         h->SetBinError(t, std::sqrt(std::max(0.0, unfolded[t - 1])));
      }
      return h;
   }

   TH1D *SVDUnfold1D(const TH1D *meas, const TH2D *respTrueReco, int kReg, const char *name)
   {
      const int nTrue = respTrueReco->GetNbinsX();
      const int nReco = respTrueReco->GetNbinsY();
      if (meas->GetNbinsX() != nReco)
      {
         Error("SVDUnfold1D", "Measurement/response binning mismatch");
         return nullptr;
      }

      TMatrixD A(nReco, nTrue);
      for (int t = 1; t <= nTrue; ++t)
      {
         double colSum = 0.0;
         for (int r = 1; r <= nReco; ++r)
            colSum += respTrueReco->GetBinContent(t, r);
         if (colSum <= 0.0)
            continue;
         for (int r = 1; r <= nReco; ++r)
            A(r - 1, t - 1) = respTrueReco->GetBinContent(t, r) / colSum;
      }

      TVectorD m(nReco);
      for (int r = 1; r <= nReco; ++r)
         m(r - 1) = meas->GetBinContent(r);

      TDecompSVD svd(A);
      if (!svd.Decompose())
      {
         Error("SVDUnfold1D", "SVD decomposition failed");
         return nullptr;
      }

      TVectorD sig = svd.GetSig();
      TMatrixD U = svd.GetU();
      TMatrixD V = svd.GetV();
      const int nSig = sig.GetNrows();
      const int k = std::max(1, std::min(kReg, nSig));

      TMatrixD Splus(nTrue, nReco);
      for (int i = 0; i < k; ++i)
      {
         if (sig(i) > 1e-12)
            Splus(i, i) = 1.0 / sig(i);
      }

      TMatrixD Ut(TMatrixD::kTransposed, U);
      TMatrixD B = V * Splus * Ut;
      TVectorD x = B * m;

      TH1D *h = new TH1D(name, name, nTrue,
                         respTrueReco->GetXaxis()->GetXmin(),
                         respTrueReco->GetXaxis()->GetXmax());
      h->Reset();
      h->SetDirectory(nullptr);
      for (int t = 1; t <= nTrue; ++t)
      {
         h->SetBinContent(t, std::max(0.0, x(t - 1)));
      }

      for (int t = 1; t <= nTrue; ++t)
      {
         double var = 0.0;
         for (int r = 1; r <= nReco; ++r)
         {
            const double merr = meas->GetBinError(r);
            const double br = B(t - 1, r - 1);
            var += br * br * merr * merr;
         }
         h->SetBinError(t, std::sqrt(std::max(0.0, var)));
      }
      return h;
   }
}

void runDNdYUnfolding_BayesSVD(int nIterBayes = 1,
                                 int kRegSVD = 8,
                                 const char *mcFile = "output/KtoPi-MC-Reco-Nominal.root",
                                 const char *dataFile = "output/KtoPi-Data-Reco-Nominal.root",
                                 const char *outRoot = "output/DNdYUnfolding_BayesSVD.root",
                                 const char *responseFile = "",
                                 bool makePlots = true,
                                 int keepBinsOverride = -1)
{
   gStyle->SetOptStat(0);
   gStyle->SetPalette(kViridis);
   auto SetReadableAxisStyle = [](TH1 *h) {
      h->GetXaxis()->SetTitleSize(0.060);
      h->GetYaxis()->SetTitleSize(0.060);
      h->GetXaxis()->SetLabelSize(0.048);
      h->GetYaxis()->SetLabelSize(0.048);
      h->GetXaxis()->SetTitleOffset(1.00);
      h->GetYaxis()->SetTitleOffset(1.00);
   };
   auto SetMatrixStyle = [](TH2 *h) {
      h->GetXaxis()->SetTitleSize(0.055);
      h->GetYaxis()->SetTitleSize(0.055);
      h->GetZaxis()->SetTitleSize(0.050);
      h->GetXaxis()->SetLabelSize(0.045);
      h->GetYaxis()->SetLabelSize(0.045);
      h->GetZaxis()->SetLabelSize(0.040);
      h->GetXaxis()->SetTitleOffset(1.05);
      h->GetYaxis()->SetTitleOffset(1.22);
      h->GetZaxis()->SetTitleOffset(1.08);
   };
   auto StyleLegend = [](TLegend &leg) {
      leg.SetBorderSize(0);
      leg.SetFillColor(kWhite);
      leg.SetFillStyle(1001);
      leg.SetLineColor(kWhite);
   };

   TFile *fMC = TFile::Open(mcFile, "READ");
   TFile *fData = TFile::Open(dataFile, "READ");
   if (fMC == nullptr || fMC->IsZombie() || fData == nullptr || fData->IsZombie())
   {
      Error("runDNdYUnfolding_BayesSVD", "Cannot open input files");
      return;
   }

   TFile *fResp = fMC;
   if (responseFile != nullptr && TString(responseFile).Length() > 0)
   {
      fResp = TFile::Open(responseFile, "READ");
      if (fResp == nullptr || fResp->IsZombie())
      {
         Error("runDNdYUnfolding_BayesSVD", "Cannot open response file");
         return;
      }
   }

   TH2D *hRespFineIn = dynamic_cast<TH2D *>(fResp->Get("hDNdYResponse"));
   TH2D *hRespKFineIn = dynamic_cast<TH2D *>(fResp->Get("hDNdYResponseK"));
   TH2D *hRespPiFineIn = dynamic_cast<TH2D *>(fResp->Get("hDNdYResponsePi"));
   TH1D *hKPriorFineIn = dynamic_cast<TH1D *>(fResp->Get("hKTruedNdY"));
   TH1D *hPiPriorFineIn = dynamic_cast<TH1D *>(fResp->Get("hPiTruedNdY"));
   TH1D *hRecoCountsFine = dynamic_cast<TH1D *>(fMC->Get("hDNdYReco"));

   auto GetRecoHist = [](TFile *f, const char *preferred, const char *fallback) -> TH1D *
   {
      TH1D *h = dynamic_cast<TH1D *>(f->Get(preferred));
      if (h != nullptr)
         return h;
      return dynamic_cast<TH1D *>(f->Get(fallback));
   };

   TH1D *hKMcReco = GetRecoHist(fMC, "hKCorrectedDNdY", "hKCorrected");
   TH1D *hPiMcReco = GetRecoHist(fMC, "hPiCorrectedDNdY", "hPiCorrected");
   TH1D *hKDataReco = GetRecoHist(fData, "hKCorrectedDNdY", "hKCorrected");
   TH1D *hPiDataReco = GetRecoHist(fData, "hPiCorrectedDNdY", "hPiCorrected");

   if (hRespFineIn == nullptr || hRespKFineIn == nullptr || hRespPiFineIn == nullptr ||
       hKPriorFineIn == nullptr || hPiPriorFineIn == nullptr || hRecoCountsFine == nullptr ||
       hKMcReco == nullptr || hPiMcReco == nullptr || hKDataReco == nullptr || hPiDataReco == nullptr)
   {
      Error("runDNdYUnfolding_BayesSVD", "Missing required histograms (did you rerun ExecuteKtoPiAnalysis after adding thrust-axis dN/dy histograms?)");
      return;
   }

   const int keepBinsAuto = DetermineOverflowKeepBins(hRecoCountsFine, 100.0);
   const int keepBins = (keepBinsOverride > 0)
      ? std::max(1, std::min(keepBinsOverride, hRecoCountsFine->GetNbinsX()))
      : keepBinsAuto;
   printf("dN/dy overflow treatment: auto keepBins=%d, using keepBins=%d, collapsing bins %d..%d into final visible bin %d\n",
          keepBinsAuto, keepBins, keepBins, hRecoCountsFine->GetNbinsX(), keepBins);

   TH1D *hKMcRecoCollapsed = CollapseTail1D(hKMcReco, keepBins, "hKMcRecoCollapsed_dNdY");
   TH1D *hPiMcRecoCollapsed = CollapseTail1D(hPiMcReco, keepBins, "hPiMcRecoCollapsed_dNdY");
   TH1D *hKDataRecoCollapsed = CollapseTail1D(hKDataReco, keepBins, "hKDataRecoCollapsed_dNdY");
   TH1D *hPiDataRecoCollapsed = CollapseTail1D(hPiDataReco, keepBins, "hPiDataRecoCollapsed_dNdY");
   TH2D *hRespFine = CollapseTail2D(hRespFineIn, keepBins, keepBins, "hDNdYResponseCollapsed");
   TH2D *hRespKFine = CollapseTail2D(hRespKFineIn, keepBins, keepBins, "hDNdYResponseKCollapsed");
   TH2D *hRespPiFine = CollapseTail2D(hRespPiFineIn, keepBins, keepBins, "hDNdYResponsePiCollapsed");
   TH1D *hKPriorFine = CollapseTail1D(hKPriorFineIn, keepBins, "hKTruedNdYCollapsed");
   TH1D *hPiPriorFine = CollapseTail1D(hPiPriorFineIn, keepBins, "hPiTruedNdYCollapsed");

   hKMcReco = hKMcRecoCollapsed;
   hPiMcReco = hPiMcRecoCollapsed;
   hKDataReco = hKDataRecoCollapsed;
   hPiDataReco = hPiDataRecoCollapsed;

   TH2D *hResp = RebinResponseToMeasurementBinning(hRespFine, hKMcReco, "hDNdYResponseRebinned");
   TH2D *hRespK = RebinResponseToMeasurementBinning(hRespKFine, hKMcReco, "hDNdYResponseKRebinned");
   TH2D *hRespPi = RebinResponseToMeasurementBinning(hRespPiFine, hPiMcReco, "hDNdYResponsePiRebinned");
   TH1D *hKPrior = RebinPriorToMeasurementBinning(hKPriorFine, hKMcReco, "hKPrior_dNdY");
   TH1D *hPiPrior = RebinPriorToMeasurementBinning(hPiPriorFine, hPiMcReco, "hPiPrior_dNdY");

   TH2D *hRespNorm = (TH2D *)hResp->Clone("hDNdYResponseNormalized");
   hRespNorm->SetDirectory(nullptr);
   for (int x = 1; x <= hRespNorm->GetNbinsX(); ++x)
   {
      double sumRow = 0.0;
      for (int y = 1; y <= hRespNorm->GetNbinsY(); ++y)
         sumRow += hRespNorm->GetBinContent(x, y);
      if (sumRow > 0.0)
      {
         for (int y = 1; y <= hRespNorm->GetNbinsY(); ++y)
            hRespNorm->SetBinContent(x, y, hRespNorm->GetBinContent(x, y) / sumRow);
      }
   }

   TH1D *hKMcBayes = IterativeBayesUnfold1D(hKMcReco, hRespK, hKPrior, nIterBayes, "hKMcBayes_dNdY");
   TH1D *hPiMcBayes = IterativeBayesUnfold1D(hPiMcReco, hRespPi, hPiPrior, nIterBayes, "hPiMcBayes_dNdY");
   TH1D *hKDataBayes = IterativeBayesUnfold1D(hKDataReco, hRespK, hKPrior, nIterBayes, "hKDataBayes_dNdY");
   TH1D *hPiDataBayes = IterativeBayesUnfold1D(hPiDataReco, hRespPi, hPiPrior, nIterBayes, "hPiDataBayes_dNdY");

   TH1D *hKPriorFlat = (TH1D *)hKPrior->Clone("hKPriorFlat_dNdY");
   TH1D *hPiPriorFlat = (TH1D *)hPiPrior->Clone("hPiPriorFlat_dNdY");
   hKPriorFlat->SetDirectory(nullptr);
   hPiPriorFlat->SetDirectory(nullptr);
   for (int i = 1; i <= hKPriorFlat->GetNbinsX(); ++i)
      hKPriorFlat->SetBinContent(i, 1.0);
   for (int i = 1; i <= hPiPriorFlat->GetNbinsX(); ++i)
      hPiPriorFlat->SetBinContent(i, 1.0);

   const int nIterVar = std::max(2, nIterBayes + 1);
   TH1D *hKMcBayesPriorVar = IterativeBayesUnfold1D(hKMcReco, hRespK, hKPriorFlat, nIterBayes, "hKMcBayesPriorVar_dNdY");
   TH1D *hPiMcBayesPriorVar = IterativeBayesUnfold1D(hPiMcReco, hRespPi, hPiPriorFlat, nIterBayes, "hPiMcBayesPriorVar_dNdY");
   TH1D *hKDataBayesPriorVar = IterativeBayesUnfold1D(hKDataReco, hRespK, hKPriorFlat, nIterBayes, "hKDataBayesPriorVar_dNdY");
   TH1D *hPiDataBayesPriorVar = IterativeBayesUnfold1D(hPiDataReco, hRespPi, hPiPriorFlat, nIterBayes, "hPiDataBayesPriorVar_dNdY");

   TH1D *hKMcBayesIterVar = IterativeBayesUnfold1D(hKMcReco, hRespK, hKPrior, nIterVar, "hKMcBayesIterVar_dNdY");
   TH1D *hPiMcBayesIterVar = IterativeBayesUnfold1D(hPiMcReco, hRespPi, hPiPrior, nIterVar, "hPiMcBayesIterVar_dNdY");
   TH1D *hKDataBayesIterVar = IterativeBayesUnfold1D(hKDataReco, hRespK, hKPrior, nIterVar, "hKDataBayesIterVar_dNdY");
   TH1D *hPiDataBayesIterVar = IterativeBayesUnfold1D(hPiDataReco, hRespPi, hPiPrior, nIterVar, "hPiDataBayesIterVar_dNdY");

   TH1D *hKMcSVD = SVDUnfold1D(hKMcReco, hRespK, kRegSVD, "hKMcSVD_dNdY");
   TH1D *hPiMcSVD = SVDUnfold1D(hPiMcReco, hRespPi, kRegSVD, "hPiMcSVD_dNdY");
   TH1D *hKDataSVD = SVDUnfold1D(hKDataReco, hRespK, kRegSVD, "hKDataSVD_dNdY");
   TH1D *hPiDataSVD = SVDUnfold1D(hPiDataReco, hRespPi, kRegSVD, "hPiDataSVD_dNdY");

   TH1D *hRatioMcTrue = BuildRatio(hKPrior, hPiPrior, "hRatioMcTrue_dNdY", "MC truth K/#pi;dN_{ch}/dy (|y_{T}|<0.5);K/#pi");
   TH1D *hRatioMcBayes = BuildRatio(hKMcBayes, hPiMcBayes, "hRatioMcBayes_dNdY", "MC K/#pi Bayes;dN_{ch}/dy (|y_{T}|<0.5);K/#pi");
   TH1D *hRatioMcSVD = BuildRatio(hKMcSVD, hPiMcSVD, "hRatioMcSVD_dNdY", "MC K/#pi SVD;dN_{ch}/dy (|y_{T}|<0.5);K/#pi");
   TH1D *hRatioDataBayes = BuildRatio(hKDataBayes, hPiDataBayes, "hRatioDataBayes_dNdY", "Data K/#pi Bayes;dN_{ch}/dy (|y_{T}|<0.5);K/#pi");
   TH1D *hRatioDataSVD = BuildRatio(hKDataSVD, hPiDataSVD, "hRatioDataSVD_dNdY", "Data K/#pi SVD;dN_{ch}/dy (|y_{T}|<0.5);K/#pi");
   TH1D *hRatioMcBayesPriorVar = BuildRatio(hKMcBayesPriorVar, hPiMcBayesPriorVar, "hRatioMcBayesPriorVar_dNdY", "MC K/#pi Bayes prior-var");
   TH1D *hRatioDataBayesPriorVar = BuildRatio(hKDataBayesPriorVar, hPiDataBayesPriorVar, "hRatioDataBayesPriorVar_dNdY", "Data K/#pi Bayes prior-var");
   TH1D *hRatioMcBayesIterVar = BuildRatio(hKMcBayesIterVar, hPiMcBayesIterVar, "hRatioMcBayesIterVar_dNdY", "MC K/#pi Bayes iter-var");
   TH1D *hRatioDataBayesIterVar = BuildRatio(hKDataBayesIterVar, hPiDataBayesIterVar, "hRatioDataBayesIterVar_dNdY", "Data K/#pi Bayes iter-var");

   TH1D *hKMcBayesRefold = FoldTruth1D(hKMcBayes, hRespK, "hKMcBayesRefold_dNdY", "MC K refolded reco");
   TH1D *hPiMcBayesRefold = FoldTruth1D(hPiMcBayes, hRespPi, "hPiMcBayesRefold_dNdY", "MC #pi refolded reco");
   TH1D *hKDataBayesRefold = FoldTruth1D(hKDataBayes, hRespK, "hKDataBayesRefold_dNdY", "Data K refolded reco");
   TH1D *hPiDataBayesRefold = FoldTruth1D(hPiDataBayes, hRespPi, "hPiDataBayesRefold_dNdY", "Data #pi refolded reco");
   TH1D *hRatioMcReco = BuildRatio(hKMcReco, hPiMcReco, "hRatioMcReco_dNdY", "MC reco K/#pi;dN_{ch}^{reco}/dy (|y_{T}|<0.5);K/#pi");
   TH1D *hRatioDataReco = BuildRatio(hKDataReco, hPiDataReco, "hRatioDataReco_dNdY", "Data reco K/#pi;dN_{ch}^{reco}/dy (|y_{T}|<0.5);K/#pi");
   TH1D *hRatioMcBayesRefold = BuildRatio(hKMcBayesRefold, hPiMcBayesRefold, "hRatioMcBayesRefold_dNdY", "MC K/#pi refolded");
   TH1D *hRatioDataBayesRefold = BuildRatio(hKDataBayesRefold, hPiDataBayesRefold, "hRatioDataBayesRefold_dNdY", "Data K/#pi refolded");
   TH1D *hRefoldRecoClosureMc = BuildRatio(hRatioMcBayesRefold, hRatioMcReco, "hRefoldRecoClosureMc_dNdY", ";dN_{ch}^{reco}/dy (|y_{T}|<0.5);Refolded / reco");
   TH1D *hRefoldRecoClosureData = BuildRatio(hRatioDataBayesRefold, hRatioDataReco, "hRefoldRecoClosureData_dNdY", ";dN_{ch}^{reco}/dy (|y_{T}|<0.5);Refolded / reco");

   TH1D *hKStressTruth = BuildWeightedTruth(hKPrior, +0.18, "hKStressTruth_dNdY", "Injected K truth");
   TH1D *hPiStressTruth = BuildWeightedTruth(hPiPrior, -0.12, "hPiStressTruth_dNdY", "Injected #pi truth");
   TH1D *hKStressReco = FoldTruth1D(hKStressTruth, hRespK, "hKStressReco_dNdY", "Injected K reco");
   TH1D *hPiStressReco = FoldTruth1D(hPiStressTruth, hRespPi, "hPiStressReco_dNdY", "Injected #pi reco");
   TH1D *hKStressUnfold = IterativeBayesUnfold1D(hKStressReco, hRespK, hKPrior, nIterBayes, "hKStressUnfold_dNdY");
   TH1D *hPiStressUnfold = IterativeBayesUnfold1D(hPiStressReco, hRespPi, hPiPrior, nIterBayes, "hPiStressUnfold_dNdY");
   TH1D *hRatioStressTruth = BuildRatio(hKStressTruth, hPiStressTruth, "hRatioStressTruth_dNdY", "Injected truth K/#pi");
   TH1D *hRatioStressReco = BuildRatio(hKStressReco, hPiStressReco, "hRatioStressReco_dNdY", "Injected reco K/#pi");
   TH1D *hRatioStressUnfold = BuildRatio(hKStressUnfold, hPiStressUnfold, "hRatioStressUnfold_dNdY", "Injected unfolded K/#pi");
   TH1D *hStressClosure = BuildRatio(hRatioStressUnfold, hRatioStressTruth, "hStressClosure_dNdY", ";dN_{ch}/dy (|y_{T}|<0.5);Unfolded / injected");

   TH1D *hClosureBayes = BuildRatio(hRatioMcBayes, hRatioMcTrue, "hClosureBayes_dNdY", "Bayes closure;dN_{ch}/dy (|y_{T}|<0.5);Unfolded / MC truth");
   TH1D *hClosureSVD = BuildRatio(hRatioMcSVD, hRatioMcTrue, "hClosureSVD_dNdY", "SVD closure;dN_{ch}/dy (|y_{T}|<0.5);Unfolded / MC truth");

   TH1D *hDataOverMcBayes = BuildRatio(hRatioDataBayes, hRatioMcBayes, "hDataOverMcBayes_dNdY", "(K/#pi)_{Data}/(K/#pi)_{MC};dN_{ch}/dy (|y_{T}|<0.5);Double ratio");
   TH1D *hDataOverMcSVD = BuildRatio(hRatioDataSVD, hRatioMcSVD, "hDataOverMcSVD_dNdY", "(K/#pi)_{Data}/(K/#pi)_{MC};dN_{ch}/dy (|y_{T}|<0.5);Double ratio");
   TH1D *hDataOverMcBayesPriorVar = BuildRatio(hRatioDataBayesPriorVar, hRatioMcBayesPriorVar, "hDataOverMcBayesPriorVar_dNdY", "(K/#pi)_{Data}/(K/#pi)_{MC};dN_{ch}/dy (|y_{T}|<0.5);Double ratio");
   TH1D *hDataOverMcBayesIterVar = BuildRatio(hRatioDataBayesIterVar, hRatioMcBayesIterVar, "hDataOverMcBayesIterVar_dNdY", "(K/#pi)_{Data}/(K/#pi)_{MC};dN_{ch}/dy (|y_{T}|<0.5);Double ratio");

   TH1D *hMethodDiff = (TH1D *)hDataOverMcBayes->Clone("hMethodDiff_BayesMinusSVD_dNdY");
   hMethodDiff->SetDirectory(nullptr);
   hMethodDiff->Add(hDataOverMcSVD, -1.0);
   hMethodDiff->SetTitle("Unfolding method difference (Bayes-SVD);dN_{ch}/dy (|y_{T}|<0.5);#Delta double ratio");
   TH1D *hPriorDiff = (TH1D *)hDataOverMcBayesPriorVar->Clone("hBayesPriorVariationDiff_dNdY");
   hPriorDiff->SetDirectory(nullptr);
   hPriorDiff->Add(hDataOverMcBayes, -1.0);
   hPriorDiff->SetTitle("Bayes prior variation (flat prior - nominal);dN_{ch}/dy (|y_{T}|<0.5);#Delta double ratio");
   TH1D *hIterDiff = (TH1D *)hDataOverMcBayesIterVar->Clone("hBayesIterVariationDiff_dNdY");
   hIterDiff->SetDirectory(nullptr);
   hIterDiff->Add(hDataOverMcBayes, -1.0);
   hIterDiff->SetTitle("Bayes iteration variation (nIter var - nominal);dN_{ch}/dy (|y_{T}|<0.5);#Delta double ratio");

   if (makePlots)
   {
      TCanvas *cResp = new TCanvas("cDNdYResp", "dN/dy response", 780, 680);
      cResp->SetLeftMargin(0.14);
      cResp->SetRightMargin(0.18);
      cResp->SetBottomMargin(0.14);
      cResp->SetTopMargin(0.08);
      SetMatrixStyle(hRespNorm);
      hRespNorm->Draw("COLZ");
      cResp->SaveAs("DNdYUnfolding_ResponseMatrix.pdf");
      cResp->SaveAs("DNdYUnfolding_ResponseMatrix.png");

      TCanvas *cC = new TCanvas("cDNdYClosure", "dN/dy closure", 860, 620);
      hClosureBayes->SetMarkerStyle(20);
      hClosureBayes->SetMarkerColor(kBlue + 1);
      hClosureBayes->SetLineColor(kBlue + 1);
      hClosureSVD->SetMarkerStyle(21);
      hClosureSVD->SetMarkerColor(kGreen + 2);
      hClosureSVD->SetLineColor(kGreen + 2);
      hClosureBayes->SetMarkerSize(1.15);
      hClosureSVD->SetMarkerSize(1.15);
      hClosureBayes->SetLineWidth(2);
      hClosureSVD->SetLineWidth(2);
      SetReadableAxisStyle(hClosureBayes);
      hClosureBayes->SetMinimum(0.6);
      hClosureBayes->SetMaximum(1.4);
      hClosureBayes->Draw("E1");
      TLine l1(hClosureBayes->GetXaxis()->GetXmin(), 1.0, hClosureBayes->GetXaxis()->GetXmax(), 1.0);
      l1.SetLineStyle(2);
      l1.SetLineWidth(2);
      l1.Draw("SAME");
      hClosureBayes->Draw("E1 SAME");
      hClosureSVD->Draw("E1 SAME");
      TLegend legC(0.58, 0.74, 0.89, 0.89);
      StyleLegend(legC);
      legC.AddEntry(hClosureBayes, "Bayes", "lep");
      legC.AddEntry(hClosureSVD, "SVD", "lep");
      legC.Draw();
      cC->SaveAs("DNdYUnfolding_MCClosure_BayesVsSVD.pdf");
      cC->SaveAs("DNdYUnfolding_MCClosure_BayesVsSVD.png");

      TCanvas *cD = new TCanvas("cDNdYDataMC", "dN/dy data/mc", 860, 620);
      hDataOverMcBayes->SetMarkerStyle(20);
      hDataOverMcBayes->SetMarkerColor(kBlue + 1);
      hDataOverMcBayes->SetLineColor(kBlue + 1);
      hDataOverMcSVD->SetMarkerStyle(21);
      hDataOverMcSVD->SetMarkerColor(kGreen + 2);
      hDataOverMcSVD->SetLineColor(kGreen + 2);
      hDataOverMcBayes->SetMarkerSize(1.15);
      hDataOverMcSVD->SetMarkerSize(1.15);
      hDataOverMcBayes->SetLineWidth(2);
      hDataOverMcSVD->SetLineWidth(2);
      SetReadableAxisStyle(hDataOverMcBayes);
      hDataOverMcBayes->SetMinimum(0.6);
      hDataOverMcBayes->SetMaximum(1.4);
      hDataOverMcBayes->Draw("E1");
      TLine l2(hDataOverMcBayes->GetXaxis()->GetXmin(), 1.0, hDataOverMcBayes->GetXaxis()->GetXmax(), 1.0);
      l2.SetLineStyle(2);
      l2.SetLineWidth(2);
      l2.Draw("SAME");
      hDataOverMcBayes->Draw("E1 SAME");
      hDataOverMcSVD->Draw("E1 SAME");
      TLegend legD(0.58, 0.74, 0.89, 0.89);
      StyleLegend(legD);
      legD.AddEntry(hDataOverMcBayes, "Bayes", "lep");
      legD.AddEntry(hDataOverMcSVD, "SVD", "lep");
      legD.Draw();
      cD->SaveAs("DNdYUnfolding_DataMC_BayesVsSVD.pdf");
      cD->SaveAs("DNdYUnfolding_DataMC_BayesVsSVD.png");

      TCanvas *cM = new TCanvas("cDNdYMethodDiff", "method diff", 860, 620);
      hMethodDiff->SetMarkerStyle(20);
      hMethodDiff->SetMarkerSize(1.15);
      hMethodDiff->SetLineWidth(2);
      SetReadableAxisStyle(hMethodDiff);
      hMethodDiff->Draw("E1");
      TLine l0(hMethodDiff->GetXaxis()->GetXmin(), 0.0, hMethodDiff->GetXaxis()->GetXmax(), 0.0);
      l0.SetLineStyle(2);
      l0.SetLineWidth(2);
      l0.Draw("SAME");
      cM->SaveAs("DNdYUnfolding_MethodDifference.pdf");
      cM->SaveAs("DNdYUnfolding_MethodDifference.png");
   }

   TFile *fout = TFile::Open(outRoot, "RECREATE");
   TParameter<int> pKeepBinsAuto("keepBinsAuto_dNdY", keepBinsAuto);
   TParameter<int> pKeepBinsUsed("keepBinsUsed_dNdY", keepBins);
   pKeepBinsAuto.Write();
   pKeepBinsUsed.Write();
   hRespNorm->Write();
   hResp->Write();
   hRespK->Write();
   hRespPi->Write();
   hKPrior->Write();
   hPiPrior->Write();
   hKMcReco->Write("hKMcReco");
   hPiMcReco->Write("hPiMcReco");
   hKDataReco->Write("hKDataReco");
   hPiDataReco->Write("hPiDataReco");

   hKMcBayes->Write();
   hPiMcBayes->Write();
   hKMcSVD->Write();
   hPiMcSVD->Write();
   hKDataBayes->Write();
   hPiDataBayes->Write();
   hKDataSVD->Write();
   hPiDataSVD->Write();
   hKMcBayesRefold->Write();
   hPiMcBayesRefold->Write();
   hKDataBayesRefold->Write();
   hPiDataBayesRefold->Write();
   hKStressTruth->Write();
   hPiStressTruth->Write();
   hKStressReco->Write();
   hPiStressReco->Write();
   hKStressUnfold->Write();
   hPiStressUnfold->Write();

   hRatioMcTrue->Write();
   hRatioMcReco->Write();
   hRatioMcBayes->Write();
   hRatioMcSVD->Write();
   hRatioMcBayesRefold->Write();
   hRatioDataReco->Write();
   hRatioDataBayes->Write();
   hRatioDataSVD->Write();
   hRatioDataBayesRefold->Write();
   hRefoldRecoClosureMc->Write();
   hRefoldRecoClosureData->Write();
   hRatioStressTruth->Write();
   hRatioStressReco->Write();
   hRatioStressUnfold->Write();
   hStressClosure->Write();

   hClosureBayes->Write();
   hClosureSVD->Write();
   hDataOverMcBayes->Write();
   hDataOverMcSVD->Write();
   hDataOverMcBayesPriorVar->Write();
   hDataOverMcBayesIterVar->Write();
   hMethodDiff->Write();
    hPriorDiff->Write();
    hIterDiff->Write();
   fout->Close();

   fMC->Close();
   fData->Close();
   if (fResp != fMC)
      fResp->Close();

   printf("Saved dN/dy unfolding output to %s\n", outRoot);
}
