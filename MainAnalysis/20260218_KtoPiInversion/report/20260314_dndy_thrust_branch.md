# 2026-03-14 thrust-axis dN/dy branch

## Scope
- Added a parallel unfolded `K/pi` branch versus `dN_ch/dy` in `|y_T|<0.5`.
- `y_T` is defined with respect to the event thrust axis:
  - `y_T = 0.5 * ln((E + p·n_thrust) / (E - p·n_thrust))`
- This branch does not replace the existing beam-axis `dN_ch/deta` result.

## Code updates
- Analyzer:
  - `KtoPiAnalysis.cpp`
  - new raw/corrected histograms, response matrices, truth priors, and pT-vs-activity histograms for `DNdY`
  - event activity counted from charged tracks or particles with `|y_T| < 0.5`
  - reco event-count histograms `hNtagReco`, `hDNdEtaReco`, and `hDNdYReco` are now filled for data as well as MC
- Unfolding macro:
  - `runDNdYUnfolding_BayesSVD.C`
- Driver:
  - `run_systematics_dndy_unfolding.py`
- Finalizer:
  - `finalize_dndy_systematics.py`
- Plot scripts:
  - `make_kpi_vs_dndy_comparison_plots.py`
  - `make_dndy_thrust_control_plot.py`

## Production outputs
- Inputs:
  - `output/systematics_20260314_dndy_inputs/`
- Unfolded outputs:
  - `output/systematics_20260314_dndy/`
- Top plots:
  - `result/20260314/top_plots_dndy/`

## Variations run
- `nominal`
- `npt10`
- `npt14`

These match the currently quoted `dN_ch/deta` branch strategy: binning envelope plus residual, unfolding, and PID-SF propagation.

## Statistical treatment
- No dedicated `dN/dy` toy-coverage calibration was run in this update.
- The current `dN/dy` note figures therefore use the nominal Bayes propagated bin errors.
- This is recorded in:
  - `output/systematics_20260314_dndy/systematics_dndy_table.txt`
  - `output/systematics_20260314_dndy/systematics_dndy_summary.root` via `dNdYStatTreatment`

## Main figures
- Control:
  - `result/20260314/top_plots_dndy/ControlPlot_DNdYThrust_DataMC.pdf`
- Final double ratio:
  - `output/systematics_20260314_dndy/KtoPi_DoubleRatio_vs_dNdY_with_Systematics.pdf`
- Unfolded `K/pi` data vs MC:
  - `result/20260314/top_plots_dndy/KtoPi_vs_dNdY_DataMC_with_Systematics.pdf`
- Unfolding validation:
  - `result/20260314/top_plots_dndy/DNdYUnfolding_ResponseMatrix.pdf`
  - `result/20260314/top_plots_dndy/DNdYUnfolding_MethodDifference.pdf`
  - `result/20260314/top_plots_dndy/DNdYUnfolding_MCClosure_BayesVsSVD.pdf`
  - `result/20260314/top_plots_dndy/DNdYUnfolding_DataMC_BayesVsSVD.pdf`

## Note integration
- Added a dedicated results subsection to:
  - `overleaf_repo_git/main.tex`
- Synced note figure assets into:
  - `overleaf_repo_git/StrangenessV3/Figures/`
- Rebuilt note:
  - `overleaf_repo_git/main.pdf`

## Current limitations
- The highest-activity `dN_ch/dy` bin is still the merged overflow bin from `keepBins=8`.
- The statistical treatment is weaker than the main `dN_ch/deta` branch until a dedicated `dN/dy` toy study is added.
- There is no external heavy-ion or pp comparison overlay for the thrust-axis variable, because those experiments do not report the same observable definition.
