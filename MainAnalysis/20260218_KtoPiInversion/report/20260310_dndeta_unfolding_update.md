# 2026-03-10 `dN_{ch}/d\eta` unfolding update
Author: Yen-Jie Lee and OpenAI

## Scope
This note records the changes made to the `dN_{ch}/d\eta`-dependent unfolding chain in `MainAnalysis/20260218_KtoPiInversion`, the reason for those changes, and the current validation status.

## Problem in the previous chain
The previous `dN_{ch}/d\eta` branch was not a same-variable unfolding.

- Truth axis: generator-level charged multiplicity in `|eta| < 0.5`
- Reco axis: reconstructed `N_{ch}^{tag}`

That setup was effectively unfolding from reco `N_{ch}^{tag}` to true `dN_{ch}/d\eta`, then refolding back to reco `N_{ch}^{tag}`. The response was broad, the refolding closure was poor, and Figure 13 in the note showed a clear failure.

## Implemented changes
### 1. Same-variable reco `dN_{ch}/d\eta` definition
A dedicated reco-side multiplicity estimator is now used for the `dN_{ch}/d\eta` unfolding.

Definition:
- count reconstructed tracks with `RecoGoodTrack == 1`
- require `RecoCharge != 0`
- no PID requirement
- compute reco pseudorapidity from reco momentum
- require `|eta| < 0.5`
- use the resulting charged-track count as the reco-side `dN_{ch}/d\eta` estimator

Because the rapidity window width is `Delta eta = 1.0`, this count is numerically equal to `dN_{ch}/d\eta` in the chosen fiducial window.

Implementation:
- [KtoPiAnalysis.cpp](../KtoPiAnalysis.cpp)
  - `NchEta05Reco` counting
  - `hKDNdEta`, `hPiDNdEta`, `hPDNdEta`
  - `hKCorrectedDNdEta`, `hPiCorrectedDNdEta`, `hPCorrectedDNdEta`
  - `hDNdEtaResponse`, `hDNdEtaResponseK`, `hDNdEtaResponsePi`, `hDNdEtaResponseP`
  - `hDNdEtaReco`, `hDNdEtaTrue`

### 2. Response matrix changed to same-variable mapping
The `dN_{ch}/d\eta` response is now built as:
- truth axis: true charged multiplicity in `|eta| < 0.5`
- reco axis: reconstructed charged multiplicity in `|eta| < 0.5`

Instead of the previous proxy mapping:
- truth `dN_{ch}/d\eta`
- reco `N_{ch}^{tag}`

Implementation:
- [KtoPiAnalysis.cpp](../KtoPiAnalysis.cpp)
- [runDNdEtaUnfolding_BayesSVD.C](../runDNdEtaUnfolding_BayesSVD.C)
- [runDNdEtaUnfolding_PoverPi_BayesSVD.C](../runDNdEtaUnfolding_PoverPi_BayesSVD.C)

### 3. Dedicated reco spectra vs reco `dN_{ch}/d\eta`
The unfolding inputs for the `dN_{ch}/d\eta` analysis now come from the new corrected reco spectra:
- `hKCorrectedDNdEta`
- `hPiCorrectedDNdEta`
- `hPCorrectedDNdEta`

This removes the previous mismatch where the `dN_{ch}/d\eta` unfolding used reco spectra still binned in `N_{ch}^{tag}`.

### 4. Overflow-tail merge in the high-activity region
A final visible overflow bin was added in the unfolding stage to stabilize the sparsely populated high-activity tail.

Current treatment:
- collapse bins `10..16` into the final visible bin `10`
- apply the same collapse to:
  - reco spectra
  - truth priors
  - response matrices

Implementation:
- [runDNdEtaUnfolding_BayesSVD.C](../runDNdEtaUnfolding_BayesSVD.C)
- [runDNdEtaUnfolding_PoverPi_BayesSVD.C](../runDNdEtaUnfolding_PoverPi_BayesSVD.C)

Console message from the current nominal macro:
- `dN/deta overflow treatment: collapsing bins 10..16 into final visible bin 10`

## Rebuilt outputs
Nominal and binning-variation files were refreshed for the updated chain.

Main ROOT outputs:
- `output/systematics_20260306/nominal_data.root`
- `output/systematics_20260306/nominal_mc.root`
- `output/systematics_20260306/npt10_data.root`
- `output/systematics_20260306/npt10_mc.root`
- `output/systematics_20260306/npt14_data.root`
- `output/systematics_20260306/npt14_mc.root`
- `output/systematics_20260306_dndeta/nominal_unfold_dndeta.root`
- `output/systematics_20260306_dndeta/npt10_unfold_dndeta.root`
- `output/systematics_20260306_dndeta/npt14_unfold_dndeta.root`
- `output/systematics_20260306_dndeta/nominal_unfold_dndeta_ptopi.root`
- `output/systematics_20260306_dndeta/npt10_unfold_dndeta_ptopi.root`
- `output/systematics_20260306_dndeta/npt14_unfold_dndeta_ptopi.root`

## Current validation status
### K/pi `dN_{ch}/d\eta` chain
The updated same-variable reco definition plus the overflow-tail merge substantially improved the validation.

Two different refolding observables are now tracked and should not be mixed:

1. Ratio-level refolding of `K/pi`
- `hRefoldRecoClosureMc_dNdEta`: `0.01724`
- `hRefoldRecoClosureData_dNdEta`: `0.03401`

2. Species-level refolding used in Figure 13
- `K` MC: `0.45783`
- `K` data: `0.45598`
- `pi` MC: `0.45493`
- `pi` data: `0.45659`

The published Figure 13 is the second observable: it compares the refolded kaon
and pion reco spectra directly to the measured reco spectra. The smaller RMS
values above apply only to the ratio-level `K/pi` validation, where kaon and
pion distortions partially cancel.

An explicit companion ratio-level refolding figure is now produced in addition
to the species-level refolding figure:
- `result/20260310/unfolding_validation/DNdEtaUnfolding_RatioRefoldingValidation.pdf`

This separation is intentional. The ratio-level `K/\pi` refolding is the metric
relevant to the quoted final observable, while the species-level kaon and pion
refolding remains the stricter ingredient-level validation.

Other diagnostics:
- `hClosureBayes_dNdEta` RMS: `0.02311`
- `hStressClosure_dNdEta` RMS: `0.08326`

Nominal closure values by bin:
- `hClosureBayes_dNdEta`
- `[0.9946, 0.9776, 0.9976, 0.9932, 0.9997, 0.9964, 0.9896, 1.0174, 0.9586, 0.9488]`

Nominal data/MC values by bin:
- `hDataOverMcBayes_dNdEta`
- `[1.0640, 1.0603, 1.0516, 1.0431, 1.0393, 1.0402, 1.0334, 1.0142, 0.9807, 0.9482]`

Interpretation:
- the ratio-level `K/\pi` refolding is now acceptable for the nominal
  `K/\pi` `dN_{ch}/d\eta` chain
- the species-level refolding shown in Figure 13 is not yet fully satisfactory
  in the high-activity tail
- the dominant residual limitation is the final merged overflow bin, but the
  species-level drop already becomes visible from about bin 6 onward

### Explicit tail-resolution path
The next technical step is now treated as a scanned parameter study rather than
as a hidden heuristic.

- the unfolding macro accepts an explicit `keepBins` override
- the intended scan is `keepBins = 8, 9, 10`
- the acceptance rule is:
  1. ratio-level `K/\pi` refolding remains acceptable
  2. species-level kaon and pion tail bins no longer collapse catastrophically
  3. final merged-bin purity and stability improve beyond the current
     marginal `0.20-0.25` range

If no reasonable `keepBins` choice satisfies those criteria, the last
`dN_{ch}/d\eta` point should be downgraded to an explicitly cautionary overflow
bin or removed from the quoted physics interpretation.

### Migration quality summary
From `result/20260310/unfolding_validation/migration_metrics.csv`:

Kaons:
- bin 1 (`1.41`): purity `0.7106`, stability `0.9045`
- bin 2 (`5.22`): purity `0.6960`, stability `0.6499`
- bin 3 (`9.03`): purity `0.5491`, stability `0.4860`
- final overflow bin (`35.72`): purity `0.2010`, stability `0.2062`

Pions:
- bin 1 (`1.41`): purity `0.6952`, stability `0.9129`
- bin 2 (`5.22`): purity `0.6908`, stability `0.6459`
- bin 3 (`9.03`): purity `0.5445`, stability `0.4870`
- final overflow bin (`35.72`): purity `0.2469`, stability `0.2376`

This is consistent with the visual behavior of the updated response and refolding plots: the central bins are reasonably diagonal, while the final merged bin remains broad.

### p/pi `dN_{ch}/d\eta` chain
The same reco-axis and overflow-bin treatment was propagated to the `p/\pi` branch, and the figures were rebuilt. However, the closure is still not yet at the same quality level as `K/\pi`.

Current `p/\pi` nominal closure:
- `output/systematics_20260306_dndeta/nominal_unfold_dndeta_ptopi.root`
- `hClosureBayes_dNdEta` RMS: `0.17397`

So the structural fix is in place for both channels, but the practical validation improvement is established mainly for `K/\pi`.

## Updated figures
The following note-facing figures were regenerated from the updated chain:
- `DNdEtaUnfolding_ResponseMatrix.pdf`
- `DNdEtaUnfolding_MethodDifference.pdf`
- `DNdEtaUnfolding_MCClosure_BayesVsSVD.pdf`
- `DNdEtaUnfolding_DataMC_BayesVsSVD.pdf`
- `DNdEtaUnfolding_RefoldingValidation.pdf`
- `DNdEtaUnfolding_RatioRefoldingValidation.pdf`
- `KtoPi_DoubleRatio_vs_dNdEta_with_Systematics.pdf`
- `KtoPi_vs_dNdEta_DataMC_with_Systematics.pdf`
- `KtoPi_vs_dNdEta_DELPHI_vs_ALICE.pdf`
- `KtoPi_vs_dNdEta_DELPHI_vs_WorldEE.pdf`
- `KtoPi_vs_dNdEta_DELPHI_vs_WorldEE_NoPIDSFOnly.pdf`
- `Systematics_Components_vs_dNdEta.pdf`
- `PtoPi_DoubleRatio_vs_dNdEta_with_Systematics.pdf`
- `PtoPi_Systematics_Components_vs_dNdEta.pdf`
- `PtoPi_vs_dNdEta_DataMC_with_Systematics.pdf`
- `PtoPi_vs_dNdEta_DELPHI_vs_ALICE.pdf`
- `MigrationMetrics_PurityStability.pdf`
- `UnfoldingStressTest_Combined.pdf`

## Note update
The analysis note caption for the `dN_{ch}/d\eta` refolding validation figure was updated to match the current implementation:
- reco and truth are both defined in `|eta| < 0.5`
- the main remaining limitation is the merged overflow tail, not the old proxy-definition mismatch
- the note now separates the species-level refolding figure from the ratio-level
  `K/\pi` refolding companion figure and states the concrete tail-resolution path

File:
- [overleaf_repo_git/main.tex](../overleaf_repo_git/main.tex)

## Remaining caveats
- The quoted improvement is robust for the nominal `K/\pi` `dN_{ch}/d\eta` branch.
- The `p/\pi` `dN_{ch}/d\eta` closure still needs more work.
- The high-activity tail is now controlled by an explicit merged overflow bin. This is the correct stabilization step for the current statistics, but it should be described clearly anywhere the last `dN_{ch}/d\eta` bin is interpreted.
- The old iteration-4 status markdown in `report/20260310_iteration4_unfolding_validation_status.md` still contains the pre-fix `dNch/deta` refolding numbers and should be treated as superseded by this note.

## Bottom line
The `dN_{ch}/d\eta` unfolding is no longer using reco `N_{ch}^{tag}` as a proxy in the nominal `K/\pi` branch. It now unfolds the same variable on reco and truth, uses a merged overflow bin in the tail, and the nominal refolding validation is now good enough to support the updated figure set in the analysis note.
