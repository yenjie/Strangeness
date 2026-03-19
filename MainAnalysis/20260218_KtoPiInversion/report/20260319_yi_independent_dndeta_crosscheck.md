# Independent Yi-style dNch/deta Cross-check

Date: `2026-03-19`

## Goal

Build a genuinely more independent `K/pi` vs `dN_ch/deta` cross-check that is
closer to Yi Chen's slide procedure than the reduced factorized stress
tests.

Target ordering:

1. fake-correct observed tagged reco yields
2. unfold the full activity/kinematics variable `mu`
3. apply the `3x3` species inversion in true `mu`
4. apply the truth-side matched/gen efficiency correction

with

- `mu = (dN_ch/deta, p, |cos(theta)|)`
- active `PassAll` event selection
- active weak-decay-daughter veto / DELPHI-like truth policy
- active PID SF convention:
  - `R_after = R_before / (SF_K / SF_PI)`

## New source files

- [tools/build_yi_independent_dndeta_inputs.cc](/raid5/data/yjlee/strangeness/Strangeness/MainAnalysis/20260218_KtoPiInversion/tools/build_yi_independent_dndeta_inputs.cc)
- [run_yi_independent_dndeta.py](/raid5/data/yjlee/strangeness/Strangeness/MainAnalysis/20260218_KtoPiInversion/run_yi_independent_dndeta.py)

## Input builder

The builder reads the active v2.5 trees and writes flat reco/truth inputs:

- fake-corrected observed reco-tag counts:
  - `hRawTagKRecoFlat`
  - `hRawTagPiRecoFlat`
  - `hRawTagPRecoFlat`
- observed-tag flat responses:
  - `hRespTagKFlat`
  - `hRespTagPiFlat`
  - `hRespTagPFlat`
- truth denominators:
  - `hGenAllKFlat`, `hGenAllPiFlat`, `hGenAllPFlat`
- matched truth denominators:
  - `hMatchedDenKFlat`, `hMatchedDenPiFlat`, `hMatchedDenPFlat`
- truth-cell species-tag numerators:
  - `hTrueKTagAsKFlat`, ...

Flat indexing is

- `mu = (dN_ch/deta, p, |cos(theta)|)`

with the original fine binning

- `dN_ch/deta`: `16` uniform bins over `[-0.5, 60.5]`
- `p`: `0.4, 0.8, 1.2, 1.6, 2.0, 3.0, 5.0, 10.0, 50.0`
- `|cos(theta)|`: `0.15, 0.25, 0.35, 0.50, 0.675`

## Important bookkeeping fix

The first implementation was wrong in one crucial way:

- I applied the fake/matching weights not only to the measured observed-tag
  yields, but also to the matched response fills and matched truth-tag
  numerators.

That was inconsistent.

Correct treatment:

- measured observed-tag yields:
  - fake-correct with the reco-match weight
- matched response:
  - fill with weight `1`
- matched truth-tag numerators:
  - fill with weight `1`

Reason:

- once the reco-gen match is explicit, those objects are already in the matched
  space
- weighting them again by the reco-match probability double-counts the fake
  correction and drives a broad nonclosure

This fix is what made the independent branch close.

## Solver definition

The standalone solver does not call the nominal unfolding macros.

It does:

1. read the flat ROOT inputs
2. aggregate the fine `p` and `|cos(theta)|` bins to a stable coarse 3D grid
3. unfold each observed tag category in the full flat reco `mu` space with a
   local iterative Bayes implementation
4. apply a truth-cell `3x3` inversion
5. divide by the truth-cell matched/gen efficiency
6. integrate over `p` and `|cos(theta)|`
7. collapse the 16 activity bins to the active 8 visible `dN_ch/deta` bins for
   direct comparison to the nominal result
8. apply the same residual nonclosure correction and PID SF convention as the
   main note comparison

Current stable coarse grouping:

- `p` groups:
  - `0.4-1.2`
  - `1.2-2.0`
  - `2.0-5.0`
  - `5.0-50.0` GeV/c
- `|cos(theta)|` groups:
  - `0.15-0.35`
  - `0.35-0.675`

So the active independent grid is:

- `16 x 4 x 2`

## Runtime products

Builder outputs:

- [output/yi_independent_dndeta/mc_inputs.root](/raid5/data/yjlee/strangeness/Strangeness/MainAnalysis/20260218_KtoPiInversion/output/yi_independent_dndeta/mc_inputs.root)
- [output/yi_independent_dndeta/data_inputs.root](/raid5/data/yjlee/strangeness/Strangeness/MainAnalysis/20260218_KtoPiInversion/output/yi_independent_dndeta/data_inputs.root)

Independent result:

- [result/20260319/yi_independent_dndeta/yi_independent_dndeta.root](/raid5/data/yjlee/strangeness/Strangeness/MainAnalysis/20260218_KtoPiInversion/result/20260319/yi_independent_dndeta/yi_independent_dndeta.root)
- [result/20260319/yi_independent_dndeta/yi_independent_dndeta.txt](/raid5/data/yjlee/strangeness/Strangeness/MainAnalysis/20260218_KtoPiInversion/result/20260319/yi_independent_dndeta/yi_independent_dndeta.txt)
- [result/20260319/yi_independent_dndeta/YiIndependent_DNdEta_Comparison.pdf](/raid5/data/yjlee/strangeness/Strangeness/MainAnalysis/20260218_KtoPiInversion/result/20260319/yi_independent_dndeta/YiIndependent_DNdEta_Comparison.pdf)
- [result/20260319/yi_independent_dndeta/YiIndependent_DNdEta_ClosureComparison.pdf](/raid5/data/yjlee/strangeness/Strangeness/MainAnalysis/20260218_KtoPiInversion/result/20260319/yi_independent_dndeta/YiIndependent_DNdEta_ClosureComparison.pdf)

## Final numbers

From [yi_independent_dndeta.txt](/raid5/data/yjlee/strangeness/Strangeness/MainAnalysis/20260218_KtoPiInversion/result/20260319/yi_independent_dndeta/yi_independent_dndeta.txt):

- `nPBins = 4`
- `nAbsCosBins = 2`
- `singularCells = 48`
- `pseudoInverseCells = 4`
- `zeroGenEffCells = 147`
- `nominalClosureRMS = 0.022939`
- `yiClosureRMS = 0.036203`
- `maxAbsShift = 0.132909`
- `maxRelShift = 0.124474`
- `maxShiftOverNominalSys = 1.758720`
- `maxShiftOverNominalSysBin = 8`

Bin-by-bin comparison:

| Bin | Center | Nominal | Yi-style | Yi/Nominal | Closure | Shift / nominal sys |
|---|---:|---:|---:|---:|---:|---:|
| 1 | 1.406 | 1.1049 | 1.0825 | 0.9797 | 0.9957 | 0.295 |
| 2 | 5.219 | 1.1212 | 1.1063 | 0.9868 | 0.9796 | 0.187 |
| 3 | 9.031 | 1.0825 | 1.0717 | 0.9900 | 1.0063 | 0.145 |
| 4 | 12.844 | 1.0939 | 1.0895 | 0.9960 | 0.9849 | 0.057 |
| 5 | 16.656 | 1.0774 | 1.0811 | 1.0035 | 0.9876 | 0.050 |
| 6 | 20.469 | 1.0630 | 1.0539 | 0.9915 | 1.0014 | 0.123 |
| 7 | 24.281 | 1.0425 | 1.0431 | 1.0006 | 1.0073 | 0.008 |
| 8 | 28.094 | 1.0678 | 0.9349 | 0.8755 | 1.0978 | 1.759 |

## Interpretation

This branch is now technically coherent.

What changed after the bookkeeping fix:

- before the fix:
  - closure RMS was about `0.124`
  - broad low-bin bias was present
- after the fix:
  - closure RMS improved to `0.036`
  - the first 7 visible bins agree with the nominal result within about `2%`

So the earlier failure was not evidence that the Yi-style ordering is
intrinsically wrong. It was a weighting bug in the matched response / matched
truth-tag inputs.

Current conclusion:

- the independent Yi-style branch supports the nominal quoted result in the
  populated low- and mid-activity bins
- the merged final tail bin remains method-sensitive
- the current nominal branch should still be kept as the quoted result because:
  - it closes better overall
  - it is already the validated production chain
  - the only significant disagreement is confined to the merged tail bin

## Note update

This study was added to the appendix of the analysis note as a separate section:

- `Independent Yi-style dNch/deta 3D Cross-check`

with figures:

- `StrangenessV3/Figures/YiIndependent_DNdEta_Comparison.pdf`
- `StrangenessV3/Figures/YiIndependent_DNdEta_ClosureComparison.pdf`
