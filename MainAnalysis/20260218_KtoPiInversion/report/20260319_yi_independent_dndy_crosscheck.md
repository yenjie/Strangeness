# 2026-03-19 independent Yi-style dN/dy cross-check

## Goal
Build an independent thrust-axis `dN_ch/dy` branch that follows the Yi-style ordering more closely than the nominal DELPHI chain:

1. fake-correct observed reco tags in full flat `mu = (dN_ch/dy, p, |cos(theta)|)`
2. unfold observed tags in full flat `mu`
3. apply truth-cell `3x3` tag inversion
4. apply truth-cell matched/gen efficiency correction
5. integrate back to `K/pi` vs thrust-axis `dN_ch/dy`

This is intended as an appendix-level cross-check of the nominal thrust-axis branch, not a replacement nominal result.

## Inputs
- MC builder output:
  - `output/yi_independent_dndy/mc_inputs.root`
- Data builder output:
  - `output/yi_independent_dndy/data_inputs.root`
- Nominal summary:
  - `output/systematics_20260314_dndy/systematics_dndy_summary.root`
- Nominal unfolding package:
  - `output/systematics_20260314_dndy/nominal_unfold_dndy.root`

## Implementation
### Builder
- source:
  - `tools/build_yi_independent_dndy_inputs.cc`
- booked fine activity axis:
  - `[-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 12.5, 15.5, 20.5, 25.5, 30.5]`
- fine kinematic binning:
  - `p`: `0.4-0.8, 0.8-1.2, 1.2-1.6, 1.6-2.0, 2.0-3.0, 3.0-5.0, 5.0-10.0, 10.0-50.0 GeV/c`
  - `|cos(theta)|`: `0.15-0.25, 0.25-0.35, 0.35-0.50, 0.50-0.675`
- selection follows the active DELPHI truth policy:
  - `PassAll`
  - active PID fiducial
  - `0.4 <= pT < 5.0 GeV/c`
  - weak-decay-daughter veto through the shared truth-counting policy

### Solver
- source:
  - `run_yi_independent_dndy.py`
- stable coarse grouping used in the solver:
  - `p`: `(0.4-1.2), (1.2-2.0), (2.0-5.0), (5.0-50.0) GeV/c`
  - `|cos(theta)|`: `(0.15-0.35), (0.35-0.675)`
- important bookkeeping rule:
  - fake/matching weights are applied only to the measured observed-tag yields
  - once reco-gen matching is explicit, the matched response and matched truth-tag numerators are unweighted
- final visible axis matches the nominal quoted `dN_ch/dy` result exactly:
  - 14 visible bins with tail merge `15.5 < dN_ch/dy < 30.5`
- same global PID scale-factor convention as the nominal result:
  - divide by `SF_K / SF_pi`

## Outputs
- ROOT summary:
  - `result/20260319/yi_independent_dndy/yi_independent_dndy.root`
- text summary:
  - `result/20260319/yi_independent_dndy/yi_independent_dndy.txt`
- comparison plot:
  - `result/20260319/yi_independent_dndy/YiIndependent_DNdY_Comparison.pdf`
- closure plot:
  - `result/20260319/yi_independent_dndy/YiIndependent_DNdY_ClosureComparison.pdf`

## Results
### Structural summary
- `nActivityBinsFine = 16`
- `nPBins = 4`
- `nAbsCosBins = 2`
- `singularCells = 3`
- `pseudoInverseCells = 6`
- `zeroGenEffCells = 16`

### Performance
- nominal closure RMS: `0.027504`
- Yi-style closure RMS: `0.014361`

So this independent thrust-axis branch closes better than the current nominal branch.

### Shift relative to nominal
- maximum absolute shift: `0.065394`
- maximum relative shift: `0.059415`
- maximum shift / nominal systematic: `0.673492`
- worst bin: `14` (`15.5 < dN_ch/dy < 30.5`)

So the Yi-style thrust-axis result stays inside the quoted nominal systematic envelope in every visible bin.

### Bin-by-bin comparison
- bins 1-12 agree with nominal within about `2.4%`
- bin 13 is lower by about `4.2%`
- bin 14 is lower by about `5.9%`
- even in the tail, the branch remains within the nominal systematic envelope

### Closure by bin
- the Yi-style closure stays close to unity in every visible bin
- bin 14 closure is `0.989228`

## Interpretation
This is materially different from the beam-axis `dN_ch/deta` case.

For thrust-axis `dN_ch/dy`, the more independent Yi-style branch:
- closes well
- is at least as stable as the nominal branch, and numerically better in closure RMS
- remains inside the quoted nominal systematic in all visible bins

So the thrust-axis result appears less sensitive to the correction ordering than the beam-axis `dN_ch/deta` branch.

## Conclusion
- keep the current validated nominal thrust-axis branch as the quoted result for consistency with the rest of the analysis
- retain the independent Yi-style `dN_ch/dy` branch as a strong appendix cross-check
- unlike the beam-axis order-swap and independent `dN_ch/deta` studies, the thrust-axis Yi-style branch does not indicate a tail-only tension beyond the quoted nominal systematic
