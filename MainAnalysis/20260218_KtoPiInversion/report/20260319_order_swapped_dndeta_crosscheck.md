# 2026-03-19 order-swapped dN/deta cross-check

## Goal

Test the impact of swapping the correction order in the `K/pi` vs `dN_ch/deta`
branch.

Nominal chain:
1. fake / reco-match correction
2. 3x3 PID inversion
3. gen-matching correction
4. activity unfolding

Cross-check chain implemented here:
1. fake-correct observed K/pi/p-tagged reco yields
2. unfold each observed tag species in `dN_ch/deta`
3. apply a truth-matched 3x3 species inversion in true `dN_ch/deta x pT`
4. apply a truth-matched gen-efficiency correction
5. form `K/pi`, divide data by MC, and apply the same global PID SF convention as the nominal result

## Scope

This is an approximate reordered cross-check within the existing factorized
framework. It is **not** Yi Chen's full 3D `(p, |cos(theta)|, activity)`
procedure.

The current implementation keeps the active DELPHI-like definitions for:
- `PassAll` event selection
- charged-particle `dN_ch/deta(|eta|<0.5)` activity
- exclusive observed PID categories with the legacy `K > pi > p` tie rule
- `0.4 < pT < 5.0 GeV/c` and the nominal PID fiducial for the counted `K/pi`
- the current global PID scale-factor convention `R_after = R_before / (SF_K / SF_PI)`

The new ingredient is only the ordering of the activity unfolding relative to
the species untangling.

## Inputs produced

A dedicated standalone tool was added:
- `tools/build_order_swapped_dndeta_inputs.cc`

It loops the nominal MC or data tree and stores:
- fake-corrected observed tagged yields vs reco `dN_ch/deta x pT`
- tag-specific activity responses `true dN_ch/deta -> reco dN_ch/deta`
- truth-side matched denominators and tagged numerators for the 3x3 inversion
- truth-side generated and matched denominators for the gen-efficiency step

Products:
- `output/order_swapped_dndeta/mc_inputs.root`
- `output/order_swapped_dndeta/data_inputs.root`

The reordered builder is:
- `run_order_swapped_dndeta_crosscheck.py`

Main outputs:
- `output/order_swapped_dndeta/order_swapped_dndeta_crosscheck.root`
- `output/order_swapped_dndeta/order_swapped_dndeta_crosscheck.txt`
- `output/order_swapped_dndeta/OrderSwapped_DNdEta_Comparison.pdf`
- `output/order_swapped_dndeta/OrderSwapped_DNdEta_ClosureComparison.pdf`
- `output/order_swapped_dndeta/OrderSwapped_DNdEta_RatioRefoldingValidation.pdf`
- `output/order_swapped_dndeta/OrderSwapped_DNdEta_MethodShiftVsSystematics.pdf`

## Result summary

The reordered cross-check remains close to the current nominal result at low and
mid activity, but diverges in the merged high-activity tail.

Bin-by-bin `order-swapped / nominal`:
- bin 1 (`1.41`): `0.995`
- bin 2 (`5.22`): `0.989`
- bin 3 (`9.03`): `0.992`
- bin 4 (`12.84`): `0.976`
- bin 5 (`16.66`): `0.954`
- bin 6 (`20.47`): `0.923`
- bin 7 (`24.28`): `0.905`
- bin 8 (`28.09`): `0.836`

Summary metrics:
- maximum absolute shift: `0.1747`
- maximum relative shift: `16.36%`
- maximum `|shift| / nominal total systematic`: `2.312` in the last visible bin

Nominal vs reordered performance metrics:
- closure RMS(`unfolded / truth - 1`):
  - nominal: `0.0229`
  - reordered: `0.0691`
- reco-space refolding RMS(`refolded / reco - 1`):
  - nominal MC: `0.0090`
  - nominal data: `0.0101`
  - reordered MC: `0.0493`
  - reordered data: `0.0533`

The reordered closure itself is not perfectly flat:
- `0.986, 0.978, 0.994, 0.992, 1.024, 1.066, 1.104, 1.148`

So the same conclusion follows from the closure and from the data/MC comparison:
changing the order matters mainly in the final merged tail bins.

The extra validation views reinforce the same point:
- low and mid bins are not very sensitive to the ordering choice
- the reordered branch becomes clearly less stable in the tail
- the tail shift is larger than the quoted nominal total systematic from bin 6
  onward

## Interpretation

This supports keeping the current chain as the nominal result.

Reason:
- the current nominal chain is the validated production factorization already
  tied to the stored reco-side PID quantities and to the quoted residual and
  unfolding checks
- the reordered implementation is useful as a method cross-check, but its own
  closure and reco-space refolding performance are both visibly worse than the
  nominal branch
- where the reordered method shift becomes larger than the nominal systematic,
  it does so in the same sparse merged tail bins where its self-validation is
  also weakest

So the reordered result should be documented as a cross-check, not promoted to
replace the nominal quoted numbers.
