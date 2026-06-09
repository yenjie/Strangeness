# Plot Final Results

This note documents the commands to reproduce the final K/pi plots from the current workflow.

## 1) Run full workflow (MC closure + nominal MC/data)

```bash
./analysis.sh
```

Outputs:
- `output/KtoPi-MC-Gen-Closure.root`
- `output/KtoPi-MC-Reco-Closure.root`
- `output/KtoPi-MC-Gen-Nominal.root`
- `output/KtoPi-MC-Reco-Nominal.root`
- `output/KtoPi-Data-Reco-Nominal.root`

## 2) Plot final K/pi comparisons

```bash
root -l -b -q plotClosure.C
root -l -b -q plotResultCor.C
root -l -b -q plotClosureVsData.C
```

Macro purpose:
- `plotClosure.C`: compares MC corrected reco `hKoverPiCorrected` from `output/KtoPi-MC-Reco-Closure.root` to generator-level `hKoverPi` from `output/KtoPi-MC-Gen-Closure.root`. The lower panel shows the closure ratio, corrected reco divided by generator truth.
- `plotResultCor.C`: overlays the nominal fully corrected data, fully corrected MC, and generator-level MC K/pi ratios. This is the main final-result plot for comparing DELPHI data with the corrected MC baseline and generator reference.
- `plotClosureVsData.C`: combines the MC closure result with the fully corrected data comparison in one figure, useful for checking whether the data/MC trend is larger than the residual MC non-closure.

Main figures:
- `KtoPiClosure_Overlay.pdf`
- `KtoPiRatio_Overlay.pdf`
- `KtoPiClosureVsData_Overlay.pdf`

## 3) Optional closure diagnostics by species

```bash
root -l -b -q plotClosureSpectra.C
root -l -b -q plotClosureByNtagSummary.C
```

Macro purpose:
- `plotClosureSpectra.C`: projects the two-dimensional species spectra over `Ntag_ch` and compares corrected reco to generator truth as a function of `pT` for kaons, pions, and protons.
- `plotClosureByNtagSummary.C`: builds species-specific closure maps in `(Ntag_ch, pT)` and integrated closure curves versus `Ntag_ch`. This is the main diagnostic for identifying where K, pi, or proton closure drives deviations in the final K/pi ratio.

Additional figures:
- `KPtClosure_Overlay.pdf`
- `PiPtClosure_Overlay.pdf`
- `PPtClosure_Overlay.pdf`
- `Closure_ByNtag_Summary.pdf`
