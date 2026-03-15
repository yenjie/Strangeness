#!/usr/bin/env python3
import math
import os
import shlex
import subprocess
from dataclasses import dataclass
from typing import Dict, List

import ROOT

ROOT.gROOT.SetBatch(True)
FORCE_RERUN = os.environ.get("FORCE_RERUN", "0") == "1"
ACTIVE_KEEPBINS = -1
KEEPBINS_OVERRIDE = int(os.environ.get("KEEPBINS_OVERRIDE", str(ACTIVE_KEEPBINS)))
OUT_DIR = os.environ.get("KPI_DNDY_DIR", "output/systematics_20260314_dndy")
PRECOMPUTED_DIR = os.environ.get("KPI_DNDY_INPUT_DIR", "output/systematics_20260314_dndy_inputs")
INPUT_DATA = os.environ.get("KPI_INPUT_DATA", "sample/Strangeness/merged_data_v2.5.root")
INPUT_MC = os.environ.get("KPI_INPUT_MC", "sample/Strangeness/merged_pythia_v2.5.root")
ANALYSIS_EXTRA_ARGS = shlex.split(os.environ.get("KPI_ANALYSIS_EXTRA_ARGS", ""))
FINALIZE_SCRIPT = os.environ.get("KPI_DNDY_FINALIZER", "finalize_dndy_systematics.py")


@dataclass
class Variation:
    name: str
    group: str
    args: Dict[str, str]


def run_cmd(cmd: List[str]) -> None:
    print("Running:", " ".join(cmd))
    subprocess.run(cmd, check=True)


def has_hist(path: str, name: str) -> bool:
    if not os.path.exists(path) or os.path.getsize(path) <= 0:
        return False
    try:
        f = ROOT.TFile.Open(path, "READ")
    except OSError:
        return False
    if not f or f.IsZombie():
        return False
    ok = f.Get(name) is not None
    f.Close()
    return ok


def run_analysis(input_file: str, output_file: str, args: Dict[str, str], require_hist: str) -> None:
    if (not FORCE_RERUN) and has_hist(output_file, require_hist):
        print("Skipping existing output:", output_file)
        return
    cmd = ["./ExecuteKtoPiAnalysis", "--Input", input_file, "--Output", output_file]
    for k, v in args.items():
        cmd.extend([f"--{k}", str(v)])
    cmd.extend(ANALYSIS_EXTRA_ARGS)
    run_cmd(cmd)


def run_unfold(
    mc_path: str,
    data_path: str,
    out_path: str,
    response_path: str = "",
    niter: int = 1,
    kreg: int = 8,
    keepbins_override: int = -1,
) -> None:
    if (not FORCE_RERUN and
        has_hist(out_path, "hDataOverMcBayes_dNdY") and
        has_hist(out_path, "hDataOverMcSVD_dNdY")):
        print("Skipping existing unfolding:", out_path)
        return
    macro_call = (
        f'runDNdYUnfolding_BayesSVD.C({niter},{kreg},"{mc_path}","{data_path}","{out_path}","{response_path}",false,{keepbins_override})'
    )
    run_cmd(["root", "-l", "-b", "-q", macro_call])


def main():
    out_dir = OUT_DIR
    os.makedirs(out_dir, exist_ok=True)
    os.makedirs(PRECOMPUTED_DIR, exist_ok=True)
    print(f"Active dN/dy status branch uses KEEPBINS_OVERRIDE={KEEPBINS_OVERRIDE}")
    variations = [
        Variation("nominal", "nominal", {"IsGen": "false", "UseMCTruthMatrix": "false", "UsePIDFiducial": "true", "PtMin": "0.4", "NtagPtMin": "0.2", "PtMax": "5.0", "MinThetaDeg": "30", "MaxThetaDeg": "150", "MinNch": "7", "NPtBins": "12"}),
        Variation("npt10", "binning", {"IsGen": "false", "UseMCTruthMatrix": "false", "UsePIDFiducial": "true", "PtMin": "0.4", "NtagPtMin": "0.2", "PtMax": "5.0", "MinThetaDeg": "30", "MaxThetaDeg": "150", "MinNch": "7", "NPtBins": "10"}),
        Variation("npt14", "binning", {"IsGen": "false", "UseMCTruthMatrix": "false", "UsePIDFiducial": "true", "PtMin": "0.4", "NtagPtMin": "0.2", "PtMax": "5.0", "MinThetaDeg": "30", "MaxThetaDeg": "150", "MinNch": "7", "NPtBins": "14"}),
    ]

    input_data = INPUT_DATA
    input_mc = INPUT_MC
    precomputed_dir = PRECOMPUTED_DIR

    unfold_roots: Dict[str, str] = {}
    for v in variations:
        data_out = os.path.join(precomputed_dir, f"{v.name}_data.root")
        mc_out = os.path.join(precomputed_dir, f"{v.name}_mc.root")
        run_analysis(input_data, data_out, v.args, "hKCorrectedDNdY")
        run_analysis(input_mc, mc_out, v.args, "hKCorrectedDNdY")

        uf_out = os.path.join(out_dir, f"{v.name}_unfold_dndy.root")
        # Use the variation-specific MC output as the response source so the
        # unfolding stays consistent with the thrust-axis dN/dy definition.
        run_unfold(mc_out, data_out, uf_out, "", keepbins_override=KEEPBINS_OVERRIDE)
        unfold_roots[v.name] = uf_out
    print("Delegating finalized dN/dy summary outputs to", FINALIZE_SCRIPT)
    run_cmd(["python3", FINALIZE_SCRIPT])
    print("dN/dy unfolding roots refreshed and finalized outputs rebuilt in", out_dir)


if __name__ == "__main__":
    main()
