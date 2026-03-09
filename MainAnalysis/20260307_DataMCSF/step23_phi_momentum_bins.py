#!/usr/bin/env python3
import csv
import math
import pathlib
import re
import subprocess

ROOT = pathlib.Path(__file__).resolve().parents[1]
REPORTS = ROOT / "Reports"


def run_root(cmd: str):
  p = subprocess.run(cmd, shell=True, cwd=ROOT, capture_output=True, text=True)
  if p.returncode != 0:
    raise RuntimeError(f"Command failed:\n{cmd}\nSTDOUT:\n{p.stdout}\nSTDERR:\n{p.stderr}")


def read_eff(summary_path: pathlib.Path):
  s = summary_path.read_text()
  m = re.search(r"tagging_efficiency,([0-9.eE+-]+),tagging_efficiency_err,([0-9.eE+-]+)", s)
  if not m:
    raise RuntimeError(f"Could not parse efficiency in {summary_path}")
  return float(m.group(1)), float(m.group(2))


def read_yields(summary_path: pathlib.Path):
  s = summary_path.read_text()
  out = {}
  for c in ["0tag", "1tag", "2tag"]:
    m = re.search(rf"{c},model,([^,]+),chi2ndf,([0-9.eE+-]+),mean,([0-9.eE+-]+),meanErr,([0-9.eE+-]+),yield,([0-9.eE+-]+),yieldErr,([0-9.eE+-]+)", s)
    if not m:
      raise RuntimeError(f"Could not parse yields for {c} in {summary_path}")
    out[c] = float(m.group(5))
  return out


def run_one_bin(pmin: float, pmax: float, label: str):
  mc_prefix = f"step2_phi_pbin_{label}_mc"
  data_prefix = f"step3_phi_pbin_{label}_data"

  cmd_mc = (
      "root -l -b -q 'Code/step2_mc_random_pairs.C("
      f"\"Samples/merged_mc_v2.2.root\",\"Reports/step1_best_fit_params.csv\",\"Reports\",\"{mc_prefix}\",\"MC\","
      f"\"auto\",\"baseline\",0.99,1.06,1.00,1.05,280,{pmin:.3f},{pmax:.3f})'"
  )
  cmd_da = (
      "root -l -b -q 'Code/step2_mc_random_pairs.C("
      f"\"Samples/merged_data_v2.2.root\",\"Reports/step1_best_fit_params.csv\",\"Reports\",\"{data_prefix}\",\"Data\","
      f"\"auto\",\"baseline\",0.99,1.06,1.00,1.05,280,{pmin:.3f},{pmax:.3f})'"
  )
  run_root(cmd_mc)
  run_root(cmd_da)

  mc_sum = REPORTS / f"{mc_prefix}_fit_summary.txt"
  da_sum = REPORTS / f"{data_prefix}_fit_summary.txt"
  emc, semc = read_eff(mc_sum)
  eda, seda = read_eff(da_sum)
  ymc = read_yields(mc_sum)
  yda = read_yields(da_sum)

  sf = eda / emc
  sf_err = sf * math.sqrt((seda / eda) ** 2 + (semc / emc) ** 2)
  return {
      "bin_label": label,
      "p_min": pmin,
      "p_max": pmax,
      "eff_mc": emc,
      "eff_mc_err": semc,
      "eff_data": eda,
      "eff_data_err": seda,
      "sf": sf,
      "sf_err": sf_err,
      "mc_y0": ymc["0tag"],
      "mc_y1": ymc["1tag"],
      "mc_y2": ymc["2tag"],
      "data_y0": yda["0tag"],
      "data_y1": yda["1tag"],
      "data_y2": yda["2tag"],
  }


def main():
  bins = [
      (0.0, 1.2, "p00_12"),
      (1.2, 2.2, "p12_22"),
      (2.2, 5.0, "p22_50"),
  ]

  rows = []
  for pmin, pmax, label in bins:
    rows.append(run_one_bin(pmin, pmax, label))

  out_csv = REPORTS / "step23_phi_momentum_binned_sf.csv"
  with out_csv.open("w", newline="") as f:
    w = csv.DictWriter(
        f,
        fieldnames=[
            "bin_label", "p_min", "p_max",
            "mc_y0", "mc_y1", "mc_y2",
            "data_y0", "data_y1", "data_y2",
            "eff_mc", "eff_mc_err", "eff_data", "eff_data_err", "sf", "sf_err",
        ],
    )
    w.writeheader()
    w.writerows(rows)

  out_txt = REPORTS / "step23_phi_momentum_binned_sf_summary.txt"
  with out_txt.open("w") as f:
    f.write("Phi->KK momentum-binned tagging SF (strict both-daughter p bin)\n")
    f.write("Seed: Reports/step1_best_fit_params.csv; model=auto; bkg=baseline\n\n")
    for r in rows:
      f.write(
          f"{r['bin_label']} [{r['p_min']:.2f},{r['p_max']:.2f}) GeV: "
          f"eff_mc={r['eff_mc']:.6f}+/-{r['eff_mc_err']:.6f}, "
          f"eff_data={r['eff_data']:.6f}+/-{r['eff_data_err']:.6f}, "
          f"SF={r['sf']:.6f}+/-{r['sf_err']:.6f}\n"
      )

  print(f"Wrote {out_csv}")
  print(f"Wrote {out_txt}")


if __name__ == "__main__":
  main()
