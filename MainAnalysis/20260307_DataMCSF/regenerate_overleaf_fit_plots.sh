#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

tex_file="OverleafProject/DataMCSF/DataMCSF.tex"
figure_dir="OverleafProject/DataMCSF/Figures"

csv_value() {
  local key="$1"
  local file="$2"
  awk -F, -v k="$key" '$1 == k { print substr($0, length($1) + 2); exit }' "$file"
}

dataset_label_for() {
  local name="$1"
  case "$name" in
    *_mc_*|step1_phi_tagging_report.pdf|step4_a0025_kshort_tagging_report.pdf) echo "MC" ;;
    *_data_*) echo "Data" ;;
    step2_phi_random_pairs_report.pdf|step5_kshort_random_pairs_report.pdf) echo "MC" ;;
    step3_phi_random_pairs_report.pdf|step6_kshort_random_pairs_report.pdf) echo "Data" ;;
    *) echo "MC" ;;
  esac
}

rerun_phi_random_pairs() {
  local pdf="$1"
  local prefix="${pdf%_phi_random_pairs_report.pdf}"
  local summary="Reports/${prefix}_fit_summary.txt"
  local input seed sig bkg mass fit mom bins pmin pmax dataset
  input="$(csv_value Input "$summary")"
  seed="$(csv_value Step1Seed "$summary")"
  sig="$(csv_value SignalOverride "$summary")"
  bkg="$(csv_value BackgroundMode "$summary")"
  mass="$(csv_value MassRange "$summary")"
  fit="$(csv_value FitRange "$summary")"
  bins="$(csv_value NBins "$summary")"
  mom="$(csv_value MomentumBin "$summary")"
  IFS=, read -r mass_min mass_max <<< "$mass"
  IFS=, read -r fit_min fit_max <<< "$fit"
  if [[ "$mom" == inclusive* ]]; then
    pmin=-1.0
    pmax=-1.0
  else
    IFS=, read -r _mom_mode pmin pmax <<< "$mom"
  fi
  dataset="$(dataset_label_for "$pdf")"
  root -l -b -q "Code/step2_mc_random_pairs.C(\"$input\",\"$seed\",\"Reports\",\"$prefix\",\"$dataset\",\"$sig\",\"$bkg\",$mass_min,$mass_max,$fit_min,$fit_max,$bins,$pmin,$pmax)" </dev/null
}

rerun_kshort_random_pairs() {
  local pdf="$1"
  local prefix="${pdf%_kshort_random_pairs_report.pdf}"
  local summary="Reports/${prefix}_fit_summary.txt"
  local input seed sig bkg mass fit mom bins pmin pmax dataset
  input="$(csv_value Input "$summary")"
  seed="$(csv_value Step4Seed "$summary")"
  sig="$(csv_value SignalOverride "$summary")"
  bkg="$(csv_value BackgroundMode "$summary")"
  mass="$(csv_value MassRange "$summary")"
  fit="$(csv_value FitRange "$summary")"
  bins="$(csv_value NBins "$summary")"
  mom="$(csv_value MomentumBin "$summary")"
  IFS=, read -r mass_min mass_max <<< "$mass"
  IFS=, read -r fit_min fit_max <<< "$fit"
  if [[ "$mom" == inclusive* ]]; then
    pmin=-1.0
    pmax=-1.0
  else
    IFS=, read -r _mom_mode pmin pmax <<< "$mom"
  fi
  dataset="$(dataset_label_for "$pdf")"
  root -l -b -q "Code/step5_mc_random_pairs_kshort.C(\"$input\",\"$seed\",\"Reports\",\"$prefix\",\"$dataset\",\"$sig\",\"$bkg\",$mass_min,$mass_max,$fit_min,$fit_max,$bins,$pmin,$pmax)" </dev/null
}

copy_figure() {
  local pdf="$1"
  cp "Reports/$pdf" "$figure_dir/$pdf"
}

root -l -b -q 'Code/step1_make_report.C("Samples/merged_mc_v2.2.root","Reports",0.01,"step1")' </dev/null
root -l -b -q 'Code/step4_kshort_make_report.C("Samples/merged_mc_v2.2.root","Reports",0.025,"step4_a0025",0.42,0.58,0.40,0.60)' </dev/null

figure_list="$(mktemp)"
perl -ne 'while(/\\includegraphics(?:\\[[^\\]]*\\])?\\{([^}]+)\\}/g){print "$1\n"}' "$tex_file" |
  grep -E '(_report\.pdf|_tagging_report\.pdf)$' |
  sort -u > "$figure_list"

while IFS= read -r pdf; do
  case "$pdf" in
    step1_phi_tagging_report.pdf|step4_a0025_kshort_tagging_report.pdf)
      copy_figure "$pdf"
      ;;
    *_phi_random_pairs_report.pdf)
      rerun_phi_random_pairs "$pdf"
      copy_figure "$pdf"
      ;;
    *_kshort_random_pairs_report.pdf)
      rerun_kshort_random_pairs "$pdf"
      copy_figure "$pdf"
      ;;
  esac
done < "$figure_list"

rm -f "$figure_list"
