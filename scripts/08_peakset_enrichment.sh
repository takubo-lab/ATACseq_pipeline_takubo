#!/usr/bin/env bash
# =============================================================================
#  Step 8: Peak-set enrichment
#
#  Substeps:
#    8a. Normalize external BED region sets onto the canonical project peak set
#    8b. Run fgsea-like enrichment on edgeR DAR rankings
#
#  Inputs:
#    - Peak_nomodel/Peaks_250bp_th0.bed
#    - Peak_nomodel/edgeR/*_DA_edgeR_Results.tsv
#    - Differential_Peak_Iwama/*.bed
#    - Differential_Peak_Kobayashi/*.bed
#
#  Outputs:
#    - Peak_nomodel/PeakSetLibrary/region_sets_manifest.tsv
#    - Peak_nomodel/PeakSetLibrary/normalized/{source}/*.bed
#    - Peak_nomodel/PeakSetEnrichment/*
# =============================================================================
set -euo pipefail

_cfg="${PIPELINE_CONFIG_FILE:-$(dirname "$0")/../config.sh}"
source "${_cfg}"

FROM_SUB="${PIPELINE_FROM_SUBSTEP:-}"
FORCE="${PIPELINE_FORCE:-false}"
ONLY_STEPS="${PIPELINE_ONLY_STEPS:-}"

should_run_sub() {
  local sub="$1"
  [[ -z "${FROM_SUB}" || ! "${sub}" < "${FROM_SUB}" ]]
}

resolve_path() {
  local path_value="$1"
  if [[ -z "${path_value}" ]]; then
    return 0
  fi
  if [[ "${path_value}" = /* ]]; then
    printf '%s\n' "${path_value}"
  else
    printf '%s\n' "${DIR}/${path_value}"
  fi
}

discover_peak_sources() {
  local -a sources=()

  if declare -p FGSEA_SET_SOURCES &>/dev/null; then
    sources=("${FGSEA_SET_SOURCES[@]}")
  else
    [[ -d "${DIR}/Differential_Peak_Iwama" ]] && sources+=("${DIR}/Differential_Peak_Iwama")
    [[ -d "${DIR}/Differential_Peak_Kobayashi" ]] && sources+=("${DIR}/Differential_Peak_Kobayashi")
  fi

  for src in "${sources[@]}"; do
    resolve_path "${src}"
  done
}

parse_set_metadata() {
  local base_name="$1"
  local source_label="$2"
  local comparison_name="${base_name}"
  local direction="mixed"
  local fdr_label=""

  if [[ "${base_name}" =~ ^(.*)_(high|low)_FDR(.+)$ ]]; then
    comparison_name="${BASH_REMATCH[1]}"
    direction="${BASH_REMATCH[2]}"
    fdr_label="${BASH_REMATCH[3]}"
  fi

  printf '%s\t%s\t%s\t%s\n' "${source_label}" "${comparison_name}" "${direction}" "${fdr_label}"
}

peak_dir="${DIR}/${DIR_PEAKS}"
reference_bed="${FGSEA_REFERENCE_BED:-${peak_dir}/Peaks_250bp_th0.bed}"
reference_bed="$(resolve_path "${reference_bed}")"

if [[ ! -f "${reference_bed}" && -f "${peak_dir}/Peak_250bp_th0.bed" ]]; then
  reference_bed="${peak_dir}/Peak_250bp_th0.bed"
fi

if [[ ! -f "${reference_bed}" ]]; then
  echo "ERROR: canonical peak BED not found: ${reference_bed}"
  exit 1
fi

library_dir="${peak_dir}/${DIR_PEAKSET_LIBRARY}"
normalized_dir="${library_dir}/normalized"
manifest_file="${library_dir}/region_sets_manifest.tsv"
results_dir="${peak_dir}/${DIR_FGSEA}"

run_prepare="false"
run_enrichment="false"

if [[ "${ONLY_STEPS}" == "8a" ]]; then
  run_prepare="true"
elif [[ "${ONLY_STEPS}" == "8b" ]]; then
  run_enrichment="true"
else
  should_run_sub a && run_prepare="true"
  should_run_sub b && run_enrichment="true"
fi

mkdir -p "${library_dir}" "${normalized_dir}" "${results_dir}"

if [[ "${run_prepare}" == "true" ]]; then
  echo "  Preparing normalized peak-set library"
  printf 'source\tcomparison\tset_name\tdirection\tfdr_label\tinput_bed\tnormalized_bed\tpeak_count\n' > "${manifest_file}"

  while IFS= read -r src_dir; do
    [[ -n "${src_dir}" ]] || continue
    if [[ ! -d "${src_dir}" ]]; then
      echo "  WARN: region-set directory not found, skipping: ${src_dir}"
      continue
    fi

    source_label="$(basename "${src_dir}")"
    source_label="${source_label#Differential_Peak_}"
    source_out_dir="${normalized_dir}/${source_label}"
    mkdir -p "${source_out_dir}"

    shopt -s nullglob
    bed_files=("${src_dir}"/*.bed)
    shopt -u nullglob

    if [[ ${#bed_files[@]} -eq 0 ]]; then
      echo "  WARN: no BED files found in ${src_dir}"
      continue
    fi

    for input_bed in "${bed_files[@]}"; do
      base_name="$(basename "${input_bed}" .bed)"
      [[ "${base_name}" == *_intersect ]] && continue

      normalized_bed="${source_out_dir}/${base_name}.bed"
      if [[ "${FORCE}" == "true" || ! -s "${normalized_bed}" ]]; then
        bedtools intersect -a "${reference_bed}" -b "${input_bed}" -wa |
          sort -k1,1 -k2,2n -k3,3n |
          awk '!seen[$0]++' > "${normalized_bed}"
      fi

      peak_count=$(wc -l < "${normalized_bed}")
      IFS=$'\t' read -r parsed_source comparison_name direction fdr_label < <(parse_set_metadata "${base_name}" "${source_label}")
      printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "${parsed_source}" \
        "${comparison_name}" \
        "${base_name}" \
        "${direction}" \
        "${fdr_label}" \
        "${input_bed}" \
        "${normalized_bed}" \
        "${peak_count}" >> "${manifest_file}"
    done
  done < <(discover_peak_sources)

  echo "  Peak-set manifest: ${manifest_file}"
fi

if [[ "${run_enrichment}" == "true" ]]; then
  if [[ ! -s "${manifest_file}" ]]; then
    echo "ERROR: peak-set manifest not found. Run Step 8a first or use --steps 8."
    exit 1
  fi

  export PIPELINE_FGSEA_REFERENCE_BED="${reference_bed}"
  export PIPELINE_PEAK_SET_LIBRARY_DIR="${library_dir}"
  export PIPELINE_PEAK_SET_MANIFEST="${manifest_file}"
  export PIPELINE_FGSEA_RESULTS_DIR="${results_dir}"

  Rscript "${SCRIPT_DIR}/08_peakset_fgsea.R"
fi