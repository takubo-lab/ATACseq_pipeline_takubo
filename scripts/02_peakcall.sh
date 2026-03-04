#!/usr/bin/env bash
# =============================================================================
#  Step 2: Peak Calling
#
#  処理内容:
#    1. 全サンプルの dedup BAM を入力として MACS3 callpeak (BAMPE, p=0.01)
#       → 全サンプル統合の単一 peak set を生成
#    2. サミットファイルからブラックリスト領域を除外
#    3. サミットを中心に ±SUMMIT_HALFWIDTH (デフォルト 125bp) に拡張
#       → Peak_250bp_noBL.bed (250bp 固定長ピーク)
#
#  入力: BAM/*.noDup.noMT.filt.sorted.bam
#  出力: Peak_nomodel/aggregated/       (MACS3 生出力)
#        Peak_nomodel/Summits_noBL.bed  (ブラックリスト除外済サミット)
#        Peak_nomodel/Peak_250bp_noBL.bed (250bp拡張ピーク, csaw入力)
# =============================================================================
set -euo pipefail
source "$(dirname "$0")/../config.sh"

# --- Substep control ---
# FROM_SUB: ""/"a"=MACS3+summit  "b"=summit only
FROM_SUB="${PIPELINE_FROM_SUBSTEP:-}"
FORCE="${PIPELINE_FORCE:-false}"

should_run_sub() {
  local sub="$1"
  [[ -z "${FROM_SUB}" || ! "${sub}" < "${FROM_SUB}" ]]
}

if [[ -n "${FROM_SUB}" ]]; then
  echo "  Substep start: ${FROM_SUB}"
fi

PEAK_DIR="${DIR}/${DIR_PEAKS}"
BAM_DIR="${DIR}/${DIR_BAM}"
AGG_DIR="${PEAK_DIR}/aggregated"

mkdir -p "${AGG_DIR}"

# --- dedup BAM ファイル収集 ---
shopt -s nullglob
BAM_FILES=( "${BAM_DIR}"/*.noDup.noMT.filt.sorted.bam )
shopt -u nullglob

if [[ ${#BAM_FILES[@]} -eq 0 ]]; then
  echo "ERROR: No dedup BAM files found in ${BAM_DIR}"
  exit 1
fi

echo "Found ${#BAM_FILES[@]} BAM file(s):"
printf '  %s\n' "${BAM_FILES[@]}"

# =============================================================================
#  [a] MACS3 peak calling — 全サンプル統合
# =============================================================================
if should_run_sub a; then
  SUMMIT_FILE="${AGG_DIR}/aggregated_narrow_summits.bed"
  if [[ "${FORCE}" != "true" && -f "${SUMMIT_FILE}" ]]; then
    echo "[a] Summit file already exists, skipping MACS3 (use --force to redo)"
  else
    echo "[a] Running MACS3 callpeak (format=${MACS3_FORMAT}, p=${MACS3_PVALUE})..."
    if [[ "${MACS3_FORMAT}" == "BAMPE" ]]; then
      macs3 callpeak -f BAMPE -t "${BAM_FILES[@]}" \
        -p "${MACS3_PVALUE}" -g "${GENOME_SIZE_MACS}" --keep-dup all \
        --outdir "${AGG_DIR}" -n aggregated_narrow \
        2>&1 | tee "${AGG_DIR}/macs3.log"
    else
      macs3 callpeak -f BAM -t "${BAM_FILES[@]}" \
        -p "${MACS3_PVALUE}" -g "${GENOME_SIZE_MACS}" --keep-dup all \
        --nomodel --shift -50 --extsize 100 \
        --outdir "${AGG_DIR}" -n aggregated_narrow \
        2>&1 | tee "${AGG_DIR}/macs3.log"
    fi
  fi
  echo "  Raw peaks called: $(wc -l < "${AGG_DIR}/aggregated_narrow_peaks.narrowPeak")"
fi

# =============================================================================
#  [b] Blacklist 除去 + summit → 250bp 拡張
# =============================================================================
if should_run_sub b; then
  SUMMIT_FILE="${AGG_DIR}/aggregated_narrow_summits.bed"
  if [[ ! -f "${SUMMIT_FILE}" ]]; then
    echo "ERROR: Summit file not found: ${SUMMIT_FILE}"
    echo "  Hint: Run substep a first:  ./run_pipeline.sh --from 2a"
    exit 1
  fi

  BL_NOCHR="${PEAK_DIR}/blacklist_noChrY.bed"
  grep -v "^chrY" "${BLACKLIST}" > "${BL_NOCHR}"

  echo "[b] Removing blacklist + extending to $((SUMMIT_HALFWIDTH * 2)) bp..."
  SUMMITS_NOBL="${PEAK_DIR}/Summits_noBL.bed"
  cat "${AGG_DIR}/aggregated_narrow_summits.bed" \
    | sort -k1,1V -k2,2n \
    | bedtools intersect -v -a stdin -b "${BL_NOCHR}" -sorted \
    > "${SUMMITS_NOBL}"
  echo "  Summits after blacklist removal: $(wc -l < "${SUMMITS_NOBL}")"

  PEAKS_250="${PEAK_DIR}/Peak_250bp_noBL.bed"
  awk -v hw="${SUMMIT_HALFWIDTH}" 'BEGIN{OFS="\t"}{
    s=$2-hw; e=$2+hw; if(s<0)s=0; print $1,s,e
  }' "${SUMMITS_NOBL}" | sort -k1,1V -k2,2n > "${PEAKS_250}"
  echo "  Fixed-size peaks ($((SUMMIT_HALFWIDTH*2)) bp): $(wc -l < "${PEAKS_250}")"
fi

echo ""
echo "Step 2 complete."
echo "  Summits (no blacklist) : ${SUMMITS_NOBL}"
echo "  Fixed-size peaks       : ${PEAKS_250}"
