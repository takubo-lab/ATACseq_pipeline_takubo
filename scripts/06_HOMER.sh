#!/usr/bin/env bash
# =============================================================================
#  Step 6: HOMER Motif Analysis
#
#  処理内容:
#    - edgeR 結果 TSV を読み込み、有意な DAR を BED 形式で抽出
#      (FDR < HOMER_FDR かつ |logFC| > HOMER_LFC)
#    - up (logFC > 0) と down (logFC < 0) を別々に HOMER へ投入
#    - findMotifsGenome.pl でモチーフ検索 (背景 = Peaks_250bp_th0.bed)
#
#  環境変数 / config.sh:
#    HOMER_FDR, HOMER_LFC, HOMER_SIZE, HOMER_GENOME, HOMER_THREADS
#
#  入力: Peak_nomodel/edgeR/*_DA_edgeR_Results.tsv
#  出力: Peak_nomodel/HOMER/
#    {G1}_vs_{G2}_up_FDR{FDR}/     (HOMER モチーフ結果)
#    {G1}_vs_{G2}_down_FDR{FDR}/
#    {G1}_vs_{G2}_up_FDR{FDR}.bed   (入力 BED)
#    {G1}_vs_{G2}_down_FDR{FDR}.bed
# =============================================================================
set -euo pipefail
_cfg="${PIPELINE_CONFIG_FILE:-$(dirname "$0")/../config.sh}"
source "${_cfg}"
FORCE="${PIPELINE_FORCE:-false}"

PEAK_DIR="${DIR}/${DIR_PEAKS}"
EDGE_DIR="${PEAK_DIR}/edgeR"
HOMER_DIR="${PEAK_DIR}/HOMER"
BG_BED="${PEAK_DIR}/Peaks_250bp_th0.bed"

mkdir -p "${HOMER_DIR}"

# 背景 BED の確認
if [[ ! -f "${BG_BED}" ]]; then
  echo "ERROR: Background peak file not found: ${BG_BED}"
  exit 1
fi

# HOMER の確認
if ! command -v findMotifsGenome.pl &>/dev/null; then
  echo "ERROR: findMotifsGenome.pl not found. Please install HOMER and add to PATH."
  echo "       See: http://homer.ucsd.edu/homer/introduction/install.html"
  exit 1
fi

# HOMER ゲノムのインストール確認
_homer_bin="$(dirname "$(command -v findMotifsGenome.pl)")"
_homer_genome_dir="${_homer_bin}/../data/genomes/${HOMER_GENOME}"
if [[ ! -d "${_homer_genome_dir}" ]]; then
  echo "ERROR: HOMER genome '${HOMER_GENOME}' is not installed."
  echo "       Install with: perl ${_homer_bin}/configureHomer.pl install ${HOMER_GENOME}"
  exit 1
fi

# edgeR 結果 TSV 一覧
shopt -s nullglob
tsv_files=( "${EDGE_DIR}"/*_DA_edgeR_Results.tsv )
shopt -u nullglob

if [[ ${#tsv_files[@]} -eq 0 ]]; then
  echo "No edgeR result TSV files found in ${EDGE_DIR}"
  exit 1
fi

echo "Found ${#tsv_files[@]} edgeR result(s)"
echo "  FDR threshold  : ${HOMER_FDR}"
echo "  |logFC| cutoff : ${HOMER_LFC}"
echo "  HOMER size     : ${HOMER_SIZE}"
echo "  Genome         : ${HOMER_GENOME}"

# =============================================================================
#  ループ: 各ペアワイズ比較について
# =============================================================================
for tsv in "${tsv_files[@]}"; do
  base=$(basename "${tsv}" _DA_edgeR_Results.tsv)   # e.g. G1_vs_G2
  echo ""
  echo "--- Processing: ${base} ---"

  # ---- up: G1 > G2 (logFC > HOMER_LFC, FDR < HOMER_FDR) ----
  UP_BED="${HOMER_DIR}/${base}_up_FDR${HOMER_FDR}.bed"
  if [[ "${FORCE}" == "true" || ! -s "${UP_BED}" ]]; then
    awk -v fdr="${HOMER_FDR}" -v lfc="${HOMER_LFC}" \
      'BEGIN{OFS="\t"; FS="\t"}
       NR==1{next}
       $4+0 > lfc+0 && $8+0 < fdr+0 {print $1,$2,$3}' \
      "${tsv}" \
      | sort -k1,1V -k2,2n \
      > "${UP_BED}"
  fi
  N_UP=$(wc -l < "${UP_BED}")
  echo "  up   DAR: ${N_UP} peaks"

  # ---- down: G1 < G2 (logFC < -HOMER_LFC, FDR < HOMER_FDR) ----
  DOWN_BED="${HOMER_DIR}/${base}_down_FDR${HOMER_FDR}.bed"
  if [[ "${FORCE}" == "true" || ! -s "${DOWN_BED}" ]]; then
    awk -v fdr="${HOMER_FDR}" -v lfc="${HOMER_LFC}" \
      'BEGIN{OFS="\t"; FS="\t"}
       NR==1{next}
       $4+0 < -lfc+0 && $8+0 < fdr+0 {print $1,$2,$3}' \
      "${tsv}" \
      | sort -k1,1V -k2,2n \
      > "${DOWN_BED}"
  fi
  N_DOWN=$(wc -l < "${DOWN_BED}")
  echo "  down DAR: ${N_DOWN} peaks"

  # ----- HOMER: up -----
  UP_OUT="${HOMER_DIR}/${base}_up_FDR${HOMER_FDR}"
  if [[ "${N_UP}" -ge 10 && ( "${FORCE}" == "true" || ! -f "${UP_OUT}/knownResults.html" ) ]]; then
    echo "  Running HOMER (up)..."
    if [[ "${FORCE}" == "true" ]]; then
      rm -rf "${UP_OUT}"
    fi
    findMotifsGenome.pl \
      "${UP_BED}" \
      "${HOMER_GENOME}" \
      "${UP_OUT}/" \
      -size "${HOMER_SIZE}" \
      -bg "${BG_BED}" \
      -p "${HOMER_THREADS:-${THREADS}}" \
      -mask \
      2>&1 | tail -5
  elif [[ "${N_UP}" -lt 10 ]]; then
    echo "  Skipping HOMER (up): too few peaks (${N_UP} < 10)"
  else
    echo "  HOMER (up) result already exists, skipping."
  fi

  # ----- HOMER: down -----
  DOWN_OUT="${HOMER_DIR}/${base}_down_FDR${HOMER_FDR}"
  if [[ "${N_DOWN}" -ge 10 && ( "${FORCE}" == "true" || ! -f "${DOWN_OUT}/knownResults.html" ) ]]; then
    echo "  Running HOMER (down)..."
    if [[ "${FORCE}" == "true" ]]; then
      rm -rf "${DOWN_OUT}"
    fi
    findMotifsGenome.pl \
      "${DOWN_BED}" \
      "${HOMER_GENOME}" \
      "${DOWN_OUT}/" \
      -size "${HOMER_SIZE}" \
      -bg "${BG_BED}" \
      -p "${HOMER_THREADS:-${THREADS}}" \
      -mask \
      2>&1 | tail -5
  elif [[ "${N_DOWN}" -lt 10 ]]; then
    echo "  Skipping HOMER (down): too few peaks (${N_DOWN} < 10)"
  else
    echo "  HOMER (down) result already exists, skipping."
  fi
done

echo ""
echo "Step 6 complete. HOMER results: ${HOMER_DIR}"
