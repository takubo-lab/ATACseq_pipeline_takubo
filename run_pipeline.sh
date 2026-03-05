#!/usr/bin/env bash
# =============================================================================
#  ATAC-seq Pipeline — Master Orchestrator
#
#  使い方:
#    ./run_pipeline.sh                   # 全ステップ実行 (1a → 7)
#    ./run_pipeline.sh --from 3          # step 3 以降を実行
#    ./run_pipeline.sh --from 4b         # step 4 の substep b 以降を実行
#    ./run_pipeline.sh --to 4            # step 1–4 を実行
#    ./run_pipeline.sh --steps 2b,4c     # 指定 substep のみ実行
#    ./run_pipeline.sh --force --from 2b # 既存出力を無視して再実行
#    ./run_pipeline.sh --config /path/to/project/config.sh  # 外部プロジェクト
#    ./run_pipeline.sh --version         # バージョン表示
#
#  ── Substep 一覧 ─────────────────────────────────────────
#  1a: trim_galore アダプタートリミング
#  1b: bowtie2 アライメント + samtools フィルタ
#  1c: picard 重複除去
#  2a: MACS3 callpeak (全サンプル統合)
#  2b: サミット → blacklist 除去 → 250bp 拡張
#  3a: csaw カウント + histogram 出力 (閾値未設定なら一時停止)
#  3b: 閾値フィルタ → BED / カウント行列出力
#  4a: グループごと BAM マージ + blacklist 除去
#  4b: スケールファクター算出 (R)
#  4c: bamCoverage → bigWig 生成
#  5 : DAR 検出 (edgeR QL)
#  6 : HOMER モチーフ解析
#  7 : PCA・可視化
#
#  ── よくある再実行例 ────────────────────────────────────
#  SUMMIT_HALFWIDTH を変えた  → ./run_pipeline.sh --from 2b
#  BINSIZE を変えた           → ./run_pipeline.sh --from 4c
#  閾値確認後に再実行         → ./run_pipeline.sh --from 3b
#  DAR/HOMER/PCA だけやり直し → ./run_pipeline.sh --steps 5,6,7
# =============================================================================
set -euo pipefail
IFS=$'\n\t'

# =============================================================================
#  引数パース (設定読み込みより前に実行)
# =============================================================================
PIPELINE_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

FROM_STEP=1
FROM_SUB=""
TO_STEP=7
ONLY_STEPS=""          # --steps 用 (カンマ区切り, 例: "2b,5,6")
FORCE="false"
CONFIG_FILE=""         # --config 用 (外部プロジェクトの config.sh)
declare -A STEP_FROM_SUB   # --steps 用: step番号 → 開始 substep

# stage spec "4b" → PARSED_STEP=4, PARSED_SUB=b
parse_stage() {
  local spec="$1"
  PARSED_STEP="${spec//[a-z]/}"
  PARSED_SUB="${spec//[0-9]/}"
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --from)
      parse_stage "$2"
      FROM_STEP="${PARSED_STEP}"
      FROM_SUB="${PARSED_SUB}"
      shift 2 ;;
    --to)
      parse_stage "$2"
      TO_STEP="${PARSED_STEP}"
      shift 2 ;;
    --steps)
      ONLY_STEPS="$2"
      IFS=',' read -ra _specs <<< "$2"
      for _s in "${_specs[@]}"; do
        parse_stage "$_s"
        STEP_FROM_SUB["${PARSED_STEP}"]="${PARSED_SUB}"
      done
      shift 2 ;;
    --config)
      CONFIG_FILE="$2"
      shift 2 ;;
    --force)
      FORCE="true"
      shift ;;
    --version)
      source "${PIPELINE_ROOT}/scripts/utils.sh"
      echo "ATACseq_pipeline_takubo v$(get_pipeline_version) ($(get_git_commit))"
      exit 0 ;;
    -h|--help)
      grep '^#' "$0" | head -50 | sed 's/^# \{0,2\}//'
      exit 0 ;;
    *)
      echo "Unknown argument: $1  (try --help)"
      exit 1 ;;
  esac
done

# --- ユーティリティ読み込み ---
source "${PIPELINE_ROOT}/scripts/utils.sh"

# --- Conda 環境チェック ---
EXPECTED_ENV="atacseq_takubo"
current_env="${CONDA_DEFAULT_ENV:-}"
if [[ "${current_env}" != "${EXPECTED_ENV}" ]]; then
  echo -e "\033[1;33m[WARN]\033[0m Conda環境 '${EXPECTED_ENV}' がアクティブではありません。"
  echo -e "       現在の環境: '${current_env:-none}'"
  echo ""
  echo "  セットアップ:   bash setup_env.sh"
  echo "  アクティベート: conda activate ${EXPECTED_ENV}"
  echo ""
  read -rp "このまま続行しますか？ [y/N]: " yn
  case "${yn}" in
    [Yy]*) echo "[INFO] 現在の環境で続行します。" ;;
    *)     echo "[INFO] 中断しました。"; exit 1 ;;
  esac
fi

# --- 設定読み込み ---
if [[ -n "$CONFIG_FILE" ]]; then
  if [[ ! -f "$CONFIG_FILE" ]]; then
    echo "ERROR: Config file not found: ${CONFIG_FILE}"
    exit 1
  fi
  source "$CONFIG_FILE"
else
  source "${PIPELINE_ROOT}/config.sh"
fi

# --- ロック取得 ---
acquire_lock "${DIR}"

# --- ログ設定 ---
mkdir -p "${DIR}/logs"
LOG="${DIR}/logs/pipeline_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG") 2>&1

# --- バナー表示 ---
echo "=============================================="
echo " ATAC-seq Pipeline"
print_version_banner
echo " Date   : $(date)"
echo " User   : $(whoami)@$(hostname)"
echo " Project : ${DIR}"
echo " Log    : ${LOG}"
echo "  --from ${FROM_STEP}${FROM_SUB}  --to ${TO_STEP}  --force=${FORCE}"
[[ -n "$ONLY_STEPS" ]] && echo "  --steps ${ONLY_STEPS}"
[[ -n "$CONFIG_FILE" ]] && echo "  --config ${CONFIG_FILE}"
echo "=============================================="

# --- Provenance 記録 ---
generate_provenance "${DIR}"

# ── step を実行するか判断 ──────────────────────────────────
should_run() {
  local step="$1"
  if [[ -n "$ONLY_STEPS" ]]; then
    echo "$ONLY_STEPS" | tr ',' '\n' | grep -qE "^${step}[a-z]?$"
  else
    [[ $step -ge $FROM_STEP && $step -le $TO_STEP ]]
  fi
}

# ── step に対応する「開始 substep」を返す ─────────────────
get_from_substep() {
  local step="$1"
  if [[ -n "$ONLY_STEPS" ]]; then
    echo "${STEP_FROM_SUB[$step]:-}"
  elif [[ $step -eq $FROM_STEP ]]; then
    echo "$FROM_SUB"
  else
    echo ""
  fi
}

# ── step ラッパー ─────────────────────────────────────────
run_step() {
  local step_num="$1" step_name="$2"
  shift 2
  echo ""
  echo "====== Step ${step_num}: ${step_name} ======"
  echo "  Started : $(date)"
  export PIPELINE_FROM_SUBSTEP="$(get_from_substep "${step_num}")"
  export PIPELINE_FORCE="${FORCE}"
  "$@"
  echo "  Finished: $(date)"
}

# --- サンプルシートの存在確認 ---
if [[ ! -f "${SAMPLES_TSV}" ]]; then
  echo "ERROR: samples.tsv not found at ${SAMPLES_TSV}"
  echo "  → samples.tsv を編集してサンプル名とグループを定義してください。"
  exit 1
fi

# --- 必須ツール確認 ---
check_tool() {
  command -v "$1" &>/dev/null || { echo "ERROR: '$1' not found. Please install or check PATH."; exit 1; }
}

echo "--- Checking required tools ---"
for tool in trim_galore bowtie2 samtools bedtools; do
  check_tool "$tool"
done
java -jar "${PICARD_JAR}" --version &>/dev/null || { echo "ERROR: picard.jar not found at ${PICARD_JAR}"; exit 1; }
echo "  Required tools: OK"

# --- 環境アクティベート ---
if [[ -n "${MACS3_ENV}" && -f "${MACS3_ENV}" ]]; then
  # shellcheck source=/dev/null
  source "${MACS3_ENV}"
fi

# --- 設定変数を R スクリプトに渡す環境変数として export ---
export PIPELINE_DIR="${DIR}"
export PIPELINE_SAMPLES="${SAMPLES_TSV}"
export PIPELINE_BLACKLIST="${BLACKLIST}"
export PIPELINE_THREADS="${THREADS}"
export PIPELINE_GENOME="${GENOME}"
export PIPELINE_HOMER_GENOME="${HOMER_GENOME}"
export PIPELINE_PEAK_DIR="${DIR}/${DIR_PEAKS}"
export PIPELINE_BAM_DIR="${DIR}/${DIR_BAM}"
export PIPELINE_MERGED_BAM_DIR="${DIR}/${DIR_MERGEDBAM}"
export PIPELINE_BW_DIR="${DIR}/${DIR_BW}"
export PIPELINE_PLOTS_DIR="${DIR}/${DIR_PLOTS}"
export PIPELINE_THRESHOLD="${PEAK_LOGCPM_THRESHOLD}"
export PIPELINE_DAR_FDR="${DAR_FDR}"
export PIPELINE_DAR_LFC="${DAR_LFC}"

# =============================================================================
#  Step 1 – Trimming & Alignment  (1a=trim  1b=align  1c=dedup)
# =============================================================================
if should_run 1; then
  run_step 1 "Trimming & Alignment" bash "${SCRIPT_DIR}/01_mapping.sh"
fi

# =============================================================================
#  Step 2 – Peak Calling          (2a=MACS3  2b=summit→250bp)
# =============================================================================
if should_run 2; then
  run_step 2 "Peak Calling" bash "${SCRIPT_DIR}/02_peakcall.sh"
fi

# =============================================================================
#  Step 3 – Peak Counts           (3a=csaw+histogram  3b=filter)
# =============================================================================
if should_run 3; then
  echo ""
  echo "====== Step 3: Peak Counts ======"
  echo "  Started : $(date)"
  export PIPELINE_FROM_SUBSTEP="$(get_from_substep 3)"
  export PIPELINE_FORCE="${FORCE}"

  set +e
  Rscript "${SCRIPT_DIR}/03_peak_counts.R"
  _exit=$?
  set -e

  if [[ $_exit -eq 99 ]]; then
    echo ""
    echo "┌──────────────────────────────────────────────────────────┐"
    echo "│  ⏸  Pipeline PAUSED — threshold 未設定                  │"
    echo "│  Peak_nomodel/Histogram_peaks.png を確認し               │"
    echo "│  config.sh の PEAK_LOGCPM_THRESHOLD に値を設定後:        │"
    echo "│    ./run_pipeline.sh --from 3b   で再実行                │"
    echo "└──────────────────────────────────────────────────────────┘"
    exit 0
  elif [[ $_exit -ne 0 ]]; then
    echo "ERROR: Step 3 failed (exit code ${_exit})"
    exit $_exit
  fi
  echo "  Finished: $(date)"
fi

# =============================================================================
#  Step 4 – Scale Factors & BigWig (4a=bammerge  4b=scalefactor  4c=bigwig)
# =============================================================================
if should_run 4; then
  run_step 4 "Scale Factors & BigWig" bash "${SCRIPT_DIR}/04_scale_deeptools.sh"
fi

# =============================================================================
#  Step 5 – DAR Detection
# =============================================================================
if should_run 5; then
  run_step 5 "DAR Detection (edgeR QL)" Rscript "${SCRIPT_DIR}/05_DAR_edgeR.R"
fi

# =============================================================================
#  Step 6 – HOMER Motif Analysis
# =============================================================================
if should_run 6; then
  run_step 6 "HOMER Motif Analysis" bash "${SCRIPT_DIR}/06_HOMER.sh"
fi

# =============================================================================
#  Step 7 – PCA & Visualization
# =============================================================================
if should_run 7; then
  run_step 7 "PCA & Visualization" Rscript "${SCRIPT_DIR}/07_PCA_plots.R"
fi

echo ""
echo "=============================================="
echo " Pipeline COMPLETE: $(date)"
echo " Log: ${LOG}"
echo " Provenance: ${DIR}/provenance.yml"
echo "=============================================="
