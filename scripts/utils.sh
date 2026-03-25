#!/usr/bin/env bash
# =============================================================================
#  ATAC-seq Pipeline — Utility Functions
#
#  バージョン管理、Provenance 記録、ロック機構などの共通ユーティリティ。
#  run_pipeline.sh から source して使用する。
# =============================================================================

# --- バージョン取得 -----------------------------------------------------------
get_pipeline_version() {
  local version_file="${PIPELINE_ROOT}/VERSION"
  if [[ -f "$version_file" ]]; then
    cat "$version_file" | tr -d '[:space:]'
  else
    echo "unknown"
  fi
}

# --- Git コミットハッシュ取得 -------------------------------------------------
get_git_commit() {
  if command -v git &>/dev/null && git -C "${PIPELINE_ROOT}" rev-parse HEAD &>/dev/null 2>&1; then
    local hash
    hash=$(git -C "${PIPELINE_ROOT}" rev-parse --short HEAD 2>/dev/null)
    local dirty=""
    if ! git -C "${PIPELINE_ROOT}" diff --quiet 2>/dev/null; then
      dirty="-dirty"
    fi
    echo "${hash}${dirty}"
  else
    echo "no-git"
  fi
}

# --- ソフトウェアバージョン取得 -----------------------------------------------
get_tool_version() {
  local tool="$1"
  case "$tool" in
    bowtie2)
      bowtie2 --version 2>/dev/null | head -1 | grep -oP '[\d.]+$' || echo "N/A"
      ;;
    samtools)
      samtools --version 2>/dev/null | head -1 | awk '{print $2}' || echo "N/A"
      ;;
    macs3)
      macs3 --version 2>/dev/null | awk '{print $2}' || echo "N/A"
      ;;
    bedtools)
      bedtools --version 2>/dev/null | awk '{print $2}' || echo "N/A"
      ;;
    trim_galore)
      trim_galore --version 2>/dev/null | grep -oP '[\d.]+' | head -1 || echo "N/A"
      ;;
    picard)
      java -jar "${PICARD_JAR}" MarkDuplicates --version 2>&1 | head -1 || echo "N/A"
      ;;
    deeptools)
      bamCoverage --version 2>/dev/null | awk '{print $2}' || echo "N/A"
      ;;
    R)
      Rscript --version 2>&1 | grep -oP '[\d.]+' | head -1 || echo "N/A"
      ;;
    *)
      echo "N/A"
      ;;
  esac
}

# --- Provenance (実行記録) 生成 -----------------------------------------------
#  パイプライン実行時に provenance.yml を出力ディレクトリに自動生成する。
#  記録内容: パイプラインバージョン、gitコミット、実行日時、ユーザー、ホスト名、
#            主要ソフトウェアバージョン、config.sh のハッシュ
generate_provenance() {
  local project_dir="$1"
  local provenance_file="${project_dir}/provenance.yml"
  local config_file="${project_dir}/config.sh"

  local pipeline_version
  pipeline_version=$(get_pipeline_version)
  local git_commit
  git_commit=$(get_git_commit)
  local config_hash="N/A"
  if [[ -f "$config_file" ]]; then
    config_hash=$(sha256sum "$config_file" 2>/dev/null | awk '{print $1}' || echo "N/A")
  fi

  cat > "$provenance_file" <<EOF
# =============================================================================
#  ATAC-seq Pipeline — Provenance Record
#  This file is auto-generated at each pipeline run. Do not edit manually.
# =============================================================================

pipeline:
  name: ATACseq_pipeline_takubo
  version: "${pipeline_version}"
  git_commit: "${git_commit}"
  repository: "https://github.com/takubo-lab/ATACseq_pipeline_takubo.git"

execution:
  date: "$(date -Iseconds)"
  user: "$(whoami)"
  hostname: "$(hostname)"
  working_directory: "${project_dir}"
  config_sha256: "${config_hash}"

software:
  bowtie2: "$(get_tool_version bowtie2)"
  samtools: "$(get_tool_version samtools)"
  macs3: "$(get_tool_version macs3)"
  bedtools: "$(get_tool_version bedtools)"
  trim_galore: "$(get_tool_version trim_galore)"
  picard: "$(get_tool_version picard)"
  deeptools: "$(get_tool_version deeptools)"
  R: "$(get_tool_version R)"

parameters:
  genome: "${GENOME}"
  threads: ${THREADS}
  bowtie2_max_insert: ${BOWTIE2_MAX_INSERT}
  filter_nfr: "${FILTER_NFR}"
  macs3_pvalue: "${MACS3_PVALUE}"
  macs3_format: "${MACS3_FORMAT}"
  macs3_control_bams: "${MACS3_CONTROL_BAMS[*]:-}"
  summit_halfwidth: ${SUMMIT_HALFWIDTH}
  peak_logcpm_threshold: "${PEAK_LOGCPM_THRESHOLD}"
  binsize: ${BINSIZE}
  dar_fdr: "${DAR_FDR}"
  dar_lfc: "${DAR_LFC}"
EOF

  echo "  Provenance recorded: ${provenance_file}"
}

# --- ロック機構 ---------------------------------------------------------------
#  同一プロジェクトでの同時実行を防止する。
#  .pipeline.lock ファイルにPID、ユーザー、ホスト名、開始時刻を記録。
LOCK_FILE=""

acquire_lock() {
  local project_dir="$1"
  LOCK_FILE="${project_dir}/.pipeline.lock"

  if [[ -f "$LOCK_FILE" ]]; then
    local lock_pid lock_user lock_host lock_date
    lock_pid=$(grep '^pid:' "$LOCK_FILE" | awk '{print $2}')
    lock_user=$(grep '^user:' "$LOCK_FILE" | awk '{print $2}')
    lock_host=$(grep '^host:' "$LOCK_FILE" | awk '{print $2}')
    lock_date=$(grep '^started:' "$LOCK_FILE" | cut -d' ' -f2-)

    # 同一ホストでプロセスがまだ生きているか確認
    if [[ "$(hostname)" == "$lock_host" ]] && kill -0 "$lock_pid" 2>/dev/null; then
      echo "ERROR: Pipeline is already running on this project."
      echo "  PID    : ${lock_pid}"
      echo "  User   : ${lock_user}"
      echo "  Host   : ${lock_host}"
      echo "  Started: ${lock_date}"
      echo ""
      echo "  If the previous run crashed, remove the lock file manually:"
      echo "    rm ${LOCK_FILE}"
      exit 1
    else
      echo "  WARNING: Stale lock file detected (PID ${lock_pid} on ${lock_host} is not running)."
      echo "  Removing stale lock and proceeding..."
      rm -f "$LOCK_FILE"
    fi
  fi

  # ロックファイルを作成
  cat > "$LOCK_FILE" <<EOF
pid: $$
user: $(whoami)
host: $(hostname)
started: $(date -Iseconds)
EOF

  # シグナルハンドラで終了時にロック解除
  trap 'release_lock' EXIT INT TERM HUP
  echo "  Lock acquired: ${LOCK_FILE} (PID $$)"
}

release_lock() {
  if [[ -n "$LOCK_FILE" && -f "$LOCK_FILE" ]]; then
    rm -f "$LOCK_FILE"
  fi
}

# --- バージョン表示バナー -----------------------------------------------------
print_version_banner() {
  local version
  version=$(get_pipeline_version)
  local commit
  commit=$(get_git_commit)
  echo " Version: v${version} (${commit})"
}
