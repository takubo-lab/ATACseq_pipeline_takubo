#!/usr/bin/env bash
# =============================================================================
#  setup_env.sh — ATAC-seq パイプライン環境セットアップ
#
#  使い方:
#    bash setup_env.sh          # 環境を新規作成
#    bash setup_env.sh --update # 既存環境を更新
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENV_NAME="atacseq_takubo"
ENV_FILE="${SCRIPT_DIR}/environment.yml"

# --- 色付きメッセージ ---
info()  { echo -e "\033[1;34m[INFO]\033[0m  $*"; }
warn()  { echo -e "\033[1;33m[WARN]\033[0m  $*"; }
error() { echo -e "\033[1;31m[ERROR]\033[0m $*"; exit 1; }

# --- conda / mamba 検出 ---
detect_installer() {
  if command -v mamba &>/dev/null; then
    echo "mamba"
  elif command -v conda &>/dev/null; then
    echo "conda"
  else
    echo ""
  fi
}

install_miniforge() {
  info "Miniforge (conda + mamba) をインストールします..."
  local installer="Miniforge3-$(uname)-$(uname -m).sh"
  wget -q "https://github.com/conda-forge/miniforge/releases/latest/download/${installer}" \
    -O /tmp/${installer}
  bash /tmp/${installer} -b -p "${HOME}/miniforge3"
  rm /tmp/${installer}

  # シェル初期化
  eval "$(${HOME}/miniforge3/bin/conda shell.bash hook)"
  conda init bash 2>/dev/null || true
  conda init zsh  2>/dev/null || true
  info "Miniforge をインストールしました。"
  info "シェルを再起動してから、再度このスクリプトを実行してください。"
  exit 0
}

# --- ツール検証 ---
verify_tools() {
  info "インストールされたツールを検証中..."
  eval "$(conda shell.bash hook)"
  conda activate "${ENV_NAME}"

  local tools=("bowtie2" "samtools" "macs3" "bedtools" "bamCoverage" "fastqc" "trim_galore" "Rscript")
  local all_ok=true
  for tool in "${tools[@]}"; do
    if command -v "${tool}" &>/dev/null; then
      local ver=""
      case "${tool}" in
        bowtie2)     ver=$(bowtie2 --version 2>&1 | head -1 | grep -oP '[\d.]+$' || echo "OK") ;;
        samtools)    ver=$(samtools --version 2>&1 | head -1 | awk '{print $2}') ;;
        macs3)       ver=$(macs3 --version 2>&1 || echo "OK") ;;
        bedtools)    ver=$(bedtools --version 2>&1 | awk '{print $2}') ;;
        bamCoverage) ver=$(bamCoverage --version 2>&1 || echo "OK") ;;
        fastqc)      ver=$(fastqc --version 2>&1 | awk '{print $2}') ;;
        trim_galore) ver=$(trim_galore --version 2>&1 | grep -i version | xargs || echo "OK") ;;
        Rscript)     ver=$(Rscript --version 2>&1 | grep -oP '[\d.]+' | head -1 || echo "OK") ;;
      esac
      info "  ✓ ${tool}: ${ver}"
    else
      warn "  ✗ ${tool}: 見つかりません"
      all_ok=false
    fi
  done

  if ${all_ok}; then
    info "全ツールの検証完了 ✓"
  else
    warn "一部ツールが見つかりません。手動で確認してください。"
  fi
}

# --- メイン ---
main() {
  info "=== ATAC-seq Pipeline 環境セットアップ ==="
  info "環境名: ${ENV_NAME}"
  info "定義ファイル: ${ENV_FILE}"
  echo ""

  # 1. conda/mamba 確認
  INSTALLER=$(detect_installer)
  if [[ -z "${INSTALLER}" ]]; then
    warn "conda/mamba が見つかりません。"
    read -rp "Miniforge を自動インストールしますか？ [y/N]: " yn
    case "${yn}" in
      [Yy]*) install_miniforge ;;
      *)     error "conda または mamba をインストールしてから再実行してください。" ;;
    esac
  fi
  info "パッケージマネージャ: ${INSTALLER}"

  # 2. 既存環境の確認
  if conda env list 2>/dev/null | grep -qw "${ENV_NAME}"; then
    warn "環境 '${ENV_NAME}' は既に存在します。"
    if [[ "${1:-}" == "--update" ]]; then
      info "環境を更新します..."
      ${INSTALLER} env update -n "${ENV_NAME}" -f "${ENV_FILE}" --prune
      verify_tools
      info "=== 環境の更新が完了しました ==="
      info "使用方法: conda activate ${ENV_NAME}"
      exit 0
    fi
    read -rp "再作成しますか（既存環境を削除）？ [y/N]: " yn
    case "${yn}" in
      [Yy]*)
        info "既存環境を削除中..."
        conda env remove -n "${ENV_NAME}" -y
        ;;
      *)
        info "既存環境を更新します..."
        ${INSTALLER} env update -n "${ENV_NAME}" -f "${ENV_FILE}" --prune
        verify_tools
        info "=== 環境の更新が完了しました ==="
        info "使用方法: conda activate ${ENV_NAME}"
        exit 0
        ;;
    esac
  fi

  # 3. 環境作成
  info "環境 '${ENV_NAME}' を作成中（10〜20分かかる場合があります）..."
  ${INSTALLER} env create -f "${ENV_FILE}"

  # 4. 検証
  verify_tools

  echo ""
  info "=== セットアップ完了 ==="
  info ""
  info "使用方法:"
  info "  conda activate ${ENV_NAME}"
  info "  bash run_pipeline.sh"
}

main "$@"
