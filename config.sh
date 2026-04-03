#!/usr/bin/env bash
# =============================================================================
#  ATAC-seq Pipeline — Global Configuration
#  【このファイルのみを編集して実験ごとにパイプラインを設定する】
#
#  使い方:
#    1. DIR, GENOME*, FILE PATHS を自分の環境に合わせて編集
#    2. samples.tsv に サンプル名とグループ名を記入
#    3. ./run_pipeline.sh で全ステップ実行
#       ./run_pipeline.sh --from 3        # step 3 以降を実行
#       ./run_pipeline.sh --steps 1,2,5  # 指定ステップのみ実行
# =============================================================================

# --- プロジェクトのルートディレクトリ (このファイルがある場所を自動検出) ---
DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# --- サンプルシート ---
# sample_name <TAB> group の TSV (ヘッダー行あり)
# sample_name は fastq/{sample_name}_{lane}_[12].fq.gz の {sample_name} に対応
SAMPLES_TSV="${DIR}/samples.tsv"

# =============================================================================
#  GENOME SETTINGS
#  Human (hg38): GENOME="hg38"  GENOME_SIZE_MACS="hs"
#                STANDARD_CHR=( $(printf 'chr%s ' {1..22} X) )
#  Mouse (mm10): GENOME="mm10"  GENOME_SIZE_MACS="mm"
#                STANDARD_CHR=( $(printf 'chr%s ' {1..19} X) )
# =============================================================================
GENOME="hg38"                        # bowtie2 インデックス名
GENOME_SIZE_MACS="hs"                # MACS3 -g フラグ (hs / mm)
HOMER_GENOME="hg38"                  # HOMER findMotifsGenome.pl 用ゲノム名

# 保持する染色体 (chrM, chrY, chrUn, _random を自動除外)
STANDARD_CHR=( $(printf 'chr%s ' {1..22} X) )

# =============================================================================
#  FILE PATHS — 自環境に合わせて変更
# =============================================================================
BOWTIE2_INDEX="/home/stemcell/Genome/hg38/hg38"
BLACKLIST="/home/stemcell/Genome/hg38/hg38_blacklist_v2.bed"
CHROM_SIZES="/home/stemcell/Genome/hg38/hg38.chrom.sizes"
PICARD_JAR="$HOME/picard/picard.jar"

# =============================================================================
#  ENVIRONMENT
#  すでに PATH に samtools/bowtie2/macs3/deeptools が入っている場合は "" に
# =============================================================================
MACS3_ENV="$HOME/venvs/venv_macs3/bin/activate"

# =============================================================================
#  PERFORMANCE
# =============================================================================
THREADS=6
MAX_RAM_RECORDS=2500000   # picard MAX_RECORDS_IN_RAM (≈ RAM(MB) / 4 を目安)

# =============================================================================
#  ASSAY TYPE — アッセイタイプ別パラメータ対応表
#
#  このパイプラインは ATAC-seq / CUT&Tag / CUT&RUN に対応しています。
#  下記の対応表を参照し、使用するアッセイに合わせて各パラメータを設定してください。
#
#  ┌──────────────┬──────────────────────────────────────┬──────────────────────────────────────────────┐
#  │              │  Bowtie2 / Trimming                  │  MACS3                                       │
#  │              │  -X    │ MAPQ  │ Trim品質/長さ        │  -f    │ nomodel │ extsize │ shift │ Peak幅  │
#  ├──────────────┼────────┼───────┼─────────────────────┼────────┼─────────┼─────────┼───────┼─────────┤
#  │ ATAC-seq     │  700   │  -    │ デフォルト           │  BAM   │  true   │   200   │ -100  │  125bp  │
#  │ CUT&Tag      │  700   │  -    │ デフォルト           │  BAMPE │  true   │    -    │   -   │  500bp  │
#  │ CUT&RUN      │ 1000   │  30   │ --quality 30 -L 15  │  BAMPE │  true   │    -    │   -   │  125bp  │
#  └──────────────┴────────┴───────┴─────────────────────┴────────┴─────────┴─────────┴───────┴─────────┘
#
ASSAY_TYPE="ATAC"   # 参考ラベル (ATAC / CUT_AND_TAG / CUT_AND_RUN)
                    # ※ 実際の動作は下記の個別パラメータで決まります

# =============================================================================
#  STEP 1 — MAPPING (trim_galore / bowtie2 / picard)
# =============================================================================

# --- Bowtie2 ---
# ATAC / CUT&Tag : 700   CUT&RUN : 1000
BOWTIE2_MAX_INSERT=700    # bowtie2 -X (最大 insert size, bp)

# samtools MAPQ フィルタ (0 = フィルタなし)
# ATAC / CUT&Tag : 0     CUT&RUN : 30
MAPQ_FILTER=0

# --- trim_galore 追加オプション ---
# ATAC / CUT&Tag : デフォルト (TRIM_QUALITY=20, TRIM_MIN_LENGTH は省略)
# CUT&RUN        : TRIM_QUALITY=30  TRIM_MIN_LENGTH=15
TRIM_QUALITY=20           # trim_galore --quality
TRIM_MIN_LENGTH=""        # trim_galore --length (空 = デフォルト、適用しない)

# Nucleosome-free region (NFR) フィルタ
# true: フラグメント長 < NFR_MAXFRAG のリードのみ使用 (NFR解析)
# false: 全フラグメントを使用 (現行パイプラインと同等)
FILTER_NFR="false"
NFR_MAXFRAG=200

# =============================================================================
#  STEP 2 — PEAK CALLING (MACS3)
# =============================================================================

# --- MACS3 オプション ---
#                    ATAC     CUT&Tag   CUT&RUN
MACS3_PVALUE="0.01"          # p-value カットオフ
MACS3_FORMAT="BAM"           # -f : BAM (ATAC) / BAMPE (CUT&Tag/CUT&RUN)
MACS3_NOMODEL="true"         # --nomodel : ATAC=true / CUT&Tag=true / CUT&RUN=true
MACS3_EXTSIZE="200"          # --extsize : ATAC=200 / CUT&Tag="" / CUT&RUN=""  (空=省略)
MACS3_SHIFT="-100"           # --shift   : ATAC=-100 / CUT&Tag="" / CUT&RUN=""  (空=省略)

# CUT&RUN などでバックグラウンド BAM を使う場合に指定
# 例:
#   MACS3_CONTROL_BAMS=(
#     "${DIR}/controls/IgG_rep1.bam"
#     "${DIR}/controls/IgG_rep2.bam"
#   )
# 空配列なら control なしで実行
MACS3_CONTROL_BAMS=()

# サミットを中心とした固定長ピークウィンドウの half-width (bp)
# 最終ピーク幅 = SUMMIT_HALFWIDTH × 2
# ATAC / CUT&RUN : 125 (→ 250bp)   CUT&Tag (H3K27ac/H3K9ac) : 500 (→ 1000bp)
SUMMIT_HALFWIDTH=125

# =============================================================================
#  STEP 3 — PEAK COUNTS (csaw)
# =============================================================================
# aveLogCPM フィルタ閾値 (この値未満のピークを除外)
#
#   "" (空)   : ヒストグラム (Peak_nomodel/Histogram_peaks.png) と
#               分布データ (Peak_nomodel/peak_logcpm.tsv) を出力して
#               一時停止する。プロットを確認して数値を入力し、
#               ./run_pipeline.sh --from 3 で再実行。
#   数値      : 例) "-1" / "0" / "1" — 指定値未満のピークを除去して続行
#   "auto"    : ヒストグラムの最初の谷を自動検出 (確認不要な場合)
PEAK_LOGCPM_THRESHOLD=""

# =============================================================================
#  STEP 4 — SCALE FACTORS & DEEPTOOLS
# =============================================================================
BINSIZE=25                # bigWig ビンサイズ (bp)

# =============================================================================
#  STEP 5 — DAR DETECTION (edgeR)
# =============================================================================
DAR_FDR="0.05"            # FDR カットオフ
DAR_LFC="1"               # 最小 |logFC|

# =============================================================================
#  STEP 6 — HOMER MOTIF ANALYSIS
# =============================================================================
HOMER_FDR="${DAR_FDR}"    # HOMER 入力ピーク選定の FDR カットオフ
HOMER_LFC="${DAR_LFC}"    # HOMER 入力ピーク選定の |logFC| カットオフ
HOMER_SIZE=250            # findMotifsGenome.pl -size パラメータ

# =============================================================================
#  OUTPUT DIRECTORIES (${DIR} からの相対パス)
# =============================================================================
DIR_TRIMMED="trimmed"
DIR_BAM="BAM"
DIR_PEAKS="Peak_nomodel"
DIR_MERGEDBAM="MergedBAM"
DIR_BW="bw_for_deeptools"
DIR_PLOTS="Plots"
SCRIPT_DIR="${DIR}/scripts"
