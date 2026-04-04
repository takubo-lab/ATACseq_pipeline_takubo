# ATAC-seq Pipeline (Takubo Lab)

ペアエンド ATAC-seq データの全工程を自動化するパイプラインです。  
Trim → Align → Peak Call → Peak Count → Scale/BigWig → DAR → HOMER → PCA → peak-set enrichment の順に実行します。

## パイプライン概要

| Step | Script | 内容 |
|------|--------|------|
| 1 | `scripts/01_mapping.sh` | trim_galore → bowtie2 → picard (トリミング・アライメント・重複除去) |
| 2 | `scripts/02_peakcall.sh` | MACS3 → 250bp 固定長ピーク |
| 3 | `scripts/03_peak_counts.R` | csaw カウント → 閾値フィルタ |
| 4 | `scripts/04_scale_deeptools.sh` | BAM マージ → scale factor → bigWig |
| 5 | `scripts/05_DAR_edgeR.R` | DAR 検出 (edgeR LRT) |
| 6 | `scripts/06_HOMER.sh` | HOMER モチーフ解析 |
| 7 | `scripts/07_PCA_plots.R` | PCA・相関ヒートマップ・Venn 図 |
| 8 | `scripts/08_peakset_enrichment.sh` | 外部 accessibility site set の正規化と peak-set enrichment |

## クイックスタート

### 1. セットアップ

```bash
# リポジトリをクローン
git clone https://github.com/takubo-lab/ATACseq_pipeline_takubo.git
cd ATACseq_pipeline_takubo

# 環境を一括セットアップ (Conda/Mamba)
bash setup_env.sh
conda activate atacseq_takubo
```

> **WSL (Windows) の場合も同じ手順で動作します。** Conda未導入の場合は `setup_env.sh` が Miniforge を自動インストールします。

### 2-A. プロジェクト初期化（推奨）

パイプラインコードとデータを分離して管理する方法:

```bash
# 新しいプロジェクトを初期化
bash init_project.sh /path/to/my_project hg38

# サンプル情報を編集
vim /path/to/my_project/samples.tsv

# FASTQファイルを配置
cp *.fq.gz /path/to/my_project/fastq/

# パイプライン実行
bash run_pipeline.sh --config /path/to/my_project/config.sh
```

### 2-B. 直接実行（シンプル）

パイプラインディレクトリ内で直接解析する方法:

**`config.sh`** — ゲノム、パス、パラメータを設定:
```bash
GENOME="hg38"              # "hg38" or "mm10"
GENOME_SIZE_MACS="hs"      # "hs" or "mm"
```

**`samples.tsv`** — サンプル情報を記載:
```
sample_name	group
Old-Mock-1	Old_Mock
Old-Mock-2	Old_Mock
Young-Ctrl-1	Young_Ctrl
```

### 3. FASTQファイルを配置

```bash
# fastq/ ディレクトリにfq.gzファイルを配置
# ファイル名: {sample_name}_{lane}_{1,2}.fq.gz
#         or: {sample_name}_R{1,2}_{anything}.fastq.gz
```

### 4. パイプライン実行

```bash
# 全ステップ実行
bash run_pipeline.sh

# 特定のステップから開始
bash run_pipeline.sh --from 4

# 特定のステップのみ実行
bash run_pipeline.sh --steps 5,6,7,8

# substep 指定
bash run_pipeline.sh --from 4b

# プロジェクト指定で実行
bash run_pipeline.sh --config /path/to/my_project/config.sh

# 既存出力を無視して強制再実行
bash run_pipeline.sh --force --from 2a

# バージョン表示
bash run_pipeline.sh --version
```

---

## バージョン管理・マルチユーザー運用

### 環境セットアップ

#### Linux フルパイプライン: Conda

Step 1〜7 をすべて実行する場合 (Linux / WSL):

```bash
bash setup_env.sh
conda activate atacseq_takubo
bash run_pipeline.sh
```

#### R 解析のみ: renv (Windows / Mac / Linux 共通)

Step 3, 5, 7 など R スクリプトのみを手元の Windows/Mac で実行する場合:

```r
# ターミナルから
Rscript setup_renv.R

# または RStudio コンソールで
renv::restore()
```

> **RStudio** でプロジェクトを開くと `.Rprofile` により renv が自動アクティベートされます。
> `renv::restore()` を実行するだけで `renv.lock` に記載された全パッケージが同一バージョンでインストールされます。

| 環境 | 用途 | セットアップ |
|------|------|-------------|
| Conda (`environment.yml`) | Linux フルパイプライン (Bash + R) | `bash setup_env.sh` |
| renv (`renv.lock`) | R 解析のみ (Windows / Mac / Linux) | `Rscript setup_renv.R` |

### バージョニング

パイプラインは `VERSION` ファイルでセマンティックバージョニングを管理:
- パッチ (1.0.x): バグ修正
- マイナー (1.x.0): 機能追加（後方互換）
- メジャー (x.0.0): 破壊的変更

リリース時は git tag を使用:
```bash
git tag -a v1.1.0 -m "Add new feature"
git push origin v1.1.0
```

### Provenance (実行記録)

パイプライン実行時に `provenance.yml` が自動生成され、以下を記録:
- パイプラインバージョン・コミット
- 実行日時・ユーザー・ホスト名
- ソフトウェアバージョン (bowtie2, samtools, MACS3, R, deeptools 等)
- config.sh の SHA256 ハッシュ
- 主要パラメータ

### ロック機構

同一プロジェクトでの同時実行を防止する `.pipeline.lock` ファイルによるロック機構を搭載。
前回の実行がクラッシュした場合は手動で削除:
```bash
rm /path/to/project/.pipeline.lock
```

---

## ディレクトリ構成

```
ATACseq_pipeline_takubo/
├── VERSION                # パイプラインバージョン
├── environment.yml        ← Conda 環境定義 (Linux フルパイプライン用)
├── setup_env.sh           ← Conda 環境セットアップスクリプト
├── renv.lock              ← R パッケージバージョンロック (全OS共通)
├── setup_renv.R           ← R 環境セットアップスクリプト
├── .Rprofile              ← renv 自動アクティベート
├── renv/                  ← renv 管理フォルダ (library/ は gitignore)
├── plot_config.default.yml ← 可視化設定テンプレート (Git管理)
├── plot_config.yml        ← 個人カスタマイズ (Git非管理, .gitignore)
├── config.sh              ← 【唯一の設定ファイル】実験ごとにここだけ編集
├── samples.tsv            ← サンプル名とグループの対応表
├── init_project.sh        ← プロジェクト初期化スクリプト
├── run_pipeline.sh        ← パイプライン実行スクリプト
├── scripts/
│   ├── utils.sh                # ユーティリティ関数 (バージョン・provenance・ロック)
│   ├── plot_utils.R            # 可視化設定読み込みユーティリティ
│   ├── 01_mapping.sh           # Step 1: trim_galore → bowtie2 → picard
│   ├── 02_peakcall.sh          # Step 2: MACS3 → 250bp 固定長ピーク
│   ├── 03_peak_counts.R        # Step 3: csaw カウント → 閾値フィルタ
│   ├── 04_scale_deeptools.sh   # Step 4: BAM マージ → scale factor → bigWig
│   ├── 04a_scale_factor.R      #         (scale factor 計算 R スクリプト)
│   ├── 05_DAR_edgeR.R          # Step 5: DAR 検出 (edgeR LRT)
│   ├── 06_HOMER.sh             # Step 6: HOMER モチーフ解析
│   ├── 07_PCA_plots.R          # Step 7: PCA・相関ヒートマップ・Venn
│   ├── 08_peakset_enrichment.sh # Step 8a/8b: peak-set enrichment wrapper
│   └── 08_peakset_fgsea.R      # Step 8b: fgsea-like peak-set enrichment
├── fastq/                 # 入力FASTQファイル (gitignore)
│   └── {sample_name}_{lane}_[1|2].fq.gz
├── trimmed/               # trim_galore 後 FASTQ (gitignore)
├── BAM/                   # 最終 BAM (gitignore)
├── Peak_nomodel/          # ピーク BED・カウント行列・DAR 結果・peak-set enrichment (gitignore)
├── MergedBAM/             # グループ別マージ BAM (gitignore)
├── bw_for_deeptools/      # bigWig ファイル (gitignore)
├── Plots/                 # PCA・ヒートマップ等 (gitignore)
└── logs/                  # 実行ログ (gitignore)
```

### 出力ディレクトリ (実行後に自動生成)

```
trimmed/          # trim_galore 後 FASTQ
BAM/              # 最終 BAM (.noDup.noMT.filt.sorted.bam)
Peak_nomodel/     # ピーク BED・カウント行列・ヒストグラム・DAR 結果・PeakSetLibrary・PeakSetEnrichment
MergedBAM/        # グループ別マージ BAM
bw_for_deeptools/ # bigWig ファイル
Plots/            # PCA・ヒートマップ・散布図
logs/             # 実行ログ (タイムスタンプ付き)
provenance.yml    # 実行記録 (バージョン・ユーザー・ソフトウェア情報)
```

---

## 前提ツール

| ツール | 用途 |
|---|---|
| trim_galore | アダプタートリミング |
| bowtie2 | アライメント |
| samtools | BAM 操作 |
| bedtools | Blacklist 除去・座標操作 |
| picard | 重複除去 (`MarkDuplicates`) |
| macs3 | ピークコール |
| bamCoverage (deeptools) | bigWig 生成 |
| HOMER | モチーフ解析 |
| R ≥ 4.1 | 統計・可視化 |

**必要な R パッケージ:**  
`GenomicRanges`, `csaw`, `edgeR`, `fgsea`, `ggplot2`, `pheatmap`, `tidyverse`, `data.table`, `BiocParallel`, `ggvenn`, `yaml`

## Step 8: Peak-set enrichment

Step 8 は project の DAR 結果を ranked peak list として扱い、外部 accessibility site BED を共通 peak universe に写像して fgsea-like enrichment を行います。

### Inputs

1. `Peak_nomodel/Peaks_250bp_th0.bed`: canonical peak universe
2. `Peak_nomodel/edgeR/*_DA_edgeR_Results.tsv`: ranked DAR results
3. `Differential_Peak_Iwama/*.bed`
4. `Differential_Peak_Kobayashi/*.bed`

### Outputs

1. `Peak_nomodel/PeakSetLibrary/region_sets_manifest.tsv`
2. `Peak_nomodel/PeakSetLibrary/normalized/{source}/*.bed`
3. `Peak_nomodel/PeakSetEnrichment/{comparison}/*_peakset_fgsea.tsv`
4. `Peak_nomodel/PeakSetEnrichment/peakset_fgsea_all_comparisons.tsv`

### Re-run examples

```bash
# Step 8 全体を再実行
bash run_pipeline.sh --steps 8

# region set 正規化だけやり直し
bash run_pipeline.sh --steps 8a

# enrichment 集計だけやり直し
bash run_pipeline.sh --steps 8b
```

---

## セットアップ手順

### 方法A: プロジェクト初期化（推奨 — マルチユーザー運用向け）

```bash
# リポジトリをクローン
git clone https://github.com/takubo-lab/ATACseq_pipeline_takubo.git
cd ATACseq_pipeline_takubo

# 新しいプロジェクトを初期化
bash init_project.sh /data/ATACseq/240101_experiment hg38

# サンプル情報を編集
vim /data/ATACseq/240101_experiment/samples.tsv

# FASTQファイルを配置
cp *.fq.gz /data/ATACseq/240101_experiment/fastq/

# パイプライン実行
bash run_pipeline.sh --config /data/ATACseq/240101_experiment/config.sh
```

### 方法B: 直接実行（シンプル）

### 1. `samples.tsv` を編集

```
sample_name	group
Old-Mock-1	Old_Mock
Old-Mock-2	Old_Mock
Young-Ctrl-1	Young_Ctrl
Young-Ctrl-2	Young_Ctrl
```

- **`sample_name`**: `fastq/{sample_name}_{lane}_[1|2].fq.gz` の `{sample_name}` に一致させる
- **`group`**: DAR 比較・bigWig 作成のグループ名（スペース不可; `_` 推奨）
- ヘッダー行 (`sample_name	group`) は必須

### 2. `config.sh` を編集

最低限、下記 5 箇所を確認・変更してください。

```bash
# ゲノム設定 (Human / Mouse)
GENOME="hg38"                   # または "mm10"
GENOME_SIZE_MACS="hs"           # または "mm"
HOMER_GENOME="hg38"             # または "mm10"
STANDARD_CHR=( $(printf 'chr%s ' {1..22} X) )   # Mouse は {1..19}

# ファイルパス
BOWTIE2_INDEX="/path/to/hg38/hg38"
BLACKLIST="/path/to/hg38_blacklist_v2.bed"
CHROM_SIZES="/path/to/hg38.chrom.sizes"
PICARD_JAR="$HOME/picard/picard.jar"

# MACS3 の仮想環境 (PATH にあれば "" に)
MACS3_ENV="$HOME/venvs/venv_macs3/bin/activate"
```

---

## 実行コマンド

```bash
# 全ステップ実行 (Step 1 → 7)
./run_pipeline.sh

# Step 3 以降を実行
./run_pipeline.sh --from 3

# Step 4 の substep b (scale factor) 以降を実行
./run_pipeline.sh --from 4b

# Step 1〜4 のみ実行
./run_pipeline.sh --to 4

# 特定ステップのみ実行 (カンマ区切り / substep 指定可)
./run_pipeline.sh --steps 5,6,7
./run_pipeline.sh --steps 2b,4c

# 既存出力を無視して強制再実行
./run_pipeline.sh --force --from 2a

# 外部プロジェクトの config.sh を指定して実行
./run_pipeline.sh --config /path/to/project/config.sh

# バージョン表示
./run_pipeline.sh --version

# ヘルプ表示
./run_pipeline.sh --help
```

---

## Substep 一覧

各ステップは内部で substep (a/b/c) に分割されており、`--from` や `--steps` で細かく再開できます。

| 指定 | 処理内容 |
|---|---|
| `1a` | trim_galore アダプタートリミング |
| `1b` | bowtie2 アライメント + samtools フィルタ (proper-pair / 標準染色体) |
| `1c` | picard MarkDuplicates (重複除去) |
| `2a` | MACS3 callpeak (全サンプル統合、`MACS3_FORMAT`/`MACS3_NOMODEL` に従う) |
| `2b` | サミット → blacklist 除去 → 250 bp 固定長ピーク |
| `3a` | csaw カウント + aveLogCPM ヒストグラム出力 ※未設定時は一時停止 |
| `3b` | 閾値フィルタ → ピーク BED / カウント行列出力 |
| `4a` | グループごと BAM マージ + blacklist 除去 |
| `4b` | スケールファクター算出 (`04a_scale_factor.R`) |
| `4c` | bamCoverage → bigWig 生成 (CPM + scale factor) |
| `5` | DAR 検出 (edgeR LRT、全ペアワイズ比較) |
| `6` | HOMER findMotifsGenome.pl (up/down DAR 別) |
| `7` | PCA / 相関ヒートマップ / 散布図 / Venn 図 |

---

## よくある再実行シナリオ

| 変更内容 | 再実行コマンド |
|---|---|
| `SUMMIT_HALFWIDTH` を変更した | `./run_pipeline.sh --from 2b` |
| `MACS3_PVALUE` を変更した | `./run_pipeline.sh --from 2a` |
| 閾値ヒストグラムを確認後に `PEAK_LOGCPM_THRESHOLD` を設定 | `./run_pipeline.sh --from 3b` |
| `BINSIZE` (bigWig) を変更した | `./run_pipeline.sh --from 4c` |
| `DAR_FDR` / `DAR_LFC` を変更した | `./run_pipeline.sh --steps 5,6,7` |
| HOMER 設定を変更した | `./run_pipeline.sh --steps 6` |
| PCA だけやり直したい | `./run_pipeline.sh --steps 7` |
| サンプル追加後に Step 1 だけ再実行 | `./run_pipeline.sh --steps 1` |

---

## `PEAK_LOGCPM_THRESHOLD` の設定手順

Step 3 の閾値は **config.sh** の `PEAK_LOGCPM_THRESHOLD` で制御します。

```bash
# config.sh の該当部分
PEAK_LOGCPM_THRESHOLD=""      # (1) まず空のまま実行 → ヒストグラム確認
PEAK_LOGCPM_THRESHOLD="-1"    # (2) 谷の値を確認して数値を設定
PEAK_LOGCPM_THRESHOLD="auto"  # 自動検出 (確認省略したい場合)
```

**操作フロー:**

```
1. config.sh を PEAK_LOGCPM_THRESHOLD="" のまま実行
   ./run_pipeline.sh

2. パイプラインが Step 3a で一時停止し、以下を出力:
     Peak_nomodel/Histogram_peaks.png   ← aveLogCPM 分布のヒストグラム
     Peak_nomodel/peak_logcpm.tsv       ← 各ピークの aveLogCPM 値

3. ヒストグラムを開いてノイズピークと真のピークの境界を確認

4. config.sh を更新:
     PEAK_LOGCPM_THRESHOLD="-1"   ← 谷の値を入力

5. Step 3b から再実行 (csaw カウントはキャッシュを使用):
   ./run_pipeline.sh --from 3b
```

---

## アッセイタイプ別パラメータ設定

このパイプラインは ATAC-seq / CUT&Tag / CUT&RUN に対応しています。
`config.sh` で下記の対応表に従って各パラメータを設定してください。

### Bowtie2 / Trimming 設定

| パラメータ | ATAC-seq | CUT&Tag | CUT&RUN |
|---|---|---|---|
| `BOWTIE2_MAX_INSERT` | `700` | `700` | `1000` |
| `MAPQ_FILTER` | `0` (フィルタなし) | `0` (フィルタなし) | `30` |
| `TRIM_QUALITY` | `20` | `20` | `30` |
| `TRIM_MIN_LENGTH` | `""` (省略) | `""` (省略) | `15` |

### MACS3 設定

| パラメータ | ATAC-seq | CUT&Tag | CUT&RUN |
|---|---|---|---|
| `MACS3_FORMAT` | `BAM` | `BAMPE` | `BAMPE` |
| `MACS3_NOMODEL` | `true` | `true` | `true` |
| `MACS3_EXTSIZE` | `200` | `""` (省略) | `""` (省略) |
| `MACS3_SHIFT` | `-100` | `""` (省略) | `""` (省略) |
| `SUMMIT_HALFWIDTH` | `125` (→250bp) | `500` (→1000bp) | `125` (→250bp) |

**設定例 (CUT&RUN の場合):**
```bash
ASSAY_TYPE="CUT_AND_RUN"
BOWTIE2_MAX_INSERT=1000
MAPQ_FILTER=30
TRIM_QUALITY=30
TRIM_MIN_LENGTH=15
MACS3_FORMAT="BAMPE"
MACS3_NOMODEL="true"
MACS3_EXTSIZE=""
MACS3_SHIFT=""
MACS3_CONTROL_BAMS=(
  "/path/to/IgG_rep1.bam"
  "/path/to/IgG_rep2.bam"
)
SUMMIT_HALFWIDTH=125
```

### CUT&RUN のバックグラウンド BAM 指定

Step 2 のピークコールは treatment 側を全サンプル aggregate して 1 回だけ MACS3 を実行します。そのため、control 側も同じ粒度で 1 つにそろえて使う前提です。

`config.sh` で `MACS3_CONTROL_BAMS` に 0 本、1 本、複数本の BAM を指定できます。

```bash
# control なし
MACS3_CONTROL_BAMS=()

# control 1本
MACS3_CONTROL_BAMS=(
  "/path/to/IgG_merged.bam"
)

# control 複数本
MACS3_CONTROL_BAMS=(
  "/path/to/IgG_rep1.bam"
  "/path/to/IgG_rep2.bam"
  "/path/to/IgG_rep3.bam"
)
```

- 1 本指定: その BAM をそのまま `macs3 callpeak -c` に渡します。
- 複数本指定: Step 2 で `Peak_nomodel/aggregated/aggregated_control.bam` に merge して index を作成し、その merged BAM を `-c` に渡します。
- 未指定: control なしで peak call します。

条件ごとに別々の background がある場合でも、このパイプラインでは treatment を pooled peak call しているため、background も pooled control としてまとめて使う設計です。条件ごとの matched control を厳密に使い分けたい場合は、条件ごとに peak call して後で union peak を作る別設計が必要です。

**設定例 (CUT&Tag の場合):**
```bash
ASSAY_TYPE="CUT_AND_TAG"
BOWTIE2_MAX_INSERT=700
MAPQ_FILTER=0
TRIM_QUALITY=20
TRIM_MIN_LENGTH=""
MACS3_FORMAT="BAMPE"
MACS3_NOMODEL="true"
MACS3_EXTSIZE=""
MACS3_SHIFT=""
SUMMIT_HALFWIDTH=500
```

---

## `config.sh` パラメータ一覧

### ゲノム設定

| 変数 | 説明 | Human | Mouse |
|---|---|---|---|
| `GENOME` | ゲノム名 (bowtie2 インデックス識別子) | `hg38` | `mm10` |
| `GENOME_SIZE_MACS` | MACS3 `-g` フラグ | `hs` | `mm` |
| `HOMER_GENOME` | HOMER ゲノム名 | `hg38` | `mm10` |
| `STANDARD_CHR` | 保持する染色体リスト | `chr1-22,X` | `chr1-19,X` |

### ファイルパス

| 変数 | 説明 |
|---|---|
| `BOWTIE2_INDEX` | bowtie2 インデックスのプレフィックスパス |
| `BLACKLIST` | Blacklist BED ファイル (ENCODE blacklist v2 推奨) |
| `CHROM_SIZES` | 染色体サイズファイル (`{genome}.chrom.sizes`) |
| `PICARD_JAR` | picard.jar のフルパス |
| `MACS3_ENV` | MACS3 仮想環境の activate スクリプト (PATH にある場合は `""`) |

### パフォーマンス

| 変数 | 説明 | デフォルト |
|---|---|---|
| `THREADS` | 並列スレッド数 | `6` |
| `MAX_RAM_RECORDS` | picard MAX_RECORDS_IN_RAM | `2500000` |

### Step 1: マッピング

| 変数 | 説明 | デフォルト (ATAC) |
|---|---|---|
| `BOWTIE2_MAX_INSERT` | bowtie2 最大 insert size (bp) | `700` |
| `MAPQ_FILTER` | samtools MAPQ フィルタ (`0`=なし, CUT&RUNは`30`) | `0` |
| `TRIM_QUALITY` | trim_galore `--quality` | `20` |
| `TRIM_MIN_LENGTH` | trim_galore `--length` (空=省略, CUT&RUNは`15`) | `""` |
| `FILTER_NFR` | NFR フィルタ (`true`/`false`) | `false` |
| `NFR_MAXFRAG` | NFR フィルタ時の最大フラグメント長 (bp) | `200` |

### Step 2: ピークコール

| 変数 | 説明 | デフォルト (ATAC) |
|---|---|---|
| `MACS3_PVALUE` | MACS3 p-value カットオフ | `0.01` |
| `MACS3_FORMAT` | `-f` オプション (`BAM`=ATAC, `BAMPE`=CUT&Tag/CUT&RUN) | `BAM` |
| `MACS3_NOMODEL` | `--nomodel` を使用 (ATAC/CUT&Tag/CUT&RUN いずれも `true`) | `true` |
| `MACS3_EXTSIZE` | `--extsize` (ATACは`200`, CUT&Tag/CUT&RUNは`""`) | `200` |
| `MACS3_SHIFT` | `--shift` (ATACは`-100`, CUT&Tag/CUT&RUNは`""`) | `-100` |
| `MACS3_CONTROL_BAMS` | MACS3 control 用 BAM 配列。複数指定時は Step 2 で merge してから使用 | `()` |
| `SUMMIT_HALFWIDTH` | 固定長ピークの半幅 (bp)。最終ピーク = ×2 | `125` (250bp) |

### Step 3: ピークカウント

| 変数 | 説明 | デフォルト |
|---|---|---|
| `PEAK_LOGCPM_THRESHOLD` | aveLogCPM フィルタ閾値 (`""` / 数値 / `"auto"`) | `""` |

### Step 4: Deeptools

| 変数 | 説明 | デフォルト |
|---|---|---|
| `BINSIZE` | bigWig ビンサイズ (bp) | `25` |

### Step 5: DAR 検出

| 変数 | 説明 | デフォルト |
|---|---|---|
| `DAR_FDR` | FDR カットオフ | `0.05` |
| `DAR_LFC` | 最小 \|logFC\| | `0.5` |

### Step 6: HOMER

| 変数 | 説明 | デフォルト |
|---|---|---|
| `HOMER_FDR` | HOMER 入力ピーク選定 FDR | `DAR_FDR` と同値 |
| `HOMER_LFC` | HOMER 入力ピーク選定 \|logFC\| | `DAR_LFC` と同値 |
| `HOMER_SIZE` | `findMotifsGenome.pl -size` | `250` |

---

## 主要な出力ファイル

| ファイル | 説明 |
|---|---|
| `BAM/{sample}.noDup.noMT.filt.sorted.bam` | 最終 BAM (dedup, フィルタ済み) |
| `Peak_nomodel/Peak_250bp_noBL.bed` | 全ピーク (250bp, blacklist 除去済み) |
| `Peak_nomodel/Histogram_peaks.png` | aveLogCPM 分布ヒストグラム |
| `Peak_nomodel/peak_logcpm.tsv` | 各ピークの aveLogCPM 値 |
| `Peak_nomodel/Peaks_250bp_th0.bed` | 閾値フィルタ後ピーク |
| `Peak_nomodel/PeakCounts_th0.txt` | フィルタ後カウント行列 |
| `Peak_nomodel/scale_factors.tsv` | グループ別スケールファクター |
| `Peak_nomodel/edgeR/{A}_vs_{B}_DAR.tsv` | DAR 結果テーブル |
| `Peak_nomodel/edgeR/{A}_vs_{B}_volcano.png` | ボルケーノプロット |
| `bw_for_deeptools/{group}.bw` | スケール補正 bigWig |
| `Plots/PCA_PC1_PC2.png` 等 | PCA・相関ヒートマップ・Venn 図 |
| `logs/pipeline_{timestamp}.log` | 実行ログ |
| `provenance.yml` | 実行記録 (パイプラインバージョン・ソフトウェア情報) |

## 可視化のカスタマイズ（個人設定）

パイプラインの解析ロジック（Git管理）と、図の見た目（個人管理）を分離しています。

| ファイル | Git管理 | 用途 |
|----------|---------|------|
| `plot_config.default.yml` | ✓ 管理 | デフォルト設定（全員共通） |
| `plot_config.yml` | ✗ 無視 | 個人カスタマイズ |
| `scripts/plot_utils.R` | ✓ 管理 | 設定読み込みユーティリティ |

### カスタマイズ手順

```bash
# テンプレートをコピー（初回のみ）
cp plot_config.default.yml plot_config.yml

# 好みに合わせて編集
vim plot_config.yml   # or RStudio で編集
```

### 設定例

```yaml
# 図のテーマを変更
theme:
  base_size: 14
  ggplot_theme: "theme_classic"

# 出力をPDFに変更
output:
  format: "pdf"

# グループに特定の色を割り当て
colors:
  group_colors:
    WT: "#4DBBD5"
    KO: "#E64B35"

# PCA の点を大きく
pca:
  point_size: 6
```

### Rスクリプト内での使用

```r
source("scripts/plot_utils.R")
cfg <- load_plot_config()

# 統一テーマ適用
p <- ggplot(df, aes(x, y, color = group)) +
  geom_point(size = cfg$pca$point_size) +
  theme_pipeline(cfg)

# グループ色を取得
cols <- get_group_colors(unique(df$group), cfg)
p <- p + scale_color_manual(values = cols)

# 図を保存（サイズ・形式は設定ファイルに従う）
save_plot(p, "PCA_result", type = "pca", outdir = "results/", cfg = cfg)
```

> **ポイント**: `plot_config.yml` は `.gitignore` に含まれるため、個人の変更がリポジトリに反映されることはありません。デフォルト設定を変更したい場合は `plot_config.default.yml` を編集してコミットしてください。

---

## 対応ゲノム

| Genome | Bowtie2 Index | Blacklist | Chrom Sizes |
|--------|--------------|-----------|-------------|
| hg38 | `~/Genome/hg38/hg38` | `~/Genome/hg38/hg38_blacklist_v2.bed` | `~/Genome/hg38/hg38.chrom.sizes` |
| mm10 | `~/Genome/mm10/mm10` | `~/Genome/mm10/mm10_blacklist_v2.bed` | `~/Genome/mm10/mm10.chrom.sizes` |

パスは `config.sh` で変更可能。

## License

Internal use — Takubo Lab
