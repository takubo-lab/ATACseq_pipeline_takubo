# ATAC-seq Pipeline

ペアエンド ATAC-seq データの全工程を自動化するパイプラインです。  
Trim → Align → Peak Call → Peak Count → Scale/BigWig → DAR → HOMER → PCA の順に実行します。

---

## ディレクトリ構成

```
./
├── config.sh            ← 【唯一の設定ファイル】実験ごとにここだけ編集
├── samples.tsv          ← サンプル名とグループの対応表
├── run_pipeline.sh      ← パイプライン実行スクリプト
├── scripts/
│   ├── 01_mapping.sh         # Step 1: trim_galore → bowtie2 → picard
│   ├── 02_peakcall.sh        # Step 2: MACS3 → 250bp 固定長ピーク
│   ├── 03_peak_counts.R      # Step 3: csaw カウント → 閾値フィルタ
│   ├── 04_scale_deeptools.sh # Step 4: BAM マージ → scale factor → bigWig
│   ├── 04a_scale_factor.R    #         (scale factor 計算 R スクリプト)
│   ├── 05_DAR_edgeR.R        # Step 5: DAR 検出 (edgeR QL)
│   ├── 06_HOMER.sh           # Step 6: HOMER モチーフ解析
│   └── 07_PCA_plots.R        # Step 7: PCA・相関ヒートマップ・Venn
└── fastq/
    └── {sample_name}_{lane}_[1|2].fq.gz
```

### 出力ディレクトリ (実行後に自動生成)

```
trimmed/          # trim_galore 後 FASTQ
BAM/              # 最終 BAM (.noDup.noMT.filt.sorted.bam)
Peak_nomodel/     # ピーク BED・カウント行列・ヒストグラム・DAR 結果
MergedBAM/        # グループ別マージ BAM
bw_for_deeptools/ # bigWig ファイル
Plots/            # PCA・ヒートマップ・散布図
logs/             # 実行ログ (タイムスタンプ付き)
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
`GenomicRanges`, `csaw`, `edgeR`, `ggplot2`, `pheatmap`, `tidyverse`, `data.table`, `BiocParallel`, `ggvenn`（任意）

---

## セットアップ手順

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
| `2a` | MACS3 callpeak (全サンプル統合、BAMPE モード) |
| `2b` | サミット → blacklist 除去 → 250 bp 固定長ピーク |
| `3a` | csaw カウント + aveLogCPM ヒストグラム出力 ※未設定時は一時停止 |
| `3b` | 閾値フィルタ → ピーク BED / カウント行列出力 |
| `4a` | グループごと BAM マージ + blacklist 除去 |
| `4b` | スケールファクター算出 (`04a_scale_factor.R`) |
| `4c` | bamCoverage → bigWig 生成 (CPM + scale factor) |
| `5` | DAR 検出 (edgeR QL、全ペアワイズ比較) |
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

| 変数 | 説明 | デフォルト |
|---|---|---|
| `BOWTIE2_MAX_INSERT` | bowtie2 最大 insert size (bp) | `700` |
| `FILTER_NFR` | NFR フィルタ (`true`/`false`) | `false` |
| `NFR_MAXFRAG` | NFR フィルタ時の最大フラグメント長 (bp) | `200` |

### Step 2: ピークコール

| 変数 | 説明 | デフォルト |
|---|---|---|
| `MACS3_PVALUE` | MACS3 p-value カットオフ | `0.01` |
| `MACS3_FORMAT` | `BAMPE` (推奨) / `BAM` | `BAMPE` |
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
| `DAR_LFC` | 最小 \|logFC\| | `1` |

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
