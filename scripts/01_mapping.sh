#!/usr/bin/env bash
# =============================================================================
#  Step 1: Trimming & Alignment
#
#  処理内容:
#    1. trim_galore --paired でアダプタートリミング
#    2. bowtie2 でゲノムアライメント (--very-sensitive, paired-end)
#    3. samtools で proper-pair (-f 3) かつ標準染色体のみ保持
#    4. picard MarkDuplicates で重複除去
#    5. 中間ファイル (.sam, 未フィルタ .bam) 削除
#
#  入力: fastq/{sample}_{lane}_{1,2}.fq.gz
#        または fastq/{sample}_R{1,2}_{anything}.fastq.gz
#  出力: BAM/{sample}.noDup.noMT.filt.sorted.bam (+.bai)
#        BAM/{sample}_picard_metrics.txt
# =============================================================================
set -euo pipefail
source "$(dirname "$0")/../config.sh"

# --- Substep control ---
# FROM_SUB: ""/"a"=trim+align+dedup  "b"=align+dedup  "c"=deduponly
FROM_SUB="${PIPELINE_FROM_SUBSTEP:-}"
FORCE="${PIPELINE_FORCE:-false}"

# returns 0 (run) if substep >= FROM_SUB
should_run_sub() {
  local sub="$1"
  [[ -z "${FROM_SUB}" || ! "${sub}" < "${FROM_SUB}" ]]
}

if [[ -n "${FROM_SUB}" ]]; then
  echo "  Substep start: ${FROM_SUB}  (skipping substeps before '${FROM_SUB}')"
fi

mkdir -p "${DIR}/${DIR_TRIMMED}" "${DIR}/${DIR_BAM}"

# --- helpers: samples.tsv から sample 一覧取得 ---
get_samples() {
  grep -v "^#" "${SAMPLES_TSV}" | grep -v "^$" | grep -v "^sample_name" | \
    awk -F'\t' '{print $1}'
}

# --- R1 ファイル検索 (複数のファイル名パターンに対応) ---
find_r1_files() {
  local sample="$1"
  local fastq_dir="${DIR}/fastq"
  shopt -s nullglob
  local files=(
    "${fastq_dir}/${sample}"*_R1_*.fastq.gz
    "${fastq_dir}/${sample}"*_R1_*.fq.gz
    "${fastq_dir}/${sample}"*_1.fq.gz
    "${fastq_dir}/${sample}"*_1.fastq.gz
  )
  shopt -u nullglob
  echo "${files[@]:-}"
}

# --- R1 から R2 パスを導出 ---
r2_from_r1() {
  local r1="$1"
  if [[ "$r1" == *"_R1_"* ]]; then
    echo "${r1/_R1_/_R2_}"
  elif [[ "$r1" == *"_1.fq.gz" ]]; then
    echo "${r1/_1.fq.gz/_2.fq.gz}"
  elif [[ "$r1" == *"_1.fastq.gz" ]]; then
    echo "${r1/_1.fastq.gz/_2.fastq.gz}"
  else
    echo "ERROR: Cannot determine R2 from R1: ${r1}" >&2
    exit 1
  fi
}

# --- trim_galore 出力の R1 val ファイルを取得 ---
trimmed_r1_from_r1() {
  local r1="$1"
  local base
  base=$(basename "$r1" | sed 's/\(\.fastq\.gz\|\.fq\.gz\)$//')
  echo "${DIR}/${DIR_TRIMMED}/${base}_val_1.fq.gz"
}
trimmed_r2_from_r2() {
  local r2="$1"
  local base
  base=$(basename "$r2" | sed 's/\(\.fastq\.gz\|\.fq\.gz\)$//')
  echo "${DIR}/${DIR_TRIMMED}/${base}_val_2.fq.gz"
}

# =============================================================================
#  MAIN LOOP: サンプルごとに処理
# =============================================================================
while IFS=$'\t' read -r sample _group; do
  [[ "$sample" =~ ^# ]] && continue
  [[ "$sample" == "sample_name" ]] && continue
  [[ -z "$sample" ]] && continue

  echo ""
  echo "---------- Processing sample: ${sample} ----------"

  # --- fastq ファイル検索 ---
  r1_list=( $(find_r1_files "$sample") )
  if [[ ${#r1_list[@]} -eq 0 ]]; then
    echo "WARNING: No R1 fastq found for sample '${sample}'. Skipping."
    continue
  fi

  # --- [a] trim_galore ---
  if should_run_sub a; then
    echo "[a] Adapter trimming..."
    for r1 in "${r1_list[@]}"; do
      r2=$(r2_from_r1 "$r1")
      if [[ ! -f "$r2" ]]; then echo "WARNING: R2 not found: ${r2}. Skipping."; continue; fi
      trim_r1=$(trimmed_r1_from_r1 "$r1")
      if [[ "${FORCE}" != "true" && -f "$trim_r1" ]]; then
        echo "  Trimmed exists, skipping: $(basename "$r1")"
      else
        echo "  Trimming: $(basename "$r1") + $(basename "$r2")"
        trim_galore --paired -o "${DIR}/${DIR_TRIMMED}/" -j 4 "$r1" "$r2"
      fi
    done
  fi

  # --- [b] bowtie2 alignment + samtools filter ---
  if should_run_sub b; then
    # 複数 lane の trimmed R1/R2 を収集
    trim_r1_files=()
    trim_r2_files=()
    for r1 in "${r1_list[@]}"; do
      r2=$(r2_from_r1 "$r1")
      tr1=$(trimmed_r1_from_r1 "$r1")
      tr2=$(trimmed_r2_from_r2 "$r2")
      [[ -f "$tr1" && -f "$tr2" ]] && trim_r1_files+=("$tr1") && trim_r2_files+=("$tr2")
    done
    if [[ ${#trim_r1_files[@]} -eq 0 ]]; then
      echo "WARNING: No trimmed files found for ${sample}. Run substep a first (or check fastq names)."
    else
      echo "[b] Aligning with bowtie2 (${#trim_r1_files[@]} lane(s))..."
      R1_ARG=$(printf ',%s' "${trim_r1_files[@]}"); R1_ARG="${R1_ARG:1}"
      R2_ARG=$(printf ',%s' "${trim_r2_files[@]}"); R2_ARG="${R2_ARG:1}"
      RAW_SAM="${DIR}/${DIR_BAM}/${sample}.sam"
      bowtie2 \
        --very-sensitive \
        -1 "${R1_ARG}" -2 "${R2_ARG}" \
        --no-mixed --no-discordant --phred33 \
        -I 10 -X "${BOWTIE2_MAX_INSERT}" \
        -x "${BOWTIE2_INDEX}" \
        -p "${THREADS}" \
        --met-file "${DIR}/${DIR_BAM}/${sample}_bowtie2.log" \
        -S "${RAW_SAM}"

      echo "  Filtering reads (standard chromosomes, proper pairs)..."
      FILT_BAM="${DIR}/${DIR_BAM}/${sample}.noMT.filt.bam"
      if [[ "${FILTER_NFR}" == "true" ]]; then
        samtools view -@ "${THREADS}" -b -f 3 "${RAW_SAM}" "${STANDARD_CHR[@]}" | \
          samtools view -b -@ "${THREADS}" -e "tlen >= 10 && tlen <= ${NFR_MAXFRAG}" \
          > "${FILT_BAM}"
      else
        samtools view -@ "${THREADS}" -b -f 3 "${RAW_SAM}" "${STANDARD_CHR[@]}" \
          > "${FILT_BAM}"
      fi
      FILT_SORTED="${DIR}/${DIR_BAM}/${sample}.noMT.filt.sorted.bam"
      samtools sort -@ "${THREADS}" -o "${FILT_SORTED}" "${FILT_BAM}"
      samtools index -@ "${THREADS}" "${FILT_SORTED}"
      rm -f "${FILT_BAM}" "${RAW_SAM}"
    fi
  fi

  # --- [c] picard MarkDuplicates ---
  if should_run_sub c; then
    FILT_SORTED="${DIR}/${DIR_BAM}/${sample}.noMT.filt.sorted.bam"
    if [[ ! -f "${FILT_SORTED}" ]]; then
      echo "WARNING: ${FILT_SORTED} not found. Run substep b first."
    else
      echo "[c] Removing duplicates with Picard..."
      BAM_FILE="${DIR}/${DIR_BAM}/${sample}.noDup.noMT.filt.sorted.bam"
      java -jar "${PICARD_JAR}" MarkDuplicates \
        REMOVE_DUPLICATES=true \
        MAX_RECORDS_IN_RAM="${MAX_RAM_RECORDS}" \
        ASSUME_SORT_ORDER=coordinate \
        I="${FILT_SORTED}" \
        O="${BAM_FILE}" \
        M="${DIR}/${DIR_BAM}/${sample}_picard_metrics.txt" \
        VALIDATION_STRINGENCY=LENIENT
      samtools index -@ "${THREADS}" "${BAM_FILE}"
      rm -f "${FILT_SORTED}" "${FILT_SORTED}.bai"
      READS=$(samtools view -c "${BAM_FILE}")
      echo "  Final read count (${sample}): ${READS}"
    fi
  fi

done < "${SAMPLES_TSV}"

echo ""
echo "Step 1 complete. Output: ${DIR}/${DIR_BAM}/"
