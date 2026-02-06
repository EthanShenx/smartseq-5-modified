#!/usr/bin/env bash
set -euo pipefail

# =====================
# User configuration
# =====================
SMARTSEQ_DIR=${SMARTSEQ_DIR:-"$HOME/AP_project/data/smartseq5/smartseq-5-modified"}
DATA_DIR=${DATA_DIR:-"$HOME/AP_project/data/smartseq5"}
OUTROOT=${OUTROOT:-"$HOME/AP_project/data/smartseq5/output_01"}
REF_DIR=${REF_DIR:-"$HOME/AP_project/ref/mm10"}
CORES=${CORES:-32}
DUPLICATE_PIXEL_DIST=${DUPLICATE_PIXEL_DIST:-2500}
PARTITION=${PARTITION:-""}
QOS=${QOS:-""}

ADAP_ILL="${SMARTSEQ_DIR}/adaptors/NexteraPE-PE.fa"
ADAP_SS2="${SMARTSEQ_DIR}/adaptors/SS2_adaptors_RCs.fa"
FASTA="${REF_DIR}/Mus_musculus.GRCm38.dna.primary_assembly.fa"
GTF="${REF_DIR}/Mus_musculus.GRCm38.102.gtf"
CHROM="${REF_DIR}/chrom.sizes"
STARIDX="${REF_DIR}/star_index"

# Optional: if you want to run a subset, set SAMPLE_GLOB or SAMPLE_LIST
SAMPLE_GLOB=${SAMPLE_GLOB:-"${DATA_DIR}/*_1.fastq.gz"}
SAMPLE_LIST=${SAMPLE_LIST:-""}  # file with one sample id per line (e.g., SRR12345)
ERROR_LOG=${ERROR_LOG:-"${OUTROOT}/batch_errors.log"}

# =====================
# Helpers
# =====================
progress_bar() {
  local current=$1
  local total=$2
  local label=${3:-""}
  local width=30
  local filled=$((current * width / total))
  local empty=$((width - filled))
  local bar space
  printf -v bar "%${filled}s" ""
  printf -v space "%${empty}s" ""
  bar=${bar// /#}
  space=${space// /-}
  printf "\r[%s%s] %d/%d %s" "$bar" "$space" "$current" "$total" "$label"
}

log_error() {
  local msg=$1
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] $msg" >> "$ERROR_LOG"
}

organize_outputs() {
  local outdir=$1
  mkdir -p "$outdir"/{logs,fastqc,qc,bam,bw,bdg,counts,fastq,other}
  shopt -s nullglob

  mv "$outdir"/*.html "$outdir"/*.zip "$outdir"/fastqc/ 2>/dev/null || true
  mv "$outdir"/QC_* "$outdir"/qc/ 2>/dev/null || true
  mv "$outdir"/*Read_numbers*.txt "$outdir"/qc/ 2>/dev/null || true
  mv "$outdir"/*.preseq.txt "$outdir"/qc/ 2>/dev/null || true
  mv "$outdir"/*.idxstat "$outdir"/qc/ 2>/dev/null || true

  mv "$outdir"/*.Log.out "$outdir"/*.Log.final.out "$outdir"/*.Log.final.out "$outdir"/logs/ 2>/dev/null || true
  mv "$outdir"/*.SJ.out.tab "$outdir"/logs/ 2>/dev/null || true
  mv "$outdir"/*Log.out "$outdir"/*Log.final.out "$outdir"/logs/ 2>/dev/null || true

  mv "$outdir"/*ReadsPerGene*.tab "$outdir"/counts/ 2>/dev/null || true

  mv "$outdir"/*.bam "$outdir"/*.bai "$outdir"/bam/ 2>/dev/null || true
  mv "$outdir"/*.bw "$outdir"/bw/ 2>/dev/null || true
  mv "$outdir"/*.bdg "$outdir"/bdg/ 2>/dev/null || true

  mv "$outdir"/*.fastq.gz "$outdir"/fastq/ 2>/dev/null || true

  shopt -u nullglob
}

check_errors_in_log() {
  local logfile=$1
  if [[ ! -f "$logfile" ]]; then
    echo "ERROR: log file not found: $logfile" >&2
    return 1
  fi

  # Patterns that indicate real failures (locale warnings are suppressed in scripts)
  local pattern='ERROR!!|FATAL|Fatal|Exception|Traceback|OutOfMemoryError|Segmentation fault|command not found|No such file or directory|EXITING:|ERROR: Invalid|Killed'
  if rg -n "$pattern" "$logfile" >/dev/null 2>&1; then
    echo "ERROR: failure detected in log: $logfile" >&2
    rg -n "$pattern" "$logfile" | head -n 20 >&2
    return 1
  fi
  return 0
}

# =====================
# Preflight checks
# =====================
for f in "$ADAP_ILL" "$ADAP_SS2" "$FASTA" "$GTF" "$CHROM"; do
  if [[ ! -f "$f" ]]; then
    echo "ERROR: missing required file: $f" >&2
    exit 1
  fi
done
if [[ ! -d "$STARIDX" ]]; then
  echo "ERROR: STAR index not found: $STARIDX" >&2
  exit 1
fi

mkdir -p "$OUTROOT"
echo "Batch run started at $(date)" > "$ERROR_LOG"

# Build sample list
sample_ids=()
if [[ -n "$SAMPLE_LIST" && -f "$SAMPLE_LIST" ]]; then
  while IFS= read -r id; do
    [[ -z "$id" ]] && continue
    sample_ids+=("$id")
  done < "$SAMPLE_LIST"
else
  for f1 in $SAMPLE_GLOB; do
    [[ -e "$f1" ]] || continue
    id=$(basename "$f1" _1.fastq.gz)
    sample_ids+=("$id")
  done
fi

if [[ ${#sample_ids[@]} -eq 0 ]]; then
  echo "ERROR: no samples found (check SAMPLE_GLOB or SAMPLE_LIST)" >&2
  exit 1
fi

total=${#sample_ids[@]}
current=0

echo "Found $total samples. Running sequentially with skip-on-error.";

# =====================
# Run samples (sequential with progress bar)
# =====================
for id in "${sample_ids[@]}"; do
  current=$((current + 1))
  progress_bar "$current" "$total" "$id"
  echo

  f1="$DATA_DIR/${id}_1.fastq.gz"
  f2="$DATA_DIR/${id}_2.fastq.gz"
  if [[ ! -f "$f1" || ! -f "$f2" ]]; then
    echo "ERROR: missing FASTQ pair for $id" >&2
    log_error "MISSING_FASTQ $id f1=$f1 f2=$f2"
    continue
  fi

  outdir="$OUTROOT/$id"
  mkdir -p "$outdir/logs"

  # Submit and wait
  sbatch_cmd=(sbatch --wait -o "$outdir/logs/slurm_%j.out" -e "$outdir/logs/slurm_%j.err")
  if [[ -n "$PARTITION" ]]; then
    sbatch_cmd+=(-p "$PARTITION")
  fi
  if [[ -n "$QOS" ]]; then
    sbatch_cmd+=(-q "$QOS")
  fi

  if ! sbatch_out=$("${sbatch_cmd[@]}" \
    "$SMARTSEQ_DIR/runRNASeqProcessingPE.sbatch" \
    "$id" \
    "$outdir" \
    "$f1" \
    "$f2" \
    "$ADAP_ILL" \
    "$ADAP_SS2" \
    "$STARIDX" \
    "$CHROM" \
    "$CORES" \
    "$DUPLICATE_PIXEL_DIST" \
    "$FASTA" \
    "$GTF" \
    "$GTF" ); then
    echo "ERROR: sbatch failed for $id" >&2
    log_error "SBATCH_FAILED $id"
    continue
  fi

  job_id=$(echo "$sbatch_out" | awk '/Submitted batch job/ {print $4}' | tail -n 1)
  if [[ -z "$job_id" ]]; then
    echo "ERROR: could not parse job id for $id" >&2
    echo "$sbatch_out" >&2
    log_error "NO_JOB_ID $id sbatch_out=$(echo "$sbatch_out" | tr '\n' ' ')"
    continue
  fi

  log_file="$outdir/logs/slurm_${job_id}.out"

  # Stop immediately if log has failure patterns
  if ! check_errors_in_log "$log_file"; then
    echo "ERROR: errors detected in $id (job $job_id)." >&2
    log_error "JOB_ERROR $id job=$job_id log=$log_file"
    continue
  fi

  # Organize outputs for this sample
  organize_outputs "$outdir"

done

progress_bar "$total" "$total" "done"
echo
echo "All samples completed. Check $ERROR_LOG for any skipped samples."
