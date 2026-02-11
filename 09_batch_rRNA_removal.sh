#!/usr/bin/env bash
set -euo pipefail

# =====================
# Batch rRNA removal for Smart-seq5 outputs
# =====================
# Uses run_rRNA_removal.sbatch and supports BAMs stored in <outdir>/bam/
# Skips samples with missing inputs and records errors to a log

SMARTSEQ_DIR=${SMARTSEQ_DIR:-"$HOME/AP_project/data/smartseq5/smartseq-5-modified"}
OUTROOT=${OUTROOT:-"$HOME/AP_project/data/smartseq5/output_downstream"}
REF_DIR=${REF_DIR:-"$HOME/AP_project/ref/mm10"}
CORES=${CORES:-32}
PARTITION=${PARTITION:-"all_gpu"}
QOS=${QOS:-""}
GRES=${GRES:-""}  # e.g., gpu:1 if required by your cluster

CHROM="${REF_DIR}/chrom.sizes"
GTF="${REF_DIR}/Mus_musculus.GRCm38.102.gtf"
RRNA="${REF_DIR}/rRNA.bed"

SAMPLE_LIST=${SAMPLE_LIST:-""}  # file with one sample id per line
SAMPLE_GLOB=${SAMPLE_GLOB:-"${OUTROOT}/SRR*"}
ERROR_LOG=${ERROR_LOG:-"${OUTROOT}/batch_errors_rRNA.log"}

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
  mv "$outdir"/QC_* "$outdir"/qc/ 2>/dev/null || true
  mv "$outdir"/*Read_numbers*.txt "$outdir"/qc/ 2>/dev/null || true
  mv "$outdir"/*.idxstat "$outdir"/qc/ 2>/dev/null || true
  mv "$outdir"/*.preseq.txt "$outdir"/qc/ 2>/dev/null || true

  mv "$outdir"/*.Log.out "$outdir"/*.Log.final.out "$outdir"/logs/ 2>/dev/null || true
  mv "$outdir"/*.SJ.out.tab "$outdir"/logs/ 2>/dev/null || true

  mv "$outdir"/*ReadsPerGene*.tab "$outdir"/counts/ 2>/dev/null || true
  mv "$outdir"/*.bam "$outdir"/*.bai "$outdir"/bam/ 2>/dev/null || true
  mv "$outdir"/*.bw "$outdir"/bw/ 2>/dev/null || true
  mv "$outdir"/*.bdg "$outdir"/bdg/ 2>/dev/null || true
  mv "$outdir"/*.fastq.gz "$outdir"/fastq/ 2>/dev/null || true
  shopt -u nullglob
}

check_required_bams() {
  local bam_dir=$1
  local id=$2
  local required=(
    "$id.FivePrime_BothStrands.Primary.Multimapping.Read2.sorted.bam"
    "$id.FivePrime_BothStrands.Primary.Unique.Read2.sorted.bam"
    "$id.Else.unique.sorted.bam"
    "$id.Else.all.multimapping.sorted.bam"
    "$id.All.unique.sorted.bam"
    "$id.All.all.multimapping.sorted.bam"
    "$id.FivePrime_BothStrands.Primary.Multimapping.BothStrands.sorted.bam"
    "$id.FivePrime_BothStrands.Primary.Unique.BothStrands.sorted.bam"
  )

  for f in "${required[@]}"; do
    if [[ ! -f "$bam_dir/$f" ]]; then
      return 1
    fi
  done
  return 0
}

# Preflight
for f in "$CHROM" "$GTF" "$RRNA"; do
  if [[ ! -f "$f" ]]; then
    echo "ERROR: missing required file: $f" >&2
    exit 1
  fi
done

mkdir -p "$OUTROOT"
echo "Batch rRNA removal started at $(date)" > "$ERROR_LOG"

echo "Using partition: ${PARTITION:-<none>}  qos: ${QOS:-<none>}  gres: ${GRES:-<none>}"

# Build sample list
sample_ids=()
if [[ -n "$SAMPLE_LIST" && -f "$SAMPLE_LIST" ]]; then
  while IFS= read -r id; do
    [[ -z "$id" ]] && continue
    sample_ids+=("$id")
  done < "$SAMPLE_LIST"
else
  for d in $SAMPLE_GLOB; do
    [[ -d "$d" ]] || continue
    sample_ids+=("$(basename "$d")")
  done
fi

if [[ ${#sample_ids[@]} -eq 0 ]]; then
  echo "ERROR: no samples found under $OUTROOT" >&2
  exit 1
fi

# Run sequentially with progress
current=0
total=${#sample_ids[@]}
for id in "${sample_ids[@]}"; do
  current=$((current + 1))
  progress_bar "$current" "$total" "$id"
  echo

  outdir="$OUTROOT/$id"
  if [[ ! -d "$outdir" ]]; then
    log_error "MISSING_OUTDIR $id $outdir"
    continue
  fi

  bam_dir="$outdir"
  if [[ -d "$outdir/bam" ]]; then
    bam_dir="$outdir/bam"
  fi

  if ! check_required_bams "$bam_dir" "$id"; then
    echo "ERROR: missing BAMs for $id (expected in $bam_dir)" >&2
    log_error "MISSING_BAM $id bam_dir=$bam_dir"
    continue
  fi

  mkdir -p "$outdir/logs"

  sbatch_cmd=(sbatch --wait -o "$outdir/logs/slurm_%j.out" -e "$outdir/logs/slurm_%j.err")
  if [[ -n "$PARTITION" ]]; then
    sbatch_cmd+=(-p "$PARTITION")
  fi
  if [[ -n "$QOS" ]]; then
    sbatch_cmd+=(-q "$QOS")
  fi
  if [[ -n "$GRES" ]]; then
    sbatch_cmd+=(--gres "$GRES")
  fi

  if ! sbatch_out=$(BAM_DIR="$bam_dir" "${sbatch_cmd[@]}" \
    "$SMARTSEQ_DIR/run_rRNA_removal.sbatch" \
    "$id" \
    "$outdir" \
    "$CHROM" \
    "$CORES" \
    "$GTF" \
    "$RRNA" ); then
    echo "ERROR: sbatch failed for $id" >&2
    log_error "SBATCH_FAILED $id"
    continue
  fi

  job_id=$(echo "$sbatch_out" | awk '/Submitted batch job/ {print $4}' | tail -n 1)
  if [[ -z "$job_id" ]]; then
    echo "ERROR: could not parse job id for $id" >&2
    log_error "NO_JOB_ID $id sbatch_out=$(echo "$sbatch_out" | tr '\n' ' ')"
    continue
  fi

  # Organize newly generated outputs
  organize_outputs "$outdir"

done

progress_bar "$total" "$total" "done"
echo
echo "rRNA removal completed. Check $ERROR_LOG for any skipped samples."
