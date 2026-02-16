# Gustafson_genomic_analysis_pipeline
This pipeline will go through the steps, taking raw Illumina sequencing output reads, and outputting read counts of individual species. It trims off adapters and conducts quality control using Trim Galore, removes the human genome using BWA, and uses a reference database to analyse reads. For this I am inside my home directory under a folder labeled Genomic_analysis, which has subdirectories for each output type, logs, and the scripts used with each step. There are three major steps. 

# Environment Setup
On the Virginia Tech ARC website, navigate to your home directory. This is where it is easiest to open a terminal from and be in the correct place. 

`interact -A introtogds -p normal_q -t 1:00:00`

# Software installation 

Software was downloaded via Conda. All of the required packages needed to run the pipeline are listed within the environment.yml file in the repository/ directory. 

Initialize Conda on ARC:
```
module load Miniconda3
```
To create the Conda environment:
```
conda env create -f environment.yml -n gustafson_analysis
```

If already created: 
```
conda activate gustafson_analysis
```

# Analysis Pipeline



# Trim Galore

Trim Galore is used to clean high-throughput sequencing reads by automatically trimming adapters and low-quality bases. It serves as a wrapper around Cutadapt and FastQC, combining adapter removal with quality control checks in a single step. By removing unwanted sequences and short or poor-quality reads, Trim Galore improves the overall accuracy and reliability of downstream analyses such as read alignment and assembly. It takes a fastq file that comes directly from Illumina as an input. 

`sbatch step1_trim_galore.sh`

<details>
  <summary>Click to expand code</summary>
  
```
#!/bin/bash

#SBATCH -t 70:00:00
#SBATCH -p normal_q
#SBATCH -A introtogds
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nicgustafson1@vt.edu 
#SBATCH --cpus-per-task=4
#SBATCH --mem=32GB
#SBATCH --output=trim_galore_%j.out
#SBATCH --error=trim_galore_%j.err

#Path to main folder (likely your home directory)
cd /home/nicgustafson1/genomic_analysis

#Set variables for loop

#create an input and output directory for trim_galore samples, set the thread count, and create a log
INPUT_DIR="/home/nicgustafson1/genomic_analysis/data"
OUTPUT_DIR="/home/nicgustafson1/genomic_analysis/trim_galore_outputs"
LOG_DIR="logs"
THREADS=8

mkdir -p "$OUTPUT_DIR" "$LOG_DIR"

#Logging function
LOGFILE="$LOG_DIR/trim_galore_${SLURM_JOB_ID}.log"
#have log set exact date and time for each iteration
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOGFILE"; }

log "Starting Trim Galore job on $(hostname)"
log "Input directory: $INPUT_DIR"
log "Output directory: $OUTPUT_DIR"

#Activate conda environment
source ~/.bashrc
module load Miniconda3
conda activate gustafson_analysis

#Main loop
#Input test data files and run trim_galore on them, outputting them to a new directory
FASTQ_FILES=("$INPUT_DIR"/sample*_test_data.fastq.gz)
[ ${#FASTQ_FILES[@]} -gt 0 ] || { log "No FASTQ files found in $INPUT_DIR"; exit 1; }

for FILE in "${FASTQ_FILES[@]}"; do
    SAMPLE="${FILE%%_test_data.fastq.gz}"
    log "Processing sample: $SAMPLE"
    
    trim_galore "$FILE" -j "$THREADS" -o "$OUTPUT_DIR"
    
    log "Finished sample: $SAMPLE"
done

log "Trim Galore completed successfully."
```
</details>

## BWA 

Burrow-Wheeler Aligner for short-read alignment. This maps DNA sequences against a large reference genome, in this case the human genome. This uses 3 algorithms- BWA-backtrack, BWA-SW and BWA-MEM. The first algorithm is designed for Illumina sequence reads up to 100bp, while the other two for longer sequences ranged from 70bp to a few megabases.

if you get an error about UNIX line endingsm use this: sed -i 's/\r$//' step2_bwa.sh then run the code again with sbatch 

`sbatch step2_bwa.sh`

<details>
  <summary>Click to expand code</summary>

```
#!/bin/bash

#SBATCH -t 70:00:00
#SBATCH -p normal_q
#SBATCH -A introtogds
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nicgustafson1@vt.edu
#SBATCH --cpus-per-task=16
#SBATCH --mem=64GB
#SBATCH --output=bwa_%j.out
#SBATCH --error=bwa_%j.err

# ---------------------------------
# Move to project directory
# ---------------------------------
cd /home/nicgustafson1/genomic_analysis || exit 1

# ---------------------------------
# Load Conda environment safely
# ---------------------------------
module load Miniconda3
source /apps/common/software/Miniconda3/24.7.1-0/etc/profile.d/conda.sh
conda activate gustafson_analysis

# Confirm tools are available
which bwa || exit 1
which samtools || exit 1

# ---------------------------------
# Parameters
# ---------------------------------
REF="/home/nicgustafson1/genomic_analysis/databases/bwa/human_ref.fna"

INPUT_DIR="/home/nicgustafson1/genomic_analysis/trim_galore_outputs"
OUTPUT_DIR="/home/nicgustafson1/genomic_analysis/bwa_outputs"
LOG_DIR="logs"
THREADS=$SLURM_CPUS_PER_TASK

mkdir -p "$OUTPUT_DIR" "$LOG_DIR"

# ---------------------------------
# Logging function
# ---------------------------------
LOGFILE="$LOG_DIR/bwa_${SLURM_JOB_ID}.log"
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOGFILE"; }

log "--------------------------------"
log "Starting BWA alignment on $(hostname)"
log "Input directory: $INPUT_DIR"
log "Output directory: $OUTPUT_DIR"
log "Reference: $REF"
log "Threads: $THREADS"
log "--------------------------------"

# ---------------------------------
# Loop through trimmed FASTQ files
# ---------------------------------
for FILE in "$INPUT_DIR"/sample*_test_data_trimmed.fq.gz; do
    [ -e "$FILE" ] || { log "No trimmed FASTQ files found in $INPUT_DIR"; exit 1; }

    SAMPLE=$(basename "$FILE" "_trimmed.fq.gz")
    log "Processing sample: $SAMPLE"

    SAM="${OUTPUT_DIR}/aln-${SAMPLE}.sam"
    SORTED="${OUTPUT_DIR}/aln-${SAMPLE}.sorted.bam"
    NONHOST="${OUTPUT_DIR}/non_host_${SAMPLE}.bam"
    FASTQ="${OUTPUT_DIR}/${SAMPLE}_nonhost.fastq"

    # ----------------------------
    # Align with BWA (single-end)
    # ----------------------------
    bwa mem -t "$THREADS" "$REF" "$FILE" > "$SAM"

    # ----------------------------
    # Check if SAM contains reads
    # ----------------------------
    READ_COUNT=$(samtools view -c "$SAM" 2>/dev/null)

    if [ "$READ_COUNT" -eq 0 ]; then
        log "WARNING: No reads aligned for $SAMPLE. Sending trimmed reads directly to Kraken."
        cp "$FILE" "${FASTQ}.gz"
        continue
    fi

    # ----------------------------
    # SAM → Sorted BAM
    # ----------------------------
    samtools view -@ "$THREADS" -Sb "$SAM" | samtools sort -@ "$THREADS" -o "$SORTED"
    samtools index "$SORTED"

    # ----------------------------
    # Extract unmapped (non-host) reads
    # ----------------------------
    samtools view -@ "$THREADS" -b -f 4 "$SORTED" > "$NONHOST"

    NONHOST_COUNT=$(samtools view -c "$NONHOST")

    if [ "$NONHOST_COUNT" -eq 0 ]; then
        log "WARNING: No non-host reads for $SAMPLE. Sending all trimmed reads to Kraken."
        cp "$FILE" "${FASTQ}.gz"
    else
        samtools fastq -@ "$THREADS" "$NONHOST" > "$FASTQ"
        gzip -f "$FASTQ"
    fi

    # ----------------------------
    # Compress intermediate files
    # ----------------------------
    gzip -f "$SAM" "$SORTED" "$NONHOST"

    log "Finished processing $SAMPLE"
    log "--------------------------------"
done

log "BWA alignment and FASTQ extraction complete for all samples."
log "--------------------------------"
```

</details>

# Kraken2

Kraken2 is a very fast way to assign taxonomic labels using k-mers to metagenomic DNA sequences. Kraken2 splits sequences into smaller fragments of DNA as "k-mers". The k-mers are then compared in a hashing table to determine similarity to reference genomes in the database. It is used for genomic reads, not protein like Diamond does. In this pipeline, the goal is not to align reads with spades because we want a direct read couunt, so I am skipping over SPAdes. 

`sbatch step3_kraken2.sh`

<details>
  <summary>Click to expand code</summary>

```
#!/bin/bash

#SBATCH -t 70:00:00
#SBATCH -p normal_q
#SBATCH -A introtogds
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nicgustafson1@vt.edu 
#SBATCH --cpus-per-task=4
#SBATCH --mem=100GB
#SBATCH --output=kraken2_%j.out
#SBATCH --error=kraken2_%j.err

# ---------------------------------
# Set working directory
# ---------------------------------
cd /home/nicgustafson1/genomic_analysis

# ---------------------------------
# Load conda environment
# ---------------------------------
source ~/.bashrc
module load Miniconda3
conda activate gustafson_analysis

# ---------------------------------
# Parameters
# ---------------------------------
DB="/home/nicgustafson1/genomic_analysis/databases/kraken2/k2_db"
INPUT_DIR="/home/nicgustafson1/genomic_analysis/bwa_outputs"
OUTPUT_BASE="/home/nicgustafson1/genomic_analysis/kraken2_outputs"
LOG_DIR="logs"

# Use the number of CPUs allocated by SLURM
THREADS=$SLURM_CPUS_PER_TASK

# Create log directory if it doesn't exist
mkdir -p "$LOG_DIR"
LOGFILE="$LOG_DIR/kraken2_${SLURM_JOB_ID:-manual}.log"

# Logging function
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOGFILE"; }

log "--------------------------------"
log "Starting Kraken2 classification job on $(hostname)"
log "Using database: $DB"
log "Scanning input assemblies in: $INPUT_DIR"
log "Threads allocated: $THREADS"
log "--------------------------------"

# ---------------------------------
# Main loop over input files
# ---------------------------------
for CONTIG_PATH in "${INPUT_DIR}"/sample*_test_data_nonhost.fastq.gz; do
    # Skip if no files found
    [ -e "$CONTIG_PATH" ] || { log "No input files found in $INPUT_DIR"; break; }

    # Extract sample name (e.g., sample1, sample2)
    SAMPLE=$(basename "$CONTIG_PATH" | sed 's/_test_data_nonhost.fastq.gz//')
    log "Processing sample: $SAMPLE"

    # Define per-sample output directory
    OUT_DIR="${OUTPUT_BASE}/${SAMPLE}"
    mkdir -p "$OUT_DIR"

    # Define Kraken2 output files
    REPORT="${OUT_DIR}/${SAMPLE}_assembly_report.txt"
    OUTPUT="${OUT_DIR}/${SAMPLE}_assembly_kraken.out"
    CLASSIFIED="${OUT_DIR}/${SAMPLE}_assembly_classified.fastq"

    # Run Kraken2 classification
    k2 classify \
        --db "$DB" \
        "$CONTIG_PATH" \
        --threads "$THREADS" \
        --report "$REPORT" \
        --output "$OUTPUT" \
        --classified-out "$CLASSIFIED" \
        2>&1 | tee -a "$LOGFILE"

    # Compress large output files safely
    [ -f "$OUTPUT" ] && gzip -f "$OUTPUT"
    [ -f "$CLASSIFIED" ] && gzip -f "$CLASSIFIED"

    log "Finished processing $SAMPLE"
    log "--------------------------------"
done

log "All samples processed successfully."
log "--------------------------------"
```

</details>

# Map Tenebrio Molitor

`sbatch step5_tenmol_mapping.sh`

<details>
  <summary>Click to expand code</summary>

```
#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH -p normal_q
#SBATCH -A introtogds
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nicgustafson1@vt.edu
#SBATCH --cpus-per-task=16
#SBATCH --mem=64GB
#SBATCH --output=midori2_blast_%j.out
#SBATCH --error=midori2_blast_%j.err

set -euo pipefail

# ---------------------------------
# Move to project directory
# ---------------------------------
cd /home/nicgustafson1/genomic_analysis || exit 1

# ---------------------------------
# Load tools
# ---------------------------------
module load Miniconda3
source /apps/common/software/Miniconda3/24.7.1-0/etc/profile.d/conda.sh
conda activate gustafson_analysis

# Confirm tools are available
which samtools || exit 1

# ---------------------------------
# BLAST (ARC-installed)
# ---------------------------------
export PATH=/common/data/ncbi/BLAST/ncbi-blast-2.17.0+/bin:$PATH
which blastn || exit 1
which makeblastdb || exit 1
blastn -version >/dev/null

# ---------------------------------
# Scratch (user-writable)
# ---------------------------------
SCRATCH_ROOT="/scratch/${USER}"
mkdir -p "${SCRATCH_ROOT}"

# ---------------------------------
# Parameters
# ---------------------------------
BWA_OUT_DIR="/home/nicgustafson1/genomic_analysis/bwa_outputs"

# LONGEST FASTA (read-only in common)
LONGEST_FASTA="/common/data/ncbi/midori2_coi/MIDORI2_LONGEST_NUC_GB269_CO1_BLAST.fasta"

# Build/use LONGEST BLAST DB in scratch (no -parse_seqids; MIDORI headers have long local IDs)
LONGEST_DB_PREFIX="${SCRATCH_ROOT}/midori2_longest/MIDORI2_LONGEST_NUC_GB269_CO1_BLAST"

# Output
OUT_DIR="/home/nicgustafson1/genomic_analysis/midori2_blast_outputs"
mkdir -p "${OUT_DIR}"

# Scratch working space
SCRATCH_BASE="${SCRATCH_ROOT}/midori2_longest_work_all"
mkdir -p "${SCRATCH_BASE}"

THREADS="${SLURM_CPUS_PER_TASK}"

# ---------------------------------
# Build BLAST DB for LONGEST (in scratch) if needed
# IMPORTANT: no -parse_seqids
# ---------------------------------
mkdir -p "$(dirname "${LONGEST_DB_PREFIX}")"

if [[ -f "${LONGEST_DB_PREFIX}.nsq" && -f "${LONGEST_DB_PREFIX}.nin" && -f "${LONGEST_DB_PREFIX}.nhr" ]]; then
  echo "[INFO] Found existing LONGEST BLAST DB in scratch: ${LONGEST_DB_PREFIX}"
else
  echo "[INFO] LONGEST BLAST DB not found in scratch. Building it now (no -parse_seqids)..."
  makeblastdb \
    -in "${LONGEST_FASTA}" \
    -dbtype nucl \
    -out "${LONGEST_DB_PREFIX}" \
    -hash_index
fi

# ---------------------------------
# Helper: summarize BLAST to read counts by species
# We BLAST with salltitles, then parse the last taxonomy field:
# ... ;Genus_XXXXX;Genus_species_TAXID
# Output: count  species_name  taxid
# ---------------------------------
summarize_counts () {
  local blast_tsv="$1"
  local out_counts="$2"

  awk -F'\t' '
    function clean_species(last,   n, parts, taxid, species, i) {
      n = split(last, parts, "_")
      if (n < 3) return last "\tNA"
      taxid = parts[n]
      species = parts[1]
      for (i=2; i<=n-1; i++) species = species "_" parts[i]
      return species "\t" taxid
    }
    {
      q=$1; title=$2; bits=$13

      split(title, hash, "###")
      taxstr = (length(hash) > 1 ? hash[2] : title)

      n = split(taxstr, fields, ";")
      last = fields[n]

      key = clean_species(last)

      if (!(q in best) || bits > best[q]) {
        best[q]=bits
        best_key[q]=key
      }
    }
    END{
      for (q in best_key) cnt[best_key[q]]++
      for (k in cnt) print cnt[k] "\t" k
    }
  ' "${blast_tsv}" | sort -nr > "${out_counts}"
}

# ---------------------------------
# Find inputs: prefer *.sorted.bam.gz if present; otherwise use *.bam and *.sam.gz
# ---------------------------------
shopt -s nullglob

BAM_GZS=( "${BWA_OUT_DIR}"/*.sorted.bam.gz "${BWA_OUT_DIR}"/*/*.sorted.bam.gz )
BAMS=( "${BWA_OUT_DIR}"/*.bam "${BWA_OUT_DIR}"/*/*.bam )
SAM_GZS=( "${BWA_OUT_DIR}"/*.sam.gz "${BWA_OUT_DIR}"/*/*.sam.gz )

if [[ ${#BAM_GZS[@]} -eq 0 && ${#BAMS[@]} -eq 0 && ${#SAM_GZS[@]} -eq 0 ]]; then
  echo "[ERROR] No input files found in ${BWA_OUT_DIR} (or its immediate subfolders)."
  echo "Looked for: *.sorted.bam.gz, *.bam, *.sam.gz"
  exit 1
fi

echo "[INFO] Found inputs:"
echo "  sorted.bam.gz: ${#BAM_GZS[@]}"
echo "  bam:           ${#BAMS[@]}"
echo "  sam.gz:         ${#SAM_GZS[@]}"

# Build a unique list of "samples" based on file basenames (strip extensions)
# Priority: sorted.bam.gz > bam > sam.gz (skip lower-priority if a higher exists)
declare -A seen
INPUTS=()

for f in "${BAM_GZS[@]}"; do
  base=$(basename "$f")
  sample="${base%.sorted.bam.gz}"
  if [[ -z "${seen[$sample]+x}" ]]; then
    seen["$sample"]=1
    INPUTS+=( "$sample|BAMGZ|$f" )
  fi
done

for f in "${BAMS[@]}"; do
  base=$(basename "$f")
  sample="${base%.bam}"
  if [[ -z "${seen[$sample]+x}" ]]; then
    seen["$sample"]=1
    INPUTS+=( "$sample|BAM|$f" )
  fi
done

for f in "${SAM_GZS[@]}"; do
  base=$(basename "$f")
  sample="${base%.sam.gz}"
  if [[ -z "${seen[$sample]+x}" ]]; then
    seen["$sample"]=1
    INPUTS+=( "$sample|SAMGZ|$f" )
  fi
done

echo "[INFO] Will process ${#INPUTS[@]} samples."

# ---------------------------------
# Main loop
# ---------------------------------
for entry in "${INPUTS[@]}"; do
  IFS='|' read -r sample ftype path <<< "${entry}"

  echo "---------------------------------"
  echo "[INFO] Processing sample: ${sample}"
  echo "[INFO] Input type: ${ftype}"
  echo "[INFO] Input path: ${path}"

  SAMPLE_DIR="${OUT_DIR}/${sample}"
  mkdir -p "${SAMPLE_DIR}"

  # Convert input -> FASTA in scratch
  FASTA_OUT="${SCRATCH_BASE}/${sample}.fasta"

  if [[ "${ftype}" == "BAMGZ" ]]; then
    echo "[INFO] Converting BAM.GZ -> name-sorted -> FASTA..."
    # head not used here; process all reads. Keep pipefail on.
    gunzip -c "${path}" \
      | samtools sort -@ "${THREADS}" -n -O BAM - \
      | samtools fasta -@ "${THREADS}" - \
      > "${FASTA_OUT}"

  elif [[ "${ftype}" == "BAM" ]]; then
    echo "[INFO] Converting BAM -> name-sorted -> FASTA..."
    samtools sort -@ "${THREADS}" -n -O BAM "${path}" \
      | samtools fasta -@ "${THREADS}" - \
      > "${FASTA_OUT}"

  elif [[ "${ftype}" == "SAMGZ" ]]; then
    echo "[INFO] Converting SAM.GZ -> name-sorted -> FASTA..."
    samtools view -@ "${THREADS}" -bh "${path}" \
      | samtools sort -@ "${THREADS}" -n -O BAM - \
      | samtools fasta -@ "${THREADS}" - \
      > "${FASTA_OUT}"
  else
    echo "[ERROR] Unknown input type: ${ftype}"
    exit 1
  fi

  # Sanity check FASTA
  if [[ ! -s "${FASTA_OUT}" ]]; then
    echo "[ERROR] FASTA output is empty for sample ${sample}: ${FASTA_OUT}"
    exit 1
  fi

  # Run blastn
  BLAST_TSV="${SAMPLE_DIR}/${sample}.midori2_longest.blast.tsv"
  echo "[INFO] Running BLAST for ${sample}..."
  blastn \
    -query "${FASTA_OUT}" \
    -db "${LONGEST_DB_PREFIX}" \
    -num_threads "${THREADS}" \
    -task megablast \
    -max_target_seqs 5 \
    -evalue 1e-20 \
    -outfmt "6 qseqid salltitles pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
    -out "${BLAST_TSV}"

  # Summarize to read counts
  COUNTS_TXT="${SAMPLE_DIR}/${sample}.species_read_counts.txt"
  summarize_counts "${BLAST_TSV}" "${COUNTS_TXT}"

  echo "[INFO] Top 10 species for ${sample}:"
  head -n 10 "${COUNTS_TXT}" || true

  # Cleanup scratch FASTA to save space
  rm -f "${FASTA_OUT}"
done

echo "[DONE] Outputs written to: ${OUT_DIR}"
```

</details>

# BLAST

`sbatch step5_tenmol_mapping.sh`

<details>
  <summary>Click to expand code</summary>

```
#!/bin/bash
#SBATCH -t 70:00:00
#SBATCH -p normal_q
#SBATCH -A introtogds
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nicgustafson1@vt.edu
#SBATCH --cpus-per-task=16
#SBATCH --mem=64GB
#SBATCH --output=tenmol_%j.out
#SBATCH --error=tenmol_%j.err

set -euo pipefail

# ---------------------------------
# Move to project directory
# ---------------------------------
cd /home/nicgustafson1/genomic_analysis || exit 1

# ---------------------------------
# Load Conda environment safely
# ---------------------------------
module load Miniconda3
source /apps/common/software/Miniconda3/24.7.1-0/etc/profile.d/conda.sh
conda activate gustafson_analysis

# Confirm tools are available
which bwa || exit 1
which samtools || exit 1
which gzip || exit 1
which zcat || exit 1

# ---------------------------------
# Parameters
# ---------------------------------
TENMOL_REF="/home/nicgustafson1/genomic_analysis/databases/tenmol/tenmol_genome.fna"
INPUT_DIR="/home/nicgustafson1/genomic_analysis/bwa_outputs"
OUTPUT_DIR="/home/nicgustafson1/genomic_analysis/tenmol_mapping_outputs"
LOG_DIR="$OUTPUT_DIR/logs"
THREADS=${SLURM_CPUS_PER_TASK:-16}

mkdir -p "$OUTPUT_DIR" "$LOG_DIR"

# ---------------------------------
# Logging function
# ---------------------------------
LOGFILE="$LOG_DIR/tenmol_map_${SLURM_JOB_ID}.log"
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOGFILE"; }

log "--------------------------------"
log "Starting Tenebrio mapping on $(hostname)"
log "Input directory: $INPUT_DIR"
log "Output directory: $OUTPUT_DIR"
log "Tenebrio reference: $TENMOL_REF"
log "Threads: $THREADS"
log "--------------------------------"

# Ensure reference + BWA index exists (basic check: .bwt)
if [ ! -f "${TENMOL_REF}.bwt" ]; then
  log "ERROR: BWA index not found next to $TENMOL_REF (missing ${TENMOL_REF}.bwt). Index it first with: bwa index $TENMOL_REF"
  exit 1
fi

# Summary TSV (one row per sample)
SUMMARY_TSV="$OUTPUT_DIR/tenmol_mapping_summary.tsv"
echo -e "sample\tinput_bam_gz\ttotal_reads\tmapped_reads\tmapped_percent\ttop_contig\ttop_contig_mapped\tpct_mapped_on_top_contig" > "$SUMMARY_TSV"

# ---------------------------------
# Loop through BWA output BAMs
# ---------------------------------
for BAMGZ in "$INPUT_DIR"/aln-sample*_test_data.sorted.bam.gz; do
  [ -e "$BAMGZ" ] || { log "No files matching $INPUT_DIR/aln-sample*_test_data.sorted.bam.gz"; exit 1; }

  BASE=$(basename "$BAMGZ")
  SAMPLE=${BASE%.sorted.bam.gz}
  SAMPLE=${SAMPLE#aln-}

  log "Processing sample: $SAMPLE"
  log "Input BAM.GZ: $BAMGZ"

  # Work files
  BAM="$OUTPUT_DIR/${SAMPLE}.sorted.bam"
  FASTQ="$OUTPUT_DIR/${SAMPLE}.all.fastq.gz"
  TENBAM="$OUTPUT_DIR/${SAMPLE}.vs_tenmol.sorted.bam"
  FLAGSTAT="$OUTPUT_DIR/${SAMPLE}.vs_tenmol.flagstat.txt"
  IDXSTATS="$OUTPUT_DIR/${SAMPLE}.vs_tenmol.idxstats.txt"

  # Decompress BAM.GZ -> BAM (samtools works reliably on plain .bam)
  log "Decompressing BAM.GZ -> BAM"
  gunzip -c "$BAMGZ" > "$BAM"

  # Sanity check: does BAM contain reads?
  TOTAL_IN_BAM=$(samtools view -c "$BAM" 2>/dev/null || echo 0)
  if [ "$TOTAL_IN_BAM" -eq 0 ]; then
    log "WARNING: $SAMPLE BAM contains 0 reads; skipping mapping."
    echo -e "${SAMPLE}\t${BAMGZ}\t0\t0\t0\tNA\t0\t0" >> "$SUMMARY_TSV"
    rm -f "$BAM"
    continue
  fi

  # BAM -> single-end FASTQ (your BAMs are flagged as unpaired; this avoids empty R1/R2 files)
  log "Converting BAM -> FASTQ (single-end)"
  samtools fastq -@ "$THREADS" "$BAM" | gzip > "$FASTQ"

  # Quick sanity check: FASTQ read count
  FASTQ_READS=$(( $(zcat "$FASTQ" | wc -l) / 4 ))
  if [ "$FASTQ_READS" -eq 0 ]; then
    log "WARNING: $SAMPLE FASTQ is empty after conversion; skipping mapping."
    echo -e "${SAMPLE}\t${BAMGZ}\t${TOTAL_IN_BAM}\t0\t0\tNA\t0\t0" >> "$SUMMARY_TSV"
    rm -f "$BAM" "$FASTQ"
    continue
  fi

  # Map to Tenebrio genome and sort
  log "Mapping reads to Tenebrio reference"
  bwa mem -t "$THREADS" "$TENMOL_REF" "$FASTQ" | samtools sort -@ "$THREADS" -o "$TENBAM"
  samtools index "$TENBAM"

  # Collect stats
  log "Computing flagstat + idxstats"
  samtools flagstat "$TENBAM" > "$FLAGSTAT"
  samtools idxstats "$TENBAM" > "$IDXSTATS"

  # Parse totals + mapped from flagstat (single-end format)
  TOTAL_READS=$(awk 'NR==1{print $1}' "$FLAGSTAT")
  MAPPED_READS=$(awk 'NR==7{print $1}' "$FLAGSTAT")
  MAPPED_PCT=$(awk 'NR==7{gsub(/[()%]/,"",$5); print $5}' "$FLAGSTAT")

  # Top contig from idxstats (by mapped reads, col3)
  TOP_LINE=$(sort -k3,3nr "$IDXSTATS" | head -n 1)
  TOP_CONTIG=$(echo "$TOP_LINE" | awk '{print $1}')
  TOP_MAPPED=$(echo "$TOP_LINE" | awk '{print $3}')

  # Percent of mapped reads that land on top contig
  TOP_PCT="0"
  if [ "$MAPPED_READS" -gt 0 ]; then
    TOP_PCT=$(awk -v top="$TOP_MAPPED" -v mapped="$MAPPED_READS" 'BEGIN{printf "%.4f", (top/mapped)*100}')
  fi

  # Write summary row
  echo -e "${SAMPLE}\t${BAMGZ}\t${TOTAL_READS}\t${MAPPED_READS}\t${MAPPED_PCT}\t${TOP_CONTIG}\t${TOP_MAPPED}\t${TOP_PCT}" >> "$SUMMARY_TSV"

  log "Finished $SAMPLE: total=${TOTAL_READS}, mapped=${MAPPED_READS} (${MAPPED_PCT}%), top_contig=${TOP_CONTIG} top_mapped=${TOP_MAPPED} top_contig_pct_of_mapped=${TOP_PCT}%"
  log "--------------------------------"

  # Optional cleanup of big intermediates (uncomment if you want to save space)
  # rm -f "$BAM" "$FASTQ"
done

log "All samples complete."
log "Summary written to: $SUMMARY_TSV"
log "--------------------------------"
```

</details>

# Pavian Shiny App
Download the report.txt file from the kraken2_outputs directory; this is what will be uploaded into the Pavion. 

https://fbreitwieser.shinyapps.io/pavian/

# Using batch jobs

In order to start a job using one of the three steps under the Genomic_analysis directory, type sbatch "script name here" and it will run that step. The three scripts used are all identical to the code nested under each analysis section in this pipeline.




# Example Output

This image was generated using this pipeline with a kraken2 standard bacteria database:
<img width="1006" height="501" alt="image" src="https://github.com/user-attachments/assets/742cd818-a327-4b15-a126-b8185168926b" />

# Extra stuff
To get to your scratch folder /scratch/"username" for me it is nicgustafson1
To access the current NCBI database do 
cd /common/data/ncbi/fcs-gx/2025-11-11/gxdb

# 10k Reads of Sample 2 testing (this worked) 

`sbatch step5_tenmol_mapping.sh`

<details>
  <summary>Click to expand code</summary>

```
#!/bin/bash
#SBATCH -t 02:00:00
#SBATCH -p normal_q
#SBATCH -A introtogds
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nicgustafson1@vt.edu
#SBATCH --cpus-per-task=16
#SBATCH --mem=64GB
#SBATCH --output=midori2_blast_test_sample2_10k_%j.out
#SBATCH --error=midori2_blast_test_sample2_10k_%j.err

set -euo pipefail

cd /home/nicgustafson1/genomic_analysis || exit 1

module load Miniconda3
source /apps/common/software/Miniconda3/24.7.1-0/etc/profile.d/conda.sh
conda activate gustafson_analysis

which samtools || exit 1

# BLAST
export PATH=/common/data/ncbi/BLAST/ncbi-blast-2.17.0+/bin:$PATH
which blastn || exit 1
which makeblastdb || exit 1

# Scratch
SCRATCH_ROOT="/scratch/${USER}"
mkdir -p "${SCRATCH_ROOT}"

# Inputs
BAM_GZ="/home/nicgustafson1/genomic_analysis/bwa_outputs/aln-sample2_test_data.sorted.bam.gz"
LONGEST_FASTA="/common/data/ncbi/midori2_coi/MIDORI2_LONGEST_NUC_GB269_CO1_BLAST.fasta"

# DB in scratch (no -parse_seqids)
LONGEST_DB_PREFIX="${SCRATCH_ROOT}/midori2_longest/MIDORI2_LONGEST_NUC_GB269_CO1_BLAST"

# Outputs
OUT_DIR="/home/nicgustafson1/genomic_analysis/midori2_blast_outputs_TEST_SAMPLE2_10K"
mkdir -p "${OUT_DIR}"

SCRATCH_BASE="${SCRATCH_ROOT}/midori2_longest_work_test_sample2_10k"
mkdir -p "${SCRATCH_BASE}"

THREADS="${SLURM_CPUS_PER_TASK}"

# Build BLAST DB (no -parse_seqids)
mkdir -p "$(dirname "${LONGEST_DB_PREFIX}")"
if [[ -f "${LONGEST_DB_PREFIX}.nsq" && -f "${LONGEST_DB_PREFIX}.nin" && -f "${LONGEST_DB_PREFIX}.nhr" ]]; then
  echo "[INFO] Found existing LONGEST BLAST DB: ${LONGEST_DB_PREFIX}"
else
  echo "[INFO] Building LONGEST BLAST DB (no -parse_seqids)..."
  makeblastdb -in "${LONGEST_FASTA}" -dbtype nucl -out "${LONGEST_DB_PREFIX}" -hash_index
fi

# Summarize by best hit per read; parse species + taxid from MIDORI title
summarize_counts () {
  local blast_tsv="$1"
  local out_counts="$2"

  awk -F'\t' '
    function clean_species(last,   n, parts, taxid, species, i) {
      n = split(last, parts, "_")
      if (n < 3) return last "\tNA"
      taxid = parts[n]
      species = parts[1]
      for (i=2; i<=n-1; i++) species = species "_" parts[i]
      return species "\t" taxid
    }
    {
      q=$1; title=$2; bits=$13

      split(title, hash, "###")
      taxstr = (length(hash) > 1 ? hash[2] : title)

      n = split(taxstr, fields, ";")
      last = fields[n]

      key = clean_species(last)

      if (!(q in best) || bits > best[q]) {
        best[q]=bits
        best_key[q]=key
      }
    }
    END{
      for (q in best_key) cnt[best_key[q]]++
      for (k in cnt) print cnt[k] "\t" k
    }
  ' "${blast_tsv}" | sort -nr > "${out_counts}"
}

# Run sample2 (first 10k reads)
if [[ ! -f "${BAM_GZ}" ]]; then
  echo "[ERROR] Input not found: ${BAM_GZ}"
  exit 1
fi

sample="aln-sample2_test_data_10k"
SAMPLE_DIR="${OUT_DIR}/${sample}"
mkdir -p "${SAMPLE_DIR}"

FASTA_OUT="${SCRATCH_BASE}/${sample}.fasta"

echo "[INFO] Converting coordinate-sorted BAM.GZ -> name-sorted -> FASTA; taking first 10,000 reads..."
# NOTE: head causes SIGPIPE upstream when it finishes; with pipefail this would fail the job.
# We temporarily disable pipefail for just this pipeline.
set +o pipefail
gunzip -c "${BAM_GZ}" \
  | samtools sort -@ "${THREADS}" -n -O BAM - \
  | samtools fasta -@ "${THREADS}" - \
  | head -n 20000 > "${FASTA_OUT}"
set -o pipefail

# Sanity check: ensure we got 10,000 reads (20,000 FASTA lines)
LINES=$(wc -l < "${FASTA_OUT}" || echo 0)
if [[ "${LINES}" -lt 20000 ]]; then
  echo "[ERROR] FASTA too short (${LINES} lines). Input may have <10,000 reads or conversion failed."
  exit 1
fi

echo "[INFO] Running BLAST..."
BLAST_TSV="${SAMPLE_DIR}/${sample}.midori2_longest.blast.tsv"
blastn \
  -query "${FASTA_OUT}" \
  -db "${LONGEST_DB_PREFIX}" \
  -num_threads "${THREADS}" \
  -task megablast \
  -max_target_seqs 5 \
  -evalue 1e-20 \
  -outfmt "6 qseqid salltitles pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
  -out "${BLAST_TSV}"

COUNTS_TSV="${SAMPLE_DIR}/${sample}.species_read_counts.tsv"
summarize_counts "${BLAST_TSV}" "${COUNTS_TSV}"

echo "[INFO] Top 20 species:"
head -n 20 "${COUNTS_TSV}" || true

rm -f "${FASTA_OUT}"
echo "[DONE] Outputs in: ${SAMPLE_DIR}"
```

</details>
