# Gustafson_genomic_analysis_pipeline
This pipeline will go through the steps, taking raw Illumina sequencing output reads, and outputting read counts of individual species. It trims off adapters and conducts quality control using Trim Galore, removes the human genome using BWA, and uses a reference database to analyse reads. For this I am inside my home directory under a folder labeled Genomic_analysis, which has subdirectories for each output type, logs, and the scripts used with each step. There are three major steps. 

# Environment Setup
On the Virginia Tech ARC website, open a terminal where you can create an environment. Use your home directory to enter the terminal. 

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

<details>
  <summary>Click to expand code</summary>

```
#!/bin/bash

#SBATCH -t 70:00:00
#SBATCH -p normal_q
#SBATCH -A introtogds
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nicgistafson1@vt.edu
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

<details>
  <summary>Click to expand code</summary>

```
#!/bin/bash

#SBATCH -t 70:00:00
#SBATCH -p normal_q
#SBATCH -A introtogds
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nicgistafson1@vt.edu 
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

# Pavian Shiny App
Download the report.txt file from the kraken2_outputs directory; this is what will be uploaded into the Pavion. 

https://fbreitwieser.shinyapps.io/pavian/

# Using batch jobs

In order to start a job using one of the three steps under the Genomic_analysis directory, type sbatch "script name here" and it will run that step. The three scripts used are all identical to the code nested under each analysis section in this pipeline.








