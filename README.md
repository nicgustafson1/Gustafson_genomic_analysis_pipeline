# Gustafson_genomic_analysis_pipeline
The goal of this project is to take short read Illumina sequencing data, trim off the adaptors and low quality bases, remove the human genome, and compare each read to the Kraken database. The end outputs will then be taken to count how many reads are from the same species.

# Environment Setup
On the virginia tech arc website, open a terminal where you can create an environment. Use your home directory to enter the terminal. 

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

The following steps will take you from Illumina output files to a read count that can be read through software. 

# Trim Galore

Trim Galore is used to clean high-throughput sequencing reads by automatically trimming adapters and low-quality bases. It serves as a wrapper around Cutadapt and FastQC, combining adapter removal with quality control checks in a single step. By removing unwanted sequences and short or poor-quality reads, Trim Galore improves the overall accuracy and reliability of downstream analyses such as read alignment and assembly. It takes a fastq file that comes directly from Illumina as an input. 

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
module load miniconda
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

Burrow-Wheeler Aligner for short-read alignment. This maps DNA sequences against a large reference genome, such as the human genome. This uses 3 algorithms- BWA-backtrack, BWA-SW and BWA-MEM. The first algorithm is designed for Illumina sequence reads up to 100bp, while the other two for longer sequences ranged from 70bp to a few megabases.

if you get an error about UNIX line endingsm use this: sed -i 's/\r$//' step3_bwa.sh then run the code again with sbatch 

<details>
  <summary>Click to expand code</summary>

```
#!/bin/bash

#SBATCH -t 70:00:00
#SBATCH -p normal_q
#SBATCH -A gustafson_analysis
#SBATCH --mail-type=ALL
#SBATCH --mail-user=###nicgustafson1@vt.edu
#SBATCH --cpus-per-task=4
#SBATCH --mem=200GB
#SBATCH --output=bwa_%j.out
#SBATCH --error=bwa_%j.err

#Path to main folder (likely your home directory)
cd /home/nicgustafson1/genomic_analysis

#Set Conda Environment
source ~/.bashrc
conda activate gustafson_analysis

#Create an input and output directory for BWA samples, set the thread count, set reference database directory, and create a log
REF="/home/nicgustafson1/genomic_analysis/databases/bwa"
INPUT_DIR="/home/nicgustafson1/genomic_analysis/trim_galore_outputs"
OUTPUT_DIR="/home/nicgustafson1/genomic_analysis/trim_galore_outputs"
LOG_DIR="logs"
THREADS=16

#have log set exact date and time for each iteration
LOGFILE="$LOG_DIR/bwa_${SLURM_JOB_ID}.log"
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOGFILE"; }

log "Starting BWA filtering on $(hostname)"
log "Reference: $REF"
mkdir -p "$OUTPUT_DIR"

#Input previous trim_galore output files and run a loop using BWA to create SAM files that will be converted to BAM files then zip them
for FILE in "$INPUT_DIR"/*_trimmed.fq.gz; do
  [ -e "$FILE" ] || { log "No trimmed FASTQ found in $INPUT_DIR"; break; }
  SAMPLE=$(basename "$FILE" _trimmed.fq.gz)
  log "Processing $SAMPLE"

  SAM="${OUTPUT_DIR}/aln-${SAMPLE}.sam"
  SORTED="${OUTPUT_DIR}/aln-${SAMPLE}.sorted.bam"
  NONHOST="${OUTPUT_DIR}/non_host_reads_${SAMPLE}.bam"
  CLEANED="${OUTPUT_DIR}/cleaned_reads_${SAMPLE}.fastq"

  bwa mem "$REF" "$FILE" > "$SAM"
  samtools view -@ "$THREADS" -Sb "$SAM" | samtools sort -@ "$THREADS" -o "$SORTED"
  samtools index "$SORTED"
  samtools view -@ "$THREADS" -b -f 4 "$SORTED" > "$NONHOST"
  samtools fastq -@ "$THREADS" -0 "$CLEANED" "$NONHOST"

  gzip -f "$SAM" "$SORTED" "$SORTED.bai" "$NONHOST" "$CLEANED"
  log "Finished $SAMPLE"
done

log "BWA complete for all samples."
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
#SBATCH -A gustafson_analysis
#SBATCH --mail-type=ALL
#SBATCH --mail-user=###nicgistafson1@vt.edu 
#SBATCH --cpus-per-task=4
#SBATCH --mem=200GB
#SBATCH --output=kraken2_%j.out
#SBATCH --error=kraken2_%j.err

#Set downloaded directory 
cd /home/nicgustafson1/genomic_analysis

#Set conda environment
source ~/.bashrc
conda activate gustafson_analysis

#Parameters
DB="/home/nicgustafson1/genomic_analysis/databases/kraken2/k2_db"
INPUT_DIR="/home/nicgustafson1/genomic_analysis/bwa_outputs"
OUTPUT_BASE="/home/nicgustafson1/genomic_analysis/kraken2_outputs"
LOG_DIR="logs"
THREADS=16

#Logging setup
#have log set exact date and time for each iteration
LOGFILE="$LOG_DIR/kraken2_${SLURM_JOB_ID:-manual}.log"
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOGFILE"; }

log "Starting Kraken2 classification job on $(hostname)"
log "Using database: $DB"
log "Scanning SPAdes assemblies in: $SPADES_DIR"

#Main loop
for CONTIG_PATH in "${INPUT_DIR}"/sample*_test_data/contigs.fasta; do
    #Skip if no files found
    [ -e "$CONTIG_PATH" ] || { log "No contigs.fasta files found in $SPADES_DIR"; break; }

    #Extract sample name (e.g. sample1_test_data)
    SAMPLE_DIR=$(basename "$(dirname "$CONTIG_PATH")")
    SAMPLE="${SAMPLE_DIR%%_test_data}"

    log "Processing sample: $SAMPLE"

    #Define per-sample output directory under kraken2_outputs
    OUT_DIR="${OUTPUT_BASE}/${SAMPLE_DIR}"
    mkdir -p "$OUT_DIR"

    #Define output file paths
    REPORT="${OUT_DIR}/${SAMPLE}_assembly_report_test_data.txt"
    OUTPUT="${OUT_DIR}/${SAMPLE}_assembly_kraken_test_data.out"
    CLASSIFIED="${OUT_DIR}/${SAMPLE}_assembly_classified_test_data.fastq"

    #Run Kraken2 classification
    k2 classify \
        --db "$DB" \
        "$CONTIG_PATH" \
        --threads "$THREADS" \
        --report "$REPORT" \
        --output "$OUTPUT" \
        --classified-out "$CLASSIFIED" \
        2>&1 | tee -a "$LOGFILE"

    #Compress large outputs
    gzip -f "$OUTPUT" "$CLASSIFIED"

    log "Finished processing $SAMPLE"
    log "--------------------------------"
done

log "All samples processed successfully."
```

</details>

# Pavian Shiny App
This website will be used ot create a graphical representation of read counts. 

https://fbreitwieser.shinyapps.io/pavian/



