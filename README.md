# Gustafson_genomic_analysis_pipeline
Taking illumina reads to analysis

# Environment Setup
On the virginia tech arc website, open a terminal where you can create an environment. Use your home directory to enter the terminal. 

`conda env create -f environment.yml -n gustafson_analysis`

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
#SBATCH -A gustafson_analysis
#SBATCH --mail-type=ALL
#SBATCH --mail-user=###nicgustafson1@vt.edu 
#SBATCH --cpus-per-task=4
#SBATCH --mem=200GB
#SBATCH --output=trim_galore_%j.out
#SBATCH --error=trim_galore_%j.err

#Path to main folder (likely your home directory)
cd /home/nicgustafson1/genomic_analysis

#Set variables for loop

#create an input and output directory for trim_galore samples, set the thread count, and create a log
INPUT_DIR="/home/nicgustafson/genomic_analysis/data"
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
conda activate g4_viruses

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
conda activate g4_viruses

#Parameters
DB="/home/nicgustafson1/genomic_analysis/databases/kraken2/k2_db"
INPUT_DIR="/home/nicgustafson1/genomic_analysis/trim_galore_outputs"
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






