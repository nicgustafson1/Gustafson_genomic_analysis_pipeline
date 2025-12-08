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

Trim Galore is used to clean high-throughput sequencing reads by automatically trimming adapters and low-quality bases. It serves as a wrapper around Cutadapt and FastQC, combining adapter removal with quality control checks in a single step. By removing unwanted sequences and short or poor-quality reads, Trim Galore improves the overall accuracy and reliability of downstream analyses such as read alignment and assembly. It takes a fastq file that comes directly from Illumina as an imput. 
