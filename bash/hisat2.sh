#!/bin/bash

#SBATCH --job-name hisat2_job
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=40
#SBATCH --time=12:00:00
#SBATCH --output=/scratch/n/nprovart/lifangy6/outputs/hisat2_output_%j.out

# Loading modules needed
module load NiaEnv/2018a
module load gcc/7.3.0
module load hisat2/2.1.0

# Define paths
FASTQ_DIR="$SLURM_SUBMIT_DIR/fastq_files"
OUTPUT_DIR="$SLURM_SUBMIT_DIR/sam_files"
BAM_DIR="$SLURM_SUBMIT_DIR/bam_files"

# Create output directories if they don't exist
mkdir -p "$OUTPUT_DIR"
mkdir -p "$BAM_DIR"

# Change directory to where the script is submitted
cd "$SLURM_SUBMIT_DIR"

# 1. Convert GTF to HISAT2's known-splice-site format
hisat2_extract_splice_sites.py Sbicolor_annotations.gtf > splice_sites.txt
hisat2_extract_exons.py Sbicolor_annotations.gtf > exons.txt

# 2. Build the HISAT2 index
hisat2-build -p 8 --ss splice_sites.txt --exon exons.txt Sbicolor_reference.fa Sbicolor_index

# 3. Align RNA-seq reads using HISAT2 (paired-end)

# Loop through FASTQ files in the directory
for fq_file in "$FASTQ_DIR"/*_1.fq; do
    # Extract file name without extension and remove the "_1" suffix
    file_name=$(basename "$fq_file" _1.fq)
    # Align reads and save output SAM file to the output directory
    hisat2 -p 8 --phred33 -x Sbicolor_index -1 "$fq_file" -2 "${fq_file/_1/_2}" -S "$OUTPUT_DIR/${file_name}.sam"
done

# Now that all SAM files are generated, load Samtools module
module load NiaEnv/2019b
module load samtools/1.13

# 4. Convert SAM files to BAM files and sort them
for sam_file in "$OUTPUT_DIR"/*.sam; do
    # Extract file name without extension
    file_name=$(basename "$sam_file" .sam)
    # Convert SAM to BAM
    samtools view -@ 8 -b -o "$BAM_DIR/${file_name}_unsorted.bam" "$sam_file"
    # Sort BAM file
    samtools sort -@ 8 -o "$BAM_DIR/${file_name}.bam" "$BAM_DIR/${file_name}_unsorted.bam"
    # Optional: Remove the unsorted BAM file to save space
    rm "$BAM_DIR/${file_name}_unsorted.bam"
done
