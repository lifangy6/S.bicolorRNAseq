#!/bin/bash

#SBATCH --job-name cufflinks_job
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=40
#SBATCH --time=12:00:00
#SBATCH --output=/scratch/n/nprovart/lifangy6/outputs/cufflinks_output_%j.out

# Load required modules
module load NiaEnv/2019b
module load cufflinks/2.2.1

# Change directory to the submission directory
cd $SLURM_SUBMIT_DIR

# Define output directory and annotation file
OUTPUT_DIR="cufflinks_output"
ANNOTATION_FILE="Sbicolor_annotations.gtf"

# Define sample labels for each group
LABELS="seed_dry,seed_plasma,seed_wt,shoot_plasma,shoot_wt"

# Define paths to BAM files for each group
GROUP1_BAMS="bam_files/seed_dry_1.bam,bam_files/seed_dry_2.bam,bam_files/seed_dry_3.bam"
GROUP2_BAMS="bam_files/seed_plasma_1.bam,bam_files/seed_plasma_2.bam,bam_files/seed_plasma_3.bam"
GROUP3_BAMS="bam_files/seed_wt_1.bam,bam_files/seed_wt_2.bam,bam_files/seed_wt_3.bam"
GROUP4_BAMS="bam_files/shoot_plasma_1.bam,bam_files/shoot_plasma_2.bam,bam_files/shoot_plasma_3.bam"
GROUP5_BAMS="bam_files/shoot_wt_1.bam,bam_files/shoot_wt_2.bam,bam_files/shoot_wt_3.bam"

# Run Cuffdiff analysis
cuffdiff -p 320 -library-type fr-firststrand -o $OUTPUT_DIR -L $LABELS -u $ANNOTATION_FILE \
$GROUP1_BAMS \
$GROUP2_BAMS \
$GROUP3_BAMS \
$GROUP4_BAMS \
$GROUP5_BAMS

