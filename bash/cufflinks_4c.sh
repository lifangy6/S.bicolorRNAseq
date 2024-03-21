#!/bin/bash

#SBATCH --job-name cufflinks_4c_job
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=40
#SBATCH --time=12:00:00
#SBATCH --output=/scratch/n/nprovart/lifangy6/outputs/cufflinks_4c_output_%j.out

# Load required modules
module load NiaEnv/2019b
module load cufflinks/2.2.1

# Change directory to the submission directory
cd $SLURM_SUBMIT_DIR

# Define output directory and annotation file
OUTPUT_DIR="cufflinks_4c_output"
ANNOTATION_FILE="Sbicolor_annotations.gtf"

# Define sample labels for each group
LABELS="seed_wt,seed_plasma,shoot_wt,shoot_plasma"

# Define paths to BAM files for each group
GROUP1_BAMS="bam_files/seed_wt_1.bam,bam_files/seed_wt_2.bam,bam_files/seed_wt_3.bam"
GROUP2_BAMS="bam_files/seed_plasma_1.bam,bam_files/seed_plasma_2.bam,bam_files/seed_plasma_3.bam"
GROUP3_BAMS="bam_files/shoot_wt_1.bam,bam_files/shoot_wt_2.bam,bam_files/shoot_wt_3.bam"
GROUP4_BAMS="bam_files/shoot_plasma_1.bam,bam_files/shoot_plasma_2.bam,bam_files/shoot_plasma_3.bam"

# Run Cuffdiff analysis
cuffdiff -p 320 -library-type fr-firststrand -o $OUTPUT_DIR -L $LABELS -u $ANNOTATION_FILE \
$GROUP1_BAMS \
$GROUP2_BAMS \
$GROUP3_BAMS \
$GROUP4_BAMS

