#!/bin/bash

#SBATCH -p cpu 
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=3:00:00
#SBATCH --mem=4G
#SBATCH -o /fast/users/a1211024/sorl1_4way_6m/00_scripts/slurm/%x_%j.out
#SBATCH -e /fast/users/a1211024/sorl1_4way_6m/00_scripts/slurm/%x_%j.err

#set up email notifications
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=karissa.barthelson@adelaide.edu.au

#load in modules
module load Java/1.8.0_121
module load fastqc/0.11.4
module load SAMtools/1.9-foss-2016b

# change directory to bams
cd /fast/users/a1211024/sorl1_4way_6m/04_star/bams/bams

# define bam files & output
BAM="/fast/users/a1211024/sorl1_4way_6m/04_star/bams/bam/*.bam"
OUTDIR="/fast/users/a1211024/sorl1_4way_6m/04_star/bams"


for BAM in *.bam
do
	samtools index ${BAM} ${BAM}.bai
	fastqc -t 8 bam_mapped -o $OUTDIR/qc $BAM
	samtools flagstat ${BAM}
done





