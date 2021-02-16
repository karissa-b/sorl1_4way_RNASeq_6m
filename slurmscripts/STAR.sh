#!/bin/bash
#SBATCH -p batch
#SBATCH --time=12:00:00
#SBATCH --mem=16GB
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -o /fast/users/a1211024/sorl1_4way_6m/04_star/logs/%x_%j.out
#SBATCH -e /fast/users/a1211024/sorl1_4way_6m/04_star/logs/%x_%j.err

# set up email notifications
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=karissa.barthelson@adelaide.edu.au

module load Java/1.8.0_121
module load fastqc/0.11.4
module load STAR/2.5.3a-foss-2016b

GENOMEDIR="/fast/users/a1211024/sorl1_W1818X_6month/referenceFiles/starIndex/"
# FIRSTREAD="*_S`${SLURM_ARRAY_TASK_ID}`_merged_R1_001_T1.fastq.gz"
OUTDIR="/fast/users/a1211024/sorl1_4way_6m/04_star/bams/"

cd /fast/users/a1211024/sorl1_4way_6m/02_trimmed_data/fastq

for FIRSTREAD in *.fastq.gz
        do
		echo "now processing" ${FIRSTREAD}

		STAR --runThreadN 16 \
		--genomeDir ${GENOMEDIR} \
		--readFilesIn ${FIRSTREAD} \
		--outFileNamePrefix ${OUTDIR}${FIRSTREAD%_merged_R1_001_T1.fastq.gz} \
		--outSAMtype BAM SortedByCoordinate \
		--readFilesCommand zcat
	 	
	done

