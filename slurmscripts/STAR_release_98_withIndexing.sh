#!/bin/bash
#SBATCH -p batch
#SBATCH --time=12:00:00
#SBATCH --mem=16GB
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -o /fast/users/a1211024/sorl1_4way_6m/00_scripts/slurm/%x_%j.out
#SBATCH -e /fast/users/a1211024/sorl1_4way_6m/00_scripts/slurm/%x_%j.err

# set up email notifications
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=karissa.barthelson@adelaide.edu.au

module load Java/1.8.0_121
module load fastqc/0.11.4
module load STAR/2.7.0d-foss-2016b
module load SAMtools/1.9-foss-2016b

GENOMEDIR="/data/biorefs/reference_genomes/ensembl-release-98/danio_rerio/star/"
OUTDIR="/fast/users/a1211024/sorl1_4way_6m/07_star_release98/"

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

cd ${OUTDIR}
for BAM in *.bam
	do
		 samtools index ${BAM} ${BAM}.bai
	
	done
