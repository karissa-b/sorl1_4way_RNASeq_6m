#!/bin/bash

#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=2:00:00
#SBATCH --mem=16GB
#SBATCH -o /fast/users/a1211024/sorl1_4way_6m/00_scripts/slurm/%x_%j.out
#SBATCH -e /fast/users/a1211024/sorl1_4way_6m/00_scripts/slurm/%x_%j.err

#set up email notifications
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=karissa.barthelson@adelaide.edu.au

#load in modules
module load Subread/1.5.2-foss-2016b
module load Java/1.8.0_121

# define paths
BAMDIR="/fast/users/a1211024/sorl1_4way_6m/07_star_release98"
OUTDIR="/fast/users/a1211024/sorl1_4way_6m/08_featureCounts_release98"
GTF="/fast/users/a1211024/sorl1_4way_6m/Danio_rerio.GRCz11.98.gtf"

# move to bam dir
cd ${BAMDIR}
pwd
ls 

featureCounts -a ${GTF} -Q 10 -s 1 --fracOverlap 1 \
	-T 8 \
	-o ${OUTDIR}/featureCounts_fracOverlap1_q10_s1_ensRelease98.txt \
	*.bam

