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

# load in modules
module load kallisto/0.43.1-foss-2017a

# make the kallisto index file only ran on first go. 
# cd /fast/users/a1211024/sorl1_4way_6m
# kallisto index -i /fast/users/a1211024/sorl1_4way_6m/kallisto_index/Danio_rerio.GRCz11.96.cdna.primary_assembly.with_unspliced_genes_sorll1_R122Pfs_V1482Afs.idx \
# Danio_rerio.GRCz11.96.cdna.primary_assembly.with_unspliced_genes_sorll1_R122Pfs_V1482Afs.fa.gz

# change directory to the trimmed fq files 
cd /fast/users/a1211024/sorl1_4way_6m/02_trimmed_data/fastq

# specify the output directory
outDir=/fast/users/a1211024/sorl1_4way_6m/06_kallisto_with_sorl1_transcripts_mean_300

# Run the kallistoQuant command
# --fr-stranded is strand specific mode, where the read aligns to the F strand from 
# https://www.nugen.com/sites/default/files/M01442v2_User_Guide%3A_Universal_Plus_mRNA-Seq_4274.pdf

	for FIRSTREAD in *fastq.gz
	do
kallisto quant -i /fast/users/a1211024/sorl1_4way_6m/kallisto_index/Danio_rerio.GRCz11.96.cdna.primary_assembly.with_unspliced_genes_sorll1_R122Pfs_V1482Afs.idx \
 -o ${outDir}/${FIRSTREAD%.fastq.gz_T1.fastq.gz} \
 --single \
 -l 300 \
 -s 60 \
 -t 16 \
 --fr-stranded \
 -b 50 \
  ${FIRSTREAD}

done

