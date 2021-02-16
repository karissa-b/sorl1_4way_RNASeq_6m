#!/bin/bash

#SBATCH -p cpu 
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=2:00:00
#SBATCH --mem=4GB
#SBATCH -o /fast/users/a1211024/sorl1_4way_6m/01_rawData/logs/%x_%j.out
#SBATCH -e /fast/users/a1211024/sorl1_4way_6m/01_rawData/logs/%x_%j.err

#set up email notifications
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=karissa.barthelson@adelaide.edu.au

#load in modules
module load Java/1.8.0_121
module load fastqc/0.11.4


OUTDIR="/fast/users/a1211024/sorl1_4way_6m/01_rawData/fastQC_raw/"

cd /fast/users/a1211024/sorl1_4way_6m/01_rawData
       
        for firstread in *.fastq.gz
        do
                fastqc -t 8 -f fastq -o ${OUTDIR} $firstread
        done

