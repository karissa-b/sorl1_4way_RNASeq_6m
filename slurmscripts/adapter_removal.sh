#! /bin/bash
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=12:00:00
#SBATCH --mem=4GB
#SBATCH -o /fast/users/a1211024/sorl1_4way_6m/02_trimmed_data/logs/%x_%j.out
#SBATCH -e /fast/users/a1211024/sorl1_4way_6m/02_trimmed_data/logs/%x_%j.er

#set up email notifications
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=karissa.barthelson@adelaide.edu.au

#load in modules
module load Java/1.8.0_121
module load fastqc/0.11.4
module load AdapterRemoval/2.2.0-foss-2016uofa

# move to the raw reads directory
cd /fast/users/a1211024/sorl1_4way_6m/01_rawData


for FIRSTREAD in *fastq.gz
do 
echo ${FIRSTREAD}

AdapterRemoval --gzip \
--file1 ${FIRSTREAD} \
--output1 /fast/users/a1211024/sorl1_4way_6m/02_trimmed_data/fastq/${FIRSTREAD%,fastq.gz}_T1.fastq.gz \
--trimns \
--trimqualities \
--minquality 20 --minlength 6 --threads 8

done
