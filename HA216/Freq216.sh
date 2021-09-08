#!/bin/bash
## This script analyzes position 216 from H9 HA. NOTE: the discrepancy in the HA numbering 
#comes from the difference between H9 (216) and H3 (226)
## The inputs of the script are a reference sequence and fastq.gz files
## This script works on Sapelo2 - for reference -> https://wiki.gacrc.uga.edu/
## Developed by Ginger Giger.
# Ginger Geiger imginger@uga.edu tested-working sapelo2 1.15.2019  
## This version was modified by Lucas M. Ferreri

#!/bin/bash

#SBATCH -J GATK-LoFreq
#SBATCH --time=48:00:00
#SBATCH --partition=highmem_p
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=100000MB

cd $SLURM_SUBMIT_DIR

NAME1=$(find . -maxdepth 1 -name "*_S*_L001_R1_001.fastq.gz" | xargs -I {} basename {} .fastq.gz)

NAME2=$(find . -maxdepth 1 -name "*_S*_L001_R2_001.fastq.gz" | xargs -I {} basename {} .fastq.gz)

NAME3=$(echo *.fastq.gz | awk -v FS='_' '{print $1}')

module load BBMap/38.90-GCC-8.3.0
module load BWA/0.7.17-GCC-8.3.0
module load cutadapt/2.8-GCCcore-8.3.0-Python-3.7.4
module load SAMtools/1.10-iccifort-2019.5.281
module load  R/4.0.0-foss-2019b

#filter reads below Q30
bbduk.sh in1=${NAME1}.fastq.gz in2=${NAME2}.fastq.gz out1=${NAME1}_clean.fastq.gz out2=${NAME2}_clean.fastq.gz qtrim=rl trimq=30 

#filter reads that do not map to the HA
bwa index HA_WF10.fasta
bwa mem -a -t 12 HA_WF10.fasta ${NAME1}_clean.fastq.gz ${NAME2}_clean.fastq.gz > ${NAME3}_bwa.sam 
samtools view -S -b ${NAME3}_bwa.sam > ${NAME3}_bwa.bam

#sort and index the .bam file
samtools sort ${NAME3}_bwa.bam -o ${NAME3}-sorted.bam
rm ${NAME3}_bwa.bam
samtools index ${NAME3}-sorted.bam

#keep only reads mapping to the 27 nucleotides surrounding p226
samtools view -H ${NAME3}-sorted.bam
samtools view -b -F 4 ${NAME3}-sorted.bam "4|HA_WF10:720-746" -o ${NAME3}_bwaKEEP.bam
samtools fastq ${NAME3}_bwaKEEP.bam > ${NAME3}_bwaKEEP.fastq

#cut the reads at the region flanking the 27 nucleotides surrounding p226
cutadapt -g CAGTGATAGGGCCAAGGCCC...GATTATTATTGGTCGGTACT -o ${NAME3}_15-27KEEP.fastq ${NAME3}_bwaKEEP.fastq -m 15 -M 27 \
--untrimmed-output=${NAME3}_untrimmedKEEP.fastq

#find the reads you need
grep "CTTGTCAA" ${NAME3}_15-27KEEP.fastq > ${NAME3}_15-27KEEP.txt

#R CMD BATCH Rvar216trunc_v2.R
# The script Rvar226trunc.R must be placed in ~/scripts
Rscript ~/scripts/Rvar226trunc.R

mv p226freq.tsv ${NAME3}_p226freq.tsv


##Add name of the sample
awk -F, -v s="$NAME3" '{$5=s; print }' ${NAME3}_p226freq.tsv > ${NAME3}_p226freq.temp.tsv

##Delete first column 

cut -f2- ${NAME3}_p226freq.temp.tsv > ${NAME3}_p226freq.tempII.tsv

##Delete header
awk '{if (NR!=1) {print}}' ${NAME3}_p226freq.tempII.tsv > ${NAME3}_p226freq.final.tsv



rm *temp*

rm *.fasta.*

mkdir temp

mv *KEEP* ./temp
mv *bam* ./temp
mv *clean* ./temp
mv *sam ./temp
