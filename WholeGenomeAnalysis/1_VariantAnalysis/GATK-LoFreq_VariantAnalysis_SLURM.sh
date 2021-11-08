#!/bin/bash

#SBATCH -J GATK-LoFreq
#SBATCH --time=48:00:00
#SBATCH --partition=highmem_p
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=100000MB

#This script prepares files for Variant Analysis using LoFreq
#Made by Lucas M. Ferreri
#Modified 20201104
cd $SLURM_SUBMIT_DIR

NAME1=$(find . -maxdepth 1 -name "*_S*_L001_R1_001.fastq.gz" | xargs -I {} basename {} .fastq.gz)

NAME2=$(find . -maxdepth 1 -name "*_S*_L001_R2_001.fastq.gz" | xargs -I {} basename {} .fastq.gz)

NAME3=$(echo *.fastq.gz | awk -v FS='_' '{print $1}')

cwd=$(pwd)

module load cutadapt/2.8-GCCcore-8.3.0-Python-3.7.4      

time cutadapt -f fastq --match-read-wildcards \
-e 0.1 -O 6 -m 32 -g ^GAGCTAGTCTG -g ^GTCGAGCTCG -g ^GTTACGCGCC -g ^GGGGGG -g ^CGGGTTATT -g ^GGTAACGCGTGATC \
-a GATCGGAAGAGCACACGTCT -b ACACTCTTTCCCTACACGACGCTCTTCCGATCT -b GCCAGAGCCGTAAGGACGACTTGGCGAGAAGGCTAGA \
-o "$NAME1"_cutadapt.fastq.gz -p "$NAME2"_cutadapt.fastq.gz \
"$NAME1".fastq.gz "$NAME2".fastq.gz 1>job.out 2>job.err \

module unload cutadapt/2.8-GCCcore-8.3.0-Python-3.7.4     

# as per LoFreq recommendation, QC is is not needed since it is handled by the GATK realingment itself.

mkdir original-data

mv *1.fastq.gz ./original-data

cat *.fastq.gz > "$NAME3"_cutadapt.fastq.gz
mv *_001_cutadapt.fastq.gz ./original-data

module load BWA/0.7.17-GCC-8.3.0
#module load samtools/1.3.1
echo `module list`

NAME=$(find . -maxdepth 1 -name "*.fastq.gz" | xargs -I {} basename {} .fastq.gz)

echo "name of the file is : $NAME"

cwd=$(pwd)

bwa index WF10_curated_varpatch.fa

time bwa mem -R '@RG\tID:S1\tSM:SAMPLE1' WF10_curated_varpatch.fa \
"$NAME".fastq.gz > "$NAME".sam

module load SAMtools/1.6-foss-2019b

samtools view -S "$NAME".sam -b -o "$NAME".bam

module load picard/2.16.0-Java-1.8.0_144

echo `module list`

time java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144" -jar \
$EBROOTPICARD/picard.jar \
SortSam I="$NAME".bam \
O="$NAME"_bwa_sorted.bam SORT_ORDER=coordinate

samtools index "$NAME"_bwa_sorted.bam

time java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144" -jar \
$EBROOTPICARD/picard.jar \
AddOrReplaceReadGroups \
I="$NAME"_bwa_sorted.bam \
O="$NAME"_bwa_sorted_Added.bam RGLB=Libreria RGPL=illumina RGPU=unit1 RGSM=20\

samtools index "$NAME"_bwa_sorted_Added.bam

time java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144" -jar \
$EBROOTPICARD/picard.jar \
ValidateSamFile \
I="$NAME"_bwa_sorted_Added.bam \
MODE=SUMMARY 1>job.out 2>job.err

rm ._*

NAME=$(find . -maxdepth 1 -name "*.fastq.gz" | xargs -I {} basename {} .fastq.gz)

NAME2=$(find . -maxdepth 1 -name "*_curated_varpatch.fa" | xargs -I {} basename {} _curated_varpatch.fa)

echo "name of the file is : $NAME"

samtools faidx "$NAME2"_curated_varpatch.fa

time java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144" -jar \
$EBROOTPICARD/picard.jar \
CreateSequenceDictionary \
REFERENCE="$NAME2"_curated_varpatch.fa \
OUTPUT="$NAME2"_curated_varpatch.dict

module load GATK/3.8-1-Java-1.8.0_144

time  java -jar $EBROOTGATK/GenomeAnalysisTK.jar \
-T RealignerTargetCreator \
-R "$NAME2"_curated_varpatch.fa \
-I "$NAME"_bwa_sorted_Added.bam \
-o target_intervals.list

module load LoFreq/2.1.3.1-foss-2019b-Python-2.7.16
#Source http://csb5.github.io/lofreq/commands/#call

time lofreq call -f "$NAME2"_curated_varpatch.fa \
-o "$NAME"_bwa_sorted_realigned.vcf "$NAME"_bwa_sorted_Added.bam

time java -jar $EBROOTGATK/GenomeAnalysisTK.jar \
-T IndelRealigner \
--maxReadsForRealignment 100000 \
-R "$NAME2"_curated_varpatch.fa \
-I "$NAME"_bwa_sorted_Added.bam \
-targetIntervals target_intervals.list \
-o "$NAME"_bwa_sorted_realigned.bam

time java -jar $EBROOTGATK/GenomeAnalysisTK.jar \
-T BaseRecalibrator \
-R "$NAME2"_curated_varpatch.fa \
-I "$NAME"_bwa_sorted_realigned.bam \
-knownSites "$NAME"_bwa_sorted_realigned.vcf \
-o recal_data.table

module unload GATK/4.1.6.0-GCCcore-8.3.0-Java-1.8

#module load LoFreq/2.1.2-foss-2016b-Python-2.7.14
#Source http://csb5.github.io/lofreq/commands/#call

time lofreq call -f "$NAME2"_curated_varpatch.fa \
-o "$NAME"_LoFreq.vcf "$NAME"_bwa_sorted_realigned.bam

##Parse columns

grep "^##" -v "$NAME"_LoFreq.vcf | awk 'BEGIN{OFS"\t"} {split($8, a, ";"); print $1, $2, $4, $5, "\t", a[2]}' | \
awk 'BEGIN{OFS"\t"} {split($5, a, "="); print $1, "\t", $2, "\t", $3, "\t", $4, "\t", a[2]}'  > "$NAME"_LoFreq-AF.vcf 
#remove header "CHROM" and keep number corresponding to segment
grep "^#" -v "$NAME"_LoFreq-AF.vcf | awk 'BEGIN{OFS"\t"} {split($1, a, "|"); print a[1], "\t", $2, "\t", $3, "\t", $4, "\t", $5}' > "$NAME"_LoFreq-nuclAF.vcf

samtools depth "$NAME"_bwa_sorted_realigned.bam | awk 'BEGIN{OFS"\t"} {split($1, a, "|"); print a[1], "\t",$2, "\t", $3}' > "$NAME"_coverage.tsv

##Organize results

mkdir 01_Reference
mv WF10_* ./01_Reference
mkdir 02_Alignments
mv *.bam ./02_Alignments
mv *.bai ./02_Alignments
mv *.sam ./02_Alignments
mv recal_data.table ./02_Alignments
mv target_intervals.list ./02_Alignments
mkdir 04_Results
mv *_cutadapt_LoFreq-nuclAF.vcf ./04_Results
mv *_cutadapt_coverage.tsv ./04_Results
mkdir 03_VCF
mv *.vcf ./03_VCF
