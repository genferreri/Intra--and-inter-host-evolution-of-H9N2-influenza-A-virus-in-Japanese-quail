#!/bin/bash

NAME=$(find . -maxdepth 1 -name "*_cutadapt-Segments_cont_LoFreq.vcf" | xargs -I {} basename {} _cutadapt-Segments_cont_LoFreq.vcf)

perl snpgenie.pl --minfreq=0.02 --vcfformat=2 --snpreport="$NAME"_cutadapt-Segments_cont_LoFreq.vcf \
--fastafile=WF10-Segments_cont.fa --gtffile=WF10-Segments_cont.gtf.txt 

cd SNPGenie_Results

cat population_summary.txt | awk 'BEGIN{OFS"\t"} {split($1, a, "-"); print a[2] "\t" $5 "\t" $10 "\t" $11}' \
| awk 'BEGIN{OFS"\t"} {split($1, a, "_"); print a[1] "\t" $2 "\t" $3 "\t" $4}' \
| awk '{if (NR!=1) {print}}' > population_summary_clean.txt

