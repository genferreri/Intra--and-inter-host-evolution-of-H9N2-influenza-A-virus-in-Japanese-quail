#!/bin/bash

NAME1=$(find . -maxdepth 1 -name "*_cutadapt_coverage.tsv" | xargs -I {} basename {} _cutadapt_coverage.tsv)

Rscript ~/scripts/iSNV/1_Auto_VCFvariant.R

mv AF.Cov.tsv "$NAME1"_AF.Cov.tsv

##Add name of the sample
awk -F, -v s="$NAME1" '{$8=s; print }' "$NAME1"_AF.Cov.tsv > "$NAME1"_AF.Cov.temp.tsv

mv "$NAME1"_AF.Cov.tsv "$NAME1"_AF.Cov
##Delete header
awk '{if (NR!=1) {print}}' "$NAME1"_AF.Cov.temp.tsv > "$NAME1"_AF.Cov.clean.temp.tsv

rm "$NAME1"_AF.Cov.temp.tsv

mv "$NAME1"_AF.Cov.clean.temp.tsv "$NAME1"_AF.Cov.clean.temp

## remove "-" from Sample's name
tr -d '-' < "$NAME1"_AF.Cov.clean.temp > "$NAME1"_AF.Cov.inter.tsv 
rm "$NAME1"_AF.Cov.clean.temp

##Split last column
cat "$NAME1"_AF.Cov.inter.tsv | awk 'BEGIN{OFS"\t"} {split($8, a, " "); print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" a[1]}'\
> "$NAME1"_AF.Cov.interm.tsv

rm "$NAME1"_AF.Cov.inter.tsv
mv "$NAME1"_cutadapt_coverage.tsv "$NAME1"_cutadapt_coverage

#Add new headers using R
Rscript ~/scripts/iSNV/colnameSample.R
##Remove "" and lead column from R
tr -d '"' < AF.Cov.clean.temp.tsv | cut -f2- > AF.Cov.clean.temp.out

mv AF.Cov.clean.temp.out "$NAME1"_AF.Cov.final.tsv
rm AF.Cov.clean.temp.tsv
rm *_AF.Cov.in*
rm *_AF.Cov

#Produce file for SNPdat
cat "$NAME1"_AF.Cov.final.tsv | awk 'NR>1' | awk 'BEGIN{OFS"\t"} {print $1 "\t" $2 "\t" $5}' > "$NAME1"_AF.Cov-ALT.vcf.txt
