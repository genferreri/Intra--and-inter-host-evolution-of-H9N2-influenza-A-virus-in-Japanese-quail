#!/bin/bash
#20200626
#This is the latest version of ParseALT_SNPdat.sh
#This script needs the reference sequence WF10_curated_varpatch.fa, SNPdat_v1.0.5.pl, 
# as well as the GTF file WF10ref.gtf.txt

#Run SNPdat
for f in *.txt; do FILENAME=${f%%.*};
perl SNPdat_v1.0.5.pl -i ${FILENAME}.AF.Cov.ALT.vcf.txt -f WF10_curated_varpatch.fa -g WF10ref.gtf.txt;
done; 

#Parse columns
for r in *.output; do FILENAME=${r%%.*};
cat ${FILENAME}.AF.Cov.ALT.vcf.txt.output | awk 'NR>1' | \
awk 'BEGIN{OFS"\t"} {split($21, a, "/"); print $1 "\t" $2 "\t" a[1] "\t" a[2] "\t" $22}' | \
awk 'BEGIN{OFS"\t"} {split($3, a, "["); split ($4, b, "]"); print $1 "\t" $2 "\t" a[2] "\t" b[1] "\t" $5}' \
> ${FILENAME}.SNPdat.temp.txt
tr -d '-' < ${FILENAME}.SNPdat.temp.txt > ${FILENAME}.SNPdat.txt


done;

rm *.summary 
rm *.SNPdat.temp.txt
mkdir SNPdat
rm WF10*.output
rm WF10ref.SNPdat.txt
mv *vcf.txt ./SNPdat
mv *.output ./SNPdat

## Add sample's name
for f in *.SNPdat.txt
do 
	NAME=${f%%.*}
	echo $NAME
	##add name of the sample
	awk -F, -v s="$NAME" '{$6="\t"s; print }' "$NAME".SNPdat.txt > "$NAME".SNPdat.temp.txt

	rm "$NAME".SNPdat.txt

	##remove "-" from name
	tr -d '-' < "$NAME".SNPdat.temp.txt > "$NAME".SNPdat.txt

	rm "$NAME".SNPdat.temp.txt

	#awk 'BEGIN { FS = OFS = "\t" } { for(i=1; i<=NF; i++) if($i ~ /^ *$/) $i = "NA" }; 1' "$NAME".SNPdat.t.txt \
	#> "$NAME".SNPdat.txt

	#rm "$NAME".SNPdat.t.txt

done
