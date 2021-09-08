#!/usr/local/bin/Rscript
## This script was made by LMF - 2019
# Modified in 20200727
## This script will output variants with coverage of 400X and frequency of 0.02.
## Takes input from LoFreq Variant analysis pipeline
library(readr)

# Import table with coverage values
file1 <-list.files(path = ".", pattern = "\\.tsv")
coverage <- read_delim(file1, "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
colnames(coverage) <- c("Segment", "Position", "Cov")
# Index Positions in Segments
coverage["length"] <- seq.int(nrow(coverage))

# Import .vcf file from GATK-Lofreq pipeline (GATK-LoFreq_VariantAnalysis_II.sh)
file2 <- list.files(path = ".", pattern = "\\.vcf")
Varfreq <- read.table(file2, header = FALSE, sep = "\t")
colnames(Varfreq) <- c("Segment", "Position", "REF", "ALT", "AF") 
# Filter variants based on frequency >= 0.02
Varfreq <- Varfreq[Varfreq$AF >= 0.02, ]
# Position 4 shows as variant in most segments because of the use of degenerated primers
# Erase these variants
Varfreq$AF[Varfreq$Position == 4] <- NA #assigns NA to position 4 since this position is detected as variant because of the primers
Varfreq$REF[Varfreq$Position == 4] <- NA
Varfreq$ALT[Varfreq$Position == 4] <- NA

# Merge *Varfreq with coverage table 
AF.Cov <- merge(coverage, Varfreq, by = c("Segment", "Position"), all = TRUE)

# Filter out variants with coverage higher than 400
AF.Cov$AF[AF.Cov$Cov < 400] <- NA 

#Delete "NA" rows
AF.Cov.vcf<-na.omit(AF.Cov)
#export the table
write_tsv(AF.Cov.vcf, "AF.Cov.tsv", quote=FALSE)