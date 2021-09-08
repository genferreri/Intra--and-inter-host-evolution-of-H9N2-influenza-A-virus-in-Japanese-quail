#!/usr/bin/env R
 
file1 <- list.files(path = ".", pattern = "\\.txt")

library(readr)
#import the dataset in tabbed delimited format and create the necessary columns to format and populate during the R analysis
hadata = read_delim(file1, "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
colnames(hadata) <- c("seq.nt")
hadata[c("nt.DNAString", "Seq.AA", "p226", "p227")]<-NA #create new columns

library(Biostrings)
hadata$nt.DNAString <- DNAStringSet(hadata$seq.nt)
hadata$Seq.AA <- Biostrings::translate(hadata$nt.DNAString, if.fuzzy.codon = "X", genetic.code = GENETIC_CODE) #tranlslate nt's to AAs. Call translate function specifically from Biostrings package.
sapply(hadata, typeof) #check that the data types are as expected, for translation of nt's you need vector class data or S4 class type
hadata$p226 <- paste(subseq(hadata$Seq.AA, 5,5))
hadata$p225 <- paste(subseq(hadata$Seq.AA, 4,4))

library(plyr)
p226.freq <- plyr::count(hadata, 'p226')
p226.freq[c("Percent")]<-NA #create new columns
p226.freq$Percent <- (p226.freq$freq/sum(p226.freq$freq))*100 #convert to % frequency

#p226.freq[[1]]<-NULL
colnames(p226.freq) <- c("AminoAcid", "ReadCount", "Freq")
p226.freq<-as.data.frame(p226.freq)

## create df with all possible AAs
AminoAcids <- c('G', 'P', 'A', 'V', 'L', 'I', 'M', 'C', 'F', 'Y', 'W', 'H', 'K', 'R', 'Q', 'N', 'E', 'D', 'S', 'T')
AminoAcids <- dplyr::as_data_frame(AminoAcids)
AminoAcids <- plyr::rename(AminoAcids, c("value"="AminoAcid"))

p226.freq<-merge(p226.freq, AminoAcids, by = c("AminoAcid"), all = TRUE)
p226.freq$ReadCount[is.na(p226.freq$ReadCount)]<-0
p226.freq$Freq[is.na(p226.freq$Freq)]<-0

#write.csv(p226.freq, "p226freq.csv")
 write.table(p226.freq, file = 'p226freq.tsv', quote=FALSE, sep='\t', col.names=NA)