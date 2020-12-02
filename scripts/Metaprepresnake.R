#Metabarcoding

#only run this the first time you do the analysis on a new machine
source('https://bioconductor.org/biocLite.R')
biocLite('DECIPHER')
install.packages('phangorn')


#start from here usually 
library(seqinr)
library(phangorn)
library(DECIPHER)
#packageVersion('dada2')

path <- "/home/fwhite/MetabarcodingCOPY"
#list.files(path)
outfile1 <- snakemake@output[["ADAPTERFILE"]]
outfile2 <- snakemake@output[["Bsamples"]]
outfile3 <- snakemake@output[["Fsamples"]]

#read in metadata and MiSeq info
MiSeqTab <- read.table("~/MetabarcodingCOPY/Mehran_Ale_2020_MiSeqReadSet_2020-08-03.csv", sep = ";", header = T, stringsAsFactors = F)
sample_metadata <- read.table("~/MetabarcodingCOPY/sample_metadata.txt", sep = "\t", header = T, stringsAsFactors = F)


#select files based on owner
SAMPS <- sample_metadata[which(sample_metadata$Owner == "Ale"),]
IDS <- c(SAMPS$Sample_id_bac,SAMPS$Sample_id_fun)


MiSeqTab <- MiSeqTab[which(MiSeqTab$Name %in% IDS),]

#split into 16S and ITS bits - different downstream analysis?


MiITS <- MiSeqTab[grep("ITS",MiSeqTab$Library.Name),]

Mi16S <- MiSeqTab[grep("16S",MiSeqTab$Library.Name),]

MakeSampleDict <- function(x){
  
  x$Alias <- as.numeric(x$Alias)
  x <- x[order(x$Alias),]
  y <- as.data.frame(cbind(x$Name, paste0(x$Filename.Prefix,"_R1.fastq.gz"), 
                           paste0(x$Filename.Prefix,"_R2.fastq.gz")))
  colnames(y) <- c("sample","fq1","fq2")
  
  return(y)
  
}

Bac <- MakeSampleDict(Mi16S)
Fung <- MakeSampleDict(MiITS)



#function to get all unique adapter and primer sequences from MiSeq output
getseq <- function(x){
  
  seqs <- unique(c(unlist(unique(strsplit(x$Reverse.Primer.Sequence, split = ";"))),
                   unlist(unique(strsplit(x$Forward.Primer.Sequence, split = ";"))),
                   unlist(unique(strsplit(x$Adaptor.Read.1..NOTE..Usage.is.bound.by.Illumina.Disclaimer.found.on.Nanuq.Project.Page., split = ";"))),
                   unlist(unique(strsplit(x$Adaptor.Read.2..NOTE..Usage.is.bound.by.Illumina.Disclaimer.found.on.Nanuq.Project.Page., split = ";")))))
}


BS <- getseq(Mi16S)
FS <- getseq(MiITS)

SEQS <- unique(c(BS,FS))


write.fasta(as.list(SEQS), names = as.character(c(1:length(SEQS))), file.out = outfile1)
write.table(Bac,outfile2, sep = "\t", row.names = F, quote = FALSE)
write.table(Fung,outfile3, sep = "\t", row.names = F, quote = FALSE)

