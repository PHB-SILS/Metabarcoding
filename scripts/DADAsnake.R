#!/usr/bin/env Rscript

# Rscript --vanilla deseq2_normalization.R raw_counts.txt scaled_counts.txt

args <- commandArgs(trailingOnly=TRUE)


#Metabarcoding

#only run this the first time you do the analysis on a new machine
#source('https://bioconductor.org/biocLite.R')
#biocLite('dada2')
#biocLite('phyloseq')
#biocLite('DECIPHER')
#install.packages('ggplot2')
#install.packages('phangorn')


#start from here usually 
library(seqinr)
library(dada2)
library(ggplot2)
library(phyloseq)
library(phangorn)
library(DECIPHER)
library(doMC)
#packageVersion('dada2')

threads <- 8

path <- paste0(getwd(),"/temp/trimmed")
#list.files(path)


fnFs <- sort(list.files(path, pattern="_R1_trimmed.fq.gz",
                        full.names=TRUE))

fnRs <- sort(list.files(path, pattern="_R2_trimmed.fq.gz",
                        full.names=TRUE))


# we also need the sample names
sample.names <- gsub(".*/","",fnFs)
sample.names <- gsub("_tr.*","",sample.names)

#need to make sample.names here the same as the {sample} wildcard
#sample.names <- snakemake@wildcards[['sample']]


pdf("test20201201.pdf")
plotQualityProfile(fnFs[1:2])
dev.off()

#####
#split the script here for inspection of quality plots and filter parameter decision
#####

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))





##this can open up limitless processes, make sure multithread option is specified!
##snakemake  can also block this
##this is set to default parameters following DADA2 1.8 tutorial. should be changed based on QC plots or better still, automated using FIGARO or equivalent

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, minLen = 20,
                     maxN=0, maxEE=c(2,2), truncQ=2, matchIDs=TRUE, rm.phix=TRUE,
                     compress=TRUE, multithread=threads)


#learn error rates
#nbases argument to increase the number of bases that the error is learned from

errF <- learnErrors(filtFs, multithread = threads)
errR <- learnErrors(filtRs, multithread = threads)

pdf(paste0(Sys.Date(),"_errorsF.pdf"))
plotErrors(errF, nominalQ=TRUE)
dev.off()

pdf(paste0(Sys.Date(),"_errorsR.pdf"))
plotErrors(errR, nominalQ=TRUE)
dev.off()

#dereplication - combine all identical sequencing reads into unique sequences with corresponding abundances


derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names


#sample inference

dadaFs <- dada(derepFs, err=errF, multithread=threads)
dadaRs <- dada(derepRs, err=errR, multithread=threads)


mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
#head(mergers[[1]])

#make seq table
seqtab <- makeSequenceTable(mergers)
#check length of seqs
#table(nchar(getSequences(seqtab)))

#remove chimeras

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=threads, verbose=TRUE)

#check chimera freq
#sum(seqtab.nochim)/sum(seqtab)



#check where reads are lost along the dada2 process
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))

# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
#head(track)

taxa <- assignTaxonomy(seqtab.nochim, args[1], multithread=threads)

taxa <- addSpecies(taxa, args[2])
#taxa <- addSpecies(taxa, args[2])

#outfile1 <- snakemake@output[["taxa"]]
#outfile2 <- snakemake@output[["asv"]]

#write taxa and ASV tables as snakemake outputs
write.table(taxa,args[3])
write.table(seqtab.nochim,args[4])

save.image("debugdada")

