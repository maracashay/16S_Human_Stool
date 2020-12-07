# Intial processing and filtering was performed following the dada2 pipeline #
setwd("~/Downloads/New Samples - set 4") #set working directory first

path<-setwd("~/Downloads/New Samples - set 4") 

list.files(path)


#### 1. Start here if primers are already trimmed####

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
#do sample.names first because we need to split the strings before sorting
sample.names
#### 2. Quality Trimming and Filtering####

plotQualityProfile(fnFs[5:9])
plotQualityProfile(fnRs[5:9])

filt_path.new <- file.path(path, "filtered")
#here we are putting the files into a folder called "filtered"

filtFs.new <- file.path(filt_path.new, paste0(sample.names, "_F_filt.fastq"))
filtRs.new <- file.path(filt_path.new, paste0(sample.names, "_R_filt.fastq"))


#for truncLen : better to be more  conservative for this data, since we
# had poorer quality samples; so later, when we track the reads through the
# pipeline, we will see less samples that ultimately made it through

out.4 <- filterAndTrim(fnFs, filtFs.new, fnRs, filtRs.new, truncLen=c(270,210), trimLeft = c(20,5),
                       maxN=1, maxEE=c(3,3), truncQ=2,
                       compress=TRUE, multithread=FALSE, matchIDs = TRUE) # On Windows set multithread=FALSE



#### 3. Learn Error Rates####
errF.new <- learnErrors(filtFs.new, multithread=TRUE)
errR.new <- learnErrors(filtRs.new, multithread=TRUE)

plotErrors(errF.new, nominalQ=TRUE)

#### 4. Dereplication####
#see which sequences are the exact same; if so, then it will count the unique sequences
#count uniqe ESVs in forward reads and reverse reads
#derepFs and derepRs are file names
derepFs.cnt <- derepFastq(filtFs.new, verbose=TRUE)
derepRs.cnt <- derepFastq(filtRs.new, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs.cnt) <- sample.names
names(derepRs.cnt) <- sample.names

dadaFs.cnt <- dada(derepFs.cnt, err=errF.new, multithread=TRUE)
dadaRs.cnt <- dada(derepRs.cnt, err=errR.new, multithread=TRUE)

dadaFs.cnt[[1]]
dadaRs.cnt[[5]]

mergers.cnt <- mergePairs(dadaFs.cnt, derepFs.cnt, dadaRs.cnt, derepRs.cnt, verbose=TRUE)
# Inspect the merger data.frame from the first sample

head(mergers.cnt[[1]])


seqtab.4 <- makeSequenceTable(mergers.cnt)
#it will say sequences being tables vary in length, but that's normal
dim(seqtab.cnt)
#this will say how many samples and total ESVs


table(nchar(getSequences(seqtab.1)))
# majority of sequences from set 1 were between 414 and 441
table(nchar(getSequences(seqtab.2)))
# majority of sequences from set 2 were betwewn 413 and 441
table(nchar(getSequences(seq.tab.3)))
# majority of sequences froms et 3 were between 414 and 441
table(nchar(getSequences(seqtab.4)))
# majority of sequences from set 4 were between 413 and 441

seqtab1.new <- seqtab.1[,nchar(colnames(seqtab.1)) %in% seq(413,441)]
seqtab2.new <- seqtab.2[,nchar(colnames(seqtab.2)) %in% seq(413,441)]
seqtab3.new <- seq.tab.3[,nchar(colnames(seq.tab.3)) %in% seq(413,441)]
seqtab4.new <- seqtab.4[,nchar(colnames(seqtab.4)) %in% seq(413,441)]



#### 5. Remove Chimeras####
#Error in PCR, a chimera is essentially 2 sequences from 2 biological samples

seqtab.nochim.1 <- removeBimeraDenovo(seqtab1.new, method="consensus", multithread=TRUE, verbose=TRUE)
seqtab.nochim.2 <- removeBimeraDenovo(seqtab2.new, method="consensus", multithread=TRUE, verbose=TRUE)
seqtab.nochim.3 <- removeBimeraDenovo(seqtab3.new, method="consensus", multithread=TRUE, verbose=TRUE)
seqtab.nochim.4 <- removeBimeraDenovo(seqtab4.new, method="consensus", multithread=TRUE, verbose=TRUE)


dim(seqtab.nochim.1)
# out of 57 samples, 5085 seqs were kept

dim(seqtab.nochim.2)
# out of 79 samples, 8573 seqs were kept
dim(seqtab.nochim.3)
# out of 70 samples, 4099 seqs were kept
dim(seqtab.nochim.4)
# out of 63 samples, 6595 seqs were kept

sum(seqtab.nochim.1)/sum(seqtab.1)
# 66.91
sum(seqtab.nochim.2)/sum(seqtab.2)
# 69.42
sum(seqtab.nochim.3)/sum(seq.tab.3)
# 78.88
sum(seqtab.nochim.4)/sum(seqtab.4)
# 59.25

#### 6. Track reads through the pipeline####

getN.new <- function(x) sum(getUniques(x))
track.4 <- cbind(out.4, sapply(dadaFs.cnt, getN.new), sapply(dadaRs.cnt, getN.new), sapply(mergers.cnt, getN.new), rowSums(seqtab.nochim.4))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track.4) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")

rownames(track.new) <- sample.names
head(track.new)

#  merge sequence tables 
seqtab.full <- mergeSequenceTables(seqtab.nochim.1, seqtab.nochim.2, seqtab.nochim.3, seqtab.nochim.42)

taxa.full<- assignTaxonomy(seqtab.full, "/Users/maraslaptop/Downloads/silva_nr_v132_train_set.fa", multithread=FALSE)

taxa.spp <- addSpecies(taxa.full, "/Users/maraslaptop/Downloads/silva_species_assignment_v132.fa")

# 7. Make Phyloseq Object ####
ps.full <- merge_phyloseq(tax_table(taxa.spp), otu_table(seqtab.full, taxa_are_rows = FALSE), sample_data(map), phy_tree(fitGTR$tree))
dna <- Biostrings::DNAStringSet(taxa_names(ps.full))
names(dna) <- taxa_names(ps.full)
ps.full <- merge_phyloseq(ps.full, dna)
taxa_names(ps.full) <- paste0("ASV", seq(ntaxa(ps.full)))
ps.full

# 8. Filter phyloseq ####
ps.stool <- subset_taxa(ps.full, Kingdom == "Bacteria" | Kingdom == "Archaea")

ps.stool.20 = filter_taxa(ps.stool, function(x) sum(x > 0) > (0.2*length(x)), TRUE)



