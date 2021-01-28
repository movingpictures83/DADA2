#First must have dada2 and Bioconductor installed: https://benjjneb.github.io/dada2/dada-installation.html 
#Must also have phyloseq installed: http://joey711.github.io/phyloseq/install.html
dyn.load(paste("RPluMA", .Platform$dynlib.ext, sep=""))
source("RPluMA.R")

library(dada2); 
packageVersion("dada2")


input <- function(inputfile) {
  pfix = prefix()
  if (length(pfix) != 0) {
     pfix <- paste(pfix, "/", sep="")
  }

  parameters <- read.table(inputfile, as.is=T);
  rownames(parameters) <- parameters[,1]
  fastqdir <<- toString(parameters["FASTQ", 2])
  DB <<- paste(pfix, toString(parameters["DB", 2]), sep="")
  specDB <<- paste(pfix, toString(parameters["species", 2]), sep="")

path <<- paste(pfix, fastqdir, sep="") # CHANGE ME to the directory containing the fastq files after unzipping.
list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <<- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <<- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <<- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

#Inspect read quality profiles------------------

# Calloc error on this
#plotQualityProfile(fnFs, n = 5e+09, aggregate = TRUE)
#plotQualityProfile(fnRs, n = 5e+09, aggregate = TRUE)

#Filter and trim--------------------------------

#Place filtered files in filtered/ subdirectory
}

run <- function() {
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(250, 225), trimLeft = c(17, 21),
                     maxN=0, maxEE=6, truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE, n = 1e+09) # On Windows set multithread=FALSE
head(out)
mean(out[,1]) #average number of input reads per sanmple
mean(out[,2]/out[,1]) #Percentage of reads passing through filterAndTrim

#Inspect quality profiles of filtered reads------------------

plotQualityProfile(filtFs, aggregate = TRUE)

plotQualityProfile(filtRs, aggregate = TRUE)

#Dereplication----------------------------------------

derep_forward <- derepFastq(filtFs, verbose=TRUE)
derep_reverse <- derepFastq(filtRs, verbose=TRUE)

#Learn the error rates--------------------------------

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

plotErrors(errF, nominalQ=TRUE)

#Sample Inference-------------------------------------

dadaFs <- dada(derep_forward, err=errF, multithread=TRUE)
dadaRs <- dada(derep_reverse, err=errR, multithread=TRUE)

dadaFs[[1]]

merger1 <- mergePairs(dadaFs, derep_forward, dadaRs, derep_reverse, minOverlap = 5, verbose=TRUE)

#Remove chimeras--------------------------------------

merger1.nochim <- removeBimeraDenovo(merger1, multithread=FALSE, verbose=TRUE)

#Inspect the merger data.frame

head(merger1.nochim[[1]])

#Construct sequence table------------------------------

seqtab.nochim <<- makeSequenceTable(merger1.nochim) #sequence table with chimeras removed
dim(seqtab.nochim) #how many ASVs you have

# Inspect distribution of sequence lengths

table(nchar(getSequences(seqtab.nochim)))

#sum(seqtab.nochim)/sum(seqtab) #percentage of merged sequences that are NOT chimeras

#Track reads through pipeline---------------------------

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(merger1.nochim, getN), rowSums(seqtab.nochim))

# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
mean(track[,5]/track[,4]) #percentage of denoised sequences that have merged

#Assign taxonomy----------------------------------------

taxa <<- assignTaxonomy(seqtab.nochim, DB, multithread=TRUE)
taxa <<- addSpecies(taxa, specDB)

taxa.print <<- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

#Extract info--------------------------------------------

asv_seqs <<- colnames(seqtab.nochim)
asv_headers <<- vector(dim(seqtab.nochim)[2], mode="character")

for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <<- paste(">ASV", i, sep="_")
}
}

output <- function(outputfile) {
# making and writing out a fasta of our final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, paste(outputfile, ".fa", sep=""))
asv_fasta

# count table:
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, paste(outputfile, ".counts.tsv", sep=""), sep="\t", quote=F, col.names=NA)
asv_tab

# tax table:
asv_tax <- taxa
row.names(asv_tax) <- sub(">", "", asv_headers)
write.table(asv_tax, paste(outputfile, ".taxonomy.tsv", sep=""), sep="\t", quote=F, col.names=NA)
asv_tax

#for export of the tax table as a csv

#asv_tax1 = as.data.frame(asv_tax) 
#write.csv(asv_tax1, file = '~/Desktop/asv_tax.csv')
}


#Bonus: Handoff to phyloseq------------------------------

#library(phyloseq); packageVersion("phyloseq")
#library(Biostrings); packageVersion("Biostrings")
#library(ggplot2); packageVersion("ggplot2")

#Creating abundance bar plots------------------------------------

#TAB = otu_table(asv_tab, taxa_are_rows=TRUE)
#TAX = tax_table(asv_tax)

#physeq = phyloseq(TAB, TAX)

#plot_bar(physeq, fill = "Phylum")
#plot_bar(physeq, fill = "Family")
#plot_bar(physeq, fill = "Genus")

#nsamples(physeq)
#sample_names(physeq)
#sample_variables(physeq)

