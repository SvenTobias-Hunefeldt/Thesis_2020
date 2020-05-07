#Enable packages
library(dada2)
library(phyloseq)  
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(plyr)
library(dplyr)
library(scales)
library(reshape)
library(reshape2)
library(RColorBrewer)
library(Rmisc)
library(grid)
library(microbiome)
library(graphics)
library(tidyr)
library('ggedit')
library(ggrepel)



#Update path inside R
path <- "/SCRATCH/Sven_Data/DADA2/Biofilm/16S/Samples/"

#List the files we have
list.files(path) #Make sure you have the files you want

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE)) #State what files are forward orientated files
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE)) #State what files are reverse orientated files
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)


#Inspect read quality
pdf("Plotquality.pdf") #Save following output as a pdf file
plotQualityProfile(fnFs[1]) #Save only the quality scores of the first sample - forwards orientation
plotQualityProfile(fnRs[1]) #Reverse is worse quality - common in illumina sequencing , and save only the quality scores of the first sample - forwards orientation
dev.off() #Turn off saving as pdf function

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz")) #Save sample names and paths and F orientation as a separate object
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz")) #Save sample names and paths and R orientation as a separate object
names(filtFs) <- sample.names #assign actual sample names to file path - F orientation
names(filtRs) <- sample.names #assign actual sample names to file path - R orientation

#Filter and Trim files according parameters set out below
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(150,150),
              maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
              compress=TRUE, multithread=TRUE, matchIDs=TRUE) # On Windows set multithread=FALSE
head(out) #Make sure its worked

#Issues during previous step can occur - usually due to a replicate sample or it randomly assigning a different name, I have no idea why it does it, but you can check for it using:
#any(duplicated(c(fnFs, fnRs)))
#any(duplicated(c(filtFs, filtRs)))
#Both should be false, if they are and previous step still fails it means that a different name has been assigned to files and you'll have to manually check file names and look for inconsistencies.
#This could also be because barcodes during demultiplexing were mismatched (try deleting one of the reverse options and re-running the whole script) and thus truncation cut off got rid of all reads.


#Learn the error rate
errF <- learnErrors(filtFs, multithread=TRUE) #Train error rate for F orientation
errR <- learnErrors(filtRs, multithread=TRUE) #Train error rate for R orientation

#Sanity check
pdf("Error_rates.pdf") #Save following output as a pdf file
plotErrors(errF, nominalQ=TRUE) #Show error rates of base to base substitutions for F orientation
plotErrors(errF, nominalQ=TRUE) #Show error rates of base to base substitutions for R orientation
dev.off() #Turn off saving as pdf function


##Dereplicate the fastq file and carry out sample Inference
exists <- file.exists(filtFs) #Save sample names of samples that passed filtering and trimming
deRepFs <- derepFastq(filtFs[exists]) #Dereplicate F samples for further processing
deRepRs <- derepFastq(filtRs[exists]) #Dereplicate R samples for further processing


#Carry out sample inference
dadaFs <- dada(deRepFs, err=errF, multithread=TRUE) #Maybe look at pooling in the future, see how it affects ASV creation?
dadaRs <- dada(deRepRs, err=errR, multithread=TRUE) #Maybe look at pooling in the future, see how it affects ASV creation?
dadaFs[[1]] #Make sure it worked


##Merge paired reads
mergers <- mergePairs(dadaFs, deRepFs, dadaRs, deRepRs, verbose=TRUE) #Merge F and R samples keeping in mind error rates
# Inspect the merger data.frame from the first sample
head(mergers[[1]]) #Make sure its worked


##Construct a sequence Table
seqtab <- makeSequenceTable(mergers) #Make the by row merged samples into a table.
dim(seqtab) #Check dimensions to see if the expected number of samples have made it into the table with the expected number of reads

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))


##Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim) #Check dimensions to see if the a good number of samples have made it into the table with a good number of reads

sum(seqtab.nochim)/sum(seqtab) #Check the proportion of reads that made it past the chimera check

saveRDS(seqtab.nochim, file=paste(path, "../seqtab.nochim.rds", sep = "")) #Save table as .rds object so we can use it later in a local version of R to create a tree


##Track reads through the pipeline
getN <- function(x) sum(getUniques(x)) #See how many unique sequences there are? Very tired, will check later.
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim)) #As it says on the box
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim") #Set column names
rownames(track) <- sample.names #Label rownames for each sample
head(track) #Check to see the proportion of samples that made it through the whole process.

write.csv(track, "Track_reads_through_pipeline.csv") #Save the tracked number of reads through the process


##Assign Taxonomy
#Only goes down to family
taxa <- assignTaxonomy(seqtab.nochim, paste(path, "taxa_DB/silva_nr_v132_train_set.fa.gz", sep = ""), multithread=TRUE)

taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL #Superflous information, get rid
head(taxa.print) #See if it looks good

#Based on exact matching and goes down to species if possible
taxa <- addSpecies(taxa, paste(path, "taxa_DB/silva_species_assignment_v132.fa.gz", sep = "")) #As on tin

taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL #Superflous information, get rid
head(taxa.print) #See if it looks good

write.csv(taxa.print, "Taxa.csv") #Save taxonomic information for later perusal



##Convert to phyloseq object
map_file <- paste(path, "Mapping_file.txt", sep = "") #Define pathway for mapping file
bmsd <- import_qiime_sample_data(map_file) #Import mapping file to later use for phyloseq object creation

myotutable = otu_table(seqtab.nochim, taxa_are_rows=FALSE) #Create otutable object
rownames(myotutable) = gsub(".fastq", "", rownames(myotutable)) #Remove last part of rownames so it matches across the three files.
mysampledata = sample_data(bmsd) #Create sample data object
mytaxatable = tax_table(taxa)#Create taxa table object
myphyloseq = phyloseq(myotutable, mysampledata, mytaxatable) #Integrate all objects into a phyloseq object


#Store DNA sequences of ASVs in the phyloseq refseq slot, renaming taxa to a short string
dna <- Biostrings::DNAStringSet(taxa_names(myphyloseq)) #Save 
names(dna) <- taxa_names(myphyloseq)
myphyloseq <- merge_phyloseq(myphyloseq, dna)
taxa_names(myphyloseq) <- paste0("ASV", seq(ntaxa(myphyloseq)))
myphyloseq #Check to see if it worked

saveRDS(myphyloseq, file=paste(path, "../Phyloseq_object.rds", sep = "")) #Save phyloseq as .rds object so we can use it later in a local version of R for data analysis
