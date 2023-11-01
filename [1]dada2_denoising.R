# imposed on https://github.com/charlie-carr/implant_microbiota 
# Carr C, Wilcox H, Burton JP, Menon S, Al KF, O'Gorman D, et al. (2021) Deciphering the low abundance microbiota of presumed aseptic hip and knee implants. PLoS ONE 16(9): e0257471. https://doi.org/10.1371/journal.pone.0257471
# and
# https://benjjneb.github.io/dada2/tutorial.html
# https://ryjohnson09.netlify.app/post/microbiome-analysis-with-dada2-and-phyloseq/

set.seed(1602000)
library(dada2)
library(dplyr)

setwd("C:/Users/ASUS/AppData/Local/Packages/CanonicalGroupLimited.Ubuntu20.04LTS_79rhkp1fndgsc/LocalState/rootfs/home/ocetiner/2209_raw_files/fastq_files")

# specify path containing demultiplexed fastq files
path <- "C:/Users/ASUS/AppData/Local/Packages/CanonicalGroupLimited.Ubuntu20.04LTS_79rhkp1fndgsc/LocalState/rootfs/home/ocetiner/2209_raw_files/fastq_files"
list.files(path)

# We have decided to use only forward reads, 
# because the quality of reverse reads were low, which was detected using MultiQC tool.

# keep fastq filenames, filename format : SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# Duplicated files
any(duplicated(sample.names))

# quality control to decided truncation and trim lengths 
plotQualityProfile(fnFs[1:2])
# plotQualityProfile(fnRs[1:2])

# take four samples and review on pdf

ids <- round(runif(4,1,length(sample.names)))

pdf("figures/final_quality_profiles.pdf")
plotQualityProfile(fnFs[ids])
dev.off()

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
names(filtFs) <- sample.names


out <- filterAndTrim(fnFs, filtFs,truncLen=c(260),
                     trimLeft = 32) 

# Learn the error rates
errF <- learnErrors(filtFs, multithread=TRUE)

# graph the error rates
pdf("figures/error_plot_final_F_3.pdf")
plotErrors(errF, nominalQ = TRUE)
dev.off()

# Use the filtered files and error rates to perform
# sample inference (without pooling)
dadaFs <- dada(filtFs, err = errF, multithread = TRUE)

# extract the read counts 
counts_raw <- makeSequenceTable(dadaFs)
dim(counts_raw)
sum(counts_raw)

# Filter chimeras
counts_nochim <- removeBimeraDenovo(counts_raw, multithread = TRUE, verbose = TRUE)

sum(counts_nochim)/sum(counts_raw)

# assign taxonomy provided by SILVA

tax_nochim <- assignTaxonomy(counts_nochim, "C:/Users/ASUS/Desktop/ferm_project/raw_data_and_taxonomy_db/silva_nr99_v138_wSpecies_train_set.fa.gz", multithread=TRUE)
tax_nochim <- as.data.frame(tax_nochim)

counts_filtered <- counts_nochim[, rownames(tax_nochim)]

# Assign readable names
tax_nochim_edited <- tax_nochim %>%
  mutate(Sequence = rownames(tax_nochim))
rownames(tax_nochim_edited) <- paste0("ASV_", 1:nrow(tax_nochim_edited))

any(colnames(counts_filtered) != tax_nochim_edited$Sequence) # FALSE, so the ASVs are in the same order
colnames(counts_filtered) <- paste0("ASV_", 1:ncol(counts_filtered))

# Construct a table to summarize the removal of reads throughout the pipeline
getN <- function(x) {sum(getUniques(x))}
track <- cbind(out, sapply(dadaFs, getN), rowSums(counts_nochim))
colnames(track) <- c("Input", "Filtered", "DenoisedF", "Non-Chimeric")
rownames(track) <- sample.names

# Output the counts, tax, and tracking tables
write.table(counts_filtered, 
            file = "C:/Users/ASUS/Desktop/prosthesis_microbiota_study/data/counts_silva.txt", 
            sep = "\t", col.names = NA, quote = F)

write.table(tax_nochim_edited, 
            file = "C:/Users/ASUS/Desktop/prosthesis_microbiota_study/data/tax_silva.txt",
            sep = "\t", col.names = NA, quote = F)

write.table(track, 
            file = "C:/Users/ASUS/Desktop/prosthesis_microbiota_study/data/track_silva.txt",
            sep = "\t", col.names = NA, quote = F)

# Create a phyloseq object from the dada2 output
library(phyloseq)
# Forward
metadata <- read.table("C:/Users/ASUS/Desktop/prosthesis_microbiota_study/data/metadata.txt",
                       header = TRUE, sep = "")

metadata$sample <- rownames(metadata)

tax_nochim_edited <- tax_nochim_edited[,!names(tax_nochim_edited) %in% c("Sequence")]

ps <- phyloseq(otu_table(counts_filtered, taxa_are_rows=FALSE), 
               sample_data(metadata), 
               tax_table(as.matrix(tax_nochim_edited)))

#Save the phyloseq object
saveRDS(ps,"C:/Users/ASUS/Desktop/prosthesis_microbiota_study/R_Scripts/PJI_forward_ps_object2.rds")
