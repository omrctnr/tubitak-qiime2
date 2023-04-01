# imposed on https://github.com/charlie-carr/implant_microbiota 
# Carr C, Wilcox H, Burton JP, Menon S, Al KF, O'Gorman D, et al. (2021) Deciphering the low abundance microbiota of presumed aseptic hip and knee implants. PLoS ONE 16(9): e0257471. https://doi.org/10.1371/journal.pone.0257471
# and
# https://benjjneb.github.io/dada2/tutorial.html
# https://ryjohnson09.netlify.app/post/microbiome-analysis-with-dada2-and-phyloseq/

library(dada2)
library(dplyr)


setwd("C:/Users/ASUS/AppData/Local/Packages/CanonicalGroupLimited.Ubuntu20.04LTS_79rhkp1fndgsc/LocalState/rootfs/home/ocetiner/2209_raw_files/fastq_files")

# specify path containing demultiplexed fastq files
path <- "C:/Users/ASUS/AppData/Local/Packages/CanonicalGroupLimited.Ubuntu20.04LTS_79rhkp1fndgsc/LocalState/rootfs/home/ocetiner/2209_raw_files/fastq_files"
list.files(path)

# keep fastq filenames, filename format : SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# Duplicated files
any(duplicated(sample.names))

# quality control to decided truncation and trim lengths 
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])

# take four samples and review on pdf

ids <- round(runif(4,1,length(sample.names)))

pdf("figures/final_quality_profiles.pdf")
plotQualityProfile(fnFs[ids])
plotQualityProfile(fnRs[ids])
dev.off()

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# filter and trim process
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(260,220),
                     trimLeft = 11) 
head(out)

# Learn the error rates
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

# graph the error rates
pdf("figures/error_plot_final_F.pdf")
plotErrors(errF, nominalQ = TRUE)
dev.off()

pdf("figures/error_plot_final_R.pdf")
plotErrors(errR, nominalQ = TRUE)
dev.off()

# Use the filtered files and error rates to perform
# sample inference (without pooling)
dadaFs <- dada(filtFs, err = errF, multithread = TRUE,)
dadaRs <- dada(filtRs, err = errR, multithread = TRUE)


# Merge the forward and reverse paired reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose = TRUE)

# Build the counts table
counts_raw <- makeSequenceTable(mergers)

# Check that all of the sequence lengths are within the expected range
table(nchar(getSequences(counts_raw)))

# Remove ASVs that are not enough lengths 
counts_trimmed <- counts_raw[, nchar(colnames(counts_raw)) %in% seq(416, 445)]

# Filter chimeras
counts_nochim <- removeBimeraDenovo(counts_trimmed, method = "consensus", 
                                    multithread = TRUE, verbose = TRUE)

sum(counts_nochim)/sum(counts_trimmed)

# assign taxonomy provided by Greengenes

tax_nochim <- assignTaxonomy(counts_nochim,
                             "C:/Users/ASUS/Desktop/prosthesis_microbiota_study/data/gg_13_8_train_set_97.fa.gz",
                             multithread = TRUE)

# Filter by taxonomy
tax_filtered <- as.data.frame(tax_nochim) %>%
  filter(!is.na(Kingdom)) %>%
  filter(Kingdom != "k__Archaea") %>%
  filter(Family != "f__mitochondria") %>%
  filter(Class != "c__Chloroplast")%>%
  mutate_all(list(~ substr(., 4, nchar(.))))%>%
  mutate_all(na_if,"")

counts_filtered <- counts_nochim[, rownames(tax_filtered)]

# Assign readable names
tax_filtered <- tax_filtered %>%
  mutate(Sequence = rownames(tax_filtered))
rownames(tax_filtered) <- paste0("ASV_", 1:nrow(tax_filtered))

any(colnames(counts_filtered) != tax_filtered$Sequence) # FALSE, so the ASVs are in the same order

colnames(counts_filtered) <- paste0("ASV_", 1:ncol(counts_filtered))

# Construct a table to summarize the removal of reads throughout the pipeline
getN <- function(x) {sum(getUniques(x))}
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(counts_trimmed), rowSums(counts_nochim), rowSums(counts_filtered))
colnames(track) <- c("Input", "Filtered", "DenoisedF", "DenoisedR", "Merged", "Trimmed", "Non-Chimeric", "Tax-Filtered")
rownames(track) <- sample.names

# Output the counts, tax, and tracking tables
write.table(counts_filtered, 
            file = "C:/Users/ASUS/Desktop/prosthesis_microbiota_study/data/counts_gg.txt", 
            sep = "\t", col.names = NA, quote = F)

write.table(tax_filtered, 
            file = "C:/Users/ASUS/Desktop/prosthesis_microbiota_study/data/tax_gg.txt",
            sep = "\t", col.names = NA, quote = F)

write.table(track, 
            file = "C:/Users/ASUS/Desktop/prosthesis_microbiota_study/data/track_gg.txt",
            sep = "\t", col.names = NA, quote = F)
