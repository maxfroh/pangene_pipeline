if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

if (!require("DESeq2", quietly = TRUE)) {
  BiocManager::install("DESeq2")
}

suppressMessages(library("DESeq2"))
library("tximport")
library("stringr")

args = commandArgs(trailingOnly = TRUE)
if (FALSE & length(args) < 10) {
  stop("You must call this file with the proper arguments.\nusage: Rscript deseq.R ...")
} else {
  annotation_file <- args[-1]
  alpha <- args[-1]
  lf2c_thresh <- args[-1]
  samples <- args[4:length(args)]
}

# Read in .gtf file and prepare it for use by tximport
gtf <- read.csv(annotation_file, sep = "\t", header = FALSE)["V9"]
gtf$V9 <- stringr::str_replace(gtf$V9, ";$", "")
gtf[c("transcript_id", "Geneid")] <- stringr::str_split_fixed(gtf$V9, "; ", 2)
gtf$transcript_id <- stringr::str_replace(gtf$transcript_id, "transcript_id ", "")
gtf$Geneid <- stringr::str_replace(gtf$Geneid, "gene_id ", "")

# tx2gene will map transcript_ids to genes
tx2gene <- unique(gtf[c("transcript_id", "Geneid")])

# read in all kallisto abundance files (one abundance file per sample/replicate)
files <- file.path(kallisto_dir, samples, "abundance.tsv")
names(files) <- samples

# tximport
txi <- tximport::tximport(files, type = "kallisto", tx2gene = tx2gene, ignoreAfterBar = TRUE, countsFromAbundance = "no")

sampleTable <- data.frame(condition = factor(rep(c("A", "B"), each = 3)))
rownames(sampleTable) <- colnames(txi$counts)

head(txi$counts)
colnames(txi$counts)
sampleTable

# DESeq2::DESeqDatasetFromTximport(txi = txi, colData = sampleTable, design = ~condition)