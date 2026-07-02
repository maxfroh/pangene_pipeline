# if (!require("BiocManager", quietly = TRUE)) {
#   suppressMessages(install.packages("BiocManager"))
# }

# if (!require("DESeq2", quietly = TRUE)) {
#   suppressMessages(BiocManager::install("DESeq2"))
# }

suppressMessages(library("DESeq2"))
suppressMessages(library("tximport"))
suppressMessages(library("stringr"))
suppressMessages(library("readr"))
suppressMessages(library("tools"))

args = commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("You must call this file with the proper arguments.\nusage: Rscript deseq.R ...")
} else {
  dge_dir <- args[1]
  column_data_file <- args[2]
  annotation_file <- args[3]
  results_file <- args[4]
  counts_file <- args[5]
  abundance_file <- args[6]
  samples <- args[7:length(args)]
}

# Read in .gtf file and prepare it for use by tximport
# if (tools::file_ext(annotation_file) == "gtf") {
#   gtf <- read.csv(annotation_file, sep = "\t", header = FALSE)["V9"]
#   gtf$V9 <- stringr::str_replace(gtf$V9, ";$", "")
#   gtf[c("transcript_id", "Geneid")] <- stringr::str_split_fixed(gtf$V9, "; ", 2)
#   gtf$transcript_id <- stringr::str_replace(gtf$transcript_id, "transcript_id ", "")
#   gtf$Geneid <- stringr::str_replace(gtf$Geneid, "gene_id ", "")
# } else if (tools::file_ext(annotation_file) == ".reduced_map") {
#   gtf <- read.csv(annotation_file, sep = "\t", header = TRUE)
# } else {
#   stop("Improper file for transcriptome:gene mapping!")
# }
if (tools::file_ext(annotation_file) == "map") {
  tx2gene_map <- read.csv(annotation_file, sep = "\t", header = TRUE)
} else {
  stop("Improper file for transcriptome:gene mapping!")
}

# tx2gene will map transcript_ids to genes
tx2gene <- unique(tx2gene_map[c("transcript_id", "Geneid")])
tx2gene$transcript_id <- as.character(tx2gene$transcript_id)
tx2gene$Geneid <- as.character(tx2gene$Geneid)

# read in all kallisto abundance files (one abundance file per sample/replicate)
files <- file.path(dge_dir, samples, "abundance.tsv")
names(files) <- samples

# tximport
txi <- tximport::tximport(files, type = "kallisto", tx2gene = tx2gene, ignoreAfterBar = TRUE, countsFromAbundance = "no")

print(counts_file)
readr::write_tsv(as.data.frame(txi$counts), file=counts_file)
readr::write_tsv(as.data.frame(txi$abundance), file=abundance_file)

# get conditions so DESeq can group
colData <- read.csv(column_data_file, sep="\t", row.names=1)
colData$condition <- factor(colData$condition)

dds <- DESeq2::DESeqDataSetFromTximport(txi = txi, colData = colData, design = ~condition)

print("Performing DESeq")
dds <- DESeq2::DESeq(dds)
# resultsNames(dds)
res <- results(dds, independentFiltering=FALSE)
res$Geneid <- rownames(res) 
# Only want false discovery rate and gene id
resTrunc <- res[, c("Geneid", "padj")]

# Write truncated results
readr::write_tsv(as.data.frame(resTrunc), file=results_file)

# summary(res)

