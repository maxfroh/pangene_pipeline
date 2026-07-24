# if (!require("BiocManager", quietly = TRUE)) {
#   suppressMessages(install.packages("BiocManager"))
# }

# if (!require("DESeq2", quietly = TRUE)) {
#   suppressMessages(BiocManager::install("DESeq2"))
# }

suppressMessages(library("DESeq2"))
suppressMessages(library("tximport"))
suppressMessages(library("stringr"))
suppressMessages(library("dplyr"))
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

print(results_file)

if (tools::file_ext(annotation_file) == "map") {
  tx2gene_map <- read.csv(annotation_file, sep = "\t", header = TRUE)
} else {
  stop("Improper file for transcriptome to gene mapping!")
}

# tx2gene will map transcript_ids to genes
tx2gene <- unique(tx2gene_map[c("transcript_id", "Geneid")])
tx2gene$transcript_id <- as.character(tx2gene$transcript_id)
tx2gene$Geneid <- as.character(tx2gene$Geneid)

# read in all kallisto abundance files (one abundance file per sample/replicate)
files <- file.path(dge_dir, samples, "abundance.tsv")
names(files) <- samples

# tximport
print("Using tximport to convert transcript-level results to gene-level.")
txi <- tximport::tximport(files, type = "kallisto", tx2gene = tx2gene, ignoreAfterBar = TRUE, countsFromAbundance = "no")

counts_data <- as.data.frame(txi$counts)
counts_data$Geneid <- rownames(counts_data)
counts_data <- counts_data %>% dplyr::relocate("Geneid")
abundance_data <- as.data.frame(txi$abundance)
abundance_data$Geneid <- rownames(abundance_data)
abundance_data <- abundance_data %>% dplyr::relocate("Geneid")

readr::write_tsv(counts_data, file=counts_file)
readr::write_tsv(abundance_data, file=abundance_file)
print("Successfully ran tximport!")

# get conditions so DESeq can group
colData <- read.csv(column_data_file, sep="\t", row.names=1)
colData$condition <- factor(colData$condition)

dds <- DESeq2::DESeqDataSetFromTximport(txi = txi, colData = colData, design = ~condition)

print("Usign DESeq to calculate false discovery rate for differential expression analysis.")
dds <- DESeq2::DESeq(dds)
# resultsNames(dds)
res <- results(dds, independentFiltering=FALSE)
res$Geneid <- rownames(res) 
# Only want false discovery rate and gene id
resTrunc <- res[, c("Geneid", "padj")]

# Write truncated results
readr::write_tsv(as.data.frame(resTrunc), file=results_file)
print("Ran DESeq2 successfully!")
# summary(res)

