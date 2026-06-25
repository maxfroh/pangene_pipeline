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

args = commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("You must call this file with the proper arguments.\nusage: Rscript deseq.R ...")
} else {
  kallisto_dir <- args[1]
  column_data_file <- args[2]
  annotation_file <- args[3]
  results_file <- args[4]
  alpha <- args[5]
  lf2c_thresh <- args[6]
  samples <- args[7:length(args)]
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

# sampleTable <- data.frame(condition = factor(rep(c("A", "B"), each = 3)))
# rownames(sampleTable) <- colnames(txi$counts)
colData <- read.csv(column_data_file, sep="\t", row.names=1)
colData$condition <- factor(colData$condition)

head(txi$counts)
colnames(txi$counts)

colData


dds <- DESeq2::DESeqDataSetFromTximport(txi = txi, colData = colData, design = ~condition)

print("Performing DESeq")
dds <- DESeq2::DESeq(dds)
resultsNames(dds)
# res <- results(dds, alpha=alpha)
res <- results(dds, independentFiltering=FALSE)
# table(res$padj == is.na(res$padj) & res$baseMean > 0)
resLFC <- lfcShrink(dds, coef=2, type="apeglm")
# resOrdered <- res[order(res$pvalue),]
# resSig <- subset(resOrdered, padj < alpha & abs(log2FoldChange) >= 1)
res$Geneid <- rownames(res) 
resTrunc <- res[, c("Geneid", "padj")]
# rownames(res)

readr::write_tsv(as.data.frame(res), file=str_replace(results_file, ".tsv", "_untrunc.tsv"))
readr::write_tsv(as.data.frame(resTrunc), file=results_file)

summary(res)

