suppressPackageStartupMessages({
  library(tximport); library(DESeq2); library(readr); library(ggrepel)
})


base_dir <- "/Applications/Document/denovo_blast_trinity_own/all_quant"
tx2gene_path <- file.path(base_dir, "Trinity_rm.Trinity.fasta.gene_trans_map")
blast_path   <- "/Applications/Document/denovo_blast_trinity_own/gene_annotation/gene_level_annotation.csv"
outdir       <- "/Applications/Document/denovo_blast_trinity_own/DEG_list"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)


sample_names <- c(
  "ALE35_35_1","ALE35_35_2","ALE35_35_3",
  "ALE35_100_1","ALE35_100_2","ALE35_100_3",
  "ALE100_35_1","ALE100_35_2","ALE100_35_3",
  "ALE100_100_1","ALE100_100_2","ALE100_100_3"
)
files <- file.path(base_dir, paste0(sample_names, "_quant.sf"))
names(files) <- sample_names

samples <- data.frame(
  row.names = sample_names,
  strain   = factor(ifelse(grepl("^ALE35", sample_names), "ALE35", "ALE100")),
  salinity = factor(ifelse(grepl("_35_",  sample_names), "35ppt", "100ppt"))
)
samples$strain   <- relevel(samples$strain, ref="ALE35")
samples$salinity <- relevel(samples$salinity, ref="35ppt")


gt <- read.delim(tx2gene_path, header=FALSE, col.names=c("gene_id","transcript_id"))
tx2gene <- gt[, c("transcript_id","gene_id")]

txi <- tximport(
  files, type="salmon", tx2gene=tx2gene,
  countsFromAbundance="lengthScaledTPM"
)

dds <- DESeqDataSetFromTximport(txi, colData=samples, design = ~ strain * salinity)
dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds)

blast <- read.csv(blast_path)


contrast_list <- list(
  ALE35_saline  = list(contrast = c("salinity", "100ppt", "35ppt"),
                       label    = "ALE35 (100ppt vs 35ppt)"),
  ALE100_saline = list(contrast = list(c("salinity_100ppt_vs_35ppt",
                                         "strainALE100.salinity100ppt")),
                       label    = "ALE100 (100ppt vs 35ppt)"),
  strain_35ppt  = list(contrast = c("strain", "ALE100", "ALE35"),
                       label    = "ALE100 vs ALE35 (35ppt)"),
  strain_100ppt = list(contrast = list(c("strain_ALE100_vs_ALE35",
                                         "strainALE100.salinity100ppt")),
                       label    = "ALE100 vs ALE35 (100ppt)")
)


for (name in names(contrast_list)) {
  message("Processing: ", name, " — ", contrast_list[[name]]$label)
  res <- results(dds, contrast = contrast_list[[name]]$contrast)
  res_df <- as.data.frame(res)
  res_df$gene_id <- rownames(res_df)
  
  
  res_annot <- merge(res_df, blast, by.x="gene_id", by.y="gene_id", all.x=TRUE)
  
  
  outfile <- file.path(outdir, paste0("DEG_", name, "_annot.csv"))
  write.csv(res_annot, file=outfile, row.names=FALSE)
  message("✅ Saved: ", outfile)
}

message( outdir)


vsd <- vst(dds, blind = FALSE)
vsd_mat <- assay(vsd)
write.csv(vsd_mat, file=file.path(outdir, "VST_expression_matrix.csv"))

