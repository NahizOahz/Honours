suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(ggrepel)
})


indir  <- "/Applications/Document/denovo_blast_trinity_own/DEG_list"
outdir <- file.path("/Applications/Document/denovo_blast_trinity_own/DEG_plot")
if (!dir.exists(outdir)) dir.create(outdir)

files <- c(
  "DEG_ALE35_saline_annot.csv",
  "DEG_ALE100_saline_annot.csv",
  "DEG_strain_35ppt_annot.csv",
  "DEG_strain_100ppt_annot.csv"
)

padj_cutoff <- 0.05
lfc_cutoff  <- 1


for (f in files) {
  message("Processing: ", f)
  
  
  deg <- read.csv(file.path(indir, f))
  deg$padj_plot <- ifelse(is.na(deg$padj), 1, deg$padj)
  
  
  deg$deg_group <- "none"
  deg$deg_group[deg$log2FoldChange >  lfc_cutoff & deg$padj_plot < padj_cutoff] <- "up"
  deg$deg_group[deg$log2FoldChange < -lfc_cutoff & deg$padj_plot < padj_cutoff] <- "down"
  
  
  deg_sig <- subset(deg, deg_group %in% c("up","down"))
  
  
  genes_up   <- deg_sig$uniprot_acc[deg_sig$deg_group == "up"]
  genes_down <- deg_sig$uniprot_acc[deg_sig$deg_group == "down"]
  
  
  comp_name <- sub("\\.csv$", "", f)
  write.table(genes_up,   file=file.path(outdir, paste0(comp_name, "_up.txt")),
              row.names=FALSE, col.names=FALSE, quote=FALSE)
  write.table(genes_down, file=file.path(outdir, paste0(comp_name, "_down.txt")),
              row.names=FALSE, col.names=FALSE, quote=FALSE)
  
  
  
  p <- ggplot(deg, aes(x=log2FoldChange, y=-log10(padj_plot), color=deg_group)) +
    geom_point(alpha=0.7, size=1) +   
    scale_color_manual(values=c(up="red", down="blue", none="grey70")) +
    geom_vline(xintercept=c(-lfc_cutoff, lfc_cutoff), linetype="dashed", color="black") +
    geom_hline(yintercept=-log10(padj_cutoff), linetype="dashed", color="black") +
    theme_minimal(base_size=14) +
    labs(title=comp_name,
         x="log2 Fold Change", y="-log10(FDR)",
         color="Regulation") 
  # diy size# + coord_cartesian(xlim = c(-25, 30), ylim = c(0, 150))
  
  ggsave(file.path(outdir, paste0(comp_name, "_volcano.png")),
         plot=p, width=10, height=10, dpi=300)
}
