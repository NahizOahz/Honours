suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(matrixStats)
})


vsd_csv <- "/Applications/Document/denovo_blast_trinity_own/DEG_list/VST_expression_matrix.csv"
outdir  <- "/Applications/Document/denovo_blast_trinity_own/DEG_list"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)


vsd_mat <- read.csv(vsd_csv, row.names = 1, check.names = FALSE)
stopifnot(ncol(vsd_mat) >= 3)
stopifnot(nrow(vsd_mat) >= 2)


keep <- apply(vsd_mat, 1, function(v) stats::var(v, na.rm = TRUE) > 0)
vsd_mat <- vsd_mat[keep, , drop = FALSE]


ntop <- 500
rv <- rowVars(as.matrix(vsd_mat))
select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
mat_sel <- vsd_mat[select, , drop = FALSE]


sample_names <- colnames(vsd_mat)
samples <- data.frame(
  sample   = sample_names,
  Strain   = ifelse(grepl("^WT", sample_names), "WT", "ALE100"),
  Salinity = ifelse(grepl("_35_",  sample_names), "35ppt", "100ppt"),
  stringsAsFactors = FALSE
)


pca <- prcomp(t(mat_sel), center = TRUE, scale. = FALSE)
pc_scores <- as.data.frame(pca$x[, 1:2, drop = FALSE])
pc_scores$sample <- rownames(pc_scores)
df <- merge(pc_scores, samples, by = "sample", all.x = TRUE)


varExpl <- 100 * (pca$sdev^2 / sum(pca$sdev^2))
xlab_txt <- sprintf("PC1: %.1f%% variance", varExpl[1])
ylab_txt <- sprintf("PC2: %.1f%% variance", varExpl[2])


theme_base <- theme_bw(base_size = 14) +
  theme(panel.grid = element_blank(),
        legend.position = "right",
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5))



lim_range <- range(c(df$PC1, df$PC2))
margin <- diff(lim_range) * 0.05             
lims <- c(lim_range[1] - margin, lim_range[2] + margin)  


p_basic <- ggplot(df, aes(PC1, PC2, color = Salinity, shape = Strain)) +
  geom_point(size = 4, alpha = 0.9, stroke = 0.4) +
  xlab(xlab_txt) + ylab(ylab_txt) +
  coord_equal(xlim = lims, ylim = lims) +    
  theme_base +
  ggtitle("PCA of RNA-seq Samples")


ggsave(file.path(outdir, "PCA_standard_basic_square.png"),
       plot = p_basic, width = 7, height = 7, dpi = 600)


set.seed(1)
p_label <- p_basic +
  geom_text_repel(aes(label = sample), size = 3.5, max.overlaps = 50)
ggsave(file.path(outdir, "PCA_standard_label.png"),
       plot = p_label, width = 7, height = 7, dpi = 600)


p_ellipse <- p_basic +
  stat_ellipse(aes(group = interaction(strain, salinity), color = strain),
               level = 0.65, linewidth = 0.5, linetype = 2, alpha = 2)
ggsave(file.path(outdir, "PCA_standard_ellipse.png"),
       plot = p_ellipse, width = 7, height = 7, dpi = 600)


write.csv(df, file = file.path(outdir, "PCA_scores_standard.csv"), row.names = FALSE)
write.csv(
  data.frame(PC = paste0("PC", seq_along(varExpl)), variance_percent = varExpl),
  file = file.path(outdir, "PCA_variance_explained_standard.csv"),
  row.names = FALSE
)

message("✅ntop=500, scale.=FALSE)")