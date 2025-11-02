library(ggplot2)
library(stringr)   


resdir <- "/Applications/Document/denovo_blast_trinity_own/GO_topGO"
setwd(resdir)


result_files <- list.files(resdir, pattern="_topGO.txt$", full.names=TRUE)


plot_topgo_dot <- function(file) {
  
  df <- read.table(file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
  
  
  df$classicFisher <- as.numeric(sub("<", "", df$classicFisher))
  df$classicFisher[is.na(df$classicFisher)] <- 1
  
  
  df$GeneRatio <- as.numeric(df$Significant) / as.numeric(df$Annotated)
  
  
  plot_data <- df[order(df$classicFisher), ][1:min(15, nrow(df)), ]
  
  
  plot_data$Term <- str_wrap(plot_data$Term, width=50)
  
  
  plotfile <- sub("_topGO.txt$", "_topGO_dotplot.png", file)
  
  
  p <- ggplot(plot_data, aes(x=GeneRatio, y=reorder(Term, GeneRatio))) +
    geom_point(aes(size=Significant, color=-log10(classicFisher))) +
    scale_color_gradient(low="blue", high="red") +
    labs(title=sub("_topGO.txt","",basename(file)),
         x="Gene Ratio", y="GO Term",
         color="-log10(p-value)", size="Gene Count") +
    theme_bw(base_size=10) +
    theme(
      axis.text.y = element_text(size=7),   
      axis.text.x = element_text(size=8)  
    ) +
    scale_y_discrete(expand = expansion(mult = c(0.01, 0.01)))  
  
  
  ggsave(plotfile, p, width=6, height=10)
}


for (f in result_files) {
  message("Dotplot for ", basename(f))
  plot_topgo_dot(f)
}
