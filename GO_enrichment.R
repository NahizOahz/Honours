suppressPackageStartupMessages({
  library(topGO)
})


workdir <- "/Applications/Document/denovo_blast_trinity_own/DEG_plot"
setwd(workdir)

bg_file <- "genelist.csv"


all_files <- list.files(workdir, pattern="^DEG_.*_(up|down)\\.txt$", full.names=TRUE)


g2go_files <- list(
  BP = "gene2go_BP.txt",
  MF = "gene2go_MF.txt",
  CC = "gene2go_CC.txt"
)


outdir <- "/Applications/Document/denovo_blast_trinity_own/GO_topGO"
if (!dir.exists(outdir)) dir.create(outdir, recursive=TRUE)


bg_genes <- read.csv(bg_file, header=FALSE, stringsAsFactors=FALSE)[,1]


run_topgo <- function(target_genes, prefix, ontology, g2go_file) {
  # geneList
  geneList <- factor(as.integer(bg_genes %in% target_genes))
  names(geneList) <- bg_genes
  
  
  gene2GO <- readMappings(file = g2go_file)
  
  
  GOdata <- new("topGOdata",
                ontology = ontology,
                allGenes = geneList,
                annot = annFUN.gene2GO,
                gene2GO = gene2GO)
  
  
  resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
  resultElim   <- runTest(GOdata, algorithm = "elim", statistic = "fisher")
  
  
  allRes <- GenTable(GOdata,
                     classicFisher = resultFisher,
                     elimFisher    = resultElim,
                     topNodes = length(score(resultElim)))
  
  
  allRes$Term <- Term(GOTERM[allRes$GO.ID])
  
  
  allRes$classicFisher <- as.numeric(allRes$classicFisher)
  allRes$elimFisher    <- as.numeric(allRes$elimFisher)
  
  
  allRes$classicFisherAdj <- p.adjust(allRes$classicFisher, method="BH")
  allRes$elimFisherAdj    <- p.adjust(allRes$elimFisher, method="BH")
  
  
  outfile_full <- file.path(outdir, paste0(prefix, "_", ontology, "_topGO_full.txt"))
  write.table(allRes, file=outfile_full, sep="\t", quote=FALSE, row.names=FALSE)
  
  
  sigRes <- subset(allRes, elimFisherAdj < 0.05)
  outfile_sig <- file.path(outdir, paste0(prefix, "_", ontology, "_topGO_sig.txt"))
  write.table(sigRes, file=outfile_sig, sep="\t", quote=FALSE, row.names=FALSE)
  
  return(list(full=allRes, sig=sigRes))
}
 


for (f in all_files) {
  fname <- basename(f)
  cond  <- sub("_(up|down)\\.txt$", "", fname)
  dir   <- sub("^DEG_.*_(up|down)\\.txt$", "\\1", fname)
  
  target_genes <- read.csv(f, header=FALSE, stringsAsFactors=FALSE)[,1]
  
  for (ont in names(g2go_files)) {
    message("Processing ", cond, "_", dir, " with ontology ", ont, "...")
    run_topgo(target_genes, paste0(cond, "_", dir), ont, g2go_files[[ont]])
  }
}