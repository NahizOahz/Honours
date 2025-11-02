suppressPackageStartupMessages({
  library(dplyr)
})


base_dir <- "/Applications/Document/denovo_blast_trinity_own"
blast_path <- file.path(base_dir, "blast_filter_transcripts/all_hits_filtered.csv")
map_path   <- file.path(base_dir, "all_quant/Trinity_rm.Trinity.fasta.gene_trans_map")


blast <- read.csv(blast_path, header=TRUE, check.names=FALSE)


gt <- read.delim(map_path, header=FALSE,
                 col.names=c("gene_id","transcript_id"),
                 stringsAsFactors=FALSE)


blast$qseqid_clean <- sub("\\.p[0-9]+$", "", blast$qseqid)


blast_gene <- merge(blast, gt, by.x="qseqid_clean", by.y="transcript_id")


blast_gene_best <- blast_gene %>%
  group_by(gene_id) %>%
  arrange(evalue, desc(bitscore), desc(pident)) %>%
  slice(1) %>%
  ungroup()


blast_gene_best <- blast_gene_best %>%
  mutate(
    uniprot_acc    = sub("^.+\\|(\\w+)\\|.*$", "\\1", sseqid),   # Q9M841
    gene_symbol    = sub("^.+\\|.+\\|(.*)$", "\\1", sseqid)      # PEX12_ARATH
  )


outdir <- file.path(base_dir, "gene_annotation")
if (!dir.exists(outdir)) dir.create(outdir)

write.csv(blast_gene_best,
          file=file.path(outdir, "gene_level_annotation.csv"),
          row.names=FALSE)