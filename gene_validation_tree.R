
library(ape)       
library(ggtree)
library(ggplot2)
library(dplyr)


tree <- read.tree("~/Downloads/all_dnastring_18StufA.fasta (1).contree")
tree$tip.label <- gsub("_", " ", tree$tip.label)


tree_rooted <- root(tree, outgroup = "Chlamydomonas reinhardtii CC-1952", resolve.root = TRUE)
tree <- tree_rooted


nd <- Ntip(tree)
bt_tbl <- tibble(node = (nd + 1):(nd + tree$Nnode),
                 bootstrap = suppressWarnings(as.numeric(tree$node.label)))
tree_data <- ggtree(tree)$data %>% left_join(bt_tbl, by = "node")
tree_data$bootstrap[is.na(tree_data$bootstrap)] <- 0


p <- ggtree(tree, layout = "rectangular", size = 0.8) +
  geom_tiplab(aes(label = label), size = 5, fontface = "bold.italic", nudge_x = 0.01) +
  geom_nodepoint(aes(subset = !isTip), shape = 21, fill = "black", size = 1.5) +
  geom_text(data = subset(tree_data, !isTip & bootstrap > 60),
            aes(x = x, y = y, label = bootstrap),
            size = 3, hjust = -0.2, nudge_x = 0.001,
            color = "black", fontface = "bold") +
  theme_tree2() +
  geom_treescale(x = 0, y = 0.3, offset = 0.2, fontsize = 3) +
  theme(legend.position = "none") +
  xlim(0, max(tree_data$x) + 0.4) + 
  ggtitle("Phylogenetic Tree of 63 Microalgae (Rooted with Chlamydomonas reinhardtii)")

ggsave("/Applications/hpc_download/phylogenetic_tree_no_group.png",
       plot = p, width = 16, height = 14, dpi = 600)

print(p)
