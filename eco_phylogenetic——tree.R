
library(ggtree)
library(ggplot2)
library(dplyr)


tree <- read.tree("~/Downloads/multigenes_63.fasta.contree")


tree$tip.label <- gsub("_", " ", tree$tip.label)


groups <- data.frame(
  label = c(
    # Chlorophyta
    "Auxenochlorella protothecoides", "Auxenochlorella pyrenoidosa",
    "Botryococcus braunii", "Chlamydomonas moewusii", "Chlamydomonas nivalis",
    "Chlamydomonas reinhardtii", "Chlamydomonas sp. ICE-L",
    "Chlorella sorokiniana", "Chlorella variabilis", "Chlorella vulgaris",
    "Chlorella zofingiensis", "Choricystis parasitica", "Coccomyxa subellipsoidea",
    "Coelastrum microporum", "Dunaliella salina", "Haematococcus lacustris", "Lobosphaera incisa",
    "Micractinium conductrix", "Micractinium pusillum", "Mychonastes homosphaera",
    "Neochloris aquatica", "Parachlorella kessleri", "Pediastrum duplex",
    "Pseudochloris wilhelmii", "Tetradesmus obliquus", "Trebouxia aggregata",
    "Tetraselmis suecica", "Volvox carteri",
    "Bathycoccus prasinos", "Ostreococcus tauri", "Pyramimonas parkeae",
    
    # Streptophyta
    "Klebsormidium flaccidum", "Klebsormidium nitens",
    
    # Bacillariophyta
    "Chaetoceros muelleri", "Cylindrotheca closterium", "Cyclotella meneghiniana",
    "Fragilariopsis cylindrus", "Fistulifera solaris", "Halamphora coffeaeformis",
    "Leptocylindrus danicus", "Nitzschia palea", "Odontella aurita",
    "Phaeodactylum tricornutum", "Pseudo-nitzschia multiseries",
    "Skeletonema marinoi", "Skeletonema costatum",
    "Thalassiosira oceanica", "Thalassiosira pseudonana",
    
    # Eustigmatophyceae
    "Nannochloropsis gaditana", "Nannochloropsis limnetica",
    "Nannochloropsis oceanica", "Nannochloropsis salina",
    "Trachydiscus minutus",
    
    # Rhodophyta
    "Porphyridium purpureum", "Cyanidioschyzon merolae",
    
    # Cryptophyta
    "Rhodomonas salina",
    
    # Haptophyta
    "Emiliania huxleyi", "Isochrysis galbana",
    "Chrysochromulina tobin", "Diacronema lutheri", "Tisochrysis lutea",
    
    # Raphidophyceae
    "Heterosigma akashiwo",
    
    # Glaucophyta
    "Cyanophora paradoxa"
  ),
  group = c(
    rep("Chlorophyta", 31),
    rep("Streptophyta", 2),
    rep("Bacillariophyta", 15),
    rep("Eustigmatophyceae", 5),
    rep("Rhodophyta", 2),
    "Cryptophyta",
    rep("Haptophyta", 5),
    "Raphidophyceae",
    "Glaucophyta"
  )
)


setdiff(tree$tip.label, groups$label)


tree_data <- ggtree(tree)$data
tree_data$bootstrap <- as.numeric(tree_data$label)
tree_data$bootstrap[is.na(tree_data$bootstrap)] <- 0


p <- ggtree(tree, layout = "rectangular", size = 0.8) %<+% groups +
  
  geom_tiplab(aes(label = label, color = group),
              size = 5, fontface = "bold.italic", nudge_x = 0.01) +
  
  geom_nodepoint(aes(subset = !isTip), shape = 21, fill = "black", size = 1.5) +
  # bootstrap > 60 
  geom_text(data = subset(tree_data, !isTip & bootstrap > 60),
            aes(x = x, y = y, label = bootstrap),
            size = 3, hjust = -0.2, nudge_x = 0.001,
            color = "black", fontface = "bold") +
  
  geom_point(aes(x = 0, y = 0, color = group), size = 8, shape = 15) +
  scale_color_manual(
    name = "Phylum",
    values = c(
      "Chlorophyta" = "#006600",
      "Streptophyta" = "#9ACD32",
      "Bacillariophyta" = "#7341D8",
      "Eustigmatophyceae" = "#FF8C00",
      "Rhodophyta" = "#FF0000",
      "Cryptophyta" = "#FF007F",
      "Haptophyta" = "#1E90FF",
      "Raphidophyceae" = "#800080",
      "Glaucophyta" = "#000000"
    )
  ) +
  theme_tree2() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12, face = "bold"),
    axis.line.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  xlim(0, max(tree_data$x) + 0.12) +
  geom_treescale(x = 0, y = -0.1, offset = 0.2, fontsize = 3) +
  ggtitle("Phylogenetic Tree of 63 Microalgae (Phylum-level)")


ggsave("/Applications/hpc_download/phylogenetic_tree_phylum_colored_bootstrap.png",
       plot = p, width = 16, height = 14, dpi = 600)

print(p)