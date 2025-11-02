library(ape)
library(ips)

folder_in <- list.files('/Users/zhao/Library/CloudStorage/OneDrive-UTS/Zihan_project/phylogenetic_tree')
folder_in
allgenes <- unique(gsub(".*__|\\.fasta", "",list.files('/Users/zhao/Library/CloudStorage/OneDrive-UTS/Zihan_project/phylogenetic_tree',
                                                       recursive = T, full.names = T)))
allgenes
####step 1 align each gene across species
for(i in 1:length(allgenes)){
  #i = 1}
  gene_goi <- allgenes[i]
  allfasta_goi <- list.files('/Users/zhao/Library/CloudStorage/OneDrive-UTS/Zihan_project/phylogenetic_tree',
                             pattern = gene_goi,
                             recursive = T, full.names = T)
  allfasta_goi
  dna_in <- c()
  for(j in 1:length(allfasta_goi)){
    fasta_in <- allfasta_goi[j]
    seqin <- ape::read.FASTA(fasta_in)
    if(is.null(dna_in)){
      dna_in <- seqin
    }else{
      dna_in <- c(dna_in, seqin)  
    }
    
    #fasta_in1 <- ape::read.FASTA("/Users/zhao/Library/CloudStorage/OneDrive-UTS/Zihan_project/phylogenetic_tree//Volvox_carteri/Volvox_carteri__18SrRNA.fasta")
    #fasta_in2 <- ape::read.FASTA("/Users/zhao/Library/CloudStorage/OneDrive-UTS/Zihan_project/phylogenetic_tree//Scenedesmus_obliquus/Scenedesmus_obliquus__18SrRNA.fasta")
    
    #dna_in <- c(fasta_in1, fasta_in2)
    
  }
  aligned_dna <- ips::mafft(dna_in)
  
  # Assuming my_dnabin is your DNAbin object
  ape::write.FASTA(aligned_dna, 
                   paste0('/Users/zhao/Library/CloudStorage/OneDrive-UTS/Zihan_project/aligned_sequences/', gene_goi, ".fasta"))
  
  aligned_in <- read.FASTA(paste0('/Users/zhao/Library/CloudStorage/OneDrive-UTS/Zihan_project/aligned_sequences/', gene_goi, ".fasta"))
  
  out_fasta <- aligned_in[["Botryococcus_braunii__18SrRNA"]]
  aligned_in$Botryococcus_braunii__18SrRNA
  
  numeric_to_base <- c("04" = "-", "f0" = "-", "88" = "A", "18" = "T", "28" = "C", "48" = "G", "90" = "W", "e0" = "V", "30" = "Y", "c0" = "R", "60" = "S", "a0" = "M", "50" = "K", "70" = "B", "d0" = "D", "b0" = "H")
  for(k in seq_along(aligned_in)){
    #k = 1}
    seq_name <- names(aligned_in)[k]  
    seq_data <- paste(sapply(aligned_in[[k]], function(num) numeric_to_base[as.character(num)]), collapse = "")
    #seq_data <- paste(sapply(aligned_in[[k]], function(num) numeric_to_base[as.character(num)]), collapse = "")
    file_name <- paste0("/Users/zhao/Library/CloudStorage/OneDrive-UTS/Zihan_project/aligned_sequences/split/",seq_name, ".fasta")
    
    cat(">", seq_name, "\n", paste(seq_data, collapse = ""), "\n", file = file_name, sep = "")
    
  }
  
}

####concatenate by each species

library(ape)
species_names <- folder_in
#species_names: "Botryococcus_braunii"

# choose one specie
for(l in 1:length(species_names)){
  #l = 1
  species_in <- species_names[l]
  
  allfasta_aligned <- list.files('/Users/zhao/Library/CloudStorage/OneDrive-UTS/Zihan_project/aligned_sequences/split',
                                 pattern = species_in, 
                                 recursive = T, full.names = T)
  allfasta_aligned
  
  # merge each genes under the chosen specie
  species_sequences <- c()
  for(m in 1:length(allfasta_aligned) ){
    #m = 1
    aligned_fasta_in <-  allfasta_aligned[m]
    aligned_seqin <- ape::read.FASTA(aligned_fasta_in)
    
    # try remove space
    aligned_seqin_char <- unlist(as.character(aligned_seqin))
    species_sequences <- c(species_sequences, aligned_seqin_char)
    
    
  }
  species_sequences <- list(species_sequences)
  names(species_sequences) <- species_in
  dnabin_data <- as.DNAbin(species_sequences)
  class(dnabin_data)
  
  #save file and path
  ape::write.FASTA(dnabin_data, paste0('/Users/zhao/Library/CloudStorage/OneDrive-UTS/Zihan_project/concatenated/', species_in, ".fasta")) 
  
}

  #merge each species
system('cat /Users/zhao/Library/CloudStorage/OneDrive-UTS/Zihan_project/concatenated/*.fasta > /Users/zhao/Library/CloudStorage/OneDrive-UTS/Zihan_project/concatenated/multigenes_65.fasta')
