# Candidate SNP overlap replication
# Create a log file to store outputs
log_file <- "Indolics_Sig.Snp.Overlap.rep_log.txt"

# Open the log file for writting 
sink(log_file, append = TRUE)

### INDOLICS
# subset candidate genes to those affecting indolic GSL biosynthesis
candidate.genes.ind = subset(candidate.genes, indolic == "yes")
nrow(candidate.genes.ind) # 17 genes

candidate.snps.ind = vector()

# Create a list of file paths
file_paths <- list(
  "output_200/200_gluc14_maf0.03_miss0.1GF.assoc.txt",
  "output_200/200_gluc15_maf0.03_miss0.1GF.assoc.txt",
  "output_200/200_gluc16_maf0.03_miss0.1GF.assoc.txt",
  "output_200/200_gluc17_maf0.03_miss0.1GF.assoc.txt"
)

Indolic_List <- list(
  "IM",
  "1moIM",
  "1hIM",
  "4moIM"
)

for (j in seq_along(file_paths)) {
  file <- file_paths[[j]]
  glucosinolate <- Indolic_List[[j]]
  
  # Read the file and filter for maf > 0.05 and p_wald <0.05
  snps <- subset(read.delim(file, sep="\t")[,c("chr","ps","af","p_wald")], af > 0.05)
  sig.snps <- subset(snps, p_wald < 0.05)
  
  # Standardize SNP identifier as "chr_position"
  sig.snps$rs <- paste(sig.snps$chr, sig.snps$ps, sep = "_")
  
  # Remove duplicates 
  sig.snps <- unique(sig.snps)
  
  # get SNPs within boundaries specified for each candidate gene
  for (i in 1:nrow(candidate.genes.ind)){
    
    hits = subset(sig.snps, chr == candidate.genes.ind[i,"chr"] & 
                    ps > candidate.genes.ind[i,"ps_left"] &
                    ps < candidate.genes.ind[i,"ps_right"]
    )
    
    candidate.snps.ind = c(candidate.snps.ind, as.vector(hits$rs))
    
    print(paste0("Glucosinolate: ", glucosinolate))
    print(paste0("gene: ", candidate.genes.ind[i,"gene_abbreviation"], " - #", i, " of ", nrow(candidate.genes.ind)))
    print(paste0("hits: ", nrow(hits), " SNPs"))
    print(head(hits))
    
  }
  
}

sink()

# Candidate SNP overlap with 2022 outputs
# Create a log file to store outputs
log_file <- "Indolics_Sig.Snp.Overlap.2022_log.txt"

# Open the log file for writting 
sink(log_file, append = TRUE)

### INDOLICS
# subset candidate genes to those affecting indolic GSL biosynthesis
candidate.genes.ind = subset(candidate.genes, indolic == "yes")
nrow(candidate.genes.ind) # 17 genes

candidate.snps.ind = vector()

# Create a list of file paths
file_paths <- list(
  "tou_miss10/blup_indolic_assoc/lmm_maf03_BlupGIM.assoc.fixed.txt",
  "tou_miss10/blup_indolic_assoc/lmm_maf03_BlupG1MOIM.assoc.fixed.txt",
  "tou_miss10/blup_indolic_assoc/lmm_maf03_BlupG4HIM.assoc.fixed.txt",
  "tou_miss10/blup_indolic_assoc/lmm_maf03_BlupG4MOIM.assoc.fixed.txt"
)

Indolic_List <- list(
  "IM",
  "1moIM",
  "1hIM",
  "4moIM"
)

for (j in seq_along(file_paths)) {
  file <- file_paths[[j]]
  glucosinolate <- Indolic_List[[j]]
  
  # Read the file and filter for maf > 0.05 and p_wald <0.05
  snps <- subset(read.delim(file, sep="\t")[,c("chr","ps","af","p_wald")], af > 0.05)
  sig.snps <- subset(snps, p_wald < 0.05)
  
  # Standardize SNP identifier as "chr_position"
  sig.snps$rs <- paste(sig.snps$chr, sig.snps$ps, sep = "_")
  
  # Remove duplicates 
  sig.snps <- unique(sig.snps)
  
  # get SNPs within boundaries specified for each candidate gene
  for (i in 1:nrow(candidate.genes.ind)){
    
    hits = subset(sig.snps, chr == candidate.genes.ind[i,"chr"] & 
                    ps > candidate.genes.ind[i,"ps_left"] &
                    ps < candidate.genes.ind[i,"ps_right"]
    )
    
    candidate.snps.ind = c(candidate.snps.ind, as.vector(hits$rs))
    
    print(paste0("Glucosinolate: ", glucosinolate))
    print(paste0("gene: ", candidate.genes.ind[i,"gene_abbreviation"], " - #", i, " of ", nrow(candidate.genes.ind)))
    print(paste0("hits: ", nrow(hits), " SNPs"))
    print(head(hits))
    
  }
  
}

sink()

# ALIPHATIC
# Candidate SNP overlap with replication outputs
# Create a log file to store outputs
log_file <- "Aliphatics_Sig.Snp.Overlap.rep_log.txt"

# Open the log file for writting 
sink(log_file, append = TRUE)

### INDOLICS
# subset candidate genes to those affecting indolic GSL biosynthesis
candidate.genes.ali = subset(candidate.genes, aliphatic == "yes")
nrow(candidate.genes.ali) # 33 genes

candidate.snps.ali = vector()

# Create a list of file paths
file_paths <- list(
  "output_200/200_gluc1_maf0.03_miss0.1GF.assoc.txt", #1
  "output_200/200_gluc2_maf0.03_miss0.1GF.assoc.txt", #2
  "output_200/200_gluc3_maf0.03_miss0.1GF.assoc.txt", #3
  "output_200/200_gluc4_maf0.03_miss0.1GF.assoc.txt", #4
  "output_200/200_gluc5_maf0.03_miss0.1GF.assoc.txt", #5
  "output_200/200_gluc6_maf0.03_miss0.1GF.assoc.txt", #6
  "output_200/200_gluc7_maf0.03_miss0.1GF.assoc.txt", #7
  "output_200/200_gluc8_maf0.03_miss0.1GF.assoc.txt",#8
  "output_200/200_gluc9_maf0.03_miss0.1GF.assoc.txt", #9
  "output_200/200_gluc10_maf0.03_miss0.1GF.assoc.txt", #10
  "output_200/200_gluc11_maf0.03_miss0.1GF.assoc.txt", #11
  "output_200/200_gluc12_maf0.03_miss0.1GF.assoc.txt", #12
  "output_200/200_gluc13_maf0.03_miss0.1GF.assoc.txt" #13
)

Aliphatic_List <- list(
  "7mSh",
  "8MTO",
  "3mSOp",
  "4mSOb",
  "5mSOp",
  "6mSOh",
  "7mSOh",
  "8mSOo",
  "Pren",
  "Buen",
  "Peen",
  "S2hBuen",
  "2hPeen"
)

for (j in seq_along(file_paths)) {
  file <- file_paths[[j]]
  glucosinolate <- Aliphatic_List[[j]]
  
  # Read the file and filter for maf > 0.05 and p_wald <0.05
  snps <- subset(read.delim(file, sep="\t")[,c("chr","ps","af","p_wald")], af > 0.05)
  sig.snps <- subset(snps, p_wald < 0.05)
  
  # Standardize SNP identifier as "chr_position"
  sig.snps$rs <- paste(sig.snps$chr, sig.snps$ps, sep = "_")
  
  # Remove duplicates 
  sig.snps <- unique(sig.snps)
  
  # get SNPs within boundaries specified for each candidate gene
  for (i in 1:nrow(candidate.genes.ali)){
    
    hits = subset(sig.snps, chr == candidate.genes.ali[i,"chr"] & 
                    ps > candidate.genes.ali[i,"ps_left"] &
                    ps < candidate.genes.ali[i,"ps_right"]
    )
    
    candidate.snps.ali = c(candidate.snps.ali, as.vector(hits$rs))
    
    print(paste0("Glucosinolate: ", glucosinolate))
    print(paste0("gene: ", candidate.genes.ali[i,"gene_abbreviation"], " - #", i, " of ", nrow(candidate.genes.ali)))
    print(paste0("hits: ", nrow(hits), " SNPs"))
    #print(head(hits))
    
  }
  
}

sink()

# ALIPHATIC
# Candidate SNP overlap with 2022 outputs
# Create a log file to store outputs
log_file <- "Aliphatics_Sig.Snp.Overlap.2022_log.txt"

# Open the log file for writting 
sink(log_file, append = TRUE)

### INDOLICS
# subset candidate genes to those affecting indolic GSL biosynthesis
candidate.genes.ali = subset(candidate.genes, aliphatic == "yes")
nrow(candidate.genes.ali) # 33 genes

candidate.snps.ali = vector()

# Create a list of file paths
file_paths <- list(
  "tou_miss10/blup_aliphatic_assoc/lmm_maf03_BlupG7MTH.assoc.fixed.txt",
  "tou_miss10/blup_aliphatic_assoc/lmm_maf03_BlupG8MTO.assoc.fixed.txt",
  "tou_miss10/blup_aliphatic_assoc/lmm_maf03_BlupG3MSP.assoc.fixed.txt",
  "tou_miss10/blup_aliphatic_assoc/lmm_maf03_BlupG4MSB.assoc.fixed.txt",
  "tou_miss10/blup_aliphatic_assoc/lmm_maf03_BlupG5MSP.assoc.fixed.txt",
  "tou_miss10/blup_aliphatic_assoc/lmm_maf03_BlupG6MSH.assoc.fixed.txt",
  "tou_miss10/blup_aliphatic_assoc/lmm_maf03_BlupG7MSH.assoc.fixed.txt",
  "tou_miss10/blup_aliphatic_assoc/lmm_maf03_BlupG8MSO.assoc.fixed.txt",
  "tou_miss10/blup_aliphatic_assoc/lmm_maf03_BlupG2P.assoc.fixed.txt",
  "tou_miss10/blup_aliphatic_assoc/lmm_maf03_BlupG3B.assoc.fixed.txt",
  "tou_miss10/blup_aliphatic_assoc/lmm_maf03_BlupG4P.assoc.fixed.txt",
  "tou_miss10/blup_aliphatic_assoc/lmm_maf03_BlupG2H3B.assoc.fixed.txt",
  "tou_miss10/blup_aliphatic_assoc/lmm_maf03_BlupG2H4P.assoc.fixed.txt"
)

Aliphatic_List <- list(
  "7mSh",
  "8MTO",
  "3mSOp",
  "4mSOb",
  "5mSOp",
  "6mSOh",
  "7mSOh",
  "8mSOo",
  "Pren",
  "Buen",
  "Peen",
  "S2hBuen",
  "2hPeen"
)

for (j in seq_along(file_paths)) {
  file <- file_paths[[j]]
  glucosinolate <- Aliphatic_List[[j]]
  
  # Read the file and filter for maf > 0.05 and p_wald <0.05
  snps <- subset(read.delim(file, sep="\t")[,c("chr","ps","af","p_wald")], af > 0.05)
  sig.snps <- subset(snps, p_wald < 0.05)
  
  # Standardize SNP identifier as "chr_position"
  sig.snps$rs <- paste(sig.snps$chr, sig.snps$ps, sep = "_")
  
  # Remove duplicates 
  sig.snps <- unique(sig.snps)
  
  # get SNPs within boundaries specified for each candidate gene
  for (i in 1:nrow(candidate.genes.ali)){
    
    hits = subset(sig.snps, chr == candidate.genes.ali[i,"chr"] & 
                    ps > candidate.genes.ali[i,"ps_left"] &
                    ps < candidate.genes.ali[i,"ps_right"]
    )
    
    candidate.snps.ali = c(candidate.snps.ali, as.vector(hits$rs))
    
    print(paste0("Glucosinolate: ", glucosinolate))
    print(paste0("gene: ", candidate.genes.ali[i,"gene_id"], " - #", i, " of ", nrow(candidate.genes.ali)))
    print(paste0("hits: ", nrow(hits), " SNPs"))
    #print(head(hits))
    
  }
  
}

sink()