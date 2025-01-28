# GWAS_200_Replication

## Description
This project is a replication of the Genome Wide Association Study of glucosinolate abundance phenotypes in Arabidopsis Thaliana carried out in Gloss et al 2022 (https://pmc.ncbi.nlm.nih.gov/articles/PMC9149790/). Specifically, GWAS of individual phenotypes, finding the intersect of significant SNPs between the replication and the original study, determining which and how many SNPs overlap with glucosinolate biosynthetic candidate gene regions, and generating Manhattan plots to visualize gene associations. 

## Table of Contents
  -Inputs
  -phen_glucosinolatesBLUPs_2021-08-02.csv:
    -phenotype input: contains best linear unbiased predictions (BLUPs) generated from raw spectrometry data. Generating BLUPs from raw data accounts for noise that can result from variability          between samples, spectrometers, time, batch effect, and randomness. The BLUPs represent abundances of glucosinolate metabolytes within ~200 accessions and are processed to be the phenotype         input to GEMMA GWAS mixed linear modeling. 
  -TOU_accessions.csv:
    -Contains ecotype IDs to map to the glucosinolate dataset for input to GEMMA GWAS mixed linear modeling
  -genotypes: 
    -the genotypes file contains three binary files genotypes.bim, genotypes.bam, and genotypes.fam (they are too large to import to this repository). The fam file contains genotype IDs that are        mapped to the ecotype IDs for each accession.
-Scripts
  -ScriptPrepforGemmaRep200.R:
    -Transforming BLUPs dataset for GEMMA formatting. Removes column headers, merges and maps BLUPs to ecotype ID (TOU_accessions.csv) and genotype ID (genotypes.fam)
  -ScriptGEMMAKin:lmm.R:
    -Generate kinship matrix to account for population structure and random effects then run GWAS mixed linear model per phenotype. 
    -Outputs:
      -____.log.txt contains information on total SNPs analyzed, percent variance explained (PVE), time to complete GWAS analysis for the given phenotype, etc.
      -____.assoc.txt contains all SNPs analyzed (rows), SNP coordinate in genome, a p_value for each SNP, etc. 
  -ScriptSigSNPIntersectValidation.R:
    -for each phenotype GEMMA .assoc.txt output from the replication and 2022 paper, significant SNPs are subsetted (p < 0.05), the total number of significant SNPs is caluclated for each               phenotype in the replication and from the 2022 paper, and the total intersecting is computed using the SNP ID (rs)
  -ScriptCandidateGeneSNPOverlapValidation.R
    -for each phenotype GEMMA .assoc.txt output from the replication and 2022 paper, significant SNPs are subsetted (p <0.05), and the number of significant SNP hits is calculated for each 
     glucosinolate metabolite candidate gene. This is done by using the genomic coordinate of each significant SNPs and determining if they fall within the range of each candidate gene.
  -ScriptSNPALLPVals.R
    -Merge all indolic and aliphatic GWAS outputs into a dataset that contains SNP information (rs, ps, chr, af, and p_wald) where there is a column for each indolic or each aliphatic phenotypes 
     p_values for each SNP. It additionally creates a column for the best p_value (the minimum p_value) accross the row. This column is used to generate manhattan plots to analyze the overall GWAS 
     results for indolic and aliphatic phenotypes
  -ScriptManhattans.R
    -Writes a function to generate manhattan plots that highlight candidate gene regions, generates log10p_values, and generates significance thresholds that account for multiple hypothesis 
     testing with Bonferonni correction
     
