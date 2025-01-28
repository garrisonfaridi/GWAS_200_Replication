# GEMMA Kinship matrix
-gk -miss 0.1 -bfile genotypes -p dat_tot_ForGemmaGF_200.csv -o dat_tot_ForGemma_kmat_200GF

# GEMMA LMM GWAS
{1..9}; do gemma -lmm 1 -maf 0.03 -miss 0.1 -bfile genotypes -k dat_tot_ForGemma_kmat_200GF.cXX.txt -p dat_tot_ForGemmaGF_200.csv -n $i -o 200_gluc${i}_maf0.03_miss0.1GF; done
{10..17}; do gemma -lmm 1 -maf 0.03 -miss 0.1 -bfile genotypes -k dat_tot_ForGemma_kmat_200GF.cXX.txt -p dat_tot_ForGemmaGF_200.csv -n $i -o 200_gluc${i}_maf0.03_miss0.1GF; done