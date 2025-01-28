library(tidyverse)
# PrepForGemma

# Load BLUPs and Accessions
dat_200 = read.csv("2022_Inputs/phen_glucosinolatesBLUPs_2021-08-02.csv")
acc_200 = read.csv("2022_Inputs/TOU_accessions.csv")

acc_200_cleaned <- acc %>%
  separate(SAMPLE.ID.ecotype_id.pop_code, into = c("SAMPLE", "ID", "ecotype_id", "pop_code"), sep = ",") 


# merge with list of genotypes, ordered as in the .fam file
# (from the set of PLINK binary files used as genotypes in GWAS with GEMMA)
fam.file.IDs.200 = read.csv("2022_Inputs/genotypes/genotypes.fam", sep = " ", header = F)[,1]

dat_200 = merge(acc_200_cleaned, dat_200, by = "ID", all.x = T)

dat_200 = dplyr::left_join(data.frame(ecotype_id= as.character(fam.file.IDs.200)), dat_200, by = "ecotype_id")

# glucosinolates to retain
dat_200 = dat_200[,c("ecotype_id","gsl.7mSh.blup","gsl.8MTO.blup","gsl.3mSOp.blup","gsl.4mSOb.blup",
                     "gsl.5mSOp.blup","gsl.6mSOh.blup","gsl.7mSOh.blup","gsl.8mSOo.blup","gsl.Pren.blup",
                     "gsl.Buen.blup","gsl.Peen.blup","gsl.S2hBuen.blup","gsl.2hPeen.blup","gsl.IM.blup",
                     "gsl.1moIM.blup","gsl.1hIM.blup","gsl.4moIM.blup")]

colnames(dat_200)[1] = "ID"

# write table prepped for GEMMA (this file is the input for generating kinship matrix and for GEMMA lmm)
write.table(dat_200[,c(2:ncol(dat_200))], "dat_tot_ForGemmaGF_200.csv",
            sep = ",", row.names = F, col.names = F, quote = F, eol = "\n")