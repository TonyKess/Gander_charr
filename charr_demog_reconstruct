setwd("~/Desktop/Charr_Reanalysis/")
library(data.table)
library(tidyverse)
library(data.table)

#plot stairways
pale_stairway <- fread("pale morph 100 to 200K years.final.summary", stringsAsFactors = F, data.table = F)
colnames(pale_stairway) <- c("mutation_per_site",     "n_estimation",          "theta_per_site_median", "theta_per_site_2.5",   "theta_per_site_97.5","year", "Ne_median" ,"Ne_2.5","Ne_97.5", "Ne_12.5",  "Ne_87.5" )

dark_stairway <- fread("dark morph 100 to 200K years.final.summary", stringsAsFactors = F, data.table = F)
colnames(dark_stairway) <- c("mutation_per_site",     "n_estimation",          "theta_per_site_median", "theta_per_site_2.5",   "theta_per_site_97.5","year", "Ne_median" ,"Ne_2.5","Ne_97.5", "Ne_12.5",  "Ne_87.5" )

min(pale_stairway$Ne_median)
min(dark_stairway$Ne_median)

#log scaled
ggplot() + geom_line(data = pale_stairway, aes(y = log10(Ne_median), x = log10(year), colour = "pale" )) + 
  geom_line(data = pale_stairway, aes(y = log10(Ne_2.5), x = log10(year), colour = "pale_bound" )) +
  geom_line(data = pale_stairway, aes(y = log10(Ne_97.5), x = log10(year), colour = "pale_bound" )) + 
  geom_line(data = dark_stairway, aes(y = log10(Ne_median), x = log10(year), colour = "dark" )) + 
  geom_line(data = dark_stairway, aes(y = log10(Ne_2.5), x = log10(year), colour = "dark_bound" )) +
  geom_line(data = dark_stairway, aes(y = log10(Ne_97.5), x = log10(year), colour = "dark_bound" )) + 
  geom_vline(xintercept = log10(5000), colour = "blue") + geom_vline(xintercept = log10(7000), colour = "blue") + 
  geom_vline(xintercept = log10(19000), colour = "green") + geom_vline(xintercept = log10(115000), colour = "green") +
  theme_classic() 

##SNEP
Charr_POP_ID_GDL <- read.delim("~/Desktop/Charr_Reanalysis/Charr_POP_ID_GDL.txt", header=FALSE)
Dark_morphs <- Charr_POP_ID_GDL[Charr_POP_ID_GDL$V2 %in% GDL_PHENO$CHIP_ID[GDL_PHENO$Phenotype %in% "Dark_morph"],]
write.table(Dark_morphs, "Dark_morphs.txt", sep = "\t", col.names = F, row.names = F, quote = F)
system("~/Desktop/Software/plink_mac_20200219/plink --file Charrchip_Poly_GDL_comp --recode --geno 0 --exclude Outlier_combo.txt --allow-extra-chr --keep Dark_morphs.txt --out Charrchip_Poly_GDL_comp_Dark")

Pale_morphs <- Charr_POP_ID_GDL[Charr_POP_ID_GDL$V2 %in% GDL_PHENO$CHIP_ID[GDL_PHENO$Phenotype %in% "Pale_morph"],]
write.table(Pale_morphs, "Pale_morphs.txt", sep = "\t", col.names = F, row.names = F, quote = F)
system("~/Desktop/Software/plink_mac_20200219/plink --file Charrchip_Poly_GDL_comp --recode --geno 0 --exclude Outlier_combo.txt --allow-extra-chr --keep Pale_morphs.txt --out Charrchip_Poly_GDL_comp_Pale")


Charr_SNep <- read.delim("Charr_SNEP_brief.txt", header = T, stringsAsFactors = F)
ggplot() + geom_point(data = Charr_SNep, aes(x = GenAgo * 4, y = Ne/max(Ne), colour = morph)) +
  geom_line(data = Charr_SNep, aes(x = GenAgo * 4, y = Ne/max(Ne), colour = morph)) + theme_classic()

