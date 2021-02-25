
charchip <- fread("~/Desktop/Charr_Reanalysis/TREEMIX/Charrchip_recall_filter01_mappedordered.ped",  data.table = F)

#Filter by population ####
GDL_PHENO <- read.delim("GDL_PHENO2.txt", header = T, stringsAsFactors = F)
DARK_MORPHS <- GDL_PHENO[GDL_PHENO$Phenotype %in% "Dark_morph",]
PALE_MORPHS <- GDL_PHENO[GDL_PHENO$Phenotype %in% "Pale_morph",]

table(charchip$V1[charchip$V1 %in% c("PBP", "MBB", "ENG", "IKI", "KAN")])
Charrchip_NFLD_Lab_pops <- charchip[charchip$V1 %in% c("PBP", "MBB", "ENG", "IKI", "KAN", "GDL"),]
Charrchip_NFLD_Lab_pops$V1[Charrchip_NFLD_Lab_pops$V2 %in% DARK_MORPHS$CHIP_ID] <- "DARK"
Charrchip_NFLD_Lab_pops$V1[Charrchip_NFLD_Lab_pops$V2 %in% PALE_MORPHS$CHIP_ID] <- "PALE"
fwrite(Charrchip_NFLD_Lab_pops , "~/Desktop/Charr_Reanalysis/TREEMIX/Charrchip_NFLD_Lab_pops.ped", sep="\t", col.names = FALSE, quote = FALSE, row.names = FALSE)


system("~/Desktop/Software/plink_mac_20200219/plink --ped TREEMIX/Charrchip_NFLD_Lab_pops.ped --map TREEMIX/Charrchip_recall_filter01_mappedordered.map --allow-extra-chr --maf .05 --recode --make-bed --out Charchip_NFLD_LAB")

#GDL Chip outliers
smnfOL <- read.delim("~/Desktop/Charr_Reanalysis/smnf_OL.txt", stringsAsFactors = F)
pcadapt_OL <- read.delim("~/Desktop/Charr_Reanalysis/PCAdapt_CHARRGDL_Poly_Outliers_K1_Q05.txt", stringsAsFactors = F)
Chip_OL <- inner_join(smnfOL, pcadapt_OL)

Charrchip_NFLD_Lab_pops$V1


#Calculate SNP PC scores
library(pcadapt)
CHARR_multi_pop <- read.pcadapt("Charchip_NFLD_LAB.ped", type = "ped")
PCs_multi_pop <- pcadapt(CHARR_multi_pop, K=3, method="mahalanobis",  min.maf = 0.01)
plot(PCs_multi_pop,option="screeplot") + geom_hline(yintercept = 0.05)


MAPPOS <- read.delim2("Charchip_NFLD_LAB.map", header = F)
colnames(MAPPOS) <- c("Chrom", "SNP", "CM", "BP")
PVALS <- PCs_multi_pop$pvalues
PCMAP <- as.data.frame(cbind(MAPPOS, PVALS))
library(qvalue)
qval <- qvalue(PCMAP$PVALS)
alpha <- 0.05
outliers <- which(qval$qvalues<alpha)
Outliers_CHARR <- PCMAP[outliers,]
ALL_OL <- unique(c(Chip_OL$SNP, Outliers_CHARR$SNP))
write.table(ALL_OL , "NFLD_LAB_GDL_OL.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
system("~/Desktop/Software/plink_mac_20200219/plink --file Charchip_NFLD_LAB --out Charchip_NFLD_LAB_NEUT --exclude NFLD_LAB_GDL_OL.txt --allow-extra-chr --recode ")
system("~/Desktop/Software/plink_mac_20200219/plink  --file Charchip_NFLD_LAB_NEUT --out char_unlink_NFLD_LAB --indep-pairwise 10 2 0.5 --allow-extra-chr")
system("~/Desktop/Software/plink_mac_20200219/plink  --file Charchip_NFLD_LAB_NEUT --extract char_unlink_NFLD_LAB.prune.in --out Charchip_NFLD_LAB_NEUT_LE --allow-extra-chr --make-bed --recode")

CLUSTS_NFLD_LAB <- data.frame(cbind(Charrchip_NFLD_Lab_pops$V1, Charrchip_NFLD_Lab_pops$V2, Charrchip_NFLD_Lab_pops$V1), stringsAsFactors = F)
write.table(CLUSTS_NFLD_LAB, "GDL_NFLD_LAB_cluster.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
system("~/Desktop/Software/plink_mac_20200219/plink --bfile Charchip_NFLD_LAB_NEUT_LE --freq --within GDL_NFLD_LAB_cluster.txt --out Charchip_NFLD_LAB_NEUT_LE  --allow-extra-chr")

