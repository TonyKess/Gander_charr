setwd("~/Desktop/Charr_Reanalysis/")
library(data.table)
library(tidyverse)
library(data.table)

#Identifying and removing CNVs, repetitive regions, then re-running popgen analyses
#Charr_LG_Chrom_Conversion
CHar_LG_Chrom_Conversion <- read.delim("~/Desktop/Charr_Final_Scripts_Data/CHar_LG_Chrom_Conversion.txt", stringsAsFactors=FALSE, header = F)
colnames(CHar_LG_Chrom_Conversion) <- c("LG", "CM_Chrom", "NC_Chrom")


#CNVNATOR
dark_cnvnator_edited <- read.delim("~/Desktop/Charr_Reanalysis/dark_cnvnator_edited.tsv", header=FALSE)
colnames(dark_cnvnator_edited) <- c("Type", "CM_Chrom", "BP", "Pos2", "Size", "Normalized_RD",
                                    "pval1", "pval2","pval3", "pval4", "q")

pale_cnvnator_edited <- read.delim("~/Desktop/Charr_Reanalysis/pale_cnvnator_edited.tsv", header=FALSE)
colnames(pale_cnvnator_edited) <- c("Type", "CM_Chrom", "BP", "Pos2", "Size", "Normalized_RD",
                                    "pval1", "pval2","pval3", "pval4", "q")


pale_cnvnator_bed<- data.frame(cbind(pale_cnvnator_edited$CM_Chrom, pale_cnvnator_edited$BP, pale_cnvnator_edited$Pos2))
write.table(pale_cnvnator_bed, "pale_cnvnator.bed", col.names = F, row.names = F, sep = "\t", quote = F)

dark_cnvnator_bed<- data.frame(cbind(dark_cnvnator_edited$CM_Chrom, dark_cnvnator_edited$BP, dark_cnvnator_edited$Pos2))
write.table(dark_cnvnator_bed, "dark_cnvnator.bed", col.names = F, row.names = F, sep = "\t", quote = F)

#filter for comparison with SI outliers
dark_cnvnator_filtered <- dark_cnvnator_edited[which(dark_cnvnator_edited$q > 0.5),]
pale_cnvnator_filtered <- pale_cnvnator_edited[which(pale_cnvnator_edited$q > 0.5),]


pale_cnvnator_filtered_bed<- data.frame(cbind(pale_cnvnator_filtered$CM_Chrom, pale_cnvnator_filtered$BP, pale_cnvnator_filtered$Pos2))
write.table(pale_cnvnator_filtered_bed, "pale_cnvnator_filter.bed", col.names = F, row.names = F, sep = "\t", quote = F)

dark_cnvnator_filtered_bed<- data.frame(cbind(dark_cnvnator_filtered$CM_Chrom, dark_cnvnator_filtered$BP, dark_cnvnator_filtered$Pos2))
write.table(dark_cnvnator_filtered_bed, "dark_cnvnator_filter.bed", col.names = F, row.names = F, sep = "\t", quote = F)

#VST Signal intensity OL
GDL_pale_GDL_dark_CNV_OL <- read.delim("~/Desktop/Charr_Reanalysis/RDA_and_VST_OL.txt", header=T, stringsAsFactors=FALSE)
colnames(GDL_pale_GDL_dark_CNV_OL)[3] <- "BP"
CHar_LG_Chrom_Conversion_VST <- read.table("CHar_LG_Chrom_Conversion_VST.txt", quote="\"", comment.char="", stringsAsFactors=FALSE)
colnames(CHar_LG_Chrom_Conversion_VST) <- c("LG", "LG_alt", "CM_Chrom", "NC_Chrom")

GDL_VST <- read.delim("~/Desktop/Charr_Reanalysis/GDL_dark_rn_GDL_pale_rn.Vst.txt", header=FALSE, stringsAsFactors=FALSE)
colnames(GDL_VST) <- c("SNP", "LG_alt", "BP", "C", "VST")
colnames(CHar_LG_Chrom_Conversion_VST) <- c("LG", "LG_alt", "CM_Chrom", "NC_Chrom")
GDL_VST <- inner_join(GDL_VST, CHar_LG_Chrom_Conversion_VST)
GDL_VST_OL <- inner_join(GDL_pale_GDL_dark_CNV_OL, GDL_VST)


VSTANDRDAOL <- fread("Charr_VSTANDRDAOL.bed", data.table = F)
colnames(VSTANDRDAOL) <- c("NC_Chrom", "Pos1", "Pos2")
VSTANDRDAOL <- inner_join(VSTANDRDAOL,CHar_LG_Chrom_Conversion)


#Removing ALL CNV and repetitive regions.
#Repeatmasker
Charr_repeat <- fread("Charr_repeatregions_simple.txt", stringsAsFactors = F, data.table = F)
Charr_repeat <- inner_join(CHar_LG_Chrom_Conversion, Charr_repeat)
Charr_repeat_bed <- data.frame(cbind(Charr_repeat$CM_Chrom, Charr_repeat$Pos1, Charr_repeat$Pos2))

colnames(Charr_repeat_bed) <- c("CM_Chrom", "Pos1", "Pos2")
colnames(dark_cnvnator_bed) <-  c("CM_Chrom", "Pos1", "Pos2")
colnames(pale_cnvnator_bed) <-  c("CM_Chrom", "Pos1", "Pos2")
VSTANDRDAOL <- data.frame(cbind(VSTANDRDAOL[5],VSTANDRDAOL[2:3]))
colnames(VSTANDRDAOL) <-  c("CM_Chrom", "Pos1", "Pos2")
All_potential_reps <- data.frame(rbind(Charr_repeat_bed, dark_cnvnator_bed, pale_cnvnator_bed, VSTANDRDAOL))

All_potential_reps$Pos1 <- as.numeric(All_potential_reps$Pos1)
All_potential_reps$Pos2 <- as.numeric(All_potential_reps$Pos2)
All_potential_reps <- All_potential_reps %>% arrange(CM_Chrom, Pos1)

write.table(All_potential_reps, "AP_reps.bed", sep = "\t", col.names = F, row.names = F, quote = F)
Charr_map <- fread("Charrchip_Poly_GDL_comp.map", data.table = F, stringsAsFactors = F)
colnames(Charr_map)[1] <- "NC_Chrom"
Charr_map_CM <- inner_join(Charr_map, CHar_LG_Chrom_Conversion)
write.table(cbind(Charr_map_CM$CM_Chrom, Charr_map_CM[2:4]), "Charrchip_Poly_GDL_comp_CM.map", 
            sep = "\t", col.names = F, row.names = F, quote = F)


###Sort it out
system("~/Desktop/Software/plink_mac_20200219/plink --ped Charrchip_Poly_GDL_comp.ped --map Charrchip_Poly_GDL_comp_CM.map --recode vcf --allow-extra-chr  --out Charrchip_Poly_GDL_comp")
system("bedtools merge -i AP_reps.bed > AP_reps_merged.bed")
system("bedtools intersect -a Charrchip_Poly_GDL_comp.vcf -b AP_rep_merge.bed -v > Charr_Poly_GDL_comp_norep.vcf")

AP_reps_merged <- read.delim("AP_reps_merged.bed", stringsAsFactors = F, header = F)
AP_reps_merged$size <- AP_reps_merged$V3 - AP_reps_merged$V2
colnames(AP_reps_merged)[1] <- "CM_Chrom"
AP_reps_merged <- inner_join(AP_reps_merged, CHar_LG_Chrom_Conversion)

ggplot() + geom_density(data = AP_reps_merged, aes(x = size)) + theme_classic()
max(AP_reps_merged$size)
#"awk '{OFS="\t"}{print $1,$2}' Charr_Poly_GDL_comp_norep.vcf > Charr_norep.tsv probably there's a nice way to do this in R terminal but just open a command line terminal 
Charrchip_norep <- read.delim("Charr_norep.tsv", header = F)
colnames(Charrchip_norep) <- c("CM_Chrom", "V4")
Charrchip_norep <- inner_join(Charrchip_norep, Charr_map_CM)
write.table(Charrchip_norep$V2, "Charrchip_norep_SNPS.txt", sep = "\t",
            col.names = F, row.names = F, quote = F)

system("~/Desktop/Software/plink_mac_20200219/plink --file Charrchip_Poly_GDL_comp --recode --maf 0.05 --extract Charrchip_norep_SNPS.txt --allow-extra-chr  --out Charrchip_Poly_GDL_comp_norep")
##Check how this performs with quick outlier test PCA etc.
library(pcadapt)
CHARR_GANDER <- read.pcadapt("Charrchip_Poly_GDL_comp_norep.ped", type = "ped")
PCs_GDL <- pcadapt(CHARR_GANDER, K=1, method="mahalanobis",  min.maf = 0.01)
PCs_GDL$singular.values ^ 2
plot(PCs_GDL,option="scores")
(PCs_GDL$singular.values)^2
GDL_IDs <-system("awk '{print $2}' Charrchip_Poly_GDL_comp_norep.ped ", intern = T)
GDL_PC_SCORES <- as.data.frame(cbind(GDL_IDs, PCs_GDL$scores))
GDL_PHENO <- read.delim("GDL_PHENO2.txt", header = T, stringsAsFactors = F)
GDL_PC_SCORES <- merge(GDL_PHENO, GDL_PC_SCORES, by.x = "CHIP_ID", by.y = "GDL_IDs")
GDL_PC_SCORES$V2 <- as.numeric(as.character(GDL_PC_SCORES$V2))
##for K=2
##GDL_PC_SCORES$V3 <- as.numeric(as.character(GDL_PC_SCORES$V3))
plot(PCs_GDL, option="scores") 
#Nicer Plot
colnames(GDL_PC_SCORES) = c("CHIP_ID","Jar_ID","Phenotype","Year","PC1")
##for K=2 
###ggplot(GDL_PC_SCORES, aes(x=PC1, y = PC2, col=Phenotype, shape=Year)) + geom_point()  + xlab("PC1") + ylab("PC2") + scale_color_manual(values = c("antiquewhite4","burlywood")) + theme_classic()
colnames(GDL_PC_SCORES) = c("CHIP_ID","Jar_ID","Phenotype","Year","PC1")
ggplot(GDL_PC_SCORES, aes(x=PC1, fill=Phenotype)) + geom_density()  + xlab("PC1")  + scale_fill_manual(values = c("antiquewhite4","burlywood")) + theme_classic()

MAPPOS <- read.delim2("Charrchip_Poly_GDL_comp_norep.map", header = F)
colnames(MAPPOS) <- c("Chrom", "SNP", "CM", "BP")
PVALS <- PCs_GDL$pvalues
PCMAP <- as.data.frame(cbind(MAPPOS, PVALS))
library(qvalue)
qval <- qvalue(PCs_GDL$pvalues)$qvalues
alpha <- 0.05
outliers <- which(qval<alpha)
Outliers_CHARR_norep <- PCMAP[outliers,]
write.table(Outliers_CHARR_norep, "PCAdapt_CHARRGDL_Poly_Outliers_K1_Q05_noreps.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

inner_join(cnvnator_filtered_distinct_genes_bed_NC, cnvnator_filtered_distinct_genes_NC_old)
cnvnator_filtered_distinct_genes_bed_NC$V2 <- as.integer(cnvnator_filtered_distinct_genes_bed_NC$V2)
cnvnator_filtered_distinct_genes_bed_NC$V3 <- as.integer(cnvnator_filtered_distinct_genes_bed_NC$V3)

colnames(cnvnator_filtered_distinct_genes_bed_NC) <- c("V1", "V2", "V3", "Symbol")

#just neutral LE SNPS
system("~/Desktop/Software/plink_mac_20200219/plink --file Charchip_GDL_Neut_Poly_LE --recode --maf 0.05 --extract Charrchip_norep_SNPS.txt --allow-extra-chr  --out Charrchip_Poly_GDL_Neut_Poly_norep")
library(pcadapt)
CHARR_GANDER <- read.pcadapt("Charrchip_Poly_GDL_Neut_Poly_norep.ped", type = "ped")
PCs_GDL <- pcadapt(CHARR_GANDER, K=1, method="mahalanobis",  min.maf = 0.01)
PCs_GDL$singular.values ^ 2
plot(PCs_GDL,option="scores")
(PCs_GDL$singular.values)^2
GDL_IDs <-system("awk '{print $2}' Charrchip_Poly_GDL_Neut_Poly_norep.ped ", intern = T)
GDL_PC_SCORES <- as.data.frame(cbind(GDL_IDs, PCs_GDL$scores))
GDL_PHENO <- read.delim("GDL_PHENO2.txt", header = T, stringsAsFactors = F)
GDL_PC_SCORES <- merge(GDL_PHENO, GDL_PC_SCORES, by.x = "CHIP_ID", by.y = "GDL_IDs")
GDL_PC_SCORES$V2 <- as.numeric(as.character(GDL_PC_SCORES$V2))
##for K=2
##GDL_PC_SCORES$V3 <- as.numeric(as.character(GDL_PC_SCORES$V3))
plot(PCs_GDL, option="scores") 
#Nicer Plot
colnames(GDL_PC_SCORES) = c("CHIP_ID","Jar_ID","Phenotype","Year","PC1")
##for K=2 
###ggplot(GDL_PC_SCORES, aes(x=PC1, y = PC2, col=Phenotype, shape=Year)) + geom_point()  + xlab("PC1") + ylab("PC2") + scale_color_manual(values = c("antiquewhite4","burlywood")) + theme_classic()
colnames(GDL_PC_SCORES) = c("CHIP_ID","Jar_ID","Phenotype","Year","PC1")
ggplot(GDL_PC_SCORES, aes(x=PC1, fill=Phenotype)) + geom_density()  + xlab("PC1")  + scale_fill_manual(values = c("antiquewhite4","burlywood")) + theme_classic()



library(ggman)
ggman(PCMAP, chrom = "Chrom", pvalue = "PVALS", snp = "SNP", bp="BP", pointSize = 1, title = "CHARR", xlabel = "Chromosome", sigLine = -log10(max(Outliers_CHARR$PVALS)),ymax = 300) + theme_classic()


#PCangsd selscan 
library(RcppCNPy)
library(qvalue)
S <- npyLoad("PC_charr_mq30.selection.npy") # Reads results from selection scan
sites <- fread("PC_charr_mq30.sites", data.table = F, stringsAsFactors = F, header = F)
pcangsd_sel <- cbind(S, sites)
colnames(pcangsd_sel) <- c("PC_S","CM_Chrom", "BP")
write.table((cbind(sites$V1, sites$V2, sites$V2 + 1)), "pcangsd_charr_sites.bed", sep = "\t", quote = F, row.names = F, col.names = F)

system("bedtools intersect -a pcangsd_charr_sites.bed -b AP_rep_merge.bed -v > pcangsd_sites_norep_jr.tsv")


#PCA and admixture with filtered 
setwd("~/Desktop/Tony/Charr/")
library(RcppCNPy)
library(ggplot2)
library(dplyr)
library(tidyverse)
C <- read.delim("PC_charr_mq30_norepmaf05.cov", sep = "", header = F)
PCA <- eigen(C)
PCA_scores <- data.frame(PCA$vectors[,1:2], stringsAsFactors = F)
plot(PCA_scores$X1, PCA_scores$X2, )
PCA_scores$MORPH <- c(rep("Pelagic", 10), rep("Deep-water", 10))
ggplot(data = PCA_scores) + geom_density(aes(x = X1, fill = MORPH)) + theme_classic() + scale_fill_manual(values = c("burlywood", "antiquewhite4"))

eigenvalues <- PCA$values
eigenvalues[1]/sum(eigenvalues)

Admix <- data.frame(npyLoad("PC_charr_mq30_norepmaf05.admix.Q.npy"), stringsAsFactors = F)
Admix$MORPH <- PCA_scores$MORPH
Admix$ID <- as.numeric(rep(1:20))
colnames(Admix)[1:2] <- c("K1", "K2")

plot_data <-  Admix %>% 
  gather('K', 'ADMIX', K1:K2) 

ggplot(plot_data, aes(ID, ADMIX, fill = K)) +
  geom_col() + theme_classic()


#check outlier overlap
CHar_LG_Chrom_Conversion <- read.delim("~/Desktop/Charr_Final_Scripts_Data/CHar_LG_Chrom_Conversion.txt", stringsAsFactors=FALSE, header = F)
colnames(CHar_LG_Chrom_Conversion) <- c("LG", "CM_Chrom", "NC_Chrom")

S <- npyLoad("PC_charr_mq30_norepmaf05.selection.npy") # Reads results from selection scan
sites <- fread("PC_charr_mq30_norepmaf05.sites", data.table = F, stringsAsFactors = F, header = F)
pcangsd_sel <- cbind(S, sites)
colnames(pcangsd_sel) <- c("PC_S","CM_Chrom", "BP")
pcangsd_sel <- inner_join(pcangsd_sel, CHar_LG_Chrom_Conversion)

pcangsd_sel_maxes <- pcangsd_sel %>% 
  group_by(LG) %>%
  filter(BP == max(BP))
pcangsd_sel_maxes <- data.frame(pcangsd_sel_maxes, stringsAsFactors = F)
pcangsd_sel_maxes$BP <- as.integer(pcangsd_sel_maxes$BP)
pcangsd_sel_maxes$BP <- pcangsd_sel_maxes$BP + 2000
pcangsd_sel_dummy <- rbind(pcangsd_sel_maxes, pcangsd_sel)
pcangsd_sel_1kwin <- winScan(x = pcangsd_sel_dummy, groups  = "LG", position = "BP", win_size = 1000, win_step = 1000,  values = "PC_S", funs = "mean")
pcangsd_sel_1kwin <- pcangsd_sel_1kwin[which(pcangsd_sel_1kwin$PC_S_n > 0),]
pcangsd_sel_1kwin$BP <- pcangsd_sel_1kwin$win_mid
#write out so that windows don't need to be recalculated (should have done this much earlier...)
write.table(pcangsd_sel_1kwin, "pcangsd_sel_1kwin.txt", col.names = T, sep = "\t", row.names = F, quote = F)
pcangsd_sel_1kwin <- fread("pcangsd_sel_1kwin.txt", data.table = F, stringsAsFactors = F)
pcangsd_sel_1k_99OL<- pcangsd_sel_1kwin[which(pcangsd_sel_1kwin$PC_S_mean > quantile(pcangsd_sel_1kwin$PC_S_mean, 0.99, na.rm=T)),]
pcangsd_sel_1k_99OL$BP <- pcangsd_sel_1k_99OL$win_mid

pcangsd_1kOL_LGBP<- data.frame(cbind("LG" = pcangsd_sel_1k_99OL$LG, "BP" = pcangsd_sel_1k_99OL$BP), stringsAsFactors = F)
pcangsd_1kOL_LGBP$BP<- as.numeric(pcangsd_1kOL_LGBP$BP)

ALL_WG_OL_combo <- read.delim("Charr_WGOL_FST_and_PCA_OL99.bed", header = F)
colnames(ALL_WG_OL_combo) <- c("CM_Chrom", "BP", "BPplsuone")
combined <- inner_join(ALL_WG_OL_combo, pcangsd_1kOL_LGBP)

#Also redo PCA output real quick for a nice density plot
C <- read.delim("PC_charr_mq30.cov", sep = "", header = F)
PCA <- eigen(C)
PCA_scores <- data.frame(PCA$vectors[,1:2], stringsAsFactors = F)
plot(PCA_scores$X1, PCA_scores$X2, )
PCA_scores$MORPH <- c(rep("Pelagic", 10), rep("Deep-water", 10))
ggplot(data = PCA_scores) + geom_density(aes(x = X1, fill = MORPH)) + theme_classic() + scale_fill_manual(values = c("burlywood", "antiquewhite4"))

eigenvalues <- PCA$values
eigenvalues[1]/sum(eigenvalues)

###Identifying distinct CNVs 
#Import annotation
charr_genome_annotation <- read.csv("~/Desktop/Charr_Reanalysis/charr_genome_annotation.csv")
colnames(charr_genome_annotation)[11] <- "LG"
charr_genome_annotation <- inner_join(charr_genome_annotation, CHar_LG_Chrom_Conversion)
write.table(cbind(charr_genome_annotation$CM_Chrom, charr_genome_annotation$start_position_on_the_genomic_accession, charr_genome_annotation$end_position_on_the_genomic_accession, charr_genome_annotation$Symbol),
            "charr_genes.bed", sep = "\t", quote = F, col.names = F, row.names = F)

system("bedtools intersect -a charr_genes.bed -b dark_cnvnator_filter.bed -wb > dark_cnvnator_charrgenes.bed")
system("bedtools intersect -a charr_genes.bed -b pale_cnvnator_filter.bed -wb  > pale_cnvnator_charrgenes.bed")
system("bedtools intersect -a charr_genes.bed -b Charr_VSTANDRDAOL.bed -wb > RDAVST_charrgenes.bed")
#import overlap
dark_genes_cnvnator <- read.delim("~/Desktop/Charr_Reanalysis/dark_cnvnator_charrgenes.bed", header = F)
colnames(dark_genes_cnvnator)<- c("CM_Chrom", "OL_start", "OL_stop", "Symbol", "Chrom2", "BP","Pos2")
dark_genes_cnvnator <- inner_join(dark_genes_cnvnator, dark_cnvnator_edited, by = c("CM_Chrom", "BP"))
dark_genes_cnvnator$morph <- "dark"
dark_genes_cnvnator$method <- "cnvnator"
dark_genes_cnvnator_unique <- unique(dark_genes_cnvnator$Symbol)

#import overlap
pale_genes_cnvnator <- read.delim("~/Desktop/Charr_Reanalysis/pale_cnvnator_charrgenes.bed", header = F)
colnames(pale_genes_cnvnator)<- c("CM_Chrom", "OL_start", "OL_stop", "Symbol", "Chrom2", "BP","Pos2")
pale_genes_cnvnator <- inner_join(pale_genes_cnvnator, pale_cnvnator_edited, by = c("CM_Chrom", "BP"))
pale_genes_cnvnator$morph <- "pale"
pale_genes_cnvnator$method <- "cnvnator"
pale_genes_cnvnator_unique <- unique(pale_genes_cnvnator$Symbol)

dark_CNVnator_distinct <- anti_join(dark_genes_cnvnator, pale_genes_cnvnator, by = c("Symbol"))
pale_CNVnator_distinct <- anti_join( pale_genes_cnvnator, dark_genes_cnvnator, by = c("Symbol"))

pale_CNVnator_distinct_annot <-  inner_join(pale_CNVnator_distinct, charr_genome_annotation, by = "Symbol")
dark_CNVnator_distinct_annot <-  inner_join(dark_CNVnator_distinct, charr_genome_annotation, by = "Symbol")

dark_CNVnator_distinct_annot_short <- dark_CNVnator_distinct_annot %>% select(c("CM_Chrom.x", "Symbol", "description", "morph", "BP" ,"Pos2.x","Type", "Size","q", "pval1", "pval2"))
pale_CNVnator_distinct_annot_short <- pale_CNVnator_distinct_annot %>% select(c("CM_Chrom.x",  "Symbol", "description","morph", "BP", "Pos2.x","Type", "Size","q", "pval1", "pval2"))

distinct_annot_short <- rbind(pale_CNVnator_distinct_annot_short, dark_CNVnator_distinct_annot_short)
write.table(distinct_annot_short, "distinct_CNVnator_annotated.txt", col.names = T, row.names = F, quote = F, sep = "\t")

#Overlap AND VST info
write.table(cbind(GDL_VST_OL$CM_Chrom, GDL_VST_OL$BP, GDL_VST_OL$BP+1, GDL_VST_OL$VST), "GDL_VSTOL.bed", 
            sep = "\t", col.names = F, row.names = F, quote = F)
VST_OL_genes <- read.table("VSTOL_gene_overlap")
colnames(VST_OL_genes) <- c("CM_Chrom", "BP", "BP2", "Symbol")
VST_OL_genes <- inner_join(VST_OL_genes, GDL_VST_OL)
VST_OL_genes_annot <- inner_join(VST_OL_genes, charr_genome_annotation, by = "Symbol")

dark_VST_WG_genes <- inner_join(dark_CNVnator_distinct_annot_short, VST_OL_genes, by = "Symbol")
pale_VST_WG_genes <- inner_join(pale_CNVnator_distinct_annot_short, VST_OL_genes, by = "Symbol")

anti_join(dark_VST_WG_genes, pale_VST_WG_genes, by = "Symbol")
anti_join(pale_VST_WG_genes, dark_VST_WG_genes, by = "Symbol")

inner_join(dark_VST_WG_genes, pale_VST_WG_genes, by = "Symbol")

VST_WG_genes<- rbind(pale_VST_WG_genes, dark_VST_WG_genes)
unique(VST_WG_genes$Symbol)

#Outlier region, fst, PCA comparison
paledark_fst <- read.delim("~/Desktop/Charr_Reanalysis/Old Charr files/Charr_Linux(1)/paledark.fst", header=FALSE)
colnames(paledark_fst) <- c("CM_Chrom", "BP", "SNPs", "Coverage_proportion", "mean_coverage", "pops_compared", "FST")
paledark_fst <- inner_join(CHar_LG_Chrom_Conversion, paledark_fst)


#Import chip map
Charrchip_map <- read.delim("~/Desktop/Charr_Reanalysis/Charrchip_Poly_GDL_comp.map", header=FALSE, stringsAsFactors=FALSE)
colnames(Charrchip_map) <- c("NC_Chrom", "SNP", "whatev", "BP")
CHIP_FST <- read.delim("~/Desktop/Charr_Reanalysis/FST_GDLCharr_diveRsity.txt", stringsAsFactors=FALSE)
CHIP_FST <- data.frame(cbind(CHIP_FST$loci, CHIP_FST$Fst_WC),stringsAsFactors = F)
colnames(CHIP_FST) <- c("SNP", "FST")
CHIP_FST <- inner_join(CHIP_FST, Charrchip_map)
CHIP_FST <- inner_join(CHIP_FST, CHar_LG_Chrom_Conversion)
CHIP_FST$FST <-  as.numeric(CHIP_FST$FST)

CHIP_VST_WIN <- winScan(x = GDL_VST_OL_dummy, groups  = "LG", position = "BP", win_size = 1000, win_step = 1000,  values = "VST", funs = "mean")
CHIP_VST_WIN <- CHIP_VST_WIN[CHIP_VST_WIN$VST_n >0,]
CHIP_VST_WIN$BP <- CHIP_VST_WIN$win_mid
write.table(CHIP_VST_WIN, "CHIP_VST_WIN.txt", sep = "\t", quote = F, col.names = T, row.names = F)
CHIP_VST_WIN <- read.delim("CHIP_VST_WIN.txt", header = T)
VSTOL_FST <- inner_join(paledark_fst, CHIP_VST_WIN)
NOTVSTOL_FST <- anti_join(paledark_fst, CHIP_VST_WIN)
mean(NOTVSTOL_FST$FST)
mean(VSTOL_FST$FST)
mean(paledark_fst$FST)
wilcox.test(NOTVSTOL_FST$FST, VSTOL_FST$FST)


#FST and gene overlap of CNVs
pale_dark_fst_bed<- data.frame(cbind(paledark_fst$CM_Chrom, paledark_fst$BP , paledark_fst$BP + 1, paledark_fst$FST))
write.table(pale_dark_fst_bed, "paledark_fst.bed", col.names = F, row.names = F, sep = "\t", quote = F)

dark_cnvnator_filtered <- dark_cnvnator_edited[which(dark_cnvnator_edited$q > 0.5),]
pale_cnvnator_filtered <- pale_cnvnator_edited[which(pale_cnvnator_edited$q > 0.5),]
cnvnator_filtered <- rbind(dark_cnvnator_filtered, pale_cnvnator_filtered)
cnvnator_filtered <- distinct_at(cnvnator_filtered, vars(Type, CM_Chrom, BP, Pos2))
cnvnator_filtered <- inner_join(cnvnator_filtered, CHar_LG_Chrom_Conversion)
cnvnator_filtered_bed <- data.frame(cbind(cnvnator_filtered$CM_Chrom, cnvnator_filtered$BP, cnvnator_filtered$Pos2))

cnvnator_filtered_distinct <- rbind(dark_CNVnator_distinct, pale_CNVnator_distinct)
cnvnator_filtered_distinct <- inner_join(cnvnator_filtered_distinct, CHar_LG_Chrom_Conversion)
charr_genes_bed <- data.frame(cbind(charr_genome_annotation$CM_Chrom, charr_genome_annotation$start_position_on_the_genomic_accession, charr_genome_annotation$end_position_on_the_genomic_accession, charr_genome_annotation$Symbol))
cnvnator_filtered_distinct_genes_bed <- charr_genes_bed[charr_genes_bed$X4 %in% cnvnator_filtered_distinct$Symbol,]
colnames(cnvnator_filtered_distinct_genes_bed)[4] <- "Symbol"
write.table(cnvnator_filtered_distinct_genes_bed, "cnvnator_filtered_distinct_genes.bed", col.names = F, row.names = F, sep = "\t", quote = F )

charr_genes_bed_NC <- data.frame(cbind(charr_genome_annotation$NC_Chrom, charr_genome_annotation$start_position_on_the_genomic_accession, charr_genome_annotation$end_position_on_the_genomic_accession, charr_genome_annotation$Symbol))
cnvnator_filtered_distinct_genes_bed_NC <- charr_genes_bed_NC[charr_genes_bed_NC$X4 %in% cnvnator_filtered_distinct$Symbol,]
colnames(cnvnator_filtered_distinct_genes_bed_NC)[4] <- "Symbol"
write.table(cnvnator_filtered_distinct_genes_bed_NC, "cnvnator_filtered_distinct_genes_NC.bed", col.names = F, row.names = F, sep = "\t", quote = F )


write.table(pale_dark_fst_bed, "paledark_fst.bed", col.names = F, row.names = F, sep = "\t", quote = F)
system("bedtools intersect -a paledark_fst.bed -b cnvnator_filtered_distinct_genes.bed -wb > paledark_fst_cnvnator_filtered_distinctgenes")
Charr_WGOL99_bed <- read.table("Charr_WGOL_FST_and_PCA_OL99.bed")
colnames(Charr_WGOL99_bed)[1] <- "NC_Chrom"
Charr_WGOL99_bed<- inner_join(Charr_WGOL99_bed,CHar_LG_Chrom_Conversion )
write.table(cbind(Charr_WGOL99_bed$CM_Chrom, Charr_WGOL99_bed$V2, Charr_WGOL99_bed$V3), "Charr_WGOL_andPC_99_CM.bed", col.names = F, row.names = F, sep = "\t", quote = F)
system("bedtools intersect -a Charr_WGOL_andPC_99_CM.bed -b cnvnator_filtered_distinct_genes.bed -wb  > WGOL99_cnvnator_filtered")

paledark_fst_cnvnator_filtered <- fread("paledark_fst_cnvnator_filtered_distinctgenes", data.table = F, stringsAsFactors = F)
paledark_WGOL99_cnvnator_filtered <- fread("WGOL99_cnvnator_filtered", data.table = F, stringsAsFactors = F)
colnames(paledark_WGOL99_cnvnator_filtered) <- c("CM_Chrom", "BP", "BPplusone", "CM_chrom_junior", "BP2", "BPsbackbaby", "Symbol")
paledark_WGOL99_cnvnator_filtered<- inner_join(paledark_fst, paledark_WGOL99_cnvnator_filtered)
paledark_WGOL99_cnvnator_filtered_annot <- inner_join(paledark_WGOL99_cnvnator_filtered, charr_genome_annotation, by = "Symbol")
colnames(paledark_fst_cnvnator_filtered) <-  c("CM_Chrom", "BP", "BPplusone", "FST","CM_chrom_junior", "BP2", "BPsbackbaby", "Symbol")
paledark_fst_cnvnator_filtered_annot <- inner_join(paledark_fst_cnvnator_filtered, charr_genome_annotation, by = "Symbol" )

colnames(ALL_WG_OL_combo)[1] <- "NC_Chrom"
ALL_WG_OL_combo <- inner_join(ALL_WG_OL_combo, CHar_LG_Chrom_Conversion)
ALL_WG_OL_combo <- inner_join(ALL_WG_OL_combo, paledark_fst, by = c("CM_Chrom", "BP"))

ggplot() + geom_density(data =paledark_fst, aes(x = FST, colour = "Genome average")) +
  geom_density(data =paledark_fst_cnvnator_filtered, aes(x = FST, colour = "CNVnator CNVs")) + 
  geom_density(data =Chip_OL, aes(x = FST, colour = "Chip FST OL")) +
  geom_density(data =ALL_WG_OL_combo , aes(x = FST, colour = "CWG OL")) +
  geom_density(data =VSTOL_FST, aes(x = FST, colour = "VST outliers")) + theme_classic() + geom_vline(xintercept = 0.11)
