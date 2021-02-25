setwd("~/Desktop/Charr_Reanalysis/")
library(RcppCNPy)
library(qvalue)
library(tidyverse)
library(windowscanr)
library(data.table)
library(ggplot2)
library(ggman)
library(data.table)
library(dplyr)
library(vegan)

#Charr_LG_Chrom_Conversion
CHar_LG_Chrom_Conversion <- read.delim("~/Desktop/Charr_Final_Scripts_Data/CHar_LG_Chrom_Conversion.txt", stringsAsFactors=FALSE, header = F)
colnames(CHar_LG_Chrom_Conversion) <- c("LG", "CM_Chrom", "NC_Chrom")

CNV_GENOS<- fread("GDL_CNV_Matrix.txt", data.table = F, stringsAsFactors = F)
CNV_MAP <- fread("GDL_CNV_mapRDA_LGs", data.table = F, stringsAsFactors = F)
colnames(CNV_MAP)
Charrchip_map <- read.delim("Charrchip_Poly_GDL_comp.map", header=FALSE, stringsAsFactors=FALSE)
colnames(Charrchip_map) <- c("NC_Chrom", "SNP", "whatev", "BP")
CHar_LG_Chrom_Conversion <- read.delim("CHar_LG_Chrom_Conversion.txt", stringsAsFactors=FALSE, header = F)
colnames(CHar_LG_Chrom_Conversion) <- c("LG", "CM_Chrom", "NC_Chrom")

CNV_MAP <- inner_join(CNV_MAP, CHar_LG_Chrom_Conversion)

GDL_PHENO <- read.delim("GDL_PHENO3.txt", header = T, stringsAsFactors = F)
colnames(GDL_PHENO) <- c("ID", "jar_ID", "Morph", "Year")
GDL_PHENOGENO_CNV<- inner_join(GDL_PHENO, CNV_GENOS)
GDL_CNV <- GDL_PHENOGENO_CNV[5:length(GDL_PHENOGENO_CNV)]

GDL_CNV_RDA <- rda(GDL_CNV ~ GDL_PHENOGENO_CNV$Morph + Condition(GDL_PHENOGENO_CNV$Year))

GDL_CNV_RDA_YEAR <- varpart(GDL_CNV, ~ GDL_PHENOGENO_CNV$Morph, ~ GDL_PHENOGENO_CNV$Year)


anova.cca(GDL_CNV_RDA, parallel=getOption("mc.cores"), by = "terms")
showvarparts(parts = 2)
plot(GDL_CNV_RDA_YEAR)


#Significance testing
RsquareAdj(GDL_CNV_RDA)
anova.cca(GDL_CNV_RDA, parallel=getOption("mc.cores")) # default is permutation=999
screeplot(GDL_CNV_RDA)

bg <- c("dodgerblue", "red")
eco <- as.factor(GDL_PHENOGENO_CNV$Morph)
plot(GDL_CNV_RDA, type="n", scaling=3)
points(GDL_CNV_RDA, display="species", pch=20, cex=0.7, col="gray32", scaling=3)           # the SNPs
points(GDL_CNV_RDA, display="sites", pch=21, cex=1.3, col="gray32", scaling=3, bg=bg[eco]) # the fish
text(GDL_CNV_RDA, scaling=3, display="bp", col="#0868ac", cex=1)                           # the predictors
legend("topleft", legend=levels(eco), bty="n", col="gray32", pch=21, cex=0.5, pt.bg=bg)
#Get SNP scores

GDL_RDA_CNVSCORES <- data.frame(GDL_CNV_RDA$CCA$v[,1], stringsAsFactors = F)

RDA_CNVSCORES_MAPPED <- cbind(CNV_MAP, GDL_RDA_CNVSCORES)
colnames(RDA_CNVSCORES_MAPPED)[6] <- "RDA1"

max(RDA_CNVSCORES_MAPPED$RDA1)
RDA_CNVSCORES_MAPPED$RDA1_abs <-  abs(as.numeric(as.character(RDA_CNVSCORES_MAPPED$RDA1)))

#Plot RDA1
#Get outliers for each axis
RDA1_CNV_GDL_OL <- RDA_CNVSCORES_MAPPED[which(RDA_CNVSCORES_MAPPED$RDA1_abs > quantile(x = RDA_CNVSCORES_MAPPED$RDA1_abs, 0.99 )),]
RDA1_CNV_GDL_OL$SNP <- rownames(RDA1_CNV_GDL_OL)

library(ggman)
CNV_RDA <- ggman(RDA_CNVSCORES_MAPPED, chrom = "LG", pvalue = "RDA1_abs", snp = "SNP", bp="POS", pointSize = 1, title = "GDLCHIP_CNV", xlabel = "Chromosome", ymax = 0.3, logTransform = F, sigLine = 0.02514899, ylabel = "RDA1" ) + theme_classic()
ggmanHighlight(ggmanPlot = CNV_RDA, highlight = GDL_VST_OL$SNP)

#VST OL
GDL_pale_GDL_dark_VST <- read.delim("GDL_dark_rn_GDL_pale_rn.Vst.txt", header=FALSE, stringsAsFactors=FALSE)
colnames(GDL_pale_GDL_dark_VST) <- c("SNP", "LG_alt", "BP", "C", "VST")
CHar_LG_Chrom_Conversion_VST <- read.table("CHar_LG_Chrom_Conversion_VST.txt", quote="\"", comment.char="", stringsAsFactors=FALSE)
colnames(CHar_LG_Chrom_Conversion_VST) <- c("LG", "LG_alt", "CM_Chrom", "NC_Chrom")
GDL_pale_GDL_dark_VST <- inner_join(CHar_LG_Chrom_Conversion_VST, GDL_pale_GDL_dark_VST)

#get vst as well
RDA_plusVST <- unique(c(GDL_VST_OL_Pos$SNP, RDA1_CNV_GDL_OL$SNP))
GDL_RDAplusVST_OL_Pos <- RDA_CNVSCORES_MAPPED[RDA_CNVSCORES_MAPPED$SNP %in% RDA_plusVST,]

#make a bed file for gene ontology of RDA outliers AND VST outliers
write.table(cbind(GDL_RDAplusVST_OL_Pos$NC_Chrom, GDL_RDAplusVST_OL_Pos$POS, (GDL_RDAplusVST_OL_Pos$POS +1)), "Charr_CNV_VSTplusRDAOL.bed", sep = "\t", quote = F, row.names = F, col.names = F)

RDA_AND_VST <- GDL_VST_OL_Pos$SNP[GDL_VST_OL_Pos$SNP  %in% RDA1_CNV_GDL_OL$SNP]
GDL_RDAANDVST_OL_Pos <- RDA_CNVSCORES_MAPPED[RDA_CNVSCORES_MAPPED$SNP %in% RDA_AND_VST,]
GDL_RDAVST_OLD <- read.delim("RDA_and_VST_OL.txt", header = T)

#make a bed file for gene ontology of RDA outliers AND VST outliers
write.table(cbind(GDL_RDAANDVST_OL_Pos$NC_Chrom, GDL_RDAANDVST_OL_Pos$POS, (GDL_RDAANDVST_OL_Pos$POS +1)), "Charr_CNV_VSTANDRDAOL.bed", sep = "\t", quote = F, row.names = F, col.names = F)
write.table(GDL_RDAANDVST_OL_Pos, "RDA_and_VST_OL.txt", col.names = T, row.names = F, sep = "\t", quote = F)

VST_Man <- ggman(GDL_pale_GDL_dark_VST, snp = "SNP", chrom = "LG", bp = "BP", pvalue = "VST", logTransform = F, ymax = 1, sigLine = quantile(GDL_pale_GDL_dark_VST$VST, 0.99), pointSize = 2) + theme_classic()
ggmanHighlight(VST_Man, highlight = RDA_AND_VST)
