#PCADAPT #####
setwd("~/Desktop/Charr_Chip_Recalled/")
library(pcadapt)
library(ggplot2)
library(ggman)
library(data.table)
library(lfmm)
library(LEA)
library(qvalue)

#Calculate SNP PC scores
CHARR_GANDER <- read.pcadapt("Charrchip_Poly_GDL_comp.ped", type = "ped")
PCs_GDL <- pcadapt(CHARR_GANDER, K=1, method="mahalanobis",  min.maf = 0.01)
plot(PCs_GDL,option="screeplot")
GDL_IDs <-system("awk '{print $2}' Charrchip_Poly_GDL_comp.ped ", intern = T)
GDL_PC_SCORES <- as.data.frame(cbind(GDL_IDs, PCs_GDL$scores))
GDL_PHENO <- read.delim("GDL_PHENO2.txt", header = T, stringsAsFactors = F)
GDL_PC_SCORES <- merge(GDL_PHENO, GDL_PC_SCORES, by.x = "CHIP_ID", by.y = "GDL_IDs")
GDL_PC_SCORES$V2 <- as.numeric(as.character(GDL_PC_SCORES$V2))
##for K=2
##GDL_PC_SCORES$V3 <- as.numeric(as.character(GDL_PC_SCORES$V3))
plot(PCs_GDL, option="scores") 
#Nicer Plot
colnames(GDL_PC_SCORES) = c("CHIP_ID","Jar_ID","Phenotype","Year","PC1","PC2")
##for K=2 
###ggplot(GDL_PC_SCORES, aes(x=PC1, y = PC2, col=Phenotype, shape=Year)) + geom_point()  + xlab("PC1") + ylab("PC2") + scale_color_manual(values = c("antiquewhite4","burlywood")) + theme_classic()
colnames(GDL_PC_SCORES) = c("CHIP_ID","Jar_ID","Phenotype","Year","PC1")
ggplot(GDL_PC_SCORES, aes(x=PC1, fill=Phenotype)) + geom_density()  + xlab("PC1")  + scale_fill_manual(values = c("antiquewhite4","burlywood")) + theme_classic()


MAPPOS <- read.delim2("Charrchip_Poly_GDL_comp.map", header = F)
colnames(MAPPOS) <- c("Chrom", "SNP", "CM", "BP")
PVALS <- PCs_GDL$pvalues
PCMAP <- as.data.frame(cbind(MAPPOS, PVALS))
library(qvalue)
qval <- qvalue(PCs_GDL$pvalues)$qvalues
alpha <- 0.05
outliers <- which(qval<alpha)
Outliers_CHARR <- PCMAP[outliers,]
write.table(Outliers_CHARR, "PCAdapt_CHARRGDL_Poly_Outliers_K1_Q05.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
class(Outliers_CHARR$PVALS)
max(Outliers_CHARR$PVALS)

library(ggman)
ggman(PCMAP, chrom = "Chrom", pvalue = "PVALS", snp = "SNP", bp="BP", pointSize = 1, title = "CHARR", xlabel = "Chromosome", sigLine = -log10(max(Outliers_CHARR$PVALS)),ymax = 300) + theme_classic()


##RDA #####
library(vegan)
library(rsed)

#genotypes
GDL_RAW <- read.delim("Charrchip_GDL_Poly_PhenoorderA.raw", sep = "", stringsAsFactors = F)
#just IDs
GDL_POPIND <- GDL_RAW[1:2]
#morphometric values, morph info
GDL_MORPHO <- read.delim("Charr_Sequenced_Photographed_RW_ID_Morph.txt", stringsAsFactors = F)

#subset SNPs
GDL_SNPS <- GDL_RAW[7:3650]
#impute SNPs
GDL_SNPS <- apply(GDL_SNPS, 2, function(x) replace(x, is.na(x), as.numeric(names(which.max(table(x))))))
#RDA by significant morphometric variables for ecotype differentiation
GDL_MORPHO_RDA <- rda(GDL_SNPS ~ RW2 + RW3 + RW5 + RW6, data = GDL_MORPHO)
RsquareAdj(GDL_MORPHO_RDA)
summary(eigenvals(GDL_MORPHO_RDA, model = "constrained"))
screeplot(GDL_MORPHO_RDA)
anova.cca(GDL_MORPHO_RDA, parallel=getOption("mc.cores")) # default is permutation=999

# Plot axes 1 & 2
bg <- c("antiquewhite4","burlywood")
eco <- as.factor(GDL_MORPHO$MORPH)
plot(GDL_MORPHO_RDA, type="n", scaling=3)
points(GDL_MORPHO_RDA, display="species", pch=20, cex=0.7, col="gray32", scaling=3)           # the SNPs
points(GDL_MORPHO_RDA, display="sites", pch=21, cex=1.3, col="gray32", scaling=3, bg=bg[eco]) # the fish
text(GDL_MORPHO_RDA, scaling=3, display="bp", col="#0868ac", cex=1)                           # the predictors
legend("bottomright", legend=levels(eco), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)

#Get SNP scores
GDL_morhpo.rda <- scores(GDL_MORPHO_RDA, choices=c(1:3), display="species")  # Species scores for the first three constrained 
GDL_RDA_SNPSCORES <- as.data.frame(cbind(SNP = rownames(GDL_morhpo.rda), GDL_morhpo.rda))
#outlier function
outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]               # locus names in these tails
}
GDL_RDA_OL <- outliers(GDL_morhpo.rda[,1],3)
library(rsed)

GDL_RDA_OL
names(GDL_RDA_OL)
GDL_RDA_OL <- cbind.data.frame(rep(1,times=length(GDL_RDA_OL)), names(GDL_RDA_OL), unname(GDL_RDA_OL))
colnames(GDL_RDA_OL) <-  c("axis","snp","loading")
#match names to map file
GDL_RDA_OL$snp <- as.character(GDL_RDA_OL$snp)
GDL_RDA_OL$snp <- sed_substitute(GDL_RDA_OL$snp, "AX.", "AX-") 
GDL_RDA_OL$snp <- sed_substitute(GDL_RDA_OL$snp, "_.*", "")

GDL_RDA_SNPSCORES$SNP <- as.character(GDL_RDA_SNPSCORES$SNP)
GDL_RDA_SNPSCORES$SNP <- sed_substitute(GDL_RDA_SNPSCORES$SNP, "AX.", "AX-") 
GDL_RDA_SNPSCORES$SNP <- sed_substitute(GDL_RDA_SNPSCORES$SNP, "_.*", "")


#compare outlier PCADAPT and RDA overlap
Outliers_GDL_PCADAPT <- read.delim("PCAdapt_CHARRGDL_Outliers_K2_Q05.txt", stringsAsFactors = F)
GDL_RDA_OL$snp %in% Outliers_GDL_PCADAPT$SNP

#Outliers not in common
GDL_RDA_OL$snp[which((GDL_RDA_OL$snp %in% Outliers_GDL_PCADAPT$SNP)==F)]

#Get RDA OL map info

MAPPOS <- read.delim2("Charrchip_GDL_Poly_Phenoorder12.map", header = F)
colnames(MAPPOS) <- c("Chrom", "SNP", "CM", "BP")
GDL_RDA_OL_MAP <- MAPPOS[MAPPOS$SNP %in% GDL_RDA_OL$snp,]
GDL_RDA_OL_MAPPED <- merge(GDL_RDA_OL_MAP, GDL_RDA_OL, by.x = "SNP", by.y = "snp")
GDL_RDA_OL_MAPPED$SNP <- as.character(GDL_RDA_OL_MAPPED$SNP)
GDL_RDA_OL_MAPPED$BP <- as.numeric(GDL_RDA_OL_MAPPED$BP)
GDL_RDA_OL_MAPPED$Chrom <- as.character(GDL_RDA_OL_MAPPED$Chrom)
write.table(GDL_RDA_OL_MAPPED, "GDL_RDA_OL_Mapped.txt", col.names = T, row.names = F, quote = F, sep = "\t")
group <- rep("RDA Outlier", length(GDL_RDA_OL_MAPPED$Chrom))
GDL_RDA_OL_GROUP <- cbind(GDL_RDA_OL_MAPPED, group)

#Get RDA total map info
GDL_RDA_MAPPED <- merge(MAPPOS, GDL_RDA_SNPSCORES, by.x = "SNP", by.y = "SNP")
GDL_RDA_MAPPED$RDA1 <- as.numeric(as.character(GDL_RDA_MAPPED$RDA1))
GDL_RDA_MAPPED$RDA2 <- as.numeric(as.character(GDL_RDA_MAPPED$RDA2))
GDL_RDA_MAPPED$RDA3 <- as.numeric(as.character(GDL_RDA_MAPPED$RDA3))


GDL_RDA_Fig <- ggman(GDL_RDA_MAPPED, chrom = "Chrom", pvalue = "RDA1", snp = "SNP", bp="BP", logTransform = F, pointSize = 1, title = "CHARR", xlabel = "Chromosome", ylabel = "RDA Loading", ymax = 0.4, ymin = -0.4) + theme_classic()
OL_Highlight <- as.character(GDL_RDA_OL$snp)
ggmanHighlight(GDL_RDA_Fig, highlight = OL_Highlight )


###snmf Outlier Test####
#now ust complete dataset
ped2lfmm(input.file = "Charrchip_Poly_GDL_comp12.ped")

project = NULL
project = snmf("Charrchip_Poly_GDL_comp12.lfmm", K = 1:10, entropy = TRUE,repetitions = 10, project = "new")
# plot cross-entropy criterion for all runs in the snmf project
plot(project, col = "blue", pch = 19, cex = 1.2)
# select the best run for K = 4
best = which.min(cross.entropy(project, K = 2))
my.colors <- c("tomato", "lightblue")
barchart(project, K = 2, run = best,
         border = NA, space = 0,
         col = my.colors,
         xlab = "Individuals",
         ylab = "Ancestry proportions",
         main = "Ancestry matrix") -> bp
axis(1, at = 1:length(bp$order),
     labels = bp$order, las=1,
     cex.axis = .3)
# Population differentiation tests
p = snmf.pvalues(project,
                 entropy = TRUE,
                 ploidy = 2,
                 K = 2)
pvalues = p$pvalues
hist(pvalues, col = "orange")
plot(-log10(pvalues), pch = 19, col = "blue", cex = .7)
MAPPOS <- read.delim2("Charrchip_Poly_GDL_comp.map", header = F, stringsAsFactors = F)
colnames(MAPPOS) <- c("Chrom", "SNP", "CM", "BP")
smnf_Qvals <- qvalue(p$pvalues)$qvalue
smnf_frame <- data.frame(cbind(MAPPOS, smnf_Qvals, p$pvalues), stringsAsFactors = F)
smnf_OL <- smnf_frame[which(smnf_frame$smnf_Qvals<0.05),]
write.table(smnf_OL, "smnf_OL.txt", col.names = T, row.names = F, quote = F, sep = "\t")
library(ggman)
ggman(smnf_frame, chrom = "Chrom", pvalue = "p.pvalues", snp = "SNP", bp="BP", pointSize = 1, title = "CHARR", xlabel = "Chromosome", sigLine = -log10(max(smnf_OL$p.pvalues))) + theme_classic()

##Split populations by morph and by year for FST test and AMOVA in Arlequin ####
library(genepopedit)
GDL_PHENO <- read.delim("GDL_PHENO2.txt", stringsAsFactors = F)
GDL_PHENO <- GDL_PHENO[GDL_PHENO$CHIP_ID %in% GDL_gp$SampleID,]
GDL_gp <- genepop_flatten("Charrchip_Poly_GDL_comp.txt")
GDL_gp$SampleID <- as.character(GDL_gp$SampleID)
#For outlier tests etc.
#Dark morph
GDL_DARK <- GDL_PHENO[GDL_PHENO$Phenotype=="Dark_morph",]
GDL_DARK_gp <- GDL_gp[GDL_gp$SampleID %in% GDL_DARK$CHIP_ID,]

GDL_DARK_gp$SampleID <- paste0("Dark_", rep(1:length(GDL_DARK_gp$SampleID)))
GDL_DARK_gp$Population  <- rep("Dark", length(GDL_DARK_gp$SampleID))
#Pale morph
GDL_PALE <- GDL_PHENO[GDL_PHENO$Phenotype=="Pale_morph",]
GDL_PALE_gp <- GDL_gp[GDL_gp$SampleID %in% GDL_PALE$CHIP_ID,]

GDL_PALE_gp$SampleID <- paste0("Pale_", rep(1:length(GDL_PALE_gp$SampleID)))
GDL_PALE_gp$Population  <- rep("Pale", length(GDL_PALE_gp$SampleID))

#Combine
Paledark <- rbind(GDL_DARK_gp, GDL_PALE_gp)
genepop_unflatten(Paledark, path = paste0("Pale_dark.txt"))
#For AMOVA
#Dark morph 2017
GDL_DARK_2017 <- GDL_DARK[GDL_DARK$Year=="2017",]
GDL_DARK_gp_2017 <- GDL_gp[GDL_gp$SampleID %in% GDL_DARK_2017$CHIP_ID,]
GDL_DARK_gp_2017$SampleID <- paste0("Dark2017_", rep(1:length(GDL_DARK_gp_2017$SampleID)))
GDL_DARK_gp_2017$Population  <- rep("Dark2017", length(GDL_DARK_gp_2017$SampleID))

#Dark morph 2005/2006
GDL_DARK_2005_2006 <- GDL_DARK[GDL_DARK$Year=="2005_2006",]
GDL_DARK_gp_2005_2006 <- GDL_gp[GDL_gp$SampleID %in% GDL_DARK_2005_2006$CHIP_ID,]
GDL_DARK_gp_2005_2006$SampleID <- paste0("Dark20056_", rep(1:length(GDL_DARK_gp_2005_2006$SampleID)))
GDL_DARK_gp_2005_2006$Population  <- rep("Dark20056", length(GDL_DARK_gp_2005_2006$SampleID))


#Pale morph 2017
GDL_PALE_2017 <- GDL_PALE[GDL_PALE$Year=="2017",]
GDL_PALE_gp_2017 <- GDL_gp[GDL_gp$SampleID %in% GDL_PALE_2017$CHIP_ID,]
GDL_PALE_gp_2017$SampleID <- paste0("Pale2017_", rep(1:length(GDL_PALE_gp_2017$SampleID)))
GDL_PALE_gp_2017$Population  <- rep("Pale2017", length(GDL_PALE_gp_2017$SampleID))

#Pale morph 2005/2006
GDL_PALE_2005_2006 <- GDL_PALE[GDL_PALE$Year=="2005_2006",]
GDL_PALE_gp_2005_2006 <- GDL_gp[GDL_gp$SampleID %in% GDL_PALE_2005_2006$CHIP_ID,]
GDL_PALE_gp_2005_2006$SampleID <- paste0("Pale20056_", rep(1:length(GDL_PALE_gp_2005_2006$SampleID)))
GDL_PALE_gp_2005_2006$Population  <- rep("Pale20056", length(GDL_PALE_gp_2005_2006$SampleID))

Paledarkyears <- rbind(GDL_DARK_gp_2017, GDL_DARK_gp_2005_2006, GDL_PALE_gp_2017, GDL_PALE_gp_2005_2006)
genepop_unflatten(Paledarkyears, path = paste0("Pale_dark_years.txt"))

##Calculate FST and HET ####
library(diveRsity)
GDL_DIVPART <- fastDivPart(infile = "Pale_dark.txt", outfile = "FST_GDL_PaleDark.txt", gp = 3, pairwise = TRUE,
               fst = TRUE, bs_locus = FALSE, bs_pairwise = FALSE,
               boots = 0, plot = FALSE, para = T)

GDL_DIVPART$estimate
GDL_FST <- data.frame(GDL_DIVPART$estimate, stringsAsFactors = F)
GDL_FST$locus <- rownames(GDL_FST)
FST_MAP <- merge( MAPPOS, GDL_FST, by.y = "locus", by.x = "SNP")
ggman(FST_MAP, chrom = "Chrom", pvalue = "Fst_WC", snp = "SNP", bp="BP", logTransform = F, pointSize = 1, title = "CHARR", xlabel = "Chromosome", ylabel = "FST", ymax = 1, sigLine = 0.64273 ) + theme_classic()
FST_95Q <- GDL_FST$locus[which(GDL_FST$Fst_WC>quantile(GDL_FST$Fst_WC, 0.95, na.rm = T))]

#Get heterozygosity
system("../plink_mac/plink --file Charrchip_Poly_GDL_comp --out Charrchip_het --hardy --allow-extra-chr")
HET <- read.delim("Charrchip_het.hwe", stringsAsFactors = F, header = T, sep = "")
HET_MAP <- merge( MAPPOS, HET, by.y = "SNP", by.x = "SNP")
ggman(HET_MAP, chrom = "Chrom", pvalue = "O.HET.", snp = "SNP", bp="BP", logTransform = F, pointSize = 1, title = "CHARR", xlabel = "Chromosome", ylabel = "FST", ymax = 1) + theme_classic()


##FDIST2 Outliers ####
GDL_FDIST <- read.delim("Charr_GDL_fdist2_ObsOut.txt", header = T, stringsAsFactors = F)
#FDR for big-time FDIST outliers
qval_fdist <- qvalue(GDL_FDIST$FST.P.value, pi0 = 1)$qvalue
GDL_FDIST <- cbind(GDL_FDIST, qval_fdist)
alpha <- 0.05
outliers_fdist <- which(qval_fdist<alpha)
GDL_FDIST$Locus[outliers_fdist]
FDIST_Frame <- merge(MAPPOS, GDL_FDIST, by.x = "SNP", by.y = "Locus" )
ggman(FDIST_Frame, chrom = "Chrom", pvalue = "FST.P.value", snp = "SNP", bp="BP", pointSize = 1, title = "CHARR", xlabel = "Chromosome", sigLine = 5.1 ) + theme_classic()
fdist_uncorr <- GDL_FDIST$Locus[which(GDL_FDIST$FST.P.value<0.05)]
fdist_uncorr


###BAYESCAN OUTLIERS ####
GDL_BS <- read.delim("GDL_FST_fst.txt", header = T, stringsAsFactors = F)
BS_Frame <- merge(MAPPOS, GDL_BS, by.x = "SNP", by.y = "Loc" )

alpha <- 0.05
outliers_BS <- BS_Frame[which(BS_Frame$qval<alpha),]
write.table(outliers_BS, "GDL_Bayesscan_outliers.txt", col.names = T, row.names = F, sep = "\t", quote = F)
ggman(BS_Frame, chrom = "Chrom", pvalue = "qval", snp = "SNP", bp="BP", pointSize = 1, title = "CHARR", xlabel = "Chromosome", ymin = 0, ymax = 20, sigLine = -log10(0.05)) + theme_classic()

##All outliers ####
PCadapt_OL <- read.delim("PCAdapt_CHARRGDL_Poly_Outliers_K1_Q05.txt", stringsAsFactors = F, header = T)
smnf_OL <- read.delim("smnf_OL.txt", stringsAsFactors = F, header = T)
fdist_OL <- GDL_FDIST$Locus[outliers_fdist]  ##there's only 3 loci and they overlap completely with other methods, not bothering writing out a file. Also fdist is semi-stochastic in what it calls an outlier.
BS_OL <- read.delim("GDL_Bayesscan_outliers.txt", stringsAsFactors = F, header = T)

#overlap
Two_method_OL<- snmf_OL$SNP[snmf_OL$SNP %in% PCadapt_OL$SNP]
Three_method_OL <- Two_method_OL[Two_method_OL %in% BS_OL$SNP]
Four_method_OL <- Three_method_OL[Three_method_OL %in% fdist_OL]
Anymethod_OL <- unique(c(smnf_OL$SNP, PCadapt_OL$SNP, fdist_OL, BS_OL$SNP))
Anymethod_OL_Frame <- MAPPOS[MAPPOS$SNP %in% Anymethod_OL,]
Important_OL_Frame <- MAPPOS[MAPPOS$SNP %in% Four_method_OL,]
Two_method_OL_frame <- MAPPOS[MAPPOS$SNP %in% Two_method_OL,]
write.table(Anymethod_OL_Frame, "Charr_GDL_Outliers", sep = "\t", col.names = T, row.names = F, quote = F)
write.table(Important_OL_Frame, "Charr_GDL_Outliers_TOP", sep = "\t", col.names = T, row.names = F, quote = F)
write.table(Two_method_OL_frame, "Charr_GDL_Outliers_2method", sep = "\t", col.names = T, row.names = F, quote = F)

#Plot FST and outliers
FST_MAN <- ggman(FST_MAP, chrom = "Chrom", pvalue = "Fst_WC", snp = "SNP", bp="BP", logTransform = F, pointSize = 1, title = "CHARR", xlabel = "Chromosome", ylabel = "FST", ymax = 1, ymin = 0, sigLine = quantile(GDL_FST$Fst_WC, 0.95, na.rm = T)) + theme_classic()
DIV_SNPS <- as.character(Div_rem$SNP)
ggmanHighlight(FST_MAN, highlight = Two_method_OL)

##Population structure ####

#prepare datasets ####
library(LEA)
library(lfmm)
library(adegenet)
library(genepopedit)
###smnf Admixture ####
###SMNF # Neutral
setwd(dir = "~/Desktop/Charr_Chip_Recalled")
ped2lfmm(input.file = "Charchip_GDL_Neut_Poly_LE12.ped")
pc = pca("Charchip_GDL_Neut_Poly_LE12.lfmm") 
tc = tracy.widom(pc)
plot(tc$percentage)
project2 = NULL
project2 = snmf(input.file = "Charchip_GDL_Neut_Poly_LE12.lfmm", K = 1:5, entropy = TRUE, project = "new")
plot(project2, col = "blue", pch = 19, cex = 1.2)

best = which.min(cross.entropy(project2, K = 2))
my.colors <- c( "burlywood4", "antiquewhite")
barchart(project2, K = 2, run = best,
         border = NA, space = 0,
         col = my.colors,
         xlab = "Individuals",
         ylab = "Ancestry proportions",
         main = "Ancestry matrix") -> bp
axis(1, at = 1:length(bp$order),
     labels = GDL_PHENO$Phenotype[bp$order], las=1,
     cex.axis = .3)

GDL_PHENO <- read.delim("GDL_PHENO2.txt", header = T, stringsAsFactors = F)
GDL_IDs <-as.data.frame(system("awk '{print $2}' Charchip_GDL_Neut_Poly_LE12.ped", intern = T), stringsAsFactors = F)
colnames(GDL_IDs) <- "CHIP_ID"
GDL_PHENO <- plyr::join(GDL_IDs, GDL_PHENO)
GDL_PHENO$Phenotype[bp$order]

###SMNF #Div
setwd(dir = "~/Desktop/Charr_Chip_Recalled")
ped2lfmm(input.file = "Charchip_GDL_Div12.ped")
pc = pca("Charchip_GDL_Div12.lfmm")
tc = tracy.widom(pc)
plot(tc$percentage)
project = NULL
project = snmf(input.file = "Charchip_GDL_Div12.lfmm", K = 1:5, entropy = TRUE, project = "new")
plot(project, col = "blue", pch = 19, cex = 1.2)

best = which.min(cross.entropy(project, K = 2))
my.colors <- c( "antiquewhite", "burlywood4")
barchart(project, K = 2, run = best,
         border = NA, space = 0,
         col = my.colors,
         xlab = "Individuals",
         ylab = "Ancestry proportions",
         main = "Ancestry matrix") -> bp
axis(1, at = 1:length(bp$order),
     labels = GDL_PHENO$Phenotype[bp$order], las=1,
     cex.axis = .3)




##PCA
setwd("~/Desktop/Charr_Chip_Recalled/")
library(adegenet)
library(genepopedit)
library(ggplot2)
#NEUTRAL

NEU <- read.genepop(file = "Pale_dark_LE.gen", ncode = 3, quiet = FALSE) 
NEU_TAB <- tab(NEU, freq = TRUE, NA.method = "mean")
NEU_ALL <- dudi.pca(NEU_TAB,scale=FALSE,scannf=FALSE,nf=3)
NEU_PCframe <- as.data.frame(NEU_ALL$l1)
NEU$pop
NEUPC_Scores <- as.data.frame(cbind(POPS = NEU$pop, NEU_PCframe))
ggplot(NEUPC_Scores, aes(x=RS1, y=RS2, col=NEUPC_Scores$POPS)) + geom_point() + theme_classic() + xlab("PC1") + ylab("PC2") + scale_colour_manual(values = c("antiquewhite4","burlywood"))
loadingplot(abs(NEU_ALL$c1$CS2), lab = NULL , threshold = 0.08, xlab = "SNP", ylab = "PC Loading", main = "PC Loading per SNP",  )
rownames(NEU_ALL$c1[which(abs(NEU_ALL$c1$CS2>0.08))])
#DAPC
grp <- find.clusters(NEU, max.n.clust=10)
#60 PCs retained (~80% cumulative variance)
#2 clusters retained
dapc1 <- dapc(NEU, grp$grp)
#60 PCs retained ~80% cumulative variance
#1 eigenvalues retained
scatter(dapc1)
scatter(dapc1, posi.da="bottomright", bg="white", pch=17:22)
#Plot DAPC by Cluster
scatter(dapc1, scree.da=FALSE, bg="white", pch=20,  cell=0, cstar=0, col=deepseasun(6), solid = .4, cex=3,clab=0, leg=TRUE, txt.leg=paste("Cluster",1:2))
#Plot DAPC by Pop
scatter(dapc1, bg="white", pch=20, leg=TRUE, grp = NEU$pop, clab=0, posi.leg = "topright",  posi.da="topleft")

which(NEU_ALL$c1$CS2>0.08)
#DIVL

ALLPOP <- read.genepop(file = "Pale_dark_DIV.gen", ncode = 3, quiet = FALSE) 
ALL_TAB <- tab(ALLPOP, freq = TRUE, NA.method = "mean")
PCA_ALL <- dudi.pca(ALL_TAB,scale=FALSE,scannf=FALSE,nf=3)
ALLframe <- as.data.frame(PCA_ALL$l1)


ALLPOP$pop


ALLPC_Scores <- as.data.frame(cbind(POPS = ALLPOP$pop, ALLframe))
ggplot(ALLPC_Scores, aes(x=RS1, y=RS2, col=ALLPC_Scores$POPS)) + geom_point() + theme_classic() + xlab("PC1") + ylab("PC2")

#DAPC
grp <- find.clusters(ALLPOP, max.n.clust=10)
#10 PCs retained (~80% cumulative variance)
#2 clusters retained
dapc1 <- dapc(ALLPOP, grp$grp)
#10 PCs retained ~80% cumulative variance
#1 eigenvalues retained
scatter(dapc1)
scatter(dapc1, posi.da="bottomright", bg="white", pch=17:22)
#Plot DAPC by Cluster
scatter(dapc1, scree.da=FALSE, bg="white", pch=20,  cell=0, cstar=0, col=deepseasun(6), solid = .4, cex=3,clab=0, leg=TRUE, txt.leg=paste("Cluster",1:6))
#Plot DAPC by Pop
scatter(dapc1, , bg="white", pch=20, leg=TRUE, grp = ALLPOP$pop, clab=0, posi.leg = "topright",  posi.da="topleft")

library(diveRsity)
divMigrate(infile = "DIVFinal.gen", outfile = NULL, boots = 0, stat = "Nm",  filter_threshold = 0.075,   plot_col = "darkblue", para = FALSE, plot_network = TRUE)
divMigrate(infile = "NeutralFinal.gen", outfile = NULL, boots = 0, stat = "Nm",  filter_threshold = 0.075,   plot_col = "darkblue", para = FALSE, plot_network = TRUE)


##tsne ####
GDL_RAW <- read.delim("Charchip_GDL_Neut_Poly_LEA.raw", sep = "", stringsAsFactors = F)
GDL_SNPS <- GDL_RAW[7:2436]
GDL_RAW$IID
#impute SNPs
library(Rtsne)
GDL_SNPS <- apply(GDL_SNPS, 2, function(x) replace(x, is.na(x), as.numeric(names(which.max(table(x))))))
testne <- Rtsne(as.matrix(GDL_SNPS), num_threads = 16, perplexity = 25 )
tsne_SNP <- data.frame(testne$Y)
GDL_PHENO <- read.delim("GDL_PHENO2.txt", header = T, stringsAsFactors = F)

tsne_frame <- cbind(tsne_SNP, GDL_RAW[2])
GDL_tsne_SCORES <- merge(GDL_PHENO, tsne_frame, by.x = "CHIP_ID", by.y = "IID" )
ggplot(GDL_tsne_SCORES, aes(x=X1, y=X2, col=Phenotype)) + geom_point() + theme_classic() + xlab("TSNE1") + ylab("TSNE2")


#tsne - div 
Charchip_GDL_DivA

GDLdiv_RAW <- read.delim("Charchip_GDL_DivA.raw", sep = "", stringsAsFactors = F)
GDLdiv_SNPS <- GDL_RAW[7:205]
GDLdiv_RAW$IID
#impute SNPs
library(Rtsne)
GDLdiv_SNPS <- apply(GDLdiv_SNPS, 2, function(x) replace(x, is.na(x), as.numeric(names(which.max(table(x))))))
testne_div <- Rtsne(as.matrix(GDLdiv_SNPS), num_threads = 16, perplexity = 25 )
tsnediv_SNP <- data.frame(testne_div$Y)
GDL_PHENO <- read.delim("GDL_PHENO2.txt", header = T, stringsAsFactors = F)

tsnediv_frame <- cbind(tsnediv_SNP, GDLdiv_RAW[2])
GDLdiv_tsne_SCORES <- merge(GDL_PHENO, tsnediv_frame, by.x = "CHIP_ID", by.y = "IID" )
ggplot(GDLdiv_tsne_SCORES, aes(x=X1, y=X2, col=Phenotype)) + geom_point() + theme_classic() + xlab("TSNE1") + ylab("TSNE2")

