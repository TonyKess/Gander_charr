library(RcppCNPy)
library(ggplot2)
library(data.table)
library(dplyr)
library(windowscanr)
setwd(dir = "~/Desktop/Tony/Charr")

#Get chromosome information
CHar_LG_Chrom_Conversion <- read.delim("~/Desktop/Tony/Charr/CHar_LG_Chrom_Conversion.txt", stringsAsFactors=FALSE, header = F)
colnames(CHar_LG_Chrom_Conversion) <- c("LG", "CM_Chrom", "NC_Chrom")

#Import and update pool WG data from 1K window scan
paledark_fst <- read.delim("~/Desktop/Tony/Charr/paledark.fst", header=FALSE, stringsAsFactors=FALSE)
colnames(paledark_fst) <- c("CM_Chrom", "BP", "SNPs", "Coverage_proportion", "mean_coverage", "pops_compared", "FST")
paledark_fst <- inner_join(CHar_LG_Chrom_Conversion, paledark_fst)

#Get outliers - top 99%
WGpool_99OL <- paledark_fst[which(paledark_fst$FST > quantile(paledark_fst$FST, 0.99)),]

#make a bed file for gene ontology of FST outliers
write.table(cbind(WGpool_99OL$NC_Chrom, WGpool_99OL$BP, (WGpool_99OL$BP +1)), "Charr_WGOL99.bed", sep = "\t", quote = F, row.names = F, col.names = F)

#Make additional dummy windows for windowscanr to work properly
WGpool_99OL_maxes <- WGpool_99OL %>% 
  group_by(LG) %>%
  filter(BP == max(BP))
WGpool_99OL_maxes <- data.frame(WGpool_99OL_maxes, stringsAsFactors = F)
WGpool_99OL_maxes$BP <- as.integer(WGpool_99OL_maxes$BP)
WGpool_99OL_maxes$BP <- WGpool_99OL_maxes$BP + 200000
class(WGpool_99OL$FST)
WGpool_99OL_dummy <- rbind(WGpool_99OL_maxes, WGpool_99OL)

#100K windows - chip SNPs and VST outliers are large distances apart (800,000 and 230,000 bp on average)
WGpool_99OL_windo100K <-  winScan(x = WGpool_99OL_dummy, groups  = "LG", position = "BP", win_size = 100000, win_step = 100000,  values = "FST", funs = "mean")
WGpool_99OL_windo100K <- WGpool_99OL_windo100K[WGpool_99OL_windo100K$FST_n >0,]
WGpool_99OL_windo100K$BP <- WGpool_99OL_windo100K$win_mid

#PCangsd selscan 
library(RcppCNPy)
library(qvalue)
S <- npyLoad("PC_charr_mq30.selection.npy") # Reads results from selection scan
sites <- fread("PC_charr_mq30.sites", data.table = F, stringsAsFactors = F, header = F)
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
pcangsd_sel_1k_99OL<- pcangsd_sel_1kwin[which(pcangsd_sel_1kwin$PC_S_mean > quantile(pcangsd_sel_1kwin$PC_S_mean, 0.99, na.rm=T)),]
pcangsd_sel_1k_99OL$BP <- pcangsd_sel_1k_99OL$win_mid

pcangsd_1kOL_LGBP<- data.frame(cbind("LG" = pcangsd_sel_1k_99OL$LG, "BP" = pcangsd_sel_1k_99OL$BP), stringsAsFactors = F)
pcangsd_1kOL_LGBP$BP<- as.numeric(pcangsd_1kOL_LGBP$BP)

WG_1kOL_LGBP<- data.frame(cbind("LG" = WGpool_99OL$LG, "BP" = WGpool_99OL$BP), stringsAsFactors = F)
WG_plus_PCA_1K_OL <- distinct(rbind(pcangsd_1kOL_LGBP, WG_1kOL_LGBP))

FST_and_PCA <- inner_join(pcangsd_sel_1kwin, paledark_fst, by = c("LG", "BP"))
cor.test(FST_and_PCA$PC_S_mean, FST_and_PCA$FST)
FST_and_PCA_OL <- inner_join(pcangsd_sel_1k_99OL, WGpool_99OL, by = c("LG", "BP"))
FST_and_PCA_OL_forwindo <- data.frame(cbind("LG" = FST_and_PCA_OL$LG, "BP" = FST_and_PCA_OL$BP, "FST" = FST_and_PCA_OL$FST), stringsAsFactors = F)
FST_and_PCA_OL_forwindo$BP <- as.numeric(FST_and_PCA_OL_forwindo$BP)
FST_and_PCA_OL_forwindo$FST <- as.numeric(FST_and_PCA_OL_forwindo$FST)

FST_and_PCA_OL_forwindo_maxes <- FST_and_PCA_OL_forwindo %>% 
  group_by(LG) %>%
  filter(BP == max(BP))
FST_and_PCA_OL_forwindo_maxes <- data.frame(FST_and_PCA_OL_forwindo_maxes, stringsAsFactors = F)
FST_and_PCA_OL_forwindo_maxes$BP <- as.integer(FST_and_PCA_OL_forwindo_maxes$BP)
FST_and_PCA_OL_forwindo_maxes$BP <- FST_and_PCA_OL_forwindo_maxes$BP + 200000
FST_and_PCA_OL_forwindo_dummy <- rbind(FST_and_PCA_OL_forwindo_maxes, FST_and_PCA_OL_forwindo)
FST_and_PCA_OL_100kwin <- winScan(x =FST_and_PCA_OL_forwindo_dummy,  groups  = "LG", position = "BP", win_size = 100000, win_step = 100000,  values = "FST", funs = "mean")
FST_and_PCA_OL_100kwin$BP <- FST_and_PCA_OL_100kwin$win_mid
FST_and_PCA_OL_100kwin<- FST_and_PCA_OL_100kwin[which(FST_and_PCA_OL_100kwin$FST_n > 0),]


#make a bed file for gene ontology
write.table(cbind(FST_and_PCA_OL$NC_Chrom, FST_and_PCA_OL$BP, (FST_and_PCA_OL$BP +1)), "Charr_WGOL_FST_and_PCA_OL99.bed", sep = "\t", quote = F, row.names = F, col.names = F)


pcangsd_1kOL_LGBP<- data.frame(cbind("LG" = pcangsd_sel_1k_99OL$LG, "BP" = pcangsd_sel_1k_99OL$BP, "PC_S_mean" = as.numeric(pcangsd_sel_1k_99OL$PC_S_mean)), stringsAsFactors = F)
pcangsd_1kOL_LGBP$PC_S_mean<- as.numeric(pcangsd_1kOL_LGBP$PC_S_mean)
pcangsd_1kOL_LGBP$BP<- as.numeric(pcangsd_1kOL_LGBP$BP)

pcangsd_sel_1kmaxes <- pcangsd_1kOL_LGBP %>% 
  group_by(LG) %>%
  filter(BP == max(BP))
pcangsd_sel_1kmaxes <- data.frame(pcangsd_sel_1kmaxes, stringsAsFactors = F)
pcangsd_sel_1kmaxes$BP <- as.integer(pcangsd_sel_1kmaxes$BP)
pcangsd_sel_1kmaxes$BP <- pcangsd_sel_1kmaxes$BP + 200000
pcangsd_sel_1k_99OLdummy <- rbind(pcangsd_sel_1kmaxes, pcangsd_1kOL_LGBP)
pcangsd_sel_100kwin <- winScan(x =pcangsd_1kOL_LGBP,  groups  = "LG", position = "BP", win_size = 100000, win_step = 100000,  values = "PC_S_mean", funs = "mean")
pcangsd_sel_100kwin$BP <- pcangsd_sel_100kwin$win_mid
pcangsd_sel_100kwin <- pcangsd_sel_100kwin[which(pcangsd_sel_100kwin$PC_S_mean_n > 0),]

#Import chip map
Charrchip_map <- read.delim("~/Desktop/Tony/Charr/Charrchip_Poly_GDL_comp.map", header=FALSE, stringsAsFactors=FALSE)
colnames(Charrchip_map) <- c("NC_Chrom", "SNP", "whatev", "BP")
CHIP_FST <- read.delim("~/Desktop/Tony/Charr/FST_GDLCharr_diveRsity.txt", stringsAsFactors=FALSE)
CHIP_FST <- data.frame(cbind(CHIP_FST$loci, CHIP_FST$Fst_WC),stringsAsFactors = F)
colnames(CHIP_FST) <- c("SNP", "FST")
CHIP_FST <- inner_join(CHIP_FST, Charrchip_map)
CHIP_FST <- inner_join(CHIP_FST, CHar_LG_Chrom_Conversion)
CHIP_FST$FST <-  as.numeric(CHIP_FST$FST)

#Chip outliers
smnfOL <- read.delim("~/Desktop/Tony/Charr/smnf_OL.txt", stringsAsFactors = F)
pcadapt_OL <- read.delim("~/Desktop/Tony/Charr/PCAdapt_CHARRGDL_Poly_Outliers_K1_Q05.txt", stringsAsFactors = F)
Chip_OL <- inner_join(smnfOL, pcadapt_OL)
Chip_OL <- CHIP_FST[CHIP_FST$SNP %in% Chip_OL$SNP,]
library(ggman)
CHIP_FSTMAN <- ggman(gwas = CHIP_FST, snp = "SNP", chrom = "LG", pvalue = "FST", bp = "BP", logTransform = FALSE, ymax = 1,  pointSize = 2) + theme_classic()
ggmanHighlight(CHIP_FSTMAN, highlight = Chip_OL$SNP)
#add fake position at 100kb for the max. 
Chip_OL_maxes <- Chip_OL %>% 
  group_by(LG) %>%
  filter(BP == max(BP))
Chip_OL_maxes <- data.frame(Chip_OL_maxes, stringsAsFactors = F)
Chip_OL_maxes$BP <- as.integer(Chip_OL_maxes$BP)
Chip_OL_maxes$BP <- Chip_OL_maxes$BP + 200000
Chip_OL_dummy <- rbind(Chip_OL_maxes, Chip_OL)


CHIP_FST_100KWIN <- winScan(x = Chip_OL_dummy, groups  = "LG", position = "BP", win_size = 100000, win_step = 100000,  values = "FST", funs = "mean")
CHIP_FST_100KWIN <- CHIP_FST_100KWIN[CHIP_FST_100KWIN$FST_n > 0,]
CHIP_FST_100KWIN$BP <- CHIP_FST_100KWIN$win_mid

#Check overlap of 100kB outlier windows on chip with WG outlier windows
chip_WG_comp <- inner_join(CHIP_FST_100KWIN, WGpool_99OL_windo100K, by = c("LG" ,"BP"))
chip_PC_comp <- inner_join(CHIP_FST_100KWIN, pcangsd_sel_100kwin, by = c("LG" ,"BP"))
chip_OL_comp <- inner_join(CHIP_FST_100KWIN, FST_and_PCA_OL_100kwin, by = c("LG" ,"BP"))



#VST OL
GDL_pale_GDL_dark_CNV_OL <- read.delim("~/Desktop/Tony/Charr/RDA_and_VST_OL.txt", header=T, stringsAsFactors=FALSE)
colnames(GDL_pale_GDL_dark_CNV_OL)[3] <- "BP"

GDL_VST <- read.delim("~/Desktop/Tony/GDL_dark_rn_GDL_pale_rn.Vst.txt", header=FALSE, stringsAsFactors=FALSE)
colnames(GDL_VST) <- c("SNP", "LG_alt", "BP", "C", "VST")
CHar_LG_Chrom_Conversion_VST <- read.table("~/Desktop/Tony/Charr/CHar_LG_Chrom_Conversion_VST.txt", quote="\"", comment.char="", stringsAsFactors=FALSE)
colnames(CHar_LG_Chrom_Conversion_VST) <- c("LG", "LG_alt", "CM_Chrom", "NC_Chrom")
GDL_VST <- inner_join(GDL_VST, CHar_LG_Chrom_Conversion_VST)
GDL_VST_OL <- inner_join(GDL_pale_GDL_dark_CNV_OL, GDL_VST)


GDL_VST_OL_maxes <- GDL_VST_OL %>% 
  group_by(LG) %>%
  filter(BP == max(BP))
GDL_VST_OL_maxes <- data.frame(GDL_VST_OL_maxes, stringsAsFactors = F)
GDL_VST_OL_maxes$BP <- as.integer(GDL_VST_OL_maxes$BP)
GDL_VST_OL_maxes$BP <- GDL_VST_OL_maxes$BP + 200000
GDL_VST_OL_dummy <- rbind(GDL_VST_OL_maxes, GDL_VST_OL)


CHIP_VST_WIN <- winScan(x = GDL_VST_OL_dummy, groups  = "LG", position = "BP", win_size = 100000, win_step = 100000,  values = "VST", funs = "mean")
CHIP_VST_WIN <- CHIP_VST_WIN[CHIP_VST_WIN$VST_n >0,]
cor.test(chip_FSTVST_comp$VST_mean, chip_FSTVST_comp$FST_mean)
chipVST_PCA_comp <- inner_join(pcangsd_sel_100kwin, CHIP_VST_WIN, by = c("LG" ,"win_mid"))
cor.test(chipVST_PCA_comp$VST_mean, chipVST_PCA_comp$PC_S_mean_mean)
chip_FSTVST_comp <- inner_join(WGpool_99OL_windo100K, CHIP_VST_WIN, by = c("LG" ,"win_mid"))
chipVST_OLWGcomp <- inner_join(FST_and_PCA_OL_100kwin, CHIP_VST_WIN, by = c("LG" ,"win_mid"))

cor.test(chip_FSTVST_comp$FST_mean, chip_FSTVST_comp$VST_mean)

library(data.table)


#Prepare for WG plot
paledark_fst$LG <- factor(paledark_fst$LG, levels = CHar_LG_Chrom_Conversion$LG)
FST_and_PCA_OL$LG <- factor(FST_and_PCA_OL$LG, levels = CHar_LG_Chrom_Conversion$LG)
Chip_OL$LG <- factor(Chip_OL$LG, levels = CHar_LG_Chrom_Conversion$LG)
GDL_VST_OL$LG <- factor(GDL_VST_OL$LG, levels = CHar_LG_Chrom_Conversion$LG)
GDL_VST_OL$FST <- GDL_VST_OL$VST


ggplot() + facet_wrap(~LG, scales = "free_x") + 
  geom_point(data = FST_and_PCA_OL, aes(x = BP, y = FST, colour = "WG_OL")) + 
  geom_point(data = Chip_OL, aes(x = BP, y = FST, colour = "Chip_OL")) + 
  geom_point(data = GDL_VST_OL, aes(x = BP, y = FST, colour = "VST_OL")) + 
  geom_smooth(data = paledark_fst,  aes(x = BP, y = FST, colour = "WG"), method = "loess", span = 0.01, se = F) +
  theme_classic() 
