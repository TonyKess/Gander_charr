library(data.table)
#REFORMAT PLINK - PED ####
##FORMAT WEIRD AXIOM OUTPUT "PLINK" FILE TO AN ACTUAL PLINK FILE

#Change character spacing from spaces "" to tabs "/t", ditch annoying header
#this can be achieved in TextWrangler/TextEdit, or via command line using sed (or gsed on mac): 
system("cd Desktop; 
       gsed  's/ /\t/g' Axiom_SalpSNP_Genotypes_05nov18.ped.txt > Charrchip_Recalled.ped; sed -i.foo '1d' Charrchip_Recalled.ped")
# Make a file of IDs that don't have extraneous axiom file information, with individual names that will work in genepopedit
system("cd Desktop; awk '{print $1}' Charrchip_Recalled.ped | sed 's/.CEL//' | sed 's/.*_.*_.*_.*_.*_//' | sed 's/./&_/3' > Charr_IDs_fixed")


#Read in using datatable
charchippy <- fread(input = '~/Desktop/Charr_Chip_Recalled/Charrchip_Recalled.ped', stringsAsFactors = F, data.table = F)

#Get IDs
IDs <- read.delim("~/Desktop/Charr_Chip_Recalled/Charr_IDs_fixed", header = F, stringsAsFactors = F)
#Get just SNPs
charchippy_justSNPs <- charchippy[c(-1)]
#Build plink prefix for ped format
OOs <- rep(0, 841)
POPS <- gsub("_.*", "", IDs$V1)
plinkpref <- as.data.frame(cbind(POPS, IDs = as.character(IDs$V1), OOs, OOs, OOs, OOs))
#Combine to make ped format
charchip_plink <- data.frame(cbind(plinkpref, charchippy_justSNPs))
#Check
charchip_plink[1:7]

##write out
fwrite(charchip_plink, "~/Desktop/Charr_Chip_Recalled/CHARR CHIP/Charrchip_Recalled.ped", col.names = F, row.names = F, sep = "\t", quote = F)



#REFORMAT PLINK - MAP ####
##Reformat with just O's first for filtering - plink can't handle the number of chromosomes present in raw file 
system("cd Desktop/Charr_Chip_Recalled; 
       sed '1,4d'  Axiom_SalpSNP_Genotypes_05nov18.map.txt > Charrchip_Recalled.map")
charmap <- fread(input = '~/Desktop/Charrchip_Recalled.map', stringsAsFactors = F, data.table = F)

charmap_0s <- data.frame(cbind(rep(0, 86504), charmap $`Marker ID`, rep(0, 86504), rep(0, 86504)), stringsAsFactors = F )
fwrite(charchip_plink, "~/Desktop/Charr_Chip_Recalled/Charrchip_Recalled_tofilter.ped", col.names = F, row.names = F, sep = "\t", quote = F)

fwrite(charmap_0s, "~/Desktop/Charr_Chip_Recalled/Charrchip_Recalled_tofilter.map", col.names = F, row.names = F, sep = "\t", quote = F)


#Use plink for to filter for missing greater than 5%, minor allele less than 1%
system("cd Desktop/Charr_Chip_Recalled; plink_mac/plink --file Charrchip_Recalled_tofilter --maf .01 --geno 0.05 --recode --out Charrchip_recall_filter01  ")

#Get physical position info
#filtered file
filtered_map <- fread("~/Desktop/Charr_Chip_Recalled/Charrchip_recall_filter01.map", data.table = F)

#mapped to contigs and chromosomes
mapped_loci <- fread("~/Desktop/Charr_Chip_Recalled/CHARR CHIP/Axiom_SalpSNP_Annotation_wMAP.r1.csv", data.table = F)

mapped_all <- data.frame(cbind(mapped_loci$Chromosome, mapped_loci$Probe.Set.ID, rep(0, 86504), mapped_loci$Physical.Position), stringsAsFactors = F )
mapped_to_chrom <- data.frame(cbind(mapped_loci$Chromosome[grep("NC",mapped_loci$Chromosome)], mapped_loci$Probe.Set.ID[grep("NC",mapped_loci$Chromosome)], rep(0, 58494), mapped_loci$Physical.Position[grep("NC",mapped_loci$Chromosome)]),stringsAsFactors = F)

##Arrange SNPs by chromosome number and physical position
library(genepopedit)
library(dplyr)

#Make a new plink map file with filtered SNPs mapped to chromosomes
Mapped_filtered <- mapped_to_chrom[mapped_to_chrom$X2 %in% filtered_map$V2,]
Mapped_filtered$X4 <- as.numeric(Mapped_filtered$X4)
Arrange_mapped2chrom <- arrange(Mapped_filtered,X1,X4)   
fwrite(Arrange_mapped2chrom, "~/Desktop/Charr_Chip_Recalled/Charrchip_recall_filter01_mappedordered.map", col.names = F, row.names = F, sep = "\t", quote = F)

##First, convert to genepop using PDGSPIDER#
#Use genepopedit to reorder SNPS

subset_genepop(genepop = "~/Desktop/Charr_Chip_Recalled/Charrchip_recall_filter01.txt", subs = Arrange_mapped2chrom$X2, keep = T, path = paste0("~/Desktop/Charr_Chip_Recalled/Charrchip_recall_filter01_mappedordered.txt"))

#Convert back to plink in PGDSPider

#Convert genotypes recognizeable by pcadapt
GP <- fread("~/Desktop/Charr_Chip_Recalled/Charrchip_recall_filter01_mappedordered.ped",  data.table = F)
JustLoci <- as.data.frame(GP[,-c(1:6)])
JustLoci <- sapply(JustLoci, as.character)
JustLoci[JustLoci=="1"] <- "A" ##genepop code 100 (switch out variables depending on allele coding)
JustLoci[JustLoci=="4"] <- "T" ##genepop code 110
JustLoci[JustLoci=="3"] <- "G" ##genepop code 120
JustLoci[JustLoci=="2"] <- "C" ##genepop code 130
GP[1] <- gsub("_.*", "", GP$V2)
catframe <- cbind(GP[,(1:6)], JustLoci)
fwrite(catframe, "~/Desktop/Charr_Chip_Recalled/Charrchip_recall_filter01_mappedordered.ped", sep="\t", col.names = FALSE, quote = FALSE, row.names = FALSE)

#Filter by population ####
Charrchip_recall_noGDL <- catframe[which(catframe$V1!="GDL"),]
fwrite(Charrchip_recall_noGDL, "~/Desktop/Charr_Chip_Recalled/Charrchip_recall_noGDL.ped", sep="\t", col.names = FALSE, quote = FALSE, row.names = FALSE)


Charrchip_recall_GDL <- catframe[which(catframe$V1=="GDL"),]
fwrite(Charrchip_recall_GDL, "~/Desktop/Charr_Chip_Recalled/Charrchip_recall_GD.ped", sep="\t", col.names = FALSE, quote = FALSE, row.names = FALSE)


#Population structure test ####
library(pcadapt)
library(ggplot2)
#working directory must be same as file locations for pcadapt
setwd("~/Desktop/Charr_Chip_Recalled/")
CHARR_ALL <- read.pcadapt("Charrchip_recall_filter01_mappedordered.ped", type = "ped")
PCs <- pcadapt(CHARR_ALL, K=5, method="mahalanobis",  min.maf = 0.01)
plot(PCs,option="screeplot")

#Plot population structure
ID <- as.character(GP$V2)
POPS <- as.character(GP$V1)
ALLframe <- as.data.frame(PCs$scores)
ALLPC_Scores <- as.data.frame(cbind(POPS, ALLframe))
ALLPC_Scores$POPS <- as.character(ALLPC_Scores$POPS)
ALLPC_Scores$POPS <- factor(ALLPC_Scores$POPS, levels=unique(ALLPC_Scores$POPS))
ggplot(ALLPC_Scores, aes(x=V1, y=V2, col=ALLPC_Scores$POPS)) + geom_point() + theme_classic() + xlab("PC1") + ylab("PC2")

MAPPOS <- read.delim2("Charrchip_recall_filter01_mappedordered.map", header = F, stringsAsFactors = F)
colnames(MAPPOS) <- c("Chrom", "SNP", "CM", "BP")
PVALS <- PCs$pvalues
PCMAP <- as.data.frame(cbind(MAPPOS, PVALS))
library(qvalue)
qval <- qvalue(PCs$pvalues)$qvalues
alpha <- 0.05
outliers <- which(qval<alpha)
Outliers_CHARR <- PCMAP[outliers,]
write.table(Outliers_CHARR, "PCAdapt_CHARR_ALL_Outliers_K5_Q05.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
class(Outliers_CHARR$PVALS)
max(Outliers_CHARR$PVALS)
#Manhattan plot of individual SNP scores
library(ggman)
ggman(PCMAP, chrom = "Chrom", pvalue = "PVALS", snp = "SNP", bp="BP", pointSize = 1, title = "CHARR", xlabel = "Chromosome", sigLine = -log10(max(Outliers_CHARR$PVALS)),ymax = 500) + theme_classic()

##PSV Filter ####

#not in function
'%!in%' <- function(x,y)!('%in%'(x,y)) # not in

#Import Amber's Axiom call summary
M_SNPstats <- read.delim("ProbeSetSummaryTable_05nov18.txt", stringsAsFactors = F)
unique(M_SNPstats$ConversionType)

#Polymorphic High Resolution SNPs
PolyHigh_M <- data.frame(M_SNPstats$probeset_id[which(M_SNPstats$ConversionType=="PolyHighResolution")])
colnames(PolyHigh_M) <- "SNP"
PolyHigh_M$SNP <- as.character(PolyHigh_M$SNP)

#Sarah's potential PSV list
S_SNPstats <- read.delim(" ", header = F, stringsAsFactors = F)

#SNPs that pass both quality checks
Good_SNPS <- data.frame(PolyHigh_M$SNP[PolyHigh_M$SNP %!in% S_SNPstats$V1])
colnames(Good_SNPS) <- "SNP"
Good_SNPS$SNP <- as.character(Good_SNPS$SNP)
SNPS_PolyHigh_Everywhere <- MAPPOS[MAPPOS$SNP %in% Good_SNPS$SNP,] 
SNPS_PolyHigh_Everywhere <- SNPS_PolyHigh_Everywhere[2]
write.table(SNPS_PolyHigh_Everywhere, "PolyHighSNPS.txt", col.names = F, row.names = F, quote = F,sep = "\t" )


#Filter for PSV in plink
system("cd Desktop/Charr_Chip_Recalled; ../plink_mac/plink --file Charrchip_recall_filter01_mappedordered --extract PolyHighSNPS.txt --recode --out Charrchip_recall_filter01_Poly --allow-extra-chr ")

# Population structure - PolyHigh SNPS only ####

#Population structure test ####
library(pcadapt)
library(ggplot2)
#working directory must be same as file locations for pcadapt
setwd("~/Desktop/Charr_Chip_Recalled/")
CHARR_ALL <- read.pcadapt("Charrchip_recall_filter01_Poly.ped", type = "ped")
PCs <- pcadapt(CHARR_ALL, K=5, method="mahalanobis",  min.maf = 0.01)
plot(PCs,option="screeplot")

#Plot population structure
ID <- as.character(GP$V2)
POPS <- as.character(GP$V1)
ALLframe <- as.data.frame(PCs$scores)
ALLPC_Scores <- as.data.frame(cbind(POPS, ALLframe))
ALLPC_Scores$POPS <- as.character(ALLPC_Scores$POPS)
ALLPC_Scores$POPS <- factor(ALLPC_Scores$POPS, levels=unique(ALLPC_Scores$POPS))
ggplot(ALLPC_Scores, aes(x=V1, y=V2, col=ALLPC_Scores$POPS)) + geom_point() + theme_classic() + xlab("PC1") + ylab("PC2")

MAPPOS <- read.delim2("Charrchip_recall_filter01_Poly.map", header = F, stringsAsFactors = F)
colnames(MAPPOS) <- c("Chrom", "SNP", "CM", "BP")
PVALS <- PCs$pvalues
PCMAP <- as.data.frame(cbind(MAPPOS, PVALS))
library(qvalue)
qval <- qvalue(PCs$pvalues)$qvalues
alpha <- 0.05
outliers <- which(qval<alpha)
Outliers_CHARR <- PCMAP[outliers,]
write.table(Outliers_CHARR, "PCAdapt_CHARR_ALL_POLY_Outliers_K5_Q05.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
class(Outliers_CHARR$PVALS)
max(Outliers_CHARR$PVALS)
#Manhattan plot of individual SNP scores
library(ggman)
ggman(PCMAP, chrom = "Chrom", pvalue = "PVALS", snp = "SNP", bp="BP", pointSize = 1, title = "CHARR", xlabel = "Chromosome", sigLine = -log10(max(Outliers_CHARR$PVALS)),ymax = 500) + theme_classic()
