library(rsed)
library(data.table)
library(tidyverse)
setwd("~/Desktop/Charr_Reanalysis/")
#turns into Obj_Info - genes in human that will be merged with positions genotyped
HUMAN_genes <- read.delim("~/Desktop/Charr_Reanalysis/Human_genes.txt")
HUMAN_genes <- HUMAN_genes[HUMAN_genes$tax_id %in% "9606",] # change silly naming
charr_genome_annotation <- read.csv("~/Desktop/Charr_Reanalysis/charr_genome_annotation.csv")
CHar_LG_Chrom_Conversion <- read.delim("~/Desktop/Charr_Final_Scripts_Data/CHar_LG_Chrom_Conversion.txt", stringsAsFactors=FALSE, header = F)
colnames(charr_genome_annotation)[11] <- "LG"
colnames(CHar_LG_Chrom_Conversion) <- c("LG", "CM_Chrom", "NC_Chrom")
charr_genome_annotation <- inner_join(charr_genome_annotation, CHar_LG_Chrom_Conversion,)
charr_genome_annotation$Symbol <- toupper(charr_genome_annotation$Symbol)
charr_human_conserved <- inner_join(charr_genome_annotation, HUMAN_genes, by = c("Symbol","description"))
##random checks of these on orthodb match too

write.table(cbind(charr_human_conserved$CM_Chrom, charr_human_conserved$start_position_on_the_genomic_accession.x, charr_human_conserved$end_position_on_the_genomic_accession.x, charr_human_conserved$Symbol),
            "charr_human_conserved_genes.bed", sep = "\t", quote = F, col.names = F, row.names = F)

Genes_only_bed <- data.frame(cbind(charr_human_conserved$CM_Chrom, charr_human_conserved$start_position_on_the_genomic_accession.x, charr_human_conserved$end_position_on_the_genomic_accession.x, charr_human_conserved$Symbol))
colnames(Genes_only_bed) <- c("CM_chrom", "start", "stop", "Symbol")


#Make sure there is no missing/nas in any of the gene sets otherwise bedtools will blow up
Genes_only_bed <- Genes_only_bed[(which(!(is.na(Genes_only_bed$X2)))),]

#Get some kind of summary statistic
#WG_FST
paledark_fst <- read.delim("~/Desktop/Charr_Reanalysis/Old Charr files/Charr_Linux(1)/paledark.fst", header=FALSE)
colnames(paledark_fst) <- c("CM_Chrom", "BP", "SNPs", "Coverage_proportion", "mean_coverage", "pops_compared", "FST")
paledark_fst <- inner_join(CHar_LG_Chrom_Conversion, paledark_fst)
pale_dark_fst_bed<- data.frame(cbind(paledark_fst$CM_Chrom, paledark_fst$BP , paledark_fst$BP + 1, paledark_fst$FST))
write.table(pale_dark_fst_bed, "paledark_fst.bed", col.names = F, row.names = F, sep = "\t", quote = F)


#Overlap using bedtools intersect
#If these files are too large intersect can be done by chromosome and then the files can be concatenated using "cat".
#Also a good way to check where the bed file is screwy
system("bedtools intersect -b paledark_fst.bed -a charr_human_conserved_genes.bed -wb > paledark_fst_human_conserved_genes.txt")
paledark_fst_conserved_genes <- fread("~/Desktop/Charr_Reanalysis/paledark_fst_human_conserved_genes.txt", stringsAsFactors = F, data.table = F )
paledark_fst_conserved_genes <- data.frame(cbind(paledark_fst_conserved_genes$V1,paledark_fst_conserved_genes$V6, 
                                                 paledark_fst_conserved_genes$V7, paledark_fst_conserved_genes$V8, paledark_fst_conserved_genes$V4))
colnames(paledark_fst_conserved_genes) <- c("Chromosome", "Position", "End", "FST", "Symbol")


#Get maximum FST

paledark_fst_conserved_genes <- paledark_fst_conserved_genes %>% 
  group_by(Symbol) %>% 
  filter(FST == max(FST))


#remove duplicates
paledark_fst_conserved_genes[which(duplicated(paledark_fst_conserved_genes$Symbol)),]
paledark_fst_conserved_genes <- as.data.frame(unique(setDT(paledark_fst_conserved_genes), by = "Symbol"))


Obj_Info <- inner_join(paledark_fst_conserved_genes, charr_human_conserved)
Obj_Info$SNPcount <- 1
colnames(Obj_Info)
Obj_Info$GeneLength <- Obj_Info$end_position_on_the_genomic_accession.x - Obj_Info$start_position_on_the_genomic_accession.x
#Make an object info file with charr positions and FST, but used humann gene IDs
Obj_Info <- data.frame(cbind(Obj_Info$GeneID.y, Obj_Info$FST, Obj_Info$Symbol, Obj_Info$GeneLength, Obj_Info$CM_Chrom, Obj_Info$start_position_on_the_genomic_accession.x, Obj_Info$end_position_on_the_genomic_accession.x, Obj_Info$orientation.x), stringsAsFactors = F)
colnames(Obj_Info) <- c("objID",	"objStat",	"objName",	"GeneLength",	"chr",	"startpos",	"endpos",	"strand")
Obj_Info$strand <- NA
Obj_Info <- Obj_Info[Obj_Info$chr %in% CHar_LG_Chrom_Conversion$CM_Chrom,]

#SetInfo - TBD
human_kegg <- read.delim("~/Desktop/Charr_Reanalysis/human_kegg.txt")
colnames(human_kegg )
human_kegg  <- data.frame(cbind(human_kegg$BSID, human_kegg$Name, human_kegg$Source),stringsAsFactors = F)

SetInfo <- human_kegg
colnames(SetInfo) <- c("setID", "setName", "setSource")


#SetObj
#turns into SetObj after filtering by GeneIDs in Salmon dataset
SetObj <- read.delim("~/Desktop/Projects/Salmon_AOM/biosystems_gene", header=FALSE)
colnames(SetObj) <- c("setID",	"objID")
SetObj <- SetObj[1:2]

SetObj <- SetObj[SetObj$setID %in% SetInfo$setID,]
colnames(SetObj)
colnames(Obj_Info)


#ensure consistency of genes, pathways across comparisons
Obj_Info <- Obj_Info[Obj_Info$objID  %in% SetObj$objID,]
SetObj <- SetObj[SetObj$objID  %in% Obj_Info$objID,]
SetInfo <- SetInfo[SetInfo$setID %in% SetObj$setID,]

#Just go through and ensure class consistency for all datasets compared to PolyLinkR example
Obj_Info$objID <- as.integer(Obj_Info$objID)
Obj_Info$objStat <- as.numeric(Obj_Info$objStat)
Obj_Info$GeneLength <- as.integer(Obj_Info$GeneLength)
Obj_Info$startpos <- as.integer(Obj_Info$startpos)
Obj_Info$endpos <- as.integer(Obj_Info$endpos)
setDT(Obj_Info)

SetInfo$setID <- as.integer(SetInfo$setID)
setDT(SetInfo)

SetObj$setID 

class(SetObj$objID)
setDT(SetObj)
class(SetObj)
colnames(Anatolia_EF_CLR)
colnames(Obj_Info)
library(PolyLinkR)
output_FST = polylinkr(obj.info = Obj_Info, set.info = SetInfo, set.obj = SetObj,
                       n.cores = 8, emp.nruns = 10000, NN = 1000)
write.table(output_FST, "polylinkR_human_FST.txt", col.names = T, row.names = F, sep = "\t", )
library(qvalue)
qvs <- qvalue(output_FST$setP)
qvs$qvalues

write.table(Obj_Info, "~/Desktop/Software/polysel/data/charmander/ObjInfo.txt", col.names = T, row.names = F, sep = "\t", quote = F)
write.table(SetInfo, "~/Desktop/Software/polysel/data/charmander/SetInfo.txt", col.names = T, row.names = F, sep = "\t", quote = F)
write.table(SetObj, "~/Desktop/Software/polysel/data/charmander/SetObj.txt", col.names = T, row.names = F, sep = "\t", quote = F)
