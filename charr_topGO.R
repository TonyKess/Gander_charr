setwd("~/Desktop/Charr_Reanalysis/gene_ontology/")

library(topGO)

#all gene ontology for salp
go.anno <- read.csv("Unique.Blast2go.Salpinus.annot", header = F, stringsAsFactors = F, sep="\t")
go.anno <- go.anno[1:2]


#Outliers (intersected)
outliers <- read.table("cnvnator_filtered_distinct_XPs", stringsAsFactors = F)[,1]
length(outliers)

#Get outliers in annotated gene ontology
outliers.w.go <- outliers[outliers %in% go.anno[,1]]

#all genes - run same intersect script for bedtools with all SNPs
all.genes <- read.table("ALL_uniq_XPs", stringsAsFactors = F)[,1]
all.genes.w.go <- all.genes[all.genes %in% go.anno[,1]]

#Create function for topGo
CreateGene2GOMapping <- function(x) {
  nr <- nrow(x)
  map <- list()
  for(r in 1:nr) {
    map[[as.character(x[r, 1])]] <- append(map[[as.character(x[r, 1])]], x[r,2])
  }
  return(map)
}


# create GO mapping, takes about 20 minutes
gene2GO.map <- CreateGene2GOMapping(go.anno)

gene2GO.map <- gene2GO.map[names(gene2GO.map) %in% all.genes.w.go]
gene.list <- factor(as.integer(all.genes.w.go %in% outliers.w.go))
names(gene.list) <- all.genes.w.go

CharrDiv.GOdata <- new("topGOdata",
                       description = "Charr_Div_GDL_Outliers", ontology = "BP",
                       allGenes = gene.list,
                       nodeSize = 5,
                       annot = annFUN.gene2GO, gene2GO = gene2GO.map)

CharrDiv.GOdata.Fisher <- runTest(CharrDiv.GOdata, algorithm = "weight01", statistic = "fisher")
geneData(CharrDiv.GOdata.Fisher)
hist(score(CharrDiv.GOdata.Fisher), 50, xlab = "p-values")

write.table(GenTable(CharrDiv.GOdata, CharrDiv.GOdata.Fisher, topNodes = length(CharrDiv.GOdata.Fisher@score)), 
            "Div_C_top5GO_CharrGDL_cnvnator_distinct_uniq.txt", quote = F, row.names = F, sep = "\t")
