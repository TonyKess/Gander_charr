##RDA #####
library(vegan)
library(rsed)

#genotypes
system("~/Desktop/Software/plink_mac_20200219/plink --file Charrchip_Poly_GDL_comp --recodeA --allow-extra-chr  --out Charrchip_Poly_GDL_comp")

GDL_PHENO <- read.delim("GDL_PHENO2.txt", header = T, stringsAsFactors = F)
GDL_IDs <-as.data.frame(system("awk '{print $2}' Charchip_GDL_Neut_Poly_LE12.ped", intern = T), stringsAsFactors = F)
colnames(GDL_IDs) <- "CHIP_ID"
GDL_PHENO <- plyr::join(GDL_IDs, GDL_PHENO)
GDL_PHENO$Phenotype[bp$order]

GDL_RAW <- read.delim("Charrchip_Poly_GDL_comp.raw", sep = "", stringsAsFactors = F)
colnames(GDL_RAW)[2] <- "CHIP_ID"
GDL_PHENO_RAW <- inner_join(GDL_PHENO, GDL_RAW)

colnames(GDL_PHENO_RAW[1:10])
#subset SNPs
GDL_SNPS <- GDL_PHENO_RAW[10:3653]
#impute SNPs
GDL_SNPS <- apply(GDL_SNPS, 2, function(x) replace(x, is.na(x), as.numeric(names(which.max(table(x))))))
#get morph and year
GDL_MORPH_YEAR <- data.frame(cbind(GDL_PHENO_RAW$Phenotype, GDL_PHENO_RAW$Year))
colnames(GDL_MORPH_YEAR) <- c("Morph", "Year")
#RDA by significant morphometric variables for ecotype differentiation
GDL_SNP_YEAR <- rda(GDL_SNPS ~ GDL_MORPH_YEAR$Morph + Condition(~GDL_MORPH_YEAR$Year))
RsquareAdj(GDL_SNP_YEAR)
VP_GDL_SNP_YEAR <- varpart(GDL_SNPS, ~ GDL_MORPH_YEAR$Morph, ~ GDL_MORPH_YEAR$Year)

showvarparts(parts = 2)
plot(VP_GDL_SNP_YEAR)

GDL_SNP_YEAR <- rda(GDL_SNPS ~ GDL_MORPH_YEAR$Morph * GDL_MORPH_YEAR$Year)
RsquareAdj(GDL_SNP_YEAR)


GDL_YEAR <- rda(GDL_SNPS ~ GDL_MORPH_YEAR$Year)
RsquareAdj(GDL_YEAR)
anova.cca(GDL_YEAR)

summary(eigenvals(GDL_SNP_YEAR, model = "constrained"))
anova.cca(GDL_SNP_YEAR, parallel=8, by = "terms" ) # default is permutation=999

#Get SNP scores
library(rsed)
GDL_SNP_YEAR_scores <- scores(GDL_SNP_YEAR, choices=c(1:3), display="species")  # Species scores for the first three constrained 
GDL_SNP_YEAR_scores <- data.frame(GDL_SNP_YEAR_scores)
SNPnames <- rownames(GDL_SNP_YEAR_scores)
SNPnames  <- sed_substitute(SNPnames , "AX.", "AX-") 
SNPnames  <- sed_substitute(SNPnames , "_.*", "")
GDL_SNP_YEAR_scores$SNPnames <- as.character(SNPnames)

outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]               # locus names in these tails
}

GDL_SNP_YEAR_scores$RDA1abs <- abs(GDL_SNP_YEAR_scores$RDA1)

GDL_SNP_YEAR_scores_OL <- GDL_SNP_YEAR_scores[which(GDL_SNP_YEAR_scores$RDA1abs > quantile(GDL_SNP_YEAR_scores$RDA1abs, 0.99)),]

GDL_SNP_YEAR_scores_OL$SNPnames %in% Chip_OL$SNP
