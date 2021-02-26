#CHIPFST_FST_COR_ALL_1Kwindo
Charrchip_map <- read.delim("~/Desktop/Charr_Reanalysis/Charrchip_Poly_GDL_comp.map", header=FALSE, stringsAsFactors=FALSE)
colnames(Charrchip_map) <- c("NC_Chrom", "SNP", "whatev", "BP")
CHIP_downsamp_FST <- read.delim("~/Desktop/Charr_Reanalysis/Poolsamps.fst", stringsAsFactors=FALSE)
CHIP_downsamp_FST <- data.frame(cbind(CHIP_downsamp_FST$SNP, CHIP_downsamp_FST$FST),stringsAsFactors = F)
colnames(CHIP_downsamp_FST) <- c("SNP", "FST")
CHIP_ds_FST <- inner_join(CHIP_downsamp_FST, Charrchip_map)
CHIP_ds_FST <- inner_join(CHIP_ds_FST, CHar_LG_Chrom_Conversion)
CHIP_ds_FST$FST <-  as.numeric(CHIP_ds_FST$FST)

CHIP_ds_FST_maxes <- CHIP_ds_FST %>% 
  group_by(LG) %>%
  filter(BP == max(BP))
CHIP_ds_FST_maxes <- data.frame(CHIP_ds_FST_maxes, stringsAsFactors = F)
CHIP_ds_FST_maxes$BP <- as.integer(CHIP_ds_FST_maxes$BP)
CHIP_ds_FST_maxes$BP <- CHIP_ds_FST_maxes$BP + 200000
CHIP_ds_FST_dummy <- rbind(CHIP_ds_FST_maxes, CHIP_ds_FST)


CHIP_ds_FST_1KWIN <- winScan(x = CHIP_ds_FST_dummy, groups  = "LG", position = "BP", win_size = 1000, win_step = 1000,  values = "FST", funs = "mean")
CHIP_ds_FST_1KWIN <- CHIP_ds_FST_1KWIN[CHIP_ds_FST_1KWIN$FST_n > 0,]
CHIP_ds_FST_1KWIN$BP <- CHIP_ds_FST_1KWIN$win_mid

CHIP_ds_PC_Overlap <- inner_join(CHIP_ds_FST_1KWIN, pcangsd_sel_1kwin, by = c("LG", "win_mid"))
sum(CHIP_ds_PC_Overlap$FST_n)
2857/(3200-39)
cor.test(CHIP_ds_PC_Overlap$PC_S_mean, CHIP_ds_PC_Overlap$FST_mean)


CHIP_ds_FST_overlap<- inner_join(CHIP_ds_FST_1KWIN, paledark_fst)
cor.test(CHIP_ds_FST_overlap$FST, CHIP_ds_FST_overlap$FST_mean)

ggplot() + geom_point(data = CHIP_ds_PC_Overlap, aes(x = PC_S_mean, y = FST_mean)) +
  geom_smooth(data = CHIP_ds_PC_Overlap, aes(x = PC_S_mean, y = FST_mean), method = "glm") + theme_classic()
#CHIPFST_PCANGSD_COR_DS_1Kwindo
ggplot() + geom_point(data = CHIP_ds_FST_overlap, aes(x = FST, y = FST_mean)) +
  geom_smooth(data = CHIP_ds_FST_overlap, aes(x = FST, y = FST_mean), method = "glm") + theme_classic()
#CHIPFST_FST_COR_DS_1Kwindo
