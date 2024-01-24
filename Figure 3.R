##################
#### Figure 3 ####
##################

#load packages
library(tidyverse)
library(ComplexHeatmap)
library(circlize)


#### Figure 3B ####

#load in data and models
models.A <- readRDS("Models_50/models_A.Rds")
models.N <- readRDS("Models_50/models_N.Rds")
models.P <- readRDS("Models_50/models_P.Rds")
models.I <- readRDS("Models_50/models_I.Rds")


data.all <- readRDS("Combined_data_Bins100_Filt10_filt90_annotated.Rds")

predict.A <- predict(models.A, data.all, na.action = na.pass)
predict.N <- predict(models.N, data.all, na.action = na.pass)
predict.P <- predict(models.P, data.all, na.action = na.pass)
predict.I <- predict(models.I, data.all, na.action = na.pass)

predict.subtypes <- rbind.data.frame(t(do.call(cbind.data.frame, predict.A)),
                                     t(do.call(cbind.data.frame, predict.N)),
                                     t(do.call(cbind.data.frame, predict.P)),
                                     t(do.call(cbind.data.frame, predict.I)))

#save to disk
predict.subtypes[predict.subtypes == "rest"] <- NA


freq <- sapply(predict.subtypes, function(x) table(factor(x, levels=unique(unlist(predict.subtypes)), ordered=TRUE)))

freq <- (t(freq) / rowSums(t(freq))) *100

consensus <- cbind.data.frame(freq)
consensus$Consensus <- ifelse(consensus$A >= 50, "A", ifelse(consensus$N >= 50, "N", ifelse(consensus$P >= 50, "P", ifelse(consensus$I >= 50, "I", "equivocal"))))
consensus$ID <- data.all$Sample
consensus$Subtype_RNA <- data.all$Subtype
consensus <- read.csv("consensus.csv")
consensus$Cohort <- ifelse(consensus$ID %in% data.train$Sample, "Training",
                           ifelse(consensus$ID %in% data.test$Sample, "Testing", "Independent"))


#create Heatmap
col.models <- circlize::colorRamp2(c(10, 50, 70), c("#375696", "yellow3", "seagreen3"))
col.subtypes <- c("A" = "darkblue", "N" = "forestgreen", "P" = "firebrick", "I" = "darkgoldenrod1", "equivocal" = "gray", "equivocal" = "gray", "A/N" = "cadetblue4","rest" = "gray90")
col.pred <- c("A" = "darkblue", "N" = "forestgreen", "P" = "firebrick", "I" = "darkgoldenrod1", "equivocal" = "gray", "A/N" = "cadetblue4","rest" = "gray90")

consensus$Subtype_RNA <- sub("SCLC-", "", consensus$Subtype_RNA)
consensus <- consensus[!c(consensus$Cohort == "Independent" & !is.na(consensus$Subtype_RNA)), ]

anno <- HeatmapAnnotation("Subtype_DMC" = consensus$Consensus, "Subtype_RNA" = consensus$Subtype_RNA,
                          col = list("Subtype_DMC" = col.pred, "Subtype_RNA" = col.pred), annotation_height = unit(c(0.5, 0.5), "cm"))


data.hm <- t(consensus[, 1:4])
colnames(data.hm) <- consensus$ID

hm <- Heatmap(data.hm, cluster_rows = FALSE, cluster_columns = TRUE, name = "consensus", show_column_names = TRUE, column_names_gp = gpar(fontsize = 4),
              show_row_names = TRUE, col = col.models, top_annotation = anno, column_split = consensus$Cohort,
              height = unit(2.4, "cm"), column_dend_height = unit(2, "cm"), row_dend_width = unit(2, "cm"), row_title = "Consensus")

pdf(paste0(Sys.Date(), "_Heatmap_consensus_FFPE.pdf"), width = 12, height = 5)
draw(hm)
dev.off()


######



#### Figure 3C ####

#load in data. The DNA methylation file as well as the ctDNA fraction data is needed
data <- read.csv("data.csv")

#create plot only for the selected
selection <- c("chr12:27974490", "chr1:7236563", "chr17:29139387", "chr19:128737209", "chr2:10401557", "chr21:34669078", "chr21:45590104")

#select sites manually
data.filt <- data[, colnames(data) %in% sites]

calc.frac <- rowMeans(data.filt[, colnames(data.filt) %in% selection], na.rm = TRUE)

data.plot <- data[, 1:2] %>% cbind.data.frame("calc.frac" = calc.frac)


#correlate the calculated fraction with actual fraction
plot <- ggplot(data.plot, aes(x = Fraction, y = calc.frac)) + geom_point() + geom_smooth(method = "lm") + theme_bw() + ggpubr::stat_cor() + 
  xlab("Calculated Fraction [%]") + ylab("ULP-WGS Fraction [%]") + xlim(0, 100) + ylim(0, 100) + 
  geom_abline(intercept =0 , slope = 1, lty = 3, size = 1, col = "gray80") +   geom_text(aes(label = Sample), hjust = 0, nudge_x = 0.05, size = 1.5)


pdf(paste0(Sys.Date(), "_Correlation_Plots_ULPWGS_Methylation_meanModel_selectedSites.pdf"), width = 7, height = 7)
plot
dev.off()


######



#### Figure 3E ####

#load in 100b bp bin data from tissue and cfDNA
tissue.data <- readRDS("100K_Bin_Tissue.Rds")
cfDNA.data <- readRDS("100k_Bin_cfDNA.Rds")

bins <- readRDS("100kBins.Rds")
bins <- bins[bins$seqnames %in% c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22"), ]


#get overlap
tissue.data <- tissue.data %>% select(one_of(colnames(cfDNA.data)))
cfDNA.data <- cfDNA.data %>% select(one_of(colnames(tissue.data)))


#combine
tissue.melt <- reshape2::melt(tissue.data, id.vars = c("bin", "seqnames"))
cfDNA.melt <- reshape2::melt(cfDNA.data, id.vars = c("bin", "seqnames"))

melt.comb <- tissue.melt %>% left_join(cfDNA.melt %>% select(bin, variable, "cfDNA" = value))
melt.comb$diff <- melt.comb$value - melt.comb$cfDNA

#recast
comb.final <- melt.comb %>% select(bin, variable, diff) %>% pivot_wider(names_from = variable, values_from = diff)

data.hm <- t(comb.final[, -1])
data.hm[!is.numeric(data.hm)] <- NA

#create means per bin
mean.df <- cbind.data.frame("bin" = bins$bin, "mean" = colMeans(data.hm, na.rm = TRUE))

#create heatmap
bins$seq_fact <- factor(bins$seqnames, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22"))

col_fun <- circlize::colorRamp2(c(-100, 0, 100), c("darkolivegreen4", "gray90", "firebrick4"))

top.anno <- HeatmapAnnotation("mean" = anno_lines(mean.df$mean, smooth = FALSE, add_points = FALSE, gp = gpar(lwd = 0.2, col = "gray60")), annotation_height = unit(c(15), "mm"))

legend.param <- list(legend_direction = "horizontal", legend_width = unit(3, "cm"))

hm <- Heatmap(data.hm, show_column_names = FALSE, column_split = bins$seq_fact, cluster_columns = FALSE, cluster_rows = FALSE, na_col = "white", cluster_column_slices = FALSE, column_order = colnames(data.hm), border = TRUE,
              top_annotation = top.anno,  heatmap_legend_param = legend.param, show_row_names = FALSE)



pdf(paste0(Sys.Date(), "_Heatmap_differencesMethylation.pdf"), width = 14, height = 4)
draw(hm, heatmap_legend_side="bottom", annotation_legend_side="right")
dev.off()


top.anno <- HeatmapAnnotation("mean" = anno_lines(mean.df$mean, smooth = TRUE, add_points = FALSE, gp = gpar(lwd = 0.8, col = "black")), annotation_height = unit(c(15), "mm"))

legend.param <- list(legend_direction = "horizontal", legend_width = unit(3, "cm"))

hm <- Heatmap(data.hm, show_column_names = FALSE, column_split = bins$seq_fact, cluster_columns = FALSE, cluster_rows = FALSE, na_col = "white", cluster_column_slices = FALSE, column_order = colnames(data.hm), border = TRUE,
              top_annotation = top.anno,  heatmap_legend_param = legend.param, show_row_names = FALSE)



pdf(paste0(Sys.Date(), "_Heatmap_differencesMethylation_smoothed.pdf"), width = 14, height = 4)
draw(hm, heatmap_legend_side="bottom", annotation_legend_side="right")
dev.off()

#create histogram
data.hist <- as.data.frame(data.hm) %>% rownames_to_column("Sample") %>% reshape2::melt()

plot.hist <- ggplot(data.hist, aes(x = value)) + geom_histogram(bins = 51, col = "black", fill = "gray80") + geom_vline(xintercept = 0, lty = "dashed", col = "red") + theme_bw()

pdf(paste0(Sys.Date(), "_histogram_methylation_diff.pdf"), width = 3, height = 5)
plot.hist
dev.off()





######


#### Figure 3F ####

#load in models
models.A <- readRDS("Models_50_cf/models_A.Rds")
models.N <- readRDS("Models_50_cf/models_N.Rds")
models.P <- readRDS("Models_50_cf/models_P.Rds")
models.I <- readRDS("Models_50_cf/models_I.Rds")


predict.A <- predict(models.A, data.cf, na.action = na.pass)
predict.N <- predict(models.N, data.cf, na.action = na.pass)
predict.P <- predict(models.P, data.cf, na.action = na.pass)
predict.I <- predict(models.I, data.cf, na.action = na.pass)

predict.subtypes <- rbind.data.frame(t(do.call(cbind.data.frame, predict.A)),
                                     t(do.call(cbind.data.frame, predict.N)),
                                     t(do.call(cbind.data.frame, predict.P)),
                                     t(do.call(cbind.data.frame, predict.I)))

#save to disk
predict.subtypes[predict.subtypes == "rest"] <- NA


freq <- sapply(predict.subtypes, function(x) table(factor(x, levels=unique(unlist(predict.subtypes)), ordered=TRUE)))

freq <- (t(freq) / rowSums(t(freq))) *100

consensus <- cbind.data.frame(freq)
consensus$Consensus <- ifelse(consensus$A >= 50, "A", ifelse(consensus$N >= 50, "N", ifelse(consensus$P >= 50, "P", ifelse(consensus$I >= 50, "I", "equivocal"))))
consensus$ID <- data.cf$Sample
consensus$ID <- sub("JHSC.", "JHSC-", consensus$ID)
consensus$Subtype <- data.cf$Subtype

#create heatmap
col.models <- circlize::colorRamp2(c(10, 50, 70), c("#375696", "yellow3", "seagreen3"))
col.subtypes <- c("A" = "darkblue", "N" = "forestgreen", "P" = "firebrick", "I" = "darkgoldenrod1", "equivocal" = "gray", "equivocal" = "gray", "A/N" = "cadetblue4","rest" = "gray90")
col.pred <- c("A" = "darkblue", "N" = "forestgreen", "P" = "firebrick", "I" = "darkgoldenrod1", "equivocal" = "gray", "A/N" = "cadetblue4","rest" = "gray90")

consensus$Subtype_RNA <- sub("SCLC.", "", consensus$Subtype)

anno <- HeatmapAnnotation("Subtype_DMC" = consensus$Consensus, "Subtype_RNA" = consensus$Subtype_RNA,
                          col = list("Subtype_DMC" = col.pred, "Subtype_RNA" = col.pred), annotation_height = unit(c(0.5, 0.5), "cm"))


data.hm <- t(consensus[, 1:4])
colnames(data.hm) <- consensus$ID

hm <- Heatmap(data.hm, cluster_rows = FALSE, cluster_columns = TRUE, name = "consensus", show_column_names = TRUE, column_names_gp = gpar(fontsize = 4),
              show_row_names = TRUE, col = col.models, top_annotation = anno, 
              height = unit(2.4, "cm"), column_dend_height = unit(2, "cm"), row_dend_width = unit(2, "cm"), row_title = "Consensus")

pdf("Models_50_cf/Heatmap_consensus_cfDNA.pdf", width = 12, height = 5)
draw(hm)
dev.off()

######




