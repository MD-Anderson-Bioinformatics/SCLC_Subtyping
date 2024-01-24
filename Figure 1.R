##################
#### Figure 1 ####
##################

#load packages
library(tidyverse)
library(ComplexHeatmap)
library(circlize)

### Figure 1B ###

# Load data from RNA expression, classification from thr GRC approach and cohort annotation
data.rna <- read.csv("RNAseq.csv")
data.class <- read.csv("classification.csv")
data.cohort <- read.csv("annotation.csv")


#define colors for heatmap
col.models <- colorRamp2(c(10, 50, 70), c("#375696", "yellow3", "seagreen3"))
col.subtypes <- c("A" = "darkblue", "N" = "forestgreen", "P" = "firebrick", "I" = "darkgoldenrod1", "equivocal" = "gray")
col.cohort <- c("C1" = "black", "C2" = "gray50")


#create data for the TF expressions
data.expression <- as.data.frame(data.rna[data.rna$Hugo_Symbol %in% c("ASCL1", "NEUROD1", "POU2F3"), ])
rownames(data.expression) <- data.expression$Hugo_Symbol
data.expression <- data.expression[, -1]

#create data for the expression of subtypes
data.models <- as.data.frame(t(data.class[, 1:4]))
colnames(data.models) <- data.class$ID


#create gene signatures
hla.vec <- c("HLA-A", "HLA-B", "HLA-C", "B2M", "HLA-DPA1", "HLA-DPB1", "HLA-DMA", "HLA-DMB", "HLA-DQA1",
             "HLA-DQB1", "HLA-DOA", "HLA-DOB", "HLA-DRA", "HLA-DRB1", "HLA-DRB5", "HLA-DQA2", "HLA-DQB2",
             "HLA-E, HLA-F, HLA-G")
tis.vec <- c("STAT1", "CXCL9", "IDO1", "PSMB10", "LAG3", "TIGIT", "CXCR6", "CCL5", "NKG7", "CD27", "HLA-E",
             "CD274", "HLA-DRB1", "HLA-DQA1", "CMKLR1", "PDCD1LG2", "CD276", "CD8A")

ne <- c("TAGLN3", "SH3GL2", "INSM1", "KIF5C", "KIF1A", "TMSB15A", "TMSB15B", "RTN1", "SYT11", "SYT4", "ASCL1",
        "TFF3", "CHGA", "CHGB", "BEX1", "SEZ6", "RUNDC3A", "SYN1", "SYP", "BSN", "MYT1", "CELF3", "CRMP1",
        "SCG3")

non.ne <- c("S100A10", "RAB27B", "TACSTD2", "ANXA1", "PTGES", "ITGB4", "EPHA2", "CCND1", "ABCC3", "IFITM2",
            "IFITM3", "YAP1", "ARHGDIB", "TGFBR2", "AHNAK", "PLAU", "TGFBI", "EMP1", "CAV2", "CAV1")



#split in annotations
data.hla <- data.rna[data.rna$Hugo_Symbol %in% hla.vec, ] %>% as.data.frame()
data.tis <- data.rna[data.rna$Hugo_Symbol %in% tis.vec, ] %>% as.data.frame()
data.ne <- data.rna[data.rna$Hugo_Symbol %in% ne, ] %>% as.data.frame()
data.nne <- data.rna[data.rna$Hugo_Symbol %in% non.ne, ] %>% as.data.frame()

rownames(data.hla) <- data.hla$Hugo_Symbol
data.hla <- t(scale(t(data.hla[, -1])))

rownames(data.tis) <- data.tis$Hugo_Symbol
data.tis <- t(scale(t(data.tis[, -1])))

rownames(data.ne) <- data.ne$Hugo_Symbol
data.ne <- t(scale(t(data.ne[, -1])))

rownames(data.nne) <- data.nne$Hugo_Symbol
data.nne <- t(scale(t(data.nne[, -1])))

#create heatmap
top.anno <- HeatmapAnnotation("GRC (RNA)" = data.class$consensus, "Cohort" = data.cohort$cohort,
                              col = list("GRC (RNA)" = col.subtypes, "Cohort" = col.cohort), annotation_height = c(unit(c(0.5, 0.5), "cm")))


hm.cons <- Heatmap(data.models, top_annotation = top.anno, cluster_rows = TRUE, cluster_columns = TRUE, name = "consensus",
                   show_row_names = TRUE, col = col.models, column_names_gp = gpar(fontsize = 4), column_split = data.class$consensus,
                   height = unit(2.4, "cm"), column_dend_height = unit(1, "cm"), row_dend_width = unit(1, "cm"), row_title = "Consensus")

hm.expr <-  Heatmap(data.expression, show_row_names = TRUE, column_names_gp = gpar(fontsize = 4), column_split = data.class$consensus, name = "Log2 + 1",
                    height = unit(2.4, "cm"), column_dend_height = unit(1, "cm"), row_dend_width = unit(1, "cm"), row_title = "RNAseq")

hm.ne <-  Heatmap(data.ne, show_row_names = FALSE, show_column_names = FALSE, column_split = data.class$consensus, name = "Z-Score",
                  height = unit(3, "cm"), column_dend_height = unit(1, "cm"), row_dend_width = unit(1, "cm"), row_title = "Neuroendocrine")

hm.nne <-  Heatmap(data.nne, show_row_names = FALSE, show_column_names = FALSE, column_split = data.class$consensus, name = "Z-Score",
                   height = unit(3, "cm"), column_dend_height = unit(1, "cm"), row_dend_width = unit(1, "cm"), row_title = "Non-neuroendocrine")

hm.tis <-  Heatmap(data.tis, show_row_names = FALSE, show_column_names = FALSE, column_split = data.class$consensus, name = "Z-Score",
                   height = unit(3, "cm"), column_dend_height = unit(1, "cm"), row_dend_width = unit(1, "cm"), row_title = "TIS")

hm.hla <-  Heatmap(data.hla, show_row_names = FALSE, column_names_gp = gpar(fontsize = 4), column_split = data.class$consensus, name = "Z-Score",
                   height = unit(3, "cm"), column_dend_height = unit(1, "cm"), row_dend_width = unit(1, "cm"), row_title = "HLA")



pdf("Heatmap_Figure1B.pdf", width = 18, height = 12)
draw(hm.cons %v% hm.expr %v% hm.ne %v% hm.nne %v% hm.tis %v% hm.hla)
dev.off()

#############



### Figure 1C ###

#analyse ESTIMATE immune infiltration

#load data 

#estimate can be found here: https://bioinformatics.mdanderson.org/estimate/index.html
estimate <- read.csv("ESTIMATE_table.csv") %>% as.data.frame()
anno <- read.csv("anno.csv")

colnames(estimate)[1] <- "Sample"

#add classification data
estimate <- estimate %>% left_join(anno %>% select(Sample, consensus))

#plot differences

estimate.filt <- estimate[estimate$consensus %in% c("A", "N", "P", "I"), ]
estimate.filt$consensus <- factor(estimate.filt$consensus, levels = c("A", "N", "P", "I"))

comp <- list(c("A", "N"), c("A", "P"), c("A", "I"), c("N", "P"), c("N", "I"), c("P", "I"))

plot <- ggplot(estimate.filt, aes(x = consensus, y = ImmuneScore)) + geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.3) + stat_compare_means(comparisons = comp) + theme_bw()

#save to disk
pdf(paste0(Sys.Date(), "_Figure_ESTIMATE.pdf"), width = 6, height = 8)
plot
dev.off()




### Figure 1D ###


#load in data. The consensus prediction is required. 
data.class <- read.csv("classification.csv")

#change direction
data$coordY <- data$A - data$P
data$coordX <- data$N - data$I


col.pred <- c("A" = "darkblue", "N" = "forestgreen", "P" = "firebrick", "I" = "darkgoldenrod1", "equivocal" = "gray")


#create figure
plot <- ggplot(data, aes(x = coordX, y = coordY, col = Consensus, fill = Consensus)) + geom_abline(intercept = 0, slope = 1, size = 0.2, lty = "dashed") +
  geom_abline(intercept = 0, slope = -1, size = 0.2, lty = "dashed") +
  geom_hline(yintercept = 0) +geom_vline(xintercept = 0) + geom_hex(data = data[data$Consensus != "equivocal", ], bins = 6, fill = NA) + 
  geom_point(size = 2, alpha = 0.5) + theme_bw() + scale_colour_manual(values = col.pred) +
  ggside::geom_ysidedensity()

pdf(paste0(Sys.Date(), "_Figure_coordinates.pdf"), width = 10, height = 10)
ggExtra::ggMarginal(plot, type="histogram", groupFill = TRUE, binwidth = 5)
dev.off()0

