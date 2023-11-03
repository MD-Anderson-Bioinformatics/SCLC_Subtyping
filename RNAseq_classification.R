
#############################################################
### This code allows to subtype SCLC based on RNAseq data ###
#############################################################

#load in packages
library(tidyverse)
library(caret)
library(ComplexHeatmap)

#load in data
data <- #load in your data

# The code assumes that the data was organized that all samples are in columns and gene names are per row in teh first column. If the format differs, the code needs to be changed
data.t <- as.data.frame(t(data[, -1]))
colnames(data.t) <- data$Hugo_Symbol

#load models
models.A <- readRDS("Models_A.Rds")
models.N <- readRDS("Models_N.Rds")
models.P <- readRDS("Models_P.Rds")
models.I <- readRDS("Models_I.Rds")



# the required gene ratios are directly taken from the models
models.filt <- c(models.A, models.N, models.P, models.I)

ratios.models <- lapply(models.filt, "[[", "coefnames")
ratios.models <- unlist(ratios.models)
ratios.models <- ratios.models[!duplicated(ratios.models)]
ratios.models <- gsub("`", "", ratios.models)

ratios.forModels <- cbind.data.frame("Ratio" = ratios.models) %>% separate(Ratio, into = c("1", "2"), sep = "/", remove = FALSE)

#create gene ratios
ratios <- c()
for(i in 1:nrow(ratios.forModels)){
  if(ratios.forModels[i, 2] %in% colnames(data.t) & ratios.forModels[i, 3] %in% colnames(data.t)) {
    j <- ratios.forModels[i, 1]
    ratios[[j]] <- data.t[, unlist(ratios.forModels[i, 2])] - data.t[, unlist(ratios.forModels[i, 3])]
  } else if (ratios.forModels[i, 2] %in% colnames(data.t)){
    k <- ratios.forModels[i, 3]
    data[, k] <- 0
    j <- ratios.forModels[i, 1]
    ratios[[j]] <- data.t[, unlist(ratios.forModels[i, 2])] - data.t[, unlist(ratios.forModels[i, 3])]
  } else {
    k <- ratios.forModels[i, 2]
    data[, k] <- 0
    j <- ratios.forModels[i, 1]
    ratios[[j]] <- data.t[, unlist(ratios.forModels[i, 2])] - data.t[, unlist(ratios.forModels[i, 3])]
  }
}

#combine data to dataframe
ratios.data <- do.call(cbind.data.frame, ratios)
colnames(ratios.data) <- names(ratios)

is.na(ratios.data) <- sapply(ratios.data, is.infinite)
ratios.data[is.na(ratios.data)] <- NA


#Predict Subtypes
predict.A <- predict(models.A, ratios.data, na.action = na.pass)
predict.N <- predict(models.N, ratios.data, na.action = na.pass)
predict.P <- predict(models.P, ratios.data, na.action = na.pass)
predict.I <- predict(models.I, ratios.data, na.action = na.pass)


predict.A.df <- t(do.call(cbind.data.frame, predict.A))
predict.N.df <- t(do.call(cbind.data.frame, predict.N))
predict.P.df <- t(do.call(cbind.data.frame, predict.P))
predict.I.df <- t(do.call(cbind.data.frame, predict.I))

predict.A.df[predict.A.df == "rest"] <- NA
predict.N.df[predict.N.df == "rest"] <- NA
predict.P.df[predict.P.df == "rest"] <- NA
predict.I.df[predict.I.df == "rest"] <- NA

predict.all <- rbind.data.frame(predict.A.df, predict.N.df, predict.P.df, predict.I.df)

#get frequency of predictions
freq.A <- colSums(predict.A.df == "A", na.rm = TRUE) / nrow(predict.A.df) * 100
freq.N <- colSums(predict.N.df == "N", na.rm = TRUE) / nrow(predict.N.df) * 100
freq.P <- colSums(predict.P.df == "P", na.rm = TRUE) / nrow(predict.P.df) * 100
freq.I <- colSums(predict.I.df == "I", na.rm = TRUE) / nrow(predict.I.df) * 100

predict.subtypes <- rbind.data.frame(t(do.call(cbind.data.frame, predict.A)),
                                     t(do.call(cbind.data.frame, predict.N)),
                                     t(do.call(cbind.data.frame, predict.P)),
                                     t(do.call(cbind.data.frame, predict.I)))

predict.subtypes[predict.subtypes == "rest"] <- NA


freq <- sapply(predict.subtypes, function(x) table(factor(x, levels=unique(unlist(predict.subtypes)), ordered=TRUE)))


freq <- t((t(freq) / rowSums(t(freq))) *100)

#create table with classification if consensus for one subtype > 50%
consensus <- as.data.frame(t(freq))
consensus$consensus <- apply(consensus,1,function(x) names(consensus)[which(x>50)])
consensus$consensus <- ifelse(consensus$consensus == "character(0)", "equivocal", consensus$consensus)
consensus$consensus <- unlist(consensus$consensus)
consensus$ID <- colnames(data)[-1]

##############################################################################################
### The consensus table now has all the information including the final subtype assignment ###
##############################################################################################


#create Heatmap with results
col.models <- circlize::colorRamp2(c(10, 50, 70), c("#375696", "yellow3", "seagreen3"))
col.subtypes <- c("A" = "darkblue", "N" = "forestgreen", "P" = "firebrick", "I" = "darkgoldenrod1", "equivocal" = "gray", "A/N" = "cadetblue4","Unknown" = "white")

data.exp <- data

data.expression <- as.data.frame(data.exp[data.exp$Hugo_Symbol %in% c("ASCL1", "NEUROD1", "POU2F3"), ])
rownames(data.expression) <- data.expression$Hugo_Symbol
data.expression <- data.expression[, -1]

data.models <- as.data.frame(t(consensus[, 1:4]))
colnames(data.models) <- consensus$ID



top.anno <- HeatmapAnnotation("Subtype_RNA" = consensus$consensus,
                              col = list("Subtype_RNA" = col.subtypes), annotation_height = c(unit(0.5, "cm")))


hm.a <- Heatmap(data.models, top_annotation = top.anno, cluster_rows = TRUE, cluster_columns = TRUE, name = "consensus",
                show_row_names = TRUE, col = col.models, column_names_gp = gpar(fontsize = 4), column_split = consensus$consensus,
                height = unit(2.4, "cm"), column_dend_height = unit(1, "cm"), row_dend_width = unit(1, "cm"), row_title = "Consensus")

hm.b <-  Heatmap(data.expression, show_row_names = TRUE, column_names_gp = gpar(fontsize = 4), column_split = consensus$consensus, name = "Log2 + 1",
                 height = unit(2.4, "cm"), column_dend_height = unit(1, "cm"), row_dend_width = unit(1, "cm"), row_title = "RNAseq")



pdf("Heatmap.pdf", width = 13, height = 4)
draw(hm.a %v% hm.b)
dev.off()


### In addition this code generates a larger heatmap as in Figure 1B

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
data.hla <- data.exp[data.exp$Hugo_Symbol %in% hla.vec, ] %>% as.data.frame()
data.tis <- data.exp[data.exp$Hugo_Symbol %in% tis.vec, ] %>% as.data.frame()
data.ne <- data.exp[data.exp$Hugo_Symbol %in% ne, ] %>% as.data.frame()
data.nne <- data.exp[data.exp$Hugo_Symbol %in% non.ne, ] %>% as.data.frame()

rownames(data.hla) <- data.hla$Hugo_Symbol
data.hla <- t(scale(t(data.hla[, -1])))



rownames(data.tis) <- data.tis$Hugo_Symbol
data.tis <- t(scale(t(data.tis[, -1])))

rownames(data.ne) <- data.ne$Hugo_Symbol
data.ne <- t(scale(t(data.ne[, -1])))

rownames(data.nne) <- data.nne$Hugo_Symbol
data.nne <- t(scale(t(data.nne[, -1])))

data.expression <- as.data.frame(data.exp[data.exp$Hugo_Symbol %in% c("ASCL1", "NEUROD1", "POU2F3"), ])
rownames(data.expression) <- data.expression$Hugo_Symbol
data.expression <- data.expression[, -1]

data.models <- as.data.frame(t(consensus[, 1:4]))
colnames(data.models) <- consensus$ID


#create heatmap
col.models <- circlize::colorRamp2(c(10, 50, 70), c("#375696", "yellow3", "seagreen3"))
col.subtypes <- c("A" = "darkblue", "N" = "forestgreen", "P" = "firebrick", "I" = "darkgoldenrod1", "equivocal" = "gray", "A/N" = "cadetblue4","Unknown" = "white")

top.anno <- HeatmapAnnotation("Subtype_RNA" = consensus$consensus,
                              col = list("Subtype_RNA" = col.subtypes),
                              annotation_height = c(unit(c(0.5, 0.5), "cm")))


hm.cons <- Heatmap(data.models, top_annotation = top.anno, cluster_rows = TRUE, cluster_columns = TRUE, name = "consensus",
                   show_row_names = TRUE, col = col.models, column_names_gp = gpar(fontsize = 4), column_split = consensus$consensus,
                   height = unit(2.4, "cm"), column_dend_height = unit(1, "cm"), row_dend_width = unit(1, "cm"), row_title = "Consensus")

hm.expr <-  Heatmap(data.expression, show_row_names = TRUE, column_names_gp = gpar(fontsize = 4), column_split = consensus$consensus, name = "Log2 + 1",
                    height = unit(2.4, "cm"), column_dend_height = unit(1, "cm"), row_dend_width = unit(1, "cm"), row_title = "RNAseq")

hm.ne <-  Heatmap(data.ne, show_row_names = FALSE, show_column_names = FALSE, column_split = consensus$consensus, name = "Z-Score",
                  height = unit(3, "cm"), column_dend_height = unit(1, "cm"), row_dend_width = unit(1, "cm"), row_title = "Neuroendocrine")

hm.nne <-  Heatmap(data.nne, show_row_names = FALSE, show_column_names = FALSE, column_split = consensus$consensus, name = "Z-Score",
                   height = unit(3, "cm"), column_dend_height = unit(1, "cm"), row_dend_width = unit(1, "cm"), row_title = "Non-neuroendocrine")

hm.tis <-  Heatmap(data.tis, show_row_names = FALSE, show_column_names = FALSE, column_split = consensus$consensus, name = "Z-Score",
                   height = unit(3, "cm"), column_dend_height = unit(1, "cm"), row_dend_width = unit(1, "cm"), row_title = "TIS")

hm.hla <-  Heatmap(data.hla, show_row_names = FALSE, column_names_gp = gpar(fontsize = 4), column_split = consensus$consensus, name = "Z-Score",
                   height = unit(3, "cm"), column_dend_height = unit(1, "cm"), row_dend_width = unit(1, "cm"), row_title = "HLA")



pdf("Heatmap_Full.pdf", width = 18, height = 12)
draw(hm.cons %v% hm.expr %v% hm.ne %v% hm.nne %v% hm.tis %v% hm.hla)
dev.off()
