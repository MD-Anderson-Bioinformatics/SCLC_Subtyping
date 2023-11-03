
############################################################################################################
### This code describes the selection of appropriate methylation sites and training of predictive models ###
############################################################################################################

library(GenomicRanges)
library(BSgenome.Hsapiens.NCBI.GRCh38)
library(ComplexHeatmap)
library(GenomicRanges)
library(annotatr)
library(tidyverse)
library(rstatix)
library(psych)
library(pROC)
library(caret)


##################################################################################
#### Create a combined dataset with cov > 10 and at least 90% of data present ####
##################################################################################

##### For Bin 100 #####


#load in data of the two cohorts with data filtered for 10 reads
data.1 <- #load in data of cohort 1
data.2 <- #load in data of cohort 2


#create 100 bins for data.1
data.1$Chr <- sub("chr", "", data.1$Chr)


bins <- as.data.frame(tileGenome(seqinfo(Hsapiens), tilewidth=100, cut.last.tile.in.chrom=TRUE))
bins$bin <- 1:nrow(bins)
bins <- as(bins, "GRanges")
bins.df <- as.data.frame(bins)



data.gr <- makeGRangesFromDataFrame(data.1, seqnames.field = "Chr", start.field = "Pos", end.field = "Pos", keep.extra.columns = TRUE)
data.gr.merge <- as.data.frame(mergeByOverlaps(bins, data.gr))
data.gr.merge <- as.data.frame(data.gr.merge)
data.gr.merge <- data.gr.merge[, c(7, 96:178)]

data.bin.mean <- data.gr.merge %>% 
  group_by(bin) %>%
  summarise_all(mean, na.rm = TRUE)

bins.df <- bins.df %>% select(bin, seqnames, start, end) %>% mutate(Position = paste0("chr", seqnames, ":", start, "-", end)) %>% select(bin, Position) %>% filter(bin %in% data.gr.merge$bin)
data.c1 <- cbind.data.frame(bins.df, data.bin.mean[, -1])


#create 100 bins for data.2
data.2$Chr <- sub("chr", "", data.2$Chr)
data.gr <- makeGRangesFromDataFrame(data.2, seqnames.field = "Chr", start.field = "Pos", end.field = "Pos", keep.extra.columns = TRUE)
data.gr.merge <- as.data.frame(mergeByOverlaps(bins, data.gr))

data.gr.merge <- as.data.frame(data.gr.merge)
data.gr.merge <- data.gr.merge[, c(7, 60:106)]

data.bin.mean <- data.gr.merge %>% 
  group_by(bin) %>%
  summarise_all(mean, na.rm = TRUE)

bins.df <- as.data.frame(bins) %>% select(bin, seqnames, start, end) %>% mutate(Position = paste0("chr", seqnames, ":", start, "-", end)) %>% select(bin, Position) %>% filter(bin %in% data.gr.merge$bin)

data.c2 <- cbind.data.frame(bins.df, data.bin.mean[, -1])


#repeat for cell lines
data.cl <- #load in cell line data

data.gr <- makeGRangesFromDataFrame(data.cl, seqnames.field = "Chr", start.field = "Pos", end.field = "Pos", keep.extra.columns = TRUE)
data.gr.merge <- as.data.frame(mergeByOverlaps(bins, data.gr))

data.gr.merge <- as.data.frame(data.gr.merge)
data.gr.merge <- data.gr.merge[, c(7, 72:130)]

data.bin.mean <- data.gr.merge %>% 
  group_by(bin) %>%
  summarise_all(mean, na.rm = TRUE)

bins.df <- as.data.frame(bins) %>% select(bin, seqnames, start, end) %>% mutate(Position = paste0("chr", seqnames, ":", start, "-", end)) %>% select(bin, Position) %>% filter(bin %in% data.gr.merge$bin)

data.cl <- cbind.data.frame(bins.df, data.bin.mean[, -1])



#combine sets
data.combined <- data.1 %>% inner_join(data.2, by = c("bin", "Position"))


#filter oust sites with at least 90% present data
nas <- rowSums(is.na(data.combined)) / (ncol(data.combined) - 2)
data.combined.filt <- data.combined[nas <= 0.1, ]

#remove X and Y chromosome for analysis and transpose data to allow annotation
data.combined.filt <- data.combined.filt[!c(grepl("chrX", data.combined.filt$Position) | grepl("chrY", data.combined.filt$Position)), ]

data.combined.filt.t <- t(data.combined.filt[, -c(1:2)])
colnames(data.combined.filt.t) <- as.character(data.combined.filt$bin)
data.combined.filt.t <- data.combined.filt.t %>% as.data.frame() %>% rownames_to_column("Sample")


#load annotation and annotate data
anno <- #load data with sample annotation data
anno.filt <- anno %>% filter(Histology == "SCLC") %>% select(Sample = ID, Subtype = Subtype)

anno.filt <- cbind.data.frame("Sample" = data.combined.filt.t$Sample) %>% left_join(anno.filt)

data.combined.filt.t <- cbind.data.frame(anno.filt,data.combined.filt.t[, -1])


#########################################
#### create training and testing set ####
#########################################

data.filt <- data.combined.filt.t[data.combined.filt.t$Subtype %in% c("SCLC-A", "SCLC-N", "SCLC-P", "SCLC-I"), ]

index <- createDataPartition(data.filt$Subtype, p = .7, 
                             list = FALSE, 
                             times = 1)

data.train <- data.filt[index, ]
data.test <- data.filt[-index, ]

###################################################
#### run pROC to define best methylation sites ####
###################################################

data.mean <- reshape2::melt(data.train, id.vars = c("Sample", "Subtype")) %>% group_by(Subtype, variable) %>% 
  summarise_at(vars(value), mean, na.rm = TRUE) %>% pivot_wider(names_from = Subtype, values_from = value)

#calculate differences
data.mean$`SCLC-A diff` <- data.mean$`SCLC-A` - rowMeans(data.mean[, 3:5], na.rm = TRUE)
data.mean$`SCLC-N diff` <- data.mean$`SCLC-N` - rowMeans(data.mean[, c(2:3,5)], na.rm = TRUE)
data.mean$`SCLC-P diff` <- data.mean$`SCLC-P` - rowMeans(data.mean[, c(2:4)], na.rm = TRUE)
data.mean$`SCLC-I diff` <- data.mean$`SCLC-I` - rowMeans(data.mean[, c(2,4,5)], na.rm = TRUE)

data.mean <- data.combined %>% mutate(bin = as.character(bin), "Chr_Range" = Position) %>% select(bin, Chr_Range) %>% inner_join(data.mean, by = c("bin" = "variable"))


##### analyse ROC for SCLC-A ####
data.roc <- data.train


#calculate for SCLC-A
data.roc$Subtype <- ifelse(data.roc$Subtype == "SCLC-A", "SCLC-A", "rest")

#prepare ROC data
data.roc <- pkgcond::suppress_messages(lapply(data.roc[,3:ncol(data.roc)], function(row){tryCatch(roc(data.roc$Subtype, row),
                                                                                                  error = function(e) NA)}))

data.roc <- data.roc[!is.na(data.roc)]

#get Youden's J
data.roc.J <- lapply(data.roc, function(row){coords(row, x="best", input="threshold", best.method="youden")})

#drop lines if there are multiple thresholds
data.roc.J <- data.roc.J[sapply(data.roc.J, nrow) == 1]

#prepare df
data.roc.J <- as.data.frame(t(sapply(data.roc.J, '[')))

#get AUC
data.roc.auc <- as.data.frame(sapply(data.roc, '[[', "auc"))

#filter overlapping names and combine both data.frames
data.roc.auc <- data.roc.auc[rownames(data.roc.auc) %in% rownames(data.roc.J), ]

data.roc.df <- cbind(data.roc.auc, data.roc.J)
colnames(data.roc.df)[1] <- "AUC"


#order by AUC
data.roc.df <- data.roc.df %>% rownames_to_column(var = "bin")

res.A <- data.roc.df %>% inner_join(data.mean) %>% arrange(desc(AUC))



##### analyse ROC for SCLC-N ####
data.roc <- data.train


#calculate for SCLC-A
data.roc$Subtype <- ifelse(data.roc$Subtype == "SCLC-N", "SCLC-N", "rest")

#prepare ROC data
data.roc <- pkgcond::suppress_messages(lapply(data.roc[,3:ncol(data.roc)], function(row){tryCatch(roc(data.roc$Subtype, row),
                                                                                                  error = function(e) NA)}))

data.roc <- data.roc[!is.na(data.roc)]

#get Youden's J
data.roc.J <- lapply(data.roc, function(row){coords(row, x="best", input="threshold", best.method="youden")})

#drop lines if there are multiple thresholds
data.roc.J <- data.roc.J[sapply(data.roc.J, nrow) == 1]

#prepare df
data.roc.J <- as.data.frame(t(sapply(data.roc.J, '[')))

#get AUC
data.roc.auc <- as.data.frame(sapply(data.roc, '[[', "auc"))

#filter overlapping names and combine both data.frames
data.roc.auc <- data.roc.auc[rownames(data.roc.auc) %in% rownames(data.roc.J), ]

data.roc.df <- cbind(data.roc.auc, data.roc.J)
colnames(data.roc.df)[1] <- "AUC"


#order by AUC
data.roc.df <- data.roc.df %>% rownames_to_column(var = "bin")

res.N <- data.roc.df %>% inner_join(data.mean) %>% arrange(desc(AUC))


##### analyse ROC for SCLC-P ####
data.roc <- data.train


#calculate for SCLC-A
data.roc$Subtype <- ifelse(data.roc$Subtype == "SCLC-P", "SCLC-P", "rest")

#prepare ROC data
data.roc <- pkgcond::suppress_messages(lapply(data.roc[,3:ncol(data.roc)], function(row){tryCatch(roc(data.roc$Subtype, row),
                                                                                                  error = function(e) NA)}))

data.roc <- data.roc[!is.na(data.roc)]

#get Youden's J
data.roc.J <- lapply(data.roc, function(row){coords(row, x="best", input="threshold", best.method="youden")})

#drop lines if there are multiple thresholds
data.roc.J <- data.roc.J[sapply(data.roc.J, nrow) == 1]

#prepare df
data.roc.J <- as.data.frame(t(sapply(data.roc.J, '[')))

#get AUC
data.roc.auc <- as.data.frame(sapply(data.roc, '[[', "auc"))

#filter overlapping names and combine both data.frames
data.roc.auc <- data.roc.auc[rownames(data.roc.auc) %in% rownames(data.roc.J), ]

data.roc.df <- cbind(data.roc.auc, data.roc.J)
colnames(data.roc.df)[1] <- "AUC"


#order by AUC
data.roc.df <- data.roc.df %>% rownames_to_column(var = "bin")

res.P <- data.roc.df %>% inner_join(data.mean) %>% arrange(desc(AUC))


##### analyse ROC for SCLC-I ####
data.roc <- data.train


#calculate for SCLC-A
data.roc$Subtype <- ifelse(data.roc$Subtype == "SCLC-I", "SCLC-I", "rest")

#prepare ROC data
data.roc <- pkgcond::suppress_messages(lapply(data.roc[,3:ncol(data.roc)], function(row){tryCatch(roc(data.roc$Subtype, row),
                                                                                                  error = function(e) NA)}))

data.roc <- data.roc[!is.na(data.roc)]

#get Youden's J
data.roc.J <- lapply(data.roc, function(row){coords(row, x="best", input="threshold", best.method="youden")})

#drop lines if there are multiple thresholds
data.roc.J <- data.roc.J[sapply(data.roc.J, nrow) == 1]

#prepare df
data.roc.J <- as.data.frame(t(sapply(data.roc.J, '[')))

#get AUC
data.roc.auc <- as.data.frame(sapply(data.roc, '[[', "auc"))

#filter overlapping names and combine both data.frames
data.roc.auc <- data.roc.auc[rownames(data.roc.auc) %in% rownames(data.roc.J), ]

data.roc.df <- cbind(data.roc.auc, data.roc.J)
colnames(data.roc.df)[1] <- "AUC"


#order by AUC
data.roc.df <- data.roc.df %>% rownames_to_column(var = "bin")

res.I <- data.roc.df %>% inner_join(data.mean) %>% arrange(desc(AUC))




#############################################
#### annotate sites with gene annotation ####
#############################################

#create annotation matrix
anno.gr <- res.A %>% separate(Chr_Range, into = c("Chr", "Pos"), sep = ":") %>% separate(Pos, into = c("Start", "End"), sep = "-") %>% 
                              dplyr::select(Chr, Start, End) %>%
                              makeGRangesFromDataFrame(keep.extra.columns=FALSE, seqnames.field= "Chr", 
                                    start.field="Start", end.field="End")

annots <- c("hg38_genes_promoters","hg38_genes_exons", "hg38_genes_introns", "hg38_genes_1to5kb", 
            "hg38_genes_5UTRs", "hg38_genes_intergenic", "hg38_genes_3UTRs", "hg38_genes_firstexons", "hg38_genes_intronexonboundaries", 
            "hg38_genes_exonintronboundaries")
annotations <- build_annotations(genome = 'hg38', annotations = annots)

dm_annotated <- annotate_regions(regions = anno.gr, annotations = annotations, ignore.strand = TRUE, quiet = FALSE)
df_dm_annotated <-  data.frame(dm_annotated)
df_dm_annotated$Position <- paste0(df_dm_annotated$seqnames, ":", df_dm_annotated$start, "-", df_dm_annotated$end)
df_dm_annotated <- df_dm_annotated %>% left_join(dplyr::select(res.A, bin, Chr_Range), by = c("Position" = "Chr_Range"))


#create annotations for CpG
annots.cpg <- c("hg38_cpgs")
annotations.cpg <- build_annotations(genome = 'hg38', annotations = annots.cpg)

df <- as.data.frame(annotations.cpg)


#combine together
df_dm_annotated_combined <- df_dm_annotated %>% group_by(bin) %>% summarise_all(funs(paste(.,collapse = ', ')))
df_anno <- df_dm_annotated_combined %>% dplyr::select(bin, annot.symbol, annot.id)

#annotate
res.A <- res.A %>% left_join(df_anno)
res.N <- res.N %>% left_join(df_anno)
res.P <- res.P %>% left_join(df_anno)
res.I <- res.I %>% left_join(df_anno)

#save to disk
data.table::fwrite(res.A, "ROC_FFPE_A_Train_filt90_annotated.csv")
data.table::fwrite(res.N, "ROC_FFPE_N_Train_filt90_annotated.csv")
data.table::fwrite(res.P, "ROC_FFPE_P_Train_filt90_annotated.csv")
data.table::fwrite(res.I, "ROC_FFPE_I_Train_filt90_annotated.csv")


###################################################################
#### run pROC to define best methylation sites from cell lines ####
###################################################################

#load in cell line data
data.cl <- #load in cell line data
nas <- rowSums(is.na(data.cl)) / (ncol(data.cl) - 2)
data.cl.filt <- data.cl[nas <= 0.3, ]

#transpose
data.cl.filt.t <- t(data.cl.filt[, -c(1:2)])
colnames(data.cl.filt.t) <- data.cl.filt$bin

anno <- #load in cell line annotation
anno.filt <- anno %>% select(Sample = `Cell line`, Subtype = Classification)
anno.filt <- cbind.data.frame("Sample" = colnames(data.cl)[-c(1:2)]) %>% left_join(anno.filt)

data.cl.filt.t <- cbind.data.frame(anno.filt, data.cl.filt.t)

data.mean <- reshape2::melt(data.cl.filt.t, id.vars = c("Sample", "Subtype")) %>% group_by(Subtype, variable) %>% 
  summarise_at(vars(value), mean, na.rm = TRUE) %>% pivot_wider(names_from = Subtype, values_from = value)

#calculate differences
data.mean$`SCLC-A diff` <- data.mean$`SCLC-A` - rowMeans(data.mean[, 3:5], na.rm = TRUE)
data.mean$`SCLC-N diff` <- data.mean$`SCLC-N` - rowMeans(data.mean[, c(2:3,5)], na.rm = TRUE)
data.mean$`SCLC-P diff` <- data.mean$`SCLC-P` - rowMeans(data.mean[, c(2:4)], na.rm = TRUE)
data.mean$`SCLC-I diff` <- data.mean$`SCLC-I` - rowMeans(data.mean[, c(2,4,5)], na.rm = TRUE)

data.mean <- data.cl %>% mutate(bin = as.character(bin), "Chr_Range" = Position) %>% select(bin, Chr_Range) %>% inner_join(data.mean, by = c("bin" = "variable"))


##### analyse ROC for SCLC-A ####
data.roc <- data.cl.filt.t


#calculate for SCLC-A
data.roc$Subtype <- ifelse(data.roc$Subtype == "SCLC-A", "SCLC-A", "rest")

#prepare ROC data
data.roc <- pkgcond::suppress_messages(lapply(data.roc[,3:ncol(data.roc)], function(row){tryCatch(roc(data.roc$Subtype, row),
                                                                                                  error = function(e) NA)}))

data.roc <- data.roc[!is.na(data.roc)]

#get Youden's J
data.roc.J <- lapply(data.roc, function(row){coords(row, x="best", input="threshold", best.method="youden")})

#drop lines if there are multiple thresholds
data.roc.J <- data.roc.J[sapply(data.roc.J, nrow) == 1]

#prepare df
data.roc.J <- as.data.frame(t(sapply(data.roc.J, '[')))

#get AUC
data.roc.auc <- as.data.frame(sapply(data.roc, '[[', "auc"))

#filter overlapping names and combine both data.frames
data.roc.auc <- data.roc.auc[rownames(data.roc.auc) %in% rownames(data.roc.J), ]

data.roc.df <- cbind(data.roc.auc, data.roc.J)
colnames(data.roc.df)[1] <- "AUC"


#order by AUC
data.roc.df <- data.roc.df %>% rownames_to_column(var = "bin")

res.A.cl <- data.roc.df %>% inner_join(data.mean) %>% arrange(desc(AUC))





##### analyse ROC for SCLC-N ####
data.roc <- data.cl.filt.t


#calculate for SCLC-A
data.roc$Subtype <- ifelse(data.roc$Subtype == "SCLC-N", "SCLC-N", "rest")

#prepare ROC data
data.roc <- pkgcond::suppress_messages(lapply(data.roc[,3:ncol(data.roc)], function(row){tryCatch(roc(data.roc$Subtype, row),
                                                                                                  error = function(e) NA)}))

data.roc <- data.roc[!is.na(data.roc)]

#get Youden's J
data.roc.J <- lapply(data.roc, function(row){coords(row, x="best", input="threshold", best.method="youden")})

#drop lines if there are multiple thresholds
data.roc.J <- data.roc.J[sapply(data.roc.J, nrow) == 1]

#prepare df
data.roc.J <- as.data.frame(t(sapply(data.roc.J, '[')))

#get AUC
data.roc.auc <- as.data.frame(sapply(data.roc, '[[', "auc"))

#filter overlapping names and combine both data.frames
data.roc.auc <- data.roc.auc[rownames(data.roc.auc) %in% rownames(data.roc.J), ]

data.roc.df <- cbind(data.roc.auc, data.roc.J)
colnames(data.roc.df)[1] <- "AUC"


#order by AUC
data.roc.df <- data.roc.df %>% rownames_to_column(var = "bin")

res.N.cl <- data.roc.df %>% inner_join(data.mean) %>% arrange(desc(AUC))






##### analyse ROC for SCLC-P ####
data.roc <- data.cl.filt.t


#calculate for SCLC-A
data.roc$Subtype <- ifelse(data.roc$Subtype == "SCLC-P", "SCLC-P", "rest")

#prepare ROC data
data.roc <- pkgcond::suppress_messages(lapply(data.roc[,3:ncol(data.roc)], function(row){tryCatch(roc(data.roc$Subtype, row),
                                                                                                  error = function(e) NA)}))

data.roc <- data.roc[!is.na(data.roc)]

#get Youden's J
data.roc.J <- lapply(data.roc, function(row){coords(row, x="best", input="threshold", best.method="youden")})

#drop lines if there are multiple thresholds
data.roc.J <- data.roc.J[sapply(data.roc.J, nrow) == 1]

#prepare df
data.roc.J <- as.data.frame(t(sapply(data.roc.J, '[')))

#get AUC
data.roc.auc <- as.data.frame(sapply(data.roc, '[[', "auc"))

#filter overlapping names and combine both data.frames
data.roc.auc <- data.roc.auc[rownames(data.roc.auc) %in% rownames(data.roc.J), ]

data.roc.df <- cbind(data.roc.auc, data.roc.J)
colnames(data.roc.df)[1] <- "AUC"


#order by AUC
data.roc.df <- data.roc.df %>% rownames_to_column(var = "bin")

res.P.cl <- data.roc.df %>% inner_join(data.mean) %>% arrange(desc(AUC))




##### analyse ROC for SCLC-I ####
data.roc <- data.cl.filt.t


#calculate for SCLC-A
data.roc$Subtype <- ifelse(data.roc$Subtype == "SCLC-I", "SCLC-I", "rest")

#prepare ROC data
data.roc <- pkgcond::suppress_messages(lapply(data.roc[,3:ncol(data.roc)], function(row){tryCatch(roc(data.roc$Subtype, row),
                                                                                                  error = function(e) NA)}))

data.roc <- data.roc[!is.na(data.roc)]

#get Youden's J
data.roc.J <- lapply(data.roc, function(row){coords(row, x="best", input="threshold", best.method="youden")})

#drop lines if there are multiple thresholds
data.roc.J <- data.roc.J[sapply(data.roc.J, nrow) == 1]

#prepare df
data.roc.J <- as.data.frame(t(sapply(data.roc.J, '[')))

#get AUC
data.roc.auc <- as.data.frame(sapply(data.roc, '[[', "auc"))

#filter overlapping names and combine both data.frames
data.roc.auc <- data.roc.auc[rownames(data.roc.auc) %in% rownames(data.roc.J), ]

data.roc.df <- cbind(data.roc.auc, data.roc.J)
colnames(data.roc.df)[1] <- "AUC"


#order by AUC
data.roc.df <- data.roc.df %>% rownames_to_column(var = "bin")

res.I.cl <- data.roc.df %>% inner_join(data.mean) %>% arrange(desc(AUC))



###################################
#### Select sites for Training ####
###################################

res.A <- data.table::fread("RRBS/Combined/Bins100_New/ROC_FFPE_A_Train_filt90.csv")
res.N <- data.table::fread("RRBS/Combined/Bins100_New/ROC_FFPE_N_Train_filt90.csv")
res.P <- data.table::fread("RRBS/Combined/Bins100_New/ROC_FFPE_P_Train_filt90.csv")
res.I <- data.table::fread("RRBS/Combined/Bins100_New/ROC_FFPE_I_Train_filt90.csv")

res.A.cl <- res.A.cl %>% filter(AUC >= 0.7)
res.N.cl <- res.N.cl %>% filter(AUC >= 0.7)
res.P.cl <- res.P.cl %>% filter(AUC >= 0.7)
res.I.cl <- res.I.cl %>% filter(AUC >= 0.7)


#filter out
res.A.filt <- res.A %>% filter(AUC >= 0.7, abs(`SCLC-A diff`) >= 25) %>% filter(bin %in% res.A.cl$bin)
res.N.filt <- res.N %>% filter(AUC >= 0.7, abs(`SCLC-N diff`) >= 30) %>% filter(bin %in% res.N.cl$bin)
res.P.filt <- res.P %>% filter(AUC >= 0.8, abs(`SCLC-P diff`) >= 35) %>% filter(bin %in% res.P.cl$bin)
res.I.filt <- res.I %>% filter(AUC >= 0.7, abs(`SCLC-I diff`) >= 30) %>% filter(bin %in% res.I.cl$bin)

#create selection
Selection.A <- as.character(res.A.filt$bin)
Selection.N <- as.character(res.N.filt$bin)
Selection.P <- as.character(res.P.filt$bin)
Selection.I <- as.character(res.I.filt$bin)




#################################################
#### Create selection for Training of Models ####
#################################################

#10 sites
sel.A <- c()

for(i in 1:500){
  set.seed(10 + i)
  sel.A[[i]] <- Selection.A[sample(c(rep(TRUE,10), rep(FALSE, length(Selection.A)-10)))]}

sel.A <- do.call(rbind.data.frame, sel.A)
saveRDS(sel.A, paste0("Models_10/Selection_A.Rds"))


sel.N <- c()

for(i in 1:500){
  set.seed(10 + i)
  sel.N[[i]] <- Selection.N[sample(c(rep(TRUE,10), rep(FALSE, length(Selection.N)-10)))]}

sel.N <- do.call(rbind.data.frame, sel.N)
saveRDS(sel.N, paste0("Models_10/Selection_N.Rds"))

sel.P <- c()

for(i in 1:500){
  set.seed(10 + i)
  sel.P[[i]] <- Selection.P[sample(c(rep(TRUE,10), rep(FALSE, length(Selection.P)-10)))]}

sel.P <- do.call(rbind.data.frame, sel.P)
saveRDS(sel.P, paste0("Models_10/Selection_P.Rds"))

sel.I <- c()

for(i in 1:500){
  set.seed(10 + i)
  sel.I[[i]] <- Selection.I[sample(c(rep(TRUE,10), rep(FALSE, length(Selection.I)-10)))]}

sel.I <- do.call(rbind.data.frame, sel.I)
saveRDS(sel.I, paste0("Models_10/Selection_I.Rds"))


#50 sites for A and N and 50 sites for P and I

sel.A <- c()
  
  for(i in 1:500){
    set.seed(50 + i)
    sel.A[[i]] <- Selection.A[sample(c(rep(TRUE,50), rep(FALSE, length(Selection.A)-50)))]}
  
  sel.A <- do.call(rbind.data.frame, sel.A)
  saveRDS(sel.A, paste0("Models_50/Selection_A.Rds"))

  
  sel.N <- c()
  
  for(i in 1:500){
    set.seed(50 + i)
    sel.N[[i]] <- Selection.N[sample(c(rep(TRUE,50), rep(FALSE, length(Selection.N)-50)))]}
  
  sel.N <- do.call(rbind.data.frame, sel.N)
  saveRDS(sel.N, paste0("Models_50/Selection_N.Rds"))
  
  sel.P <- c()
  
  for(i in 1:500){
    set.seed(50 + i)
    sel.P[[i]] <- Selection.P[sample(c(rep(TRUE,50), rep(FALSE, length(Selection.P)-50)))]}
  
  sel.P <- do.call(rbind.data.frame, sel.P)
  saveRDS(sel.P, paste0("Models_50/Selection_P.Rds"))
  
  sel.I <- c()
  
  for(i in 1:500){
    set.seed(50 + i)
    sel.I[[i]] <- Selection.I[sample(c(rep(TRUE,50), rep(FALSE, length(Selection.I)-50)))]}
  
  sel.I <- do.call(rbind.data.frame, sel.I)
  saveRDS(sel.I, paste0("Models_50/Selection_I.Rds"))
  
# for 100 sites
sel.A <- c()
  
for(i in 1:500){
 set.seed(100 + i)
 sel.A[[i]] <- Selection.A[sample(c(rep(TRUE,100), rep(FALSE, length(Selection.A)-100)))]}
  
sel.A <- do.call(rbind.data.frame, sel.A)
saveRDS(sel.A, paste0("Models_100/Selection_A.Rds"))
  
  
sel.N <- c()
  
for(i in 1:500){
  set.seed(100 + i)
  sel.N[[i]] <- Selection.N[sample(c(rep(TRUE,100), rep(FALSE, length(Selection.N)-100)))]}
  
sel.N <- do.call(rbind.data.frame, sel.N)
saveRDS(sel.N, paste0("Models_100/Selection_N.Rds"))
  
sel.P <- c()
  
for(i in 1:500){
  set.seed(100 + i)
  sel.P[[i]] <- Selection.P[sample(c(rep(TRUE,100), rep(FALSE, length(Selection.P)-100)))]}
  
sel.P <- do.call(rbind.data.frame, sel.P)
saveRDS(sel.P, paste0("Models_100/Selection_P.Rds"))
  
sel.I <- c()
  
for(i in 1:500){
  set.seed(100 + i)
  sel.I[[i]] <- Selection.I[sample(c(rep(TRUE,100), rep(FALSE, length(Selection.I)-100)))]}
  
sel.I <- do.call(rbind.data.frame, sel.I)
saveRDS(sel.I, paste0("Models_100/Selection_I.Rds"))


##### Train Models #####

# The training has been carried out on a high-performance cluster.
# The code below highlights the snippet that has been used to train the individual models 


#load libraries
library(tidyverse)
library(caret)


#load variables
index <- commandArgs(trailingOnly = TRUE)

i <- as.numeric(index)


#load data
data <- readRDS("Combined_data_Bins100_Filt10_filt90_annotated_train.Rds")
selection <- readRDS("Models_10/Selection_A.Rds")


ctrlMod <- trainControl(method = "LOOCV",
                        number = 1,
                        classProbs = TRUE,
                        verboseIter = FALSE)



data.model <- data[, colnames(data) %in% c("Subtype", selection[i, ])]
data.model$Subtype <- ifelse(data.model$Subtype == "SCLC-A", "A", "rest")
models <- train(Subtype ~ .,
                data = data.model, 
                method = "xgbDART",
                trControl = ctrlMod,
                na.action = na.pass)


saveRDS(models, paste0("Models_10/Models/Model_A_", i, ".Rds"))


##### load in models #####

#For 10 sites
for(i in c("A", "N", "P", "I")){
  files <- list.files(paste0("Models_10/Models"), pattern = paste0("Model_", i), recursive = TRUE)
  models <- c()
  
  for(m in 1:length(files)){
    models[[m]] <- readRDS(paste0("Models_10/Models/", files[m]))}
  
  saveRDS(models, paste0("Models_10/models_", i, ".Rds"))
}

#For 50 sites
for(i in c("A", "N", "P", "I")){
  files <- list.files(paste0("Models_50/Models"), pattern = paste0("Model_", i), recursive = TRUE)
  models <- c()
  
  for(m in 1:length(files)){
    models[[m]] <- readRDS(paste0("Models_50/Models/", files[m]))}
  
  saveRDS(models, paste0("Models_50/models_", i, ".Rds"))
}


#For 100 sites
for(i in c("A", "N", "P", "I")){
  files <- list.files(paste0("Models_100/Models"), pattern = paste0("Model_", i), recursive = TRUE)
  models <- c()
  
  for(m in 1:length(files)){
    models[[m]] <- readRDS(paste0("Models_100/Models/", files[m]))}
  
  saveRDS(models, paste0("Models_100/models_", i, ".Rds"))
}


#########################################
#### Predict Subtypes on FFPE tissue ####
#########################################

data.train <- readRDS("Combined_data_Bins100_Filt10_filt90_annotated_train.Rds")
data.test <- readRDS("Combined_data_Bins100_Filt10_filt90_annotated_test.Rds")

#run prediction for 10 sites/model
models.A <- readRDS("Models_10/models_A.Rds")
models.N <- readRDS("Models_10/models_N.Rds")
models.P <- readRDS("Models_10/models_P.Rds")
models.I <- readRDS("Models_10/models_I.Rds")
  

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
  
#save to disk
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
  
pdf("Models_10/Heatmap_consensus_FFPE.pdf", width = 12, height = 5)
draw(hm)
dev.off()

#run prediction for 50 sites/model
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
consensus <- read.csv("RRBS/Combined/Bins100_New/Models_50/consensus.csv")
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

pdf("Models_50/Heatmap_consensus_FFPE.pdf", width = 12, height = 5)
draw(hm)
dev.off()

#run prediction for 100 sites/model
models.A <- readRDS("Models_100/models_A.Rds")
models.N <- readRDS("Models_100/models_N.Rds")
models.P <- readRDS("Models_100/models_P.Rds")
models.I <- readRDS("Models_100/models_I.Rds")


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
saveRDS(predict.subtypes, "Models_100/Predictions_all.Rds")

predict.subtypes[predict.subtypes == "rest"] <- NA


freq <- sapply(predict.subtypes, function(x) table(factor(x, levels=unique(unlist(predict.subtypes)), ordered=TRUE)))

freq <- (t(freq) / rowSums(t(freq))) *100

consensus <- cbind.data.frame(freq)
consensus$Consensus <- ifelse(consensus$A >= 50, "A", ifelse(consensus$N >= 50, "N", ifelse(consensus$P >= 50, "P", ifelse(consensus$I >= 50, "I", "equivocal"))))
consensus$ID <- data.all$Sample
consensus$Subtype_RNA <- data.all$Subtype
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

pdf("Models_100/Heatmap_consensus_FFPE.pdf", width = 12, height = 5)
draw(hm)
dev.off()



########################################
#### Train new models on cfDNA data ####
########################################

#load in selection
res.A <- data.table::fread("ROC_FFPE_A_Train_filt90.csv")
res.N <- data.table::fread("ROC_FFPE_N_Train_filt90.csv")
res.P <- data.table::fread("ROC_FFPE_P_Train_filt90.csv")
res.I <- data.table::fread("ROC_FFPE_I_Train_filt90.csv")

res.A.cl <- data.table::fread("ROC_CL_100Bin_A_Train_filt70.csv") %>% filter(AUC >= 0.7)
res.N.cl <- data.table::fread("ROC_CL_100Bin_N_Train_filt70.csv") %>% filter(AUC >= 0.7)
res.P.cl <- data.table::fread("ROC_CL_100Bin_P_Train_filt70.csv") %>% filter(AUC >= 0.7)
res.I.cl <- data.table::fread("ROC_CL_100Bin_I_Train_filt70.csv") %>% filter(AUC >= 0.7)


#filter out
res.A.filt <- res.A %>% filter(AUC >= 0.7, abs(`SCLC-A diff`) >= 25) %>% filter(bin %in% res.A.cl$bin)
res.N.filt <- res.N %>% filter(AUC >= 0.7, abs(`SCLC-N diff`) >= 30) %>% filter(bin %in% res.N.cl$bin)
res.P.filt <- res.P %>% filter(AUC >= 0.8, abs(`SCLC-P diff`) >= 35) %>% filter(bin %in% res.P.cl$bin)
res.I.filt <- res.I %>% filter(AUC >= 0.7, abs(`SCLC-I diff`) >= 30) %>% filter(bin %in% res.I.cl$bin)

Selection.A <- as.character(res.A.filt$bin)
Selection.N <- as.character(res.N.filt$bin)
Selection.P <- as.character(res.P.filt$bin)
Selection.I <- as.character(res.I.filt$bin)



#### Load in cfDNA ####
data.cf <- #load in cfDNA data
data.cf.filt <- data.cf[data.cf$bin %in% c(Selection.A, Selection.N, Selection.P, Selection.I), ]
data.cf.filt.t <- as.data.frame(t(data.cf.filt[, -c(1:2)]))
colnames(data.cf.filt.t) <- data.cf.filt$bin
data.cf.filt.t <- data.cf.filt.t %>% rownames_to_column("Sample")
data.cf.filt.t <- data.cf.filt.t %>% separate(Sample, into = c("Sample", "Timepoint"), sep = "\\.2$", remove = FALSE)
data.cf.filt.t$Timepoint <- ifelse(is.na(data.cf.filt.t$Timepoint), "Baseline", "Progression")


anno <- #read in annotation data
anno.filt <- anno %>% filter(Histology == "SCLC") %>% select(Sample = ID, Subtype = RNA_New)
anno.filt <- cbind.data.frame("Sample" = data.cf.filt.t$Sample) %>% left_join(anno.filt)
data.cf.filt.t <- cbind.data.frame(anno.filt, data.cf.filt.t[, -1])
data.cf.filt.t$Sample <- ifelse(data.cf.filt.t$Timepoint == "Progression", paste0(data.cf.filt.t$Sample, "-2"), data.cf.filt.t$Sample)

#### create models


#10 sites
sel.A <- c()

for(i in 1:500){
  set.seed(10 + i)
  sel.A[[i]] <- Selection.A[sample(c(rep(TRUE,10), rep(FALSE, length(Selection.A)-10)))]}

sel.A <- do.call(rbind.data.frame, sel.A)
saveRDS(sel.A, paste0("Models_10_cf/Selection_A.Rds"))


sel.N <- c()

for(i in 1:500){
  set.seed(10 + i)
  sel.N[[i]] <- Selection.N[sample(c(rep(TRUE,10), rep(FALSE, length(Selection.N)-10)))]}

sel.N <- do.call(rbind.data.frame, sel.N)
saveRDS(sel.N, paste0("Models_10_cf/Selection_N.Rds"))

sel.P <- c()

for(i in 1:500){
  set.seed(10 + i)
  sel.P[[i]] <- Selection.P[sample(c(rep(TRUE,10), rep(FALSE, length(Selection.P)-10)))]}

sel.P <- do.call(rbind.data.frame, sel.P)
saveRDS(sel.P, paste0("Models_10_cf/Selection_P.Rds"))

sel.I <- c()

for(i in 1:500){
  set.seed(10 + i)
  sel.I[[i]] <- Selection.I[sample(c(rep(TRUE,10), rep(FALSE, length(Selection.I)-10)))]}

sel.I <- do.call(rbind.data.frame, sel.I)
saveRDS(sel.I, paste0("Models_10_cf/Selection_I.Rds"))


#50 sites

sel.A <- c()

for(i in 1:500){
  set.seed(50 + i)
  sel.A[[i]] <- Selection.A[sample(c(rep(TRUE,50), rep(FALSE, length(Selection.A)-50)))]}

sel.A <- do.call(rbind.data.frame, sel.A)
saveRDS(sel.A, paste0("Models_50_cf/Selection_A.Rds"))


sel.N <- c()

for(i in 1:500){
  set.seed(50 + i)
  sel.N[[i]] <- Selection.N[sample(c(rep(TRUE,50), rep(FALSE, length(Selection.N)-50)))]}

sel.N <- do.call(rbind.data.frame, sel.N)
saveRDS(sel.N, paste0("Models_50_cf/Selection_N.Rds"))

sel.P <- c()

for(i in 1:500){
  set.seed(50 + i)
  sel.P[[i]] <- Selection.P[sample(c(rep(TRUE,50), rep(FALSE, length(Selection.P)-50)))]}

sel.P <- do.call(rbind.data.frame, sel.P)
saveRDS(sel.P, paste0("Models_50_cf/Selection_P.Rds"))

sel.I <- c()

for(i in 1:500){
  set.seed(50 + i)
  sel.I[[i]] <- Selection.I[sample(c(rep(TRUE,50), rep(FALSE, length(Selection.I)-50)))]}

sel.I <- do.call(rbind.data.frame, sel.I)
saveRDS(sel.I, paste0("Models_50_cf/Selection_I.Rds"))

# for 100 sites
sel.A <- c()

for(i in 1:500){
  set.seed(100 + i)
  sel.A[[i]] <- Selection.A[sample(c(rep(TRUE,100), rep(FALSE, length(Selection.A)-100)))]}

sel.A <- do.call(rbind.data.frame, sel.A)
saveRDS(sel.A, paste0("Models_100_cf/Selection_A.Rds"))


sel.N <- c()

for(i in 1:500){
  set.seed(100 + i)
  sel.N[[i]] <- Selection.N[sample(c(rep(TRUE,100), rep(FALSE, length(Selection.N)-100)))]}

sel.N <- do.call(rbind.data.frame, sel.N)
saveRDS(sel.N, paste0("Models_100_cf/Selection_N.Rds"))

sel.P <- c()

for(i in 1:500){
  set.seed(100 + i)
  sel.P[[i]] <- Selection.P[sample(c(rep(TRUE,100), rep(FALSE, length(Selection.P)-100)))]}

sel.P <- do.call(rbind.data.frame, sel.P)
saveRDS(sel.P, paste0("Models_100_cf/Selection_P.Rds"))

sel.I <- c()

for(i in 1:500){
  set.seed(100 + i)
  sel.I[[i]] <- Selection.I[sample(c(rep(TRUE,100), rep(FALSE, length(Selection.I)-100)))]}

sel.I <- do.call(rbind.data.frame, sel.I)
saveRDS(sel.I, paste0("Models_100_cf/Selection_I.Rds"))



##### Train Models #####

# The training has been carried out on a high-performance cluster.
# The code below highlights the snippet that has been used to train the individual models 


#load libraries
library(tidyverse)
library(caret)


#load variables
index <- commandArgs(trailingOnly = TRUE)

i <- as.numeric(index)


#load data
data <- readRDS("data.Rds")
selection <- readRDS("Models_10_cf/Selection_A.Rds")


ctrlMod <- trainControl(method = "LOOCV",
                        number = 1,
                        classProbs = TRUE,
                        verboseIter = FALSE)



data.model <- data[, colnames(data) %in% c("Subtype", selection[i, ])]
data.model$Subtype <- ifelse(data.model$Subtype == "SCLC-A", "A", "rest")
models <- train(Subtype ~ .,
                data = data.model, 
                method = "xgbDART",
                trControl = ctrlMod,
                na.action = na.pass)


saveRDS(models, paste0("Models_10_cf/Models/Model_A_", i, ".Rds"))






#### Load in Trained Models ####
for(i in c("A", "N", "P", "I")){
   files <- list.files(paste0("Models_10_cf/Models"), pattern = paste0("Model_", i), recursive = TRUE)
   models <- c()
  
   for(m in 1:length(files)){
     models[[m]] <- readRDS(paste0("Models_10_cf/Models/", files[m]))}
   
   saveRDS(models, paste0("Models_10_cf/models_", i, ".Rds"))
}


for(i in c("A", "N", "P", "I")){
  files <- list.files(paste0("Models_50_cf/Models"), pattern = paste0("Model_", i), recursive = TRUE)
  models <- c()
  
  for(m in 1:length(files)){
    models[[m]] <- readRDS(paste0("Models_50_cf/Models/", files[m]))}
  
  saveRDS(models, paste0("Models_50_cf/models_", i, ".Rds"))
}


for(i in c("A", "N", "P", "I")){
  files <- list.files(paste0("Models_100_cf/Models"), pattern = paste0("Model_", i), recursive = TRUE)
  models <- c()
  
  for(m in 1:length(files)){
    models[[m]] <- readRDS(paste0("Models_100_cf/Models/", files[m]))}
  
  saveRDS(models, paste0("Models_100_cf/models_", i, ".Rds"))
}

############################
##### Predict on cfDNA #####
############################

#load in data
data.cf <- #load in data of cfDNA


### for 10 methylation sites ###

#load in models
models.A <- readRDS("Models_10_cf/models_A.Rds")
models.N <- readRDS("Models_10_cf/models_N.Rds")
models.P <- readRDS("Models_10_cf/models_P.Rds")
models.I <- readRDS("Models_10_cf/models_I.Rds")


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

pdf("Models_10_cf/Heatmap_consensus_cfDNA.pdf", width = 12, height = 5)
draw(hm)
dev.off()
  
### for 50 methylation sites ###

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



### for 100 methylation sites ###

#load in models
models.A <- readRDS("Models_100_cf/models_A.Rds")
models.N <- readRDS("Models_100_cf/models_N.Rds")
models.P <- readRDS("Models_100_cf/models_P.Rds")
models.I <- readRDS("Models_100_cf/models_I.Rds")


predict.A <- predict(models.A, data.cf, na.action = na.pass)
predict.N <- predict(models.N, data.cf, na.action = na.pass)
predict.P <- predict(models.P, data.cf, na.action = na.pass)
predict.I <- predict(models.I, data.cf, na.action = na.pass)

predict.subtypes <- rbind.data.frame(t(do.call(cbind.data.frame, predict.A)),
                                     t(do.call(cbind.data.frame, predict.N)),
                                     t(do.call(cbind.data.frame, predict.P)),
                                     t(do.call(cbind.data.frame, predict.I)))

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

pdf("Models_100_cf/Heatmap_consensus_cfDNA.pdf", width = 12, height = 5)
draw(hm)
dev.off()




###############################
#### Predict on Cell lines ####
###############################

data.cl <- #read in cell line data

#run prediction for 10 sites/model
models.A <- readRDS("Models_10/models_A.Rds")
models.N <- readRDS("Models_10/models_N.Rds")
models.P <- readRDS("Models_10/models_P.Rds")
models.I <- readRDS("Models_10/models_I.Rds")

predict.A <- predict(models.A, data.cl, na.action = na.pass)
predict.N <- predict(models.N, data.cl, na.action = na.pass)
predict.P <- predict(models.P, data.cl, na.action = na.pass)
predict.I <- predict(models.I, data.cl, na.action = na.pass)

predict.subtypes <- rbind.data.frame(t(do.call(cbind.data.frame, predict.A)),
                                     t(do.call(cbind.data.frame, predict.N)),
                                     t(do.call(cbind.data.frame, predict.P)),
                                     t(do.call(cbind.data.frame, predict.I)))

predict.subtypes[predict.subtypes == "rest"] <- NA


freq <- sapply(predict.subtypes, function(x) table(factor(x, levels=unique(unlist(predict.subtypes)), ordered=TRUE)))

freq <- (t(freq) / rowSums(t(freq))) *100

consensus <- cbind.data.frame(freq)
consensus$Consensus <- ifelse(consensus$A >= 50, "A", ifelse(consensus$N >= 50, "N", ifelse(consensus$P >= 50, "P", ifelse(consensus$I >= 50, "I", "equivocal"))))
consensus$ID <- data.cl$Sample
consensus$Subtype_RNA <- data.cl$Subtype

#create Heatmap
col.models <- circlize::colorRamp2(c(10, 50, 70), c("#375696", "yellow3", "seagreen3"))
col.subtypes <- c("A" = "darkblue", "N" = "forestgreen", "P" = "firebrick", "I" = "darkgoldenrod1", "equivocal" = "gray", "equivocal" = "gray", "A/N" = "cadetblue4","rest" = "gray90")
col.pred <- c("A" = "darkblue", "N" = "forestgreen", "P" = "firebrick", "I" = "darkgoldenrod1", "equivocal" = "gray", "A/N" = "cadetblue4","rest" = "gray90")

consensus$Subtype_RNA <- sub("SCLC-", "", consensus$Subtype_RNA)

anno <- HeatmapAnnotation("Subtype_DMC" = consensus$Consensus, "Subtype_RNA" = consensus$Subtype_RNA,
                          col = list("Subtype_DMC" = col.pred, "Subtype_RNA" = col.pred), annotation_height = unit(c(0.5, 0.5), "cm"))


data.hm <- t(consensus[, 1:4])
colnames(data.hm) <- consensus$ID

hm <- Heatmap(data.hm, cluster_rows = FALSE, cluster_columns = TRUE, name = "consensus", show_column_names = TRUE, column_names_gp = gpar(fontsize = 4),
              show_row_names = TRUE, col = col.models, top_annotation = anno, column_split = consensus$Cohort,
              height = unit(2.4, "cm"), column_dend_height = unit(2, "cm"), row_dend_width = unit(2, "cm"), row_title = "Consensus")

pdf("Models_10/Heatmap_consensus_FFPE_onCellLines.pdf", width = 12, height = 5)
draw(hm)
dev.off()

#run prediction for 50 sites/model
models.A <- readRDS("Models_50/models_A.Rds")
models.N <- readRDS("Models_50/models_N.Rds")
models.P <- readRDS("Models_50/models_P.Rds")
models.I <- readRDS("Models_50/models_I.Rds")

predict.A <- predict(models.A, data.cl, na.action = na.pass)
predict.N <- predict(models.N, data.cl, na.action = na.pass)
predict.P <- predict(models.P, data.cl, na.action = na.pass)
predict.I <- predict(models.I, data.cl, na.action = na.pass)

predict.subtypes <- rbind.data.frame(t(do.call(cbind.data.frame, predict.A)),
                                     t(do.call(cbind.data.frame, predict.N)),
                                     t(do.call(cbind.data.frame, predict.P)),
                                     t(do.call(cbind.data.frame, predict.I)))

predict.subtypes[predict.subtypes == "rest"] <- NA


freq <- sapply(predict.subtypes, function(x) table(factor(x, levels=unique(unlist(predict.subtypes)), ordered=TRUE)))

freq <- (t(freq) / rowSums(t(freq))) *100

consensus <- cbind.data.frame(freq)
consensus$Consensus <- ifelse(consensus$A >= 50, "A", ifelse(consensus$N >= 50, "N", ifelse(consensus$P >= 50, "P", ifelse(consensus$I >= 50, "I", "equivocal"))))
consensus$ID <- data.cl$Sample
consensus$Subtype_RNA <- data.cl$Subtype

#create Heatmap
col.models <- circlize::colorRamp2(c(10, 50, 70), c("#375696", "yellow3", "seagreen3"))
col.subtypes <- c("A" = "darkblue", "N" = "forestgreen", "P" = "firebrick", "I" = "darkgoldenrod1", "equivocal" = "gray", "equivocal" = "gray", "A/N" = "cadetblue4","rest" = "gray90")
col.pred <- c("A" = "darkblue", "N" = "forestgreen", "P" = "firebrick", "I" = "darkgoldenrod1", "equivocal" = "gray", "A/N" = "cadetblue4","rest" = "gray90")

consensus$Subtype_RNA <- sub("SCLC-", "", consensus$Subtype_RNA)

anno <- HeatmapAnnotation("Subtype_DMC" = consensus$Consensus, "Subtype_RNA" = consensus$Subtype_RNA,
                          col = list("Subtype_DMC" = col.pred, "Subtype_RNA" = col.pred), annotation_height = unit(c(0.5, 0.5), "cm"))


data.hm <- t(consensus[, 1:4])
colnames(data.hm) <- consensus$ID

hm <- Heatmap(data.hm, cluster_rows = FALSE, cluster_columns = TRUE, name = "consensus", show_column_names = TRUE, column_names_gp = gpar(fontsize = 4),
              show_row_names = TRUE, col = col.models, top_annotation = anno, column_split = consensus$Cohort,
              height = unit(2.4, "cm"), column_dend_height = unit(2, "cm"), row_dend_width = unit(2, "cm"), row_title = "Consensus")

pdf("Models_50/Heatmap_consensus_FFPE_onCellLines.pdf", width = 12, height = 5)
draw(hm)
dev.off()

#run prediction for 100 sites/model
models.A <- readRDS("Models_100/models_A.Rds")
models.N <- readRDS("Models_100/models_N.Rds")
models.P <- readRDS("Models_100/models_P.Rds")
models.I <- readRDS("Models_100/models_I.Rds")

predict.A <- predict(models.A, data.cl, na.action = na.pass)
predict.N <- predict(models.N, data.cl, na.action = na.pass)
predict.P <- predict(models.P, data.cl, na.action = na.pass)
predict.I <- predict(models.I, data.cl, na.action = na.pass)

predict.subtypes <- rbind.data.frame(t(do.call(cbind.data.frame, predict.A)),
                                     t(do.call(cbind.data.frame, predict.N)),
                                     t(do.call(cbind.data.frame, predict.P)),
                                     t(do.call(cbind.data.frame, predict.I)))

predict.subtypes[predict.subtypes == "rest"] <- NA


freq <- sapply(predict.subtypes, function(x) table(factor(x, levels=unique(unlist(predict.subtypes)), ordered=TRUE)))

freq <- (t(freq) / rowSums(t(freq))) *100

consensus <- cbind.data.frame(freq)
consensus$Consensus <- ifelse(consensus$A >= 50, "A", ifelse(consensus$N >= 50, "N", ifelse(consensus$P >= 50, "P", ifelse(consensus$I >= 50, "I", "equivocal"))))
consensus$ID <- data.cl$Sample
consensus$Subtype_RNA <- data.cl$Subtype

#create Heatmap
col.models <- circlize::colorRamp2(c(10, 50, 70), c("#375696", "yellow3", "seagreen3"))
col.subtypes <- c("A" = "darkblue", "N" = "forestgreen", "P" = "firebrick", "I" = "darkgoldenrod1", "equivocal" = "gray", "equivocal" = "gray", "A/N" = "cadetblue4","rest" = "gray90")
col.pred <- c("A" = "darkblue", "N" = "forestgreen", "P" = "firebrick", "I" = "darkgoldenrod1", "equivocal" = "gray", "A/N" = "cadetblue4","rest" = "gray90")

consensus$Subtype_RNA <- sub("SCLC-", "", consensus$Subtype_RNA)

anno <- HeatmapAnnotation("Subtype_DMC" = consensus$Consensus, "Subtype_RNA" = consensus$Subtype_RNA,
                          col = list("Subtype_DMC" = col.pred, "Subtype_RNA" = col.pred), annotation_height = unit(c(0.5, 0.5), "cm"))


data.hm <- t(consensus[,c("A", "N", "P", "I")])
colnames(data.hm) <- consensus$ID

hm <- Heatmap(data.hm, cluster_rows = FALSE, cluster_columns = TRUE, name = "consensus", show_column_names = TRUE, column_names_gp = gpar(fontsize = 4),
              show_row_names = TRUE, col = col.models, top_annotation = anno, column_split = consensus$Cohort,
              height = unit(2.4, "cm"), column_dend_height = unit(2, "cm"), row_dend_width = unit(2, "cm"), row_title = "Consensus")

pdf("Models_100/Heatmap_consensus_FFPE_onCellLines.pdf", width = 12, height = 5)
draw(hm)
dev.off()