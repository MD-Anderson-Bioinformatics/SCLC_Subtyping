################################################################################################
### This code describes the selection of appropriate genes and training of predictive models ###
################################################################################################

library(ComplexHeatmap)
library(tidyverse)
library(rstatix)
library(pROC)
library(caret)



###############################################################
### Define top genes that are associated with SCLC subtypes ###
###############################################################

#load in data (data can be of various input such as TPM or FPKM and log2 or linear)
data <- read.csv("data.csv")
  
#the code expects that the frist column have information of sample name and teh second column teh classification while each additional column is a gene with gene expression data. 
#each sample is in one row. If th eformat is different you need to adjust accordingly.


# for each subtype, the highest association is computed using ROC analysis


### For SCLC-A
data.A <- data
data.A$Classification <- ifelse(data.A$Classification == "A", "A", "rest")

#prepare ROC data
data.A.roc <- pkgcond::suppress_messages(lapply(data.A[,3:ncol(data.A)], function(row){roc(data.A$Classification, row)}))

#get Youden's J
data.A.roc.J <- lapply(data.A.roc, function(row){coords(row, x="best", input="threshold", best.method="youden")})

#drop lines if there are multiple thresholds
data.A.roc.J <- data.A.roc.J[sapply(data.A.roc.J, nrow) == 1]

#prepare df
data.A.roc.J <- as.data.frame(t(sapply(data.A.roc.J, '[')))

#get AUC
data.A.roc.auc <- as.data.frame(sapply(data.A.roc, '[[', "auc"))

#filter overlapping names and combine both data.frames
data.A.roc.auc <- data.A.roc.auc[rownames(data.A.roc.auc) %in% rownames(data.A.roc.J), ]

data.A.roc.df <- cbind(data.A.roc.auc, data.A.roc.J)
colnames(data.A.roc.df)[1] <- "AUC"


#order by AUC
data.A.roc.df <- data.A.roc.df %>% rownames_to_column(var = "Gene") %>% 
  arrange(desc(AUC))


#safe to disk
saveRDS(data.A.roc.df, file = "Data_A_ROC.Rds")


### For SCLC-N
data.N <- data
data.N$Classification <- ifelse(data.N$Classification == "N", "N", "rest")

#prepare ROC data
data.N.roc <- pkgcond::suppress_messages(lapply(data.N[,3:ncol(data.N)], function(row){roc(data.N$Classification, row)}))

#get Youden's J
data.N.roc.J <- lapply(data.N.roc, function(row){coords(row, x="best", input="threshold", best.method="youden")})

#drop lines if there are multiple thresholds
data.N.roc.J <- data.N.roc.J[sapply(data.N.roc.J, nrow) == 1]

#prepare df
data.N.roc.J <- as.data.frame(t(sapply(data.N.roc.J, '[')))

#get AUC
data.N.roc.auc <- as.data.frame(sapply(data.N.roc, '[[', "auc"))

#filter overlapping names and combine both data.frames
data.N.roc.auc <- data.N.roc.auc[rownames(data.N.roc.auc) %in% rownames(data.N.roc.J), ]

data.N.roc.df <- cbind(data.N.roc.auc, data.N.roc.J)
colnames(data.N.roc.df)[1] <- "AUC"


#order by AUC
data.N.roc.df <- data.N.roc.df %>% rownames_to_column(var = "Gene") %>% 
  arrange(desc(AUC))

#safe to disk
saveRDS(data.N.roc.df, file = "Data_N_ROC.Rds")


### For SCLC-P
data.P <- data
data.P$Classification <- ifelse(data.P$Classification == "P", "P", "rest")

#prepare ROC data
data.P.roc <- pkgcond::suppress_messages(lapply(data.P[,3:ncol(data.P)], function(row){roc(data.P$Classification, row)}))

#get Youden's J
data.P.roc.J <- lapply(data.P.roc, function(row){coords(row, x="best", input="threshold", best.method="youden")})

#drop lines if there are multiple thresholds
data.P.roc.J <- data.P.roc.J[sapply(data.P.roc.J, nrow) == 1]

#prepare df
data.P.roc.J <- as.data.frame(t(sapply(data.P.roc.J, '[')))

#get AUC
data.P.roc.auc <- as.data.frame(sapply(data.P.roc, '[[', "auc"))

#filter overlapping names and combine both data.frames
data.P.roc.auc <- data.P.roc.auc[rownames(data.P.roc.auc) %in% rownames(data.P.roc.J), ]

data.P.roc.df <- cbind(data.P.roc.auc, data.P.roc.J)
colnames(data.P.roc.df)[1] <- "AUC"


#order by AUC
data.P.roc.df <- data.P.roc.df %>% rownames_to_column(var = "Gene") %>% 
  arrange(desc(AUC))


#safe to disk
saveRDS(data.P.roc.df, file = "Data_P_ROC.Rds")


### For SCLC-I
data.I <- data
data.I$Classification <- ifelse(data.I$Classification == "I", "I", "rest")

#prepare ROC data
data.I.roc <- pkgcond::suppress_messages(lapply(data.I[,3:ncol(data.I)], function(row){roc(data.I$Classification, row)}))

#get Youden's J
data.I.roc.J <- lapply(data.I.roc, function(row){coords(row, x="best", input="threshold", best.method="youden")})

#drop lines if there are multiple thresholds
data.I.roc.J <- data.I.roc.J[sapply(data.I.roc.J, nrow) == 1]

#prepare df
data.I.roc.J <- as.data.frame(t(sapply(data.I.roc.J, '[')))

#get AUC
data.I.roc.auc <- as.data.frame(sapply(data.I.roc, '[[', "auc"))

#filter overlapping names and combine both data.frames
data.I.roc.auc <- data.I.roc.auc[rownames(data.I.roc.auc) %in% rownames(data.I.roc.J), ]

data.I.roc.df <- cbind(data.I.roc.auc, data.I.roc.J)
colnames(data.I.roc.df)[1] <- "AUC"


#order by AUC
data.I.roc.df <- data.I.roc.df %>% rownames_to_column(var = "Gene") %>% 
  arrange(desc(AUC))

#safe to disk
saveRDS(data.I.roc.df, file = "Data_I_ROC.Rds")



#select genes for training
selection <- c(data.A.roc.df$Gene[1:50], 
               data.N.roc.df$Gene[1:50],
               data.P.roc.df$Gene[1:50],
               data.I.roc.df$Gene[1:50])


#########################################
### Create gene ratios from selection ###
#########################################


combinations <- as.data.frame(t(combn(selection, 2)))

combinations$formula <- paste0(combinations$V1, "/", combinations$V2)

#create new table
ratios <- c()
for(i in 1:nrow(combinations)){
  j <- combinations[i, 3]
  ratios[[j]] <- select(data, combinations[i, 1]) / select(data, combinations[i, 2])
  
}

#combine data to dataframe. If data is log2 transformed the formula needs to be adjusted to value = 2^(geneA - geneB)
ratios.df <- do.call(cbind.data.frame, ratios)
colnames(ratios.df) <- names(ratios)

ratios.df <- cbind.data.frame("Sample" = data$Sample,
                              "Classification" = data$Classification,
                              ratios.df)

  
#save to disk
saveRDS(ratios.df, "ratios_df.Rds")

#################################################
#### Create selection for Training of Models ####
#################################################


#set seed as you prefer to any numeric value (replace A, N, P, I) but additiona of number ensures that each time new combinations are selected

sel.A <- c()

for(i in 1:500){
  set.seed(A + i)
  sel.A[[i]] <- combinations[sample(c(rep(TRUE,20), rep(FALSE, length(combinations)-20)))]}

sel.A <- do.call(rbind.data.frame, sel.A)
saveRDS(sel.A, paste0("Models_20/Selection_A.Rds"))

sel.N <- c()

for(i in 1:500){
  set.seed(N + i)
  sel.N[[i]] <- combinations[sample(c(rep(TRUE,20), rep(FALSE, length(combinations)-20)))]}

sel.N <- do.call(rbind.data.frame, sel.N)
saveRDS(sel.N, paste0("Models_20/Selection_N.Rds"))

sel.P <- c()

for(i in 1:500){
  set.seed(P + i)
  sel.P[[i]] <- combinations[sample(c(rep(TRUE,20), rep(FALSE, length(combinations)-20)))]}

sel.P <- do.call(rbind.data.frame, sel.P)
saveRDS(sel.P, paste0("Models_20/Selection_P.Rds"))

sel.I <- c()

for(i in 1:500){
  set.seed(I + i)
  sel.I[[i]] <- combinations[sample(c(rep(TRUE,10), rep(FALSE, length(combinations)-10)))]}

sel.I <- do.call(rbind.data.frame, sel.I)
saveRDS(sel.I, paste0("Models_20/Selection_I.Rds"))



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
ratios <- readRDS("ratios_df.Rds")
selection <- readRDS("Models_20/Selection_A.Rds")



ctrlMod <- trainControl(method = "repeatedcv",
                        repeats = 20,
                        number = 5,
                        classProbs = TRUE,
                        verboseIter = TRUE)



data.model <- ratios[, colnames(ratios) %in% c("Subtype", selection[i, ])]
models <- train(Subtype ~ .,
                data = data.model, 
                method = "xgbDART",
                trControl = ctrlMod)


saveRDS(models, paste0("Models_20/Model_A_", i, ".Rds"))


##########################
##### load in models #####
##########################

for(i in c("A", "N", "P", "I")){
  files <- list.files(paste0("Models_20/Models"), pattern = paste0("Model_", i), recursive = TRUE)
  models <- c()
  
  for(m in 1:length(files)){
    models[[m]] <- readRDS(paste0("Models_20/Models/", files[m]))}
  
  saveRDS(models, paste0("Models_20/models_", i, ".Rds"))
}

