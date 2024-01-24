##################
#### Figure 2 ####
##################

#load packages
library(tidyverse)
library(karyoploteR)
library(GenomicRanges)


#### Figure 2A ####

# load in annotation
anno <- read.csv("anno.csv")

#load in data of 100kbp bins of DNA methylation data
data.bin <- readRDS("100kBins.Rds")

#load in a data frame with the region annotation for each of the bins
bins.df <- readRDS("100kBins_anno.Rds")

#filter for only the autosomes
bins.df <- bins.df[bins.df$seqnames %in% 1:22, ]

#combine data and calculate rolling average
data.bin.comb <- cbind.data.frame("Sample" = data.bin$Sample) %>% left_join(select(data.bin, Sample, consensus)) %>% cbind.data.frame(data.bin[, -1])
data.bin.comb.mean <- data.bin.comb %>% filter(consensus %in% c("A", "N", "P", "I")) 
data.bin.comb.mean <- aggregate(data.bin.comb.mean[, -c(1:3)], by = list(c(data.bin.comb.mean$consensus)), mean, na.rm = TRUE)
data.bin.comb.mean <- as.data.frame(t(data.bin.comb.mean[, -1]))
colnames(data.bin.comb.mean) <- c("A", "I", "N", "P")
data.bin.comb.mean$bin <- data.bin$bin

data.bin.comb.mean$RMean_A <- zoo::rollapply(data.bin.comb.mean$A, width=500, FUN=mean, na.rm = TRUE, fill=NA, align="center")
data.bin.comb.mean$RMean_N <- zoo::rollapply(data.bin.comb.mean$N, width=500, FUN=mean, na.rm = TRUE, fill=NA, align="center")
data.bin.comb.mean$RMean_P <- zoo::rollapply(data.bin.comb.mean$P, width=500, FUN=mean, na.rm = TRUE, fill=NA, align="center")
data.bin.comb.mean$RMean_I <- zoo::rollapply(data.bin.comb.mean$I, width=500, FUN=mean, na.rm = TRUE, fill=NA, align="center")

#create Figure
data.plot <- reshape2::melt(data.bin.comb.mean, id.vars = c("bin", "RMean_A", "RMean_N", "RMean_P", "RMean_I"))
data.plot$variable <- as.numeric(as.character(data.plot$variable))
data.plot <- arrange(data.plot, bin, variable)


data.plot <- data.plot %>% filter(bin <= max(as.numeric(bins.df$bin)))

colnames(data.plot)[6] <- "Subtype"

col.subtypes <- c("A" = "darkblue", "N" = "forestgreen", "P" = "firebrick", "I" = "darkgoldenrod1")

plot <- ggplot(data.plot, aes(x = bin, y = value, col = Subtype)) + 
  geom_vline(xintercept = max(bins.df[bins.df$seqnames == "1", 6]), lty = "dashed", col = "gray80") +
  geom_vline(xintercept = max(bins.df[bins.df$seqnames == "2", 6]), lty = "dashed", col = "gray80") +
  geom_vline(xintercept = max(bins.df[bins.df$seqnames == "3", 6]), lty = "dashed", col = "gray80") +
  geom_vline(xintercept = max(bins.df[bins.df$seqnames == "4", 6]), lty = "dashed", col = "gray80") +
  geom_vline(xintercept = max(bins.df[bins.df$seqnames == "5", 6]), lty = "dashed", col = "gray80") +
  geom_vline(xintercept = max(bins.df[bins.df$seqnames == "6", 6]), lty = "dashed", col = "gray80") +
  geom_vline(xintercept = max(bins.df[bins.df$seqnames == "7", 6]), lty = "dashed", col = "gray80") +
  geom_vline(xintercept = max(bins.df[bins.df$seqnames == "8", 6]), lty = "dashed", col = "gray80") +
  geom_vline(xintercept = max(bins.df[bins.df$seqnames == "9", 6]), lty = "dashed", col = "gray80") +
  geom_vline(xintercept = max(bins.df[bins.df$seqnames == "10", 6]), lty = "dashed", col = "gray80") +
  geom_vline(xintercept = max(bins.df[bins.df$seqnames == "11", 6]), lty = "dashed", col = "gray80") +
  geom_vline(xintercept = max(bins.df[bins.df$seqnames == "12", 6]), lty = "dashed", col = "gray80") +
  geom_vline(xintercept = max(bins.df[bins.df$seqnames == "13", 6]), lty = "dashed", col = "gray80") +
  geom_vline(xintercept = max(bins.df[bins.df$seqnames == "14", 6]), lty = "dashed", col = "gray80") +
  geom_vline(xintercept = max(bins.df[bins.df$seqnames == "15", 6]), lty = "dashed", col = "gray80") +
  geom_vline(xintercept = max(bins.df[bins.df$seqnames == "16", 6]), lty = "dashed", col = "gray80") +
  geom_vline(xintercept = max(bins.df[bins.df$seqnames == "17", 6]), lty = "dashed", col = "gray80") +
  geom_vline(xintercept = max(bins.df[bins.df$seqnames == "18", 6]), lty = "dashed", col = "gray80") +
  geom_vline(xintercept = max(bins.df[bins.df$seqnames == "19", 6]), lty = "dashed", col = "gray80") +
  geom_vline(xintercept = max(bins.df[bins.df$seqnames == "20", 6]), lty = "dashed", col = "gray80") +
  geom_vline(xintercept = max(bins.df[bins.df$seqnames == "21", 6]), lty = "dashed", col = "gray80") + theme_bw() + 
  geom_line(data = data.plot[!duplicated(data.plot$bin) & !is.na(data.plot$RMean_A), ], aes(x = bin, y = RMean_A), col = "darkblue", size = 1) +
  geom_line(data = data.plot[!duplicated(data.plot$bin) & !is.na(data.plot$RMean_N), ], aes(x = bin, y = RMean_N), col = "forestgreen", size = 1) +
  geom_line(data = data.plot[!duplicated(data.plot$bin) & !is.na(data.plot$RMean_P), ], aes(x = bin, y = RMean_P), col = "firebrick", size = 1) +
  geom_line(data = data.plot[!duplicated(data.plot$bin) & !is.na(data.plot$RMean_I), ], aes(x = bin, y = RMean_I), col = "darkgoldenrod1", size = 1) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_color_manual(values = col.subtypes) + ylim(0,100)

pdf(paste0(Sys.Date(), "_Average_Methylation_100kBins.pdf"), width = 25, height = 10)
plot
dev.off()


######


#### Figure 2B-F ####

#Read in data
data <- read.csv("expression_data.csv")
data <- data[data$classification %in% c("SCLC-A", "SCLC-N", "SCLC-P", "SCLC-I"), ]

#Specify x-axis order to match desired A/N/P/I ouptut
x_order <- c('SCLC-A', 'SCLC-N', 'SCLC-P', 'SCLC-I')

#Define comparisons for statistical analysis
my_comparisons <- list(c("SCLC-P", "SCLC-N"), c("SCLC-P", "SCLC-I"), c("SCLC-P", "SCLC-A"),
                       c("SCLC-N", "SCLC-I"), c("SCLC-N", "SCLC-A"), c("SCLC-I", "SCLC-A"))

#Draw boxplot plot for each of the respective genes


plot <- ggplot(data, aes(x = factor(data$classification, levels = x_order), y = data$GOI, #replace with Gene you want to analyse
                                                fill = data$classification)) +
  geom_boxplot(outlier.shape = 1,
               width = 0.55) +
  stat_boxplot(geom = "errorbar", width = 0.1) +
  theme_bw() +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.text = element_text(size = 10, color = "black"),
        axis.title.y = element_text(size = 12, color = "black"),
        plot.title = element_text(hjust = 0.5),
        legend.position = "None") +
  stat_compare_means(comparisons = my_comparisons) +
  labs(title = "GOI") + #replace with Gene of interest
  scale_fill_manual(values = c("SCLC-A" = "darkblue", 
                               "SCLC-N" = "forestgreen",
                               "SCLC-I" = "darkgoldenrod1",
                               "SCLC-P" = "firebrick")) +
  scale_y_continuous("mRNA expression (TPM)")


pdf(paste0(Sys.Date(), "_plot.pdf", width = 8, height = 8)
plot
dev.off()

######


#### Figure 2I-L ####



##### For SCLC-A

#load in SCLC-A data
data.roc.A  <- read.csv("ROC_A_Train_annotated_AUC80.csv")
data.roc.A <- data.roc.A %>% separate(Chr_Range, into = c("chr", "range"), sep = ":", remove = TRUE) %>% separate(range, into = c("start", "end"), remove = TRUE)

#create genomic ranges data
data.a.gr <- makeGRangesListFromDataFrame(select(data.roc.A, chr, start, end, AUC))

#create plot

pdf(paste0(Sys.Date(), "_Karyogram_AUC80_A.pdf"), width = 8, height = 4)

plot <- plotKaryotype(chromosomes=c("autosomal"))
kpPlotRegions(plot, unlist(data.a.gr), col = "darkblue")

dev.off()




##### For SCLC-N

#load in SCLC-N data
data.roc.N  <- read.csv("ROC_N_Train_annotated_AUC80.csv")
data.roc.N <- data.roc.N %>% separate(Chr_Range, into = c("chr", "range"), sep = ":", remove = TRUE) %>% separate(range, into = c("start", "end"), remove = TRUE)

#create genomic ranges data
data.n.gr <- makeGRangesListFromDataFrame(select(data.roc.N, chr, start, end, AUC))

#create plot

pdf(paste0(Sys.Date(), "_Karyogram_AUC80_N.pdf"), width = 8, height = 4)

plot <- plotKaryotype(chromosomes=c("autosomal"))
kpPlotRegions(plot, unlist(data.n.gr), col = "darkblue")

dev.off()




##### For SCLC-P

#load in SCLC-P data
data.roc.P  <- read.csv("ROC_P_Train_annotated_AUC80.csv")
data.roc.P <- data.roc.P %>% separate(Chr_Range, into = c("chr", "range"), sep = ":", remove = TRUE) %>% separate(range, into = c("start", "end"), remove = TRUE)

#create genomic ranges data
data.a.gr <- makeGRangesListFromDataFrame(select(data.roc.P, chr, start, end, AUC))

#create plot

pdf(paste0(Sys.Date(), "_Karyogram_AUC80_P.pdf"), width = 8, height = 4)

plot <- plotKaryotype(chromosomes=c("autosomal"))
kpPlotRegions(plot, unlist(data.p.gr), col = "darkblue")

dev.off()



##### For SCLC-I

#load in SCLC-A data
data.roc.I  <- read.csv("ROC_I_Train_annotated_AUC80.csv")
data.roc.I <- data.roc.I %>% separate(Chr_Range, into = c("chr", "range"), sep = ":", remove = TRUE) %>% separate(range, into = c("start", "end"), remove = TRUE)

#create genomic ranges data
data.i.gr <- makeGRangesListFromDataFrame(select(data.roc.I, chr, start, end, AUC))

#create plot

pdf(paste0(Sys.Date(), "_Karyogram_AUC80_I.pdf"), width = 8, height = 4)

plot <- plotKaryotype(chromosomes=c("autosomal"))
kpPlotRegions(plot, unlist(data.i.gr), col = "darkblue")

dev.off()








