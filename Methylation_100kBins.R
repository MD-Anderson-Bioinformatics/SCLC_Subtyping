
###########################################################
### This code generates the global methylation analysis ###
###########################################################

#load in packages
library(BSgenome.Hsapiens.NCBI.GRCh38)
library(GenomicRanges)
library(tidyverse)


#create 100K bins

data.1 <- #load in methylation data for cohort 1
data.2 <- #load in methylation data for cohort 2


#create 100k bins for data.1
data.1$Chr <- sub("chr", "", data.1$Chr)
bins <- as.data.frame(tileGenome(seqinfo(Hsapiens), tilewidth=100000, cut.last.tile.in.chrom=TRUE))
bins$bin <- 1:nrow(bins)
bins <- as(bins, "GRanges")
bins.df <- as.data.frame(bins)


#create average methylaton per bin for cohort 1
data.gr <- makeGRangesFromDataFrame(data.1, seqnames.field = "Chr", start.field = "Pos", end.field = "Pos", keep.extra.columns = TRUE)
data.gr.merge <- as.data.frame(mergeByOverlaps(bins, data.gr))
data.gr.merge <- as.data.frame(data.gr.merge)
data.gr.merge <- data.gr.merge[, c(7, 96:178)]

data.bin.mean <- data.gr.merge %>% 
  group_by(bin) %>%
  summarise_all(mean, na.rm = TRUE)

bins.df <- bins.df %>% select(bin, seqnames, start, end) %>% mutate(Position = paste0("chr", seqnames, ":", start, "-", end)) %>% select(bin, Position) %>% filter(bin %in% data.gr.merge$bin)
data.bin.c1 <- cbind.data.frame(bins.df, data.bin.mean[, -1])


#create average methylaton per bin for cohort 2
data.2$Chr <- sub("chr", "", data.2$Chr)
data.gr <- makeGRangesFromDataFrame(data.2, seqnames.field = "Chr", start.field = "Pos", end.field = "Pos", keep.extra.columns = TRUE)
data.gr.merge <- as.data.frame(mergeByOverlaps(bins, data.gr))

data.gr.merge <- as.data.frame(data.gr.merge)
data.gr.merge <- data.gr.merge[, c(7, 60:106)]

data.bin.mean <- data.gr.merge %>% 
  group_by(bin) %>%
  summarise_all(mean, na.rm = TRUE)

bins.df <- as.data.frame(bins) %>% select(bin, seqnames, start, end) %>% mutate(Position = paste0("chr", seqnames, ":", start, "-", end)) %>% select(bin, Position) %>% filter(bin %in% data.gr.merge$bin)
data.bin.c2 <- cbind.data.frame(bins.df, data.bin.mean[, -1])


### create mean per subtype and Figure

consensus <- #load in consensus classification for cohort

#calculate for each cohort individually
data.bin.c1.t <- t(data.bin.c1[, -c(1:2)])
colnames(data.bin.c1.t) <- data.bin.c1$bin
data.bin.c1.t <- as.data.frame(data.bin.c1.t) %>% rownames_to_column("Sample")

data.bin.c1.comb <- cbind.data.frame("Sample" = data.bin.c1.t$Sample) %>% left_join(select(data.combined, Sample, consensus, Cohort)) %>% cbind.data.frame(data.bin.c1.t[, -1])
data.bin.c1.comb.mean <- data.bin.c1.comb %>% filter(consensus %in% c("A", "N", "P", "I")) 
data.bin.c1.comb.mean <- aggregate(data.bin.c1.comb.mean[, -c(1:3)], by = list(c(data.bin.c1.comb.mean$consensus)), mean, na.rm = TRUE)
data.bin.c1.comb.mean <- as.data.frame(t(data.bin.c1.comb.mean[, -1]))
colnames(data.bin.c1.comb.mean) <- c("A", "I", "N", "P")
data.bin.c1.comb.mean$bin <- data.bin.c1$bin

data.bin.c1.comb.mean$RMean_A <- zoo::rollapply(data.bin.c1.comb.mean$A, width=500, FUN=mean, na.rm = TRUE, fill=NA, align="center")
data.bin.c1.comb.mean$RMean_N <- zoo::rollapply(data.bin.c1.comb.mean$N, width=500, FUN=mean, na.rm = TRUE, fill=NA, align="center")
data.bin.c1.comb.mean$RMean_P <- zoo::rollapply(data.bin.c1.comb.mean$P, width=500, FUN=mean, na.rm = TRUE, fill=NA, align="center")
data.bin.c1.comb.mean$RMean_I <- zoo::rollapply(data.bin.c1.comb.mean$I, width=500, FUN=mean, na.rm = TRUE, fill=NA, align="center")



#create Figure
data.plot <- reshape2::melt(data.bin.c1.comb.mean, id.vars = c("bin", "RMean_A", "RMean_N", "RMean_P", "RMean_I"))
data.plot$variable <- as.numeric(as.character(data.plot$variable))
data.plot <- arrange(data.plot, bin, variable)


data.plot <- data.plot %>% filter(bin <= max(as.numeric(bins.df$bin)))

colnames(data.plot)[6] <- "Subtype"

col.subtypes <- c("A" = "darkblue", "N" = "forestgreen", "P" = "firebrick", "I" = "darkgoldenrod1", "equivocal" = "gray")

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
  geom_vline(xintercept = max(bins.df[bins.df$seqnames == "21", 6]), lty = "dashed", col = "gray80") +
  theme_bw() + 
  geom_line(data = data.plot[!duplicated(data.plot$bin) & !is.na(data.plot$RMean_A), ], aes(x = bin, y = RMean_A), col = "darkblue", size = 1) +
  geom_line(data = data.plot[!duplicated(data.plot$bin) & !is.na(data.plot$RMean_N), ], aes(x = bin, y = RMean_N), col = "forestgreen", size = 1) +
  geom_line(data = data.plot[!duplicated(data.plot$bin) & !is.na(data.plot$RMean_P), ], aes(x = bin, y = RMean_P), col = "firebrick", size = 1) +
  geom_line(data = data.plot[!duplicated(data.plot$bin) & !is.na(data.plot$RMean_I), ], aes(x = bin, y = RMean_I), col = "darkgoldenrod1", size = 1) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_color_manual(values = col.subtypes)

pdf("Average_Methylation_100kBins_C1.pdf", width = 25, height = 10)
plot
dev.off()



#cohort 2
data.bin.c2.t <- t(data.bin.c2[, -c(1:2)])
colnames(data.bin.c2.t) <- data.bin.c2$bin
data.bin.c2.t <- as.data.frame(data.bin.c2.t) %>% rownames_to_column("Sample")

data.bin.c2.comb <- cbind.data.frame("Sample" = data.bin.c2.t$Sample) %>% left_join(select(data.combined, Sample, consensus, Cohort)) %>% cbind.data.frame(data.bin.c2.t[, -1])
data.bin.c2.comb.mean <- data.bin.c2.comb %>% filter(consensus %in% c("A", "N", "P", "I")) 
data.bin.c2.comb.mean <- aggregate(data.bin.c2.comb.mean[, -c(1:3)], by = list(c(data.bin.c2.comb.mean$consensus)), mean, na.rm = TRUE)
data.bin.c2.comb.mean <- as.data.frame(t(data.bin.c2.comb.mean[, -1]))
colnames(data.bin.c2.comb.mean) <- c("A", "I", "N", "P")
data.bin.c2.comb.mean$bin <- data.bin.c2$bin

data.bin.c2.comb.mean$RMean_A <- zoo::rollapply(data.bin.c2.comb.mean$A, width=500, FUN=mean, na.rm = TRUE, fill=NA, align="center")
data.bin.c2.comb.mean$RMean_N <- zoo::rollapply(data.bin.c2.comb.mean$N, width=500, FUN=mean, na.rm = TRUE, fill=NA, align="center")
data.bin.c2.comb.mean$RMean_P <- zoo::rollapply(data.bin.c2.comb.mean$P, width=500, FUN=mean, na.rm = TRUE, fill=NA, align="center")
data.bin.c2.comb.mean$RMean_I <- zoo::rollapply(data.bin.c2.comb.mean$I, width=500, FUN=mean, na.rm = TRUE, fill=NA, align="center")



#create Figure
data.plot <- reshape2::melt(data.bin.c2.comb.mean, id.vars = c("bin", "RMean_A", "RMean_N", "RMean_P", "RMean_I"))
data.plot$variable <- as.numeric(as.character(data.plot$variable))
data.plot <- arrange(data.plot, bin, variable)


data.plot <- data.plot %>% filter(bin <= max(as.numeric(bins.df$bin)))

colnames(data.plot)[6] <- "Subtype"

col.subtypes <- c("A" = "darkblue", "N" = "forestgreen", "P" = "firebrick", "I" = "darkgoldenrod1", "equivocal" = "gray")

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
  geom_vline(xintercept = max(bins.df[bins.df$seqnames == "21", 6]), lty = "dashed", col = "gray80") +
  theme_bw() + 
  geom_line(data = data.plot[!duplicated(data.plot$bin) & !is.na(data.plot$RMean_A), ], aes(x = bin, y = RMean_A), col = "darkblue", size = 1) +
  geom_line(data = data.plot[!duplicated(data.plot$bin) & !is.na(data.plot$RMean_N), ], aes(x = bin, y = RMean_N), col = "forestgreen", size = 1) +
  geom_line(data = data.plot[!duplicated(data.plot$bin) & !is.na(data.plot$RMean_P), ], aes(x = bin, y = RMean_P), col = "firebrick", size = 1) +
  geom_line(data = data.plot[!duplicated(data.plot$bin) & !is.na(data.plot$RMean_I), ], aes(x = bin, y = RMean_I), col = "darkgoldenrod1", size = 1) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_color_manual(values = col.subtypes)

pdf("Average_Methylation_100kBins_c2.pdf", width = 25, height = 10)
plot
dev.off()

