
#########################################################################################################################################
### This script allows the clinical analysis of survival differences between RNA-based and DNA-based classification as in Figure 4C/D ###
#########################################################################################################################################

#load in packages
library(survival)
library(survminer)
library(tidyverse)


####################################
#### Run clinical data analysis ####
####################################

#compare results from RNAseq and Methylation

#load in data
data <- #load in data with clinical information and subtyping information 

fit.A <- survfit(Surv(OS , Deceased_State) ~ Class, data = data[data$Subtype == "SCLC-A", ])
fit.N <- survfit(Surv(OS , Deceased_State) ~ Class, data = data[data$Subtype == "SCLC-N", ])

#create median OS
med_OS.A <- surv_median(fit.A)
med_OS.A[is.na(med_OS.A)] <- "NR"

med_OS.N <- surv_median(fit.N)
med_OS.N[is.na(med_OS.N)] <- "NR"


#create cox proprotional hazard ratio 
cox.A <- coxph(Surv(OS , Deceased_State) ~ Class, data = data[data$Subtype == "SCLC-A", ])
cox.N <- coxph(Surv(OS , Deceased_State) ~ Class, data = data[data$Subtype == "SCLC-N", ])

ggsplot <- ggsurvplot(fit.A, data = data, size = 1, palette = c("blue4", "darkgoldenrod1"),
                      pval = FALSE, risk.table = TRUE, ggtheme = theme_bw(), surv.median.line = "hv")

ggsplot$plot <- ggsplot$plot + ggplot2::annotate("text", x = Inf, y = 0.85, vjust = 1, hjust = 1,
                                                 label = paste0("HR = ", round(summary(cox.A)$conf.int[1], 2), 
                                                                " (", round(summary(cox.A)$conf.int[3], 2), " - ", 
                                                                round(summary(cox.A)$conf.int[4], 2), ")")) +
  ggplot2::annotate("text", x = Inf, y = 0.95, vjust = 1, hjust = 1,
                    label = paste0("median OS (95% CI) = ", 
                                   round(med_OS.A$median[1], 1), " (",
                                   round(med_OS.A$lower[1], 1), " - ",
                                   med_OS.A$upper[1], ") vs ",
                                   round(med_OS.A$median[2], 1), " (",
                                   round(med_OS.A$lower[2], 1), " - ",
                                   med_OS.A$upper[2], "); ",
                                   surv_pvalue(fit.A)$pval.txt)) +
  ggplot2::xlab("Months")

pdf("Differenc_Classification_A_KM.pdf", width = 5, height = 5)
ggsplot
dev.off()


ggsplot <- ggsurvplot(fit.N, data = data, size = 1, palette = c("blue4", "darkgoldenrod1"),
                      pval = FALSE, risk.table = TRUE, ggtheme = theme_bw(), surv.median.line = "hv")

ggsplot$plot <- ggsplot$plot + ggplot2::annotate("text", x = Inf, y = 0.85, vjust = 1, hjust = 1,
                                                 label = paste0("HR = ", round(summary(cox.N)$conf.int[1], 2), 
                                                                " (", round(summary(cox.N)$conf.int[3], 2), " - ", 
                                                                round(summary(cox.N)$conf.int[4], 2), ")")) +
  ggplot2::annotate("text", x = Inf, y = 0.95, vjust = 1, hjust = 1,
                    label = paste0("median OS (95% CI) = ", 
                                   round(med_OS.N$median[1], 1), " (",
                                   round(med_OS.N$lower[1], 1), " - ",
                                   med_OS.N$upper[1], ") vs ",
                                   round(med_OS.N$median[2], 1), " (",
                                   round(med_OS.N$lower[2], 1), " - ",
                                   med_OS.N$upper[2], "); ",
                                   surv_pvalue(fit.N)$pval.txt)) +
  ggplot2::xlab("Months")

pdf("Differenc_Classification_N_KM.pdf", width = 5, height = 5)
ggsplot
dev.off()