# This code is to plot the difference in the rank of biomarkers. 
# 1st section is for CCLE individual model vs GDSC-CCLE model. (V1)
# 2nd section is for GDSC individual model vs GDSC-CCLE model. (v2)
# The calculation is: V1 - V2. This value must be greater than 0. 
# Value > 0 indicated that the GDSC-CCLE model identifies signals better than
# individual models. 

# e - ccle exp
# h - ccle mutation (hybrid)
# se - sanger exp 
# sh - sanger mutation (hybrid)

# libraries 
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(cowplot)
library(ggrepel)
library(viridis)
library(plyr)
library(lattice)
library(hrbrthemes)

# loading the csv files with the ranks. 

FPATH1 <- "aml"
FPATH2 <- "pancancer"

## Common
drug <- c("X17AAG", "AZD0530", "AZD6244", "ERLOTINIB", "LAPATINIB", "NILOTINIB", 
          "NUTLIN3", "PD0332991", "PF2341066", "PHA665752", "PLX4720", "SORAFENIB")

## 1. CCLE - 

ccle_df <- data.frame(matrix(nrow = 0, ncol = 5))
longer_ccle <- data.frame(matrix(nrow = 0, ncol = 5))

for (i in drug) {
  
  #temp <- data.frame(matrix(nrow = 0, ncol = 5))
  #i = "X17AAG"
  
  path_ind <- paste("Spring'23-CancerSystemsBiology/Project/trained_models/", FPATH2, "/ccle Ranks/", sep="")
  path_dual <- paste("Spring'23-CancerSystemsBiology/Project/trained_models/", FPATH2, "/ccle Ranks/", sep = "")
  
  #fname_all <- paste(path_dual, i, "_cdf_e_h.csv", sep = "")
  fname_ind <- paste(path_ind,  i, "_cdf_e_h.csv", sep = "")
  fname_dual <- paste(path_dual, i, "_cdf_e_h_c.csv", sep = "")
  #fname_triple <- paste(path_dual, i, "_cdf_e_h_mgidb.csv", sep = "")
  
  
  #ccle_all <- read_csv(fname_all, col_names = TRUE)
  #colnames(ccle_all) <- c("gene", "e_p", "h_p")
  #ccle_all$max_p <- apply(ccle_all[c(2, 3)], 1, FUN = max, na.rm = TRUE)
  #ccle_all$max_p <- round(ccle_all$max_p, digits = 2)
  #ccle_all <- ccle_all[, c(1, 4)]
  #ccle_all["analysis"] <- c(rep("All", length(rownames(ccle_all)))) 
  #ccle_all["dr"] <- c(rep(i, length(rownames(ccle_all))))
  

  ccle_ind <- read_csv(fname_ind, col_names = TRUE)
  colnames(ccle_ind) <- c("gene", "e_p", "h_p")
  
  ccle_ind$max_p <- apply(ccle_ind[c(2, 3)], 1, FUN = max, na.rm = TRUE)
  ccle_ind$max_p <- round(ccle_ind$max_p, digits = 2)
  ccle_ind <- ccle_ind[, c(1, 4)]
  
  ccle_ind["analysis"] <- c(rep("e,h", length(rownames(ccle_ind)))) 
  ccle_ind["dr"] <- c(rep(i, length(rownames(ccle_ind))))
  
  ccle_dual <- read_csv(fname_dual, col_names = TRUE)
  #colnames(ccle_dual) <- c("rn", "max_score", "abs_max_score", "ccle_dual") #"max_score", "abs_max_score",
  colnames(ccle_dual) <- c("gene", "e_p", "h_p", "c_p")
  
  ccle_dual$max_p <- apply(ccle_dual[c(2, 3, 4)], 1, FUN = max, na.rm = TRUE)
  ccle_dual$max_p <- round(ccle_dual$max_p, digits = 2)
  ccle_dual <- ccle_dual[, c(1, 5)]
  
  ccle_dual["analysis"] <- c(rep("e,h,c", length(rownames(ccle_dual)))) 
  ccle_dual["dr"] <- c(rep(i, length(rownames(ccle_dual))))
  
  #ccle_triple <- read_csv(fname_triple, col_names = TRUE)
  #colnames(ccle_triple) <- c("gene", "e_p", "h_p")
  #ccle_triple$max_p <- apply(ccle_triple[c(2, 3)], 1, FUN = max, na.rm = TRUE)
  #ccle_triple$max_p <- round(ccle_triple$max_p, digits = 2)
  #ccle_triple <- ccle_triple[, c(1, 4)]
  #ccle_triple["analysis"] <- c(rep("MSigDB", length(rownames(ccle_triple)))) 
  #ccle_triple["dr"] <- c(rep(i, length(rownames(ccle_triple))))
  
  # p-values
  ccle_join <- left_join(ccle_ind, ccle_dual, by="gene")
  
  ccle_join["drug"] <- c(rep(i, length(rownames(ccle_join)))) 
  ccle_df <- rbind(ccle_df, ccle_join)
  
  longer_ccle <- rbind(longer_ccle, ccle_ind, ccle_dual)
}

##### molecular features comparison #########

## 2D scatter plot - with line
png(paste("Spring'23-CancerSystemsBiology/Project/trained_models/", FPATH2, "/ccle Ranks/plots/CCLE-2dScatter_ehc.jpg", sep=""))
ggplot(ccle_df, aes(x=max_p.x, y=max_p.y, color=drug, label=gene)) +
  geom_point() + #aes(size=delta)) +
  #ylim(-1000, 2000) +
  #geom_text(label= ccle_df$rn) +
  #coord_flip() +
  xlab("Probability of biomarkers (e,h)") +
  ylab("Probability of biomarkers (e,h,c)") +
  #geom_hline(yintercept=0, linetype="dashed", color = "black") +
  #geom_text_repel(max.overlaps = 100) +
  geom_abline(intercept=0, slope=1) +
  #geom_smooth(method=lm, se=FALSE, fullrange=TRUE, formula = y~x-1) +
  #theme(legend.position = "none") + 
  #facet_wrap(~drug) +
  ggtitle(paste("CCLE: Pancancer - effect of curation", sep = ""))
dev.off()
  
  ############################### p-values #####################################
  
png(paste("Spring'23-CancerSystemsBiology/Project/trained_models/", 
          FPATH2, "/ccle Ranks/plots/CCLE-boxplot_ehc.jpg", sep=""), width = 700, height = 700)
ggplot(longer_ccle, aes(dr, max_p, fill=analysis)) +
  #geom_violin() +
  geom_boxplot() +
  #geom_point() +
  xlab("Drugs") +
  ylab("Probability") +
  guides(fill=guide_legend(title="Features")) +
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  scale_fill_manual(values = c("#00BFC4",  "#C77CFF")) + # "#F8766D",  "#7CAE00"
  ggtitle(paste("CCLE: Pancancer - Drugs vs Probability ", sep="")) 
dev.off()

################################# 2. GDSC ######################################

## 2. uncommon GDSC - 

#drugs <- c("X17AAG", "AZD6244", "NILOTINIB", "NUTLIN3", "PD0325901", "PD0332991", "PLX4720")

gdsc_df <- data.frame(matrix(nrow = 0, ncol = 5))
longer_gdsc <- data.frame(matrix(nrow = 0, ncol = 5))

for (i in drug) {
  
  #i <- "X17AAG"
  
  path_ind <- paste("Spring'23-CancerSystemsBiology/Project/trained_models/", FPATH2, "/gdsc Ranks/", sep="")
  path_dual <- paste("Spring'23-CancerSystemsBiology/Project/trained_models/", FPATH2, "/gdsc Ranks/", sep = "")
  
  #fname_all <- paste(path_ind,  i, "_cdf_e_h.csv", sep = "")
  fname_ind <- paste(path_ind,  i, "_cdf_e_h.csv", sep = "")
  fname_dual <- paste(path_dual, i, "_cdf_e_h_c.csv", sep = "")
  #fname_triple <- paste(path_dual, i, "_cdf_e_h_mgidb.csv", sep = "")
  
  #gdsc_all <- read_csv(fname_all, col_names = TRUE)
  #colnames(gdsc_all) <- c("gene", "e_p", "h_p")
  #gdsc_all$max_p <- apply(gdsc_all[c(2, 3)], 1, FUN = max, na.rm = TRUE)
  #gdsc_all$max_p <- round(gdsc_all$max_p, digits = 2)
  #gdsc_all <- gdsc_all[, c(1, 4)]
  #gdsc_all["analysis"] <- c(rep("All", length(rownames(gdsc_all)))) 
  #gdsc_all["dr"] <- c(rep(i, length(rownames(gdsc_all))))
  
  
  gdsc_ind <- read_csv(fname_ind, col_names = TRUE)
  colnames(gdsc_ind) <- c("gene", "e_p", "h_p")
  
  gdsc_ind$max_p <- apply(gdsc_ind[c(2, 3)], 1, FUN = max, na.rm = TRUE)
  gdsc_ind$max_p <- round(gdsc_ind$max_p, digits = 2)
  gdsc_ind <- gdsc_ind[, c(1, 4)]
  
  gdsc_ind["analysis"] <- c(rep("e,h", length(rownames(gdsc_ind)))) 
  gdsc_ind["dr"] <- c(rep(i, length(rownames(gdsc_ind))))
  
  
  gdsc_dual <- read_csv(fname_dual, col_names = TRUE)
  colnames(gdsc_dual) <- c("gene", "e_p", "h_p", "c_p")
  
  gdsc_dual$max_p <- apply(gdsc_dual[c(2, 3, 4)], 1, FUN = max, na.rm = TRUE)
  gdsc_dual$max_p <- round(gdsc_dual$max_p, digits = 2)
  gdsc_dual <- gdsc_dual[, c(1, 5)]
  
  gdsc_dual["analysis"] <- c(rep("e,h,c", length(rownames(gdsc_dual)))) 
  gdsc_dual["dr"] <- c(rep(i, length(rownames(gdsc_dual))))
  
  #gdsc_triple <- read_csv(fname_triple, col_names = TRUE)
  #colnames(gdsc_triple) <- c("gene", "e_p", "h_p")
  #gdsc_triple$max_p <- apply(gdsc_triple[c(2, 3)], 1, FUN = max, na.rm = TRUE)
  #gdsc_triple$max_p <- round(gdsc_triple$max_p, digits = 2)
  #gdsc_triple <- gdsc_triple[, c(1, 4)]
  #gdsc_triple["analysis"] <- c(rep("MSigDB", length(rownames(gdsc_triple)))) 
  #gdsc_triple["dr"] <- c(rep(i, length(rownames(gdsc_triple))))
  
  gdsc_join <- left_join(gdsc_ind, gdsc_dual, by="gene")
  gdsc_join["drug"] <- c(rep(i, length(rownames(gdsc_join))))
  gdsc_df <- rbind(gdsc_df, gdsc_join)
  
  longer_gdsc <- rbind(longer_gdsc, gdsc_ind, gdsc_dual)
}

## 2D scatter plot - with line
png(paste("Spring'23-CancerSystemsBiology/Project/trained_models/", FPATH2, "/gdsc Ranks/plots/GDSC-2dScatter_ehc.jpg", sep=""))
ggplot(gdsc_df, aes(x=max_p.x, y=max_p.y, color=drug, label=gene)) +
  geom_point() + #aes(size=delta)) +
  #coord_flip() +
  xlab("Probability of biomarkers (e,h)") +
  ylab("Probability of biomarkers (e,h,c)") +
  #geom_hline(yintercept=0, linetype="dashed", color = "black") +
  #geom_text_repel(max.overlaps = 100) +
  geom_abline(intercept=0, slope=1) +
  #geom_smooth(method=lm, se=FALSE, fullrange=TRUE, formula = y~x-1) +
  #theme(legend.position = "none") + 
  #facet_wrap(~drug) +
  ggtitle(paste("GDSC: Pancancer - effect of curation", sep = ""))
dev.off()
  
  
############################### p-values #####################################
  
png(paste("Spring'23-CancerSystemsBiology/Project/trained_models/", 
          FPATH2, "/gdsc Ranks/plots/GDSC-boxplot_ehc.jpg", sep=""), width = 700, height = 700)
ggplot(longer_gdsc, aes(dr, max_p, fill=analysis)) +
  #geom_violin() +
  geom_boxplot() +
  #geom_point() +  
  xlab("Drugs") +
  ylab("Probability") +
  guides(fill=guide_legend(title="Features")) +
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  scale_fill_manual(values = c("#00BFC4",  "#C77CFF")) + # "#F8766D",  "#7CAE00"
  ggtitle(paste("GDSC: Pancancer - Drugs vs Probability", sep="")) 
dev.off()

####### pan cancer vs aml biomarkers comparison #############

comp_cancer <- data.frame(matrix(nrow = 0, ncol = 4))

for (i in drug) {
  
  #temp <- data.frame(matrix(nrow = 0, ncol = 5))
  #i = "X17AAG"
  
  path_aml_ccle <- paste("Spring'23-CancerSystemsBiology/Project/trained_models/", FPATH1, "/ccle Ranks/", sep="")
  path_pan_ccle <- paste("Spring'23-CancerSystemsBiology/Project/trained_models/", FPATH2, "/ccle Ranks/", sep = "")
  path_aml_gdsc <- paste("Spring'23-CancerSystemsBiology/Project/trained_models/", FPATH1, "/gdsc Ranks/", sep="")
  path_pan_gdsc <- paste("Spring'23-CancerSystemsBiology/Project/trained_models/", FPATH2, "/gdsc Ranks/", sep = "")
  path_baml <- paste("Spring'23-CancerSystemsBiology/Project/trained_models/", FPATH1, "/baml Ranks/", sep="")
  
  fname_aml_ccle <- paste(path_aml_ccle,  i, "_cdf_e_h_mgidb.csv", sep = "")
  fname_pan_ccle <- paste(path_pan_ccle, i, "_cdf_e_h_mgidb.csv", sep = "")
  fname_aml_gdsc <- paste(path_aml_gdsc,  i, "_cdf_e_h_mgidb.csv", sep = "")
  fname_pan_gdsc <- paste(path_pan_gdsc, i, "_cdf_e_h_mgidb.csv", sep = "")
  fname_baml <- paste(path_baml,  i, "_cdf_e_mgidb.csv", sep = "")
  
  
  ccle_aml <- read_csv(fname_aml_ccle, col_names = TRUE)
  colnames(ccle_aml) <- c("gene", "e_p", "h_p")
  
  ccle_aml$max_p <- apply(ccle_aml[c(2, 3)], 1, FUN = max, na.rm = TRUE)
  ccle_aml$max_p <- round(ccle_aml$max_p, digits = 2)
  ccle_aml <- ccle_aml[, c(1, 4)]
  ccle_aml_fin <- ccle_aml %>% filter(max_p >= 0.6)
  
  # rank 
  #ccle_aml <- ccle_aml[, c(1, 4)]
  #ccle_aml_fin <- ccle_aml #%>% filter(rank <= 10000)
  
  ccle_aml_fin["analysis"] <- c(rep("AML", length(rownames(ccle_aml_fin)))) 
  ccle_aml_fin["dr"] <- c(rep(i, length(rownames(ccle_aml_fin))))
  
  ##############################################################################
  
  ccle_pan <- read_csv(fname_pan_ccle, col_names = TRUE)
  colnames(ccle_pan) <- c("gene", "e_p", "h_p")
  
  ccle_pan$max_p <- apply(ccle_pan[c(2, 3)], 1, FUN = max, na.rm = TRUE)
  ccle_pan$max_p <- round(ccle_pan$max_p, digits = 2)
  ccle_pan <- ccle_pan[, c(1, 4)]
  ccle_pan_fin <- ccle_pan %>% filter(max_p >= 0.6)
  
  #rank
  #ccle_pan <- ccle_pan[, c(1, 4)]
  #ccle_pan_fin <- ccle_pan #%>% filter(max_p >= 0.6)
  
  ccle_pan_fin["analysis"] <- c(rep("Pancancer", length(rownames(ccle_pan_fin)))) 
  ccle_pan_fin["dr"] <- c(rep(i, length(rownames(ccle_pan_fin))))
  
  ##############################################################################
  
  gdsc_aml <- read_csv(fname_aml_gdsc, col_names = TRUE)
  colnames(gdsc_aml) <- c("gene", "e_p", "h_p")
  
  gdsc_aml$max_p <- apply(gdsc_aml[c(2, 3)], 1, FUN = max, na.rm = TRUE)
  gdsc_aml$max_p <- round(gdsc_aml$max_p, digits = 2)
  gdsc_aml <- gdsc_aml[, c(1, 4)]
  gdsc_aml_fin <- gdsc_aml %>% filter(max_p >= 0.6)
  
  #rank
  #gdsc_aml <- gdsc_aml[, c(1, 4)]
  #gdsc_aml_fin <- gdsc_aml #%>% filter(max_p >= 0.6)
  
  gdsc_aml_fin["analysis"] <- c(rep("AML", length(rownames(gdsc_aml_fin)))) 
  gdsc_aml_fin["dr"] <- c(rep(i, length(rownames(gdsc_aml_fin))))
  
  ##############################################################################
  
  gdsc_pan <- read_csv(fname_pan_gdsc, col_names = TRUE)
  colnames(gdsc_pan) <- c("gene", "e_p", "h_p")
  
  gdsc_pan$max_p <- apply(gdsc_pan[c(2, 3)], 1, FUN = max, na.rm = TRUE)
  gdsc_pan$max_p <- round(gdsc_pan$max_p, digits = 2)
  gdsc_pan <- gdsc_pan[, c(1, 4)]
  gdsc_pan_fin <- gdsc_pan %>% filter(max_p >= 0.6)
  
  #rank
  #gdsc_pan <- gdsc_pan[, c(1, 4)]
  #gdsc_pan_fin <- gdsc_pan #%>% filter(max_p >= 0.6)
  
  gdsc_pan_fin["analysis"] <- c(rep("Pancancer", length(rownames(gdsc_pan_fin)))) 
  gdsc_pan_fin["dr"] <- c(rep(i, length(rownames(gdsc_pan_fin))))
  
  ##############################################################################

  baml <- read_csv(fname_baml, col_names = TRUE)
  colnames(baml) <- c("gene", "max_p")
  baml$max_p <- round(baml$max_p, digits = 2)
  baml_fin <- baml %>% filter(max_p >= 0.6)
  
  #rank
  #baml <- baml[, c(1, 4)]
  #baml_fin <- baml #%>% filter(max_p >= 0.6)
  
  baml_fin["analysis"] <- c(rep("AML", length(rownames(baml_fin)))) 
  baml_fin["dr"] <- c(rep(i, length(rownames(baml_fin))))
  
  
  # p-values
  #ccle_join <- left_join(ccle_ind, ccle_dual, by="gene")
  #ccle_join["drug"] <- c(rep(i, length(rownames(ccle_join)))) 
  #ccle_df <- rbind(ccle_df, ccle_join)
  
  comp_cancer <- rbind(comp_cancer, ccle_aml_fin, gdsc_aml_fin, ccle_pan_fin, gdsc_pan_fin, baml_fin)
}

#aml <- ccle_df %>% filter(max_p.x >= 0.5 & drug == "X17AAG")

#png(paste("Spring'23-CancerSystemsBiology/Project/trained_models/", FPATH, "/gdsc Ranks/plots/GDSC-boxplot.jpg", sep=""))
#ggplot(aml, aes(gene, max_p.x)) + 
#  geom_bar(stat = "identity", fill="purple") +
#  geom_text(aes(label=max_p.x), vjust=1.6, color="white", size=3.5) +
#  xlab("Biomarkers") +
#  ylab("Probability") +
#  ggtitle(paste("AML biomarkers for drug ", drug, sep=""))

#dotplot(max_p~gene|analysis, data = (comp_cancer %>% filter(dr == "X17AAG")))

#data = (comp_cancer %>% filter(dr == "X17AAG"))

#dotchart(data$max_p, labels = data$gene, pch = 21, bg = "green", pt.cex = 1.5)#, groups = rev(data$analysis))

#ggplot(data) +
#  geom_point(aes(gene, analysis, size=max_p), color="darkred") +
#  geom_point(aes(gene, analysis, size=max_p), color="black",pch=5) + #, color=max_p)) +
#  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 

data_filter_analysis = (comp_cancer %>% filter(analysis == "AML"))
aggr_aml <- data_filter_analysis %>% group_by(gene, dr) %>% summarise_at(vars(max_p), list(max_p = mean))
#aggr_aml$rank <- round(aggr_aml$rank) 

aml_geneset <- aggr_aml %>% group_by(gene) %>% tally()
aml_geneset <- aml_geneset %>% filter(n > 1)

df_aml <- aggr_aml[aggr_aml$gene %in% aml_geneset$gene,]

gene_order_aml <- aml_geneset[order(aml_geneset$n),]

png(paste("Spring'23-CancerSystemsBiology/Project/trained_models/aml_heatmap_cdf_mgidb.jpg", sep=""), width = 800, height =800)
ggplot(df_aml, aes(dr, gene, fill= max_p)) + 
  geom_tile() +
  scale_fill_gradient(low="red", high="green") +
  theme_ipsum() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position = "none") +
  scale_x_discrete(limits=c("SORAFENIB","ERLOTINIB","LAPATINIB", "AZD6244", 
                            "PD0332991", "AZD0530", "X17AAG", "PF2341066",
                            "PLX4720", "NILOTINIB", "NUTLIN3")) +
  scale_y_discrete(limits=c(gene_order_aml$gene)) +
  #guides(fill=guide_legend(title="Probability", reverse = TRUE)) +
  #theme(legend.position = "none")     
  xlab("Drugs") +
  ylab("Biomarkers/Targets") +
  ggtitle("AML")
dev.off()

# pancancer
data_filter_analysis_pan = (comp_cancer %>% filter(analysis == "Pancancer"))
aggr_pan <- data_filter_analysis_pan %>% group_by(gene, dr) %>% summarise_at(vars(max_p), list(max_p = mean))
#aggr_pan$rank <- round(aggr_pan$rank) 

pan_geneset <- aggr_pan %>% group_by(gene) %>% tally()
pan_geneset <- pan_geneset %>% filter(n > 1)

df_pan <- aggr_pan[aggr_pan$gene %in% pan_geneset$gene,]

gene_order_pan <- pan_geneset[order(pan_geneset$n),]

png(paste("Spring'23-CancerSystemsBiology/Project/trained_models/pancancer_heatmap_cdf_mgidb.jpg", sep=""))
ggplot(df_pan, aes(dr, gene, fill= max_p)) + 
  geom_tile() +
  #scale_fill_distiller(palette = "RdPu") +
  scale_fill_gradient(low="red", high="green") +
  theme_ipsum() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_x_discrete(limits=c("SORAFENIB","ERLOTINIB","LAPATINIB", "AZD6244", 
                            "PD0332991", "AZD0530", "X17AAG", "PF2341066",
                            "PLX4720", "NILOTINIB", "NUTLIN3")) +
  scale_y_discrete(limits=c(gene_order_pan$gene)) +
  guides(fill=guide_legend(title="Probability", reverse = TRUE)) +
  xlab("Drugs") +
  ylab("Biomarkers/Targets") +
  ggtitle("Pancancer")
dev.off()


