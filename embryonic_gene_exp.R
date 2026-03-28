library('biomaRt')
library('ggplot2')
library('tidyr')
library('dplyr')
library('data.table')
library('DESeq2')
library('ggpubr')
library('readxl')

xol1sex1_N2_ee_genes <- read.table("salmon_analysis/xol1_sex1/xol1sex1_N2_ee_genes.csv")
xol1sex1_xol1_ee_genes <- read.table("salmon_analysis/xol1_sex1/xol1sex1_xol1_ee_genes.csv")
xol1_N2_ee_genes <- read.table("salmon_analysis/xol1_sex1/xol1_N2_ee_genes.csv")

EE_genes <- read.table("Gene Sets/EE_EGs.txt", header = TRUE)
LE_genes <- read.table("Gene Sets/LE_EGs.txt", header = TRUE)
L1_genes <- read.table("Gene Sets/L1_EGs.txt", header = TRUE)
L2_genes <- read.table("Gene Sets/L2_EGs.txt", header = TRUE)
L3_genes <- read.table("Gene Sets/L3_EGs.txt", header = TRUE)
L4_genes <- read.table("Eshna/Gene Sets/L4_EGs.txt", header = TRUE)
YA_gonad <- read.table("Gene Sets/YA-gonad_EGs.txt", header = TRUE)
YA_genes <- read.table("/Gene Sets/YA_EGs.txt", header = TRUE)

common_embryo_genes <- EE_genes$gene_id[EE_genes$gene_id %in% LE_genes$gene_id] 
#5243
common_embryo_larval_genes <- EE_genes$gene_id[EE_genes$gene_id %in% LE_genes$gene_id &
                                                 EE_genes$gene_id %in% L1_genes$gene_id &
                                                 EE_genes$gene_id %in% L2_genes$gene_id &
                                                 EE_genes$gene_id %in% L3_genes$gene_id &
                                                 EE_genes$gene_id %in% L4_genes$gene_id]
#4368 genes
common_larval_genes <- L1_genes$gene_id[L1_genes$gene_id %in% L2_genes$gene_id &
                                          L1_genes$gene_id %in% L3_genes$gene_id &
                                          L1_genes$gene_id %in% L4_genes$gene_id]
#4807 genes
common_all_stages <- EE_genes$gene_id[EE_genes$gene_id %in% common_embryo_genes &
                                        EE_genes$gene_id %in% common_embryo_larval_genes &
                                        EE_genes$gene_id %in% common_larval_genes &
                                        EE_genes$gene_id %in% YA_genes]
#none
LE_genes_core <- LE_genes[!LE_genes$gene_id %in% common_embryo_larval_genes &
                            !LE_genes$gene_id %in% common_embryo_genes,]
EE_genes_core <- EE_genes[!EE_genes$gene_id %in% common_embryo_larval_genes &
                            !EE_genes$gene_id %in% common_embryo_genes,]

#### late embryonic genes ##

colnames(xol1sex1_N2_ee_genes)[1] <- "EnsemblID"
xol1sex1_N2_ee_genes <- drop_na(xol1sex1_N2_ee_genes)  

xol1sex1_N2_ee_legenes <- xol1sex1_N2_ee_genes[xol1sex1_N2_ee_genes$EnsemblID %in% LE_genes_core$gene_id,]
xol1sex1_N2_ee_legenes$genotype <- rep("late embryo enriched genes", nrow(xol1sex1_N2_ee_legenes))
xol1sex1_N2_ee_notlegenes <- xol1sex1_N2_ee_genes[!xol1sex1_N2_ee_genes$EnsemblID %in% LE_genes_core$gene_id,]
xol1sex1_N2_ee_notlegenes$genotype <- rep("not late embryo enriched genes", nrow(xol1sex1_N2_ee_notlegenes))
xol1sex1_N2_ee_legenes_labeled <- rbind(xol1sex1_N2_ee_legenes, xol1sex1_N2_ee_notlegenes)

plot_LE_xol1sex1_WT <- ggplot(xol1sex1_N2_ee_legenes_labeled, aes(x = as.factor(genotype) ,y = log2FoldChange)) + geom_boxplot(outlier.shape = NA) + 
  coord_cartesian(ylim = c(-2,2)) +
  geom_hline(linetype = 2, alpha = 0.3, yintercept = median(xol1sex1_N2_ee_legenes_labeled$log2FoldChange[xol1sex1_N2_ee_legenes_labeled$genotype == "not late embryo enriched genes"])) +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(size = 15), axis.title.y = element_text(size = 18)) +
  labs(title = "xol-1 sex-1 / WT")
plot_LE_xol1sex1_WT

wilcox.test(xol1sex1_N2_ee_legenes_labeled$log2FoldChange[xol1sex1_N2_ee_legenes_labeled$genotype == "late embryo enriched genes"],
            xol1sex1_N2_ee_legenes_labeled$log2FoldChange[xol1sex1_N2_ee_legenes_labeled$genotype == "not late embryo enriched genes"])
#p-value = 0.1679

###moods median test
moodtest_legenes_xol1sex1_wt_legenes <- xol1sex1_N2_ee_legenes_labeled[xol1sex1_N2_ee_legenes_labeled$genotype == "late embryo enriched genes", c("log2FoldChange", "genotype")]
moodtest_legenes_xol1sex1_wt_notlegenes <- xol1sex1_N2_ee_legenes_labeled[xol1sex1_N2_ee_legenes_labeled$genotype == "not late embryo enriched genes", c("log2FoldChange", "genotype")]
moodtest_legenes_xol1sex1_wt <- rbind(moodtest_legenes_xol1sex1_wt_legenes,
                                      moodtest_legenes_xol1sex1_wt_notlegenes)
mood.medtest(log2FoldChange ~ genotype,
             data = moodtest_legenes_xol1sex1_wt)
#pvalue=1

ggsave(plot_LE_xol1sex1_WT, filename = "boxplot_xs_wt_LEgenes.png" , path = "salmon_analysis/xol1_sex1")

##
colnames(xol1sex1_xol1_ee_genes)[1] <- "EnsemblID"
xol1sex1_xol1_ee_genes <- drop_na(xol1sex1_xol1_ee_genes)

xol1sex1_xol1_ee_legenes <- xol1sex1_xol1_ee_genes[xol1sex1_xol1_ee_genes$EnsemblID %in% LE_genes_core$gene_id,]
xol1sex1_xol1_ee_legenes$genotype <- rep("late embryo enriched genes", nrow(xol1sex1_xol1_ee_legenes))
xol1sex1_xol1_ee_notlegenes <- xol1sex1_xol1_ee_genes[!xol1sex1_xol1_ee_genes$EnsemblID %in% LE_genes_core$gene_id,]
xol1sex1_xol1_ee_notlegenes$genotype <- rep("not late embryo enriched genes", nrow(xol1sex1_xol1_ee_notlegenes))
xol1sex1_xol1_ee_legenes_labeled <- rbind(xol1sex1_xol1_ee_legenes, xol1sex1_xol1_ee_notlegenes)

plot_LE_xol1sex1_xol1 <- ggplot(xol1sex1_xol1_ee_legenes_labeled, aes(x = as.factor(genotype) ,y = log2FoldChange)) + geom_boxplot(outlier.shape = NA) + 
  coord_cartesian(ylim = c(-2.3,2)) +
  geom_hline(linetype = 2, alpha = 0.3, yintercept = median(xol1sex1_xol1_ee_legenes_labeled$log2FoldChange[xol1sex1_xol1_ee_legenes_labeled$genotype == "not late embryo enriched genes"])) +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(size = 15), axis.title.y = element_text(size = 18)) +
  labs(title = "xol-1 sex-1/ xol-1")
plot_LE_xol1sex1_xol1

wilcox.test(xol1sex1_xol1_ee_legenes_labeled$log2FoldChange[xol1sex1_xol1_ee_legenes_labeled$genotype == "late embryo enriched genes"],
            xol1sex1_xol1_ee_legenes_labeled$log2FoldChange[xol1sex1_xol1_ee_legenes_labeled$genotype == "not late embryo enriched genes"])
#p-value = 1.391e-08

###moods median test
moodtest_legenes_xol1sex1_xol1_legenes <- xol1sex1_xol1_ee_legenes_labeled[xol1sex1_xol1_ee_legenes_labeled$genotype == "late embryo enriched genes", c("log2FoldChange", "genotype")]
moodtest_legenes_xol1sex1_xol1_notlegenes <- xol1sex1_xol1_ee_legenes_labeled[xol1sex1_xol1_ee_legenes_labeled$genotype == "not late embryo enriched genes", c("log2FoldChange", "genotype")]
moodtest_legenes_xol1sex1_xol1 <- rbind(moodtest_legenes_xol1sex1_xol1_legenes,
                                        moodtest_legenes_xol1sex1_xol1_notlegenes)
mood.medtest(log2FoldChange ~ genotype,
             data = moodtest_legenes_xol1sex1_xol1)
#p-value = 3.567e-06

ggsave(plot_LE_xol1sex1_xol1, filename = "boxplot_xs_x_LEgenes.png" , path = "salmon_analysis/xol1_sex1")

##
colnames(xol1_N2_ee_genes)[1] <- "EnsemblID"
xol1_N2_ee_genes <- drop_na(xol1_N2_ee_genes)

xol1_N2_ee_legenes <- xol1_N2_ee_genes[xol1_N2_ee_genes$EnsemblID %in% LE_genes_core$gene_id,]
xol1_N2_ee_legenes$genotype <- rep("late embryo enriched genes", nrow(xol1_N2_ee_legenes))
xol1_N2_ee_notlegenes <- xol1_N2_ee_genes[!xol1_N2_ee_genes$EnsemblID %in% LE_genes_core$gene_id,]
xol1_N2_ee_notlegenes$genotype <- rep("not late embryo enriched genes", nrow(xol1_N2_ee_notlegenes))
xol1_N2_ee_legenes_labeled <- rbind(xol1_N2_ee_legenes, xol1_N2_ee_notlegenes)

plot_LE_xol1_N2 <- ggplot(xol1_N2_ee_legenes_labeled, aes(x = as.factor(genotype) ,y = log2FoldChange)) + geom_boxplot(outlier.shape = NA) + 
  coord_cartesian(ylim = c(-2.5,2.5)) +
  geom_hline(linetype = 2, alpha = 0.3, yintercept = median(xol1_N2_ee_legenes_labeled$log2FoldChange[xol1_N2_ee_legenes_labeled$genotype == "not late embryo enriched genes"])) +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(size = 15), axis.title.y = element_text(size = 18)) +
  labs(title = "xol-1 sex-1/ xol-1")
plot_LE_xol1_N2

wilcox.test(xol1_N2_ee_legenes_labeled$log2FoldChange[xol1_N2_ee_legenes_labeled$genotype == "late embryo enriched genes"],
            xol1_N2_ee_legenes_labeled$log2FoldChange[xol1_N2_ee_legenes_labeled$genotype == "not late embryo enriched genes"])
#p-value = 8.039e-08

###moods median test
moodtest_legenes_xol1_wt_legenes <- xol1_N2_ee_legenes_labeled[xol1_N2_ee_legenes_labeled$genotype == "late embryo enriched genes", c("log2FoldChange", "genotype")]
moodtest_legenes_xol1_wt_notlegenes <- xol1_N2_ee_legenes_labeled[xol1_N2_ee_legenes_labeled$genotype == "not late embryo enriched genes", c("log2FoldChange", "genotype")]
moodtest_legenes_xol1_wt <- rbind(moodtest_legenes_xol1_wt_legenes,
                                  moodtest_legenes_xol1_wt_notlegenes)
mood.medtest(log2FoldChange ~ genotype,
             data = moodtest_legenes_xol1_wt)
#p-value = 0.0005415

ggsave(plot_LE_xol1_N2, filename = "boxplot_x_wt_LEgenes.png" , path = "salmon_analysis/xol1_sex1")

##combined

xol1_N2_ee_legenes_labeled$dataset <- rep("xol1/WT", nrow(xol1_N2_ee_legenes_labeled))
xol1sex1_N2_ee_legenes_labeled$dataset <- rep("xol1_sex1/WT", nrow(xol1sex1_N2_ee_legenes_labeled))
xol1sex1_xol1_ee_legenes_labeled$dataset <- rep("xol1_sex1/xol1", nrow(xol1sex1_xol1_ee_legenes_labeled))

combined_le <- rbind(xol1_N2_ee_legenes_labeled,
                     xol1sex1_N2_ee_legenes_labeled,
                     xol1sex1_xol1_ee_legenes_labeled)
combined_le$dataset <- factor(combined_le$dataset, levels = c("xol1/WT",
                                                              "xol1_sex1/WT",
                                                              "xol1_sex1/xol1"))

ggplot(combined_le, aes(x = as.factor(genotype), y = log2FoldChange, fill = genotype)) + 
  stat_boxplot(geom = "errorbar", width = 0.25) +  
  geom_boxplot(outlier.shape = NA) + 
  facet_wrap(~ dataset) +
  coord_cartesian(ylim = c(-2.5,2.5)) +
  geom_hline(linetype = 2, alpha = 0.3, yintercept = 0) +
  scale_fill_manual(values = c("#008080","lightgrey")) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 11, color = "black"), 
        axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 18), 
        title = element_text(size = 8),
        legend.position = "none") +
  theme(axis.title.y = element_blank()) +
  scale_x_discrete(labels = c("Late-embryonic \n Genes", " All Other \n Genes"))

ggsave(filename = "boxplot_lateemb_all.png",
       path = "sex-1 paper/Data/",
       height = 4,
       width = 8)

### early embryonic genes ###

xol1sex1_N2_ee_eegenes <- xol1sex1_N2_ee_genes[xol1sex1_N2_ee_genes$EnsemblID %in% EE_genes_core$gene_id,]
xol1sex1_N2_ee_eegenes$genotype <- rep("early embryo enriched genes", nrow(xol1sex1_N2_ee_eegenes))
xol1sex1_N2_ee_noteegenes <- xol1sex1_N2_ee_genes[!xol1sex1_N2_ee_genes$EnsemblID %in% EE_genes_core$gene_id,]
xol1sex1_N2_ee_noteegenes$genotype <- rep("not early embryo enriched genes", nrow(xol1sex1_N2_ee_noteegenes))
xol1sex1_N2_ee_eegenes_labeled <- rbind(xol1sex1_N2_ee_eegenes, xol1sex1_N2_ee_noteegenes)

plot_EE_xol1sex1_WT <- ggplot(xol1sex1_N2_ee_eegenes_labeled, aes(x = as.factor(genotype) ,y = log2FoldChange)) + geom_boxplot(outlier.shape = NA) + 
  coord_cartesian(ylim = c(-2,2)) +
  geom_hline(linetype = 2, alpha = 0.3, yintercept = median(xol1sex1_N2_ee_eegenes_labeled$log2FoldChange[xol1sex1_N2_ee_eegenes_labeled$genotype == "not early embryo enriched genes"])) +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(size = 15), axis.title.y = element_text(size = 18)) +
  labs(title = "xol-1 sex-1 / WT")
plot_EE_xol1sex1_WT

wilcox.test(xol1sex1_N2_ee_eegenes_labeled$log2FoldChange[xol1sex1_N2_ee_eegenes_labeled$genotype == "early embryo enriched genes"],
            xol1sex1_N2_ee_eegenes_labeled$log2FoldChange[xol1sex1_N2_ee_eegenes_labeled$genotype == "not early embryo enriched genes"])
#p-value = 0.6111

###moods median test
moodtest_eegenes_xol1sex1_wt_eegenes <- xol1sex1_N2_ee_eegenes_labeled[xol1sex1_N2_ee_eegenes_labeled$genotype == "early embryo enriched genes", c("log2FoldChange", "genotype")]
moodtest_eegenes_xol1sex1_wt_noteegenes <- xol1sex1_N2_ee_eegenes_labeled[xol1sex1_N2_ee_eegenes_labeled$genotype == "not early embryo enriched genes", c("log2FoldChange", "genotype")]
moodtest_eegenes_xol1sex1_wt <- rbind(moodtest_eegenes_xol1sex1_wt_eegenes,
                                      moodtest_eegenes_xol1sex1_wt_noteegenes)
mood.medtest(log2FoldChange ~ genotype,
             data = moodtest_eegenes_xol1sex1_wt)
#p-value = 0.3558

ggsave(plot_EE_xol1sex1_WT, filename = "boxplot_xs_wt_EEgenes.png" , path = "salmon_analysis/xol1_sex1")

##

xol1sex1_xol1_ee_eegenes <- xol1sex1_xol1_ee_genes[xol1sex1_xol1_ee_genes$EnsemblID %in% EE_genes_core$gene_id,]
xol1sex1_xol1_ee_eegenes$genotype <- rep("early embryo enriched genes", nrow(xol1sex1_xol1_ee_eegenes))
xol1sex1_xol1_ee_noteegenes <- xol1sex1_xol1_ee_genes[!xol1sex1_xol1_ee_genes$EnsemblID %in% EE_genes_core$gene_id,]
xol1sex1_xol1_ee_noteegenes$genotype <- rep("not early embryo enriched genes", nrow(xol1sex1_xol1_ee_noteegenes))
xol1sex1_xol1_ee_eegenes_labeled <- rbind(xol1sex1_xol1_ee_eegenes, xol1sex1_xol1_ee_noteegenes)

plot_EE_xol1sex1_xol1 <- ggplot(xol1sex1_xol1_ee_eegenes_labeled, aes(x = as.factor(genotype) ,y = log2FoldChange)) + geom_boxplot(outlier.shape = NA) + 
  coord_cartesian(ylim = c(-2.3,3)) +
  geom_hline(linetype = 2, alpha = 0.3, yintercept = median(xol1sex1_xol1_ee_eegenes_labeled$log2FoldChange[xol1sex1_xol1_ee_eegenes_labeled$genotype == "not early embryo enriched genes"])) +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(size = 15), axis.title.y = element_text(size = 18)) +
  labs(title = "xol-1 sex-1/ xol-1")
plot_EE_xol1sex1_xol1

wilcox.test(xol1sex1_xol1_ee_eegenes_labeled$log2FoldChange[xol1sex1_xol1_ee_eegenes_labeled$genotype == "early embryo enriched genes"],
            xol1sex1_xol1_ee_eegenes_labeled$log2FoldChange[xol1sex1_xol1_ee_eegenes_labeled$genotype == "not early embryo enriched genes"])
#p-value = 1.818e-07

###moods median test
moodtest_eegenes_xol1sex1_xol1_eegenes <- xol1sex1_xol1_ee_eegenes_labeled[xol1sex1_xol1_ee_eegenes_labeled$genotype == "early embryo enriched genes", c("log2FoldChange", "genotype")]
moodtest_eegenes_xol1sex1_xol1_noteegenes <- xol1sex1_xol1_ee_eegenes_labeled[xol1sex1_xol1_ee_eegenes_labeled$genotype == "not early embryo enriched genes", c("log2FoldChange", "genotype")]
moodtest_eegenes_xol1sex1_xol1 <- rbind(moodtest_eegenes_xol1sex1_xol1_eegenes,
                                        moodtest_eegenes_xol1sex1_xol1_noteegenes)
mood.medtest(log2FoldChange ~ genotype,
             data = moodtest_eegenes_xol1sex1_xol1)
#p-value = 4.615e-05

ggsave(plot_EE_xol1sex1_xol1, filename = "boxplot_xs_x_EEgenes.png" , path = "salmon_analysis/xol1_sex1")

##
xol1_N2_ee_eegenes <- xol1_N2_ee_genes[xol1_N2_ee_genes$EnsemblID %in% EE_genes_core$gene_id,]
xol1_N2_ee_eegenes$genotype <- rep("early embryo enriched genes", nrow(xol1_N2_ee_eegenes))
xol1_N2_ee_noteegenes <- xol1_N2_ee_genes[!xol1_N2_ee_genes$EnsemblID %in% EE_genes_core$gene_id,]
xol1_N2_ee_noteegenes$genotype <- rep("not early embryo enriched genes", nrow(xol1_N2_ee_noteegenes))
xol1_N2_ee_eegenes_labeled <- rbind(xol1_N2_ee_eegenes, xol1_N2_ee_noteegenes)

plot_EE_xol1_N2 <- ggplot(xol1_N2_ee_eegenes_labeled, aes(x = as.factor(genotype),y = log2FoldChange)) + 
  geom_boxplot(outlier.shape = NA) + 
  coord_cartesian(ylim = c(-5,3)) +
  geom_hline(linetype = 2, alpha = 0.3, yintercept = median(xol1_N2_ee_eegenes_labeled$log2FoldChange[xol1_N2_ee_eegenes_labeled$genotype == "not early embryo enriched genes"])) +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(size = 15), axis.title.y = element_text(size = 18)) +
  labs(title = "xol-1 sex-1/ xol-1")
plot_EE_xol1_N2

wilcox.test(xol1_N2_ee_eegenes_labeled$log2FoldChange[xol1_N2_ee_eegenes_labeled$genotype == "early embryo enriched genes"],
            xol1_N2_ee_eegenes_labeled$log2FoldChange[xol1_N2_ee_eegenes_labeled$genotype == "not early embryo enriched genes"])
#p-value = 1.381e-06

###moods median test
moodtest_eegenes_xol1_wt_eegenes <- xol1_N2_ee_eegenes_labeled[xol1_N2_ee_eegenes_labeled$genotype == "early embryo enriched genes", c("log2FoldChange", "genotype")]
moodtest_eegenes_xol1_wt_noteegenes <- xol1_N2_ee_eegenes_labeled[xol1_N2_ee_eegenes_labeled$genotype == "not early embryo enriched genes", c("log2FoldChange", "genotype")]
moodtest_eegenes_xol1_wt <- rbind(moodtest_eegenes_xol1_wt_eegenes,
                                  moodtest_eegenes_xol1_wt_noteegenes)
mood.medtest(log2FoldChange ~ genotype,
             data = moodtest_eegenes_xol1_wt)
#p-value = 4.803e-06

ggsave(plot_EE_xol1_N2, filename = "boxplot_x_wt_EEgenes.png" , path = "salmon_analysis/xol1_sex1")

##combined

xol1_N2_ee_eegenes_labeled$dataset <- rep("xol1/WT", nrow(xol1_N2_ee_eegenes_labeled))
xol1sex1_N2_ee_eegenes_labeled$dataset <- rep("xol1_sex1/WT", nrow(xol1sex1_N2_ee_eegenes_labeled))
xol1sex1_xol1_ee_eegenes_labeled$dataset <- rep("xol1_sex1/xol1", nrow(xol1sex1_xol1_ee_eegenes_labeled))

combined_ee <- rbind(xol1_N2_ee_eegenes_labeled,
                     xol1sex1_N2_ee_eegenes_labeled,
                     xol1sex1_xol1_ee_eegenes_labeled)
combined_ee$dataset <- factor(combined_ee$dataset, levels = c("xol1/WT",
                                                              "xol1_sex1/WT",
                                                              "xol1_sex1/xol1"))

ggplot(combined_ee, aes(x = as.factor(genotype), y = log2FoldChange, fill = genotype)) + 
  stat_boxplot(geom = "errorbar", width = 0.25) +  
  geom_boxplot(outlier.shape = NA) + 
  facet_wrap(~ dataset) +
  coord_cartesian(ylim = c(-4.8,3.5)) +
  geom_hline(linetype = 2, alpha = 0.3, yintercept = 0) +
  scale_fill_manual(values = c("#b2d8d8","lightgrey")) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 11, color = "black"), 
        axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 18), 
        title = element_text(size = 8),
        legend.position = "none") +
  theme(axis.title.y = element_blank()) +
  scale_x_discrete(labels = c("Early-embryonic \n Genes", " All Other \n Genes"))

ggsave(filename = "boxplot_earlyemb_all.png",
       path = "sex-1 paper/Data/",
       height = 4,
       width = 8)
