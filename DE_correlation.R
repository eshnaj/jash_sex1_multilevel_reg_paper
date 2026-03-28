library('biomaRt')
library('ggplot2')
library('tidyr')
library('dplyr')
library('data.table')
library('DESeq2')
library('ggpubr')
library('readxl')
library("VennDiagram")

### raw data analysis ########
### volcano ##
xol1sex1_xol1_ee_raw <- read.table("salmon_analysis/xol1_sex1/xol1sex1_xol1_ee_raw.csv")
xol1_N2_ee_raw <- read.table("salmon_analysis/xol1_sex1/xol1_N2_ee_raw.csv")

volcano2 <- ggplot(data = xol1sex1_xol1_ee_raw, 
                   aes(x = log2FoldChange, y = -log10(padj), 
                       color = abs(log2FoldChange) > 1.5 & padj<0.05)) +
  geom_point(show.legend = FALSE, alpha = 0.3, size = 1.5) +
  scale_color_manual(values = c("lightgrey", "#006666")) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, alpha = 0.7) +
  geom_vline(xintercept = c(-1.5,1.5), linetype = 2, alpha = 0.7) +
  coord_cartesian(xlim = c(-15,15)) +
  labs(title = "xol-1 sex-1 vs xol-1", x = "Log2 Fold Change (Early Embryos)", y = "-log10 (adjusted p-value)") +
  theme(axis.title.x = element_text(size = 17), 
        axis.title.y = element_text(size = 17), 
        axis.text = element_text(size = 12))
volcano2

ggsave(plot = volcano2, 
       filename = "volcano_xol1sex1_xol1.tiff",
       path = "sex-1 paper/",
       height = 5,
       width = 6)

### venn ##

venn2_1 <- subset(xol1_N2_ee_raw, padj < 0.05)
venn2_2 <- subset(xol1sex1_N2_ee_raw, padj < 0.05)
venn2_3 <- subset(xol1sex1_xol1_ee_raw, padj < 0.05)

venn_4 <- list(xol1_WT = row.names(venn2_1),
               xol1sex1_xol1 = row.names(venn2_3))

venn.diagram(venn_4, 
             "sex-1 paper/Data/venn.png",
             scaled = 3, compression = "lzw",
             fill = c("#b2d8d8", "#006666"),
             category.names = "",
             fontface = 2)

### correlation ##

xol1_N2_ee_scatter <- xol1_N2_ee_raw[,c("log2FoldChange","padj")]
xol1_N2_ee_scatter$dataset <- rep("xol1_N2", nrow(xol1_N2_ee_scatter))
xol1sex1_xol1_ee_scatter <- xol1sex1_xol1_ee_raw[,c("log2FoldChange","padj")]
xol1sex1_xol1_ee_scatter$dataset <- rep("xol1sex1_xol1", nrow(xol1sex1_xol1_ee_scatter))


scatter1 <- cbind(xol1_N2_ee_scatter, xol1sex1_xol1_ee_scatter)
scatter1 <- na.omit(scatter1)
colnames(scatter1) <- c("log2FC_xol1_N2", "padj_xol1_N2", "dataset_1",
                        "log2FC_xol1sex1_xol1", "padj_xol1sex1_xol1", "dataset_2")

sig_correlation <- function(inputdataframe, padj1_col, padj2_col) {
  x = 1
  inputdataframe$sig <- rep("PH", nrow(inputdataframe))
  while(TRUE) {
    if(x > nrow(inputdataframe)) 
      break
    if(inputdataframe[x,padj1_col] < 0.05 & inputdataframe[x,padj2_col] < 0.05) {
      inputdataframe[x,"sig"] <- "yes"
    }
    else {
      inputdataframe[x,"sig"] <- "no" 
    }
    x = x+1
  }
  return(inputdataframe)
}

scatter1 <- sig_correlation(scatter1, "padj_xol1_N2", "padj_xol1sex1_xol1") 

# plot correlation
scatter1_plot <- ggplot(subset(scatter1, abs(log2FC_xol1_N2) > 0.5 & abs(log2FC_xol1sex1_xol1) > 0.5), 
                        aes(x = log2FC_xol1_N2, 
                            y = log2FC_xol1sex1_xol1, 
                            color = sig, 
                            alpha = sig)) +
  geom_point() +
  scale_color_manual(values = c("salmon", "darkblue")) +
  scale_alpha_manual(values = c(0.2, 0.9)) +
  theme_bw() +
  theme(legend.position = "none", axis.text = element_text(size = 10)) +
  stat_cor(method = "pearson")

scatter1_plot

# figure for paper

scatter1$sig <- as.factor(scatter1$sig)
scatter1_plot <- ggplot() +
  geom_point(data = subset(scatter1, 
                           abs(log2FC_xol1_N2) > 0.5 & 
                             abs(log2FC_xol1sex1_xol1) > 0.5 &
                             sig == "no"),
             aes(x = log2FC_xol1_N2, 
                 y = log2FC_xol1sex1_xol1),
             color = "salmon", alpha = 0.3,
             size = 2) +
  geom_point(data = subset(scatter1, 
                           abs(log2FC_xol1_N2) > 0.5 & 
                             abs(log2FC_xol1sex1_xol1) > 0.5 &
                             sig =="yes"), 
             aes(x = log2FC_xol1_N2, 
                 y = log2FC_xol1sex1_xol1),
             color = "darkblue", alpha = 0.6,
             size = 2) +
  theme_bw() +
  theme(legend.position = "none", axis.text = element_text(size = 14), axis.title = element_blank()) +
  annotate(geom = "text", 
           label = substitute(paste(italic('R '), '= - 0.91, ', italic('p '), '< 2.2e-16')), 
           x = 18, y = 22,
           color = "darkblue",
           size = 6)
scatter1_plot

ggsave(filename = "Fig3C.png",
       path = "sex-1 paper/Data/")


scatter_anticorr <- scatter3

scatter_anticorr$corr <- rep(NA, nrow(scatter_anticorr))
for (i in 1:nrow(scatter_anticorr)) {
  if (scatter_anticorr$log2FC_xol1_N2[i] > 0 & 
      scatter_anticorr$log2FC_xol1sex1_xol1[i] > 0) {
    scatter_anticorr$corr[i] <- "corr"
  }
  else if (scatter_anticorr$log2FC_xol1_N2[i] < 0 & 
           scatter_anticorr$log2FC_xol1sex1_xol1[i] < 0) {
    scatter_anticorr$corr[i] <- "corr"
  } else {
    scatter_anticorr$corr[i] <- "anti_corr"
  }
}

### linear modeling ##

scatter_anticorr <- subset(scatter_anticorr, sig == "yes")
scatter_anticorr <- subset(scatter_anticorr, corr == "anti_corr")

anticorr_plot_conf <- ggscatter(subset(scatter_anticorr, 
                                       abs(log2FC_xol1_N2) > 0.5 &
                                         abs(log2FC_xol1sex1_xol1) > 0.5),
                                x = "log2FC_xol1_N2", 
                                y = "log2FC_xol1sex1_xol1",
                                color = "darkblue", 
                                fill = "black",
                                alpha = 0.6,
                                size = 2,
                                add = "reg.line",
                                add.params = list(color = "black", 
                                                  fill = "darkgray", 
                                                  size = 0.3),
                                conf.int = TRUE,
                                cor.coef = TRUE,
                                conf.int.level = 0.99,
                                cor.coeff.args = list(method = "pearson", 
                                                      label.x = -21.5, 
                                                      label.y = -6,
                                                      label.sep = "\n")) +
  stat_regline_equation(label.y.npc = 'bottom', label.x.npc = 'left', size = 4) +
  geom_hline(yintercept = 0, linetype = 2, color = "lightgray") +
  geom_vline(xintercept = 0, linetype = 2, color = "lightgray") +
  theme_bw() +
  theme(axis.text = element_text(size = 14))

anticorr_plot_conf

#figure for paper
anticorr_plot_conf <- ggscatter(subset(scatter_anticorr, 
                                       abs(log2FC_xol1_N2) > 0.5 &
                                         abs(log2FC_xol1sex1_xol1) > 0.5),
                                x = "log2FC_xol1_N2", 
                                y = "log2FC_xol1sex1_xol1",
                                color = "darkblue", 
                                fill = "black",
                                alpha = 0.6,
                                size = 2,
                                add = "reg.line",
                                add.params = list(color = "black", 
                                                  fill = "darkgray", 
                                                  size = 0.3),
                                conf.int = TRUE,
                                cor.coef = FALSE,
                                conf.int.level = 0.99) +
  # stat_regline_equation(label.y.npc = 'bottom', label.x.npc = 'left', size = 4) +
  annotate(geom = "text", 
           label = "italic(R) == -0.96 * ',' ~ italic(p) < 2.2e-16",
           parse = TRUE,
           x = 5, y = 21,
           color = "darkblue",
           size = 5)  +
  annotate(geom = "text", 
           label = "italic(y) == 0.19 + 0.88 * italic(x)",
           parse = TRUE,
           x = 5, y = 23,
           color = "black",
           size = 6)  +
  geom_hline(yintercept = 0, linetype = 2, color = "lightgray") +
  geom_vline(xintercept = 0, linetype = 2, color = "lightgray") +
  theme_bw() +
  theme(axis.title = element_blank(), axis.text = element_text(size = 14))
anticorr_plot_conf

ggsave(anticorr_plot_conf, filename = "Fig3D.png", path = "sex-1 paper/Data/")

### correction for developmental timing ####

#importing gene expression datasets from https://www.vanderbilt.edu/wormdoc/wormmap/Welcome.html
EE_genes <- read.table("Genesets/EE_EGs.txt", header = TRUE)
LE_genes <- read.table("Genesets/LE_EGs.txt", header = TRUE)

common_embryo_genes <- EE_genes$gene_id[EE_genes$gene_id %in% LE_genes$gene_id] 

LE_genes_core <- LE_genes[!LE_genes$gene_id %in% common_embryo_genes,]
EE_genes_core <- EE_genes[!EE_genes$gene_id %in% common_embryo_genes,]

LE_EE_0.1FC <- EE_genes %>% 
  inner_join(LE_genes, by = "gene_id", suffix = c("_EE", "_LE")) %>%
  mutate(log2FC = log2(avg_expr_LE/avg_expr_EE)) %>%
  filter(abs(log2FC) > 0.1)

## Load datasets
ee <- readRDS("/Users/eshna/Desktop/Work/Csankovszki Lab/2026/early_emb/DE_run.rds")
xol1sex1_xol1 <- results(ee, contrast=c("condition","xol1sex1","xol1"))
xol1sex1_xol1 <- data.frame(xol1sex1_xol1)

xol1_N2 <- results(ee, contrast=c("condition","xol1","N2"))
xol1_N2 <- data.frame(xol1_N2)

## FILTER OUT GENES

xol1sex1_xol1$EnsemblID <- row.names(xol1sex1_xol1)
xol1sex1_xol1 <- xol1sex1_xol1[!is.na(xol1sex1_xol1$log2FoldChange),]
xol1sex1_xol1 <- subset(xol1sex1_xol1, !EnsemblID %in% LE_genes_core$gene_id)
xol1sex1_xol1 <- subset(xol1sex1_xol1, !EnsemblID %in% EE_genes_core$gene_id)

xol1_N2$EnsemblID <- row.names(xol1_N2)
xol1_N2 <- xol1_N2[!is.na(xol1_N2$log2FoldChange),]
xol1_N2 <- subset(xol1_N2, !EnsemblID %in% LE_genes_core$gene_id)
xol1_N2 <- subset(xol1_N2, !EnsemblID %in% EE_genes_core$gene_id)

## PLOT

xol1sex1_xol1_0.1 <- subset(xol1sex1_xol1, !EnsemblID %in% LE_EE_0.1FC$gene_id)
xol1_N2_0.1 <- subset(xol1_N2, !EnsemblID %in% LE_EE_0.1FC$gene_id)

xol1_N2_ee_scatter3 <- xol1_N2_0.1[,c("log2FoldChange","padj")]
xol1_N2_ee_scatter3$dataset <- rep("xol1_N2", nrow(xol1_N2_ee_scatter3))
xol1sex1_xol1_ee_scatter3 <- xol1sex1_xol1_0.1[,c("log2FoldChange","padj")]
xol1sex1_xol1_ee_scatter3$dataset <- rep("xol1sex1_xol1", nrow(xol1sex1_xol1_ee_scatter3))


scatter3 <- cbind(xol1_N2_ee_scatter3, xol1sex1_xol1_ee_scatter3)
scatter3 <- na.omit(scatter3)
colnames(scatter3) <- c("log2FC_xol1_N2", "padj_xol1_N2", "dataset_1",
                        "log2FC_xol1sex1_xol1", "padj_xol1sex1_xol1", "dataset_2")

scatter3 <- sig_correlation(scatter3, "padj_xol1_N2", "padj_xol1sex1_xol1") 

scatter3_plot <- ggplot(subset(scatter3, abs(log2FC_xol1_N2) > 0.5 & abs(log2FC_xol1sex1_xol1) > 0.5), 
                        aes(x = log2FC_xol1_N2, 
                            y = log2FC_xol1sex1_xol1, 
                            color = sig, 
                            alpha = sig)) +
  geom_point() +
  scale_color_manual(values = c("salmon", "darkblue")) +
  scale_alpha_manual(values = c(0.2, 0.9)) +
  theme_bw() +
  theme(legend.position = "none", axis.text = element_text(size = 10)) +
  stat_cor(method = "pearson")

scatter3_plot

#figure for paper

scatter3$sig <- as.factor(scatter3$sig)
scatter3_plot <- ggplot() +
  geom_point(data = subset(scatter3, 
                           abs(log2FC_xol1_N2) > 0.5 & 
                             abs(log2FC_xol1sex1_xol1) > 0.5 &
                             sig == "no"),
             aes(x = log2FC_xol1_N2, 
                 y = log2FC_xol1sex1_xol1),
             color = "salmon", alpha = 0.3,
             size = 2) +
  geom_point(data = subset(scatter3, 
                           abs(log2FC_xol1_N2) > 0.5 & 
                             abs(log2FC_xol1sex1_xol1) > 0.5 &
                             sig =="yes"), 
             aes(x = log2FC_xol1_N2, 
                 y = log2FC_xol1sex1_xol1),
             color = "darkblue", alpha = 0.6,
             size = 2) +
  theme_bw() +
  theme(legend.position = "none", axis.text = element_text(size = 14), axis.title = element_blank()) +
  annotate(geom = "text", 
           label = "italic(R) == -0.92 * ',' ~ italic(p) < 2.2e-16",
           parse = TRUE,
           x = 16, y = 22,
           color = "darkblue",
           size = 6)
scatter3_plot

ggsave(filename = "Fig3C.png",
       path = "sex-1 paper/Data/2026")

## linear modeling ###

scatter_anticorr <- scatter3

scatter_anticorr$corr <- rep(NA, nrow(scatter_anticorr))
for (i in 1:nrow(scatter_anticorr)) {
  if (scatter_anticorr$log2FC_xol1_N2[i] > 0 & 
      scatter_anticorr$log2FC_xol1sex1_xol1[i] > 0) {
    scatter_anticorr$corr[i] <- "corr"
  }
  else if (scatter_anticorr$log2FC_xol1_N2[i] < 0 & 
           scatter_anticorr$log2FC_xol1sex1_xol1[i] < 0) {
    scatter_anticorr$corr[i] <- "corr"
  } else {
    scatter_anticorr$corr[i] <- "anti_corr"
  }
}

scatter_anticorr <- subset(scatter_anticorr, sig == "yes")
scatter_anticorr <- subset(scatter_anticorr, corr == "anti_corr")

anticorr_plot <- ggplot(subset(scatter_anticorr, 
                               abs(log2FC_xol1_N2) > 0.5 & abs(log2FC_xol1sex1_xol1) > 0.5), 
                        aes(x = log2FC_xol1_N2, 
                            y = log2FC_xol1sex1_xol1)) +
  geom_point() +
  #  scale_color_manual(values = c("salmon", "darkblue")) +
  #  scale_alpha_manual(values = c(0.2, 0.9)) +
  theme_bw() +
  stat_cor(method = "pearson")
anticorr_plot

anticorr_plot_conf <- ggscatter(subset(scatter_anticorr, 
                                       abs(log2FC_xol1_N2) > 0.5 &
                                         abs(log2FC_xol1sex1_xol1) > 0.5),
                                x = "log2FC_xol1_N2", 
                                y = "log2FC_xol1sex1_xol1",
                                color = "darkblue", 
                                fill = "black",
                                alpha = 0.6,
                                size = 2,
                                add = "reg.line",
                                add.params = list(color = "black", 
                                                  fill = "darkgray", 
                                                  size = 0.3),
                                conf.int = TRUE,
                                cor.coef = TRUE,
                                conf.int.level = 0.99,
                                cor.coeff.args = list(method = "pearson", 
                                                      label.x = -21.5, 
                                                      label.y = -6,
                                                      label.sep = "\n")) +
  stat_regline_equation(label.y.npc = 'bottom', label.x.npc = 'left', size = 4) +
  geom_hline(yintercept = 0, linetype = 2, color = "lightgray") +
  geom_vline(xintercept = 0, linetype = 2, color = "lightgray") +
  theme_bw() +
  theme(axis.text = element_text(size = 14))
anticorr_plot_conf

#figure for paper
anticorr_plot_conf <- ggscatter(subset(scatter_anticorr, 
                                       abs(log2FC_xol1_N2) > 0.5 &
                                         abs(log2FC_xol1sex1_xol1) > 0.5),
                                x = "log2FC_xol1_N2", 
                                y = "log2FC_xol1sex1_xol1",
                                color = "darkblue", 
                                fill = "black",
                                alpha = 0.6,
                                size = 2,
                                add = "reg.line",
                                add.params = list(color = "black", 
                                                  fill = "darkgray", 
                                                  size = 0.3),
                                conf.int = TRUE,
                                cor.coef = FALSE,
                                conf.int.level = 0.99) +
  # stat_regline_equation(label.y.npc = 'bottom', label.x.npc = 'left', size = 4) +
  annotate(geom = "text", 
           label = "italic(R) == -0.96 * ',' ~ italic(p) < 2.2e-16",
           parse = TRUE,
           x = 5, y = 21,
           color = "darkblue",
           size = 5)  +
  annotate(geom = "text", 
           label = "italic(y) == 0.19 + 0.88 * italic(x)",
           parse = TRUE,
           x = 5, y = 23,
           color = "black",
           size = 6)  +
  geom_hline(yintercept = 0, linetype = 2, color = "lightgray") +
  geom_vline(xintercept = 0, linetype = 2, color = "lightgray") +
  theme_bw() +
  theme(axis.title = element_blank(), axis.text = element_text(size = 14))
anticorr_plot_conf

ggsave(anticorr_plot_conf, filename = "Fig3D.png", path = "sex-1 paper/Data/2026/")

## volcano ###

volcano2 <- ggplot(data = xol1sex1_xol1_0.1, 
                   aes(x = log2FoldChange, y = -log10(padj), 
                       color = abs(log2FoldChange) > 1.5 & padj<0.05)) +
  geom_point(show.legend = FALSE, alpha = 0.3, size = 1.5) +
  scale_color_manual(values = c("lightgrey", "#006666")) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, alpha = 0.7) +
  geom_vline(xintercept = c(-1.5,1.5), linetype = 2, alpha = 0.7) +
  coord_cartesian(xlim = c(-15,15)) +
  labs(title = "xol-1 sex-1 vs xol-1", x = "Log2 Fold Change (Early Embryos)", y = "-log10 (adjusted p-value)") +
  theme(axis.title.x = element_text(size = 17), 
        axis.title.y = element_text(size = 17), 
        axis.text = element_text(size = 12))
volcano2

ggsave(plot = volcano2, 
       filename = "Fig3A.png",
       path = "sex-1 paper/Data/2026/",
       height = 5,
       width = 6)

### venn diagram ####

venn2_1 <- subset(xol1_N2_0.1, padj < 0.05)
venn2_3 <- subset(xol1sex1_xol1_0.1, padj < 0.05)

venn_4 <- list(xol1_WT = row.names(venn2_1),
               xol1sex1_xol1 = row.names(venn2_3))

venn.diagram(venn_4, 
             "/Users/eshna/Desktop/Work/Csankovszki Lab/2026/Fig3B.png",
             scaled = 3, compression = "lzw",
             fill = c("#b2d8d8", "#006666"),
             category.names = "",
             fontface = 2)

