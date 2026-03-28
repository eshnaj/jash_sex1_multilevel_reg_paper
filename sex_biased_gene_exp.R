library('biomaRt')
library('ggplot2')
library('tidyr')
library('dplyr')
library('AnnotationHub')
library('ensembldb')
library('data.table')
library('DESeq2')

xol1sex1_N2_ee_genes <- read.table("salmon_analysis/xol1_sex1/xol1sex1_N2_ee_genes.csv")
xol1sex1_xol1_ee_genes <- read.table("salmon_analysis/xol1_sex1/xol1sex1_xol1_ee_genes.csv")
xol1_N2_ee_genes <- read.table("salmon_analysis/xol1_sex1/xol1_N2_ee_genes.csv")

time_resolved <- read.table("Boeck_2016/time_resolved_embryo_RNAseq_dcpm.txt", header = TRUE)

time_resolved_embryo <- time_resolved[,1:19]
time_resolved_embryo <- drop_na(time_resolved_embryo)
time_resolved_embryo <- time_resolved_embryo[rowSums(time_resolved_embryo[2:19]) > 0,]
time_resolved_direction_bin2 <- data.frame("WormbaseID" = time_resolved_embryo$WormbaseName,
                                           "direction" = rep(NA, nrow(time_resolved_embryo)))

for (i in 1:nrow(time_resolved_embryo)) {
  if ( rowSums(time_resolved_embryo[i,2:10]) > rowSums(time_resolved_embryo[i,11:19]) ) {
    time_resolved_direction_bin2[i,2] <- "up" 
  } 
  else if ( rowSums(time_resolved_embryo[i,2:10]) < rowSums(time_resolved_embryo[i,11:19]) ) {
    time_resolved_direction_bin2[i,2] <- "down"
  } else {
    time_resolved_direction_bin2[i,2] <- "unchanged"
  }
}

## male-biased genesets with correction
## since we expect expression to be moved higher in xol-1 sex-1/xol1
## and xol-1 sex-1 is younger than xol-1
## remove the genes that have higher expression in late embryos compared to early embryos

###loading him8 datasets

him8_data_ee <- read.table("DEraw_him8.csv")
him8_male_bias <- subset(him8_data_ee, baseMean > 1 & padj < 0.05 & log2FoldChange > 0)

#importing wormbase ID
mart<-useDataset("celegans_gene_ensembl", 
                 useMart("ENSEMBL_MART_ENSEMBL", host="https://www.ensembl.org"))
attributes <- as.data.frame(mart@attributes[["name"]])

genes<-rownames(him8_male_bias)
gene_list<- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name", "wormbase_cds"),
  filters = "ensembl_gene_id",
  values = genes,
  mart = mart, 
  useCache = FALSE)
head(gene_list)
gene_list<-data.frame(gene_list)

gene_list_cds <- gene_list[,"wormbase_cds"]
gene_list_cds <- gsub('[abcdefghijklmnopqrstuvwxyz]', '', gene_list_cds)

gene_list$wormbase_cds <- gene_list_cds
gene_list <- gene_list[!duplicated(gene_list[,1:2]),]

him8_male_bias <-merge.data.frame(him8_male_bias, gene_list, by.x="row.names" , by.y = "ensembl_gene_id")
head(him8_male_bias)

## xol-1 sex-1/ WT with log2FC > 2 cutoff
him8_male_bias_filter_bin2 <- him8_male_bias[!him8_male_bias$wormbase_cds %in% time_resolved_direction_bin2$WormbaseID[time_resolved_direction_bin2$direction == "down"],] 
him8_male_bias_filter_bin2 <- subset(him8_male_bias_filter_bin2, log2FoldChange > 2)

add_sex_bias_column <- function(inputdataframe, bias_dataframe, bias_name) {
  sex_bias <- inputdataframe[rownames(inputdataframe) %in% bias_dataframe$Row.names,]
  sex_bias$geneset <- rep(bias_name, nrow(sex_bias))
  
  sex_bias_rest <- inputdataframe[!rownames(inputdataframe) %in% bias_dataframe$Row.names,]
  sex_bias_rest$geneset <- rep("All Other Genes", nrow(sex_bias_rest))
  
  combined_ggplot <- rbind(sex_bias, sex_bias_rest)
  combined_ggplot$geneset <- as.factor(combined_ggplot$geneset)
  combined_ggplot$geneset <- factor(combined_ggplot$geneset, 
                                    levels = c(bias_name, "All Other Genes"))
  return(combined_ggplot)
}
plot_sex_bias <- function(inputdataframe){
  ggplot_sex_bias <- ggplot(inputdataframe, aes(x = geneset, 
                                                y = log2FoldChange, 
                                                fill = geneset)) + 
    geom_boxplot(outlier.shape = NA, show.legend = FALSE) +
    geom_hline(linetype = 2, 
               alpha = 0.3, 
               yintercept = median(inputdataframe$log2FoldChange[inputdataframe$geneset == "All Other Genes"])) +
    theme(axis.text.x = element_text(size = 18), axis.title.x = element_blank(), axis.title.y = element_text(size = 18), title = element_text(size = 18))
  return(ggplot_sex_bias)
}

xol1sex1_WT_ee <- read.table("salmon_analysis/xol1_sex1/xol1sex1_N2_ee_raw.csv")
xol1sex1_WT_ee <- drop_na(xol1sex1_WT_ee)

xol1sex1_WT_male_bias2 <- add_sex_bias_column(xol1sex1_WT_ee, him8_male_bias_filter_bin2, "Male-biased Genes")
plot_xol1sex1_WT_male_bias2 <- plot_sex_bias(xol1sex1_WT_male_bias2) + 
  coord_cartesian(ylim = c(-3,2))

plot_xol1sex1_WT_male_bias2 

ggsave(plot_xol1sex1_WT_male_bias2, 
       filename = "boxplot_malebias2_corrected_xs_wt_ee.png",
       path = "salmon_analysis/xol1_sex1/")

## xol-1 sex-1/ WT with log2FC > 1 cutoff

him8_male_bias_filter_bin1 <- him8_male_bias[!him8_male_bias$wormbase_cds %in% time_resolved_direction_bin2$WormbaseID[time_resolved_direction_bin2$direction == "down"],] 
him8_male_bias_filter_bin1 <- subset(him8_male_bias_filter_bin1, log2FoldChange > 1)

xol1sex1_WT_male_bias1 <- add_sex_bias_column(xol1sex1_WT_ee, him8_male_bias_filter_bin1, "Male-biased Genes")
plot_xol1sex1_WT_male_bias1 <- plot_sex_bias(xol1sex1_WT_male_bias1) +
  coord_cartesian(ylim = c(-3,2))

plot_xol1sex1_WT_male_bias1 

ggsave(plot_xol1sex1_WT_male_bias1, 
       filename = "boxplot_malebias1_corrected_xs_wt_ee.png",
       path = "salmon_analysis/xol1_sex1/")

wilcox.test(xol1sex1_WT_male_bias1$log2FoldChange[xol1sex1_WT_male_bias1$geneset == "Male-biased Genes"],
            xol1sex1_WT_male_bias1$log2FoldChange[xol1sex1_WT_male_bias1$geneset == "All Other Genes"])
#p-value < 2.2e-16

##for paper
ggplot(xol1sex1_WT_male_bias1, aes(x = geneset, 
                                   y = log2FoldChange, 
                                   fill = geneset)) +
  stat_boxplot(geom = "errorbar", width = 0.25) +  
  geom_boxplot(outlier.shape = NA, width = 0.5, show.legend = FALSE) +
  coord_cartesian(ylim = c(-3,2)) +
  geom_hline(yintercept = median(xol1sex1_WT_male_bias1$log2FoldChange[xol1sex1_WT_male_bias1$geneset == "All Other Genes"]), 
             linetype = 2, 
             alpha = 0.3) +
  scale_fill_manual(values = c("#5065A7","lightgrey")) +
  theme(axis.text.x = element_text(size = 16, color = "black"), 
        axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 18), 
        title = element_text(size = 8)) +
  theme(axis.title.y = element_blank()) +
  scale_x_discrete(labels = c("Male-biased \n Genes", " All Other \n Genes"))

ggsave(filename = "boxplot_malebias1_corrected_xs_wt_ee.png",
       path = "sex-1 paper/Data/")

## xol-1 sex-1/ xol-1 with log2FC > 2 cutoff

xol1sex1_xol1_ee <- read.table("salmon_analysis/xol1_sex1/xol1sex1_xol1_ee_raw.csv")
xol1sex1_xol1_ee <- drop_na(xol1sex1_xol1_ee)

xol1sex1_xol1_male_bias2 <- add_sex_bias_column(xol1sex1_xol1_ee, him8_male_bias_filter_bin2, "Male-biased Genes")
plot_xol1sex1_xol1_male_bias2 <- plot_sex_bias(xol1sex1_xol1_male_bias2) +
  coord_cartesian(ylim = c(-3,3.5))

plot_xol1sex1_xol1_male_bias2 

ggsave(plot_xol1sex1_xol1_male_bias2, 
       filename = "boxplot_malebias2_corrected_xs_x_ee.png",
       path = "salmon_analysis/xol1_sex1/")

## xol-1 sex-1/ xol-1 with log2FC > 1 cutoff

him8_male_bias_filter_bin1 <- him8_male_bias[!him8_male_bias$wormbase_cds %in% time_resolved_direction_bin2$WormbaseID[time_resolved_direction_bin2$direction == "down"],] 
him8_male_bias_filter_bin1 <- subset(him8_male_bias_filter_bin1, log2FoldChange > 1)

xol1sex1_xol1_male_bias1 <- add_sex_bias_column(xol1sex1_xol1_ee, him8_male_bias_filter_bin1, "Male-biased Genes")
plot_xol1sex1_xol1_male_bias1 <- plot_sex_bias(xol1sex1_xol1_male_bias1) +
  coord_cartesian(ylim = c(-2,2.3))

plot_xol1sex1_xol1_male_bias1 

ggsave(plot_xol1sex1_xol1_male_bias1, 
       filename = "boxplot_malebias1_corrected_xs_x_ee.png",
       path = "salmon_analysis/xol1_sex1/")

wilcox.test(xol1sex1_xol1_male_bias1$log2FoldChange[xol1sex1_xol1_male_bias1$geneset == "Male-biased Genes"],
            xol1sex1_xol1_male_bias1$log2FoldChange[xol1sex1_xol1_male_bias1$geneset == "All Other Genes"])
#p-value = 1.417e-11

##for paper
ggplot(xol1sex1_xol1_male_bias1, aes(x = geneset, 
                                     y = log2FoldChange, 
                                     fill = geneset)) +
  stat_boxplot(geom = "errorbar", width = 0.25) +  
  geom_boxplot(outlier.shape = NA, width = 0.5, show.legend = FALSE) +
  coord_cartesian(ylim = c(-1.9,2.3)) +
  geom_hline(yintercept = median(xol1sex1_xol1_male_bias1$log2FoldChange[xol1sex1_xol1_male_bias1$geneset == "All Other Genes"]), 
             linetype = 2, 
             alpha = 0.3) +
  scale_fill_manual(values = c("#5065A7","lightgrey")) +
  theme(axis.text.x = element_text(size = 16, color = "black"), 
        axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 18), 
        title = element_text(size = 8)) +
  theme(axis.title.y = element_blank()) +
  scale_x_discrete(labels = c("Male-biased \n Genes", " All Other \n Genes"))

ggsave(filename = "boxplot_malebias1_corrected_xs_x_ee.png",
       path = "sex-1 paper/Data/")

## xol-1/ WT with log2FC > 1 cutoff

him8_male_bias_filter_bin1 <- him8_male_bias[!him8_male_bias$wormbase_cds %in% time_resolved_direction_bin2$WormbaseID[time_resolved_direction_bin2$direction == "down"],] 
him8_male_bias_filter_bin1 <- subset(him8_male_bias_filter_bin1, log2FoldChange > 1)

xol1_WT_male_bias1 <- add_sex_bias_column(xol1_WT_ee, him8_male_bias_filter_bin1, "Male-biased Genes")
plot_xol1_WT_male_bias1 <- plot_sex_bias(xol1_WT_male_bias1) +
  coord_cartesian(ylim = c(-4.3,3))

plot_xol1_WT_male_bias1 

ggsave(plot_xol1_WT_male_bias1, 
       filename = "boxplot_malebias1_corrected_x_wt_ee.png",
       path = "salmon_analysis/xol1_sex1/")

wilcox.test(xol1_WT_male_bias1$log2FoldChange[xol1_WT_male_bias1$geneset == "Male-biased Genes"],
            xol1_WT_male_bias1$log2FoldChange[xol1_WT_male_bias1$geneset == "All Other Genes"])
#p-value < 2.2e-16

##for paper
ggplot(xol1_WT_male_bias1, aes(x = geneset, 
                               y = log2FoldChange, 
                               fill = geneset)) +
  stat_boxplot(geom = "errorbar", width = 0.25) +  
  geom_boxplot(outlier.shape = NA, width = 0.5, show.legend = FALSE) +
  coord_cartesian(ylim = c(-4.2,3)) +
  geom_hline(yintercept = median(xol1_WT_male_bias1$log2FoldChange[xol1_WT_male_bias1$geneset == "All Other Genes"]), 
             linetype = 2, 
             alpha = 0.3) +
  scale_fill_manual(values = c("#5065A7","lightgrey")) +
  theme(axis.text.x = element_text(size = 16, color = "black"), 
        axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 18), 
        title = element_text(size = 8)) +
  theme(axis.title.y = element_blank()) +
  scale_x_discrete(labels = c("Male-biased \n Genes", " All Other \n Genes"))

ggsave(filename = "boxplot_malebias1_corrected_x_wt_ee.png",
       path = "sex-1 paper/Data/")

####combined

xol1_WT_male_bias1$dataset <- rep("xol1/WT", nrow(xol1_WT_male_bias1))
xol1sex1_WT_male_bias1$dataset <- rep("xol1_sex1/WT", nrow(xol1sex1_WT_male_bias1))
xol1sex1_xol1_male_bias1$dataset <- rep("xol1_sex1/xol1", nrow(xol1sex1_xol1_male_bias1))

combined_male_bias1 <- rbind(xol1_WT_male_bias1,
                             xol1sex1_WT_male_bias1,
                             xol1sex1_xol1_male_bias1)
combined_male_bias1$dataset <- factor(combined_male_bias1$dataset, 
                                      levels = c("xol1/WT",
                                                 "xol1_sex1/WT",
                                                 "xol1_sex1/xol1"))

ggplot(combined_male_bias1, aes(x = geneset, 
                                y = log2FoldChange, 
                                fill = geneset)) +
  stat_boxplot(geom = "errorbar", width = 0.25) +  
  geom_boxplot(outlier.shape = NA, width = 0.5, show.legend = FALSE) +
  facet_wrap(~ dataset) +
  coord_cartesian(ylim = c(-4.2,3)) +
  geom_hline(yintercept = 0, 
             linetype = 2, 
             alpha = 0.3) +
  scale_fill_manual(values = c("#5065A7","lightgrey")) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10, color = "black"), 
        axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 18), 
        title = element_text(size = 8)) +
  theme(axis.title.y = element_blank()) +
  scale_x_discrete(labels = c("Male-biased \n Genes", " All Other \n Genes"))

ggsave(filename = "boxplos_malebias1_corrected_all.png",
       path = "sex-1 paper/Data/",
       height = 4,
       width = 8)
time_resolved_direction_bin2$direction <- as.factor(time_resolved_direction_bin2$direction)
summary(time_resolved_direction_bin2$direction)

######### herm-biased genes with boeck correction #######
## since we expect expression to be lower in xol-1 sex-1/xol1
## and xol-1 sex-1 is younger than xol-1
## remove the genes that have lower expression in late embryos compared to early embryos

him8_herm_bias <- subset(him8_data_ee, padj < 0.05 & log2FoldChange < 0)
#importing wormbase ID
mart<-useDataset("celegans_gene_ensembl", 
                 useMart("ENSEMBL_MART_ENSEMBL", host="https://www.ensembl.org"))
attributes <- as.data.frame(mart@attributes[["name"]])

genes<-rownames(him8_herm_bias)
gene_list<- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name", "wormbase_cds"),
  filters = "ensembl_gene_id",
  values = genes,
  mart = mart, 
  useCache = FALSE)
head(gene_list)
gene_list<-data.frame(gene_list)

gene_list_cds <- gene_list[,"wormbase_cds"]
gene_list_cds <- gsub('[abcdefghijklmnopqrstuvwxyz]', '', gene_list_cds)

gene_list$wormbase_cds <- gene_list_cds
gene_list <- gene_list[!duplicated(gene_list[,1:2]),]

him8_herm_bias <-merge.data.frame(him8_herm_bias, gene_list, by.x="row.names" , by.y = "ensembl_gene_id")
head(him8_herm_bias)

## xol-1 sex-1/WT with log2FC > 2 cutoff

him8_herm_bias_filter_bin2 <- him8_herm_bias[!him8_herm_bias$wormbase_cds %in% time_resolved_direction_bin2$WormbaseID[time_resolved_direction_bin2$direction == "up"],]
him8_herm_bias_filter_bin2 <- subset(him8_herm_bias_filter_bin2, log2FoldChange < -2)

xol1sex1_WT_herm_bias2 <- add_sex_bias_column(xol1sex1_WT_ee, him8_herm_bias_filter_bin2, "Hermaphrodite-biased Genes")
plot_xol1sex1_WT_herm_bias2 <- plot_sex_bias(xol1sex1_WT_herm_bias2) +
  coord_cartesian(ylim = c(-2,2))

plot_xol1sex1_WT_herm_bias2

ggsave(plot_xol1sex1_WT_herm_bias2, 
       filename = "boxplot_hermbias2_corrected_xs_wt_ee.png",
       path = "salmon_analysis/xol1_sex1/")

## xol-1 sex-1/WT with log2FC > 1 cutoff

him8_herm_bias_filter_bin1 <- him8_herm_bias[!him8_herm_bias$wormbase_cds %in% time_resolved_direction_bin2$WormbaseID[time_resolved_direction_bin2$direction == "up"],]
him8_herm_bias_filter_bin1 <- subset(him8_herm_bias_filter_bin1, log2FoldChange < -1)

xol1sex1_WT_herm_bias1 <- add_sex_bias_column(xol1sex1_WT_ee, him8_herm_bias_filter_bin1, "Hermaphrodite-biased Genes")
plot_xol1sex1_WT_herm_bias1 <- plot_sex_bias(xol1sex1_WT_herm_bias1) +
  coord_cartesian(ylim = c(-2,2))

plot_xol1sex1_WT_herm_bias1

ggsave(plot_xol1sex1_WT_herm_bias1, 
       filename = "boxplot_hermbias1_corrected_xs_wt_ee.png",
       path = "salmon_analysis/xol1_sex1/")

wilcox.test(xol1sex1_WT_herm_bias1$log2FoldChange[xol1sex1_WT_herm_bias1$geneset == "Hermaphrodite-biased Genes"],
            xol1sex1_WT_herm_bias1$log2FoldChange[xol1sex1_WT_herm_bias1$geneset == "All Other Genes"])
# p.value = 0.0001208

###moods median test
moodtest_herm1_xol1sex1_wt_hermgenes <- xol1sex1_WT_herm_bias1[xol1sex1_WT_herm_bias1$geneset == "Hermaphrodite-biased Genes", c("log2FoldChange", "geneset")]
moodtest_herm1_xol1sex1_wt_allother <- xol1sex1_WT_herm_bias1[xol1sex1_WT_herm_bias1$geneset == "All Other Genes", c("log2FoldChange", "geneset")]
moodtest_herm1_xol1sex1_wt <- rbind(moodtest_herm1_xol1sex1_wt_hermgenes,
                                    moodtest_herm1_xol1sex1_wt_allother)

mood.medtest(log2FoldChange ~ geneset,
             data = moodtest_herm1_xol1sex1_wt)

#for paper
ggplot(xol1sex1_WT_herm_bias1, aes(x = geneset, 
                                   y = log2FoldChange, 
                                   fill = geneset)) +
  stat_boxplot(geom = "errorbar", width = 0.25) +  
  geom_boxplot(outlier.shape = NA, width = 0.5, show.legend = FALSE) +
  coord_cartesian(ylim = c(-1.9,2.3)) +
  geom_hline(yintercept = median(xol1sex1_WT_herm_bias1$log2FoldChange[xol1sex1_WT_herm_bias1$geneset == "All Other Genes"]), 
             linetype = 2, 
             alpha = 0.3) +
  scale_fill_manual(values = c("#A367B1","lightgrey")) +
  theme(axis.text.x = element_text(size = 12, color = "black"), 
        axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 18), 
        title = element_text(size = 8)) +
  theme(axis.title.y = element_blank()) +
  scale_x_discrete(labels = c("Hermaphrodite-biased \n Genes", " All Other \n Genes"))

ggsave(filename = "boxplot_hermbias1_corrected_xs_wt_ee.png",
       path = "sex-1 paper/Data/",
       height = 4,
       width = 4)

## xol-1 sex-1/xol-1 with log2FC > 2 cutoff

him8_herm_bias_filter_bin2 <- him8_herm_bias[!him8_herm_bias$wormbase_cds %in% time_resolved_direction_bin2$WormbaseID[time_resolved_direction_bin2$direction == "up"],]
him8_herm_bias_filter_bin2 <- subset(him8_herm_bias_filter_bin2, log2FoldChange < -2)

xol1sex1_xol1_herm_bias2 <- add_sex_bias_column(xol1sex1_xol1_ee, him8_herm_bias_filter_bin2, "Hermaphrodite-biased Genes")
plot_xol1sex1_xol1_herm_bias2 <- plot_sex_bias(xol1sex1_xol1_herm_bias2) +
  coord_cartesian(ylim = c(-2,2))

plot_xol1sex1_xol1_herm_bias2

ggsave(plot_xol1sex1_xol1_herm_bias2, 
       filename = "boxplot_hermbias2_corrected_xs_x_ee.png",
       path = "salmon_analysis/xol1_sex1/")

## xol-1 sex-1/xol1 with log2FC > 1 cutoff

him8_herm_bias_filter_bin1 <- him8_herm_bias[!him8_herm_bias$wormbase_cds %in% time_resolved_direction_bin2$WormbaseID[time_resolved_direction_bin2$direction == "up"],]
him8_herm_bias_filter_bin1 <- subset(him8_herm_bias_filter_bin1, log2FoldChange < -1)

xol1sex1_xol1_herm_bias1 <- add_sex_bias_column(xol1sex1_xol1_ee, him8_herm_bias_filter_bin1, "Hermaphrodite-biased Genes")
plot_xol1sex1_xol1_herm_bias1 <- plot_sex_bias(xol1sex1_xol1_herm_bias1) +
  coord_cartesian(ylim = c(-2,2))

plot_xol1sex1_xol1_herm_bias1

ggsave(plot_xol1sex1_xol1_herm_bias1, 
       filename = "boxplot_hermbias1_corrected_xs_x_ee.png",
       path = "salmon_analysis/xol1_sex1/")

wilcox.test(xol1sex1_xol1_herm_bias1$log2FoldChange[xol1sex1_xol1_herm_bias1$geneset == "Hermaphrodite-biased Genes"],
            xol1sex1_xol1_herm_bias1$log2FoldChange[xol1sex1_xol1_herm_bias1$geneset == "All Other Genes"])
#p-value < 2.2e-16

###moods median test
moodtest_herm1_xol1sex1_xol1_hermgenes <- xol1sex1_xol1_herm_bias1[xol1sex1_xol1_herm_bias1$geneset == "Hermaphrodite-biased Genes", c("log2FoldChange", "geneset")]
moodtest_herm1_xol1sex1_xol1_allother <- xol1sex1_xol1_herm_bias1[xol1sex1_xol1_herm_bias1$geneset == "All Other Genes", c("log2FoldChange", "geneset")]
moodtest_herm1_xol1sex1_xol1 <- rbind(moodtest_herm1_xol1sex1_xol1_hermgenes,
                                      moodtest_herm1_xol1sex1_xol1_allother)
mood.medtest(log2FoldChange ~ geneset,
             data = moodtest_herm1_xol1sex1_xol1)
#p-value < 2.2e-16

ggplot(xol1sex1_xol1_herm_bias1, aes(x = geneset, 
                                     y = log2FoldChange, 
                                     fill = geneset)) +
  stat_boxplot(geom = "errorbar", width = 0.25) +  
  geom_boxplot(outlier.shape = NA, width = 0.5, show.legend = FALSE) +
  coord_cartesian(ylim = c(-2.1,1.7)) +
  geom_hline(yintercept = median(xol1sex1_xol1_herm_bias1$log2FoldChange[xol1sex1_xol1_herm_bias1$geneset == "All Other Genes"]), 
             linetype = 2, 
             alpha = 0.3) +
  scale_fill_manual(values = c("#A367B1","lightgrey")) +
  theme(axis.text.x = element_text(size = 12, color = "black"), 
        axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 18), 
        title = element_text(size = 8)) +
  theme(axis.title.y = element_blank()) +
  scale_x_discrete(labels = c("Hermaphrodite-biased \n Genes", " All Other \n Genes"))

ggsave(filename = "boxplot_hermbias1_corrected_xs_x_ee.png",
       path = "sex-1 paper/Data/",
       height = 4,
       width = 4)

## xol-1/ WT with log2FC > 1 cutoff

him8_herm_bias_filter_bin1 <- him8_herm_bias[!him8_herm_bias$wormbase_cds %in% time_resolved_direction_bin2$WormbaseID[time_resolved_direction_bin2$direction == "up"],]
him8_herm_bias_filter_bin1 <- subset(him8_herm_bias_filter_bin1, log2FoldChange < -1)

xol1_WT_ee <- read.table("salmon_analysis/xol1_sex1/xol1_N2_ee_raw.csv")
xol1_WT_ee <- drop_na(xol1_WT_ee)

xol1_WT_herm_bias1 <- add_sex_bias_column(xol1_WT_ee, him8_herm_bias_filter_bin1, "Hermaphrodite-biased Genes")
plot_xol1_WT_herm_bias1 <- plot_sex_bias(xol1_WT_herm_bias1) +
  coord_cartesian(ylim = c(-3,3))

plot_xol1_WT_herm_bias1

ggsave(plot_xol1_WT_herm_bias1, 
       filename = "boxplot_hermbias1_corrected_x_wt_ee.png",
       path = "salmon_analysis/xol1_sex1/")

wilcox.test(xol1_WT_herm_bias1$log2FoldChange[xol1_WT_herm_bias1$geneset == "Hermaphrodite-biased Genes"],
            xol1_WT_herm_bias1$log2FoldChange[xol1_WT_herm_bias1$geneset == "All Other Genes"])
# p.value < 2.2e-16

median(xol1_WT_herm_bias1$log2FoldChange[xol1_WT_herm_bias1$geneset == "Hermaphrodite-biased Genes"])
#[1] 0.5310543
median(xol1_WT_herm_bias1$log2FoldChange[xol1_WT_herm_bias1$geneset == "All Other Genes"])
#[1] -0.05881747

#for paper
ggplot(xol1_WT_herm_bias1, aes(x = geneset, 
                               y = log2FoldChange, 
                               fill = geneset)) +
  stat_boxplot(geom = "errorbar", width = 0.25) +  
  geom_boxplot(outlier.shape = NA, width = 0.5, show.legend = FALSE) +
  coord_cartesian(ylim = c(-2.4,2.4)) +
  geom_hline(yintercept = median(xol1sex1_WT_herm_bias1$log2FoldChange[xol1sex1_WT_herm_bias1$geneset == "All Other Genes"]), 
             linetype = 2, 
             alpha = 0.3) +
  scale_fill_manual(values = c("#A367B1","lightgrey")) +
  theme(axis.text.x = element_text(size = 12, color = "black"), 
        axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 18), 
        title = element_text(size = 8)) +
  theme(axis.title.y = element_blank()) +
  scale_x_discrete(labels = c("Hermaphrodite-biased \n Genes", " All Other \n Genes"))

ggsave(filename = "boxplot_hermbias1_corrected_x_wt_ee.png",
       path = "sex-1 paper/Data/",
       height = 4,
       width = 4)

####combined

xol1_WT_herm_bias1$dataset <- rep("xol1/WT", nrow(xol1_WT_herm_bias1))
xol1sex1_WT_herm_bias1$dataset <- rep("xol1_sex1/WT", nrow(xol1sex1_WT_herm_bias1))
xol1sex1_xol1_herm_bias1$dataset <- rep("xol1_sex1/xol1", nrow(xol1sex1_xol1_herm_bias1))

combined_herm_bias1 <- rbind(xol1_WT_herm_bias1,
                             xol1sex1_WT_herm_bias1,
                             xol1sex1_xol1_herm_bias1)
combined_herm_bias1$dataset <- factor(combined_herm_bias1$dataset, 
                                      levels = c("xol1/WT",
                                                 "xol1_sex1/WT",
                                                 "xol1_sex1/xol1"))

ggplot(combined_herm_bias1, aes(x = geneset, 
                                y = log2FoldChange, 
                                fill = geneset)) +
  stat_boxplot(geom = "errorbar", width = 0.25) +  
  geom_boxplot(outlier.shape = NA, width = 0.5, show.legend = FALSE) +
  facet_wrap(~ dataset) +
  coord_cartesian(ylim = c(-2.7,2.5)) +
  geom_hline(yintercept = 0, 
             linetype = 2, 
             alpha = 0.3) +
  scale_fill_manual(values = c("#A367B1","lightgrey")) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10, color = "black"), 
        axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 18), 
        title = element_text(size = 8)) +
  theme(axis.title.y = element_blank()) +
  scale_x_discrete(labels = c("Hermaphrodite-biased \n Genes", " All Other \n Genes"))

ggsave(filename = "boxplos_hermbias1_corrected_all.png",
       path = "sex-1 paper/Data/",
       height = 4,
       width = 8)
