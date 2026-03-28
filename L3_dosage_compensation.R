## analyzing xol-1 sex-1, xol-1 and WT samples
library('biomaRt')
library('ggplot2')
library('tidyr')
library('dplyr')
library('data.table')
library('tximport')
library('DESeq2')
library("ggpubr")
library('AnnotationHub')
library('ensembldb')


############# deseq2 ##########
data_xol1sex1_L3 <- data.frame(sampleID=rep("PH", 12),
                               condition=rep("PH", 12),
                               genotype_rep=rep("PH", 12))
#setting genotypes
data_xol1sex1_L3$genotype_rep <- c("N2_rep1",
                                   "N2_rep2",
                                   "N2_rep3",
                                   "N2_rep4",
                                   "xol1_rep1",
                                   "xol1_rep2",
                                   "xol1_rep3",
                                   "xol1_rep4",
                                   "xol1sex1_rep1",
                                   "xol1sex1_rep2",
                                   "xol1sex1_rep3",
                                   "xol1sex1_rep4")
#setting sampleIDs
for (i in c(1:12)) {
  y=paste0("11600-EJ-",i)
  data_xol1sex1_L3[i,1]= y
}

#setting condition
data_xol1sex1_L3$condition <- c(rep("N2", 4), rep("xol1", 4), rep("xol1sex1", 4))

#creating file path for DESeq2
dir_xol1sex1_L3 <- ""
list.files(dir_xol1sex1_L3)
files <- file.path(dir_xol1sex1_L3, 
                   paste0("salmon_counts_", data_xol1sex1_L3$sampleID), 
                   "quant.sf")
names(files) <- data_xol1sex1_L3$sampleID
files
all(file.exists(files))
#returned TRUE

#creating transcript ID to gene ID reference
celegansdb_formalclassobject <- query(AnnotationHub(), pattern = c("Caenorhabditis elegans", "EnsDb", "109"))
celegansdb <- celegansdb_formalclassobject[[1]]
genes <- genes(celegansdb)
tx2gene <- data.frame(TXNAME=genes$canonical_transcript,
                      GENEID=genes$gene_id)
head(tx2gene)

#using tximport to compile salmon output
salmon_xol1sex1_L3 <- tximport(files, type = "salmon", tx2gene = tx2gene)
names(salmon_xol1sex1_L3)
head(salmon_xol1sex1_L3$counts)

#piping into DESeq2
ddsxol1sex1L3 <- DESeqDataSetFromTximport(salmon_xol1sex1_L3,
                                          colData = data_xol1sex1_L3,
                                          design = ~ condition)
ddsxol1sex1L3$condition <- relevel(ddsxol1sex1L3$condition, ref = "N2")
DErunxol1sex1L3 <- DESeq(ddsxol1sex1L3)
DErun_xol1sex1L3_res <- results(DErunxol1sex1L3)
head(DErun_xol1sex1L3_res)

#saving DESeq2 results
saveRDS(DErunxol1sex1L3, 
        file = "results/xol1sex1_L3.rds")

#running PCA
rld <- rlog(DErunxol1sex1L3, blind = TRUE)
plotPCA(rld)
plotPCA(rld, intgroup="condition", returnData = TRUE)

#fpkm
DErunxol1sex1L3 <- readRDS("results/xol1sex1_L3.rds")
fpkm_l3 <- fpkm(DErunxol1sex1L3)
colnames(fpkm_l3) <- c("WT_rep1", "WT_rep2", "WT_rep3", "WT_rep4", 
                       "xol1_rep1", "xol1_rep2", "xol1_rep3", "xol1_rep4", 
                       "xol1sex1_rep1", "xol1sex1_rep2", "xol1sex1_rep3", "xol1sex1_rep4")
write.table(fpkm_l3, file = "fpkm_L3.txt", sep = "\t",quote = FALSE, na = "NA", row.names = TRUE, col.names = TRUE)

##### extracting comparisons ######

DErunxol1sex1L3 <- readRDS("results/xol1sex1_L3.rds")

xol1sex1_N2 <- results(DErunxol1sex1L3, contrast=c("condition","xol1sex1","N2"))
write.table(xol1sex1_N2,"results/xol1sex1_N2_raw.csv")

xol1sex1_xol1 <- results(DErunxol1sex1L3, contrast=c("condition","xol1sex1","xol1"))
write.table(xol1sex1_xol1,"/results/xol1sex1_xol1_raw.csv")

xol1_N2 <- results(DErunxol1sex1L3, contrast=c("condition","xol1","N2"))
write.table(xol1_N2,"results/xol1_N2_raw.csv")

xol1sex1_N2_anno <- annotate_deseq_res(xol1sex1_N2)
xol1sex1_xol1_anno <- annotate_deseq_res(xol1sex1_xol1)
xol1_N2_anno <- annotate_deseq_res(xol1_N2)

write.table(xol1sex1_N2_anno,"results/xol1sex1_N2_annotated.csv")
write.table(xol1sex1_xol1_anno,"results/xol1sex1_xol1_annotated.csv")
write.table(xol1_N2_anno,"results/xol1_N2_annotated.csv")

###### X chromosome derepression #####

xol1sex1_N2_anno <- annotate_XorA(xol1sex1_N2_anno)
xol1sex1_xol1_anno <- annotate_XorA(xol1sex1_xol1_anno)
xol1_N2_anno <- annotate_XorA(xol1_N2_anno)

xol1sex1_N2_anno <- drop_na(xol1sex1_N2_anno)
xol1sex1_xol1_anno <- drop_na(xol1sex1_xol1_anno)
xol1_N2_anno <- drop_na(xol1_N2_anno)

plot_xs_wt <- create_XA_boxplot(xol1sex1_N2_anno)
plot_xs_wt <- plot_xs_wt +
  coord_cartesian(ylim = c(-2,2))
plot_xs_wt
ggsave(plot=plot_xs_wt, filename="boxplot_xol1sex1_wt_x_derep.png", path="results/")

plot_xs_x <- create_XA_boxplot(xol1sex1_xol1_anno)
plot_xs_x <- plot_xs_x +
  coord_cartesian(ylim = c(-2,2))
plot_xs_x
ggsave(plot=plot_xs_x, filename="boxplot_xol1sex1_xol1_x_derep.png", path="results/")

##for paper
theme_set(theme_bw())
plot_xs_x <- ggplot(data = subset(xol1sex1_xol1_anno, baseMean>1), 
                    aes(x = XorA, y = log2FoldChange, fill = XorA)) +
  stat_boxplot(geom = "errorbar", width = 0.25) +  
  geom_boxplot(outlier.shape = NA, show.legend = FALSE, width = 0.5) +
  coord_cartesian(ylim = c(-1.7,2)) +
  geom_hline(yintercept = median(xol1sex1_xol1_anno$log2FoldChange[xol1sex1_xol1_anno$XorA == "Autosome"]), 
             linetype = 2, 
             alpha = 0.6) +
  scale_fill_manual(values = c("lightgrey","#008080")) +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(size = 12, color = "black"), 
        axis.title.y = element_blank(), 
        title = element_text(size = 16),
        axis.text.y = element_text(size = 10, color = "black")) +
  scale_x_discrete(labels = c("Autosomal Genes", " X-linked Genes"))
plot_xs_x

ggsave(plot = plot_xs_x,
       filename = "boxplot_XA_global_xs_x_l3.png",
       path = "sex-1 paper/Data/",
       height = 3,
       width = 4)

plot_x_wt <- create_XA_boxplot(xol1_N2_anno)
plot_x_wt <- plot_x_wt +
  coord_cartesian(ylim = c(-1.5,1.5))
plot_x_wt
ggsave(plot=plot_x_wt, filename="boxplot_xol1_wt_x_derep.png", path="results/")

##for paper
theme_set(theme_bw())
plot_x_wt <- ggplot(data = subset(xol1_N2_anno, baseMean>1), 
                    aes(x = XorA, y = log2FoldChange, fill = XorA)) +
  stat_boxplot(geom = "errorbar", width = 0.25) +  
  geom_boxplot(outlier.shape = NA, show.legend = FALSE, width = 0.5) +
  coord_cartesian(ylim = c(-1.2,1.4)) +
  geom_hline(yintercept = median(xol1_N2_anno$log2FoldChange[xol1_N2_anno$XorA == "Autosome"]), 
             linetype = 2, 
             alpha = 0.6) +
  scale_fill_manual(values = c("lightgrey","#008080")) +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(size = 12, color = "black"), 
        axis.title.y = element_blank(), 
        title = element_text(size = 16),
        axis.text.y = element_text(size = 10, color = "black")) +
  scale_x_discrete(labels = c("Autosomal Genes", " X-linked Genes"))
plot_x_wt

ggsave(plot = plot_x_wt,
       filename = "boxplot_XA_global_x_wt_l3.png",
       path = "sex-1 paper/Data/",
       height = 3,
       width = 4)

wilcox.test(xol1sex1_N2_anno$log2FoldChange[xol1sex1_N2_anno$XorA == "X chromosome"],
            xol1sex1_N2_anno$log2FoldChange[xol1sex1_N2_anno$XorA == "Autosome"])
#p-value = 1.655e-08
median(xol1sex1_N2_anno$log2FoldChange[xol1sex1_N2_anno$XorA == "X chromosome"])
#[1] -0.03836155
median(xol1sex1_N2_anno$log2FoldChange[xol1sex1_N2_anno$XorA == "Autosome"])
#[1] 0.04413783
median(xol1sex1_N2_anno$log2FoldChange[xol1sex1_N2_anno$XorA == "X chromosome"]) - median(xol1sex1_N2_anno$log2FoldChange[xol1sex1_N2_anno$XorA == "Autosome"])
#[1] -0.08249938


wilcox.test(xol1sex1_xol1_anno$log2FoldChange[xol1sex1_xol1_anno$XorA == "X chromosome"],
            xol1sex1_xol1_anno$log2FoldChange[xol1sex1_xol1_anno$XorA == "Autosome"])
#p-value = 0.08433
median(xol1sex1_xol1_anno$log2FoldChange[xol1sex1_xol1_anno$XorA == "X chromosome"])
#[1] 0.02258512
median(xol1sex1_xol1_anno$log2FoldChange[xol1sex1_xol1_anno$XorA == "Autosome"])
#[1] 0.04052173
median(xol1sex1_xol1_anno$log2FoldChange[xol1sex1_xol1_anno$XorA == "X chromosome"]) - median(xol1sex1_xol1_anno$log2FoldChange[xol1sex1_xol1_anno$XorA == "Autosome"])
#[1] -0.01793661


wilcox.test(xol1_N2_anno$log2FoldChange[xol1_N2_anno$XorA == "X chromosome"],
            xol1_N2_anno$log2FoldChange[xol1_N2_anno$XorA == "Autosome"])
#p-value = 2.134e-10
median(xol1_N2_anno$log2FoldChange[xol1_N2_anno$XorA == "X chromosome"])
#[1] -0.04948222
median(xol1_N2_anno$log2FoldChange[xol1_N2_anno$XorA == "Autosome"])
#[1] 0.01075838
median(xol1_N2_anno$log2FoldChange[xol1_N2_anno$XorA == "X chromosome"]) - median(xol1_N2_anno$log2FoldChange[xol1_N2_anno$XorA == "Autosome"])
#[1] -0.0602406


###### comparison to DPY-27::AID dataset #######

###### dpy-27 #######

##xol-1 sex-1/ xol-1 ####
dpy27AID <- read.table("Trombley_2024/DE_lists/tir1_dpy27AID_aux_vs_tir1_dpy27AID_noaux.csv")
dpy27_X <- subset(dpy27AID, 
                  log2FoldChange >= 1 & 
                    chromosome_name == "X")
dpy27_sensitive <- dpy27_X$EnsemblID

xol1sex1_xol1_L3 <- read.table("results/xol1sex1_xol1_annotated.csv")

X_derep_xs_x_dpy27 <- define_X_sensitive(xol1sex1_xol1_L3, dpy27_sensitive)

X_derep_xs_x_dpy27_plot <- make_X_sensitive_plot(X_derep_xs_x_dpy27, 
                                                 c(-1.5,2), 
                                                 "xol-1 sex-1/ xol-1 (L3 stage)",
                                                 "X-sensitive defined as genes with log2FC > 1 in dpy-27::AID depletion")
X_derep_xs_x_dpy27_plot

ggsave(X_derep_xs_x_dpy27_plot, filename = "boxplot_dpy27_sensitive_xs_x.png", path = "results/")

#for paper
X_derep_xs_x_dpy27$derep <- factor(X_derep_xs_x_dpy27$derep, 
                                   levels = c("autosome",
                                              "X_sensitive",
                                              "X_nonsensitive"))
ggplot(data=subset(X_derep_xs_x_dpy27, baseMean>1), 
       aes(x=derep, y=log2FoldChange, fill = derep)) + 
  stat_boxplot(geom = "errorbar", width = 0.25) +  
  geom_boxplot(outlier.shape = NA, width = 0.5, show.legend = FALSE) +
  coord_cartesian(ylim = c(-1.6, 2.7)) +
  geom_hline(yintercept = median(X_derep_xs_x_dpy27$log2FoldChange[X_derep_xs_x_dpy27$derep == "autosome"]), 
             linetype = 2, 
             alpha = 0.6) +
  scale_fill_manual(values = c("lightgrey","#008080", "#b2d8d8")) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10, color = "black"), 
        axis.title.x = element_blank(), 
        axis.text.y = element_text(size = 10), 
        title = element_text(size = 8),
        axis.title.y = element_blank()) +
  scale_x_discrete(labels = c("Autosomal Genes", 
                              "X Genes\nDPY-27 sensitive",
                              "X Genes\nOther"))

ggsave(filename = "boxplot_dpy27sensitive_xs_x_l3.png",
       path = "sex-1 paper/Data/",
       height = 3.5,
       width = 4.5)


run_wilcox_X_derep_sensitive(X_derep_xs_x_dpy27)
#[1] "X_sensitive vs X_nonsensitive = 1.09542821227579e-08"
#1] "X_sensitive vs autosome = 5.86899151085706e-06"
#[1] "X_nonsensitive vs autosome = 0.000655759134635008"

median(X_derep_xs_x_dpy27$log2FoldChange[X_derep_xs_x_dpy27$derep == "X_sensitive"]) - median(X_derep_xs_x_dpy27$log2FoldChange[X_derep_xs_x_dpy27$derep == "autosome"])
#[1] 0.1902943
median(X_derep_xs_x_dpy27$log2FoldChange[X_derep_xs_x_dpy27$derep == "X_sensitive"])
#[1] 0.2331744

###xol-1/ WT ####

xol1_N2_L3 <- read.table("results/xol1_N2_annotated.csv")

X_derep_x_wt_dpy27 <- define_X_sensitive(xol1_N2_L3, dpy27_sensitive)

X_derep_x_wt_dpy27_plot <- make_X_sensitive_plot(X_derep_x_wt_dpy27, 
                                                 c(-2,1.5), 
                                                 "xol-1 / WT (L3 stage)",
                                                 "X-sensitive defined as genes with log2FC > 1 in dpy-27::AID depletion")
X_derep_x_wt_dpy27_plot

ggsave(X_derep_x_wt_dpy27_plot, filename = "boxplot_dpy27_sensitive_x_wt.png", path = "results/")

run_wilcox_X_derep_sensitive(X_derep_x_wt_dpy27)
#[1] "X_sensitive vs X_nonsensitive = 6.42975723580618e-06"
#[1] "X_sensitive vs autosome = 1.15159332726006e-08"
#[1] "X_nonsensitive vs autosome = 4.10219816824093e-05"

#for paper
X_derep_x_wt_dpy27$derep <- factor(X_derep_x_wt_dpy27$derep, 
                                   levels = c("autosome",
                                              "X_sensitive",
                                              "X_nonsensitive"))
ggplot(data=subset(X_derep_x_wt_dpy27, baseMean>1), 
       aes(x=derep, y=log2FoldChange, fill = derep)) + 
  stat_boxplot(geom = "errorbar", width = 0.25) +  
  geom_boxplot(outlier.shape = NA, width = 0.5, show.legend = FALSE) +
  coord_cartesian(ylim = c(-1.8, 2.2)) +
  geom_hline(yintercept = median(X_derep_x_wt_dpy27$log2FoldChange[X_derep_x_wt_dpy27$derep == "autosome"]), 
             linetype = 2, 
             alpha = 0.6) +
  scale_fill_manual(values = c("lightgrey","#008080", "#b2d8d8")) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10, color = "black"), 
        axis.title.x = element_blank(), 
        axis.text.y = element_text(size = 10), 
        title = element_text(size = 8),
        axis.title.y = element_blank()) +
  scale_x_discrete(labels = c("Autosomal Genes", 
                              "X Genes\nDPY-27 sensitive",
                              "X Genes\nOther"))

ggsave(filename = "boxplot_dpy27sensitive_x_wt_l3.png",
       path = "sex-1 paper/Data/",
       height = 3.5,
       width = 4.5)

###### dpy-21 #######

## xol-1 sex-1/ xol-1 #####
dpy21 <- read.table("Trombley_2024/DE_lists/dpy21_WT.csv")
dpy21_X <- subset(dpy21, 
                  log2FoldChange >= 1 & 
                    chromosome_name == "X")
dpy21_sensitive <- dpy21_X$EnsemblID

#xol1sex1_xol1_L3 <- read.table("/Volumes/lsa-gyorgyi/Eshna/Eshna_xol1sex1_L3/Eshna_analysis/results/xol1sex1_xol1_annotated.csv")

X_derep_xs_x_dpy21 <- define_X_sensitive(xol1sex1_xol1_L3, dpy21_sensitive)

X_derep_xs_x_dpy21_plot <- make_X_sensitive_plot(X_derep_xs_x_dpy21, 
                                                 c(-2,2), 
                                                 "xol-1 sex-1/ xol-1",
                                                 "X-sensitive defined as genes with log2FC > 1 in dpy-21 mutant")
X_derep_xs_x_dpy21_plot

ggsave(X_derep_xs_x_dpy21_plot, filename = "boxplot_dpy21_sensitive_xs_x.png", path = "results/")

run_wilcox_X_derep_sensitive(X_derep_xs_x_dpy21)
#[1] "X_sensitive vs X_nonsensitive = 2.42514056362189e-16"
#[1] "X_sensitive vs autosome = 2.29247130091049e-08"
#[1] "X_nonsensitive vs autosome = 4.40334704875394e-08"

#for paper
X_derep_xs_x_dpy21$derep <- factor(X_derep_xs_x_dpy21$derep, 
                                   levels = c("autosome",
                                              "X_sensitive",
                                              "X_nonsensitive"))
ggplot(data=subset(X_derep_xs_x_dpy21, baseMean>1), 
       aes(x=derep, y=log2FoldChange, fill = derep)) + 
  stat_boxplot(geom = "errorbar", width = 0.25) +  
  geom_boxplot(outlier.shape = NA, width = 0.5, show.legend = FALSE) +
  coord_cartesian(ylim = c(-1.6, 2.7)) +
  geom_hline(yintercept = median(X_derep_xs_x_dpy21$log2FoldChange[X_derep_xs_x_dpy21$derep == "autosome"]), 
             linetype = 2, 
             alpha = 0.6) +
  scale_fill_manual(values = c("lightgrey","#008080", "#b2d8d8")) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10, color = "black"), 
        axis.title.x = element_blank(), 
        axis.text.y = element_text(size = 10), 
        title = element_text(size = 8),
        axis.title.y = element_blank()) +
  scale_x_discrete(labels = c("Autosomal Genes", 
                              "X Genes\nDPY-21 sensitive",
                              "X Genes\nOther"))

ggsave(filename = "boxplot_dpy21sensitive_xs_x_l3.png",
       path = "sex-1 paper/Data/",
       height = 3.5,
       width = 4.5)

median(X_derep_xs_x_dpy21$log2FoldChange[X_derep_xs_x_dpy21$derep == "X_sensitive"]) - median(X_derep_xs_x_dpy21$log2FoldChange[X_derep_xs_x_dpy21$derep == "autosome"])
#[1] 0.1251845
median(X_derep_xs_x_dpy21$log2FoldChange[X_derep_xs_x_dpy21$derep == "X_sensitive"])
#[1] 0.1680646


### xol-1 / WT ####

#xol1_N2_L3 <- read.table("/Volumes/lsa-gyorgyi/Eshna/Eshna_xol1sex1_L3/Eshna_analysis/results/xol1_N2_annotated.csv")

X_derep_x_wt_dpy21 <- define_X_sensitive(xol1_N2_L3, dpy21_sensitive)

X_derep_x_wt_dpy21_plot <- make_X_sensitive_plot(X_derep_x_wt_dpy21, 
                                                 c(-2,1.5), 
                                                 "xol-1 / WT (L3 stage)",
                                                 "X-sensitive defined as genes with log2FC > 1 in dpy-21 mutation")
X_derep_x_wt_dpy21_plot

ggsave(X_derep_x_wt_dpy21_plot, filename = "boxplot_dpy21_sensitive_x_wt.png", path = "results/")

run_wilcox_X_derep_sensitive(X_derep_x_wt_dpy21)
#[1] "X_sensitive vs X_nonsensitive = 3.71532579271909e-14"
#[1] "X_sensitive vs autosome = 9.00102616511374e-19"
#[1] "X_nonsensitive vs autosome = 0.151978775164334"

#for paper
X_derep_x_wt_dpy21$derep <- factor(X_derep_x_wt_dpy21$derep, 
                                   levels = c("autosome",
                                              "X_sensitive",
                                              "X_nonsensitive"))
ggplot(data=subset(X_derep_x_wt_dpy21, baseMean>1), 
       aes(x=derep, y=log2FoldChange, fill = derep)) + 
  stat_boxplot(geom = "errorbar", width = 0.25) +  
  geom_boxplot(outlier.shape = NA, width = 0.5, show.legend = FALSE) +
  coord_cartesian(ylim = c(-1.6, 2)) +
  geom_hline(yintercept = median(X_derep_x_wt_dpy21$log2FoldChange[X_derep_x_wt_dpy21$derep == "autosome"]), 
             linetype = 2, 
             alpha = 0.6) +
  scale_fill_manual(values = c("lightgrey","#008080", "#b2d8d8")) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10, color = "black"), 
        axis.title.x = element_blank(), 
        axis.text.y = element_text(size = 10), 
        title = element_text(size = 8),
        axis.title.y = element_blank()) +
  scale_x_discrete(labels = c("Autosomal Genes", 
                              "X Genes\nDPY-21 sensitive",
                              "X Genes\nOther"))

ggsave(filename = "boxplot_dpy21sensitive_x_wt_l3.png",
       path = "sex-1 paper/Data/",
       height = 3.5,
       width = 4.5)

##### functions #####

annotate_deseq_res <- function(input_dataframe){
  mart <- useDataset("celegans_gene_ensembl", 
                     useMart("ENSEMBL_MART_ENSEMBL", host="https://www.ensembl.org"))
  genes <- rownames(input_dataframe)
  gene_list <- getBM(
    attributes = c("ensembl_gene_id", "external_gene_name", "chromosome_name", "description"),
    filters = "ensembl_gene_id",
    values = genes,
    mart = mart, 
    useCache = FALSE)
  gene_list <- data.frame(gene_list)
  out_dataframe <- merge.data.frame(input_dataframe, gene_list, by.x= 0 , by.y = "ensembl_gene_id")
  colnames(out_dataframe)[1] <- "EnsemblID"
  return(out_dataframe)
}

annotate_XorA <- function(input_dataframe){
  XorA <- data.frame(XorA = character(), stringsAsFactors = FALSE)
  for (i in 1:nrow(input_dataframe)) {
    if (input_dataframe[i,"chromosome_name"] == "X"){
      XorA[i,1] <- "X chromosome"
    } else {
      XorA[i,1] <- "Autosome"
    }
  }
  XorA$EnsemblID <- input_dataframe$EnsemblID
  out_dataframe <- merge.data.frame(input_dataframe, XorA, by = "EnsemblID")
  return(out_dataframe)
}

create_XA_boxplot <- function(dataframe) {
  theme_set(theme_bw())
  boxplot <- ggplot(data = subset(dataframe, baseMean>1), 
                    aes(x = XorA, y = log2FoldChange, fill = XorA)) +
    geom_boxplot(outlier.shape = NA, show.legend = FALSE, width = 0.5) +
    geom_hline(yintercept = median(dataframe$log2FoldChange[dataframe$XorA == "Autosome"]), linetype = 2, alpha = 0.6) +
    theme(axis.title.x = element_blank(), axis.text.x = element_text(size = 18, color = "black"), axis.title.y = element_text(size = 18), title = element_text(size = 16), axis.text.y = element_text(size = 16, color = "black")) +
    labs(y = "Log2 Fold Change")
  return(boxplot)
}

define_X_unaffected <- function(inputdataframe, X_unaffected_list){
  
  df_X_unaffected <- subset(inputdataframe,
                            EnsemblID %in% X_unaffected_list &
                              chromosome_name == "X")
  df_X_unaffected$derep <- rep("X_unaffected", nrow(df_X_unaffected))
  df_X_affected <- subset(inputdataframe,
                          !EnsemblID %in% X_unaffected_list &
                            chromosome_name == "X")
  df_X_affected$derep <- rep("X_affected", nrow(df_X_affected))
  df_A <- subset(inputdataframe,
                 chromosome_name != "X")
  df_A$derep <- rep("autosome", nrow(df_A))
  
  out_df <- rbind(df_X_unaffected, df_X_affected, df_A)
  out_df <- subset(out_df, baseMean > 1)
  
  return(out_df)
}

make_X_sensitive_plot <- function(inputdataframe, coord, plot_title, plot_subtitle) {
  plot <- ggplot(inputdataframe, aes(x = derep, y = log2FoldChange)) + 
    geom_boxplot(outlier.shape = NA) +
    coord_cartesian(ylim = coord) +
    geom_hline(yintercept = median(inputdataframe$log2FoldChange[inputdataframe$derep == "autosome"]),
               linetype = 2,
               alpha = 0.5) +
    theme(axis.title.x = element_blank()) +
    labs(title = plot_title, subtitle = plot_subtitle)
  return(plot)
}

run_wilcox_X_derep_unaffected <- function(inputdataframe) {
  a <- wilcox.test(inputdataframe$log2FoldChange[inputdataframe$derep == "X_affected"],
                   inputdataframe$log2FoldChange[inputdataframe$derep == "X_unaffected"])[3]
  b <- wilcox.test(inputdataframe$log2FoldChange[inputdataframe$derep == "X_affected"],
                   inputdataframe$log2FoldChange[inputdataframe$derep == "autosome"])[3]
  c <- wilcox.test(inputdataframe$log2FoldChange[inputdataframe$derep == "X_unaffected"],
                   inputdataframe$log2FoldChange[inputdataframe$derep == "autosome"])[3]
  d <- list("X_affected vs X_unaffected" = a,
            "X_affected vs autosome" = b,
            "X_unaffected vs autosome" = c)
  print(paste0("X_affected vs X_unaffected = ", a))
  print(paste0("X_affected vs autosome = ", b))
  print(paste0("X_unaffected vs autosome = ", c))
}

run_wilcox_X_derep_sensitive <- function(inputdataframe) {
  a <- wilcox.test(inputdataframe$log2FoldChange[inputdataframe$derep == "X_sensitive"],
                   inputdataframe$log2FoldChange[inputdataframe$derep == "X_nonsensitive"])[3]
  b <- wilcox.test(inputdataframe$log2FoldChange[inputdataframe$derep == "X_sensitive"],
                   inputdataframe$log2FoldChange[inputdataframe$derep == "autosome"])[3]
  c <- wilcox.test(inputdataframe$log2FoldChange[inputdataframe$derep == "X_nonsensitive"],
                   inputdataframe$log2FoldChange[inputdataframe$derep == "autosome"])[3]
  d <- list("X_sensitive vs X_nonsensitive" = a,
            "X_sensitive vs autosome" = b,
            "X_nonsensitive vs autosome" = c)
  print(paste0("X_sensitive vs X_nonsensitive = ", a))
  print(paste0("X_sensitive vs autosome = ", b))
  print(paste0("X_nonsensitive vs autosome = ", c))
}

define_X_sensitive <- function(inputdataframe, X_sensitive_list){
  
  df_X_sensitive <- subset(inputdataframe,
                           EnsemblID %in% X_sensitive_list &
                             chromosome_name == "X")
  df_X_sensitive$derep <- rep("X_sensitive", nrow(df_X_sensitive))
  df_X_nonsensitive <- subset(inputdataframe,
                              !EnsemblID %in% X_sensitive_list &
                                chromosome_name == "X")
  df_X_nonsensitive$derep <- rep("X_nonsensitive", nrow(df_X_nonsensitive))
  df_A <- subset(inputdataframe,
                 chromosome_name != "X")
  df_A$derep <- rep("autosome", nrow(df_A))
  
  out_df <- rbind(df_X_sensitive, df_X_nonsensitive, df_A)
  out_df <- subset(out_df, baseMean > 1)
  
  return(out_df)
}

