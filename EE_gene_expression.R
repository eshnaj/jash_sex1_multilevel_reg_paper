library('biomaRt')
library('ggplot2')
library('tidyr')
library('dplyr')
library('AnnotationHub')
library('ensembldb')
library('data.table')
library('ggprism')
library('tximport')
library('DESeq2')
library("ggVennDiagram")
library("ggrepel")
library("VennDiagram")
library("ggpubr")
library("pheatmap")

data_xol1sex1_ee <- data.frame(sampleID=rep("PH", 9),
                               condition=rep("PH", 9),
                               genotype_rep=rep("PH", 9))
#setting genotypes
data_xol1sex1_ee$genotype_rep <- c("N2_rep1",
                                   "N2_rep2",
                                   "N2_rep3",
                                   "xol1_rep1",
                                   "xol1_rep2",
                                   "xol1_rep3",
                                   "xol1sex1_rep1",
                                   "xol1sex1_rep2",
                                   "xol1sex1_rep3")
#setting sampleIDs
for (i in c(1:9)) {
  y=paste0("3489-EJ-",i)
  data_xol1sex1_ee[i,1]= y
}

#setting condition
data_xol1sex1_ee$condition <- c(rep("N2", 3), rep("xol1", 3), rep("xol1sex1", 3))

#creating file path for DESeq2
dir_xol1sex1_ee <- "salmon_analysis"
list.files(dir_xol1sex1_ee)
files <- file.path(dir_xol1sex1_ee, 
                   paste0("counts_salmon_Sample_", data_xol1sex1_ee$sampleID), 
                   "quant.sf")
names(files) <- data_xol1sex1_ee$sampleID
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
salmon_xol1sex1_ee <- tximport(files, type = "salmon", tx2gene = tx2gene)
names(salmon_xol1sex1_ee)
head(salmon_xol1sex1_ee$counts)

#piping into DESeq2
ddsxol1sex1ee <- DESeqDataSetFromTximport(salmon_xol1sex1_ee,
                                          colData = data_xol1sex1_ee,
                                          design = ~ condition)
ddsxol1sex1ee$condition <- relevel(ddsxol1sex1ee$condition, ref = "N2")
DErunxol1sex1ee <- DESeq(ddsxol1sex1ee)
DErun_xol1sex1ee_res <- results(DErunxol1sex1ee)
head(DErun_xol1sex1ee_res)

#saving DESeq2 results
saveRDS(DErunxol1sex1ee, 
        file = "DE_run.rds")

#running PCA
rld <- rlog(DErunxol1sex1ee, blind = TRUE)
plotPCA(rld)
plotPCA(rld, intgroup="condition", returnData = TRUE)

#fpkm
DErun_ee_nxxs <- readRDS("DE_run.rds")
fpkm_ee <- fpkm(DErun_ee_nxxs)
colnames(fpkm_ee) <- c("WT_rep1", "WT_rep2", "WT_rep3", "xol1_rep1", "xol1_rep2", "xol1_rep3", "xol1sex1_rep1", "xol1sex1_rep2", "xol1sex1_rep3")
write.table(fpkm_ee, file = "fpkm_earlyemb.txt", sep = "\t",quote = FALSE, na = "NA", row.names = TRUE, col.names = TRUE)

############# extracting comparisons ########

DErunxol1sex1ee <- readRDS("DE_run.rds")

write.table(DErun_xol1sex1ee_res,"salmon_analysis/xol1_sex1/xol1sex1_N2_ee_raw.csv")

xol1sex1_xol1 <- results(DErunxol1sex1ee, contrast=c("condition","xol1sex1","xol1"))
xol1sex1_xol1 <- data.frame(xol1sex1_xol1)

write.table(xol1sex1_xol1, "salmon_analysis/xol1_sex1/xol1sex1_xol1_ee_raw.csv")

xol1_N2 <- results(DErunxol1sex1ee, contrast=c("condition","xol1","N2"))
xol1_N2 <- data.frame(xol1_N2)

write.table(xol1_N2, "salmon_analysis/xol1_sex1/xol1_N2_ee_raw.csv")


############# appending gene names ########

xol1sex1_N2_ee_data <- read.table("salmon_analysis/xol1_sex1/xol1sex1_N2_ee_raw.csv")

mart <- useDataset("celegans_gene_ensembl", useMart("ENSEMBL_MART_ENSEMBL", host="www.ensembl.org"))
genes <- row.names(xol1sex1_N2_ee_data)
gene_list<- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name", "chromosome_name", "description"),
  filters = "ensembl_gene_id",
  values = genes,
  mart = mart, 
  useCache = FALSE)
head(gene_list)

gene_list <- data.frame(gene_list)
xol1sex1_N2_ee_genes <- merge.data.frame(xol1sex1_N2_ee_data, gene_list, by.x=0 , by.y = "ensembl_gene_id")

XorA <- data.frame(XorA = character(), stringsAsFactors = FALSE)
for (i in 1:nrow(xol1sex1_N2_ee_genes)) {
  if (xol1sex1_N2_ee_genes[i,"chromosome_name"] == "X"){
    XorA[i,1] <- "X"
  } else {
    XorA[i,1] <- "A"
  }
}
row.names(XorA) <- xol1sex1_N2_ee_genes$Row.names
head(XorA)

xol1sex1_N2_ee_genes <- merge.data.frame(xol1sex1_N2_ee_genes, XorA, by.x = "Row.names", by.y = 0)
str(xol1sex1_N2_ee_genes)
xol1sex1_N2_ee_genes$XorA <- as.factor(xol1sex1_N2_ee_genes$XorA)
xol1sex1_N2_ee_genes <- drop_na(xol1sex1_N2_ee_genes)

###

xol1_N2_ee_data <- read.table("salmon_analysis/xol1_sex1/xol1_N2_ee_raw.csv")

mart <- useDataset("celegans_gene_ensembl", useMart("ENSEMBL_MART_ENSEMBL", host="www.ensembl.org"))
genes <- row.names(xol1_N2_ee_data)
gene_list<- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name", "chromosome_name", "description"),
  filters = "ensembl_gene_id",
  values = genes,
  mart = mart, 
  useCache = FALSE)
head(gene_list)

gene_list <- data.frame(gene_list)
xol1_N2_ee_genes <- merge.data.frame(xol1_N2_ee_data, gene_list, by.x=0 , by.y = "ensembl_gene_id")

XorA <- data.frame(XorA = character(), stringsAsFactors = FALSE)
for (i in 1:nrow(xol1_N2_ee_genes)) {
  if (xol1_N2_ee_genes[i,"chromosome_name"] == "X"){
    XorA[i,1] <- "X"
  } else {
    XorA[i,1] <- "A"
  }
}
row.names(XorA) <- xol1_N2_ee_genes$Row.names
head(XorA)

xol1_N2_ee_genes <- merge.data.frame(xol1_N2_ee_genes, XorA, by.x = "Row.names", by.y = 0)
str(xol1_N2_ee_genes)
xol1_N2_ee_genes$XorA <- as.factor(xol1_N2_ee_genes$XorA)
xol1_N2_ee_genes <- drop_na(xol1_N2_ee_genes)

write.table(xol1_N2_ee_genes, "salmon_analysis/xol1_sex1/xol1_N2_ee_genes.csv")

