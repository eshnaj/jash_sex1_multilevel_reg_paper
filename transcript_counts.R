library('biomaRt')
library('ggtext')
library('stringr')
library('ggprism')
library('ggplot2')

time_resolved <- read.table("time_resolved_embryo_RNAseq.txt", 
                            header = TRUE)

dim(time_resolved)
time_resolved_embryo <- time_resolved[,-20:-32]

mart<-useDataset("celegans_gene_ensembl", 
                   useMart("ENSEMBL_MART_ENSEMBL", 
                   host="https://www.ensembl.org"))
genes<-time_resolved_embryo$WormbaseName
gene_list<- getBM(
attributes = c("ensembl_gene_id", "external_gene_name", "wormbase_cds"),
filters = "wormbase_cds",
values = genes,
mart = mart, 
useCache = FALSE)

time_resolved_embryo <- merge.data.frame(time_resolved_embryo, 
                                         gene_list, 
                                         by.x = "WormbaseName", 
                                         by.y = "wormbase_cds", 
                                         all.x = TRUE)
  
time_resolved_embryo <- time_resolved_embryo[rowSums(time_resolved_embryo[,c(-1, -20:-21)]) > 0,]

#sex-1
sex1 <- time_resolved_embryo[time_resolved_embryo$WormbaseName == "F44A6.2",]
sex1 <- sex1[1,]
sex1 <- as.data.frame(t(sex1))
sex1$embryo_time <- rownames(sex1)
colnames(sex1) <- c("sex1", "embryo_time")
sex1 <- sex1[-1,]
sex1 <- sex1[-19:-20,]

sex1$embryo_time <- gsub("X", "", sex1$embryo_time)
sex1$embryo_time <- as.factor(sex1$embryo_time)
sex1$sex1 <- as.numeric(sex1$sex1)
sex1$sex1 <- round(sex1$sex1, digits = 2)
sex1$embryo_time <- gsub("min", " min", sex1$embryo_time)
sex1$embryo_time <- gsub("cell", " cell", sex1$embryo_time)
sex1$embryo_time = factor(sex1$embryo_time, levels = c("4 cell", 
                                                       "44 min", 
                                                       "83 min", 
                                                       "122 min", 
                                                       "161 min", 
                                                       "199 min", 
                                                       "238 min", 
                                                       "277 min", 
                                                       "316 min", 
                                                       "355 min", 
                                                       "393 min", 
                                                       "432 min", 
                                                       "471 min", 
                                                       "510 min", 
                                                       "548 min", 
                                                       "587 min", 
                                                       "626 min", 
                                                       "665 min"))

sex1_counts <- ggplot(sex1, aes(y = sex1, x = embryo_time, group = 1)) + 
  geom_line(color = "#5065A7", alpha = 0.5) + 
  geom_point(color = "#5065A7", size = 3, alpha = 0.9) +
  theme_prism(base_line_size = 0.5) +
  labs(y =  substitute(paste(italic('xol-1  '), 'transcripts','\n')), 
       x = " Time In Embryo Development") +
  theme(axis.text.x = element_text(angle = 90, 
                                   vjust = 0.5, 
                                   hjust=1, 
                                   size = 15),
        axis.text.y = element_text(size = 12), 
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18, family = "Arial Bold")) 

sex1_counts

ggsave(plot = sex1_counts, 
       filename = "sex1_transcripts.tiff", 
       path = "sex-1 paper/Data/", 
       height = 5, 
       width = 7)
