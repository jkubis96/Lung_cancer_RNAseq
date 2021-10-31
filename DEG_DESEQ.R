
##Download files
# rnaSEQ_manifest <- read.csv("data/HTseq/Metadata/gdc_manifest_20211023_160225.txt", sep = "\t", header = T)
# rnaSEQ_manifest <- rnaSEQ_manifest[grep('counts', rnaSEQ_manifest$filename),]
# rnaSEQ_manifest <- rnaSEQ_manifest[1:200,]
# 
# write.table(rnaSEQ_manifest, 'data/manifest.txt', sep =  '\t', row.names = F, quote = F)
# 


##Read data

#BiocManager::install("DESeq2")

list_files <- list.files(path="data/HTseq/Mutual/", pattern="*counts", full.names=TRUE, recursive=FALSE)

n = 0
for (files in list_files) {
  if (n == 0) {
    df <- read.csv(files, sep = '\t', col.names = c('gen', gsub('*.(.+?/.+?/.+?)/', '', files)), check.names = F)
    n = n + 1
    print(n)
  } else if (n > 0) {
    n = n + 1
    print(n)
    df <- merge(df,read.csv(files, sep = '\t', col.names = c('gen', gsub('*.(.+?/.+?/.+?)/', '', files)), check.names = F), by  = 'gen')
  }
}  

 df <- df[-grep('__', df$gen),]

##
# Change name
HGNC <- read.csv("data/genes/gencode.gene.info.v22.tsv", sep = "\t")[,1:2]
genes <- as.data.frame(df$gen)
colnames(genes)[1] <- 'gene_id'

library(dplyr)

genes <- left_join(genes, HGNC)
df$gen <- genes$gene_name

df <- aggregate(x = df[2:length(df)], list(df$gen), FUN = sum)
rownames(df) <- df$Group.1
df <- df[,-1]
##
## Metadata

sample_sheet <- read.csv("data/HTseq/Metadata/gdc_sample_sheet.2021-10-25.tsv", sep = "\t", header = T)
sample_sheet <- sample_sheet[sample_sheet$File.Name %in% colnames(df),]
colnames(sample_sheet)[colnames(sample_sheet) %in% 'Case.ID'] <- 'id'
sample_sheet <- sample_sheet[sample_sheet$Sample.Type %in% 'Primary Tumor',]
sample_sheet <- sample_sheet[!duplicated(sample_sheet$Sample.ID),]

df <- df[, colnames(df) %in% sample_sheet$File.Name]

for (col in 1:length(colnames(df))) {
  for (row in 1:length(rownames(sample_sheet))) {
    if (colnames(df)[col] == sample_sheet$File.Name[row]) {
      colnames(df)[col] <- sample_sheet$Sample.ID[row]
    }
      
  }
  
}


# DSEQ - metadata prepare

clinical <- read.csv("data/HTseq/Metadata/clinical.tsv", sep = "\t", header = T)
colnames(clinical)[colnames(clinical) %in% 'case_submitter_id'] <- 'id'


meta <- as.data.frame(colnames(df))
colnames(meta)[1] <- 'Sample.ID'
meta <- left_join(meta, sample_sheet[,c('Sample.Type','id', 'Sample.ID')], by  = 'Sample.ID')
meta <- left_join(meta, clinical[,c('ajcc_pathologic_stage','vital_status','age_at_index','race','gender','id', 'treatment_or_therapy')], by  = 'id')
meta <- meta[meta$Sample.ID %in% colnames(df),]
meta <- meta[!duplicated(meta$Sample.ID),]
rownames(meta) <- meta$Sample.ID


#male
meta_male <- meta[meta$gender %in% 'male',]
meta_male <- meta_male[order(rownames(meta_male)), , drop = F]
df_male <- df[,colnames(df) %in% meta_male$Sample.ID]
df_male <- df_male[, order(colnames(df_male)) , drop = F]


# Data normalization
#male

library(DESeq2)
library(tidyverse)

dds <- DESeqDataSetFromMatrix(countData = df_male, colData = meta_male, design = ~  treatment_or_therapy + age_at_index)
dds$treatment_or_therapy <- relevel(dds$treatment_or_therapy, ref = "no")
dds <- estimateSizeFactors(dds)
normalized_matrix_male <- counts(dds, normalized=TRUE)



jpeg(file.path("results/PCA_plot_man.jpeg") , units="in", width = 15, height = 10, res=600)
pcs <- vst(dds)
plotPCA(pcs, intgroup = c('treatment_or_therapy', 'vital_status'))
dev.off()




dds <- DESeq(dds, test = 'Wald')



write.csv(normalized_matrix_male, 'results/normalized_data_male.csv', row.names =  T)



results_male <- as.data.frame(results(dds))
results_male <- drop_na(results_male)
results_male$log_pval <- -log(results_male$padj)

#Vulcano

results_male$labels <- ''
results_male$labels[head(order(results_male$log_pval, decreasing = T), n = 100)] <- rownames(results_male)[head(order(results_male$log_pval, decreasing = T), n = 100)]
results_male$regulation <- NA
results_male$regulation[results_male$log2FoldChange < 0] <- "Downregulated"
results_male$regulation[results_male$log2FoldChange > 0] <- "Upregulated"
results_male$regulation[results_male$padj > 0.05] <- "No_Significant"

jpeg(file.path("results/Vulcano_plot_man.jpeg") , units="in", width = 15, height = 10, res=600)
ggplot(results_male, aes(x = log2FoldChange, y = log_pval), na.rm = T) + 
  geom_point(aes(color = regulation)) +
  ggtitle('Difference genes') +
  labs(color='Regulation') 

dev.off()
#############

results_male <- results_male[results_male$padj <= 0.05,]
results_male$labels <- rownames(results_male)
write.csv(results_male, 'results/DESEQ_male.csv')

jpeg(file.path("results/MA_plot_man.jpeg") , units="in", width = 15, height = 10, res=600)
plotMA(dds)
dev.off()

#heatmap
markers_up <- rownames(results_male)[order(results_male$log2FoldChange, decreasing = T) ,drop = F][1:10]
markers_down <- rownames(results_male)[order(results_male$log2FoldChange, decreasing = F) ,drop = F][1:10]

plotCounts(dds, gene="LYPD4", intgroup="treatment_or_therapy")

sample <- DataFrame(meta_male$treatment_or_therapy, meta_male$vital_status)
colnames(sample)[1:2] <- c('Treatment','Status')
gen <- as.data.frame(c(markers_up, markers_down))
gen$reg <- 'DOWN'
gen$reg[gen$`c(markers_up, markers_down)` %in% markers_up] <- 'UP'
rownames(gen) <- gen$`c(markers_up, markers_down)`
gen <- as.data.frame(gen[,-1, drop = F])
colnames(gen)[1] <- 'Regulation'

rownames(sample) <- meta_male$Sample.ID
sample <- as.data.frame(sample)
normalized_matrix_heat <- (normalized_matrix_male[rownames(normalized_matrix_male) %in% rownames(gen),])
jpeg(file.path("results/pheatmap.jpeg"),units="in", width = 30, height = 10 ,  res=300)
pheatmap::pheatmap(log(normalized_matrix_heat + 1), 
                   clustering_method = 'ward.D',
                   fontsize = 8,
                   annotation_col = sample,
                   annotation_row = gen,
                   cluster_rows = F,
                   angle_col = 270, fontsize_row = 6, fontsize_col = 6)
dev.off()



check <- as.data.frame(t(normalized_matrix_heat))
check$treat <- meta_male$treatment_or_therapy



check <- pivot_longer(check, cols = unique(rownames(gen)))
check$name <- factor(check$name, levels = unique(rownames(gen)))
check <- unique(check)
check$value <- log(check$value, 2)


library(ggpubr)

g <- ggboxplot(check, x = 'treat', y = 'value', color = 'name',
               add = "jitter")+
  stat_compare_means(method = 'kruskal.test', hide.ns = T) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 12, face = "bold", color = "black"),
        legend.position = "none",
        axis.text.y = element_text(angle = 90, hjust = 0.5, size = 12, face = "bold", color = "black"),
        axis.text=element_text(size=14),
        axis.title=element_text(size=14,face="bold")) +
  xlab('Related genes') +
  ylab('Expression')+
  theme(strip.text.x = element_text(
    size = 14, color = "red", face = "bold.italic"
  )) +
  facet_grid(~name) 


ggsave(g, filename = paste("results/man_genes.png", sep = ''), dpi = 300, units = "in", height = 10, width = 50, limitsize = F)







#female
meta_female <- meta[meta$gender %in% 'female',]
meta_female <- meta_female[order(rownames(meta_female)), , drop = F]
df_female <- df[,colnames(df) %in% meta_female$Sample.ID]
df_female <- df_female[, order(colnames(df_female)) , drop = F]


# Data normalization
#df_female



dds <- DESeqDataSetFromMatrix(countData = df_female, colData = meta_female, design = ~  treatment_or_therapy + age_at_index)
dds$treatment_or_therapy <- relevel(dds$treatment_or_therapy, ref = "no")
dds <- estimateSizeFactors(dds)
normalized_matrix_female <- counts(dds, normalized=TRUE)



jpeg(file.path("results/PCA_plot_woman.jpeg") , units="in", width = 15, height = 10, res=600)
pcs <- vst(dds)
plotPCA(pcs, intgroup = c('treatment_or_therapy', 'vital_status'))
dev.off()




dds <- DESeq(dds, test = 'Wald')



write.csv(normalized_matrix_female, 'results/normalized_data_female.csv', row.names =  T)



results_female <- as.data.frame(results(dds))
results_female <- drop_na(results_female)
results_female$log_pval <- -log(results_female$padj)

#Vulcano

results_female$labels <- ''
results_female$labels[head(order(results_female$log_pval, decreasing = T), n = 100)] <- rownames(results_female)[head(order(results_female$log_pval, decreasing = T), n = 100)]
results_female$regulation <- NA
results_female$regulation[results_female$log2FoldChange < 0] <- "Downregulated"
results_female$regulation[results_female$log2FoldChange > 0] <- "Upregulated"
results_female$regulation[results_female$padj > 0.05] <- "No_Significant"

jpeg(file.path("results/Vulcano_plot_woman.jpeg") , units="in", width = 15, height = 10, res=600)
ggplot(results_female, aes(x = log2FoldChange, y = log_pval), na.rm = T) + 
  geom_point(aes(color = regulation)) +
  ggtitle('Difference genes') +
  labs(color='Regulation') 

dev.off()
#############

results_female <- results_female[results_female$padj <= 0.05,]
results_female$labels <- rownames(results_female)
write.csv(results_female, 'results/DESEQ_female.csv')

jpeg(file.path("results/MA_plot_woman.jpeg") , units="in", width = 15, height = 10, res=600)
plotMA(dds)
dev.off()

#heatmap
markers_up <- rownames(results_female)[order(results_female$log2FoldChange, decreasing = T) ,drop = F][1:10]
markers_down <- rownames(results_female)[order(results_female$log2FoldChange, decreasing = F) ,drop = F][1:10]

#plotCounts(dds, gene="LYPD4", intgroup="treatment_or_therapy")

sample <- DataFrame(meta_female$treatment_or_therapy, meta_female$vital_status)
colnames(sample)[1:2] <- c('Treatment','Status')
gen <- as.data.frame(c(markers_up, markers_down))
gen$reg <- 'DOWN'
gen$reg[gen$`c(markers_up, markers_down)` %in% markers_up] <- 'UP'
rownames(gen) <- gen$`c(markers_up, markers_down)`
gen <- as.data.frame(gen[,-1, drop = F])
colnames(gen)[1] <- 'Regulation'

rownames(sample) <- meta_female$Sample.ID
sample <- as.data.frame(sample)
normalized_matrix_heat <- (normalized_matrix_female[rownames(normalized_matrix_female) %in% rownames(gen),])
jpeg(file.path("results/pheatmap_woman.jpeg"),units="in", width = 30, height = 10 ,  res=300)
pheatmap::pheatmap(log(normalized_matrix_heat + 1), 
                   clustering_method = 'ward.D',
                   fontsize = 8,
                   annotation_col = sample,
                   annotation_row = gen,
                   cluster_rows = F,
                   angle_col = 270, fontsize_row = 6, fontsize_col = 6)
dev.off()



check <- as.data.frame(t(normalized_matrix_heat))
check$treat <- meta_female$treatment_or_therapy



check <- pivot_longer(check, cols = unique(rownames(gen)))
check$name <- factor(check$name, levels = unique(rownames(gen)))
check <- unique(check)
check$value <- log(check$value, 2)


library(ggpubr)

g <- ggboxplot(check, x = 'treat', y = 'value', color = 'name',
               add = "jitter")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 12, face = "bold", color = "black"),
        legend.position = "none",
        axis.text.y = element_text(angle = 90, hjust = 0.5, size = 12, face = "bold", color = "black"),
        axis.text=element_text(size=14),
        axis.title=element_text(size=14,face="bold")) +
  xlab('Related genes') +
  ylab('Expression')+
  theme(strip.text.x = element_text(
    size = 14, color = "red", face = "bold.italic"
  )) +
  facet_grid(~name) 


ggsave(g, filename = paste("results/woman_genes.png", sep = ''), dpi = 300, units = "in", height = 10, width = 50, limitsize = F)


###### H2O

#man

#heatmap
h2oman <- read.csv('results/summary_man_treatment.csv')

markers_up <- h2oman$X[order(h2oman$log.FC., decreasing = T) ,drop = F][1:10]
markers_down <- h2oman$X[order(h2oman$log.FC., decreasing = F) ,drop = F][1:10]
markers_down <- markers_down[!is.na(markers_down)]

sample <- DataFrame(meta_male$treatment_or_therapy, meta_male$vital_status)
colnames(sample)[1:2] <- c('Treatment','Status')
gen <- as.data.frame(c(markers_up, markers_down))
gen$reg <- 'DOWN'
gen$reg[gen$`c(markers_up, markers_down)` %in% markers_up] <- 'UP'
rownames(gen) <- gen$`c(markers_up, markers_down)`
gen <- as.data.frame(gen[,-1, drop = F])
colnames(gen)[1] <- 'Regulation'

rownames(sample) <- meta_male$Sample.ID
sample <- as.data.frame(sample)
normalized_matrix_heat <- (normalized_matrix_male[rownames(normalized_matrix_male) %in% rownames(gen),])
jpeg(file.path("results/h20pheatmap.jpeg"),units="in", width = 30, height = 10 ,  res=300)
pheatmap::pheatmap(log(normalized_matrix_heat + 1), 
                   clustering_method = 'ward.D',
                   fontsize = 8,
                   annotation_col = sample,
                   annotation_row = gen,
                   cluster_rows = F,
                   angle_col = 270, fontsize_row = 6, fontsize_col = 6)
dev.off()



check <- as.data.frame(t(normalized_matrix_heat))
check$treat <- meta_male$treatment_or_therapy



check <- pivot_longer(check, cols = unique(rownames(gen)))
check$name <- factor(check$name, levels = unique(rownames(gen)))
check <- unique(check)
check$value <- log(check$value, 2)


library(ggpubr)

g <- ggboxplot(check, x = 'treat', y = 'value', color = 'name',
               add = "jitter")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 12, face = "bold", color = "black"),
        legend.position = "none",
        axis.text.y = element_text(angle = 90, hjust = 0.5, size = 12, face = "bold", color = "black"),
        axis.text=element_text(size=14),
        axis.title=element_text(size=14,face="bold")) +
  xlab('Related genes') +
  ylab('Expression')+
  theme(strip.text.x = element_text(
    size = 14, color = "red", face = "bold.italic"
  )) +
  facet_grid(~name) 


ggsave(g, filename = paste("results/h2pman_genes.png", sep = ''), dpi = 300, units = "in", height = 10, width = 50, limitsize = F)


#woman

#heatmap
h2owoman <- read.csv('results/summary_woman_treatment.csv')

markers_up <- h2owoman$X[order(h2owoman$log.FC., decreasing = T) ,drop = F][1:10]
markers_down <- h2owoman$X[order(h2owoman$log.FC., decreasing = F) ,drop = F][1:10]
markers_down <- markers_down[!is.na(markers_down)]

sample <- DataFrame(meta_female$treatment_or_therapy, meta_female$vital_status)
colnames(sample)[1:2] <- c('Treatment','Status')
gen <- as.data.frame(c(markers_up, markers_down))
gen$reg <- 'DOWN'
gen$reg[gen$`c(markers_up, markers_down)` %in% markers_up] <- 'UP'
rownames(gen) <- gen$`c(markers_up, markers_down)`
gen <- as.data.frame(gen[,-1, drop = F])
colnames(gen)[1] <- 'Regulation'

rownames(sample) <- meta_female$Sample.ID
sample <- as.data.frame(sample)
normalized_matrix_heat <- (normalized_matrix_female[rownames(normalized_matrix_female) %in% rownames(gen),])
jpeg(file.path("results/h20pheatmapfemale.jpeg"),units="in", width = 30, height = 10 ,  res=300)
pheatmap::pheatmap(log(normalized_matrix_heat + 1), 
                   clustering_method = 'ward.D',
                   fontsize = 8,
                   annotation_col = sample,
                   annotation_row = gen,
                   cluster_rows = F,
                   angle_col = 270, fontsize_row = 6, fontsize_col = 6)
dev.off()



check <- as.data.frame(t(normalized_matrix_heat))
check$treat <- meta_female$treatment_or_therapy



check <- pivot_longer(check, cols = unique(rownames(gen)))
check$name <- factor(check$name, levels = unique(rownames(gen)))
check <- unique(check)
check$value <- log(check$value, 2)


library(ggpubr)

g <- ggboxplot(check, x = 'treat', y = 'value', color = 'name',
               add = "jitter")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 12, face = "bold", color = "black"),
        legend.position = "none",
        axis.text.y = element_text(angle = 90, hjust = 0.5, size = 12, face = "bold", color = "black"),
        axis.text=element_text(size=14),
        axis.title=element_text(size=14,face="bold")) +
  xlab('Related genes') +
  ylab('Expression  ')+
  theme(strip.text.x = element_text(
    size = 14, color = "red", face = "bold.italic"
  )) +
  facet_grid(~name) 


ggsave(g, filename = paste("results/femLEh2pman_genes.png", sep = ''), dpi = 300, units = "in", height = 10, width = 50, limitsize = F)


