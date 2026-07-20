library(dplyr)
library(tidyr)
library(ggplot2)
library(ggsci)
library(ggpubr)
library(tidyr)
library(ggrepel)
library(DESeq2)
library(stringr)

#Figure 3A
count <- read.table("Figure3/PRJDB15555.count.txt", sep = "\t", header = T, stringsAsFactors = FALSE)
PAS<-read.table("PRJDB15555.PAS.header.bed", header = T,sep="\t")
PAS <- PAS %>% dplyr::filter(PAU>=0.05 & PAU <= 0.99 & Count >=100)
count <- count %>% filter(PAS_ID %in% PAS$PAS_ID)
rownames(count) <- count$PAS_ID
samples <-read.table("sample.txt", header = T,sep=",")
sample_id <- samples$sample_id
counts=count[,c("gene_name", "PAS_ID", sample_id)]
colnames(counts)[1] <- "gene_id"
colnames(counts)[2] <- "feature_id"
counts <- as.data.frame(counts)
PAU <- counts[,1:2]
for (sample in sample_id) {
  total_reads_per_gene <- tapply(counts[,sample], counts$gene_id, sum)
  PAU[,sample] <- counts[,sample] / total_reads_per_gene[counts$gene_id]
}
PAU <- PAU[complete.cases(PAU[,3:14]),]
transposed_combat <- t(PAU[,3:14])
pca_proc <- prcomp(transposed_combat[,apply(transposed_combat, 2, var, na.rm=TRUE) != 0],scale=TRUE,center=TRUE)
summary(pca_proc)

plotData = samples[,c("sample_id","Group")]
plotData$PC1 <- pca_proc$x[,1]
plotData$PC2 <- pca_proc$x[,2]
colors <- c("#E64B35FF","#4DBBD5FF","#00A087FF")
colors <- colors[as.numeric(as.factor(plotData$Group))]
ggplot(plotData, aes(x = PC1, y = PC2, color = Group, shape = Group)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = colors) +
  labs(
    x = "PC1 (24.2% variability)",
    y = "PC2 (15.6% variability)",
    color = "Group",
    shape = "Group"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank()
  )


#Figure 3B
WARM <- read.delim("Figure3/WARM.TUTR.txt",header = T)
rownames(WARM) <- WARM$gene_id
names(WARM) <- gsub("\\.HQ\\.bam$", "", names(WARM))
Ce <- names(WARM)[grep("_C", names(WARM))]
Te <- names(WARM)[grep("_T", names(WARM))]
Hi <- names(WARM)[grep("_H", names(WARM))]
WARM$C_WARM <-  rowMeans(WARM[,Ce],na.rm = TRUE)
WARM$T_WARM <-  rowMeans(WARM[,Te],na.rm = TRUE)
WARM$H_WARM <-  rowMeans(WARM[,Hi],na.rm = TRUE)
WARM <- WARM %>%  dplyr::select(gene_id,C_WARM,T_WARM,H_WARM) 
WARM <- na.omit(WARM)
long_data <- WARM %>%
  pivot_longer(cols = c(C_WARM, T_WARM,H_WARM), values_to = "WARM")
ggplot(long_data, aes(x = WARM, color = name)) +
  stat_ecdf(geom = "step", size = 1) +
  labs(title = "Cumulative Density Plot of gWARM values",
       x = "gWARM value",
       y = "Cumulative Density") +
  theme_minimal() +  scale_color_npg("nrc", alpha = 1) +  
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)  # Add border around the plot
  )

##Figure 3C
## Cb vs. Tc
CT <- read.table("Figure3/C_T.TUTR.diff.DRIMSeq.txt",sep="\t",header=T)
CT$gene <- sapply(strsplit(as.character(CT$gene_id), ":"), function(x) x[1])
CT$logP <- -log10(CT$adj_pvalue)
CT$type <- "NS"
CT$type[which((CT$adj_pvalue <= 0.05) & (CT$MPRO >= 0.2))] = "lengthen"
CT$type[which((CT$adj_pvalue <= 0.05) & (CT$MPRO <= -0.2))] = "shorten"
top_genes <- CT %>% dplyr::filter(adj_pvalue < 0.05 & abs(MPRO) >= 0.2) %>% 
  arrange(desc(logP)) %>% head(10)
CT_TUTR <- ggscatter(CT, x="MPRO", y="logP", color="type", 
                       palette=c("#E64B35FF","gray","#00A087FF"),size=2,
                       ylab="-log10 (adjusted P-value)",xlab="MPRO (Cb - Tc)",xlim=c(-1,1),title = "cerebellum vs temporal cortex",
                       font.label = 12)+
  geom_hline(yintercept=1.30,linetype="dashed")+
  geom_vline(xintercept=c(-0.2,0.2
  ),linetype="dashed")+
  geom_text_repel(data = top_genes, aes(label = gene), vjust = -1, hjust = 1) 

## Cb vs. Hy
CH <- read.table("Figure3/C_H.TUTR.diff.DRIMSeq.txt",sep="\t",header=T)
CH$gene <- sapply(strsplit(as.character(CH$gene_id), ":"), function(x) x[1])
CH$logP <- -log10(CH$adj_pvalue)
CH$type <- "NS"
CH$type[which((CH$adj_pvalue <= 0.05) & (CH$MPRO >= 0.2))] = "lengthen"
CH$type[which((CH$adj_pvalue <= 0.05) & (CH$MPRO <= -0.2))] = "shorten"
top_genes <- CH %>% dplyr::filter(adj_pvalue < 0.05 & abs(MPRO) >= 0.2) %>% 
  arrange(desc(logP)) %>% head(10)
CH_TUTR <- ggscatter(CH, x="MPRO", y="logP", color="type", 
                       palette=c("#E64B35FF","gray","#4DBBD5FF"),size=2,
                       ylab="-log10 (adjusted P-value)",xlab="MPRO",xlim=c(-1,1),title = "cerebellum vs hypothalamus",
                       font.label = 12)+
  geom_hline(yintercept=1.30,linetype="dashed")+
  geom_vline(xintercept=c(-0.2,0.2
  ),linetype="dashed")+
  geom_text_repel(data = top_genes, aes(label = gene), vjust = -1, hjust = 1) 

## Hy vs. Tc
HT <- read.table("Figure3/H_T.TUTR.diff.DRIMSeq.txt",sep="\t",header=T)
HT$gene <- sapply(strsplit(as.character(HT$gene_id), ":"), function(x) x[1])
HT$logP <- -log10(HT$adj_pvalue)
HT$type <- "NS"
HT$type[which((HT$adj_pvalue <= 0.05) & (HT$MPRO >= 0.2))] = "lengthen"
HT$type[which((HT$adj_pvalue <= 0.05) & (HT$MPRO <= -0.2))] = "shorten"
top_genes <- HT %>% dplyr::filter(adj_pvalue < 0.05 & abs(MPRO) >= 0.2) %>%  
  arrange(desc(logP)) %>% head(10)
HT_TUTR <- ggscatter(HT, x="MPRO", y="logP", color="type", 
                       palette=c("#4DBBD5FF","gray","#00A087FF"),size=2,
                       ylab="-log10 (adjusted P-value)",xlab="MPRO",xlim=c(-1,1),title = "hypothalamus vs temporal cortex",
                       font.label = 12)+
  geom_hline(yintercept=1.30,linetype="dashed")+
  geom_vline(xintercept=c(-0.2,0.2
  ),linetype="dashed")+
  geom_text_repel(data = top_genes, aes(label = gene), vjust = -1, hjust = 1) 

arrange2 <- ggarrange(CT_TUTR,CH_TUTR, HT_TUTR, ncol = 3, nrow = 1)

##Figure 3D
CT <- CT %>% dplyr::filter(adj_pvalue <= 0.05 & abs(MPRO) >= 0.2)
CT$gene_name <- sapply(strsplit(as.character(CT$gene_id), ":"), function(x) x[1])
CH <- CH %>% dplyr::filter(adj_pvalue <= 0.05 & abs(MPRO) >= 0.2) 
CH$gene_name <- sapply(strsplit(as.character(CH$gene_id), ":"), function(x) x[1])
HT <- HT %>% dplyr::filter(adj_pvalue <= 0.05 & abs(MPRO) >= 0.2) 
HT$gene_name <- sapply(strsplit(as.character(HT$gene_id), ":"), function(x) x[1])
diff_event <- unique(c(CT$gene_id,CH$gene_id,HT$gene_id))
WARM <- read.delim("Figure3/WARM.TUTR.txt",header = T)
rownames(WARM) <- WARM$gene_id
names(WARM) <- gsub("\\.HQ\\.bam$", "", names(WARM))
WARM_diff <- WARM[WARM$gene_id %in% diff_event, ]

WARM_diff <- WARM_diff %>% dplyr::select(N1_T,N28_T,R2_T,R8_T,N1_C,N28_C,R2_C,R8_C,N1_H,N28_H,R2_H,R8_H)
library(pheatmap)
library(RColorBrewer)
pheatmap(WARM_diff,cluster_rows=T,scale="row",cluster_col=F,gaps_col =c(4,8),show_rownames=F,show_colnames =T, 
         cellwidth = 16, cellheight = 1, filename = "different_APA.TUTR.pdf", col=rev(colorRampPalette(brewer.pal(10, "RdBu"))(20)))


##Figure 3G
count <- read.table("Figure3/OUT.gene_grouped_counts.tsv", header=TRUE, row.names=1)
sample_info <- read.table("Figure3/metadata.txt", header=TRUE,sep = "\t")
dds <- DESeqDataSetFromMatrix(countData = count,
                              colData = sample_info,
                              design = ~ condition)
smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]
dds <- DESeq(dds)
results_C_vs_T <- results(dds, contrast = c("condition", "C", "T"))
sorted_results_C_vs_T <- results_C_vs_T[order(results_C_vs_T$padj), ]
sorted_results_C_vs_T$sig <- ifelse(
  sorted_results_C_vs_T$log2FoldChange > 1 & as.numeric(sorted_results_C_vs_T$padj) < 0.05, "Up",
  ifelse(sorted_results_C_vs_T$log2FoldChange < -1 & sorted_results_C_vs_T$padj < 0.05, "Down", "nonsig")
)
sorted_results_C_vs_T$gene_id <- rownames(sorted_results_C_vs_T)
gtf_file <- "reference/gencode.v40.annotation_utr.gtf"
gtf <- read.delim(gtf_file, header = FALSE, comment.char = "#",
                  col.names = c("seqname", "source", "feature", "start", "end", 
                                "score", "strand", "frame", "attribute"))
gtf_genes <- gtf %>%
  filter(feature == "gene") %>%
  mutate(
    gene_id = str_extract(attribute, "(?<=gene_id )[^;]+") %>% str_remove_all('"'),
    gene_name = str_extract(attribute, "(?<=gene_name )[^;]+") %>% str_remove_all('"')
  ) %>%
  select(gene_id, gene_name)
sorted_results_C_vs_T <- as.data.frame(sorted_results_C_vs_T) %>% left_join(gtf_genes, by = "gene_id")
rownames(sorted_results_C_vs_T) <- sorted_results_C_vs_T$gene_id
C_vs_T_diff.sig <- sorted_results_C_vs_T %>% filter(abs(log2FoldChange) > 1 & padj < 0.05)
C_vs_T_diff.sig$group <- "C_vs_T"


## correlation between APA change and gene expression
CT <- read.table("C_T.gene.diff.DRIMSeq.txt",sep="\t",header=T)
CT$type <- ifelse(
  (abs(as.numeric(CT$MPRO)) > 0.1 & as.numeric(CT$adj_pvalue) < 0.05), "DAG", "nonDAG")
CT$gene_name <- CT$gene_id
CT$gene_name <- sapply(strsplit(as.character(CT$gene_id), ":"), function(x) x[1])

CT <- CT %>% select(gene_name, MPRO,sig, type) 
sorted_results_C_vs_T$sig <- ifelse(abs(sorted_results_C_vs_T$log2FoldChange) > 1 & as.numeric(sorted_results_C_vs_T$padj) < 0.05, "DEG", "nonDEG")
CT_exp <- sorted_results_C_vs_T  %>% select (gene_name, log2FoldChange, sig)
interg_CT <- CT %>% left_join(CT_exp, by = "gene_name")
interg_CT <- na.omit(interg_CT)
interg_CT$label <- paste0(interg_CT$type,"_",interg_CT$sig.y)

ggscatter(interg_CT, x="MPRO", y="log2FoldChange", color="label", alpha = 0.7,
          size=2,
          ylab="log2 fold change (Cb/Tc)",xlab="MPRO (Cb-Tc)",title = "Cb vs Tc",
          font.label = 12)+
  geom_hline(yintercept=c(-1,1),linetype="dashed")+
  stat_cor(method = "pearson", label.x = 0.2, label.y = 6) +
  geom_vline(xintercept=c(-0.2,0.2
  ),linetype="dashed") + 
  scale_color_manual(values = c("DAG_DEG" = "#DC0000",  
                                "DAG_nonDEG" = "#4DBBD5","nonDAG_nonDEG" = "grey", 
                                "nonDAG_DEG" = "#00A087"))+
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1)) 


