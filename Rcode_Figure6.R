library(dplyr)
library(tidyr)
library(ggplot2)
library(ggsci)
library(tibble)
library(stringr)
library(pheatmap)
library(RColorBrewer)

#Figure 6B
count <- read.table("Figure6/iNGN.count.txt",sep="\t",header=T)
sum_counts <- rowSums(count[, 2:8])
count <- count[sum_counts > 40, ]
samples <- read.table("Figure6/sample.txt", sep = ",", header = T, stringsAsFactors = FALSE)
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
PAU <- PAU[complete.cases(PAU[,3:9]),]

PAU <- PAU %>%
  mutate(mean_value = rowMeans(select(., -gene_id, -feature_id), na.rm = TRUE)) %>%
  filter(mean_value >= 0.05 & mean_value <= 0.95) %>%
  select(-mean_value)

transposed_combat <- t(PAU[,3:9])
pca_proc <- prcomp(transposed_combat[,apply(transposed_combat, 2, var, na.rm=TRUE) != 0],scale=TRUE,center=TRUE)
summary(pca_proc)
plotData = samples[,c("sample_id","Group")]
plotData$PC1 <- pca_proc$x[,1]
plotData$PC2 <- pca_proc$x[,2]
pdf("PCA12_based_PA_count.pdf", 4,3.5)
ggplot(plotData, aes(PC1, PC2, color = Group)) +
  geom_point(size = 5) +
  scale_color_manual(values = c("#E64B35", "#4DBBD5", "#00A087")) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 13),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    panel.grid.major = element_line(color = "grey85"),
    panel.grid.minor = element_blank()
  )
dev.off()

#Figure 6C
WARM <- read.delim("Figure6/WARM.TUTR.txt",header = T)
rownames(WARM) <- WARM$gene_id
WARM$D0_WARM <-  rowMeans(WARM[,2:3],na.rm = TRUE)
WARM$D7_WARM <-  rowMeans(WARM[,4:5],na.rm = TRUE)
WARM$D12_WARM <-  rowMeans(WARM[,6:8],na.rm = TRUE)
WARM <- WARM %>%  dplyr::select(gene_id,D0_WARM,D7_WARM, D12_WARM) 
WARM <- na.omit(WARM)
long_data <- WARM %>%
  pivot_longer(cols = c(D0_WARM, D7_WARM,D12_WARM), values_to = "WARM")

ggplot(long_data, aes(x = WARM, color = name)) +
  stat_ecdf(geom = "step", size = 1) +
  labs(title = "Cumulative Density Plot of WARM values",
       x = "WARM",
       y = "Cumulative Density") +
  theme_minimal() + 
  scale_color_manual(values = c("D0_WARM" = nrc_colors[1], 
                                "D7_WARM" = nrc_colors[3], 
                                "D12_WARM" = nrc_colors[2])) + 
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)  # Add border around the plot
  )


#Figure 6E
D7_D0 <- read.table("Figure6/D7_D0.TUTR.diff.DRIMSeq.txt",sep="\t",header=T)
D12_D0 <- read.table("Figure6/D12_D0.TUTR.diff.DRIMSeq.txt",sep="\t",header=T)
D12_D7 <- read.table("Figure6/D12_D7.TUTR.diff.DRIMSeq.txt",sep="\t",header=T)
D7_D0 <- D7_D0 %>% dplyr::filter(sig == "TRUE")
D12_D0 <- D12_D0 %>% dplyr::filter(sig == "TRUE") 
D12_D7 <- D12_D7 %>% dplyr::filter(sig == "TRUE") 

diff_gene <- unique(c(D7_D0$gene_id,D12_D0$gene_id,D12_D7$gene_id))
WARM_diff <- WARM[WARM$gene_id %in% diff_gene, ]
pheatmap(WARM_diff[,2:8],cluster_rows=T,scale="row",cluster_col=F,gaps_col =c(2,4),show_rownames=F,show_colnames =T, 
         cellwidth = 16, cellheight = 0.2, filename = "different_APA.TUTR.pdf",col=rev(colorRampPalette(brewer.pal(10, "RdBu"))(20)))


#Figure 6G
## SE-TSS coupling
all <- read.table("Figure6/iNGN_merge_exon-PA.coordination.chiqtest.txt", sep = "\t", header = T, stringsAsFactors = FALSE)
all <- all %>% dplyr::filter(FDR < 0.05)
human_fetal_SE <- read.table("Figure5/human_fetal_SE_PA_coupling_0622.txt", sep = "\t", header = T, stringsAsFactors = FALSE)
human_fetal_SE_sig  <- human_fetal_SE %>% dplyr::filter(FDR < 0.05)
human_adult_SE  <- read.table("Figure5/human_encode_SE_PA_coupling_0622.txt", sep = "\t", header = T, stringsAsFactors = FALSE)
human_adult_SE_sig  <- human_adult_SE %>% dplyr::filter(FDR < 0.05)
x3 <- list(iNGN = unique(all$gene_id), fetal = unique(human_fetal_SE_sig$gene_id), adult = unique(human_adult_SE_sig$gene_id))
plot(Venn(x3)) 

## PA-TSS coupling
all <- read.table("Figure6/iNGN_merge_TSS-PA.coordination.chiqtest.txt", sep = "\t", header = T, stringsAsFactors = FALSE)
all <- all %>% dplyr::filter(FDR < 0.05)
human_fetal_TSS <- read.table("Figure5/human_fetal_TSS_PA_coupling_0622.txt", sep = "\t", header = T, stringsAsFactors = FALSE)
human_fetal_TSS_sig <- human_fetal_TSS %>% dplyr::filter(FDR < 0.05)
human_adult_TSS <- read.table("Figure5/human_encode_TSS_PA_coupling_0622.txt", sep = "\t", header = T, stringsAsFactors = FALSE)
human_adult_TSS_sig <- human_adult_TSS %>% dplyr::filter(FDR < 0.05)
x1 <- list(iNGN = all$gene_id, fetal = unique(human_fetal_TSS_sig$gene_id), adult = unique(human_adult_TSS_sig$gene_id))
pdf("Venn.coupling_TSS.gene_level.pdf", 4,4)
plot(Venn(x1)) 

