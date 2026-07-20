library(tidyr)
library(stringr)
library(tibble)
library(dplyr)
library(ggpubr)
library(pheatmap)
library(ggrepel)
library(ggvenn)
## Figure 4A
listdir <- list.files(path = "Figure4/gene/",pattern="*.gene.diff.txt") 
diff_DAG <- list()
for (k in 1:length(listdir)){
  dag<- read.table(paste0(listdir[k]),header=T,sep='\t')
  group1 <- strsplit(listdir[k],"[.]")[[1]][1]
  dag <- dag %>% dplyr::filter(abs(MPRO) >= 0.4 & adj_pvalue <= 0.05) %>% select(gene_id, pvalue, adj_pvalue,dWARM, MPRO) %>%
    mutate(group = group1)
  diff_DAG[[k]] <- dag
}
results_df <- do.call(rbind, diff_DAG)
results_df <- results_df %>% select(gene_id,group)
binary_matrix <- results_df %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = group, values_from = value, values_fill = 0) %>%
  column_to_rownames("gene_id")
overlap_matrix <- as.matrix(t(binary_matrix)) %*% as.matrix(binary_matrix)
log_overlap_matrix <- log2(overlap_matrix)  
log_overlap_matrix[lower.tri(log_overlap_matrix)] <- NA  
order <- c("Ex_Pro", "Ex_In", "In_Pro", "ExDp_ExM", "ExDp_ExMU", "ExDp_ExN1", "ExDp_ExN2", "ExDp_ExN3", "ExDp_ExNExM", "ExM_ExMU",  
           "ExM_ExN1", "ExM_ExN2", "ExM_ExN3", "ExMU_ExN1", "ExMU_ExN2", "ExMU_ExN3", "ExN1_ExN2", "ExN1_ExN3", "ExN2_ExN3", 
          "ExNExM_ExM", "ExNExM_ExMU", "ExNExM_ExN1", "ExNExM_ExN2", "ExNExM_ExN3", "IP_oRG", "IP_PgG2M", "IP_PgS",    
         "IP_vRG","oRG_PgG2M", "oRG_PgS", "oRG_vRG", "PgG2M_PgS", "PgG2M_vRG", "PgS_vRG" )
overlap_matrix <- overlap_matrix[order, order]
log_overlap_matrix <- overlap_matrix[order, order]
pheatmap(log2(log_overlap_matrix),
         scale = 'none',
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         gaps_row=c(3,24),
         gaps_col=c(3,24),
         legend = TRUE,  
         display_numbers = round(overlap_matrix, 2),
         color = colorRampPalette(c("#b7d1ac", "white", "firebrick3"))(10),
         number_color = "black",
         na_col = "white",
         cellwidth = 14,cellheight = 14,
         filename = "Figure4A.DAG_overlap.pdf"
) 

## Figure 4B
fetal <- read.table("Figure4/gene/Ex_Pro.gene.diff.txt",sep="\t",header=T)
fetal$logP <- -log10(fetal$adj_pvalue)
fetal$type <- "NS"
fetal$type[which((fetal$adj_pvalue <= 0.05) & (fetal$MPRO >= 0.4) )] = "Ex > Pro"
fetal$type[which((fetal$adj_pvalue <= 0.05) & (fetal$MPRO <= -0.4))] = "Pro > Ex"
top_genes <- fetal %>% dplyr::filter(adj_pvalue <= 0.05 & abs(fetal$MPRO) >= 0.4) %>% 
  arrange(desc(logP)) %>% head(10)
ggscatter(fetal, x="MPRO", y="logP", color="type", 
          palette=c("#E64B35FF","gray","#00A087FF"),size=1,
          ylab="-log10 (adjusted P-value)",xlab="MPRO (Ex-Pro)",title = "Ex vs Pro",
          font.label = 12)+
  geom_hline(yintercept=1.30,linetype="dashed")+
  geom_vline(xintercept=c(-0.2,0.2
  ),linetype="dashed")+
  geom_text_repel(data = top_genes, aes(label = gene_id), vjust = -1, hjust = 1) 


##Figure4D
fetal <- read.table("Figure4/CP_GZ.TUTR.diff.DRIMSeq.txt",sep="\t",header=T)
fetal$logP <- -log10(fetal$adj_pvalue)
fetal$type <- "NS"
fetal$type[which((fetal$adj_pvalue <= 0.05) & (fetal$MPRO >= 0.2))] = "CP > GZ"
fetal$type[which((fetal$adj_pvalue <= 0.05) & (fetal$MPRO <= -0.2))] = "GZ > CP"
fetal$gene_name <- sapply(strsplit(as.character(fetal$gene_id), ":"), function(x) x[1])
top_genes <- fetal %>% dplyr::filter(adj_pvalue <= 0.05) %>% 
  arrange(desc(logP)) %>% head(10)
TUTR <- ggscatter(fetal, x="MPRO", y="logP", color="type", 
                  palette=c("#EE667799","#4477AA99","gray"),size=2,
                  ylab="-log10 (adjusted P-value)",xlab="MPRO  (CP-GZ)",
                  font.label = 12)+
  geom_hline(yintercept=1.30,linetype="dashed")+
  geom_vline(xintercept=c(-0.2,0.2
  ),linetype="dashed")+
  geom_text_repel(data = top_genes, aes(label = gene_name), vjust = -1, hjust = 1) 

#Figure 4E
DAG <- read.table("Figure4/CP_GZ.gene.diff.DRIMSeq.txt", sep = "\t", header = T, stringsAsFactors = FALSE)
DAG <- DAG %>% dplyr::filter(abs(MPRO)>= 0.2 & adj_pvalue < 0.05)
DEG <- read.table("Figure4/Exp/DEG_CP_vs_GZ.csv", sep = ",", header = T, stringsAsFactors = FALSE)
DEG <- DEG %>% filter(abs(log2FoldChange) > 1 & padj < 0.05)
write.table(unique(DEG$gene_name), "Figure4E.DEG.txt", sep = "\t", row.names = F, quote = FALSE)

DET <- read.table("Figure4/Exp/DET_CP_vs_GZ.csv", sep = ",", header = T, stringsAsFactors = FALSE)
DET <- DET %>% filter(abs(log2FoldChange) > 1 & padj < 0.05)
write.table(unique(DET$gene_name), "Figure4E.DET.txt", sep = "\t", row.names = F, quote = FALSE)

DIU <- read.table("Figure4/Exp/DIU_CP_vs_GZ.csv", sep = "\t", header = T, stringsAsFactors = FALSE)
DIU <- DIU %>% filter(adj_pvalue < 0.05)
listInput <- list(DEG=unique(DEG$gene_name),DET=unique(DET$gene_name),DIU=unique(DIU$gene_id),DAG=unique(DAG$gene_id))
ggvenn(
  listInput, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4
)

#Figure 4F
bulk <- read.table("Figure4/CP_GZ.gene.diff.DRIMSeq.txt",sep="\t",header=T)
bulk <- bulk %>% dplyr::filter(abs(MPRO) >= 0.2 & adj_pvalue <= 0.05)
bulk <- bulk %>% select(gene_id, MPRO,Note)
colnames(bulk) <- c("gene_id","MPRO_bulk","type")
sciso <- read.table("Figure4/gene/Ex_Pro.gene.diff.txt",sep="\t",header=T)
sciso <- sciso %>% select(gene_id, MPRO,Note)
sciso <- sciso %>% filter(gene_id %in% bulk$gene_id)
colnames(sciso) <- c("gene_id","MPRO_sc","Note")
sciso <- sciso %>% left_join(bulk, by="gene_id")
ggscatter(sciso, y = "MPRO_sc", x = "MPRO_bulk",fill="type", add = "reg.line",add.params = list(color = "#4477AA",fill = "lightgray"),
          xlab="MPRO (CP-GZ)",ylab="MPRO (Excitatory-Progenitor)",color="type",font.label = 12)+ 
  stat_cor(method = "pearson", label.x = -1, label.y = 1.5) +scale_color_manual(values = c("#EE667799", "#4477AA99"))+
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1)) 

