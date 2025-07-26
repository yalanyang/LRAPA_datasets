library(tidyr)
library(stringr)
library(tibble)
library(dplyr)
library(ggpubr)
library(GenomicRanges)
library(pheatmap)
## Figure 4A
setwd("/Users/yangyalan/OneDrive - The University of Chicago/Chicago/Long-read-APA/delatorre_human_scisorseq/isoseurat/Diff/gene/")
listdir <- list.files(path = "Figure4A/",pattern="*.diff.betweenPAS_gene_level.txt") 

diff_DAG <- list()

for (k in 1:length(listdir)){
  dag<- read.table(paste0(listdir[k]),header=T,sep='\t')
  group1 <- strsplit(listdir[k],"[.]")[[1]][1]
  dag <- dag %>% dplyr::filter(abs(dgDPAU) > 0.2 & adj_pvalue < 0.05)
  dag$group <- group1
  diff_DAG[[k]] <- dag
}

results_df <- do.call(rbind, diff_DAG)
results_df <- results_df %>% select(gene_name,group)
binary_matrix <- results_df %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = group, values_from = value, values_fill = 0) %>%
  column_to_rownames("gene_name")
overlap_matrix <- as.matrix(t(binary_matrix)) %*% as.matrix(binary_matrix)

log_overlap_matrix <- log2(overlap_matrix)  

log_overlap_matrix[lower.tri(log_overlap_matrix)] <- NA  

order <- c("Ex_Pro", "Ex_In", "In_Pro", "ExDp_ExNExM", "ExM_ExDp","ExM_ExMU", "ExM_ExNExM", "ExMU_ExDp", "ExMU_ExNExM",
           "ExN1_ExDp", "ExN1_ExM", "ExN1_ExMU", "ExN1_ExN2", "ExN1_ExN3", "ExN1_ExNExM",
           "ExN2_ExDp", "ExN2_ExM",  "ExN2_ExMU", "ExN2_ExN3", "ExN2_ExNExM", "ExN3_ExDp", "ExN3_ExM",
           "ExN3_ExMU", "ExN3_ExNExM",
           "IP_oRG", "IP_PgG2M", "IP_PgS", "IP_vRG","oRG_PgG2M", "oRG_PgS", "oRG_vRG", "PgG2M_PgS",
           "PgG2M_vRG", "PgS_vRG")
overlap_matrix <- overlap_matrix[order, order]
log_overlap_matrix <- overlap_matrix[order, order]
pheatmap(log_overlap_matrix,
         scale = 'none',
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         gaps_row=c(3,24),
         gaps_col=c(3,24),
         legend = TRUE, 
         display_numbers = round(overlap_matrix, 2),
         color = colorRampPalette(c("#b7d1ac", "white", "firebrick3"))(10),
         number_color = "black",
         na_col = "white")

write.table(overlap_matrix, "overlap_matrix.txt", sep = "\t", row.names = T, quote = FALSE)

## Figure 4B
fetal <- read.table("Figure4A/Ex_Pro.diff.betweenPAS_gene_level.txt",sep="\t",header=T)
fetal$logP <- -log10(fetal$adj_pvalue)
fetal$type <- "NS"
fetal$type[which((fetal$adj_pvalue < 0.05) & (fetal$dgDPAU > 0.2) )] = "Ex > Pro"
fetal$type[which((fetal$adj_pvalue < 0.05) & (fetal$dgDPAU < -0.2))] = "Pro > Ex"
top_genes <- fetal %>% dplyr::filter(adj_pvalue <= 0.05 & abs(fetal$dgDPAU) > 0.2) %>% 
  arrange(desc(logP)) %>% head(20)
ggscatter(fetal, x="dgDPAU", y="logP", color="type", 
          palette=c("#E64B35FF","gray","#00A087FF"),size=1,
          ylab="-log10 (adjusted P-value)",xlab="dgDPAU (Ex-Pro)",xlim=c(-1,1),title = "Ex vs In",
          font.label = 12)+
  geom_hline(yintercept=1.30,linetype="dashed")+
  geom_vline(xintercept=c(-0.1,0.1
  ),linetype="dashed")+
  geom_text_repel(data = top_genes, aes(label = gene_name), vjust = -1, hjust = 1) 



#Figure 4D
ex.length <- read.table("Figure4/excitatory_lengthening_GO.txt",sep="\t",header=T)
ex.length$log <- -log10(ex.length$FDR)
ggbarplot(ex.length[1:5,], x = "Term", y = "log",
          fill = "#E64B35FF", 
          color = "white",
          ylab = "-log10 (P_value)", xlab="GO Biologcial Process",
          rotate = TRUE,
          ggtheme =  theme_pubr(),sort.val = "asc"
) + theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1)) 

#Figure 4E
ex.short <- read.table("Figure4/excitatory_shortening_GO.txt",sep="\t",header=T)
ex.short$log <- -log10(ex.short$FDR)
ggbarplot(ex.short[1:5,], x = "Term", y = "log",
          fill = "#00A087FF", 
          color = "white",
          ylab = "-log10 (P_value)", xlab="GO Biologcial Process",
          rotate = TRUE,
          ggtheme =  theme_pubr(),sort.val = "asc"
) + theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1)) 

#Figure 4F
bulk <- read.table("Figure4/CP_GZ.gene.diff.DRIMSeq.txt",sep="\t",header=T)
bulk <- bulk %>% dplyr::filter(abs(dgDPAU) > 0.1 & adj_pvalue < 0.05)
sciso <- read.table("fetal.brain.sc.pseudobulk.all.dpau_matrix.txt",sep="\t",header=T)
sciso <- sciso %>% filter(gene_id %in% bulk$gene_id)
bulk <- bulk %>% select(gene_id, gDPAU1, gDPAU2, dgDPAU,Note)
colnames(bulk) <- c("gene_id", "CP", "GZ", "dgDPAU_bulk","Type")
sciso <- sciso %>% left_join(bulk, by="gene_id")

sciso$dgDPAU_sc <- sciso$Ex-sciso$Pro
ggscatter(sciso, y = "dgDPAU_sc", x = "dgDPAU_bulk",fill="Type", add = "reg.line",add.params = list(color = "#4477AA",fill = "lightgray"),
          xlab="dgDPAU (CP-GZ)",ylab="dgDPAU (Excitatory-Progenitor)",color="Type",font.label = 12)+ 
  stat_cor(method = "pearson", label.x = 0, label.y = -0.5) +scale_color_manual(values = c("#EE667799", "#4477AA99"))+
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1)) 


