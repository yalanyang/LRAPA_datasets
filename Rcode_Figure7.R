library(dplyr)
library(readr)
library(stringr)
library(ggpubr)
library(tidyr)
library(data.table)
library(ggsci)
library(pheatmap)
library(tidyverse)
library(readxl)

##Figure 7A
file_list <- list.files(path = "Figure7/count", pattern = "\\.count$")
process_file <- function(file) {
  df <- read_tsv(file, col_types = cols())
  df <- df %>%
    mutate(ID = str_c(gene_id, TSS, PAS, sep = "_")) %>%
    select(ID, read_counts)
  names <- sub("\\.TSS_PA.count$", "", file)
  colnames(df) <- c("ID",names)
  return(df)
}

data_list <- lapply(file_list, process_file)
merged_data <- Reduce(function(x, y) full_join(x, y, by = "ID"), data_list)

merged_data <- merged_data %>%
  mutate(gene_name = str_extract(ID, "^[^_]+"))

merged_data <- merged_data %>% select(ID, gene_name, vRG,oRG,PgS,PgG2M, IP, `ExN-1`, `ExN-2`, `ExN-3`, `ExN-ExM`, `vRG-ExN`, ExDp, ExM, `ExM-U`, In, -Mic_Per_End) 
merged_data2 <- merged_data[apply(merged_data[,c(3:16)], 1, function(x) max(x, na.rm = TRUE)) >= 5, ]
merged_data2 <- merged_data2[apply(merged_data2[,c(3:16)], 1, function(x) sum(x, na.rm = TRUE)) >= 10, ]
merged_data2$transcript_count <- rowSums(merged_data2[, 3:16], na.rm = TRUE)

##normalized transcript count
numeric_cols <- merged_data2 %>% select(-ID, -gene_name,-transcript_count)
col_sums <- colSums(numeric_cols, na.rm=TRUE)
cpm_transcript <- (numeric_cols / col_sums) * 1e6
row_sums <- rowSums(cpm_transcript, na.rm=TRUE)
row_normalized_transcript <- cpm_transcript / row_sums
row_normalized_transcript[is.na(row_normalized_transcript)] <- 0
final_transcript <- cbind(merged_data2[, c("ID", "gene_name","transcript_count")], row_normalized_transcript)


##normalized gene count
summed_gene <- merged_data %>% select(-ID) %>%
  group_by(gene_name) %>%
  summarise(across(where(is.numeric), sum, na.rm=TRUE)) %>%
  ungroup()
numeric_cols <- summed_gene %>% select(-gene_name)
col_sums <- colSums(numeric_cols, na.rm=TRUE)
cpm_gene <- (numeric_cols / col_sums) * 1e6
row_sums <- rowSums(cpm_gene, na.rm=TRUE)
gene_normalized <- cpm_gene / row_sums
gene_normalized[is.na(gene_normalized)] <- 0
summed_gene$sum_gene <- rowSums(summed_gene[, 2:15], na.rm = TRUE)
final_gene <- cbind(gene_name = summed_gene$gene_name, gene_count = summed_gene$sum_gene, gene_normalized)
merged_df <- inner_join(final_transcript, final_gene, by = "gene_name", suffix = c("_transcript", "_gene"))


cor_results <- apply(merged_df[,c(4:32)], 1, function(x) {
  transcript_expr <- x[1:14]
  gene_expr <- x[16:29]
  cor(transcript_expr, gene_expr, method = "spearman", use = "complete.obs")
})

p_values <- apply(merged_df[,c(4:32)], 1, function(x) {
  transcript_expr <- x[1:14]
  gene_expr <- x[16:29]
  #t.test(transcript_expr, gene_expr, paired = TRUE)$p.value
  wilcox.test(transcript_expr, gene_expr, paired = F)$p.value
})

adjusted_p_values <- p.adjust(p_values, method = "fdr")

results_df <- data.frame(
  correlation = cor_results,
  p_value = p_values,
  adjusted_p_value = adjusted_p_values
)

results_df <- cbind(merged_df,results_df)
results_df <- results_df[order(results_df$p_value), ]

##coupling_genes
sc_coupling <- read.table("Figure5/fetal_sc.tss-pas.chiqtest.txt", sep = "\t", header = T, stringsAsFactors = FALSE)
sig_coupling <- sc_coupling %>% dplyr::filter(FDR < 0.05)
non_sig_coupling <- sc_coupling %>% dplyr::filter(FDR >= 0.05)

##select the TSS-PA with the smallest P value.
results_df_unique <- results_df %>%
  group_by(gene_name) %>%
  slice_min(order_by = p_value, n = 1, with_ties = FALSE) %>%
  ungroup()

results_df_unique <- results_df_unique[order(results_df_unique$p_value), ]
results_df_sig <- results_df_unique[results_df_unique$gene_name %in% sig_coupling$gene_id,]
results_df_nosig <- results_df_unique[results_df_unique$gene_name %in% non_sig_coupling$gene_id,]

results_df_sig$group <- "sig"
results_df_nosig$group <- "nonsig"
stat <- rbind(results_df_sig,results_df_nosig)
my_comparisons <- list(c("sig", "nonsig"))

corr<- ggboxplot(stat, x = "group", y = "correlation", 
                 ylab = "Spearman correlation coefficient", 
                 color = "group",
                 add.params = list(size = 1, jitter = 0.2)) + stat_compare_means(comparisons = my_comparisons) + rotate_x_text(45)

pvalue <- ggboxplot(stat, x = "group", y = "p_value", 
                    ylab = " Wilcoxon test P-value", 
                    color = "group",
                    add.params = list(size = 1, jitter = 0.2)) + stat_compare_means(comparisons = my_comparisons) + rotate_x_text(45)

ggarrange(corr, pvalue, ncol = 2, nrow = 1)


##Figure 7B
gene <- "TPM1"
transcript_exp <- results_df[results_df$gene_name == gene, c(1,3:17)] 
transcript_exp <- transcript_exp %>% dplyr::rename(count=transcript_count)
colnames(transcript_exp) <- gsub("_transcript", "", colnames(transcript_exp))
gene_exp <- results_df[results_df$gene_name == gene, c(2,18:32)][1,] %>%
  dplyr::rename(ID = gene_name,count=gene_count)
colnames(gene_exp) <- gsub("_gene", "", colnames(gene_exp))
combined_exp <- bind_rows(transcript_exp, gene_exp)
rownames(combined_exp) <- paste0(combined_exp$ID," (",combined_exp$count,")")
annotation_col <- data.frame(
  Category = c("RG", "RG", "MP", "MP", "IP", 
               "EN", "EN", "EN", "EN", "EN", "EN", "EN", "EN",
               "IN")
)
annotation_colors <- list(
  Category = c("RG" = "#E64B35",   # Red
               "MP" = "#4DBBD5",   # Blue
               "IP" = "#00A087",   # Teal
               "EN" = "#F39B7F",   # Orange
               "IN" = "#3C5488"   # Dark Blue
  )
)
rownames(annotation_col) <- colnames(combined_exp[,c(3:16)])


pheatmap(as.matrix(combined_exp[,c(3:16)]),
         color = viridis::viridis(50),
         cluster_rows = FALSE, 
         cluster_cols = FALSE,  
         show_rownames = TRUE, 
         show_colnames = TRUE,
         main = as.character(gene),annotation_col = annotation_col,  
         annotation_colors = annotation_colors)



##Figure 7D
gene <- "SCN2A"
transcript_exp <- results_df[results_df$gene_name == gene, c(1,3:17)] 
transcript_exp <- transcript_exp %>% dplyr::rename(count=transcript_count)
colnames(transcript_exp) <- gsub("_transcript", "", colnames(transcript_exp))

gene_exp <- results_df[results_df$gene_name == gene, c(2,18:32)][1,] %>%
  dplyr::rename(ID = gene_name,count=gene_count)
colnames(gene_exp) <- gsub("_gene", "", colnames(gene_exp))
combined_exp <- bind_rows(transcript_exp, gene_exp)
write.table(combined_exp,"combined_exp.TPM1.txt",sep="\t",quote=F)

rownames(combined_exp) <- paste0(combined_exp$ID," (",combined_exp$count,")")

annotation_col <- data.frame(
  Category = c("RG", "RG", "MP", "MP", "IP", 
               "EN", "EN", "EN", "EN", "EN", "EN", "EN", "EN",
               "IN")
)
annotation_colors <- list(
  Category = c("RG" = "#E64B35",   # Red
               "MP" = "#4DBBD5",   # Blue
               "IP" = "#00A087",   # Teal
               "EN" = "#F39B7F",   # Orange
               "IN" = "#3C5488"   # Dark Blue
  )
)
rownames(annotation_col) <- colnames(combined_exp[,c(3:16)])

pheatmap(as.matrix(combined_exp[,c(3:16)]),
         color = viridis::viridis(50),
         cluster_rows = FALSE, 
         cluster_cols = FALSE,  
         show_rownames = TRUE, 
         show_colnames = TRUE,
         main = as.character(gene),annotation_col = annotation_col,  
         annotation_colors = annotation_colors)


TFregulon <- read.table("Figure7/TF_regulon_2019_neuron.txt", header=TRUE, fill=TRUE)
TFregulon_list <- lapply(TFregulon, function(x) x[x != "BLANK"])


perform_fisher_test <- function(disease_genes, sig_genes, nonsig_genes) {
  in_sig <- sum(disease_genes %in% sig_genes)
  in_nonsig <- sum(disease_genes %in% nonsig_genes)
  not_in_sig <- length(sig_genes) - in_sig
  not_in_nonsig <- length(nonsig_genes) - in_nonsig

  contingency_table <- matrix(c(in_sig, in_nonsig, not_in_sig, not_in_nonsig),
                              nrow = 2, 
                              dimnames = list(c("In TF", "Not in TF"),
                                              c("adult_tss_sig", "adult_tss_nonsig")))
  fisher_result <- fisher.test(contingency_table)
  p_value <- fisher_result$p.value
  odds_ratio <- fisher_result$estimate
  conf_int <- fisher_result$conf.int

  return(c(in_sig=in_sig, in_nonsig=in_nonsig, not_in_sig=not_in_sig, not_in_nonsig=not_in_nonsig, p_value = p_value, odds_ratio = odds_ratio, conf_lower = conf_int[1], conf_upper = conf_int[2]))
}

fetal_tss <- read.table("Figure5/human_fetal_TSS-PAS.coordination.chiqtest.txt", sep = "\t", header = T, stringsAsFactors = FALSE)
fetal_tss$type <- "fetal_tss_sig"
fetal_tss$type[which(fetal_tss$FDR >= 0.05)] = "fetal_tss_nonsig"
fetal_tss <- fetal_tss %>% select("gene_id","type") %>% unique()


tf_fetal_tss <- lapply(names(TFregulon_list), function(d) {
 perform_fisher_test(TFregulon_list[[d]], fetal_tss[fetal_tss$type=="fetal_tss_sig",]$gene_id, fetal_tss[fetal_tss$type=="fetal_tss_nonsig",]$gene_id)
})
tf_fetal_tss <- as.data.frame(do.call(rbind, tf_fetal_tss), row.names = names(TFregulon_list))
tf_fetal_tss$TF <- rownames(tf_fetal_tss)
tf_fetal_tss$fdr = p.adjust(as.numeric(tf_fetal_tss$p_value),'fdr')
tf_fetal_tss  <- tf_fetal_tss %>% dplyr::filter(fdr < 0.05)

tf_fetal_tss_sig <- tf_fetal_tss %>% inner_join(cellTFregu, by="TF")
tf_fetal_tss_sig <- tf_fetal_tss_sig[order(tf_fetal_tss_sig$fdr, decreasing = FALSE), ] %>% mutate(TF = factor(TF, levels = unique(TF))) 

ggplot(tf_fetal_tss_sig, aes(x = TF, y = -log10(fdr), fill = `odds_ratio.odds ratio`)) + geom_bar(stat = 'identity', position = position_dodge2()) +
  theme_bw() +
  geom_hline(yintercept = 1.30, lty = 'dashed', size = 0.5, color = 'red') +
  labs(y = '-log10 (adjusted P value)', x = '') +
  # Nature-style color scale (muted tones)
  scale_fill_gradientn(colors = c("#1f78b4", "#a6cee3", "#33a02c", "#fb9a99", "#e31a1c"),
                       name = "Odds Ratio") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 10),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        legend.key.size = unit(0.3, 'cm'),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        plot.margin = unit(c(5, 5, 0, 5), "pt"),
        legend.position = "right",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black", size = 0.2))


##Figure 7F
encodeRBPs = read.csv("reference/RBP_Data/RBP_targets_ENCODE.csv", header=TRUE);
encodeRBPs = encodeRBPs %>% filter(cell.type=="HepG2") 
encodeRBPs = encodeRBPs %>% select(hgnc_symbol,RBP)
encodeRBPs$source <- paste0(encodeRBPs$RBP,"_encode")
brainRBPs = read.csv("reference/RBP_Data/RBP_targets_v5.csv", header=TRUE);
brainRBPs = dplyr::select(brainRBPs, -c(MGI.symbol, ENSMUSG))
brainRBPs = brainRBPs %>% select(HGNC.symbol,RBP)
colnames(brainRBPs) <- c("hgnc_symbol","RBP")
FMRP <- read.table('reference/genelist/FMR_targets.txt',header=TRUE,sep='\t')
FMRP$RBP <- "FMRP"
colnames(FMRP) <- c("hgnc_symbol","RBP")
brainRBPs <- rbind(brainRBPs,FMRP)
brainRBPs$source <- paste0(brainRBPs$RBP,"_brain")
RBPs <- rbind(encodeRBPs,brainRBPs)
RBPs_list <- split(RBPs$hgnc_symbol, RBPs$source)

rbp_fetal_tss <- lapply(names(RBPs_list), function(d) {
 perform_fisher_test(RBPs_list[[d]], fetal_tss[fetal_tss$type=="fetal_tss_sig",]$gene_id, fetal_tss[fetal_tss$type=="fetal_tss_nonsig",]$gene_id)
})
rbp_fetal_tss <- as.data.frame(do.call(rbind, rbp_fetal_tss), row.names = names(RBPs_list))
rbp_fetal_tss$id <- rownames(rbp_fetal_tss)
rbp_fetal_tss <- rbp_fetal_tss %>%
  separate(id, into = c("RBP", "source"), sep = "_")
rbp_fetal_tss$fdr = p.adjust(as.numeric(rbp_fetal_tss$p_value),'fdr')
rbp_fetal_tss  <- rbp_fetal_tss %>% dplyr::filter(fdr < 0.05)
rbp_fetal_tss <- rbp_fetal_tss[order(rbp_fetal_tss$fdr, decreasing = FALSE), ] %>% mutate(RBP = factor(RBP, levels = unique(RBP))) 

ggplot(rbp_fetal_tss, aes(x = RBP, y = -log10(fdr), fill = `odds_ratio.odds ratio`)) + geom_bar(stat = 'identity', position = position_dodge2()) +
  theme_bw() +
  geom_hline(yintercept = 1.30, lty = 'dashed', size = 0.5, color = 'red') +
  labs(y = '-log10 (adjusted P value)', x = '') +
  # Nature-style color scale (muted tones)
  scale_fill_gradientn(colors = c("#1f78b4", "#a6cee3", "#33a02c", "#fb9a99", "#e31a1c"),
                       name = "Odds Ratio") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 10),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        legend.key.size = unit(0.3, 'cm'),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        plot.margin = unit(c(5, 5, 0, 5), "pt"),
        legend.position = "right",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black", size = 0.2))




##Figure 7G
matching_RBP_names <- match(rownames(rbp_fetal_tss), names(RBPs_list))
RBPs_list_sig <- RBPs_list[matching_RBP_names]

matching_TF_names <- match(tf_fetal_tss_sig$TF, names(TFregulon_list))
TFregulon_list_sig <- TFregulon_list[matching_TF_names]

fetal_tss <- read.table("Figure5/human_fetal_TSS-PAS.coordination.chiqtest.txt", sep = "\t", header = T, stringsAsFactors = FALSE)
fetal_tss <- fetal_tss[fetal_tss$FDR <= 0.05,]
fetal_tss_list <- fetal_tss$gene_id

result_list <- list()
for (i in 1:length(TFregulon_list_sig)) {
  for (j in 1:length(RBPs_list_sig)) {

    overlap_B_A <- intersect(fetal_tss_list, TFregulon_list_sig[[i]])  
    overlap_C_A <- intersect(fetal_tss_list, RBPs_list_sig[[j]])  
    overlap_B_C_A <- intersect(overlap_B_A, overlap_C_A)
    
    N <- length(unique(c(overlap_B_A, overlap_C_A)))  
    M <- length(overlap_B_A)  # Number of genes in overlap_B_A
    n <- length(overlap_C_A)  # Number of genes in overlap_C_A
    X_B_C_A <- length(overlap_B_C_A)  # Number of overlapping genes between B-A and C-A
    
    # Hypergeometric test for the overlap between overlap_B_A and overlap_C_A
    p_value_B_C_A <- phyper(X_B_C_A - 1, M, 1211 - M, n, lower.tail = FALSE)
    
    # Store the result for this pair in the result list
  result_list[[paste(names(TFregulon_list_sig)[i], names(RBPs_list_sig)[j], sep = "-")]] <- list(
      TF_B = names(TFregulon_list_sig)[i],
      RBP_C = names(RBPs_list_sig)[j],
      P_value = p_value_B_C_A,
      N = N,
      M = M,
      n = n,
      overlap = X_B_C_A,
      Total_TF = length(TFregulon_list_sig[[i]]),
      Total_RBP = length(RBPs_list_sig[[j]]),
      gene = paste(overlap_B_C_A, collapse = ", ")
    )
  }
}

result_df <- do.call(rbind, lapply(result_list, function(x) {
  data.frame(TF = x$TF_B,
             RBP = x$RBP_C,
             p_value = x$P_value,
             total_coupling = 1211,
             TF_coupling = x$M,
             RBP_coupling = x$n,
             Total_TF <- x$Total_TF,
             Total_RBP <- x$Total_RBP,
             overlap = x$overlap,
             overlap_gene = x$gene,
             stringsAsFactors = FALSE)
})) 
result_df$fdr = p.adjust(as.numeric(result_df$p_value),'fdr')
result_df <- result_df[order(result_df$fdr, decreasing = FALSE), ] 
result_df <- result_df[result_df$fdr <= 0.05,]
result_df <- result_df %>%
  tidyr::separate(RBP, into = c("RBP", "source"), sep = "_")
result_df$logP <- -log10(result_df$fdr)


result_df$RBP <- reorder(result_df$RBP, table(result_df$RBP)[result_df$RBP])

ggplot(result_df, aes(x = TF, y =RBP , color = logP, size = overlap)) +
  geom_point() +                          # Dots for each TF-RBP pair
  scale_color_viridis_c(trans = "log10") + # Use a color scale based on log-transformed p_value
  scale_size_continuous(range = c(1, 3)) + # Adjust dot size based on overlap
  theme_minimal() +
  labs(title = "Hypergeometric test of overlapping genes between RBP targets and TF targets",
       x = "Transcriptional factors", 
       y = "RNA binding ptoteins",
       color = "Hypergeometric test\n P-value (log10)",
       size = "Overlap") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))  

