library(tidyverse)
library(readxl)
library(ggpubr)

##Figure 6B
fetal_SE <- read.table("Figure5/human_fetal_exon-PAS.full-length.coordination.chiqtest.txt", sep = "\t", header = T, stringsAsFactors = FALSE)
fetal_SE$type <- "fetal_se_sig"
fetal_SE$type[which(fetal_SE$FDR >= 0.05)] = "fetal_se_nonsig"
adult_SE  <- read.table("Figure5/human_encode_exon-PAS.coordination.chiqtest0203.txt", sep = "\t", header = T, stringsAsFactors = FALSE)
adult_SE$type <- "adult_se_sig"
adult_SE$type[which(adult_SE$FDR >= 0.05)] = "adult_se_nonsig"

perform_fisher_test <- function(disease_genes, sig_genes, nonsig_genes) {
  in_sig <- sum(disease_genes %in% sig_genes)
  in_nonsig <- sum(disease_genes %in% nonsig_genes)
  not_in_sig <- length(sig_genes) - in_sig
  not_in_nonsig <- length(nonsig_genes) - in_nonsig

  contingency_table <- matrix(c(in_sig, in_nonsig, not_in_sig, not_in_nonsig),
                              nrow = 2, 
                              dimnames = list(c("In Disease", "Not in Disease"),
                                              c("adult_tss_sig", "adult_tss_nonsig")))
  fisher_result <- fisher.test(contingency_table)
  p_value <- fisher_result$p.value
  odds_ratio <- fisher_result$estimate
  conf_int <- fisher_result$conf.int

  return(c(in_sig=in_sig, in_nonsig=in_nonsig, not_in_sig=not_in_sig, not_in_nonsig=not_in_nonsig, p_value = p_value, odds_ratio = odds_ratio, conf_lower = conf_int[1], conf_upper = conf_int[2]))
}

rmats <- read.table("Figure6/SE.MATS.JCEC.txt", sep = "\t", header = T, stringsAsFactors = FALSE)
rmats <- rmats %>% filter(FDR <= 0.05 & abs(IncLevelDifference) > 0.05) %>% arrange(FDR)
rmats$exonStart_0base=rmats$exonStart_0base+1
rmats$exon <- paste0(rmats$chr,":",rmats$exonStart_0base,"-",rmats$exonEnd,":",rmats$strand)
adult_SE_sig <- adult_SE %>% filter(FDR <= 0.05) %>% arrange(FDR)
adult_SE_rmats <- inner_join(rmats, adult_SE_sig, by="exon")
SE_rmats_adult <- data.frame(perform_fisher_test(unique(rmats$exon), adult_SE[adult_SE$type=="adult_se_sig",]$exon, adult_SE[adult_SE$type=="adult_se_nonsig",]$exon))
SE_rmats_fetal <- data.frame(perform_fisher_test(unique(rmats$exon), fetal_SE[fetal_SE$type=="fetal_se_sig",]$exon, fetal_SE[fetal_SE$type=="fetal_se_nonsig",]$exon))
 
SE_rmats <- cbind(SE_rmats_adult,SE_rmats_fetal)
colnames(SE_rmats) <- c("adult","fetal")
SE_rmats <- data.frame(t(SE_rmats))
SE_rmats$stage <- rownames(SE_rmats)

SE_rmats$p_label <- sprintf("p = %.2g", SE_rmats$p_value) 
ggplot(SE_rmats, aes(x = odds_ratio.odds.ratio, y = stage)) +
  geom_point(size = 3, position = position_dodge(width = 0.5),color = "#E64B35FF") +  # Dots for OR
  geom_errorbarh(aes(xmin = conf_lower, xmax = conf_upper), height = 0.2, position = position_dodge(width = 0.5)) +  # 95% CI
  geom_text(aes(label = p_label), hjust = -0.2, size = 3.5, position = position_dodge(width = 0.5)) +  # Show p-values next to points
  theme_minimal() +
  labs(
       x = "Odds Ratio (95% confidence interval)") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black") +  # Reference line at OR = 1
  theme(axis.text.y = element_text(size = 10), 
        axis.text.x = element_text(size = 10),
        plot.title = element_text(size = 12, face = "bold"))



##Figure 6C
apalyzer_UTR <- read.table("Figure6/APAlyzer/3UTRAPA.diff.txt", sep = "\t", header = T, stringsAsFactors = FALSE)
apalyzer_IPA <- read.table("Figure6/APAlyzer/IPAmuti.diff.txt", sep = "\t", header = T, stringsAsFactors = FALSE)
apalyzer_UTR <- apalyzer_UTR %>% filter(pvalue <= 0.05) %>% arrange(pvalue)
apalyzer_IPA <- apalyzer_IPA %>% filter(pvalue <= 0.05) %>% arrange(pvalue)
apalyzer_diff <- unique(c(apalyzer_UTR$gene_symbol,apalyzer_IPA$gene_symbol))
SE_apalyzer_adult <- data.frame(perform_fisher_test(unique(apalyzer_diff), unique(adult_SE[adult_SE$type=="adult_se_sig",]$gene_id), unique(adult_SE[adult_SE$type=="adult_se_nonsig",]$gene_id)))
SE_apalyzer_fetal <- data.frame(perform_fisher_test(unique(apalyzer_diff), unique(fetal_SE[fetal_SE$type=="fetal_se_sig",]$gene_id), unique(fetal_SE[fetal_SE$type=="fetal_se_nonsig",]$gene_id)))

SE_apalyzer <- cbind(SE_apalyzer_adult,SE_apalyzer_fetal)
colnames(SE_apalyzer) <- c("adult","fetal")
SE_apalyzer <- data.frame(t(SE_apalyzer))
SE_apalyzer$stage <- rownames(SE_apalyzer)
SE_apalyzer$p_label <- sprintf("p = %.2g", SE_apalyzer$p_value) 
ggplot(SE_apalyzer, aes(x = odds_ratio.odds.ratio, y = stage)) +
  geom_point(size = 3, position = position_dodge(width = 0.5),color = "#E64B35FF") +  # Dots for OR
  geom_errorbarh(aes(xmin = conf_lower, xmax = conf_upper), height = 0.2, position = position_dodge(width = 0.5)) +  # 95% CI
  geom_text(aes(label = p_label), hjust = -0.2, size = 3.5, position = position_dodge(width = 0.5)) +  # Show p-values next to points
  theme_minimal() +
  labs(
    x = "Odds Ratio (95% confidence interval)") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black") +  # Reference line at OR = 1
  theme(axis.text.y = element_text(size = 10), 
        axis.text.x = element_text(size = 10),
        plot.title = element_text(size = 12, face = "bold"))


