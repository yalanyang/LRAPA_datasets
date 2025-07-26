library(tidyverse)
library(readxl)
library(ggpubr)
library(readr)
library(UpSetR)


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

gtf_genes$gene_id <- sub("\\..*", "", gtf_genes$gene_id)

##coupling genes
fetal_tss <- read.table("Figure5/human_fetal_TSS-PAS.coordination.chiqtest.txt", sep = "\t", header = T, stringsAsFactors = FALSE)
fetal_tss$type <- "fetal_tss_sig"
fetal_tss$type[which(fetal_tss$FDR >= 0.05)] = "fetal_tss_nonsig"

fetal_tss_anno <- fetal_tss %>% left_join(SFARI, by = c("gene_id"="gene.symbol"))
fetal_tss <- fetal_tss %>% select("gene_id","type") %>% unique()

adult_tss <- read.table("Figure5/human_encode_TSS-PAS.coordination.chiqtest.txt", sep = "\t", header = T, stringsAsFactors = FALSE)
adult_tss$type <- "adult_tss_sig"
adult_tss$type[which(adult_tss$FDR >= 0.05)] = "adult_tss_nonsig"
adult_tss_anno <- adult_tss %>% left_join(SFARI, by = c("gene_id"="gene.symbol"))
adult_tss <- adult_tss %>% select("gene_id","type") %>% unique()


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


##Figure 5B
enrichment <- read.table("Figure5/Figure5B.enrichment.txt", sep = "\t", header = T)
library(ggplot2)
enrichment$Type <- factor(enrichment$Type, levels = c("DAG", "DIU", "DET", "DEG"))
ggplot(enrichment, aes(x = odds_ratio, y = Type, color = group)) +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  geom_errorbarh(aes(xmin = conf_lower, xmax = conf_upper), height = 0.2, position = position_dodge(width = 0.5)) +
  geom_text(aes(label = paste0("P = ", format(p_value, scientific = TRUE))), 
            hjust = -0.2, size = 2) +  
  theme_minimal() +
  labs(x = "Odds Ratio",  color = "Coupling type") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black") +  
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) + 
  theme(axis.text.y = element_text(size = 10), 
        axis.text.x = element_text(size = 10),
        plot.title = element_text(size = 12, face = "bold"))


##Figure 5C
fetal_TSS_sig <- fetal_TSS %>% dplyr::filter(FDR < 0.05)
adult_TSS_sig <- adult_TSS %>% dplyr::filter(FDR < 0.05)
fetal_SE <- read.table("Figure5/human_fetal_exon-PAS.full-length.coordination.chiqtest.txt", sep = "\t", header = T, stringsAsFactors = FALSE)
fetal_SE_sig  <- fetal_SE %>% dplyr::filter(FDR < 0.05)
adult_SE  <- read.table("Figure5/human_encode_exon-PAS.coordination.chiqtest0203.txt", sep = "\t", header = T, stringsAsFactors = FALSE)
adult_SE_sig  <- adult_SE %>% dplyr::filter(FDR < 0.05)
listInput <- list(fetal_TSS=unique(fetal_TSS_sig$gene_id),adult_TSS=unique(adult_TSS_sig$gene_id),
                  fetal_SE=unique(fetal_SE_sig$gene_id),adult_SE=unique(adult_SE_sig$gene_id))
upset(fromList(listInput), order.by = "freq")


##Figure5D
splicesome <- read.table('reference/genelist/Splicesome',header=F,sep='\t')
splicesome <- as.character(splicesome$V1)
constrain <- read.table('reference/genelist/Constrained_gene.txt',header=T,sep='\t')
constrain <- as.character(constrain$gene)
genes <- read.table('reference/genelist/gene.list.txt',header=TRUE,sep='\t')
PSD<- genes[which(genes$set_PSD == "1") ,]
PSD<- as.character(PSD$gene)
Embryonic<- genes[which(genes$set_Embryonic == "1") ,]
Embryonic<- as.character(Embryonic$gene)
ChromatinModifiers <- genes[which(genes$set_ChromatinModifiers == "1") ,]
ChromatinModifiers<- as.character(ChromatinModifiers$gene)
EssentialGenes <- genes[which(genes$set_EssentialGenes == "1") ,]
EssentialGenes<- as.character(EssentialGenes$gene)
target <- list(PSD=PSD, Embryonic=Embryonic, Constrainedgenes=constrain,Splicesome=splicesome,
               Chromatinmodifiers=ChromatinModifiers,EssentialGenes=EssentialGenes)

target_adult_tss <- lapply(names(target), function(d) {
  perform_fisher_test(target[[d]], adult_tss[adult_tss$type=="adult_tss_sig",]$gene_id, adult_tss[adult_tss$type=="adult_tss_nonsig",]$gene_id)
})
target_adult_tss <- as.data.frame(do.call(rbind, target_adult_tss), row.names = names(target))

target_fetal_tss <- lapply(names(target), function(d) {
  perform_fisher_test(target[[d]], fetal_tss[fetal_tss$type=="fetal_tss_sig",]$gene_id, fetal_tss[fetal_tss$type=="fetal_tss_nonsig",]$gene_id)
})
target_fetal_tss <- as.data.frame(do.call(rbind, target_fetal_tss), row.names = names(target))

target_fetal_tss$type <- "fetal_tss"
target_fetal_tss$disease <- rownames(target_fetal_tss)
target_adult_tss$type <- "adult_tss"
target_adult_tss$disease <- rownames(target_adult_tss)
target_tss_combined <- rbind(target_fetal_tss, target_adult_tss)
target_tss_combined$fdr = p.adjust(as.numeric(target_tss_combined$p_value),'fdr')
target_tss_combined$p_label <- sprintf("p=%.2g", target_tss_combined$fdr)  # Format to 2 significant digits
ggplot(target_tss_combined, aes(x = `odds_ratio.odds ratio`, y = disease, color = type)) +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +  # Dots for OR
  geom_errorbarh(aes(xmin = conf_lower, xmax = conf_upper), height = 0.2, position = position_dodge(width = 0.5)) +  # 95% CI
  geom_text(aes(label = p_label), hjust = -0.2, size = 3.5, position = position_dodge(width = 0.5)) +  # Show p-values next to points
  theme_minimal() +
  labs(title = "Odds Ratios for Disease Gene Enrichment in Adult vs. Fetal TSS",
       x = "Odds Ratio", y = "Disease", color = "TSS Type") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black") +  # Reference line at OR = 1
  theme(axis.text.y = element_text(size = 10), 
        axis.text.x = element_text(size = 10),
        plot.title = element_text(size = 12, face = "bold"))



##Figure 5E
##SFARI, https://www.sfari.org/resource/sfari-gene/
SFARI <- read.csv("reference/genelist/SFARI-Gene_genes_01-13-2025release_02-14-2025export.csv")
SFARI <- filter(SFARI, gene.score == 1 | syndromic == 1)
SFARI$disease <- "ASD"
SFARI <- SFARI %>% select("gene.symbol","disease")

##reference: https://www.nature.com/articles/s41588-022-01148-2
ASD_2022NG<- read.csv("reference/genelist/Autism_185_2022NG.txt",header=F)
ASD_2022NG$disease <- "ASD"
colnames(ASD_2022NG) <- c("gene.symbol","disease")
ASD <- unique(rbind(SFARI,ASD_2022NG))

##DDD(developmental disorder, 285 genes), https://www.nature.com/articles/s41586-020-2832-5
DDD <- read_excel("reference/genelist/41586_2020_2832_MOESM4_ESM.xlsx")
DDD <- filter(DDD, significant == TRUE)
DDD$disease <- "DDD"
DDD <- DDD %>% select("symbol","disease")
colnames(DDD) <- c("gene.symbol","disease")

##Schizophrenia exome meta-analysis consortium: https://schema.broadinstitute.org/results
SCZ_SCHEMA <- read.csv("reference/genelist/SCHEMA_meta_results_2025_02_14_16_15_53.csv",header=T)
SCZ_SCHEMA <- filter(SCZ_SCHEMA, Q.meta < 0.1)
SCZ_SCHEMA$disease <- "SCZ"
SCZ_SCHEMA <- SCZ_SCHEMA %>% select("Gene","disease")
nrow(SCZ_SCHEMA)  # 34 genes
SCZ_SCHEMA <-  as.data.frame(SCZ_SCHEMA) %>% left_join(gtf_genes, by = c("Gene"="gene_id"))
SCZ_SCHEMA <- SCZ_SCHEMA %>% select("gene_name","disease")
colnames(SCZ_SCHEMA) <- c("gene.symbol","disease")

##http://resource.psychencode.org
SCZ_PsychENCODE <- read.csv("reference/genelist/INT-17_SCZ_High_Confidence_Gene_List.csv",header=T)
SCZ_PsychENCODE$disease <- "SCZ"
SCZ_PsychENCODE <- SCZ_PsychENCODE %>% select("sczgenenames","disease")
colnames(SCZ_PsychENCODE) <- c("gene.symbol","disease")
SCZ <- unique(rbind(SCZ_SCHEMA,SCZ_PsychENCODE))

##https://epi25.broadinstitute.org/downloads
epilepsy = read_tsv('reference/genelist/Epi25_gene_results.tsv') %>% filter(group=="EPI" & pval < 0.05)
epilepsy <-  as.data.frame(epilepsy) %>% left_join(gtf_genes, by = "gene_id")
epilepsy$disease <- "epilepsy"
epilepsy <- epilepsy %>% select("gene_name","disease")
colnames(epilepsy) <- c("gene.symbol","disease")

##https://epi25.broadinstitute.org/downloads
BIP <-  read.table('reference/genelist/BipEx_gene_results.tsv',header=T,sep='\t')
BIP <- filter(BIP, damaging_missense_fisher_gnom_non_psych_pval < 0.01)
BIP$disease <- "BIP"
BIP <-  as.data.frame(BIP) %>% left_join(gtf_genes, by = "gene_id")
BIP <- BIP %>% select("gene_name","disease") %>% unique()
colnames(BIP) <- c("gene.symbol","disease")

##AD (Alzheimer’s disease)_associated_genes.txt: AD_associated_genes were downloaded from ADSP: https://adsp.niagads.org/index.php/gvc-top-hits-list/
AD <-  read.table('reference/genelist/AD_associated_genes.txt',header=F,sep='\t')
AD$disease <- "AD"
colnames(AD) <- c("gene.symbol","disease")

disease <- list(BIP=BIP$gene.symbol,ASD=ASD$gene.symbol,epilepsy=epilepsy$gene.symbol, SCZ=SCZ$gene.symbol,
                DDD=DDD$gene.symbol,AD=AD$gene.symbol)



#Enrichment in different diseases
disease_adult_tss <- lapply(names(disease), function(d) {
 perform_fisher_test(disease[[d]], adult_tss[adult_tss$type=="adult_tss_sig",]$gene_id, adult_tss[adult_tss$type=="adult_tss_nonsig",]$gene_id)
})
disease_adult_tss <- as.data.frame(do.call(rbind, disease_adult_tss), row.names = names(disease))

disease_fetal_tss <- lapply(names(disease), function(d) {
 perform_fisher_test(disease[[d]], fetal_tss[fetal_tss$type=="fetal_tss_sig",]$gene_id, fetal_tss[fetal_tss$type=="fetal_tss_nonsig",]$gene_id)
})
disease_fetal_tss <- as.data.frame(do.call(rbind, disease_fetal_tss), row.names = names(disease))


disease_adult_tss$type <- "adult_tss"
disease_adult_tss$disease <- rownames(disease_adult_tss)
disease_fetal_tss$type <- "fetal_tss"
disease_fetal_tss$disease <- rownames(disease_fetal_tss)
disease_tss_combined <- rbind(disease_adult_tss, disease_fetal_tss)
disease_tss_combined$fdr = p.adjust(as.numeric(disease_tss_combined$p_value),'fdr')
disease_tss_combined$p_label <- sprintf("p=%.2g", disease_tss_combined$fdr)  # Format to 2 significant digits
##Figure 5E plot
ggplot(disease_tss_combined, aes(x = `odds_ratio.odds ratio`, y = disease, color = type)) +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +  # Dots for OR
  geom_errorbarh(aes(xmin = conf_lower, xmax = conf_upper), height = 0.2, position = position_dodge(width = 0.5)) +  # 95% CI
  geom_text(aes(label = p_label), hjust = -0.2, size = 3.5, position = position_dodge(width = 0.5)) +  # Show p-values next to points
  theme_minimal() +
  labs(title = "Odds Ratios for Disease Gene Enrichment in Adult vs. Fetal TSS",
       x = "Odds Ratio", y = "Disease", color = "TSS Type") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black") +  # Reference line at OR = 1
  theme(axis.text.y = element_text(size = 10), 
        axis.text.x = element_text(size = 10),
        plot.title = element_text(size = 12, face = "bold"))
