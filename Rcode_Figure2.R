library(Biostrings)
library(optparse)
library(rtracklayer)
library(ggplot2)
library(reshape2)
library(data.table)
library(dplyr)
library(bedtoolsr)
library(tidyr)
library(GenomicRanges)
library(ggpubr)
library(ggsci)

reference_genome <- "reference/GRCh38.primary_assembly.genome.fa"
reference <- readDNAStringSet(reference_genome)
# Adjust chromosome names in reference if necessary
names(reference) <- sub(" +.*", "", names(reference))
Encode <- fread("Figure2/Encode.brain.new.PAS.bed", sep = "\t", header = T, stringsAsFactors = FALSE)
Encode <- Encode %>% filter(Chr %in% paste0("chr", c(1:22, "X", "Y")))
Encode <- Encode %>% dplyr::filter(gene_count>=10)
Encode <- Encode %>% dplyr::filter(PAU > 0.05)


###Figure 2B
feature <- data.frame(table(Encode$Feature))
feature$ratio <- feature$Freq / sum(feature$Freq)
sum_counts <- aggregate(Count ~ Feature, data = Encode, FUN = sum)
sum_counts$ratio <- sum_counts$Count / sum(sum_counts$Count)
g <- ggplot(feature, aes(x = "", y = ratio, fill = Var1)) +
  geom_bar(stat = "identity")  +
  theme_minimal() +  scale_fill_npg("nrc", alpha = 0.7) +
  labs(
    x = "",
    y = "Percent of locations in genome") +
  theme_classic() 
h <- ggplot(sum_counts, aes(x = "", y = ratio, fill = Feature)) +
  geom_bar(stat = "identity")  +
  theme_minimal() +  scale_fill_npg("nrc", alpha = 0.7) +
  labs(
    x = "",
    y = "Percent of PAS reads in genome") +
  theme_classic()

arrange <- ggarrange(g,h, ncol = 2, nrow = 1)
ggsave("PAS distribution.pdf", arrange, width = 6, height = 3)


###Figure 2C
hexamer_stat <- data.frame(table(Encode$Hexamer))
ggpie(data = hexamer_stat, "Freq", label ="Var1",,fill="Var1") +  scale_fill_npg("nrc", alpha = 0.7) + theme(legend.position = "none") 



##Figure 2D
PAS_bed <- Encode %>% separate(PAS_ID_new,  c("chr", "end", "strand"), sep = '[:]')
PAS_bed$end <- as.numeric(PAS_bed$end)
PAS_bed$start <- PAS_bed$end-1

bed_df <- PAS_bed %>% select(chr, start, end, Feature, Hexamer,Strand, anno, Count)

# Function to extract sequences around the end site considering strand information
extract_sequences <- function(row, genome) {
  chrom <- as.character(row['chr'])
  strand <- as.character(row['Strand'])
  
  # Convert row coordinates to numeric
  end <- as.numeric(row['end'])
  
  # Calculate start and end positions
  start <- if (strand == "-") {
    end - 100
  } else {
    end - 50
  }
  end <- if (strand == "-") {
    end + 50
  } else {
    end + 100
  }
  
  # Ensure chromosome is present
  if (!chrom %in% names(genome)) {
    stop(paste("Chromosome", chrom, "not found in genome"))
  }

  # Extract and potentially reverse complement the sequence
  seq <- subseq(genome[[chrom]], start = start, end = end)
  if (strand == "-") {
    seq <- reverseComplement(seq)
  }
  
  return(as.character(seq))
}

# Extract sequences around end sites considering strand information
window_size <- 100

# Function to calculate nucleotide frequencies at each position
calculate_frequencies <- function(sequences,window_size = 100) {
  pos_range <- -50:window_size
  freq_matrix <- matrix(0, nrow = length(pos_range), ncol = 4)
  colnames(freq_matrix) <- c("A", "C", "G", "T")
  for (i in seq_along(pos_range)) {
    pos_nucleotides <- substr(sequences, i, i)
    freq_matrix[i,] <- table(factor(pos_nucleotides, levels = c("A", "C", "G", "T")))
  }
  
  freq_df <- as.data.frame(freq_matrix)
  freq_df$Position <- pos_range
  return(freq_df)
}

sequences <- apply(bed_df, 1, extract_sequences, genome=reference)
  
freq_df <- calculate_frequencies(sequences, window_size = window_size)
  # Melt the data frame for plotting
 
  freq_df_normalized <- freq_df
  freq_df_normalized[1:4] <- t(apply(freq_df[1:4], 1, function(x) x / sum(x))) 
  freq_melted <- reshape2::melt(freq_df_normalized, id.vars = "Position", variable.name = "Nucleotide", value.name = "Ratio")
  
ggplot(freq_melted, aes(x = Position, y = Ratio, color = Nucleotide)) +
    geom_line(size = 1) +
    labs(
      title = "Nucleotide Ratios Around Cleavage Sites (-50 bp to +100 bp)",
      x = "Relative Position to Poly(A) Cleavage Site (nt)",
      y = "Nucleotide Ratio"
    ) +
    scale_color_npg("nrc", alpha = 0.7) +
    theme_minimal(base_size = 15) +
    theme(
      legend.position = "top",
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12)
    )

##Figure 2E
##Intergrative public database
chromosomes_to_keep <- c(paste0("chr", 1:22), "chrX", "chrY", "chrM")
PolyASite <- fread("reference/polyA_ref/atlas.clusters.2.0.GRCh38.96.bed",header=F,sep='\t')
PolyASite <- PolyASite[V1 %in% chromosomes_to_keep]
PolyASite <- unique(PolyASite[,1:6])
colnames(PolyASite) <- c('chr', 'start', 'end', 'name', 'score', 'strand')

PolyA_DB <- fread("reference/polyA_ref/human.polyDB.PAS.hg19to38.bed",header=F,sep='\t')
PolyA_DB <- PolyA_DB[V1 %in% chromosomes_to_keep]
PolyA_DB <- PolyA_DB %>% separate(V4,  c("name", "score", "strand"), sep = '[_]')
PolyA_DB$name="."
PolyA_DB$score="."
colnames(PolyA_DB) <- c('chr', 'start', 'end', 'name', 'score', 'strand')

gr <- fread("reference/hg38.refGene.gtf", sep = "\t", header = FALSE, 
            col.names = c("Chromosome", "Source", "Feature", "Start", "End", "Score", "Strand", "Frame", "Attributes"))
refgene <- gr[!grepl("*random|*alt|*fix|chrUn*|chrM",gr$Chromosome),]
refgene <- refgene[grepl("transcript",refgene$Feature),]
refgene$PAS <- ifelse(refgene$Strand == "+", refgene$End, refgene$Start)
refgene$Start <- refgene$PAS-1
refgene <- refgene %>% select(Chromosome,Start,PAS,Score,Frame,Strand)
colnames(refgene) <- c('chr', 'start', 'end', 'name', 'score', 'strand')
refgene <- unique(refgene)


PolyA_DB$database <- "PolyA_DB"
PolyASite$database <- "PolyASite"
refgene$database <- "RefGene"
annotated_PAS <- rbind(PolyA_DB, PolyASite, refgene)
annotated_PAS <- annotated_PAS[order(annotated_PAS$chr, annotated_PAS$start, annotated_PAS$end), ]
##Distribution of PAs in different genomics regions and annnotated
Encode_3UTR <- bed_df[grepl("3UTR|stop_codon",bed_df$Feature),]
Encode_nonUTR <- bed_df[!grepl("3UTR|stop_codon",bed_df$Feature),]
UTR3_anno <- unique(bt.intersect(a=Encode_3UTR, b=annotated_PAS,wa=T, s=T))
UTR3_noanno <- unique(bt.intersect(a=Encode_3UTR, b=annotated_PAS,wa=T, s=T, v=T))
nonUTR3_anno <- unique(bt.intersect(a=Encode_nonUTR, b=annotated_PAS,wa=T, s=T))
nonUTR3_noanno <- unique(bt.intersect(a=Encode_nonUTR, b=annotated_PAS,wa=T, s=T, v=T))

overlap <- data.frame(
  UTR3 = c(nrow(UTR3_anno),nrow(UTR3_noanno)),
  nonUTR3 = c(nrow(nonUTR3_anno),nrow(nonUTR3_noanno)))
rownames(overlap) <- c("Annotated", "Unannotated")


overlap$Type <- rownames(overlap)

overlap_melted <- reshape2::melt(overlap, id.vars = "Type")
colnames(overlap_melted) <- c("Type", "Feature", "PAS_num")

ggplot(overlap_melted, aes(x = Feature, y = PAS_num, fill = Type)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(x = "Feature",
       y = "PAS_num") +  theme(axis.text.x = element_text(hjust = 1, size = 14),
                               axis.text.y = element_text(size = 14),
                               axis.title.x = element_text(size = 16),
                               axis.title.y = element_text(size = 16),
                               plot.title = element_text(size = 18, hjust = 0.5),
                               legend.text = element_text(size = 14),
                               legend.title = element_text(size = 16)) + scale_fill_npg("nrc", alpha = 0.7)

##Figure 2F
UTR3_anno$feature <- "UTR3_anno"
UTR3_noanno$feature <- "UTR3_noanno"
nonUTR3_anno$feature <- "nonUTR3_anno"
nonUTR3_noanno$feature <- "nonUTR3_noanno"
anno <- rbind(UTR3_anno,UTR3_noanno,nonUTR3_anno,nonUTR3_noanno)
anno$read <- log10(anno$V8)
ggplot(anno, aes(x = read, color = feature)) +
  stat_ecdf(geom = "step", size = 1) + scale_color_npg("nrc", alpha = 0.7)+
  labs(x = "log10 (read count)", y = "Cumulative density") +
  theme_minimal()


##Figure 2G-I
temp=list.files(path="subsample/", pattern="*.bed")
stat <- data.frame()
names <- c("file","peak","peak_polyDB")
for (i in 1:length(temp)){
  Isoseq1 <- fread(paste0("subsample/",temp[i]), sep = "\t", header = TRUE, stringsAsFactors = FALSE,
                   col.names = c("chr", "start", "end", "Hexamer", "PAS_ID_old", "strand", "Count", "UTR_id", "Feature",
                                 "gene_name", "gene_count", "PAU", "annotated_PAS",  "anno","PAS_ID_new"))
  Isoseq1 <- Isoseq1[chr %in% chromosomes_to_keep]
  Isoseq1 <- Isoseq1 %>% dplyr::filter(PAU > 0.05)
  
  for (cover in c(1,5,10,20,40)){
    Isoseq1 <- Isoseq1 %>% dplyr::filter(gene_count >= cover)
    polyDB_Isoseq <- nrow(unique(bt.intersect(a=Isoseq1, b=annotated_PAS,wa=T,wb=F)))
    data_sum <- c(temp[i],cover,nrow(Isoseq1),polyDB_Isoseq)
    stat <-rbind(stat,data_sum)
  }
}

colnames(stat) <- c("file","coverage", "peak","peak_polyDB")
stat$peak_polyDB <- as.numeric(stat$peak_polyDB) 
stat$peak <- as.numeric(stat$peak)      
stat$ratio <- stat$peak_polyDB/stat$peak*100


stat$ID <-c(rep("10%",5),rep("20%",5), rep("30%",5), rep("40%",5), rep("50%",5), rep("60%",5),rep("70%",5),rep("80%",5),
            rep("90%",5), rep("All",5))

#Figure 2G
ggline(stat, x = "coverage", y = "peak",title="peak_number_diff_coverage", color = "ID",ylim=c(0,60000),legend = "right") + scale_colour_npg("nrc", alpha = 0.7)
#Figure 2H
ggline(stat, x = "coverage", y = "ratio",title="peak_in_PolyDB",color = "ID",ylim=c(60,100),legend = "right") + scale_colour_npg("nrc", alpha = 0.7)



#Figure 2I
stat2 <- data.frame()
for (i in 1:length(temp)){
  Isoseq1 <- fread(paste0("subsample/",temp[i]), sep = "\t", header = TRUE, stringsAsFactors = FALSE,
                   col.names = c("chr", "start", "end", "Hexamer", "PAS_ID_old", "strand", "Count", "UTR_id", "Feature",
                                 "gene_name", "gene_count", "PAU", "annotated_PAS",  "anno","PAS_ID_new"))
  Isoseq1$Feature <- ifelse(Isoseq1$Feature == "stop_codon", "3UTR", Isoseq1$Feature)
  Isoseq1 <- Isoseq1[chr %in% chromosomes_to_keep]
  Isoseq1 <- Isoseq1 %>% dplyr::filter(PAU > 0.05 & gene_count >= 10)
  data_sum <- data.frame(table(Isoseq1$Feature))
  stat2 <-rbind(stat2,data_sum)
}

colnames(stat2) <- c("Feature","Num")
stat2$group <- c(rep("10%",5), rep("20%",5), rep("30%",5), rep("40%",5), rep("50%",5), rep("60%",5),rep("70%",5),
                 rep("80%",5), rep("90%",5),rep("All",5))
ggplot(stat2, aes(group, weight = Num, fill = Feature,title="Distribution of PASs")) +
  geom_bar(color = "black", width = .7, position = 'fill') +
  scale_fill_npg("nrc", alpha = 0.7)+coord_flip()+
  labs( y = 'ratio') + 
  scale_y_continuous(expand = c(0,0)) +
  theme_classic()
