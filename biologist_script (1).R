library(ggplot2)
library(tidyverse)

#import rnaseq results 
rnaseq_results <- read.delim("GSE75070_MCF7_shRUNX1_shNS_RNAseq_log2_foldchange.txt.gz", sep = "\t")
diff_genes <- subset(rnaseq_results, padj <= 0.01 & abs(log2FoldChange) > 1)

#up-regulated
up_reg <- diff_genes %>% 
  dplyr::filter(log2FoldChange > 1)
#down-regulated
down_reg <- diff_genes %>% 
  dplyr::filter(log2FoldChange < -1)

#combining up and down-reg
data_filtered <- rbind(up_reg, down_reg)

#annotated peak file and diff genes
anno_results <- read.delim("/projectnb/bf528/users/group_5/project_3/programmer/annotated_peaks.txt", sep = "\t")
colnames(anno_results)[16] <- "genename"

#getting peaks
tss_peaks <- subset(anno_results, abs(`Distance.to.TSS`) <= 5000)

data_filtered <- dplyr::mutate(data_filtered,
                               bind_status = ifelse(data_filtered$genename %in% anno_results$genename, 'RUNX1_bound', 'Not_bound'),
                               reg_status = ifelse(data_filtered$genename %in% up_reg$genename, 'Up-regulated', 'Down-regulated'))


# creating a new data frame with imp columns
plot_df1 <- group_by(data_filtered, reg_status, bind_status) %>%
  summarize(count = n()) %>%
  mutate(percent = count / sum(count) * 100)

#factoring to get up before down in the plot
plot_df1$reg_status <- factor(plot_df1$reg_status, levels = c('Up-regulated', 'Down-regulated'))

#Create the barplot using ggplot2
ggplot(plot_df1, aes(x = reg_status, y = percent, fill = bind_status)) +
  geom_bar(stat = "identity",position = "stack") +
  scale_fill_manual(values = c("grey", "red")) + 
  labs(x = "+/- 5kb of TSS", y = "Percentage of Genes") +
  ggtitle("RUNX1 peak binding +/− 5kb of transcriptional start site (TSS)") +
  theme_minimal() +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  geom_text(aes(label=count, y=percent) ,position = position_stack(vjust = 0.5)) +
  theme(plot.title = element_text(size = 10))

#######Heatmap part########

#Load the HI-C data
hic_data <- read.delim("GSE75070_HiCStein-MCF7-shGFP_hg19_chr10_C-40000-iced.matrix.gz", sep = "\t")
subset_hic_data <- hic_data[1:134, 1:134] #first 5.3 mega base
subset_hic_data[is.na(subset_hic_data)] <- 0 #convert NaN values to 0
subset_without_col_1 <- subset_hic_data[, -1]
subset_log <- log2(subset_without_col_1 + 1) #log normalization
complete_df <- cbind(subset_hic_data[1], subset_log)
df_numeric <- data.frame(lapply(complete_df[,-1], function(x) as.numeric(x, na.rm = TRUE)))
heatmap(as.matrix(complete_df[,-1]), Rowv = NA, Colv = NA, labCol= FALSE, labRow = FALSE) #heatmap without clustering










































