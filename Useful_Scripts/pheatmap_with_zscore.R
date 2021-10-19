

##> Function of Z-score ,where each value of gene-counts is given as input
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}


signif.co = 0.05
l2fc.co = 1

#> Choose only signif from DESeq-results (Complete Table)
D.sig <- subset(DESeq_table,(padj <= signif.co &  abs(log2FoldChange) >= l2fc.co) )

#> in counts:
##> Use table "gene_counts" : normalised_counts per sample( in columns ) on all genes (in rows)
D_counts.sig <- D.sig[,colnames(gene_counts)]

##> Z-transform Counts to better illustrate in Heatmap :
D_counts.sig.z <- t(apply(D_counts.sig, 1, cal_z_score)) ## apply with "1" : applies function per row

#> remove NAs if any
##DEcounts.sig_norm <- DEcounts.sig_norm[-which(is.na(DEcounts.sig_norm[,1])),]


##> Apply Clustering methods for Genes:
genes_hclust <- hclust(dist(D_counts.sig.z), method = "complete")

##> Cut in groups as many as your Conditions:
clust_assign <- cutree(tree = genes_hclust, k = 3)

clust_name   <- data.frame(cluster = ifelse(clust_assign == 1,yes="cluster 1",
                                            ifelse(clust_assign == 2,yes="cluster 2",
                                                   no = "cluster 3")))

conditions.dds <- colnames(rawdata)
sample_annot <- data.frame(conditions=conditions.dds )
rownames(sample_annot) <-  colnames(rawdata)

#>> sample_annot
#             conditions
# CTRL_1        CTRL
# CTRL_2        CTRL
# CTRL_3        CTRL
# DM1_1          DM1
# DM1_2          DM1
# DM1_3          DM1



pheatmap(D_counts.sig.z,
         color = colorRampPalette(rev(brewer.pal(n=10,name ="RdYlBu")))(100),
         annotation_col= sample_annot,
         clustering_distance_rows="euclidean",
         cluster_cols=FALSE,
         annotation_legend = TRUE,
         annotation_names_col = TRUE,
         show_rownames=FALSE,
         labels_row = NULL,
         angle_col = 45,
         main =paste0(nrow(D_counts.sig.z)," Signif.Diff/ly expressed Genes in DM1-vs-CTRL"))
