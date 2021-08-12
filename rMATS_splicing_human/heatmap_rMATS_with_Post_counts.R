#!/usr/bin/env R


###-- Author  : Maria Kondili
###-- Date    : 11 March 2021

####
# Subject: After merging Post+Stats tables,merging all event-types, calculating Correction,
# it's time to create a Heatmap and filter Significant events, to see Deregulated events side-by-side
####


cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

plot_Pheatmap <- function(ase, signif.co, psi.co, title,conditions=c("CTRL", "DM1","MBNLd"), plot_name ){

  suppressPackageStartupMessages(library(DESeq2))
  suppressPackageStartupMessages(library(dendextend))
  suppressPackageStartupMessages(library(pheatmap))
  suppressPackageStartupMessages(library(RColorBrewer))

  ase.filt <- subset(ase, (FDR_Ctrl_vs_DM1 <= signif.co &  abs(dPSI_Ctrl_vs_DM1) >= psi.co) )
  # (FDR_CTRL448_vs_MBNLd <= signif.co & abs(dPSI_CTRL448_vs_MBNLd) >= psi.co)

  cat("\nafter filtering,significant events are : ", nrow(ase.filt), "\n")

  ###  z.gene =  (E(g) - Mean(g.replicates) ) / M.A.D (g.replicates), mad=median.absolute.deviation
  x <- as.matrix(ase.filt[,8:19])

  s1 = conditions[1]
  s2 = conditions[2]
  s3 = conditions[3]

  x_norm <- t(apply(x, 1, cal_z_score))

  x_norm <- x_norm[-which(is.na(x_norm[,1])),]

  colnames(x_norm) <- c( paste(s1,"_",c("r1","r2","r3","r4","r5","r6"),sep=""),
                         paste(s2,"_",c("r1","r2","r3"),sep=""),
                         paste(s3,"_",c("r1","r2","r3"),sep=""))

  my_hclust_gene <- hclust(dist(x_norm), method = "complete")

  clust_assign <- cutree(tree = my_hclust_gene, k = 3)
  clust_name   <- data.frame(cluster = ifelse(clust_assign == 1,yes="cluster 1",
                                       ifelse(clust_assign == 2,yes="cluster 2",
                                             no = "cluster 3")))


  sample_annot <- data.frame(sample = rep(conditions, c(6,3,3)) )
  row.names(sample_annot) <- colnames(x_norm)

  #png(plot_name)
  pheatmap(x_norm,
           color = colorRampPalette(rev(brewer.pal(n=10,name ="RdYlBu")))(100),
           annotation_col=sample_annot,
           clustering_distance_rows="euclidean",
           cluster_cols=FALSE,
           annotation_legend = TRUE,
           annotation_names_col = TRUE,
           show_rownames=FALSE,
           labels_row =NULL,
           angle_col = 45)

  # print(p)
  # dev.off()
}
