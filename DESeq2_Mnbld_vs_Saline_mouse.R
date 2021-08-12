#!/usr/bin/env R

### R used : R version 4.0.2 (2020-06-22)
###~~~ Run DESeq2 for genes by HTSeq-counts, after STAR.Alignment ~~~~###
##~~~~ Compare Diff.expression of genes AAV-MBNLdelta -vs-Saline & AAV-GFP-vs-Saline ~~~~
### Author: Maria Kondili


suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(org.Mm.eg.db))
library(AnnotationDbi)
suppressPackageStartupMessages(library(devtools))
suppressPackageStartupMessages(library(dendextend))
suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(RColorBrewer))


data_dir <- "/mbnl_dct/mouse_seq/Nextflow_pipeline/"
deseq_dir <- "/mbnl_dct/mouse_seq/DESeq2_analysis/"


x <- read.table(paste0(deseq_dir, "metadata_MBNLdelta_samples.txt"),
        header=T,stringsAsFactors=F, as.is=TRUE)

## Create a table of counts,merging all samples
## genes in lines // samples in columns

counts_dir <- paste0(data_dir,"HTSeq_counts/")


for(i in 1:nrow(x)) {
    fileName = x$Sample[i]
    sampleID = paste(x$treatment[i],x$replicate[i],sep="_")
    ## remove other '_' in gene.names: (grep "_" $f)
    ## for f in ./*.tsv ; do  sed -i 's/Y_RNA/YRNA/g' $f;done
    ## for f in ./*.tsv ; do  sed -i 's/_.A/-A/g' $f;done
    tmp <- read.table(paste0(counts_dir,fileName), comment.char = "_", stringsAsFactors=F,header=F)
    colnames(tmp) <- c("geneID",sampleID)

    if(i==1){
      counts_data <- tmp; print(dim(counts_data))
    }else{
      counts_data <- merge(counts_data, tmp , by="geneID",all=T)}
}


genes_IDs <- rownames(counts_data) <- gsub("\\..*" ,"", counts_data$geneID)

write.table(counts_data, file=paste0(deseq_dir,"Merged_rawCountData_from_STAR_Saline_AAVGFP_MBNLd.txt"),
            row.names=F, sep="\t")

counts_data <- counts_data[,-1] # keep rownames but ignore 1st col with geneID


###>  Filter Low-count Reads
nbreads <- rowSums(counts_data)
counts_data.filt <- counts_data[-which(nbreads<30),]

suppressPackageStartupMessages(library(DESeq2))

###> Design of Counts for D.E.genes analysis
conditions  <- gsub("_.*", "",colnames(counts_data.filt))

dctDesign <- data.frame( row.names=colnames(counts_data.filt), conditions )
# no colname ="condition" needed
#It should be TRUE : all(rownames(dctDesign) ==  colnames(counts_data)), for dds

#> dctDesign
#              condition
# Saline_r1    Saline
# Saline_r2    Saline
# Saline_r3    Saline
# AAVGFP_r1    AAVGFP
# AAVGFP_r2    AAVGFP
# AAVGFP_r3    AAVGFP
# MBNLd_r1      MBNLd
# MBNLd_r2      MBNLd
# MBNLd_r4      MBNLd


dds <- DESeqDataSetFromMatrix(countData = counts_data.filt, colData=dctDesign, design=~conditions)

dds <- estimateSizeFactors(dds)
norm_counts <- counts(dds,normalized=TRUE)

# Saline_r1   Saline_r2   Saline_r3   AAVGFP_r1   AAVGFP_r2    AAVGFP_r3     MBNLd_r1    MBNLd_r2     MBNLd_r4
# ENSMUSG00000000001  1660.61646  1637.18008  1728.08191  2139.66388  2965.33620  1843.656913  2091.259233  1947.53384  2142.635157
# ENSMUSG00000000028   164.36136   177.73970   165.91628   214.35261   281.50416   186.242554   158.967298   141.94502   147.226419
# ENSMUSG00000000031 64265.29011 56776.97080 94978.13847 61640.85821 36046.84634 88894.004280 52505.614839 51405.75128 51695.121895
# ENSMUSG00000000037    25.50435    22.71118    13.61364    22.20770    23.85628     7.218704     6.911622    12.02924     9.815095
# ENSMUSG00000000049    49.11949    40.48515    29.77985    25.10436    23.85628    31.762296    42.457104    57.74035    63.798115

des <- DESeq(dds)

res.treat1 <- results(des, contrast=c("conditions","AAVGFP" ,"Saline"))
#res.treat1@elementMetadata$description
colnames(res.treat1)[-1] <- paste0(colnames(res.treat1)[-1],"_AAVGFPvsCTRL")

res.treat2 <- results(des, contrast=c("conditions","MBNLd","AAVGFP"))
colnames(res.treat2)[-1] <- paste0(colnames(res.treat2)[-1],"_MBNLdvsAAVGFP")

res.treat3 <- results(des, contrast=c("conditions","MBNLd", "Saline"))
colnames(res.treat3)[-1] <- paste0(colnames(res.treat3)[-1],"_MBNLdVsCTRL")

# rownames(res.treat1) == rownames(res.treat2) == rownames(res.treat3)

#! baseMean is same for both res
#! Results contain geneIDs in rownames !


##
## Merge Diff.analysis of 3 conditions
##

#y <- cbind(as.data.frame(res.treat1),as.data.frame(res.treat2)[,-1])
y <- cbind(as.data.frame(res.treat1),
            as.data.frame(res.treat2)[,-1],
            as.data.frame(res.treat3)[,-1])


### Find out Highly Regulated per comparison:
subset(y, abs(log2FoldChange_AAVGFPvsCTRL)>= 1 & padj_AAVGFPvsCTRL  < 0.05) %>% dim
subset(y, abs(log2FoldChange_MBNLdVsCTRL) >= 1 & padj_MBNLdVsCTRL   < 0.05) %>% dim
subset(y, abs(log2FoldChange_MBNLdvsAAVGFP) >=1 & padj_MBNLdvsAAVGFP< 0.05) %>% dim


# y=merge(res,res2)
nums <- which(sapply(y, is.numeric))
for(i in nums) y[,i] <- round(y[,i],3)

y$GeneID <- rownames(y)
#DE_tab <- cbind(y, counts_data.filt)

#---> or: with Normalised counts :
DE_tab <- cbind(y, norm_counts)

###> Annotate further with gene-name and Biotype
convert_ID2Name <- function(id_list){

  library(org.Mm.eg.db) # installed via biocLite. Cannot call library from an R-variable/object
  require(AnnotationDbi)
  gene_symbols <- mapIds(x=org.Mm.eg.db, keys=id_list,column="SYMBOL", keytype="ENSEMBL", multiVals = "first")
  return(gene_symbols)
}

DE_tab$geneName <- convert_ID2Name(DE_tab$GeneID)

out_dir <- paste0(deseq_dir,"Results/")
dir.create(out_dir,showWarnings = FALSE)

de_rawCounts_filename ="Diff_Expression_and_rawCounts_MBNLdecoy-AAVGFP-CTRL.tsv"
de_normCounts_filename = "Diff_Expression_and_normCounts_MBNLdecoy-AAVGFP-CTRL.tsv"

write.table(DE_tab, paste0(out_dir,de_normCounts_filename),
            quote=F,sep="\t",row.names=F,col.names=T)


oCols <- c("GeneID","geneName",
           "baseMean",
           "log2FoldChange_MBNLdVsCTRL",
           "pvalue_MBNLdVsCTRL", "padj_MBNLdVsCTRL",
           "log2FoldChange_AAVGFPvsCTRL",
           "pvalue_AAVGFPvsCTRL","padj_AAVGFPvsCTRL",
            colnames(norm_counts))



if ( all(oCols %in% colnames(DE_tab)) ) {
  D <- DE_tab[order(DE_tab$padj_MBNLdVsCTRL,decreasing=F), oCols]
  saveRDS(D,"DE_data_padj_sorted.rds")
  sigD <- subset(D, padj_MBNLdVsCTRL < 0.05 & abs(log2FoldChange_MBNLdVsCTRL)>=1 )
  #> 863 genes

  write.table(sigD,
              paste0(out_dir,"Signif_DataExpression_log2fc=1_MBNLd_vs_Saline.txt"),sep="\t",
              row.names=F,quote=F,na="NA")

  ## or : Choose the diff.genes by AAVGFP-vs-CTRL higher contrast,so we can see the MBNLdct & CTRL in common:
  sigD <- subset(D, padj_AAVGFPvsCTRL < 0.05 & abs(log2FoldChange_AAVGFPvsCTRL )>=1 )

  write.table(sigD,
              paste0(out_dir,"Signif_DataExpression_log2fc=1_AAVGFP_vs_Saline.txt"),sep="\t",
              row.names=F,quote=F,na="NA")
}



##------- PCA -----------##
plots_dir <- paste0(deseq_dir, "Plots/")
dir.create(plots_dir,showWarnings = F)

library(pheatmap)
suppressPackageStartupMessages(library(gplots))

png(paste0(plots_dir,"PCA_of_samples.png"))
    rld <- rlog(des,blind=F)
    DESeq2::plotPCA(rld,intgroup="conditions")
dev.off()


####-------------------- HEATMAP ------------------####

##> Keep the lines of DESeq with Significant genes
# I want FC > 1.5
# log(2* 1.5) > x
# 2^x > 1.5
# log2(2^x) > log2(1.5)
# x > 0.58

## Keep the highly expr. and Signif. genes to move on..

# sigD <- subset(D, padj_AAVGFPvsCTRL < 0.05 & abs(log2FoldChange_AAVGFPvsCTRL )>=1 )

sigCounts.norm <- norm_counts[which(rownames(norm_counts) %in% sigD$GeneID),]

##> From raw-Counts :
#sigCounts <- counts_data.filt[which(rownames(counts_data.filt) %in% sigD$GeneID),]


create_heatmap <- function(expr_counts, design,plot_name) {

  calc_z_score <- function(x){
    (x - mean(x)) / sd(x)
  }

  Counts.z <- t(apply(expr_counts, 1, calc_z_score))

  ###> PRETTY-HEATMAP
  genes_hclust <- hclust(dist(Counts.z), method = "complete")

  clust_assign <- cutree(tree = genes_hclust, k = 3)
  clust_name   <- data.frame(cluster = ifelse(clust_assign == 1,yes="cluster 1",
                                              ifelse(clust_assign == 2,yes="cluster 2",
                                                     no = "cluster 3")))

  sample_annot <- design

  #png(plot_name)

  pheatmap(sigCounts.z,
           color = colorRampPalette(rev(brewer.pal(n=10,name ="RdYlBu")))(100),
           annotation_col=sample_annot,
           clustering_distance_rows="euclidean",
           cluster_cols=FALSE,
           annotation_legend = TRUE,
           annotation_names_col = TRUE,
           show_rownames=FALSE,
           labels_row = NULL,
           angle_col = 45,
           main = paste(nrow(sigCounts.z),plot_name,sep=" "))

  #dev.off()
}


###
### Call Heatmap-function for MBNLd-vs-Saline diff.expr. genes :
###

sigD <- subset(D, padj_MBNLdVsCTRL < 0.05 & abs(log2FoldChange_MBNLdVsCTRL)>=1 )

create_heatmap(sigCounts.norm, dctDesign,plot_name = "Signif.Diff/ly expressed Genes MBNLdct-vs-Saline")



####------ CROSS-CORRELATION OF SAMPLES -----------------####

sigCounts.rcorr = Hmisc::rcorr(as.matrix(sigCounts.norm),type="spearman")
#> nothing changes if raw counts_data.filt are used here !

#pdf("Plots/Cross-Correlation_Heatmap_10samples.pdf")
heatmap(x = sigCounts.rcorr$r,
        col = colorRampPalette(brewer.pal(8, "YlOrRd"))(25),
        symm = TRUE,cexCol=0.9,cexRow=0.9,cex=3,
        labCol=rownames(sigCounts.rcorr$r),
        margins = c(7,7),
        main="Spearman Cross-Correlation of Samples")
#dev.off()


###
### without AAVGFP:
###

sigCounts.rcorr = Hmisc::rcorr(as.matrix(sigCounts.norm[,c(1:3,7:9)]),type="spearman")

heatmap(x = sigCounts.rcorr$r,
        col = colorRampPalette(brewer.pal(8, "YlOrRd"))(25),
        symm = TRUE,cexCol=0.9,cexRow=0.9,cex=3,
        labCol=rownames(sigCounts.rcorr$r),
        margins = c(7,7),
        main="Spearman Cross-Correlation of Samples")
