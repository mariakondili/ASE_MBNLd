#!/usr/bin/env R

###~~~ Run DESeq2 for genes by HTSeq-counts, after STAR.Alignment ~~~~###
### Author: Maria Kondili

library(biomaRt)
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(devtools))
suppressPackageStartupMessages(library(dendextend))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(reshape2))


data_dir <- "/projects/mbnl_dct/mouse_seq/"

x <- read.table(paste0(data_dir, "metadata_MBNLdelta_samples.txt"),
        header=T,stringsAsFactors=F, as.is=TRUE)

##> metadata_MBNLdelta_samples:
# Sample	                 mouse_id     replicate	  treatment	  group
# 1051_htseq_counts.txt	           1051	        r1	    Saline        CTRL
# ...                                           r1          AAVGFP        AAVGFP
# ...                                           r1          MBNLd         AAVGFP-decoy


#### --> Read each htseq-file and Merge them all to a Data.Frame  ####
####    (genes in lines // samples in columns)

counts_dir <- paste0(data_dir,"Cutadapt/STAR_Aligned/HTSeq_counts/")

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


genes_IDs <- rownames(counts_data) <- counts_data$geneID

write.table(counts_data, file=paste0(data_dir,"DESeq2_analysis/Merged_CountData_AAVGFP_MBNLd_trimmed.txt"),
            row.names=F, sep="\t")

counts_data <- counts_data[,-1] # keep rownames but ignore 1st col with geneID
rownames(counts_data) <- genes_IDs


###-->  Filter Low-count Reads
nbreads <- rowSums(counts_data)
counts_data.filt <- counts_data[-which(nbreads<30),] # 20140 genes left


mart <- useMart("ensembl", dataset="mmusculus_gene_ensembl");
## doc = https://rdrr.io/bioc/biomaRt/f/vignettes/biomaRt.Rmd



condition  <- gsub("_.*", "",colnames(counts_data.filt))

###---> Create Design for Diff.expression Analysis
dctDesign <- data.frame( row.names=colnames(counts_data.filt), condition )
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


dds <- DESeqDataSetFromMatrix(countData = counts_data.filt, colData=dctDesign, design=~condition)
des <- DESeq(dds)


res.treat1 <- results(des, contrast=c("condition","AAVGFP" ,"Saline"))
#res.treat1@elementMetadata$description
colnames(res.treat1)[-1] <- paste0(colnames(res.treat1)[-1],"_AAVGFPvsCTRL")

res.treat2 <- results(des, contrast=c("condition","MBNLd","AAVGFP"))
colnames(res.treat2)[-1] <- paste0(colnames(res.treat2)[-1],"_MBNLdvsAAVGFP")

res.treat3 <- results(des, contrast=c("condition","MBNLd", "Saline"))
colnames(res.treat3)[-1] <- paste0(colnames(res.treat3)[-1],"_MBNLdVsCTRL")

# rownames(res.treat1) == rownames(res.treat2) == rownames(res.treat3)

#! baseMean is same for both res
#! Results contain geneIDs in rownames !


###---->  Merge Diff.analysis of 3 conditions

#y <- cbind(as.data.frame(res.treat1),as.data.frame(res.treat2)[,-1])
y <- cbind(as.data.frame(res.treat1),
            as.data.frame(res.treat2)[,-1],
            as.data.frame(res.treat3)[,-1])


# y=merge(res,res2)
nums <- which(sapply(y, is.numeric))
for(i in nums) y[,i] <- round(y[,i],3)


y$GeneID <- rownames(y)
DE_tab <- cbind(y, counts_data.filt)

###--> Add more Annotation Information
df.anns <- getBM(attributes = c("ensembl_gene_id","external_gene_name","description","gene_biotype"),
                 filters    = "ensembl_gene_id",
                 values     = DE_tab$GeneID,
                 mart       = mart)

colnames(df.anns)[1:2] <- c("GeneID","geneName")
DE_tab.ann <- merge(DE_tab, df.anns,by="GeneID", all.x=T)

###--> Save new data.frame with annot.
out_dir <- paste0(data_dir,"DESeq2_analysis/Results")
dir.create(out_dir,showWarnings = FALSE)
write.table(DE_tab.ann, paste0(out_dir,"DCT_Diff_Expression_MBNLdecoy-AAVGFP-CTRL.tsv"),
            quote=F,sep="\t",row.names=F,col.names=T)


oCols <- c("GeneID","geneName",
           "description","gene_biotype","baseMean",
           "log2FoldChange_AAVGFPvsCTRL","log2FoldChange_MBNLdVsCTRL",
           "pvalue_AAVGFPvsCTRL","padj_AAVGFPvsCTRL",
           "pvalue_MBNLdVsCTRL", "padj_MBNLdVsCTRL",
            colnames(counts_data.filt))


if ( all(oCols %in% colnames(DE_tab.ann)) ) {
  # D <- DE_tab.ann.filt[order(abs(DE_tab.ann.filt$log2FoldChange_MBNLdVsCTRL),decreasing=T), oCols]
  # ! attention: for large log2FC values ,the p-values and p.adj are NA !
        
  D <- DE_tab.ann[order(DE_tab.ann$padj_MBNLdVsCTRL,decreasing=F), oCols]
  sigD <- subset(D, padj_MBNLdVsCTRL < 0.05 & abs(log2FoldChange_MBNLdVsCTRL)>=1 )

  write.table(sigD,
              paste0(out_dir,"DataExpression_signif_log2fc=1_MBNLd_vs_Ctrl.txt"),sep="\t",
              row.names=F,quote=F,na="NA")
}


####---VOLCANO_plot to show Expression of MBNL1-dependent genes, among all Upregulated----####

#D <- D[- which(is.na(D$padj_MBNLdVsCTRL)),]

sigD[,"mbnl_dependent"] <- FALSE
sigD$mbnl_dependent[which(sigD$geneName %in% unique(hits_clip$gene_symbol) )] <- TRUE
cols <- c("red", "grey")
names(cols) <- c(TRUE, FALSE)

# Remove Genes that have padj==0,bcs log10(0)= Infinity --> shown on top of y-axis always !

sigD.pmod <- sigD %>% mutate(padj_modif = case_when(padj_MBNLdVsCTRL == 0 & mbnl_dependent == FALSE ~ 0.000001,
                                                    padj_MBNLdVsCTRL == 0 & mbnl_dependent ==TRUE ~0.00001,
                                                    padj_MBNLdVsCTRL !=0  ~ padj_MBNLdVsCTRL) )

##log10(0.000001) = 6


ggplot(sigD.pmod, aes(x=log2FoldChange_MBNLdVsCTRL,
                      y= -log10(padj_modif),
                      color = mbnl_dependent)) +
  scale_colour_manual(values=cols) +
  geom_point(size = 2, alpha = 1, na.rm = T) +
  ggtitle(label = "Volcano Plot of MBNLdecoy-signif.expressed genes")  +
  xlab("log2 fold change") +
  ylab("-log10 adjusted p-value") +
  theme_bw(base_size = 14) +
  theme(legend.position = "right") +
  scale_y_continuous(trans = "log1p")


### LABEL More signif +highly expr genes
  ## Create a column to indicate which genes to label
  # res_tableOE_ordered$genelabels <- ""
  # res_tableOE_ordered$genelabels[1:10] <- rownames(res_tableOE_ordered)[1:10]
  # geom_text_repel(aes(x = log2FoldChange, y = -log10(padj), label = ifelse(genelabels == T, rownames(res_tableOE_ordered),"")))



####------- PCA -----------####
plots_dir <- paste0(out_dir, "Plots/")

# png(paste0(plots_dir,"PCA_of_samples.png"))
rld <- rlog(des,blind=F)
DESeq2::plotPCA(rld, intgroup="conditions")
# dev.off()


####-------------------- HEATMAP ------------------####

# I want FC > 1.5
# log(2* 1.5) > x
# 2^x > 1.5
# log2(2^x) > log2(1.5)
# x > 0.58

## Keep the highly expr. and Signif. genes to move on..
sigCounts <- counts_data.filt[which(rownames(counts_data.filt) %in% sigD$GeneID),]


cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}


sigCounts.z <- t(apply(sigCounts, 1, cal_z_score))

###> PRETTY-HEATMAP
genes_hclust <- hclust(dist(sigCounts.z), method = "complete")

## Split in 3 Clusters, for 3 conditions
clust_assign <- cutree(tree = genes_hclust, k = 3)
clust_name   <- data.frame(cluster = ifelse(clust_assign == 1,yes="cluster 1",
                                            ifelse(clust_assign == 2,yes="cluster 2",
                                                   no = "cluster 3")))


sample_annot <- dctDesign
plot_name <- paste0(plots_dir,"Expression_PHeatmap_log2fc1_padj0.05_zscoredCounts_",nrow(sigCounts.z) ,"genes.png")
png(plot_name)

  main_title = "Signif.Diff/ly expressed Genes in MBNLd-vs-CTRL"
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
           main = paste0(nrow(sigCounts.z),main_title ))

dev.off()

�������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������
