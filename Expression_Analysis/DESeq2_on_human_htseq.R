#!/usr/bin/env R

###~~~ Run DESeq2 for HTSeq-counts, after STAR.Alignment ~~~~###
### Author: Naira Naouar
### Updated by : Maria Kondili


suppressPackageStartupMessages(library(DESeq2))

suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(devtools))
suppressPackageStartupMessages(library(ggplot2))

suppressPackageStartupMessages(library(dendextend))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(RColorBrewer))


work_dir <- "/projects/mbnl_dct/human_cellline/DESeq2_analysis/"
getwd(work_dir)
data_dir <- paste0( work_dir,"HtS_counts/")

samplesheet <- read.table(paste0(data_dir, "DESeq2.metadata.txt"),
        header=T,stringsAsFactors=F, as.is=TRUE)

###
### Create Counts-table,merging all samples in columns
###

for(i in 1:nrow(samplesheet)){

    sample=samplesheet$sampleName[i]
    fileName=samplesheet$fileName[i]
    ## remove other '_' in gene.names: (grep "_" $f)
    ## for f in ./*.tsv ; do  sed -i 's/Y_RNA/YRNA/g' $f;done
    ## for f in ./*.tsv ; do  sed -i 's/_.A/-A/g' $f;done
    tmp = read.table(fileName,comment.char = "_", stringsAsFactors=F,h=F)
    colnames(tmp)=c("geneID","gene_name",sample)

    if(i==1) {rawdata <- tmp[,c(1,3)]; print(dim(rawdata))}

    else rawdata <- merge(rawdata,tmp[,c(1,3)] , by="geneID",all=T)
}


genes_IDs <- rownames(rawdata) <- rawdata$geneID
rawdata <- rawdata[,-1] # keep rownames but ignore 1st col with geneID

###> Keep in an RDS object to easily retrieve:
# saveRDS(rawdata,"RawData_Expr_all_samples.rds")
# rawdata <- readRDS("RawData_Expr_all_samples.rds")


nb_reps = 6 # For CTRL ,we have 2 C25-clones x 3 replicates =6 CTRL replicates

###-------------#

colnames(rawdata) <- c( paste(rep("CTRL",6),c(1:6),sep="_"),
                        paste(rep( c("DM1","MBNLd"),c(3,3)),c("1","2","3"),sep="_" ))

#colnames(rawdata) <- paste( rep( c("CTRL","DM1","MBNLd") , c(3,3,3)), c("1","2","3"),sep="_" )

write.table(rawdata, file=paste0(work_dir,"merged_countData_dct.txt"), row.names=T, sep="\t")

#####   Filter Low-count Reads & Create DESEQ dataset #####
genes2keep <- which(rowSums(rawdata) >=30)
rawdata.filt <- rawdata[genes2keep,]

saveRDS(rawdata.filt, "rawdata_filt.rds")

#(i) rowSums(matrix >10) = gives you how many samples have count >30 for a gene.
#(i) rowSums(matrix >10 ) >3 : if found more than 3 samples with counts >10 == TRUE



###
###  Create DESeq Object :
###


conditions.dds <- gsub("_.*", "",colnames(rawdata.filt))
colData.dds <- data.frame(conditions=conditions.dds )
rownames(colData.dds) <- colnames(rawdata.filt)

#$ colData.dds
#             conditions
# CTRL_1        CTRL
# CTRL_2        CTRL
# CTRL_3        CTRL
# CTRL_4        CTRL
# CTRL_5        CTRL
# CTRL_6        CTRL
# DM1_1          DM1
# DM1_2          DM1
# DM1_3          DM1
# MBNLd_1      MBNLd
# MBNLd_2      MBNLd
# MBNLd_3      MBNLd

dds.mat <- DESeqDataSetFromMatrix(rawdata.filt, colData=colData.dds, design=~conditions)
# colData= Rows of colData correspond to columns of countData
# design formula: from colname  "~conditions" of colData

dds.mat <- estimateSizeFactors(dds.mat)
norm_dds_counts <- counts(dds.mat, normalized=TRUE)
## DESeq is always applied in DDS object. Norm_counts will be used later
des <- DESeq(dds.mat)

res.dm1 <- results(des,contrast=c("conditions","DM1" ,"CTRL"))
colnames(res.dm1)[-1] <- paste0(colnames(res.dm1)[-1],"_DM1vsCTRL")

res.dct <- results(des,contrast=c("conditions","MBNLd","CTRL"))
colnames(res.dct)[-1] <- paste0(colnames(res.dct)[-1],"_MBNLdvsCTRL")

## baseMean is same for both res
## Results contain geneIDs in rownames !

y <- cbind(as.data.frame(res.dm1),as.data.frame(res.dct)[,-1])


nums <- which(sapply(y, is.numeric))
for(i in nums) y[,i]=round(y[,i],3)

y$GeneID <- rownames(y)


###> Join Statistics of DESeq with Norm.counts table
DE_tab <- as.data.frame(cbind(y, norm_dds_counts))


## Change GeneID to match it to Ensembl-dB on GeneNames :
DE_tab$GeneID <- gsub("\\..*",  "", DE_tab$GeneID)
# DE_tab$geneIdx <- sub(".*\\.",  "",DE_tab$GeneID)

convert_ID2Name <- function(id_list){
    library(org.Hs.eg.db) # installed via biocLite. Cannot call library from an R-variable/object
    require(AnnotationDbi)
    gene_symbols <- mapIds(x=org.Hs.eg.db,
                          keys=id_list,column="SYMBOL",
                          keytype="ENSEMBL",
                          multiVals = "first")
    return(gene_symbols)
}

DE_tab$geneName <- convert_ID2Name(DE_tab$GeneID)


oCols <- c("GeneID","geneName",
           "log2FoldChange_DM1vsCTRL",
           "log2FoldChange_MBNLdvsCTRL",
            colnames(norm_dds_counts),"baseMean",
           "lfcSE_DM1vsCTRL","stat_DM1vsCTRL",
           "pvalue_DM1vsCTRL","padj_DM1vsCTRL",
           "lfcSE_MBNLdvsCTRL","stat_MBNLdvsCTRL",
           "pvalue_MBNLdvsCTRL","padj_MBNLdvsCTRL")

if ( all(oCols %in% colnames(DE_tab))) {
  DE_tab <- DE_tab[order(DE_tab$padj_DM1vsCTRL,decreasing=F),oCols]
}

##> save to a new file
write.table(DE_tab,paste0(pathData,"All_Counts_and_DESeq_Expression.txt"),sep="\t",row.names=F,quote=F,na="NA")


####---- CORRECTION Calculation -----#####

DE_tab$mean_Ctrl_expr  <- rowMeans(DE_tab[, c("CTRL_1","CTRL_2","CTRL_3", "CTRL_4", "CTRL_5","CTRL_6")])
# nb_reps ==3 :
# DE_tab$mean_Ctrl_expr  <- rowMeans(DE_tab[, c("CTRL_1","CTRL_2","CTRL_3")])
DE_tab$mean_DM1_expr   <- rowMeans(DE_tab[, c("DM1_1", "DM1_2","DM1_3")])
DE_tab$mean_MBNLd_expr <- rowMeans(DE_tab[, c("MBNLd_1","MBNLd_2","MBNLd_3")])

DE_tab$Correction_Pcnt <-  with(DE_tab,
                               ( (mean_MBNLd_expr - mean_DM1_expr ) / (mean_Ctrl_expr - mean_DM1_expr ) *100) )

# C[1] = DE_tab.ord$mean_MBNLd_expr[1] - DE_tab.ord$mean_DM1_expr[1] / (DE_tab.ord$mean_Ctrl_expr[1] - DE_tab.ord$mean_DM1_expr[1] )

# C  = (65.5 - 75.7) / (33.5 - 75.7) = -10.2 / -42.2 =  0.24

DE_tab$Correction_Pcnt  <- round(DE_tab$Correction_Pcnt,2)

write.table(DE_tab,paste0(pathData,"AllDataExpression_6CTRL+Correction.txt"),sep="\t",row.names=F,quote=F,na="NA")

###----------------------------- PLOTS -------------------------------------###

plots_dir <- paste0(work_dir,"Plots/")
dir.create(plots_dir,showWarnings = F)


#####-----PCA -----#####

rld <- rlog(dds.mat,blind=F)
pca_plot <- DESeq2::plotPCA(rld, intgroup="conditions")

png(paste0(plots_dir,"PCA_3conditions.png"))
  pca_plot + ggtitle("PCA on 3 conditions based on D.E.genes")
dev.off()


#####-------------------- HEATMAP ------------------######

# Necessary function for heatplot algo using heatmap
distCor <- function(x) as.dist(1-cor(x))


cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}


## Keep only Significant :
signif.co = 0.05
l2fc.co = 1

#> In Complete Table
D.sig <- subset(DE_tab,(padj_DM1vsCTRL <= signif.co &  abs(log2FoldChange_DM1vsCTRL) >= l2fc.co) )

## On these, subset the signif. of MBNLd-vs-CTRL:
D.sig.dct <- subset(D.sig, (padj_MBNLdvsCTRL <= signif.co &  abs(log2FoldChange_MBNLdvsCTRL) >= l2fc.co))

#> in counts
D_counts.sig <- D.sig[,colnames(norm_dds_counts)]

##> Z-transform Counts to better illustrate in Heatmap :
D_counts.sig.z <- t(apply(D_counts.sig, 1, cal_z_score))

#> remove NAs if any
##DEcounts.sig_norm <- DEcounts.sig_norm[-which(is.na(DEcounts.sig_norm[,1])),]


##> Prepare Hierarchical clustering of Rows/Genes
genes_hclust <- hclust(dist(D_counts.sig.z), method = "complete")

clust_assign <- cutree(tree = genes_hclust, k = 3)
clust_name   <- data.frame(cluster = ifelse(clust_assign == 1,yes="cluster 1",
                                            ifelse(clust_assign == 2,yes="cluster 2",
                                                   no = "cluster 3")))

sample_annot <- colData.dds

plot_name <- paste0("ExpressionPHeatmap_log2fc1_padj0.05_",nrow(D_counts.sig) ,"genes.png")

png(paste0(plots_dir,plot_name))
    pheatmap(D_counts.sig.z,
             color = colorRampPalette(rev(brewer.pal(n=10,name ="RdYlBu")))(100),
             annotation_col=sample_annot,
             clustering_distance_rows="euclidean",
             cluster_cols=FALSE,
             annotation_legend = TRUE,
             annotation_names_col = TRUE,
             show_rownames=FALSE,
             labels_row = NULL,
             angle_col = 45,
             main =paste0(nrow(D_counts.sig.z)," Signif.Diff/ly expressed Genes in DM1-vs-CTRL"))

dev.off()
