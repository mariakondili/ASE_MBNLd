#!/usr/bin/env R

###~~~ Run DESeq2 for gene-HTSeq-counts, after STAR.Alignment ~~~~###
### Author: Maria Kondili 

## Paths @IFB.core.cluster

data_dir <- "/shared/projects/mbnl_dct/mouse_seq/"

x <- read.table(paste0(data_dir, "metadata_MBNLdelta_samples.txt"),
        header=T,stringsAsFactors=F, as.is=TRUE)

## Create a table of counts,merging all samples
## genes in lines // samples in columns

counts_dir <- paste0(data_dir,"Galaxy/Cutadapt/HISAT2_aligned/HTSeq_counts/") 


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

write.table(counts_data, file=paste0(data_dir,"DESeq2_analysis/Merged_CountData_AAVGFP_MBNLd_Hisat2_trimmed.txt"),
            row.names=F, sep="\t")
## L=  55401

#counts_data <- read.table(paste0(data_dir,"DESeq2_analysis/Merged_CountData_AAVGFP_MBNLd_Hisat2_trimmed.txt"), h=T)
counts_data <- counts_data[,-1] # keep rownames but ignore 1st col with geneID
rownames(counts_data) <- genes_IDs
projet <- "dct_mm"
nbreads <- rowSums(counts_data)

# Filter Low-count Reads
counts_data.filt <- counts_data[-which(nbreads<30),] # 20140 genes left

library(DESeq2)
library(biomaRt)

mart <- useMart("ensembl", dataset="mmusculus_gene_ensembl");
## doc = https://rdrr.io/bioc/biomaRt/f/vignettes/biomaRt.Rmd


## matrice de design
condition  <- gsub("_.*", "", colnames(counts_data.filt))


dctDesign <- data.frame(row.names=colnames(counts_data.filt), condition ) 

# no colname ="condition" needed 
#It should be TRUE : all(rownames(dctDesign) ==  colnames(counts_data)), for dds

dds <- DESeqDataSetFromMatrix(countData = counts_data.filt, colData=dctDesign, design=~condition)
des <- DESeq(dds)


# deseq_results_per_cond <-  function(dds, condition, treat, ctrl) { 
# 
#     res <- results(dds, contrast=c("condition",treat, ctrl))
#     colnames(res)[-1] <- paste0(colnames(res)[-1], paste0("_", ctrl,"vs",treat ))
#     return(res)
# }
# res.dm1 <- deseq_results_per_cond(data, cond.dm1, "DM1","WT")
# res.mbnld <- deseq_results_per_cond(data, cond.mnbld , "MBNLd","WT")


res.treat1 <- results(des, contrast=c("condition","AAVGFP" ,"Saline"))
#res.treat1@elementMetadata$description
colnames(res.treat1)[-1] <- paste0(colnames(res.treat1)[-1],"_AAVGFPvsCTRL")

res.treat2 <- results(des, contrast=c("condition","MBNLd","AAVGFP"))
colnames(res.treat2)[-1] <- paste0(colnames(res.treat2)[-1],"_MBNLdvsAAVGFP")

res.treat3 <- results(des, contrast=c("condition","MBNLd", "Saline"))
colnames(res.treat3)[-1] <- paste0(colnames(res.treat3)[-1],"_MBNLdVsCTRL")


#! baseMean is same for both res
#! Results contain geneIDs in rownames ! 


##
## Merge Diff.analysis of 3 conditions 
## 

#y <- cbind(as.data.frame(res.treat1),as.data.frame(res.treat2)[,-1])
y <- cbind(as.data.frame(res.treat1),
            as.data.frame(res.treat2)[,-1],
            as.data.frame(res.treat3)[,-1])



# y=merge(res,res2)
# nums <- which(sapply(y, is.numeric))
# for(i in nums) y[,i] <- round(y[,i],3)


y$GeneID <- rownames(y)
DE_tab <- cbind(y, counts_data.filt)

df.anns <- getBM(attributes = c("ensembl_gene_id","external_gene_name","description","gene_biotype"),
                 filters    = "ensembl_gene_id", 
                 values     = DE_tab$GeneID, 
                 mart       = mart)


colnames(df.anns)[1:2] <- c("GeneID","geneName")

DE_tab.ann <- merge(DE_tab, df.anns,by="GeneID", all.x=T)

DE_tab.ann$significant_AAVGFPvsCTRL  <- DE_tab.ann$padj_AAVGFPvsCTRL < 0.1
# 1973 TRUE
DE_tab.ann$significant_MBNLdvsAAVGFP <- DE_tab.ann$padj_MBNLdvsAAVGFP < 0.1
# 5424 TRUE : Mbnl_delta vs AAV.GFP

DE_tab.ann$significant_MBNLdvsCTRL   <- DE_tab.ann$padj_MBNLdVsCTRL < 0.1
# 4184 TRUE: Mbnl_delta vs Saline 


common_aavgfp_delta<- with(DE_tab.ann,
                           which(significant_MBNLdvsCTRL ==TRUE & 
                                   significant_AAVGFPvsCTRL==TRUE))

common_signif_genam <- DE_tab.ann$geneName[common_aavgfp_delta]

## Calculated Significant and UP/DOWN Regulated :
signif_Down_aavgfp_vs_ctrl <- with (DE_tab.ann, 
                                    which(log2FoldChange_AAVGFPvsCTRL <0 & padj_AAVGFPvsCTRL <0.1))
signif_Up_aavgfp_vs_ctrl <- with (DE_tab.ann, 
                                    which(log2FoldChange_AAVGFPvsCTRL >0 & padj_AAVGFPvsCTRL <0.1))


signif_Down_Mnbld_vs_aavgfp <- with (DE_tab.ann, 
                                    which(log2FoldChange_MBNLdvsAAVGFP <0 & padj_MBNLdvsAAVGFP <0.1))
signif_Up_Mnbld_vs_aavgfp <- with (DE_tab.ann, 
                                  which(log2FoldChange_MBNLdvsAAVGFP  >0 & padj_MBNLdvsAAVGFP <0.1))

##> Change sign of Log2FC values --> ?? WHY?
# y$log2FoldChange_WTvsDM1=-y$log2FoldChange_WTvsDM1
# y$log2FoldChange_WTvsMBNL1d=-y$log2FoldChange_WTvsMBNL1d


# oCols <- c("GeneID","geneName",
#            "significant_AAVGFPvsCTRL", "significant_MBNLdvsAAVGFP",
#            "log2FoldChange_AAVGFPvsCTRL", "log2FoldChange_MBNLdvsAAVGFP",
#            "description", "gene_biotype",
#             colnames(counts_data.filt),"baseMean",
#            "lfcSE_AAVGFPvsCTRL","stat_AAVGFPvsCTRL", 
#            "pvalue_AAVGFPvsCTRL","padj_AAVGFPvsCTRL",
#            "lfcSE_MBNLdvsAAVGFP","stat_MBNLdvsAAVGFP",
#            "pvalue_MBNLdvsAAVGFP","padj_MBNLdvsAAVGFP")


oCols <- c("GeneID","geneName",
           "significant_AAVGFPvsCTRL","significant_MBNLdvsCTRL", "significant_MBNLdvsAAVGFP",
           "log2FoldChange_AAVGFPvsCTRL","log2FoldChange_MBNLdVsCTRL", "log2FoldChange_MBNLdvsAAVGFP",
           "description", "gene_biotype",
           colnames(counts_data.filt),"baseMean",
           "lfcSE_AAVGFPvsCTRL","stat_AAVGFPvsCTRL","lfcSE_MBNLdVsCTRL", "stat_MBNLdVsCTRL", 
           "lfcSE_MBNLdvsAAVGFP","stat_MBNLdvsAAVGFP","pvalue_AAVGFPvsCTRL","padj_AAVGFPvsCTRL",
           "pvalue_MBNLdVsCTRL", "padj_MBNLdVsCTRL","pvalue_MBNLdvsAAVGFP","padj_MBNLdvsAAVGFP")

out_dir <- paste0(data_dir,"DESeq2_analysis/")

if (all(oCols %in% colnames(DE_tab.ann))) {
  D <- DE_tab.ann[order(abs(DE_tab.ann$log2FoldChange_MBNLdvsAAVGFP),decreasing=T), oCols]
  # ! attention: for large log2FC values ,the p-values and p.adj are NA ! 
  # 19402 genes 
  
  
  write.table(subset(D, abs(D$log2FoldChange_MBNLdvsAAVGFP)>=0),
              paste0(out_dir,"AllDataExpression_log2fc=0_MBNLd_vs_AAVGFP+allPairs.txt"),sep="\t",
              row.names=F,quote=F,na="NA")

  sigD <- D[which(D$significant_MBNLdvsAAVGFP==T),] ## 1857 are TRUE
  sigD.l2fc <- subset(sigD,abs(sigD$log2FoldChange_MBNLdvsAAVGFP)>= l2fc.co)
  
  write.table(sigD.l2fc,
              paste0(out_dir,"SignificantDataExpression_log2fc=",l2fc.co,"_MBNLd_vs_AAVGFP.txt"),
              sep="\t",row.names=F,quote=F,na="NA")
  
  write.table(sigD.l2fc$GeneID, 
              paste0(out_dir,"EnsemblID_SignifGenes_log2fc=",l2fc.co,"_MBNLd_vs_AAVGFP.txt"),
              sep="\t",row.names=F,quote=F,na="NA")
  
}



##----- Find DE.Genes of CTRL.base , that are D.E. in Decoy-vs-AAVGFP----------##

genes_l2fc2.ctrl <- D$geneName[intersect(which(abs(D$log2FoldChange_MBNLdVsCTRL)  > 2), 
                                       which(abs(D$log2FoldChange_AAVGFPvsCTRL) > 2) )]

genes_l2fc2.decoy <- intersect ( genes_l2fc2.ctrl,
                               D$geneName[which(abs(D$log2FoldChange_MBNLdvsAAVGFP) > 2)])
#186 genes 



##------- PCA -----------##
plots_dir <-paste0(out_dir, "Plots/")

library(pheatmap)
library(gplots)

# png(paste0(out_dir,"PCA_of_samples.png"))
#   rld <- rlog(des,blind=F)
#   DESeq2::plotPCA(rld)
# dev.off()


##-------------------- HEATMAP ------------------##
## muscle="Gas" ?? 


# Necessary function for heatplot algo using heatmap
distCor <- function(x) as.dist(1-cor(x));


# Necessary function for heatplot algo using heatmap
zClust <- function(x, scale="row", zlim=c(-3,3), method="average") {

    if (scale=="row"){
      z <- t(scale(t(x)))
      z <- pmin(pmax(z, zlim[1]), zlim[2])
      hcl_row <- hclust(distCor(t(z)), method=method)
      return(list(data=z, Rowv=as.dendrogram(hcl_row)))
    }
  
    if (scale=="col"){
      z <- scale(x)
      z <- pmin(pmax(z, zlim[1]), zlim[2])
      # pmin/pmax : parallel min-max of 2 vectors(here table and digit) in 1
      # pmin, pmax give same z table, so confirm it is in {-3,3}
      hcl_col <- hclust(distCor(z), method=method)
      return(list(data=z,Colv=as.dendrogram(hcl_col)))
    }
}


norm_10to6 <- function(data){
  N=NULL
  for(i in 1: ncol(data)) {
    N <- cbind(N,data[,i]/sum(data[,i])*10000000)
  }
  colnames(N)=colnames(data)
  return(N)
}



## Second Filtering of Low.counts < 5 
# lowReadsPos = c()
# for(i in 1:ncol(counts_data.filt)) {
#   lowReadsPos <- c(lowReadsPos,which(counts_data.filt[,i]<5))
# }
# counts_data.filt <-  counts_data.filt[-unique(lowReadsPos),]
## lowReadsPos = 3562


##> Keep the lines of DESeq with Significant genes
# I want FC > 1.5
# log(2* 1.5) > x 
# 2^x > 1.5
# log2(2^x) > log2(1.5)
# x > 0.58

# posSig <- intersect(which(D$significant_MBNLdecoyVsCTRL == T),
#                     which(abs(D$log2FoldChange_MBNLdecoyVsCTRL) > 0.57) )

l2fc.co  <- 2
posSig.d <- intersect(which(D$significant_MBNLdvsAAVGFP == T),
                      which(abs(D$log2FoldChange_MBNLdvsAAVGFP) > l2fc.co ))

## Detect signif.genes in counts.data :
regGenes <- intersect(rownames(counts_data.filt),D$GeneID[posSig.d]) 

# 1875 genes, with lowReadsPos
# 2435 , without lowReadsPos

## Keep the highly expr. and Signif. genes to move on..
data.sig <- counts_data.filt[rownames(counts_data.filt) %in% regGenes,]

genes.sig <- D$GeneID[posSig.d]
nb.bins <- dim(data.sig)[1]

data.sig.norm <- norm_10to6(data.sig)
##> REMOVE SAMPLE " MBNLd_r3" 
data.sig.norm <- data.sig.norm[,-9]

# select <- order(rowMeans(as.matrix(data.sig.norm )), decreasing=TRUE)[1:nb.bins];
# data.sort <- as.matrix(data.sig.norm )[select, ]; #--> same results w/o sorting !
z <- zClust(data.sig.norm, scale="row")


##------------- PCA on significant genes --------------## 

data.pca    <- prcomp(t(na.omit(as.matrix(data.sig.norm))));
PC.variance <- round(data.pca$sdev^2/sum(data.pca$sdev^2)*100, 1);

library(factoextra)

png(paste0(plots_dir,"PCA_on_SignifGenesExpression_l2fc=",l2fc.co,"samplesLabeled.png"))
    fviz_pca_ind(data.pca, habillage=condition[-9], 
                invisible="quali", geom="auto",
                pointsize=2.5, axes.linetype=0,
                xlab=paste0("PC1 = ", PC.variance[1], "%"),
                ylab=paste0("PC2 = ", PC.variance[2], "%"),
    	          title=paste0("PCA on Significant Expression Data (|log2FC| >",l2fc.co,") in murine model")) +
                theme_gray()
dev.off()


# pdf(paste0(out_dir,"PCA_on_Signif_Genes_Expression_noLabels.pdf"))
#     fviz_pca_ind(data.pca , habillage=condition,invisible="quali",
#                  geom="point",pointsize=4,axes.linetype=0,
#                  xlab=paste0("PC1 = ", PC.variance[1], "%"),
#                  ylab=paste0("PC2 = ", PC.variance[2], "%"),
#     	           title=paste0("PCA on ",datatype," expression data in murine model")) +
#       theme_gray()
# dev.off()


library(gplots)
library(devtools)

png(paste0(plots_dir,"ExpressionHeatmap_on_SignifGenes_l2fc=",l2fc.co,"MBNLdecoy_vs_AAVGFP.png"))
  #? without sorting:  z <- zClust(as.matrix(data.sig.norm ),scale="row");
  heatmap.2(z$data, dendrogram = "row", trace='none', 
            Rowv=z$Rowv, Colv=F, labRow="",
            labCol=colnames(z$data), col="greenred", 
            key.title="Color range", key.xlab="",key.ylab="Counts",
            margins=c(8,2), cexCol=1.3, key.par=list(mar=c(3,4,3,0)),
            main=paste(nrow(z$data),"Diff/ly Expressed Genes "))
dev.off()

# 'RowSideColors' must be a character vector of length nrow(x) 



## CLUSTER BY COLUMNS / SAMPLES --> Show that AAV similar to Decoy

col.hc <- hclust(as.dist(1-cor(data.sig.norm)),method="median")
row.hc <- hclust(as.dist(1-cor(t(data.sig.norm) )),method="median")

#method="ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), 
#       "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC)

myhclust <- function(x) hclust(x, method="median")
mydist <- function(x) dist(x,method="euclidean") # euclidean, "maximum", "manhattan", "canberra", "binary" or "minkowski"

png(paste0(plots_dir,"Heatmap_by_column_NoDendogr_data_sig_norm.png"))
library(gplots)
heatmap(as.matrix(data.sig.norm),
        Colv= as.dendrogram(col.hc),
        Rowv = NA, labRow="",
        distfun = mydist,
        labCol = colnames(data.sig.norm),
        keep.dendro=FALSE,
        col=greenred(256),
        margins=c(8,2),cexCol=1.3 )

dev.off()


png(paste0(plots_dir,"Heatmap_by_column+RowDendogr_data_sig_norm.png"))
heatmap(as.matrix(data.sig.norm),
        Colv= as.dendrogram(col.hc),
        Rowv = as.dendrogram(row.hc),
        labRow = "",
        distfun = mydist,
        labCol = colnames(data.sig.norm),
        col=greenred(256),
        margins=c(8,2),cexCol=1.3 )

dev.off()


###------ CROSS-CORRELATION OF SAMPLES -----------------###
library(Hmisc)
library(RColorBrewer)

data.sig.rcorr = rcorr(as.matrix(data.sig.norm),type="spearman")


pdf("Plots/Cross-Correlation_Heatmap_10samples.pdf")
heatmap(x = data.sig.rcorr$r, 
        col = colorRampPalette(brewer.pal(8, "YlOrRd"))(25), 
        symm = TRUE,cexCol=0.9,cexRow=0.9,
        labCol=rownames(data.sig.rcorr$r),
        main="Spearman Cross-Correlation of Samples")
 dev.off()

 

 ## data: AAVGFP & MBNLdecoy only : 
 # data.sig.2c <- data.sig.norm[,-c(1:3,9)]
 # data.2c.rcorr <- rcorr(as.matrix(data.sig.2c),type="spearman")
 
heatmap(x = data.2c.rcorr$r, 
         col = colorRampPalette(brewer.pal(8, "YlOrRd"))(25), 
         symm = TRUE,cexCol=0.9,cexRow=0.9,
         labCol=rownames(data.2c.rcorr$r))
 ## find colourPalette : https://www.r-graph-gallery.com/215-the-heatmap-function.html
 


###--------- MBNL1-regulated_genes (Batras et al, supp.data) ---------------###

batras_diff_splice_genes <- read.delim("../diff_Exon_Alt_splicing_Genes_Batras_suppTab2.txt",
                                      sep="\t",header=F,as.is=T)

batras_diff_splice_genes<- sapply(batras_diff_splice_genes[,1], function(g) 
                                  gsub(".+\\//" ,"",g),
                                  simplify = "array")


