#!/usr/bin/env R

### run in R 3.6.3
### updated in 17/12/20 : R 4.0.2
### Author : Naira Naouar (2017)
### Edited by : Maria Kondili (2020)

### Subject :Detect Alternative Splicing Exons in mouse skeletal muscles of MBNL.Delta vs AAV-GFP vs Saline (CTRL)


library(yaml)
#suppressPackageStartupMessages(library(made4)) #> multiv.analysis for microarrays
#suppressPackageStartupMessages(library(genefilter))
#library(pheatmap)
library(gplots)
library(RColorBrewer)
library(BiocParallel)
source("FunctionsAltSplice_Mbnldelta.R")
source("custom_functions_mouse.R")
suppressPackageStartupMessages(library(DEXSeq)) #v1.32
##> tutorial with 1.37 :
## https://bioconductor.org/packages/devel/bioc/vignettes/DEXSeq/inst/doc/DEXSeq.html#5_Testing_for_differential_exon_usage


l2fc_cutoff = 1.5
muscle = "dct"

##-------- Paths of data files to process
data_dir = "/shared/projects/mbnl_dct/mouse_seq/"
# "D:/XENA/myWorkspace/PROJECTS/MBNL_DCT/mouse_analysis/"
work_dir = paste0(data_dir,"DEXSeq_analysis/")
counts_dir = paste0(data_dir,"Galaxy/Cutadapt/HISAT2_aligned/DEXSeq_counts/")

## @ counts_dir ,renamed in bash :
# for f in ./*.fastq.tabular; do rename .fastq.tabular _dexseq_counts.txt  $f ; done

# Samples files
sampleInfoFile <- paste0(work_dir,"metadata_MBNLdelta_dexseq.txt")
sampleInfo  <- read.table(sampleInfoFile, header=TRUE,as.is=TRUE)
sampleFiles <- paste0(counts_dir, sampleInfo$Sample)


sampleTable <- data.frame(row.names = sampleFiles,
                          condition = sampleInfo$group,
                          group = paste(sampleInfo$treatment,
                                        sampleInfo$replicate,sep="_"),
                          color = sampleInfo$color,
                          treatment= sampleInfo$treatment)
# condition == group


# Parallelization des calculs
# cores = MulticoreParam(workers=7);


# Comparison variable
compa.cond = "condition"

# Stat formulas (dispersion, null model, alternative model)
formulaD = ~ sample + exon + treatment:exon + condition:exon
formula0 = ~ sample + exon
formula1 = ~ sample + exon + condition:exon
formulas = c(formulaD, formula0, formula1)


gff ="/shared/projects/mbnl_dct/mouse_seq/Galaxy_produced_DEXSeq_annot_mm10.gtf"

# Loading data in DEXSeqDataSet object
dxd <- DEXSeqDataSetFromHTSeq(countfiles = sampleFiles,
                              sampleData = sampleTable,
                              flattenedfile = gff,
                              design = formula1)

# saveRDS(dxd, "DEXseq_DataSet_FromHTSeq_mm.rds")
# dxd <- readRDS("DEXseq_DataSet_FromHTSeq_mm.rds")

## Filter low-count genes/exons :
#cnts_per_gn <- rowSums(dxd@assays@data$counts)
# dxd@assays@data$counts <- dxd@assays@data$counts[-which(cnts_per_gn < 10), ]
# dim  = 375205  x   20


## Normalization = measure relative sequencing depth

dxd <-  estimateSizeFactors(dxd)
#> dxd@modelFrameBM$sizeFactor ---> 10 diff.values
# [1] "0.69" "0.77"  "0.82" "1.073"  "1.076"
# [6] "1.096"  "1.11"  "1.11"  "1.23"  "1.25"


## Dispersion estimation = distinguish technical and biological variation
## from real effects on exon usage (dispersion formula)

dxd <- estimateDispersions(dxd, formula=formulaD)


#############################################
##                 Analysis                 #
#############################################

## DEU (Differential Exon Usage) analysis :
## test for difference in exon bin inclusion
## compared to other bins of the gene following formulas

#cores = MulticoreParam(workers = 6 )

dxd <- testForDEU(dxd, fullModel=formula1, reducedModel=formula0) # , BPPARAM=cores)
# dxd@colData@listData$condition

##> Calculate SI (remove gene expression effect)
dxd.ctrl <-  estimateExonFoldChanges(dxd,
                                    fitExpToVar = compa.cond,
                                    denominator = "CTRL") #,BPPARAM = cores)
##! 3h for cmd to finish ...!

# dxd.gfp <-  estimateExonFoldChanges(dxd,
#                                      fitExpToVar = compa.cond,
#                                      denominator = "MBNLdecoy")
#--> analysis for AAVGFP -vs- Decoy @ ase_AAVGFP_vs_decoy.R


##ERRORS/WARNINGS :
# /2/ FOR 29 gene/exons:  threw the next warning(s): the matrix is either rank-deficient or indefinite



##################
##    Results   ##
##################


dxres = DEXSeqResults(dxd.ctrl)
#dxres.gfp = DEXSeqResults(dxd.gfp)


###
### check out on results :
###

mcols(dxres)$description

# [1] "group/gene identifier"
# [2] "feature/exon identifier"
# [3] "mean of the counts across samples in each feature/exon"
# [4] "exon dispersion estimate"
# [5] "LRT statistic: full vs reduced"
# [6] "LRT p-value: full vs reduced"
# [7] "BH adjusted p-values"
# [8] "exon usage coefficient"
# [9] "exon usage coefficient"
# [10] "exon usage coefficient"
# [11] "relative exon usage fold change"
# [12] "relative exon usage fold change"
# [13] "GRanges object of the coordinates of the exon/feature"
# [14] "matrix of integer counts, of each column containing a sample"
# [15] "list of transcripts overlapping with the exon"


table ( dxres$padj < 0.05 )
# < 0.1
# FALSE     TRUE
# 265976  14576

## < 0.05 :
# FALSE     TRUE
# 268915  11637

table ( tapply( dxres$padj < 0.05, dxres$groupID, any ) )
# dxres$groupID = geneID
# FALSE  TRUE
# 3692  3765

out_dir = paste0(data_dir,"DEXSeq_analysis/Plots/")
plotMA( dxres, cex=0.8 )

## MBNL2= ENSMUSG00000022139
## MBNL1= ENSMUSG00000027763
png(paste0(out_dir,"Expression_per_Exon_Mbnl1.png"))
  plotDEXSeq( object=dxres, geneID="ENSMUSG00000027763", legend=TRUE, cex.axis=1.2, cex=0.6, lwd=2 )
dev.off()



####
#### FILTER MBNLdecoy/CTRL  - MBNLd/AAVGFP ####
###

#> keep only those signif & Up-spliced in MBNLdecoy
dxres.gfp.sig.up <- subset( dxres.gfp, (padj < 0.05 & log2FC_MBNLdecoy_vs_AAVGFP > 0) )

dxres.ctrl.sig <- subset( dxres, padj < 0.05)
#> L = 11637
dxres.ctrl.sig.filtAAV <- dxres.ctrl.sig[ rownames(dxres.gfp.sig.up),]

# length(which(rownames(dxres.ctrl.sig) %in% rownames(dxres.gfp.sig.up)))
# [1] 8884


###
### ANNOTATION
###


# source("FunctionsAltSplice_Mbnldelta.R")-> reannotate function
l2fc_cutoff = 1
padj_cutoff = 0.05
res.reann <- reannotate(dxres.ctrl.sig.filtAAV, sampleInfo,seuil_l2fc = l2fc_cutoff, seuil_padj= padj_cutoff ,muscle)
str(res.reann)

y <- res.reann[which(res.reann$significant==T),]
y <- as.data.frame(y)
#y$Gene <- gsub("\\..*","",y$ensembl_ID)


## Elimination des potentiel faux positifs :
## lignes pour lesquelless on a moins de 30 reads


## Remove genes for which Mean-reads < 30, per condition
low_reads_pos <- unique(unlist( sapply( grep("mean",colnames(y)), function(k) which(y[,k] < 30) )))
y <- y[-low_reads_pos,]
#> L = 5 648


## On élimine les bins avec -3bases
y <- y[-which(y$genomicData.width<3),] # L = 4 430

### On elimine les lignes -1 read/position => max.read.count / width
# y$maxCov <- apply(y[,grep("mean.reads",colnames(y))],1,max) = NA !!!
# y$valCov <- round(y$maxCov/y$genomicData.width,2)
# y <- y[-which(y$valCov<1),]


## On arrondit les valeurs chiffrees
nums <- which(sapply(y, is.numeric))
for(i in nums) y[,i]=round(y[,i],3)


y <- y[,-grep("dispersion|stat|rep|id$",colnames(y))]

## Ajouter le sens du splicing inclus/exclus
y$Splicing_Direction = NA
y$Splicing_Direction[y$log2FC_MBNLdecoy_vs_CTRL < 0] = "exclus"
y$Splicing_Direction[y$log2FC_MBNLdecoy_vs_CTRL > 0] = "inclus"

species="mm10"
##Ajouter les 150 premieres bases de la séquence du bin

# y$sequence=NA
# y$sequence=mapply(getSeq,as.character(y$genomicData.seqnames),y$genomicData.start,y$genomicData.end,species)


library(biomaRt)
##### Ajouter les annotations supplémentaires presentes dans Biomart
## for mm9:
# useMart('ENSEMBL_MART_ENSEMBL',dataset='mmusculus_gene_ensembl',
#         host="may2012.archive.ensembl.org")

# for mm10 :
mm10_mart <- useMart('ensembl',dataset='mmusculus_gene_ensembl',
                     host="https://dec2017.archive.ensembl.org")

# listEnsemblArchives()-> find host=

df.anns <- getBM(attributes = c("ensembl_gene_id", "external_gene_name","description", "gene_biotype"),
                 filters = "ensembl_gene_id",
                 values = y$ensembl_ID,
                 mart = mm10_mart )

colnames(df.anns) = c("ensembl_ID", "geneName", "description", "gene_biotype")


###
###  Merge DEXSeq data.frame=y  with mart.annotation of Gene.IDs
###


# df.anns$ensembl_ID = 899
# y$ensembl_ID = 1040
y.annot.ctrl = merge(y, df.anns, by="ensembl_ID")

#-> keep only genes with annot filled.in
# y.annot$ensembl_ID = 4876
# all(y.annot$ensembl_ID %in% df.anns$ensembl_ID) = TRUE !

#> previous : y.annot <- merge(y, df.anns, all.x=T) #--> creates Millions of new lines !

## Simplify the counts-columns (filename as colname)
colnames(y.annot.ctrl)[16:25] <- colnames(y)[16:25] <- paste("dexseq_counts", sampleInfo$mouse_id,sep="_")


## Add links of ENSEMBL and NCBI per Gene
y.annot.ctrl$ensembl_link <- paste0("=lien_hypertexte(\"http://www.ensembl.org/Mus_musculus/Gene/Summary?g=",
                               y.annot.ctrl$ensembl_ID,"\";\"",
                               y.annot.ctrl$geneName,"\")")

y.annot.ctrl$fasterdb_link <- paste0("=lien_hypertexte(\"http://fasterdb.ens-lyon.fr/faster/main.pl?ensembl_initial=",
                              y.annot.ctrl$ensembl_ID,"&stable_id_ensembl=", y.annot.ctrl$ensembl_ID,
	                            "&organism=mouse&bio_mol=cDNA_Only&positions=",
	                            y.annot.ctrl$genomicData.start,"-",y.annot.ctrl$genomicData.end,"\";\"",
                              y.annot.ctrl$geneName,"-",y.annot.ctrl$exon_bin,"\")")


### Ne conserver que les 10 premiers transcrits associes au gene (sinon bug d'excel)
y.annot.ctrl$transcripts <- as.character(sapply(y.annot.ctrl$transcripts, head,n=10L))


##Conserver les colonnes dans l'ordre d'interet

oCols <- c("ensembl_ID","geneName","description", "gene_biotype",
           "exon_bin","Splicing_Direction",
           "fasterdb_link","ensembl_link",
            grep("log2FC",colnames(y.annot.ctrl),value=T),
           "genomicData.width","exonBaseMean",
            grep("EUC|mean",colnames(y.annot.ctrl),value=T),
           "correction_prcnt","pvalue","padj",
           "genomicData.seqnames","genomicData.start",
           "genomicData.end","genomicData.strand","transcripts")

if ( all(oCols %in% colnames(y.annot.ctrl))  )  y.annot.ctrl <- y.annot.ctrl[,oCols]


######################################################
## Ajouter le numero de l'exon annotation fasterdb  ##
######################################################


## Read exons & genes .csv Files ?
annot_dir ="/shared/projects/mbnl_dct/annotation_files/"
species = "mouse"
exonFile <- paste0(annot_dir, species,"_exons_genomiques_bis.csv")
geneFile <- paste0(annot_dir, species,"_genes.csv")

exons_mm9 <- read.csv2(exonFile,comment.char="#",header=F,as.is=T)
colnames(exons_mm9) <- c("id_exon","id_gene","pos_on_gene",
                          "start_on_gene","end_on_gene",
                          "exon_length", "chromosome",
                          "start_chr","end_chr","strand")


## Correction des start > end
inv_start <- which(exons_mm9$start_chr > exons_mm9$end_chr)
#! No inversions in start-end !

# for(i in inv_start){
#     tmp <- exons_mm10$start_chr[i]
#     exons_mm9$start_chr[i] <- exons_mm9$end_chr[i]
#     exons_mm9$end_chr[i] <- tmp
#
# }


genes_mm9 <- read.csv2(geneFile, comment.char="#", header=F,as.is=TRUE)
colnames(genes_mm9) <- c("id_gene","stable_id_ensembl","official_symbol",
                          "synonyms","description","chromosome",
                          "strand","start_chr","end_chr","X","gene_width")

genes_mm9 <- genes_mm9[, -10] # whole column is NA


## Convert coordinates with function "mm10Tomm9"
coord_mm10 <- as.data.frame(cbind(
                    paste0("chr",y.annot.ctrl$genomicData.seqnames),
                    y.annot.ctrl$genomicData.start,
                    y.annot.ctrl$genomicData.end))

coord_mm9 <- mm10Tomm9(coord_mm10)


##
## MERGE mm10 Annotated Data with mm9 coordinates
##

# Remove lines with "unmapped" Coordinates after Conversion mm10->9
# cnvrt_dir <- "/shared/projects/mbnl_dct/Software/liftOver/"
# unmapped_out <- paste0(cnvrt_dir,"unmapped_mm10To9.bed")

## verify if all coordinates converted :
coord_mm10$id <- paste0(coord_mm10$V1,":",coord_mm10$V2, "-",coord_mm10$V3 )

if (all(coord_mm10$id %in% coord_mm9$mm10.id )) {
  y.annot.ctrl[,c("genomicData.start.mm9", "genomicData.end.mm9")] <- coord_mm9[,2:3]

} #else {
#   unm <- which(coord_mm10$id %in% coord_mm9$mm10.id)
#   y.annot.ctrl <- y.annot.ctrl[-unm, ]
#   y.annot.ctrl[,c("genomicData.start.mm9", "genomicData.end.mm9")] <- coord_mm9[,2:3]
# }



##
## Add Exon-Number in y-mm9-Table
##


if (length(grep("exon_number",colnames(y.annot.mm9)))==0)  {
  # chunks <- seq(1,nrow(y.annot)-99 )
  # exonNumberChunks<-c()
  # for (a in chunks){ b=a+99; exonNumber(gene.id[a:b],name[a:b],..) }

  y.annot.ctrl$exon_number <- unlist(sapply(1:nrow(y.annot.ctrl),
                                           function(i) {
                                             exonNumber(db_genes = genes_mm9,
                                                        db_exons = exons_mm9,
                                                        gene_id=y.annot.ctrl$ensembl_ID[i],
                                                        name=y.annot.ctrl$geneName[i],
                                                        chromosome=y.annot.ctrl$genomicData.seqnames[i],
                                                        posStart=y.annot.ctrl$genomicData.start.mm9[i],
                                                        posEnd=y.annot.ctrl$genomicData.end.mm9[i])
                                           },USE.NAMES=FALSE))

  # add new column
  #y.annot.mm9$exon_number <- exonNumbers
  #rm(exonNumbers)
}


##
## Annotate exons with position : 3',5' UTR, Constitutive, TSS
##

library(AnnotationDbi)
library(biomaRt)
# Reminder : mm10_mart <- useMart('ensembl',dataset='mmusculus_gene_ensembl',
#                    host="https://dec2017.archive.ensembl.org")

# if(length(grep("annotation",colnames(y.annot)))==0) {
#   # y.annot <- as.data.frame(append(y.annot,list(annotation=NA),after = 3))
#   y.annot.ctrl$annotation <- sapply(1:nrow(y.annot),function(i) {
#                                               getGeneAnnotation(
#                                                 myMart = mm10_mart,
#                                                 gene = y.annot.ctrl$ensembl_ID[i],
#                                                 binStart = y.annot.ctrl$genomicData.start[i],
#                                                 binEnd = y.annot.ctrl$genomicData.end[i] )
#   })
#
# }


posLog <- grep("log",colnames(y.annot.ctrl))
# 9 10

##Ordonner les lignes en fonction du Fold Change MBNLdecoy_vs_Ctrl

y.annot.ctrl <- y.annot.ctrl[order(abs(y.annot.ctrl[,posLog[2] ]),decreasing=T),]

write.table(y.annot.ctrl, "DEXSeq_Results_MBNLΔvsCTRL_FiltAAV_.txt",
             na="NA", sep="\t", row.names=F, quote=F,dec=".")


############################
## Add EXPRESSION  ##
############################


#expressionFile <- paste0(getwd(),"../DESeq2_analysis/AllDataExpression_log2fc0_MBNLdecoy_vs_AAVGFP_vs_CTRL.txt")


# cores = MulticoreParam(workers=2)
# custom.restau.DEXSeqHTML(project, res.reannotated, dxres,
#                          genes=head(as.character(unique(y.annot$ensembl_ID))),
#                          expressionFile=expressionFile,
#                          formulas=formulas, mart = mm10_mart,
#                          species="mm10", filter="ensembl_gene_id",
#                          attributes=c("external_gene_name",
#                                       "description", "gene_biotype"),
#                          fitExpToVar=compa.cond, BPPARAM=cores)



#########################################
##                PLOTS                ##
#########################################

## Splicing Heatmap

l2fc_cutoff <- 1
regALL <- which(abs(y.annot.ctrl[,"log2FC_MBNLdecoy_vs_CTRL"]) >=l2fc_cutoff)

upreg <- which(y.annot.ctrl[,"log2FC_MBNLdecoy_vs_CTRL"] >= l2fc_cutoff)
downreg <- which(y.annot.ctrl[,"log2FC_MBNLdecoy_vs_CTRL"] <= (-l2fc_cutoff) )

##.... if both log2FC columns are considered .....
# regALL <- unique(which(abs(y.annot[,9]) >= log2(l2fc_cutoff) ) ,
#                  which(abs(y.annot[,10]) >= log2(l2fc_cutoff) ))



DFsplices <- paste0(y.annot.ctrl[regALL,1], ":", y.annot.ctrl$exon_bin[regALL])

# On enleve les MBNL !
# DFsplices <- DFsplices[-grep("ENSMUSG000000152601*",DFsplices)] --> not found

dfCounts.ctrl <- counts(dxres.ctrl.sig.filtAAV, normalized=T)
colnames(dfCounts.ctrl) <- colnames(y)[16:25]

spliced_counts <- dfCounts.ctrl[DFsplices,]

## W/O MBNLd_r3 :
spliced_counts <- spliced_counts[,-9]

## W/O AAVGFP :
spliced_counts_2cond <- spliced_counts[,-c(4:6)]
samples_labels <- sampleTable$group[-c(4:6,9)]

z_row <- zClust(x=spliced_counts_2cond, scale="row",method="average")
#order.clust=rownames(z$data)[rev(z$row)]
mycolhc <- c("dodgerblue3","red4","purple")
mycolhc <- mycolhc[as.vector(z_row$tree)]
## col = rev(colorRampPalette(brewer.pal(11, "RdYlBu"))(256)), or "RdBu"

suppressPackageStartupMessages(library(gplots))
pair <- "MBNLd_vs_CTRL"
#png(paste0(out_dir,"Heatmap_Splicing_ExonBins_",nrow(z_row$data),"_",pair,"_l2fc_",l2fc_cutoff,".png"))
  heatmap.2(z_row$data,
            dendrogram = 'row',trace='none',
            Rowv=z_row$Rowv, Colv=F, labRow="",
            labCol= samples_labels,
            col=rev(colorRampPalette(brewer.pal(8, "RdBu"))(256)),
            RowSideColors=mycolhc, notecex=3,
            key=TRUE, key.title="Color range",key.xlab="",
            key.ylab="Counts",density.info="density",
            key.par=list(mar=c(4,4,4,1)),
            offsetCol=0,margins=c(9,2),cexCol=1.3,
            main=paste(nrow(z_row$data),"Diff/ly Included Exon.bins \nin " ,
                       pair, "\n(|log2fc|>",l2fc_cutoff, ")"))

#dev.off()



z_col <- zClust(x=spliced_counts_2cond, scale="col",method="average")

#png(paste0(out_dir,"Heatmap_Splicing_Exons_bins_",nrow(z_col$data),pair,"_l2fc_",l2fc_cutoff,".png"))
heatmap.2(z_col$data,
          dendrogram ='col',trace='none',
          Rowv=F, Colv=z_col$Colv, labRow="",
          labCol=samples_labels,
          col=rev(colorRampPalette(brewer.pal(8, "RdBu"))(125)),
          notecex=3, key=TRUE, key.title="Color range",
          key.xlab="",key.ylab="Counts",density.info="density",
          key.par=list(mar=c(3,4,3,1)),
          offsetCol=0,margins=c(9,2),cexCol=1.3,
          main=paste(nrow(z_col$data),
               "Diff/ly Included Exon.bins \nin " ,pair))

dev.off()




#######
####### MATCH EXON-BINS SPLICING TO DIFF.EXPRESSED GENES #####
#######



expr_data <- read.delim("../DESeq2_analysis/AllDataExpression_log2fc=0_MBNLd_vs_AAVGFP+allPairs.txt",
                        sep="\t",header=T,as.is=TRUE)

signif_decoy_vs_ctrl_expr <- expr_data[which(expr_data$padj_MBNLdVsCTRL < padj_cutoff ),]

geneID_spliced <- sapply(rownames(spliced_counts_2cond) , function(nm) strsplit(nm,":")[[1]][1] )


which(signif_decoy_vs_ctrl_expr$GeneID %in%  geneID_spliced) %>% length




####--------- Genes with CTG-repeats :are they Diff.Spliced ? -----------####

up_decoy_ctg_genes <- c("Ikzf3" ,  "Spon1"  ,"Wdr17", "Col6a6", "Nell1",  "Mab21l3" ,"Ly9",  "Cdh22",   "Pglyrp2")
down_decoy_ctg_genes <- c("Pygl" ,"Col24a1" ,"Rasgrf1","Greb1","Fer1l6","Ptger3" , "Cyp2s1","Cdh4" ,"Aox2", "Ranbp3l")

which(y.annot.ctrl$geneName %in% c(up_decoy_ctg_genes, down_decoy_ctg_genes))
# integer(0)


## List Corrected events

# 1/ High AAVGFP/CTRL, Low Mbnld/CTRL
# lc=which(abs(y.annot$log2FC_Treat_AAVGFP_CTRL)>=0.58 & abs(y.annot$log2FC_Treat_AAVGFP_decoy_CTRL)<=0.58)
# id_lc=paste0(y.annot$ensembl_ID[lc],"_",y.annot$exon_bin[lc])
# # L = 57
#
# # 2/ High AAVGFP/CTRL, High Mbnld/CTRL
# lnc=which(abs(y.annot$log2FC_Treat_AAVGFP_CTRL)>=0.58 & abs(y.annot$log2FC_Treat_AAVGFP_decoy_CTRL)>0.58)
# id_lnc=paste0(y.annot$ensembl_ID[lnc],"_",y.annot$exon_bin[lnc])
# # L = 63
#
#
# # 3/ Low AAVGFP/CTRL, High Mbnld/CTRL
# ln=which(abs(y.annot$log2FC_Treat_AAVGFP_CTRL)<0.58 & abs(y.annot$log2FC_Treat_AAVGFP_decoy_CTRL)>=0.58)
# id_ln=paste0(y.annot$ensembl_ID[ln],"_",y.annot$exon_bin[ln])
#L = 2161



correctionTable <- data.frame(Gene=y.annot.ctrl[,1],
                              ExonBin=y.annot.ctrl$exon_bin,
                              log2D=y.annot.ctrl[,"log2FC_AAVGFP_vs_CTRL"],
                              log2T=y.annot.ctrl[,"log2FC_MBNLdecoy_vs_CTRL"],
                              correction_pcnt=y.annot.ctrl[,grep("corr",colnames(y.annot.ctrl))])

correctionTable$percent <- round(100*(correctionTable$log2D - correctionTable$log2T)/
                                   correctionTable$log2D,2)

# just serial.numbers of lines= rownames
rownames(correctionTable) <- rownames(y.annot.ctrl)



exonbin_CorrectionStatus <-  function(exon_bin, fc.co) {

  ##' Run for each line= exon_bin
  ##' Input considered given: correctionTable

  minCorrec=20
  maxCorrec=50

  log2D  <- abs(correctionTable[exon_bin,"log2D"]) # log2FC_AAVGFP_vs_CTRL
  log2T  <- abs(correctionTable[exon_bin,"log2T"]) # log2FC_MBNLdecoy_vs_CTRL
  correc <- abs(correctionTable[exon_bin,"correction_pcnt"] )

	if (fc.co >1) {
	  l2fc.co <- abs(log2(1/fc.co))
	  print( c(log2D,log2T,correc) )
   }
	if (log2D < l2fc.co & log2T >= l2fc.co)
	    return("New")

	if(log2D>=l2fc.co & log2T>=l2fc.co & log2D<log2T){
	    return("No_Correction")
	}
	if(log2D>=l2fc.co & log2T>=l2fc.co & correc <= minCorrec){
	    return("Bad_Correction")
	}
	if(log2D>=l2fc.co & log2T>=l2fc.co & correc > minCorrec & correc <= maxCorrec){
	    return("Mild_Correction")
	}
	if(log2D>=l2fc.co & log2T>=l2fc.co & correc > maxCorrec) {
	  return("Good_Correction")
	}
	if(log2D>=l2fc.co & log2T<l2fc.co){
	  return("Corrected")
	}
	if(log2D<l2fc.co & log2T<l2fc.co){
    return("Under_Threshold")
	}
	return("Problem_not_affected")
}


##> Apply function CorrectionStatus for each line = each rowname
fc_cutoff=2
minCorrec=20
maxCorrec=50


correctionTable$Status    <-  sapply(rownames(correctionTable),
                                     exonbin_CorrectionStatus,
                                     fc.co = fc_cutoff,
                                     simplify = "array")


cam <- data.frame(table(correctionTable$Status))
# cam$percent <- round(cam$Freq/(sum(cam$Freq))*100,2)


## On fait le camembert en regroupant No+Bad / Mild+Good / Corrected

experiment="Splicing"
treatment="dct"
groups <- c("Complete","Partial","Bad")
# groups <- c("Complete Correction",(paste0(" Correction >",minCorrec,"%")),
#            (paste0("Correction <=",minCorrec,"%")))


#piedata<-data.frame(Group=c("Corrected","Partial_Correction","Bad_Correction"),Freq=c(cam$Freq[cam[,1]=="Corrected"],
piedata <- data.frame(Group=groups,Freq=c(cam$Freq[cam[,1]=="Corrected"],
	                    sum(cam$Freq[cam$Var1=="Good_Correction"],
	                        cam$Freq[cam$Var1=="Mild_Correction"]),
	                    sum(cam$Freq[cam$Var1=="No_Correction"],
	                        cam$Freq[cam$Var1=="Bad_Correction"])))


#colors <- c('rgb(211,94,96)', 'rgb(128,133,133)', 'rgb(144,103,167)')
#colors <- c('rgb(93,71,139)','rgb(171,130,255)','rgb(218,165,32)')
colors <- c('rgb(153,51,255)','rgb(178,102,255)','rgb(229,204,255)')


# createCam <- function(piedata){
#     p <- plot_ly(piedata, labels = ~Group, values = ~Freq, type = 'pie',
#                 textposition = 'inside',
#                 textinfo = 'label+text',
#               	text=paste0(piedata$percent,"%"),
#                 insidetextfont = list(color = '#000000',size=18),
#                 # hoverinfo = 'text',
#                 # text = ~paste(Freq, ' exon bins'),
#                 marker = list(colors = colors,
#                 line = list(color = '#FFFFFF', width = 0)),
#                 showlegend = FALSE) %>%
#                 #The 'pull' attribute can also be used to create space between the sectors
#                 layout(title = paste0('Correction rates for ',experiment,' in ',muscle),
#                        xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
#                        yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))
#     return(p)
# }
#
# #p <- createCam(piedata2)
# # export(p, file = paste0("Pie",experiment,"_",minCorrec,"-",maxCorrec,"_",muscle,".png"))
# p <- createCam(piedata)
# export(p, file = paste0("Pie",experiment,"_",muscle,".png"))


## NEW FUNCTION by M.K.
create_pie <- function(piedata) {

  suppressPackageStartupMessages(library(ggplot2))
  pie <- ggplot(piedata, aes(x="",y=Freq,fill=factor(Group))) +
          geom_bar(width=1,stat="identity") +
          theme_classic() +
          theme(axis.line = element_blank(),axis.text = element_blank(),axis.ticks = element_blank()) +
          geom_text(aes(label = paste(Group,":", round(Freq / sum(Freq) * 100, 1), "%"),
                        x = 1.1), position = position_stack(vjust = 0.6)) +
    labs(fill="Correction", x=NULL,y=NULL,
               title="Exon Correction in DCT") +
          coord_polar(theta = "y", start=0) +
          scale_fill_brewer(palette="Dark2")
  return(pie)

}


## call function:
p <- create_pie(piedata)

png("Correction_PieChart_MBNLd_vs_Ctrl.png")
print(p)
dev.off()
