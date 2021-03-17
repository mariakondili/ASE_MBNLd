

###############################################################
#                          Functions                          #
###############################################################

options(scipen = 999)

getSeq <- function(chr,start,end,database){

  if(database=="mm10") db="/media/naira/Data/Genomes/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa"
  
  if(start>end){ tmp=start;start=end;end=tmp;}
  
  if(end-start>150){ end=start+150;}
  return(paste(system(paste0("blastdbcmd -db ",db," -dbtype nucl -entry ",chr," -range ",start,"-",end),intern=T)[-1],collapse=""))

}


################################
### Convert coordinates to mm9 
### create BED file from coordinates and return coordinates converted
################################

mm10Tomm9 <- function(coordinates){
  
  #' Input: coordinates should be from annot.genes of mm10 in data.frame of 3 columns : chr , start, end 
  #' Usage:
  #' $ ./liftOver oldFile map.chain newFile unMapped
  #' To get map.chain : 
  #' wget --timestamping ftp://hgdownload.cse.ucsc.edu/goldenPath/mm10/liftOver/mm10ToMm9.over.chain.gz
 
  cnvrt_dir <- "/shared/projects/mbnl_dct/Software/liftOver/"
  mm10_file  <- paste0(cnvrt_dir,"coordmm10_mart.bed")
  
  if ( ! grepl("chr",coordinates[1,1])) {  
      coordinates[,1] <- paste0("chr", coordinates[,1])
  }
  
  #Add a 4th column as id 
  coordinates$id <- paste0(coordinates[,1],":",coordinates[,2],"-",coordinates[,3])
  
  #> save in file to use in tool "liftOver"
  write.table(coordinates,paste0(cnvrt_dir,"coordmm10_mart.bed"), 
              sep="\t", col.names=F, row.names=F, quote=F)
  
  map.chain  <- paste0(cnvrt_dir,"mm10ToMm9.over.chain.gz")
  mm9_out    <- paste0(cnvrt_dir,"coordmm9_lift.bed")
  unmapped_out <- paste0(cnvrt_dir,"unmapped_mm10To9.bed")
  
  
  #launch bash command : 
  system(paste(paste0(cnvrt_dir,"liftOver"),mm10_file, map.chain, mm9_out, unmapped_out, 
               sep=" ")) #,ignore.stderr = T)
  ## Online version to run: http://genome.ucsc.edu/cgi-bin/hgLiftOver 
  cat(mm9_out, "was created !\n")
  coordmm9 <- read.delim(mm9_out, header=FALSE)
  colnames(coordmm9) <- c("chr", "start", "end", "mm10.id")
  
  # if (file.info(unmapped_out)$size == 0) {
  #       coordinates <- coordmm9
  # } else {
  #     unm <- read.delim(unmapped_out, comment.char="#", h=F)
  #     colnames(unm) <- c("chr", "start", "end", "id")
  #     pos <- which(coordinates$id %in% unm$id)
  #     ## Create new data.frame "coordinates" with mm9 start-end
  #     coordinates <- coordinates[-pos,]
  #     coordinates <- coordmm9
  # }
  
  
  ##> Return the converted coordinates of mm9
  coordinates <- coordmm9  # or : coordmm9[,1:3], if I don't want the ID column
  return(coordinates)
}


###
### Get exonNumber from genomic positions
###


## call : 
# db_genes = genes_mm9; db_exons=exons_mm9 ;
# gene_id=y.annot.mm9$ensembl_ID[199]
# name=y.annot.mm9$geneName[199]
# chromosome=y.annot.mm9$genomicData.seqnames[199]
# posStart=y.annot.mm9$genomicData.start[199]
# posEnd=y.annot.mm9$genomicData.end[199]

## exonNumber( genes_mm9,exons_mm9, gene_id ,name, chromosome,posStart,posEnd)

exonNumber <- function(db_genes, db_exons, 
                       gene_id=gene_id, name=gene_name,
                       chromosome=gene_chr,
                       posStart=gene_start,
                       posEnd=gene_end)
{
  ## Requires:
  ##      1 . genes_mm9  & exons_mm9: from files "genes.csv" & "exons_genomiques_bis.csv"
  ##      2. The genomic coordinates from the y.data.frame,after merging with converted mm10Tomm9().
  

  db_genes <- as.data.frame(db_genes) 
  db_exons <- as.data.frame(db_exons)

	if ((is.na(posStart)) & (is.na(posEnd))) {
	  cat("\n Start or End position is NA !")
	  return(NA)
	}
  #> from "synonyms" only ensembl_IDs will be matched 
  pos_found <- unique(c(which(db_genes$synonyms  == gene_id ) , 
                        which(db_genes$stable_id_ensembl == gene_id)))
  
  cat("\nGene_ID found in pos: ", pos_found)
  ## Case 1 : Gene_Id Matches FasterDB ids
  if  (length(pos_found) != 0) 	{
      tmpGeneId <- db_genes$id_gene[pos_found]
      cat("\ntmpGeneId:", tmpGeneId)
  ## Case 2: Gene_Id Doesn't match -> search GeneName/symbol
  }else{ 
    tmpGeneId <- db_genes$id_gene[unique(which(db_genes$official_symbol==name))]
  }
  
	## Case 1 @ tmpGeneId 
	if ( length(tmpGeneId) >= 1 ) {
	      ## Case 1.a @ tmpGeneId : Multiple GeneIDs matching 
	      if (length(tmpGeneId) > 1 ) {
	          cat("\ntmpGeneIDs found are:", tmpGeneId)
	          cat("\nMore than 1 geneID for ",gene_id, ".Only First will be considered!")
	          tmpGeneId <- tmpGeneId[1]
	      }
	      ## Case 1.b @ tmpGeneId : Just the one Ideal Gene.ID found || The first was used !
    	  cat("\ntmp GeneID is unique : ", tmpGeneId)
    	  tmpExon <- db_exons[db_exons$id_gene==tmpGeneId, ]
    	  #head(tmpExon)
    	  ## remove "chr" prefix -> IRanges doesn't accept it 
    	  chromosome <- ifelse(grepl("chr",chromosome), gsub("chr", "", chromosome), chromosome )
    	  rangesBin <- split(IRanges(start=posStart,end=posEnd),chromosome) 
    	  rangesExons <- split(IRanges(start=tmpExon$start_chr, end=tmpExon$end_chr),chromosome) 
    	  ov <- findOverlaps(rangesBin, rangesExons, type="any")
    	  targets <- as.matrix(ov)[,2]
    	  
    	  if (length(targets) == 1)  { return(tmpExon[targets,"pos_on_gene"]) }
    	  if (length(targets) >  1)  { return(paste0(unique(tmpExon[targets,"pos_on_gene"]),collapse=","))
    	  } else { cat("\nNo Targets overlapping !") ; return(NA) }
    	    
  ## Case 3 @ tmpGeneId : No Gene.ID found 
	} else { # length(tmpGeneId) < 1 ==0 
	      cat("\n", gene_id," not retrieved in FasterDB ")
    		return(NA)
		    }
}


###
### Get gene annotation from genomic positions 
###


# call : 
# gene = as.character(y.annot.mm9$ensembl_ID)[i]
# binStart = y.annot.mm9$genomicData.start[i]
# binEnd = y.annot.mm9$genomicData.end[i]

getGeneAnnotation <- function(myMart, gene, binStart, binEnd){
  
	## Extract Info of constit.genes & UTRs 
	attributes=c("ensembl_gene_id",
              "exon_chrom_start","exon_chrom_end",
              "is_constitutive","5_utr_start","5_utr_end",
              "3_utr_start","3_utr_end")
	## to see possible attributes: listAttributes(mart)
	geneExons <- getBM(attributes=attributes,
	                   filter="ensembl_gene_id",
	                   value=gene, 
	                   mart = myMart)  
	
	# colnames(geneExons) <- attributes
	rangesBin   <- IRanges(binStart,binEnd)
	rangesExons <- IRanges(geneExons[,'exon_chrom_start'],geneExons[,'exon_chrom_end'])
	ov <- as.matrix(findOverlaps(rangesBin, rangesExons, type="any"))
  cat( "\nOverlap of gene-region and db_exons was researched..")
  
  if (length(ov) == 0) cat("\nNo Overlap was found !")
  
	annotation <- c("")
	if(nrow(ov)>0) {
	    cat("\nThere are Overlaps!")
    	annot <- geneExons[ov[,2],]
    
    	annotation <- ifelse(length(which(annot[,'5_utr_start']>0))>0, '5utr', '')
    	annotation <- paste0(annotation,ifelse(length(which(annot[,'3_utr_start']>0))>0, '3utr', ''))
    	annotation <- paste0(annotation,ifelse(sum(annot[,'is_constitutive'])<1, '' , 'constitutive'))
    
    	## extract TSS info per gene
    	attributes <- c("transcription_start_site")
    	
    	tss <- as.vector(t(getBM(attributes=attributes,
    	                         filter="ensembl_gene_id",
    	                         value=gene ,
    	                         mart= myMart)))
    	
    	if(length(tss)>0) {
    		for(i in 1:length(tss)){
    			if(binStart <= tss[i] & binEnd >= tss[i]){
    				annotation <- paste0(annotation,",tss")
    				break
    	 }}}
    
    	annotation=gsub("^,","",annotation)
    	annotation=gsub(",$","",annotation)
	}
  return(annotation)
}

#####
### Remove NAs in specific column


removeNA<- function(obj,nbCol){
    posNA <- which(is.na(obj[,nbCol]))
    #print(length(posNA))
    if(length(posNA>0)) return(obj[-posNA,])
    else return(obj)
}



reannotate <- function (res,sampleInfo,seuil_l2fc = 1,seuil_padj=0.05 , muscle)
{

    mean.reads.min <- 30   # min read coverage for an exon bin in one condition

    #res=res.reannotated
    # Rename columns for clarity
    colnames(res)[1]            = "ensembl_ID";
    colnames(res)[2]            = "exon_bin";
    colnames(res)[8:10]         = paste("EUC",colnames(res)[8:10],sep="_");
    colnames(res)               = gsub("log2fold","log2FC",colnames(res));
    counts_data                 = as.matrix(res$countData)
    colnames(counts_data)       = paste0(sampleInfo$treatment,"_",sampleInfo$replicate)
    
    # for (i in 1:ncol(counts_data)){
    #     nco <- colnames(counts_data)[i]
    #     pos <- grep(substr(colnames(counts_data)[i],1,6),sampleInfo$treatment)
    # 		  
    #     colnames(counts_data)[i] <- paste0(substr(nco,1,3),"_","nb.reads_",
    #                                      sampleInfo$treatment[pos],"_",
    #                                      sampleInfo$replicate[pos])
    # }

    # Add a column "mean.reads" per condition
	  for(i in (unique(sampleInfo$treatment))) {
		  #for(j in muscle){
		  res[,paste0(i,"_mean.reads")] <- rowMeans(counts_data[,grep(i,colnames(counts_data))])
	  }

    # Add column "id" as "ensembl_ID:exon_bin"
    res[,"ExonID"] <- with(res, paste(ensembl_ID, exon_bin, sep = ":"))

    # Add regulation columns
    # res$reg_DM1_C25    <- as.factor(res$id %in% df.exons.sel$id);
    # if ( project == "u7" ){
    #      res$reg_DM1_U7_C25 <- as.factor(abs(res$log2FC_DM1_U7_C25) >= seuil.fc & res$padj <= seuil_padj & ( res$mean.reads.C25 >= mean.reads.min | res$mean.reads.DM1_U7 >= mean.reads.min ) );
    # } else if ( project == "dct" ){
    #     res$reg_DM1_dCT_C25 <- as.factor(abs(res$log2FC_DM1_dCT_C25) >= seuil.fc & res$padj <= seuil_padj & ( res$mean.reads.C25 >= mean.reads.min | res$mean.reads.DM1_dCT >= mean.reads.min ) );


	  res[,"significant"] <- ((res$padj <= seuil_padj) &
                              !(is.na(res$log2FC_AAVGFP_vs_CTRL))  & 
                              !(is.na(res$log2FC_MBNLdecoy_vs_CTRL)))
    
	  res[ which(is.na(res$significant)), "significant"] <- FALSE
    res[,"reg_Mbnld_Ctrl"]  <- abs(res$log2FC_MBNLdecoy_vs_CTRL)>= seuil_l2fc & res$significant
    res[,"reg_AAVGFP_Ctrl"] <- abs(res$log2FC_AAVGFP_vs_CTRL)>= seuil_l2fc & res$significant
    res[,"regALL"]          <- (res$reg_Mbnld_Ctrl | res$reg_AAVGFP_Ctrl)
    # Add correction %
    res[,"correction_prcnt"] <-  -(res$EUC_MBNLdecoy - res$EUC_AAVGFP) / 
                                  (res$EUC_AAVGFP - res$EUC_CTRL) * 100;


    # Dealing with NA values (necessary for graphics)
    #res$reg_DM1_C25[is.na(res$reg_DM1_C25)] <- FALSE;
    #if ( project == "u7" ){ res$reg_DM1_U7_C25[is.na(res$reg_DM1_U7_C25)] <- FALSE; }else if ( project == "dct" ){ res$reg_DM1_dCT_C25[is.na(res$reg_DM1_dCT_C25)] <- FALSE; }

    for (i in grep("EUC|val",colnames(res))) {
	    #print(dim(res))
	    res <- removeNA(res,i)
	  }

    return(res)

}

reannotate2 <- function (res,sampleInfo, seuil.l2fc= 1, seuil.padj=0.05, muscle)
{

  # Selection threshold
  mean.reads.min <- 30   # min read coverage for an exon bin in one condition
  
  #res=res.reannotated
  # Rename columns for clarity
  colnames(res)[1]            = "ensembl_ID";
  colnames(res)[2]            = "exon_bin";
  colnames(res)[8:10]         = paste("EUC",colnames(res)[8:10],sep="_");
  counts_data                 = as.matrix(res$countData)
  colnames(counts_data)       = paste0(sampleInfo$treatment,"_",sampleInfo$replicate)
  
  # Add a column "mean.reads" per condition
  for(i in (unique(sampleInfo$treatment))) {
    #for(j in muscle){
    res[,paste0(i,"_mean.reads")] <- rowMeans(counts_data[,grep(i,colnames(counts_data))])
  }
  
  # Add column "id" as "ensembl_ID:exon_bin"
  res[,"ExonID"] <- with(res, paste(ensembl_ID, exon_bin, sep = ":"))
  
  res[,"significant"] <- ((res$padj <= seuil.padj) &
                            !(is.na(res$log2FC_CTRL_vs_AAVGFP))  & 
                            !(is.na(res$log2FC_MBNLdecoy_vs_AAVGFP)))
  
  res[ which(is.na(res$significant)), "significant"] <- FALSE
  
  ## Add regulation Columns : 
  res[,"reg_Ctrl_AAVGFP"]  <- abs(res$log2FC_CTRL_vs_AAVGFP) >= seuil.l2fc & res$significant
  res[,"reg_Mbnld_AAVGFP"] <- abs(res$log2FC_MBNLdecoy_vs_AAVGFP) >= seuil.l2fc & res$significant
  res[,"regALL"]           <- (res$reg_Mbnld_AAVGFP | res$reg_Ctrl_AAVGFP)
  # Add correction %
  res[,"correction_prcnt"] <-  -(res$EUC_MBNLdecoy - res$EUC_AAVGFP) / 
                                (res$EUC_AAVGFP - res$EUC_CTRL) * 100;
  
  
  # Dealing with NA values (necessary for graphics)
  #res$reg_DM1_C25[is.na(res$reg_DM1_C25)] <- FALSE;
  #if ( project == "u7" ){ res$reg_DM1_U7_C25[is.na(res$reg_DM1_U7_C25)] <- FALSE; }else if ( project == "dct" ){ res$reg_DM1_dCT_C25[is.na(res$reg_DM1_dCT_C25)] <- FALSE; }
  
  for (i in grep("EUC|val",colnames(res))) {
    #print(dim(res))
    res <- removeNA(res,i)
  }
  
  return(res)
  
}

# Generate graphs for all bins and regulated bins (sample distance matrix, heatmap/heatplot and PCA)
#    Input  : type of graph to generate, new reannotated res object
#    Output : 2 png files (all bins and regulated bins)
plotting <- function(function.to.apply, res.obj,sampleInfo)
{

    # All bins/regulated bins
#    if ( project == "u7_all" ){
#       df.all = as.data.frame(res.obj)[, c(8:15)];
#        df.reg = subset(as.data.frame(res.obj), reg_HSA_NT_WT == TRUE)[, c(8:15)];
#    }else if (project == "dct_all") {
#        df.all = as.data.frame(res.obj)[, c(8:16)];
#        df.reg = subset(as.data.frame(res.obj), reg_HSA_NT_WT == TRUE)[, c(8:16)];
#    }

dfCounts=as.data.frame(res.obj$countData)
df.all=as.data.frame(res.obj)
df.reg=as.data.frame(res.obj[res.obj$significant==T,])
dfCounts.reg=dfCounts[res.obj$significant==T,]
#df.all=cbind(df.all[,grep("EUC",colnames(df.all))],dfCounts)
#df.reg=cbind(df.reg[,grep("EUC",colnames(df.reg))],dfCounts.reg)

df.all=dfCounts
df.sig=dfCounts.reg
df.regHSA=dfCounts[res.obj$reg_HSA_control==T,]
df.regPip6=dfCounts[res.obj$reg_pip6_control==T,]
df.regALL=dfCounts[res.obj$reg_HSA_control==T | res.obj$reg_pip6_control==T,]
df.regHSAPip=dfCounts[res.obj$reg_HSA_control==T & res.obj$reg_pip6_control==T,]


    # Graphs
    if ( function.to.apply == "distance" ){
	get.sample.matrix(df.all, "all",sampleInfo);
        get.sample.matrix(df.sig, "significant",sampleInfo);
        get.sample.matrix(df.regALL, "reg",sampleInfo);


    } else if ( function.to.apply == "heatmap") {
#        get.heatmap(df.all, "all", 50);
#        get.heatmap(df.all, "all", 100);
#        get.heatmap(df.all, "all", 200);
#        get.heatmap(df.all, "all", 500);
#        get.heatmap(df.all, "all", 1000);
#        get.heatmap(df.all, "all", 2000);
        get.heatmap(df.all, "all", 5000);
        get.heatmap(df.regALL, "reg");
        get.heatmap(df.regPip6, "regPip6");
        get.heatmap(df.regHSA, "regHSA");
        get.heatmap(df.sig, "significant");
        get.heatmap(df.regHSAPip, "regHSAPip");
    } else if ( function.to.apply == "pca" ) {

        get.pca(df.all, "all",sampleInfo);
        get.pca(df.sig, "significant",sampleInfo);
        get.pca(df.regHSA, "regHSA",sampleInfo);
        get.pca(df.regPip6, "regPip6",sampleInfo);
        get.pca(df.regALL, "regALL",sampleInfo);
        get.pca(df.regHSAPip, "regHSAPip",sampleInfo);


    }

}


# Generate png files showing sample distance matrix
#    Input  : df with EUC values, list type (all bins or reg bins)
#    Output : png file
get.sample.matrix <- function (df, type,sampleInfo)
{

    sampleDists <- dist(t(as.matrix(df)));
    sampleDistMatrix <- as.matrix(sampleDists);
    colors <- sampleInfo$Color

    png(paste0(output.dir, "sample_distance_matrix_", type, ".png"), units="px", width=3000, height=3000, res=450);
    pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, col=colors);
    dev.off();

}


# Generate png files showing heatmap and heatplot graphs
#    Input  : df with EUC values, list type (all bins or reg bins), nb bins to process (optional)
#    Output : png file
get.heatmap <- function (df, type, nb.bins = NULL)
{

    if( is.null(nb.bins) ){ nb.bins = dim(df)[1]; }
   # df=df[,-grep("EUC",colnames(df))]
    # Heatmap
   # topVar <- head(order(-rowVars(as.matrix(df))), nb.bins);
    #mat <- as.matrix(df)[topVar, ];
    #mat <- mat - rowMeans(mat);
    #png(paste0(output.dir, "heatmap_", type, "_", nb.bins, ".png"), units="px", width=3000, height=3000, res=450);
    #pheatmap(mat);
    #dev.off();

    # Heatplot (heatmap with heatplot algo)
    select <- order(rowMeans(as.matrix(df)), decreasing=TRUE)[1:nb.bins];
    mat.data <- as.matrix(df)[select, ];
    z <- zClust(mat.data);
    cols <- colorRampPalette(brewer.pal(10, "RdBu"))(256);
    png(paste0(output.dir, "heatplot_", type, "_", nb.bins, ".png"), units="px", width=3000, height=3000, res=450);
    #heatmap.2(z$data, trace='none', col=rev(cols), Rowv=z$Rowv, Colv=z$Colv, cexCol=0.7);
    heatmap.2(z$data, dendogram="row",trace='none', Rowv=z$Rowv, cexCol=0.7);
    dev.off();


}


## Necessary function for heatplot algo ussing heatmap
#zClust <- function(x, scale="row", zlim=c(-3,3), method="average") {

#    if (scale=="row") z <- t(scale(t(x)))
#    if (scale=="col") z <- scale(x)
#    z <- pmin(pmax(z, zlim[1]), zlim[2])
#    hcl_row <- hclust(distCor(t(z)), method=method)
#    ct<- cutree(hcl_row, k=4)
#   # hcl_col <- hclust(distCor(z), method=method)
#   # return(list(data=z, Rowv=as.dendrogram(hcl_row), Colv=as.dendrogram(hcl_col)))
#   return(list(data=z, Rowv=as.dendrogram(hcl_row),tree=ct))
#
#}


# Necessary function for heatplot algo using heatmap
distCor <- function(x) as.dist(1-cor(x));


# Generate PCA png file (prcomp function)
#    Input  : df with EUC values
#    Output : png file
get.pca <- function (df, type,sampleInfo)
{
#	df=df[,-grep("EUC",colnames(df))]
    # Colors
    colors <- as.character(sampleInfo$Color)
    # PCA
    data.pca    <- prcomp(t(as.matrix(df)));
    PC.variance <- round(data.pca$sdev^2/sum(data.pca$sdev^2)*100, 1);
    nb.bins=dim(df)[1]
    # 2 Principal Components
    png(paste0(output.dir, "pca2D_", type, "_",nb.bins,".png"), units="px", width=3000, height=3000, res=450);
    p<-plot(data.pca$x[,1], data.pca$x[,2], pch=19, cex=2, col=colors, xlab=paste0("PC1 = ", PC.variance[1], "%"), ylab=paste0("PC2 = ", PC.variance[2], "%"))+theme_gray(12)
;
    abline(h=0, col="grey", lty=2);
    abline(v=0, col="grey", lty=2);
    #text(data.pca$x[,1], data.pca$x[,2], labels=rownames(data.pca$x), col=colors);
    text(data.pca$x[,1], data.pca$x[,2], labels=sampleInfo$SampleName, col=colors,adj=c(0.2,1.2));
    dev.off();

    # 3 Principal Components
    #png(paste0(output.dir, "pca3D_", type, ".png"), units="px", width=3000, height=3000, res=450);
    plot3d(data.pca$x, xlab=paste0("PC1 = ", PC.variance[1], "%"), ylab=paste0("PC2 = ", PC.variance[2], "%"), zlab=paste0("PC3 = ", PC.variance[3], "%"), type = "s", col = colors);
    #text3d(data.pca$x, text=rownames(data.pca$x), adj=1.3);
    text3d(data.pca$x, text=sampleInfo$SampleName, adj=1.3);
    rgl.postscript(paste0(output.dir,"pca3D_",type,"_",nb.bins,".pdf"),"pdf")
    #dev.off();

}

norm<-function(data){

for(i in 1:dim(data)[2]){

data[,i]=data[,i]/sum(data[,i])*1000000

}
return(data)
}



zClust <- function(x, scale="row", zlim=c(-3,3), method="average") {
  # other Clustering methods=
  # c("ward.D", "single", "complete", "average", "mcquitty","median", "centroid", "ward.D2")
  if (scale=="row") {
    z <- t(scale(t(x)))
    z <- pmin(pmax(z, zlim[1]), zlim[2])
    hcl_row <- hclust(distCor(t(z)), method=method)
    ct      <- cutree(hcl_row, k=2)
    hcl_col <- hclust(distCor(z), method=method)
    return(list(data=z, hc=hcl_row, Rowv=as.dendrogram(hcl_row),tree=ct))
  }
  if (scale=="col") {
    z <- scale(x)
    z <- pmin(pmax(z, zlim[1]), zlim[2])
    hcl_col <- hclust(distCor(z), method=method)
    ct      <- cutree(hcl_col, k=2)
    return(list(data=z, hc=hcl_col, Colv=as.dendrogram(hcl_col),tree=ct ))
  }
}
