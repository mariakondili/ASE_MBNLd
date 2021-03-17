######################################################################################
#   Ré-écriture de fonctions du package DEXSeq(corrections bug + personalisation)    #
######################################################################################

#custom.restau.DEXSeqHTML(project, res.reannotated, res, genes=NULL, formulas=formulas, mart=mart, filter="ensembl_gene_id", attributes=c("ensembl_gene_id", "external_gene_name", "description", "gene_biotype"), color=NULL, color.samples=NULL, fitExpToVar=compa.cond);

#custom.restau.DEXSeqHTML(project, res.reannotated, res, genes=as.character(unique(y$geneID)), formulas=formulas, mart=mart, filter="ensembl_gene_id", attributes=c("ensembl_gene_id", "external_gene_name", "description", "gene_biotype"), color=NULL, color.samples=NULL, fitExpToVar=compa.cond);

# Custom function to generate html ouput (bug correction, less information, more clarity)
#    Input  : project, reannotated DEXSeqResults object, original DEXSeqResults object, etc 
#    Output : html file with 4 sheets per gene (counts, expression, splicing and results)
# Remark : necessity of original DEXSeqResults object for plotDEXSeq function call

library("yaml")
library("DEXSeq")
library("biomaRt")
library("hwriter")
library("ggplot2")
library("reshape2")
library("BiocParallel")
library("stringi")

################## Couleurs pour <=0.01, <=0.05, <=0.1, >0.1 

matCol <- function(padj,m2col){
    if(is.null(padj)) return(m2col[4]);
    if(!is.numeric(padj)) return(m2col[4]);
    if(padj<=0.01) return(m2col[1]);
    if(padj<=0.05) return(m2col[2]);
    if(padj<=0.1) return(m2col[3]);
    return(m2col[4]);

}


custom.restau.DEXSeqHTML <- function (projet, object, original.object, 
                                      expressionFile=NULL, genes = NULL, formulas = NULL, 
                                      path = "DEXSeqReport", file = "testForDEU.html", 
                                      fitExpToVar = "condition", FC=0.58, FDR = 0.1, 
                                      color = NULL, color.samples = NULL, mart = mart,
                                      species=species, filter = "", attributes = "", 
                                      extraCols = NULL, BPPARAM = MulticoreParam(1))
{

    ensPage=ifelse(species=="hg19","http://grch37.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=","https://www.ensembl.org/Mus_musculus/Gene/Summary?db=core;g=")
    stopifnot(is(object, "DEXSeqResults"));
  
    if(!(is.null(genes))){
	    original.object <- original.object[original.object$groupID %in% genes,]
	    object <- object[object$ensembl_ID %in% genes,]  
    }

    # -- Information of interest -- #
    
    # genomicData
    genomicData           <- as.data.frame(object$genomicData);
    colnames(genomicData) <- c("chr", "start", "end", "length", "strand");

    # DEU analysis
   
    results <- data.frame(object[, c(grep("ID|bin|geneName|description|biotype|log2FC|padj",colnames(object)),grep("EUC|correction_per|reg_|mean.reads",colnames(object)))], stringsAsFactors = TRUE);
    # Merge infos of interest
    results <- cbind(results, genomicData);
    #results <- results[, c(1,2,14,18,15,16,17,3,4,5,6,7,8,9,10,11,12,13,14)]; # genomicData at the begining
    
    # Data format
   
    if(length(grep("correction_per",colnames(results)))>0) results[, c("correction_per")] <- round(results[, c("correction_per")], 0);
    results[, grep("EUC",colnames(results))]<- round( results[,  grep("EUC",colnames(results))], 2 );
    results[, grep("log2FC",colnames(results))] <- round(results[, grep("log2FC",colnames(results))], 3);
    results[, grep("mean",colnames(results))]  <- floor(results[, grep("mean",colnames(results))]);
    results[, grep("padj",colnames(results))] <- round(results[, grep("padj",colnames(results))], 3);
    
    #results=results[order(abs(results[,grep("log2FC",colnames(results))[1]]),decreasing=T)  ,] 

    ## Colonnes plus petites => changer les noms

	#   if(projet=="pipao"){
	# 	  tmp=colnames(results)
	# 	  tmp=gsub("HSA_NT","HSA",tmp)
	# 	  tmp=gsub("control","ctrl",tmp)
	# 	  tmp=gsub("pip6a_treated_HSA","pip6a",tmp)
	# 	  colnames(results)=tmp
	#   }
    if(projet=="dct"){
        tmp=colnames(results)
    }
    #Pos de reg_TmT_control => MUST BE THE 1ST REG
    posReg=grep("reg",colnames(results))[1]

    # Colors per condition
    rownames(results) <- NULL;
    sampleData <- object@sampleData;
    
    if (is.null(color)) {
        numcond <- length(unique(sampleData[[fitExpToVar]]));
        #color <- rgb(colorRamp(c("#2B83BA", "#00CC00", "#FF0000"))(seq(0, 1, length.out = numcond)), maxColorValue = 255, alpha = 175);
        #color <- rgb(colorRamp(c("blue", "red", "green"))(seq(0, 1, length.out = numcond)), maxColorValue = 255, alpha = 175);
        sampleData$colors=gsub("bled","blue",sampleData$colors)
        color<-rgb(colorRamp(sampleData$colors)(seq(0, 1, length.out = numcond)), maxColorValue = 255, alpha = 175);
        names(color) <- sort(levels(sampleData[[fitExpToVar]]));
        # print(color)
    }
    

    ############### Adding expression data
    expTable <- NULL
    if(!(is.null(expressionFile))){
	      expData=read.delim(expressionFile);
     }
    else{expData=NULL;}

    # Genes to process (defaut : with TmT/Ctl regulated exon bins) 
    if (is.null(genes)) {
        #gns <- gsub("\\..*","",as.character(unique(results$ensembl_ID[which(results[,posReg] == TRUE)])));
	      gns <- unique(as.character(results$ensembl_ID[which(results[,posReg] == TRUE)]));
	      #gns=unique(as.character(original.object$groupID))
    } else {
        gns <- genes;
    }
    # print(head(gns))
    
    # Opening and writing of output directory + main sheet 
    if(!dir.exists(file.path(path, "files"))) dir.create(file.path(path, "files"), recursive = TRUE);
	    # print(path);
    p <- openPage(file.path(path, file));
    hwrite("DEXSeq differential exon usage test", p, heading = 1);
    hwrite("Experimental design", p, heading = 2);
 
    # Colors per sample
    cond <- as.matrix(as.data.frame(sampleData[, !colnames(sampleData) %in% "sizeFactor"]));
    rownames(cond) <- NULL;
    condcolor <- matrix(rep("white", nrow(cond) * ncol(cond)), nrow(cond));
    condcolor[, which(colnames(cond) %in% fitExpToVar)] <- color[as.character(sampleData[[fitExpToVar]])];
    if (!is.null(color.samples)) {
        condcolor[, 1] <- color.samples;
    }
    hwrite(cond, bgcolor = condcolor, p);
    hwrite(paste0("\n\nformulaDispersion = ", formulas[1]), p, heading = 3);
    hwrite(paste0("\nformula0 = ", formulas[2]), p, heading = 3);
    hwrite(paste0("\nformula1 = ", formulas[3]), p, heading = 3);
    hwrite("\nTestForDEU result table", p, heading = 2);
 
    # results sheet
    ptowrite <- file.path(path, "files/");
    global.legend <- hwrite(rbind(c("EUC", " : ", "Exon Usage Coefficient in condition (exon bin expression normalized by global gene expression in log2)"), c("log2FC", " : ", "Fold Change of 1st condition EUC against 2nd condition EUC (in log2)"), c("mean.reads", " : ", "mean number of reads mapped against exon bin in condition")), border=0);
    m2col <- colorRampPalette(c("#FF706B", "#FEE08B", "white"), space = "rgb")(4);
    

    colMat=rep(m2col[4],nrow(results))
    j <- 3
    for (i in c(0.1, 0.05, 0.01)) {
        colMat[which(results$padj <= i & results[,posReg] == TRUE)] <- m2col[j];
        j <- j - 1;
    }
    # matcol <- matrix(rep(m2col[4], (ncol(results)+1) * nrow(results)), nrow(results));
    # ## j=nomber of conditions!
    # j <- 3;
    # for (i in c(0.1, 0.05, 0.01)) {
    #     matcol[which(results$padj <= i & results[,posReg] == TRUE), ] <- m2col[j];
    #     j <- j - 1;
    # }
    pval.legend <- hwrite(c("Color code for selected differentially expressed genes with fold change >=2 AND padj :", "<= 0.01", "<= 0.05", "<= 0.1", "> 0.1"), bgcolor = c("#FFFFFF", m2col), border=0);
	  # print("hwrite")
    
    # Function to generate visualization sheets
    #    Input  : gene ensembl ID
    #    Output : 3 visu sheets (normalized counts, expression and splicing)
    makePagesForGene <- function(gene) {
      # print(gene)
  	  # write("prints to stdout", stdout())
  	  # write("prints to stderr", stderr())
      ## For overlapping genes
      nameforlinks <- sapply(strsplit(as.character(gene), "\\+"), "[[", 1);
  	  # print(c("name",nameforlinks))
      # Creation of the sheets
      back <- hwrite("back", link = file.path("..", file));
      otherlinks <- hwrite(c("counts", "expression", "splicing", "results"), link = c(paste(nameforlinks, "counts.html", sep=""), paste(nameforlinks, "expression.html", sep=""), paste(nameforlinks, "splicing.html", sep=""), paste(nameforlinks, "results.html", sep="")), table=FALSE);
          
      # Information of intrest
      loc <- as.character(results$ensembl_ID) %in% as.character(gene);
      subres <- results[loc, ];
      # print(subres)
      ############################## On ajoute l'annotation des positions dans la table des genes
  	  geneBins=subres[,c('start','end')]
  
  	  annotation=rep(NA,nrow(subres))
  
      ### exons constitutifs
  
  	  attributes=c("ensembl_gene_id","exon_chrom_start","exon_chrom_end","is_constitutive")
  	  geneExons=getBM(attributes=attributes,filter="ensembl_gene_id",value=gsub("\\..*","",gene),mart=mart)
  	  colnames(geneExons)=attributes
  	  rangesBins=IRanges(subres[,'start'],subres[,'end'])
  	  rangesExons=IRanges(geneExons[,'exon_chrom_start'],geneExons[,'exon_chrom_end'])
  	  ov=as.matrix(findOverlaps(rangesBins, rangesExons, type="any"))
  	  if(nrow(ov)>0){
  		    sumConst=rep(NA,nrow(subres))
  		    for(i in 1:nrow(ov)){
  		      if(length(grep(i,ov[,1]))>0) sumConst[i]=sum(geneExons$is_constitutive[ov[which(ov[,1]==i),2]])
  		    }
  	      annotation=ifelse(sumConst>0,'constitutive,',NA)
  	  }
      ##5 UTR
    	attributes=c("5_utr_start","5_utr_end")
    	geneUTRs=removeNA(getBM(attributes=attributes,filter="ensembl_gene_id",value=gsub("\\..*","",gene),mart=mart))
    	colnames(geneUTRs)=attributes
    	
    	rangesExons=IRanges(geneUTRs[,'5_utr_start'],geneUTRs[,'5_utr_end'])
    	ov=as.matrix(findOverlaps(rangesBins, rangesExons, type="any"))
    	annotation[unique(ov[,1])]=paste0(annotation[unique(ov[,1])],",5utr")
    
    
    	attributes=c("3_utr_start","3_utr_end")
    	geneUTRs=removeNA(getBM(attributes=attributes,filter="ensembl_gene_id",value=gsub("\\..*","",gene),mart=mart))
    	colnames(geneUTRs)=attributes
    	
    	rangesExons=IRanges(geneUTRs[,'3_utr_start'],geneUTRs[,'3_utr_end'])
    	ov=as.matrix(findOverlaps(rangesBins, rangesExons, type="any"))
    	annotation[unique(ov[,1])]=paste0(annotation[unique(ov[,1])],",3utr")
    
    	attributes=c("transcription_start_site")
    	tss=as.vector(t(getBM(attributes=attributes,filter="ensembl_gene_id",value=gsub("\\..*","",gene),mart=mart)))
    
    	posTSS=NULL
    	for(i in tss){
    	  posTSS=c(posTSS,which(subres[,'start']<=i & subres[,'end']>=i))
    	}
     	annotation[unique(posTSS)]=paste0(annotation[unique(posTSS)],",tss")
     	annotation=gsub("NA,","",annotation)
     	annotation=gsub("^,","",annotation)
     	annotation=gsub(",$","",annotation)
     	annotation=gsub(",,",",",annotation)
    	subres=as.data.frame(append(subres,list(annotation=annotation),after = 2))
    	### print(dim(subres))	
    
    
      ### il faut ajouter les coordonées mm9 convertir mm10->mm9
      # tmp=mm10Tomm9(subres[,c("chr","start","end")])
      # print(c("yes converting",gsub("\\..*","",subres$ensembl_ID),as.character(subres$ensembl_ID),tmp$chr,tmp$start.mm9,tmp$end.mm9))
    	# Ajouter le numero des exons!!!
    	subres=as.data.frame(append(subres,list(exon_nb=NA),after = 2))
    	subres$exon_nb=mapply(exonNumber,gsub("\\..*","",subres$ensembl_ID),as.character(subres$ensembl_ID),subres$chr,subres$start,subres$end)
    
      matcol=t(as.data.frame(lapply(colMat,rep,ncol(subres))))
    	submatcol <- matcol[loc, ];
      # print(subres)
      if(subres[1,'strand']=="-"){ subres=subres[nrow(subres):1,];submatcol=submatcol[nrow(submatcol):1,]} 
      # rownames(subres) <- NULL;
      # results sheet
    	# print(c("name",nameforlinks))
      genpage <- openPage(paste0(ptowrite, nameforlinks, "results.html"));
    	# print("open genpage 1")
      hwrite(c(back, otherlinks), table = TRUE, border = 0, genpage);
      hwrite("= Legend =", page=genpage);
      hwrite(global.legend, page=genpage, border=0);
      hwrite(pval.legend, page=genpage, border=0);
      # http://fasterdb.ens-lyon.fr/faster/main.pl?ensembl_initial=ENSG00000152601&stable_id_ensembl=ENSG00000152601&organism=human&bio_mol=cDNA_Only#
      ## FasterDB link
      # fasterdbLink="toto"
      # ensemblLink="toto"
      # print(gene)
      # print(species)
    
      ensemblLink= paste0("Ensembl link : <a href ='",ensPage,gsub("\\..*","",gene),"' target='_blank'>", gsub("\\..*","",gene)," </a>")
      fasterdbLink=paste0("Fasterdb link : <a href ='http://fasterdb.ens-lyon.fr/faster/main.pl?ensembl_initial=",gsub("\\..*","",gene),"&stable_id_ensembl=",gsub("\\..*","",gene),"&organism=",species,"&bio_mol=cDNA_Only#' target='_blank'>", gsub("\\..*","",gene)," </a>")
      # print(paste0("Fasterdb link : <a href ='http://fasterdb.ens-lyon.fr/faster/main.pl?ensembl_initial=",gene,"&stable_id_ensembl=",gene,"&organism=",species,"&bio_mol=cDNA_Only#' target='_blank'>", gene," </a>"))
      hwrite(fasterdbLink,heading=3,border=0,page=genpage)
      hwrite(ensemblLink,heading=3,border=0,page=genpage)
       
    	hwrite("Splicing information for all bins of gene", border=0,heading=3, genpage);
      # print(dim(subres))
      # print(dim(submatcol))
      hwrite(as.matrix(subres), bgcolor = submatcol, table.class = "sortable", row.names=F,row.bgcolor="#d8d6d6",row.style=list('text-align:center;'),style = "margin:16px; border:0px solid black; border-width:1px; width:200px; text-align:center", table = TRUE, page = genpage);
      close(genpage, splash = TRUE);
      # print("close genpage 1")
      # Image size regarding number of transcripts
      transcripts <- object$transcripts[loc];
      if (length(unlist(transcripts)) > 0) {
        trans <- Reduce(union, transcripts);
        h <- ifelse(length(trans) > 10, 7 + (length(trans) * 0.3), 7);
        if (sum(loc) > 30) {
            h <- h + (sum(loc) * 0.1);
        }
      }
	# print("TOTO")

	if(!is.null(expData)){	
		posGene=grep(sub("\\..*","",gene),expData[,1])
	  if(!is.null(posGene)){ 
			  expTable=expData[posGene,]
		}
	}

  # Generate images
  custom.makePlotPage(original.object, ptowrite=ptowrite, gene=gene, whichtag="expression", links=c(back, otherlinks), color=color, color.samples=color.samples, FDR=FDR, fitExpToVar=fitExpToVar, width=1400, height=1400, h=h,m2col=m2col,expTable=expTable);
  custom.makePlotPage(original.object, ptowrite=ptowrite, gene=gene, whichtag="counts", links=c(back, otherlinks), color=color, color.samples=color.samples, FDR=FDR, fitExpToVar=fitExpToVar, width=1400, height=1400, h=h,m2col=m2col,expTable=expTable);
  custom.makePlotPage(original.object, ptowrite=ptowrite, gene=gene, whichtag="splicing", links=c(back, otherlinks), color=color, color.samples=color.samples, FDR=FDR, fitExpToVar=fitExpToVar, width=1400, height=1400, h=h,m2col=m2col,expTable=expTable);
                
  return()
    
    }
    
    # print(head(gns))
    # Generate visualization sheets

    ### Don't make Gene pages!!!
    # bplapply(gns, makePagesForGene, BPPARAM = BPPARAM);
    # print(head(gns))
    # Gene table main sheet
    results <- results[as.character(results$ensembl_ID) %in% gns, ];
    # print("YOUPALAS")
    # print(head(results))
    #  print("YOUPALAF")
    
    splitCols <- split(seq_len(nrow(results)), results$ensembl_ID);
    genetable <- lapply(splitCols, function(x) {
        data.frame( chr = unique(results$chr[x]), start = min(results$start[x]), end = max(results$end[x]), total_bins = length(x), bin_changes = sum((results[,posReg][x] == TRUE), na.rm=TRUE))
    })	
    genetable <- do.call(rbind, genetable);
    genetable <- cbind(geneID = rownames(genetable), genetable);
    genetable=genetable[genes,]
    
    # print(genetable)
    #print(attributes(mart)$dataset)
    # Gene annotation info
    
    if (class(mart) == "Mart") {
        if (attributes(mart)$dataset != "") {
            forvalues <- strsplit(as.character(genetable$geneID), "\\+");
            names(forvalues) <- gsub("\\..*E","E",genetable$geneID);
            #	print(forvalues)
            if (length(filter) > 1) {
                warning("length(filter) > 2, only first element will be taken");
                filter <- filter[1];
            }
            #print(gsub("\\..*","",forvalues))
            extra <- getBM(attributes = c(filter, attributes), filters = filter, values = gsub("\\..*","",forvalues), mart = mart);
	          # print(extra)
            fromart <- lapply(gsub("\\..*","",genetable$geneID), function(x) {
                sep <- do.call(c, strsplit(as.character(x), "\\+"))
                extra[which(extra[, filter] %in% sep), ]
            })
            extra <- sapply(attributes, function(r) {
                      unlist(lapply(fromart, function(x) {
                              paste(x[, r], collapse = " ")
                        }))
                      })
           genetable <- cbind(geneID = genetable$geneID, extra, genetable[, 2:length(genetable)]);
        }
        else {
            warning("No dataset in biomart specified");
        }
    }
    else if (mart != "") {
        warning("Please provide a Mart class object for parameter mart");
    }
    if (!is.null(extraCols)) {
        genetable <- cbind(extraCols[match(genetable$geneID, rownames(extraCols)), ], genetable);
    }
    
    #################################
    # Main sheet
    genetable$geneID <- sapply(as.character(genetable$geneID), function(m) {
            w <- strsplit(m, "\\+")
            ns <- sapply(w, "[[", 1)
            hwrite(paste(unlist(w), collapse = " "), link = paste("files/", ns, "splicing.html", sep = ""))
    })
    
    rownames(genetable) <- NULL;
    hwrite(genetable, page = p, table = TRUE, table.class = "table-layout:fixed", row.bgcolor="#d8d6d6",row.style=list('text-align:center;'),style = "margin:16px; border:0px solid black; border-width:1px; width:20%");
    close(p, splash = TRUE);
    
}



# Custom function to generate gene visualization
#    Input  : original DEXSeqResults object 
#    Output : one of visu sheet (counts, expression or splicing) 
# Remark : necessity of original DEXSeqResults object for plotDEXSeq function call (impossible to rewrite...)

custom.makePlotPage <- function(ecs, ptowrite, gene, whichtag, links, color, color.samples, FDR, fitExpToVar, width, height, h,m2col,expTable)
{

    # Récupération des options
    allopts <- c("expression", "splicing", "counts");
    opts <- allopts %in% whichtag;
    
    # File creation
    onlytag <- allopts[max(which(opts))];
    pagename <- sapply(strsplit(as.character(gene), "\\+"), "[[", 1);   
    genpage <- openPage(paste(ptowrite, pagename, onlytag, ".html", sep=""));

    # Links
    hwrite(links, table=TRUE, border=0, genpage);
    
    # Caption
    if ( whichtag == "counts" ){
         hwrite("Counts visualization : Plotting of count values for each sample", border=0, genpage);
    }else if( whichtag == "expression" ){
        
    ##On rajoute une ligne pour les valeurs d'expression si on a le fichier d'expression
		if(!is.null(expTable)){
			rownames(expTable)<-NULL;
			hwrite("Expression data based on global gene expression with DESeq", border=0, heading=3,genpage);
			padj=expTable[1,grep("padj",colnames(expTable))[1]]			
			if(!is.na(padj)){ mt=matCol(padj,m2col)
		
			hwrite(expTable, bgcolor = mt, table.class = "sortable", row.bgcolor="#d8d6d6",row.style=list('text-align:center;'),style = "margin:16px; border:1px solid black; border-width:1px; width:200px", table = TRUE, page = genpage);
			}
			hwrite("<hr>",page=genpage)
			hwrite("Expression visualization based on splicing analysis with DEXSeq : Plotting of fitted expression estimates for each condition", border=0, heading=3,page=genpage);
		}
    	
    }else{
         hwrite("Splicing visualization : Plotting of Exon bin Usage values (exon bin expression normalized by global gene expression) = Best visualization for AS events", border=0, genpage);
    }
    
    # print(gene)
    # print(dim(ecs[ecs[,1]==gene,])[1]); 
    if(dim(ecs[ecs[,1]==gene,])[1]<30){ wd=14;ps=8}
    #  else{ wd=16*2;ps=24}
    else {wd=24;ps=18}
    # print(wd)
    # Image
    # svg(paste(ptowrite, pagename, onlytag, ".svg", sep=""), height=h, width=16, pointsize=14);
    
    svg(paste(ptowrite, pagename, onlytag, ".svg", sep=""), height=h, width=wd, pointsize=ps);
      # print(gene)
      # plotDEXSeq(ecs, geneID=gene, FDR=FDR, lwd=2, expression=opts[1], splicing=opts[2], norCounts=opts[3], displayTranscripts=TRUE, fitExpToVar=fitExpToVar, names=TRUE, legend=TRUE, color=color, color.samples=color.samples, cex.axis=1.5); 
      plotDEXSeq(ecs, geneID=gene, FDR=FDR, expression=opts[1], splicing=opts[2], norCounts=opts[3], displayTranscripts=TRUE, fitExpToVar=fitExpToVar, names=TRUE, legend=TRUE, color=color, color.samples=color.samples); 
    dev.off();
   
    hwrite(hmakeTag("iframe", src=paste(pagename, onlytag, ".svg", sep=""), width=width, height=1000, border=0), p=genpage);
    # if(dim(ecs[ecs[,1]==gene,])[1]<30){ 
    #   hwrite(hmakeTag("iframe", src=paste(pagename, onlytag, ".svg", sep=""), width=width, height=1000, border=0), p=genpage);
    #}
    #else{
    #   hwrite(hmakeTag("iframe", src=paste(pagename, onlytag, ".svg", sep=""), width=width, height=height*1.5, border=0), p=genpage);

    #}

    close(genpage, splash=TRUE);

}

