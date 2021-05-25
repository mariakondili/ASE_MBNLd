#!/usr/bin/env R 

add_psi_counts_per_repl <- function(event, rmats_dir, suffix) {
  #splice_events <- c("A3SS", "A5SS","MXE", "RI", "SE" )
  
  t = read.delim(paste0(rmats_dir, event, suffix),as.is=TRUE,header=TRUE)
  for (li in 1:nrow(t) ){
    ## Create new columns for each replicate's PSI-value
    t$IncLevel1_r1[li]  <- as.numeric(unlist(strsplit(t$IncLevel1[li],","))[1])
    t$IncLevel1_r2[li]  <- as.numeric(unlist(strsplit(t$IncLevel1[li],","))[2])
    t$IncLevel1_r3[li]  <- as.numeric(unlist(strsplit(t$IncLevel1[li],","))[3])
    
    t$IncLevel2_r1[li]  <- as.numeric(unlist(strsplit(t$IncLevel1[li],","))[4])
    t$IncLevel2_r2[li]  <- as.numeric(unlist(strsplit(t$IncLevel1[li],","))[5])
    t$IncLevel2_r3[li]  <- as.numeric(unlist(strsplit(t$IncLevel1[li],","))[6])
    
    t$IncLevel3_r1[li]  <- as.numeric(unlist(strsplit(t$IncLevel1[li],","))[7])
    t$IncLevel3_r2[li]  <- as.numeric(unlist(strsplit(t$IncLevel1[li],","))[8])
    t$IncLevel3_r3[li]  <- as.numeric(unlist(strsplit(t$IncLevel1[li],","))[9])
    
    t$splice_event[li]  <- event
  }
  # remove useless columns after new ones are created
  to_rmv <- c( str_extract(t, ".JC_SAMPLE_([1-2])"), str_extract(t,"IncLevel([1-2])$") )
  t <- dplyr::select(t, -all_of(to_rmv))

  # write.table(t, file=paste0(out_dir, "/",s,new_suffix), sep="\t",
  #                              quote=F,col.names=TRUE,row.names=F,na = "NA")
  return(t)
}


add_exonLen_and_coord <- function(splice_event, tab_suffix=".MATS.JCEC.txt", in_dir,out_dir) {
  #""" Receive the newly created tables with one column per replicate
  #    Create 2 new columns: ExonLength & upstream-alternative_exon-downstream coordinates """ 
  
  t <- add_psi_counts_per_repl(splice_event,in_dir,tab_suffix)
  
  library(stringr)
  # t = read.delim(paste0(in_dir, "/",splice_event, tab_suffix),as.is=TRUE,header=TRUE)
  
  for (i in 1:nrow(t)){ 
    if ( splice_event == "A3SS" || splice_event =="A5SS" ){
      # alt 5' 
      # when strand=="+" : t[,"flankingES"] > t[,"longExonEnd"]
      # when strand=="-" : t[,"flankingES"] < t[,"longExonEnd"]
      if ( t[i,"flankingES"] > t[i,"longExonEnd"] ) { # flanking on right of exon
        t$ExonLength[i] <-  t[i,"longExonEnd"] - t[i,"shortES"] #short length on left of exon
        alt_exon <- paste0(t[i,"shortEE"],"-",t[i,"longExonEnd"])
        downstream   <-  paste0(t[i,"flankingES"],"-",t[i,"flankingEE"])
        upstream  <- paste0(t[i,"shortES"], "-", t[i,"shortEE"]) # the short part of "exon-e"
      } 
      ## alt 3'
      # when t$strand[1] =="+" : t[,"flankingEE"] < t[,"longExonStart_0base"]
      # when t$strand[1] =="-" : t[,"flankingEE"] > t[,"longExonStart_0base"]
      
      if ( t[i,"flankingEE"] < t[i,"longExonStart_0base"] ) { # flanking on left of exon
        t$ExonLength[i] <- (t[i,"shortES"] - t[i,"longExonStart_0base"]) # short-length on right of exon
        alt_exon <- paste0(t[i,"longExonStart_0base"],"-",t[i,"shortES"])
        upstream <- paste0(t[i,"flankingES"],"-",t[i,"flankingEE"])
        downstream <- paste0(t[i,"shortES"], "-", t[i,"shortEE"])
      } 
      t$alternative_exon_coordinates[i]  <- paste(upstream, 
                                                  alt_exon, 
                                                  downstream,
                                                  sep="~")
    }
    
    if (splice_event == "RI" || splice_event =="SE" ){ 
      t$ExonLength[i] <- (t[i,which(str_detect(colnames(t),"xonEnd"))] - t[i,which(str_detect(colnames(t),"xonStart"))])
      upstream   <-  paste0(t[i,"upstreamES"],"-",t[i,"upstreamEE"])
      downstream <- paste0(t[i,"downstreamES"],"-",t[i,"downstreamEE"])
      alt_exon   <- paste0(t[i,which(str_detect(colnames(t),"xonStart"))],
                           "-",	
                           t[i,which(str_detect(colnames(t),"xonEnd"))])
      
      
      t$alternative_exon_coordinates[i]  <- paste(upstream, 
                                                  alt_exon, 
                                                  downstream,
                                                  sep="~")
    }
    
    if (splice_event == "MXE") {
      # X1stExonStart_0base	X1stExonEnd	
      # X2ndExonStart_0base	X2ndExonEnd	
      # upstreamES	upstreamEE	
      # downstreamES	downstreamEE
      exonLength_1 <- t[i,which(str_detect(colnames(t),"xonEnd"))[1]] - t[i,which(str_detect(colnames(t),"xonStart"))[1]]
      exonLength_2 <- t[i,which(str_detect(colnames(t),"xonEnd"))[2]] - t[i,which(str_detect(colnames(t),"xonStart"))[2]]
      t$ExonLength[i] <- paste(exonLength_1,exonLength_2,sep=";")
      upstream   <-  paste0(t[i,"upstreamES"],"-",t[i,"upstreamEE"])
      downstream <-  paste0(t[i,"downstreamES"],"-",t[i,"downstreamEE"])
      alt_exon   <-  paste( paste0(t[i,which(str_detect(colnames(t),"xonStart"))[1]],"-", t[i,which(str_detect(colnames(t),"xonEnd"))[1]]),
                            paste0(t[i,which(str_detect(colnames(t),"xonStart"))[2]],"-", t[i,which(str_detect(colnames(t),"xonEnd"))[2]]),
                            sep=";")
      
      t$alternative_exon_coordinates[i]  <- paste(upstream, 
                                                  alt_exon, 
                                                  downstream,
                                                  sep="~")
    }
  }
  
  # Write to new file : 
  # new_suffix=".JCEC_PerRepl_exonCoord.txt"
  # write.table(t, file=paste0(out_dir, "/", splice_event , new_suffix), sep="\t",
  #             quote=F,col.names=TRUE,row.names=F,na = "NA")
  # 
  
  t <- dplyr::select(t, -c("ID.1","PValue", "FDR","IncLevelDifference"))
  return(t)
}
