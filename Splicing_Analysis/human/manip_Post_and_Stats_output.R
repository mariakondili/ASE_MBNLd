#!/usr/bin/env R

###-- Author  : Maria Kondili
###-- Date    : 11 March 2021

source("add_exonLength_and_Coord.R")

##> RUN function for Replicates.per.column, Exon-Length & up/down-stream Coordinates of splicing-event
## Add P-value, FDR, IncLevelDifference
## Read the Stat-results ,Add them to the Post-results-table with All Conditions.

add_stat_values <- function(modif_t, stats_pairs, s, tab_suffix) {

  #""" Extract statistics-values for each table (Ctrl-vs-DM1, Ctrl-vs-Decoy , Dm1-vs-Decoy)
  #""" and save them in list. Incorporate them all to the Post-results -modified table .

  stats_cols <- list()
  for (cond_path in stats_pairs) {
        stats_tab <- read.delim(paste0(cond_path,"/",s, tab_suffix))
        k = basename(cond_path) # = "Ctrl_vs_DM1"
        # Find corresponding events of Post-modified table with Stats-results-table, with ID.
        idx <- match(modif_t$ID, stats_tab$ID)
        ## modif_t$ID[52] == stats_tab$ID[1]
        names(idx) <- modif_t$ID+1 # ID is 0-based,but for index in R should be 1-based.
        stats_cols[[k]]<- as.data.frame(cbind(stats_tab$PValue[idx],
                                              stats_tab$FDR[idx],
                                              stats_tab$IncLevelDifference[idx]))

         colnames(stats_cols[[k]]) <- c(paste0("PValue_",k), paste0("FDR_",k),paste0("dPSI_",k))
         # Add each part of stats to modif.table with new columns
         modif_t[names(idx),colnames(stats_cols[[k]])] <- stats_cols[[k]][,1:3]
   }


  ## modif_t$ID[8] == stats_tab$ID[idx[8]]
  # all_stats_cols <- do.call(cbind, stats_cols)
  # colnames(all_stats_cols) <- sapply( 1:length(all_stats_cols) ,
  #                                     function(i) strsplit(colnames(all_stats_cols), "*\\.")[[i]][2])


  ##> add to Big-table of All columns :
  #modif_t[,colnames(all_stats_cols)] <- all_stats_cols

  return(modif_t)

}


## test : read a Post-results table :
#a3ss <- read.delim(paste0(post_dir,"A3SS.MATS.JCEC.txt"), header=T)

###> Merge columns to create a Large consensus table with ALL-EVENTS ####
merge_splice_events <- function(splice_events,stats_pairs, tab_suffix, post_dir){

  ##""" Take modified tables with Replicates-columns & Exon-coordinates
  ##""" and harmonize columns for all splicing-events
  listSpliceTabs <- list()

  for (s in splice_events) {

    cat("\n~~event parsed:", s,"\n")
    # Read events of Post-step output,
    # add new columns : IncLevel-1-2-3, exonCoord,Length ,splice_event
    modif_t <- add_exonLen_and_coord(splice_event=s,tab_suffix, in_dir=post_dir)
    #tab_suffix =".MATS.JCEC.txt"

    #modif_t : Is the Post-results table with All events,but diff.files per type.
    # Results after stats-step do not have ALL events.

    ## Append the new modif_t with new Pvalue,FDR for each condition-pair.
    modifstats_tab <- add_stat_values(modif_t, stats_pairs, s, tab_suffix)
    #modifstats_tab <- read.delim(paste0(indir, s, tab_suffix),header=T, as.is=T)

    ## arrange columns Of Interest:
    colsOI <- c("GeneID",	"geneSymbol","chr","strand",
                "IncFormLen", "SkipFormLen","ExonLength",
                grep ("IncLevel*", colnames(modifstats_tab)),
                "splice_event","alternative_exon_coordinates",
                grep("dPSI_",colnames(modifstats_tab), value=T),
                grep("PValue_",colnames(modifstats_tab),value=T),
                grep("FDR_",colnames(modifstats_tab),value=T)
               )

    listSpliceTabs[[s]] <-  modifstats_tab[,colsOI]
  }

  ##---SAVE ----##
  merged_splice_events <- do.call(rbind, listSpliceTabs)

  return(merged_splice_events)
}


## InclLevelDifference = Average(Inc_level_group1) - Average(Inc_level_group2)
## InclLevelDIfference  < 0  ==> group2 Larger incl.level
## InclLevelDIfference  > 0  ==> group1 Larger incl.level
