#!/usr/bin/env R

###-- Subject : -Merge rMATS output tables in one ,add a column for the type of splicing event.
###             -Create Heatmaps & Pie Chart of splic.events
###-- Author  : Maria Kondili, Centre de recherche de Myologie
###-- Date    : 11 March 2021

work_dir="/shared/projects/mbnl_dct/human_cellline/rMATS_SplicingCounts/tasks_Prep_Post_Stat/"


splice_events <- c("A3SS","A5SS","MXE","RI","SE")

source("add_exonLength_and_Coord.R")


##--->RUN function for Replicates.per.column, Exon-Length & up/down-stream Coordinates of splicing-event 

splice_events <- c("A3SS", "A5SS","MXE", "RI", "SE" )
## output in same folder but with suffix "_exonCoord"

post_dir <- paste0(work_dir,"Post_Ctrl_vs_DM1_vs_MBNLdecoy/" )

stats_pairs <-  c(paste0(stats_dir,"Ctrl_vs_DM1"), 
                  paste0(stats_dir,"Ctrl_vs_MBNLd"),
                  paste0(stats_dir,"DM1_vs_MBNLd") )


## Create new folder to put modified tables
outdir_stats  <- paste0(work_dir,"modif_Stats_tables/")
dir.create(outdir_stats,showWarnings = F)


####----> Add P-value, FDR, IncLevelDifference 
####----> Read the Stat-results ,Add them to the Post-results-table with All Conditions.

add_stat_values <- function(modif_t, cond_pairs, s, tab_suffix) { 
  
  #""" Extract statistics-values for each table (Ctrl-vs-DM1, Ctrl-vs-Decoy , Dm1-vs-Decoy)
  #""" and save them in list. Incorporate them all to the Post-results -modified table .
  
  stats_cols<- list()
  for (cond_path in cond_pairs) { 
        stats_tab <- read.delim(paste0(cond_path,"/",s, tab_suffix))
        k = strsplit(cond_path,"/")[[1]][9] # = "Ctrl_vs_DM1"
        idx <- match(modif_t$ID, stats_tab$ID)
        names(idx) <- modif_t$ID
        stats_cols[[k]] <- as.data.frame(cbind(stats_tab$PValue[idx],
                            stats_tab$FDR[idx],
                            stats_tab$IncLevelDifference[idx]))
         colnames(stats_cols[[k]]) <- c(paste0("PValue_",k), paste0("FDR_",k),paste0("dPSI_",k))
   }
 
                  
  ## modif_t$ID[8] == stats_tab$ID[idx[8]]
  all_stats_cols <- do.call(cbind, stats_cols)
  colnames(all_stats_cols) <- sapply( 1:length(all_stats_cols) ,
                                      function(i) strsplit(colnames(all_stats_cols), "*\\.")[[i]][2])
  
  ##> add to Big-table of All columns :
  modif_t[,colnames(all_stats_cols)] <- all_stats_cols
 
  return(modif_t)
  
}


## test : read a Post-results table : 
#a3ss <- read.delim(paste0(post_dir,"A3SS.MATS.JCEC.txt"), header=T)
tab_suffix = ".MATS.JCEC.txt"


##--> PARALLELISE IN CORES ???
# library(foreach)
# library(doParallel)
# 
# #setup parallel backend to use many processors
# cores=detectCores()
# cl <- makeCluster(cores[1]-1) #not to overload your computer
# registerDoParallel(cl)
# 
# finalMatrix <- foreach(i=1:150000, .combine=cbind) %dopar% {
#   tempMatrix = functionThatDoesSomething() #calling a function
#   #do other things if you want
#   
#   tempMatrix #Equivalent to finalMatrix = cbind(finalMatrix, tempMatrix)
# }
# #stop cluster
# stopCluster(cl)

for (s in splice_events) { 
  
    cat("\n~~event parsed:", s,"\n")
    # add new columns : IncLevel-1-2-3, exonCoord,Length ,splice_event
    modif_t <- add_exonLen_and_coord(splice_event=s,
                                     tab_suffix=".MATS.JCEC.txt",
                                     in_dir=post_dir,
                                     outdir_stats)
    #modif_t : Is the Post-results tables with All-events
    # stats_tab doesn't have ALL events.
    
    #for (k in stats_pairs){
    ## Append the new modif_t with new Pvalue,fdr  to next loop to add the next pvalue,fdr aside.
    final_tab <- add_stat_values(modif_t, stats_pairs, s, tab_suffix)
      
    new_suffix=".JCEC_3conds_exonCoord.txt"
    write.table(final_tab, paste0(outdir_stats, s, new_suffix),
                sep="\t", col.names = T, quote=F )
}



####-----> Merge columns to create a Large consensus table with ALL-EVENTS ####

merge_splice_events <- function(splice_events, tab_suffix, indir, filename){

  ##""" Take modified tables with Replicates-columns & Exon-coordinates
  ##""" and harmonize columns for all splicing-events
  listSpliceTabs <- list()
  
  for (s in splice_events) {
    modifstats_tab <- read.delim(paste0(indir, s, tab_suffix),header=T, as.is=T)
    colsOI <- c("GeneID",	"geneSymbol","chr","strand",
                "IncFormLen", "SkipFormLen","ExonLength",
                "IncLevel1_r1","IncLevel1_r2","IncLevel1_r3",	
                "IncLevel2_r1","IncLevel2_r2","IncLevel2_r3",	
                "IncLevel3_r1","IncLevel3_r2","IncLevel3_r3",
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


##> call function 
new_suffix=".JCEC_3conds_exonCoord.txt"
outdir.ase <- paste0(post_dir,"Final_Results/")
dir.create(outdir.ase,showWarnings = F)
all_ase <- merge_splice_events(splice_events,
                               new_suffix,
                               indir=outdir_stats, 
                               outdir=outdir.ase,
                               filename="Ctrl+DM1+MBNLd")

lev1_r1 <- which(is.na(all_ase$IncLevel1_r1))
lev1_r2 <- which(is.na(all_ase$IncLevel1_r2))
lev1_r3 <- which(is.na(all_ase$IncLevel1_r3))

lev2_r2 <- which(is.na(all_ase$IncLevel2_r2))
lev2_r3 <- which(is.na(all_ase$IncLevel2_r3))
lev2_r1 <- which(is.na(all_ase$IncLevel2_r1))

lev3_r1 <- which(is.na(all_ase$IncLevel3_r1))
lev3_r2 <- which(is.na(all_ase$IncLevel3_r2))
lev3_r3 <- which(is.na(all_ase$IncLevel3_r3))


nas <- intersect( c(c(lev1_r1,lev1_r2, lev1_r3), 
                  c(lev2_r3, lev2_r2,lev2_r1)), 
                  c(lev3_r1,lev3_r2,lev3_r3))
all_ase <- all_ase[-nas,]

write.table(all_ase[order(all_ase$FDR_Ctrl_vs_MBNLd, decreasing = F),],
            paste0(outdir.ase,"merged_modif_ASE_",filename,"_FDRsorted.tsv"),
            sep="\t",col.names = T,quote=F)


