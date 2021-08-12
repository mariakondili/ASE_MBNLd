
###-- Subject : -> Merge rMATS output tables in one ,add a column for the type of splicing event.
###             -> Create Heatmaps & Pie Chart of splic.events
###-- Author  : Maria Kondili
###-- Date    : 11 March 2021

work_dir="/projects/mbnl_dct/human_cellline/rMATS_SplicingCounts/tasks_Prep_Post_Stat/"
post_dir <- paste0(work_dir,"Post_Ctrl_vs_DM1_vs_MBNLdecoy/")
stats_pairs <-  c(paste0(post_dir,"Ctrl_vs_DM1"),
                  paste0(post_dir,"Ctrl_vs_MBNLd") )

source("manip_Post_and_Stats_output.R")
##> call function
# new_suffix=".JCEC_3conds_exonCoord.txt" #--> didn't write intermed tables with stats, use R-object modif_t instead .
tab_suffix = ".MATS.JCEC.txt"
splice_events <- c("A3SS","A5SS","MXE","RI","SE")

outdir.ase <- paste0(post_dir,"Final_Results/")
dir.create(outdir.ase,showWarnings = F)

#> call function for all steps:
#1 / add exonCoord column
#2 / add p-value,FDR, dPSI per pair of conditions
# 3/ Merge tables for all splicing-types (A3SS,A5SS,SE,MXE,RI)
all_ase <- merge_splice_events(splice_events,
                                 stats_pairs,
                                 tab_suffix,
                                 post_dir)
##> 112186 events !

## $ID = idx -1 => to find original events in Stats.table should look for ID= pos-1
## all_ase["SE.4666",] == SE.MATS.JCEC.["4665",]

#> CTRL (clone48 = 3 replicates)
lev1_r1 <- which(is.na(all_ase$IncLevel1_r1))
lev1_r2 <- which(is.na(all_ase$IncLevel1_r2))
lev1_r3 <- which(is.na(all_ase$IncLevel1_r3))
lev1_r4 <- which(is.na(all_ase$IncLevel1_r4))
lev1_r5 <- which(is.na(all_ase$IncLevel1_r5))
lev1_r6 <- which(is.na(all_ase$IncLevel1_r6))

#> DM1 samples
lev2_r1 <- which(is.na(all_ase$IncLevel2_r1))
lev2_r2 <- which(is.na(all_ase$IncLevel2_r2))
lev2_r3 <- which(is.na(all_ase$IncLevel2_r3))

#> MBNLdct samples
lev3_r1 <- which(is.na(all_ase$IncLevel3_r1))
lev3_r2 <- which(is.na(all_ase$IncLevel3_r2))
lev3_r3 <- which(is.na(all_ase$IncLevel3_r3))


nas <- intersect( c(c(lev1_r1,lev1_r2, lev1_r3,
                     lev1_r4,lev1_r5,lev1_r6),
                  c(lev2_r1, lev2_r2,lev2_r3)),
                  c(lev3_r1,lev3_r2,lev3_r3)  )

all_ase.clean <- all_ase[-nas,] #91 176

filename="Ctrl_DM1_MBNLd"
write.table(all_ase[order(all_ase$FDR_CTRL_vs_MBNLd, decreasing = F),],
            paste0(outdir.ase,"merged_modif_ASE_",filename,"_FDRsorted.tsv"),
            sep="\t",col.names = T,quote=F)

source("heatmap_rMATS_with_Post_counts.R")
## FILTER NON-SIGNIFICANT

signif.co <- 0.05
psi.co <- 0.20

title="CTRL_vs_DM1_vs_MBNLdecoy"
plot_name="Pheatmap_SignifAltSplicing_dPSI=0.20_Ctrl-DM1-MBNLdecoy.png"

plot_Pheatmap(ase, signif.co, psi.co, title, conditions=c("CTRL","DM1","MBNLd"), plot_name)
