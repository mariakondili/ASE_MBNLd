#!/usr/bin/env R
library(dplyr)
library(ggplot2)



####---Calc MEAN PSI per Condition ----####

ctrl  <- grep("IncLevel1_*", colnames(all_ase))
dm1   <- grep("IncLevel2_*", colnames(all_ase))
mbnld <- grep("IncLevel3_*", colnames(all_ase))

mean_ctrl  <- rowMeans(all_ase[,ctrl])
mean_dm1   <- rowMeans(all_ase[, dm1])
mean_mbnld <- rowMeans(all_ase[,mbnld])

all_ase$Mean_Psi_CTRL <- mean_ctrl
all_ase$Mean_Psi_DM1  <- mean_dm1
all_ase$Mean_Psi_MBNLd <- mean_mbnld


##---- Correction-----##  =  ((MNBLd - DM1 ) / (CTRL - DM1) )*100

all_ase$correction_Pcnt <- with(all_ase,
                                 ((Mean_Psi_MBNLd - Mean_Psi_DM1)/ (Mean_Psi_CTRL -Mean_Psi_DM1)) *100)

 ## OR: 
#all_ase$correction_Pcnt_v2 <-  with(all_ase, (1- (dPSI_Ctrl_vs_MBNLd / dPSI_Ctrl_vs_DM1)) *100)
# is.na =  27499

all_ase$correction_Pcnt <- round(all_ase$correction_Pcnt,2)

outdir.ase <- paste0(post_dir,"Final_Results/")
write.table(all_ase[order(all_ase$FDR_Ctrl_vs_MBNLd, decreasing = F),],
            paste0(outdir.ase,"All_ASE_Ctrl+DM1+MBNLd_FDRsorted+Correction.tsv"),
            sep="\t",col.names = T,quote=F)


###--- PIE CHART OF CORRECTION LEVELS ----###
## create function to call with different cut-offs 

piechart_for_correction <- function(ase, psi.co, fdr.co, full_corr=80, plot_name="") { 
  
  filt.ase <- subset(ase, abs(dPSI_Ctrl_vs_DM1) >= psi.co & FDR_Ctrl_vs_DM1 <= fdr.co )


  full_correction <- ((filt.ase$correction_Pcnt > full_corr)  %>% which) %>% length
  partial_correction  <-  (((filt.ase$correction_Pcnt <= full_corr) & (filt.ase$correction_Pcnt >= 20)) %>% which) %>% length
  bad_correction <-   ((filt.ase$correction_Pcnt < 20) %>% which %>% length )


  corr_levels <- data.frame(category=c( paste0("full(>",full_corr ,"%)"), paste0("partial(20-", full_corr ,"%)"), "bad(<20%)"),
                            count=c(full_correction,partial_correction,bad_correction))

  ## calculate percentage of each count : 
  corr_levels$fraction = round(corr_levels$count / sum(corr_levels$count),2) * 100
  # Compute the cumulative percentages (top of each rectangle)
  corr_levels$ymax = cumsum(corr_levels$fraction)
  # Compute the bottom of each rectangle
  corr_levels$ymin = c(0, head(corr_levels$ymax, n=-1))
  # Compute label position
  # pos <- (splice_counts$ymax + splice_counts$ymin) / 2 # = 49.0 98.5 99.5
  # 98 & 99 are very close, I ll replace them with 97 and 10, left and right of the end of circle.
  # corr_levels$labelPosition <- c(5.0, 25.0,80.0)
  corr_levels$labelPosition <- (corr_levels$ymax + corr_levels$ymin) / 2

  # Compute a good label
  #corr_levels$label <- paste0(corr_levels$category,"\n", corr_levels$fraction, " %")
  # only fraction : 
  corr_levels$label <- paste0(corr_levels$fraction, " %")

  png(paste0(plot_name, ".png") )
  
    ggplot(corr_levels, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=category)) +
      geom_rect() +
      geom_label( x=3.5, aes(y=labelPosition, label=label), size=3) +
      coord_polar(theta="y") + 
      xlim(c(2, 4)) + 
      theme_void() +
      ggtitle( paste0("Splicing Correction MBNLΔ,\n|dPSI|>",dpsi.co,",FDR< ",fdr.co) )
  dev.off()      
 
 } 


## 1st Filter BY | dPSI | > 0.15 or |dPSI| > 0.20 :
my_psi_co= 0.15 
#my_psi_co=0.20
piechart_for_correction(all_ase, psi.co=my_psi_co, fdr.co=0.05,full_corr=80,
                        plot_name=paste0("PieChart_Correction_filt_ASE_CTRL-DM1-MBNLd_dPSI=",my_psi_co,".png") )
