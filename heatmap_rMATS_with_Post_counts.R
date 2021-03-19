#!/usr/bin/env R 

# source("manip_rmats_output.R")
distCor <- function(x) as.dist(1-cor(x));

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

cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

plot_heatmap_signif_splice <- function(ase, signif.co, psi.co, title, plot_name){
    # ase.filt <- subset(ase, (FDR < signif.co &  abs(IncLevelDifference) > psi.co) )
    ase.filt <- subset(ase, (FDR_Ctrl_vs_MBNLd < signif.co &  abs(dPSI_Ctrl_vs_MBNLd) > psi.co) )
    
    ###  z.gene =  (E(g) - Mean(g.replicates) ) / M.A.D (g.replicates), mad=median.absolute.deviation
    x <- ase.filt[,8:16]
    s1 = strsplit(title, "_vs_")[[1]][1]
    s2 = strsplit(title, "_vs_")[[1]][2]
    s3 = strsplit(title, "_vs_")[[1]][3]
    colnames(x) <- c( paste(s1,"_", c("r1","r2","r3"),sep=""),
                      paste(s2,"_",c("r1","r2","r3"),sep=""),
                      paste(s3,"_",c("r1","r2","r3"),sep=""))
    
    
    x_norm <- t(apply(x, 1, cal_z_score))
    # x_norm <- x_norm[-which(apply(x_norm, 2, is.na)) ,]
    
    x_norm <- x_norm[-which(is.na(x_norm[,1])),]
    ## or zClust: 
    x_z <- zClust(x,scale="row", zlim=c(-2,2), method="average")
    
    #Rowv=as.dendrogram(hcl_row),tree=ct))
    library(gplots)
    library(RColorBrewer)
    
   # png(plot_name)
    print( heatmap.2(x_norm, 
              Rowv=F,
              Colv=F,
              distfun = dist, 
              hclustfun = hclust,
              dendrogram ='none',
              scale='row', na.rm=TRUE,
              trace='none',
              labCol=colnames(x),
              col=rev(colorRampPalette(brewer.pal(10, "RdBu"))(256)),
              margins=c(9,2),
              cexCol=1.3,
              srtCol=45, # declination of col.names 45degrees
              labRow="",
              offsetCol=0,
              key=TRUE,
              keysize = 1,
              density.info="density",
              key.title="Color range",
              key.xlab="norm.counts",
              key.ylab="Counts",
              key.par=list(mar=c(3,4,3,1)),
              notecex=2, 
              main=paste(nrow(x_norm),"signif. Diff/ly Spliced Exons\nin ",title))  )
    # + distfun = dist, hclustfun = hclust
    #print(hmp)
    #dev.off()

}


plot_Pheatmap <- function(){ 
  
  library(DESeq2)
  library(dendextend)
  library(pheatmap)
  
  ase.filt <- subset(ase, (FDR_Ctrl_vs_MBNLd < signif.co &  abs(dPSI_Ctrl_vs_MBNLd) > psi.co) )
  
  ###  z.gene =  (E(g) - Mean(g.replicates) ) / M.A.D (g.replicates), mad=median.absolute.deviation
  x <- ase.filt[,8:16]
  s1 = strsplit(title, "_vs_")[[1]][1]
  s2 = strsplit(title, "_vs_")[[1]][2]
  s3 = strsplit(title, "_vs_")[[1]][3]
  colnames(x) <- c( paste(s1,"_", c("r1","r2","r3"),sep=""),
                    paste(s2,"_",c("r1","r2","r3"),sep=""),
                    paste(s3,"_",c("r1","r2","r3"),sep=""))
  
  x_norm <- t(apply(x, 1, cal_z_score))
  # x_norm <- x_norm[-which(apply(x_norm, 2, is.na)) ,]
  
  x_norm <- x_norm[-which(is.na(x_norm[,1])),]
  my_hclust_gene <- hclust(dist(x_norm), method = "complete")
  pheatmap(x_norm,annotation_col = colnames(x_norm))
  
  
  }


## FILTER NON-SIGNIFICANT 

signif.co <- 0.05
psi.co <- 0.15

title="CTRL_vs_DM1_vs_MBNLdecoy"
plot_name="Heatmap_SignifAltSplicing_3conditions_rMATScounts.png"

plot_heatmap_signif_splice(all_ase, signif.co, psi.co, title, plot_name)


####--> Plot MBNLd-vs-Saline , after removing Common.Genes,spliced in Saline-vs-AAVGFP

# gfp.ase.genes <-  with (ase.gfp,
#                         geneSymbol[(abs(IncLevelDifference) > 0.15 & FDR < 0.05) %>% which])
# ase.saline_wo_gfp <- subset(ase.saline, ! (geneSymbol %in% gfp.ase.genes) ) # 2718 s.e Left !

plot_name = "Heatmap_SignifAltSplicing_Saline_vs_MBNLd_filtAav_rMATScounts.png"
plot_heatmap_signif_splice(ase.saline_wo_gfp, signif.co, psi.co , pair, plot_name)


pair = "Ctrl_vs_MBNLdecoy"
plot_name = "Heatmap_SignifAltSplicing_Ctrl_vs_MBNLdecoy_rMATScounts.png"
plot_heatmap_signif_splice(ase.decoy, signif.co, psi.co, pair, plot_name)

#---->> 2070 Signif.Diff/ly Spliced Exons 


####----- PIE-CHART of SIGNIF-SPLICED EVENTS and Others -------####
#### Show that differential splicing events in MBNLdecoy compaired to Saline are very few 

create_pieChart_Splicing <- function(ase,pair,signif.co,psi.co){
  ##""" Test with : ase = ase.saline; pair="Saline_vs_MBNLdecoy" ; signif.co =0.05 ; psi.co=0.15 """
  
  library(ggplot2)

  ## Count events Signif. splicing : 
  ## count InclDifference(dPSI) = mean(Saline) -  mean(MBNLdecoy)
  ## INCLUDED : PSI < -0.15 
  ## EXCLUDED : PSI > 0.15
  
  signif_exclu <- nrow(subset(ase, (FDR < signif.co  &  IncLevelDifference >= psi.co )))
  signif_inclu <- nrow(subset(ase, (FDR < signif.co &  IncLevelDifference <= -psi.co )))
  all_splice <-  nrow(ase)
  
  splice_counts <- data.frame(category=c("All", "Signif_Inclu", "Signif_Exclu"),
                              count=c(all_splice,signif_inclu,signif_exclu))
  
  ## calculate percentage of each count : 
  splice_counts$fraction = round(splice_counts$count / sum(splice_counts$count),2) * 100
  
  # Compute the cumulative percentages (top of each rectangle)
  splice_counts$ymax = cumsum(splice_counts$fraction)
  
  # Compute the bottom of each rectangle
  splice_counts$ymin = c(0, head(splice_counts$ymax, n=-1))
  
  # Compute label position
  # pos <- (splice_counts$ymax + splice_counts$ymin) / 2 # = 49.0 98.5 99.5
  # 98 & 99 are very close, I ll replace them with 97 and 10, left and right of the end of circle.
  splice_counts$labelPosition <- c(50.0, 97.0, 10.0)
  
  # Compute a good label
  splice_counts$label <- paste0(splice_counts$category,"\n", splice_counts$fraction, " %")
  
  # Make the plot
  dir.create("Plots/",showWarnings = F)
  png(paste0("Plots/Donut_Plot_Signif_SplicingEvents_",pair,".png"))
    ggplot(splice_counts, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=category)) +
      geom_rect() +
      geom_label( x=3.5, aes(y=labelPosition, label=label), size=3) +
      coord_polar(theta="y") + 
      xlim(c(2, 4)) + 
      theme_void()
 dev.off()
    
}

## call : 
create_pieChart_Splicing(ase.saline, pair="Saline_vs_MBNLdecoy",signif.co=0.05, psi.co=0.15)

create_pieChart_Splicing(ase.gfp, pair="AAVGFP_vs_Saline",signif.co=0.05, psi.co=0.15)



