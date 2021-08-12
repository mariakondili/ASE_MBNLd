#!/usr/bin/env R

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
  corr_levels$labelPosition <- (corr_levels$ymax + corr_levels$ymin) / 2

  # Compute a good label
  #corr_levels$label <- paste0(corr_levels$category,"\n", corr_levels$fraction, " %")
  # only fraction :
  corr_levels$label <- paste0(corr_levels$fraction, " %")

  #png(paste0(plot_name, ".png") )

    ggplot(corr_levels, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=category)) +
      geom_rect() +
      geom_label( x=3.5, aes(y=labelPosition, label=label), size=3) +
      coord_polar(theta="y") +
      xlim(c(2, 4)) +
      theme_void() +
      ggtitle( paste0("Splicing Correction MBNLΔ,\n|dPSI|>",dpsi.co,",FDR< ",fdr.co) )
  #dev.off()

 }
