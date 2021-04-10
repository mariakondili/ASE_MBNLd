piechart_for_correction <- function(D, l2fc.co, padj.co, plot_name=""){

  suppressPackageStartupMessages(library(dplyr))
  D <- subset(ase, abs(dPSI_CTRL_vs_DM1) > l2fc.co & FDR_CTRL_vs_DM1 < padj.co )

  full <- ((D$Correction_Pcnt > 80)  %>% which) %>% length #874
  partial  <-  (((D$Correction_Pcnt <= 80) & (D$Correction_Pcnt >= 20)) %>% which) %>% length
  bad <-   ((D$Correction_Pcnt < 20) %>% which %>% length )

  corr_levels <- data.frame(category=c("full(>80%)", "partial(20-80%)", "bad(<20%)"),
                            count=c(full,partial,bad))

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

  #png("PieChart_CorrectionPercent_filtASE.png")
  ggplot(corr_levels, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=category)) +
    geom_rect() +
    geom_label( x=3.5, aes(y=labelPosition, label=label), size=3) +
    coord_polar(theta="y") +
    xlim(c(2, 4)) +
    theme_void() +
    ggtitle( paste0("Expr.Correction MBNLΔ,\n|log2fc|>",l2fc.co, ", p.adj < " ,padj.co) )
  #dev.off()

}


piechart_for_correction(D.sig, l2fc.co=1, padj.co=0.05, plot_name="")
