make_piechart <- function(mbnl_dct_genes, clip_genes, main_title){

  suppressPackageStartupMessages(library(dplyr))

  # keep only unique genes of set A=mbnl_dct,remove the subset B="clip" genes
  my_DEgenes <- (mbnl_dct_genes  - clip_genes)

  groups <- data.frame(category=c("DE_MBNLdct", "CLIP_MBNL1"),
                       count=c(my_DEgenes ,clip_genes ))
  ## calculate percentage of each count :
  groups$fraction = round(groups$count / sum(groups$count),2) * 100

  # label only fraction
  groups$label <- paste0(groups$fraction, " %")

  suppressPackageStartupMessages(library(ggplot2))

  ggplot(groups, aes(x="",y=fraction, fill=category) ) +
    geom_bar(stat="identity",width=1) +
    coord_polar(theta="y",start=0) +
    theme_void() +
    scale_fill_manual(values=c("#E69F00", "#56B4E9")) +
    geom_text(aes(y = fraction/2 + c(0, cumsum(fraction)[-length(fraction)]),
                  label = label), size=5)+
    ggtitle(main_title )

}

##> call
## use sigD object from main script "DESeq2_analysis"
# make_piechart( length(sigD$geneName), nb_clip_genes,"MBNL1-CLIP-genes found in D.E.genes MBNLdct-vs-Saline")
