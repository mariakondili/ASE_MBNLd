#!/usr/bin/env R
## Author: Maria Kondili
## Subject : See degerulated pathways and processes from Diff.expressed genes in 3 conditions : saline-AAVGFP- MBNL-decoy

## ClusterProfiler :
# doc : https://bioconductor.org/packages/release/bioc/vignettes/clusterProfiler/inst/doc/clusterProfiler.html
# and : http://yulab-smu.top/clusterProfiler-book/index.html
suppressPackageStartupMessages(library(DOSE))
suppressPackageStartupMessages(library(enrichplot))
suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(org.Mm.eg.db))
suppressPackageStartupMessages(library(AnnotationDbi))

create_dotplot_GOterms <- function(genenames,gene_format = "SYMBOL", my_ontol="ALL",main=""){
  EntrezID <- mapIds(x=org.Mm.eg.db, keys=genenames,
                     column="ENTREZID", keytype=gene_format,
                     multiVals = "first")
  ##to see accepted keytypes :  keytypes(org.Mm.eg.db)
  ## column = what I want to obtain
  ## keytype = what I give in keys

  go.enrichm <- enrichGO(EntrezID, ont=my_ontol,
                            OrgDb=org.Mm.eg.db,
                            keyType = "ENTREZID",
                            pvalueCutoff=0.1,
                            qvalueCutoff=0.1,
                            pAdjustMethod = "BH",
                            readable=TRUE)

  # go.enrichm.mf <- enrichGO(sig.data.d$Entrez_ID, ont="MF",OrgDb=org.Mm.eg.db,keyType = "ENTREZID",
  #                           pvalueCutoff=0.05, qvalueCutoff=0.1, pAdjustMethod = "BH", readable=TRUE)
  #
  if (my_ontol=="BP") {msg = "Biol.Processes" }
  if (my_ontol=="MF") {msg = "Molec.Functions"}
  if (my_ontol=="CC") {msg = "Cellular Components" }
  if (my_ontol=="ALL") {msg = "All GO-terms"}

  # barplot(go.enrichm,showCategory = 20,
  #         x="GeneRatio",
  #         color="p.adjust", font.size=9,
  #         title=paste(msg," of ", main))

  dotplot(go.enrichm,showCategory = 30,
          x="GeneRatio",color="p.adjust",
          font.size=9,title=paste(msg," of ", main))

}


####------- signif & Highly-expressed Genes MBNLdecoy vs CTRL -----------####
l2fc.co  <- 1
# sigD : created in DESeq_analysis script, contains geneIDs and stat.meausures of DESeq2 
high_sig_genes <- subset(sigD$GeneID, abs(sigD$log2FoldChange_MBNLdVsCTRL)>= l2fc.co ) #874 genes

create_dotplot_GOterms(high_sig_genes,gene_format = "ENSEMBL",my_ontol="ALL",
                       main="Highly Regul.Genes in\n MBNL.Decoy/CTRL,(|log2FC|>1)")

###----- Investigate  enrichGO results to manipulate terms to plot ------###
myEntrez_genes <- mapIds(x=org.Mm.eg.db,keys=high_sig_genes,column="ENTREZID",keytype="ENSEMBL",multiVals = "first")

go_terms_sigHighGenes <- enrichGO(myEntrez_genes, ont="ALL",
                       OrgDb=org.Mm.eg.db,
                       keyType = "ENTREZID",
                       pvalueCutoff=0.1,
                       qvalueCutoff=0.1,
                       pAdjustMethod = "BH",
                       readable=TRUE)

#go_res <- go_terms_sigHighGenes@result ## $GeneRatio , $Count , $p.adjust , $pvalue
#goTerms_HighCounts <- subset(go_res$Description, go_res$Count > 5)

#order All from High-->Low Count
#go_res.ord <- go_res[order(go_res$Count,decreasing = T),]

# Update DOSE object with Ordered by Count Terms
go_terms_sigHighGenes@result <- go_terms_sigHighGenes@result[order(go_terms_sigHighGenes@result$Count,decreasing = T),]

## Write Down the Ordered GO-Terms Description :
##( Hand-pick which to plot )

write.table (as.data.frame(cbind(go_terms_sigHighGenes@result$Description,
                    go_terms_sigHighGenes@result$Count,
                    go_terms_sigHighGenes@result$p.adjust )),
                    "Ordered_by_Count_GO_terms.tsv",
                    row.names = F, col.names = c("Description","Count","p_adj"),
                    sep="\t", quote=F)

dotplot(go_terms_sigHighGenes,showCategory = 30,
        x="GeneRatio",color="p.adjust",
        font.size=9)


#--> BP.terms :
# go.enrichm.bp@result$Description
# [1]"negative regulation of locomotion"                  "negative regulation of cellular component movement"
# [3]"brown fat cell differentiation"                     "negative regulation of cell migration"
# [5] "negative regulation of cell motility"               "axonogenesis"
# [7] "renal system development"                           "urogenital system development"
# [9] "fat cell differentiation"                           "regulation of developmental growth"


##------------ Genes signif & UPregul .MBNLdecoy ----------##

sigUp_decoy <- sigD$geneName[which( sigD$log2FoldChange_MBNLdVsCTRL > 0)]
create_barplot_GOterms(sigUp_decoy, gene_format = "SYMBOL",my_ontol="BP",
                       main="signif.Upreg.Genes in\n Decoy-vs-Ctrl,(log2fc>0)")


###--------------- HITS-CLIP Genes found in UPregulated MBNLdecoy ---------------##

len_upreg_decoy <- with(sigD, subset(geneName, log2FoldChange_MBNLdVsCTRL > 0)) %>% length
upreg_decoy_names <- with(sigD, subset(geneName, log2FoldChange_MBNLdVsCTRL > 0))

##>> signif + Upreg @ CLIP-MBNL
clip_genes_in_UPdecoy <- (subset(hits_clip$gene_symbol, hits_clip$gene_symbol %in% upreg_decoy_names)
                        %>% unique)

clip_genes_in_UPdecoy %>% length
##> 17 genes when l2fc >1
##> 197  genes when l2fc >0

create_barplot_GOterms(clip_genes_in_UPdecoy,gene_format = "SYMBOL",my_ontol="BP",
                       main="MBNLdecoy-Upreg.genes \nfound in CLIP-MBNL.KO")


#####
###-------- HITS-CLIP Genes found in DOWN-regulated MBNLdecoy  ------------------##
downreg_decoy_names <- with(sigD, subset(geneName, log2FoldChange_MBNLdVsCTRL < 0))

len_downreg_decoy <- downreg_decoy_names %>% length

##>> signif + Downreg @ CLIP-MBNL
clip_genes_in_down_decoy <- (subset(hits_clip$gene_symbol,
                       hits_clip$gene_symbol %in% downreg_decoy_names)
                        %>% unique)

clip_genes_in_down_decoy %>% length
##> 182  genes when l2fc >0

create_barplot_GOterms(clip_genes_in_down_decoy,gene_format = "SYMBOL",my_ontol="BP",
                       main="MBNLdecoy-Downreg.genes \nfound in CLIP-MBNL.KO")



##-----------DOT PLOT ---------------##
#png("ClusterProfiler_Graphs/Dotplot_40BP_decoy_vs_AAVGFP_log2FC=2.png",width = 1700 )
dotplot(go.enrichm.bp , x="GeneRatio",
        color = "p.adjust", showCategory=40,
        title = "Most affected Biol.Processes)",
        font.size = 10)
#dev.off()
##--> or : go.enrichm.mf



## Biolog.Processes related to Muscle :
## found by go.enrichm.bp$Description
# "angiogenesis involved in wound healing"
# "positive regulation of smooth muscle cell proliferation"
# "regulation of muscle system process"
# "vascular endothelial growth factor production"
# "regulation of muscle contraction"
# "regulation of smooth muscle contraction
# "regulation of myelination"
# "regulation of angiogenesis"
# "regulation of muscle organ development" ?


muscle_pos_terms <- which( go.enrichm.bp$Description %in%
         c("positive regulation of smooth muscle cell proliferation",
           "vascular endothelial growth factor production",
           "regulation of muscle contraction",
           "regulation of angiogenesis","angiogenesis",
           "cell adhesion",  "cell matrix adhesion") )


## Find "geneID":
go.enrichm.bp[muscle_pos_terms,c("Description" ,"p.adjust","geneID")]
#go.enrichm.bp@gene2Symbol[muscle_pos_terms]

enrichplot::dotplot(go.enrichm.bp.common , x="GeneRatio",
                    color = "p.adjust", showCategory=40,
                    title = "Most affected Biol.Processes (|log2FC|>1)",
                    font.size = 10)



##--------- CNET PLOT -----------##
## create Data.frame with EntrezID <> FC,
# fc_df <- rbind(sig.data$gene, sig.data$log2FoldChange_MBNLdecoyVsCTRL)

## create Data.frame with gene.symbol <> FC.value
match.gn <- match(go.enrichm.bp@gene2Symbol, sig.data$geneName)
# length(sig.data$log2FoldChange_MBNLdecoyVsCTRL[match.gn]) ==
# length(go.enrichm.bp@gene2Symbol)

fc_genes_df <- data.frame(go.enrichm.bp@gene2Symbol[-which(is.na(go.enrichm.bp@gene2Symbol))],
                        "foldChange"=sig.data$log2FoldChange_MBNLdecoyVsCTRL[match.gn][-which(is.na(go.enrichm.bp@gene2Symbol))])

fc = sig.data$log2FoldChange_MBNLdecoyVsCTRL[match.gn]

library(ggplot2)
library(enrichplot)
library(DOSE)
library(ggnewscale)


cnetplot(go.enrichm.mf,circular=T,
         colorEdge=TRUE,showCategory = 10,
         foldChange = fc_genes ) +
         ggtitle("Gene-Concept Network\nGO Molecular Functions\nD.E genes |log2FC|>1 ") +
         theme(text=element_text(size=9),aspect.ratio = 0.40)
