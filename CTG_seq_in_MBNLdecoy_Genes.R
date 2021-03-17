
ctg_data <- read.delim("../CTG_sequences_in_genes.txt",sep="\t", header=TRUE,as.is=TRUE)
## Obtained by Naira,given to Arnaud, since project of Pip6-therapy


length(which(ctg_data$external_gene_name %in% sigD$geneName))

## in UPREGULATED 
reg_decoy_names <- with(sigD, subset(geneName, abs(log2FoldChange_MBNLdVsCTRL) > 1))


ctg_genes.decoy <- subset(reg_decoy_names, reg_decoy_names  %in% ctg_data$external_gene_name )
length(ctg_genes.decoy)

# create_barplot_GOterms(ctg_genes.decoy, gene_format = "SYMBOL",main="Genes with CTG,(|log2fc|>1)")
library(biomaRt)
mart <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
ctg_genes_Descr <- getBM(attributes=c("external_gene_name","description"),
                          filters="external_gene_name", 
                          values=ctg_genes.decoy,
                          mart=mart)


