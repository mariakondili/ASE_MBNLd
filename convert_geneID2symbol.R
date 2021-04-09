#!/bin/env R
convert_ID2Name <- function(id_list){

	library(org.Hs.eg.db) # installed via biocLite. Cannot call library from an R-variable/object
	require(AnnotationDbi)
	gene_symbols <- mapIds(x=org.Hs.eg.db, keys=id_list,column="SYMBOL", keytype="ENSEMBL", multiVals = "first")
	return(gene_symbols)
}
