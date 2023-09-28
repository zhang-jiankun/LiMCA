#!/usr/env/bin Rscript

## See https://bioconductor.org/packages/release/bioc/vignettes/nullranges/inst/doc/matchRanges.html

library(GenomicRanges)
library(nullranges)
library(tidyrverse)
library(dplyr)

match_cell <- function(input_table, sample_size) {

	table <- read.table(input_table, sep="\t", header=FALSE)
	colnames(table) <- c("Cell", "Group", "nFeatures", "Contacts", "nCount")

	# arbitrary chrom, start, end, just placeholder
	table$chrom <- "chr1"
	table$start <- seq(1, dim(table)[1])*100
	table$end   <- table$start + 99

	set.seed(1)

	focal_table <- sample_n(table[table$Group=="high", ], sample_size)
	focal_gr <- GRanges(
		seqnames = focal_table$chrom,
		ranges = IRanges(start=focal_table$start, end=focal_table$end),
		Cell = focal_table$Cell,
		Group = focal_table$Group,
		nFeatures = focal_table$nFeatures,
		Contacts = focal_table$Contacts,
		nCount = focal_table$nCount
	)

	pool_table <- table[table$Group=="low", ]
	pool_gr <- GRanges(
		seqnames = pool_table$chrom,
		ranges = IRanges(start=pool_table$start, end=pool_table$end),
		Cell = pool_table$Cell,
		Group = pool_table$Group,
		nFeatures = pool_table$nFeatures,
		Contacts = pool_table$Contacts,
		nCount = pool_table$nCount
	)

	mgr <- matchRanges(focal=focal_gr, pool=pool_gr, covar = ~ nFeatures, method="stratified", replace=FALSE)
	output <- rbind(data.frame(focal_gr), data.frame(mgr))
	return(output)
}





