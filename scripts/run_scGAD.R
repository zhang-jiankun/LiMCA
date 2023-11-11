#!/usr/bin/env Rscript

options(future.globals.maxSize= 15*850*1024^2 )
suppressMessages(lapply(c( "tidyverse", "ggplot2", "umap", "cowplot", "patchwork"), library, character.only=TRUE))
library(BandNorm)
library(future.apply)

# cat /share/home/zhangjk/data/hsp_genome/genecode_v34_gene_coords.bed \
#    | awk -F"\t" '{OFS="\t"}{print $1, $2, $3, $6, $4}' | sort -k1,1V -k2,2n > scGAD/hg38_gene_coords.txt 

# scgad_df <- read.table('scGAD/GM12878_data_10000.txt', sep="\t", header=FALSE)
# colnames(scgad_df) <- c("chrom", "binA", "binB", "count", "cell")
# saveRDS(scgad_df, file="scGAD/GM12878_data_10000.rds")

scgad_df <- readRDS("scGAD/GM12878_data_10000.rds")
hg38 <- read.table("scGAD/hg38_gene_coords.uniq.txt", sep="\t", header=FALSE)
colnames(hg38) <- c("chr", "s1", "s2", "strand", "gene_name")
geneAnno <- hg38 %>% filter(s2 - s1 >= 100000)

gad_score = scGAD(hic_df = scgad_df, genes = geneAnno, depthNorm = TRUE, cores=36, threads=72)
write.table(gad_score, file="scGAD/GM12878_scgad_score_10000.100000.txt", sep="\t", quote=FALSE)

custom_scGAD = function(path = NULL, hic_df = NULL, genes, depthNorm = TRUE, cores = 4, threads = 8, binPair = TRUE, format = "short", res = 10000){
  setDTthreads(threads)
  discardCounts = max(genes$s2 - genes$s1)
  genes$s1 <- ifelse(genes$strand == "+", genes$s1 - 1000, genes$s1)
  genes$s2 <- ifelse(genes$strand == "-", genes$s2, genes$s2 + 1000)
  colnames(genes) = c("chr", "start", "end", "strand", "names")
  genes = makeGRangesFromDataFrame(genes, keep.extra.columns=TRUE)
  g_names = genes$names
  if (is.null(hic_df)){
    if (binPair){
      names = basename(list.files(path, recursive = TRUE))
      paths = list.files(path, full.names = TRUE, recursive = TRUE)
      getCount = function(k){
        cell = fread(paths[k], select = c(1, 2, 4, 5))
        colnames(cell) = c("V1", "V2", "V4", "V5")
        cell = cell[abs(V4 - V2) <= discardCounts]
        GInt = GenomicInteractions(GRanges(cell$V1,
                                           IRanges(cell$V2, width = res)),
                                   GRanges(cell$V1,
                                           IRanges(cell$V4, width = res)),
                                   counts = cell$V5)
        one <- overlapsAny(anchorOne(GInt), genes)
        two <- overlapsAny(anchorTwo(GInt), genes)
        x.valid <- GInt[one & two]
        hits <- list()
        hits$one <- findOverlaps(anchorOne(x.valid), genes, select = "first")
        hits$two <- findOverlaps(anchorTwo(x.valid), genes, select = "first")
        counts = data.table(reads = x.valid$counts[hits[[1]] == hits[[2]]], pos = hits$one[hits$one == hits$two])
        tabulated <- unique(counts$pos)
        counts <- setDT(counts)[,.(reads = sum(reads)), by = 'pos']$reads
        dat = data.table(names = c(g_names[unique(tabulated)], g_names[-unique(tabulated)]), 
                         counts = c(counts, rep(0, length(g_names) - length(unique(tabulated)))))
        dat[match(g_names, dat$names), ]$counts
      }
      plan(multicore, workers = cores)
      output <- future_sapply(1:length(names), getCount)
      output[is.na(output)] = 0
    }
    else {
      names = basename(list.files(path, recursive = TRUE))
      paths = list.files(path, full.names = TRUE, recursive = TRUE)
      getCount = function(k){
        cell = fread(paths[k])
        if (format == "short"){
          cell = cell[cell$V2 == cell$V6, ]
          cell = cell[, c(2, 3, 7)]
        }else if (format == "medium"){
          cell = cell[cell$V3 == cell$V7, ]
          cell = cell[, c(3, 4, 8)]
        }else if (format == "long"){
          cell = cell[cell$V2 == cell$V6, ]
          cell = cell[, c(2, 3, 7)]
        }else if (format == "4DN"){
          cell = cell[cell$V2 == cell$V4, ]
          cell = cell[, c(2, 3, 5)]
        }
        colnames(cell) = c("V1", "V2", "V4")
        cell = cell[abs(V4 - V2) <= discardCounts]
        GInt = GenomicInteractions(GRanges(cell$V1,
                                           IRanges(cell$V2, width = res)),
                                   GRanges(cell$V1,
                                           IRanges(cell$V4, width = res)),
                                   counts = 1)
        one <- overlapsAny(anchorOne(GInt), genes)
        two <- overlapsAny(anchorTwo(GInt), genes)
        x.valid <- GInt[one & two]
        hits <- list()
        hits$one <- findOverlaps(anchorOne(x.valid), genes, select = "first")
        hits$two <- findOverlaps(anchorTwo(x.valid), genes, select = "first")
        counts = data.table(reads = x.valid$counts[hits[[1]] == hits[[2]]], pos = hits$one[hits$one == hits$two])
        tabulated <- unique(counts$pos)
        counts <- setDT(counts)[,.(reads = sum(reads)), by = 'pos']$reads
        dat = data.table(names = c(g_names[unique(tabulated)], g_names[-unique(tabulated)]), 
                         counts = c(counts, rep(0, length(g_names) - length(unique(tabulated)))))
        dat[match(g_names, dat$names), ]$counts
      }
      plan(multicore, workers = cores)
      output <- future_sapply(1:length(names), getCount)
      output[is.na(output)] = 0
    }
  } else{
    hic_df = setDT(hic_df)
    names = unique(hic_df$cell)
    getCount = function(k){
      cell = hic_df[cell == names[k]]
      cell = cell[abs(binB - binA) <= discardCounts]
      GInt = GenomicInteractions(GRanges(cell$chrom,
                                         IRanges(cell$binA, width = res)),
                                 GRanges(cell$chrom,
                                         IRanges(cell$binB, width = res)),
                                 counts = cell$count)
      one <- overlapsAny(anchorOne(GInt), genes)
      two <- overlapsAny(anchorTwo(GInt), genes)
      x.valid <- GInt[one & two]
      hits <- list()
      hits$one <- findOverlaps(anchorOne(x.valid), genes, select = "first")
      hits$two <- findOverlaps(anchorTwo(x.valid), genes, select = "first")
      counts = data.table(reads = x.valid$counts[hits[[1]] == hits[[2]]], pos = hits$one[hits$one == hits$two])
      tabulated <- unique(counts$pos)
      counts <- aggregate(reads ~ pos, data = counts, FUN = sum)$reads
      dat = data.table(names = c(g_names[unique(tabulated)], g_names[-unique(tabulated)]), 
                       counts = c(counts, rep(0, length(g_names) - length(unique(tabulated)))))
      dat[match(g_names, dat$names), ]$counts
    }
    plan(multicore, workers = cores)
    output <- future_sapply(1:length(names), getCount)
    output[is.na(output)] = 0
  }
  rownames(output) = g_names
  colnames(output) = names
  output = output[rowSums(output) > 0, ]
  output = output[!is.na(rowSums(output)), ]
  if (depthNorm) {
    output = t(t(output) / colSums(output)) * 1e06
  }
  # GAD = (output - rowMeans(output))/sqrt(rowVars(output))
  output
}

# raw GAD score
gad_score = custom_scGAD(hic_df = scgad_df, genes = geneAnno, depthNorm = FALSE, cores=36, threads=72)
write.table(gad_score, file="scGAD/GM12878_scgad_score_10000.100000.raw.txt", sep="\t", quote=FALSE)
