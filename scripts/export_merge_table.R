#!/usr/bin/env Rscript

################################################################################
## Mock `snakemake` preamble
################################################################################

if (interactive()) {
  library(methods)
  Snakemake <- setClass(
    "Snakemake", 
    slots=c(
      input='list', 
      output='list',
      params='list',
      wildcards='list',
      threads='numeric'
    )
  )
  snakemake <- Snakemake(
      input=list(gtf="data/gff/utrome.e3.t200.gc39.pas3.f0.9999.w500.gtf.gz"),
      output=list(tsv="/fscratch/fanslerm/utrome.test.m200.tsv"),
      params=list(genome="hg38"),
      wildcards=list(merge="200"),
      threads=1
  )
}

################################################################################
## Libraries and Parameters
################################################################################

library(txcutr)
library(BSgenome)
library(GenomicFeatures)

## convert arguments
minDistance <- as.integer(snakemake@wildcards$merge)

## load genome
message("Loading genome...")
bsg <- getBSgenome(snakemake@params$genome)

################################################################################
## Load Data
################################################################################

message("Converting GTF to TxDb...")
txdb <- makeTxDbFromGFF(file=snakemake@input$gtf, organism=organism(bsg))
txdb <- keepStandardChromosomes(txdb, pruning.mode="coarse")
seqlevelsStyle(txdb) <- "UCSC"

exportMergeTable(txdb, snakemake@output$tsv, minDistance=minDistance)
