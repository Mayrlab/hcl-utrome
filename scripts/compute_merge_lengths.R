#!/usr/bin/env Rscript

library(tidyverse)
library(plyranges)
library(magrittr)
library(BiocParallel)

################################################################################
                                        # Fake Arguments
################################################################################

if (interactive()) {
    Snakemake <- setClass("Snakemake", slots=c(input='list', output='list',
                                               wildcards='list',
                                               threads='numeric'))
    snakemake <- Snakemake(
        input=list(gtf="data/gff/utrome.e20.t5.gc39.pas3.f0.9999.w500.gtf.gz",
                   tsv="data/gff/utrome.e20.t5.gc39.pas3.f0.9999.w500.m200.tsv"),
        output=list(qc="qc/gff/utrome.merged_lengths.e20.t5.gc39.pas3.f0.9999.w500.m200.tsv"),
        wildcards=list(epsilon="20", threshold="5", version="39", tpm="3",
                       likelihood="0.9999", width="500", merge="200"),
        threads=3
        )
}

################################################################################
                                        # Load Data
################################################################################

## parameters
epsilon <- as.integer(snakemake@wildcards$epsilon)
min_tpm <- as.integer(snakemake@wildcards$threshold)

## set parallelization defaults
register(MulticoreParam(as.integer(snakemake@threads)))

## Load data
message("Loading merge table...")
df_merges <- read_tsv(snakemake@input$tsv, col_select=c("tx_in", "tx_out")) %>%
    dplyr::rename(merged_id=tx_out) %>%
    group_by(merged_id) %>%
    summarize(merging_ids=list(tx_in))

message("Loading UTRome annotation...")
gr_exons <- read_gff(snakemake@input$gtf) %>%
    filter(type == 'exon')

## Compute lengths
message("Computing measured region lengths...")
df_lengths <- df_merges %>%
    deframe() %>%
    bplapply(function (tx_ids) {
        gr_exons %>%
            filter(transcript_id %in% tx_ids) %>%
            reduce_ranges_directed() %>%
            width() %>%
            sum()
    }) %>%
    unlist() %>%    
    enframe(name="merged_id", value="nts_measured")

message("Finalizing table...")
df_out <- df_merges %>%
    left_join(df_lengths, by="merged_id") %>%
    mutate(n_merged=map_int(merging_ids, length),
           merging_ids=map_chr(merging_ids, str_flatten, ";"),
           epsilon=epsilon, min_tpm=min_tpm)

message("Exporting...")
write_tsv(df_out, snakemake@output$tsv)
