#!/usr/bin/env Rscript

library(rtracklayer)
library(readr)
library(dplyr)
library(stringr)
library(plyranges)
library(Matrix)

################################################################################
                                        # Fake Arguments
################################################################################

if (interactive()) {
    Snakemake <- setClass("Snakemake", slots=c(input='list', output='list',
                                               wildcards='list',
                                               params='list', threads='numeric'))
    snakemake <- Snakemake(
        input=list(unvalidated="data/bed/cleavage-sites/utrome.unvalidated.e3.t200.gc39.bed.gz",
                   atlas="data/cleavage-sites/polyAsite.atlas.tsv.gz"),
        output=list(supported="data/bed/cleavage-sites/utrome.supported.e3.t200.gc39.pas3.bed.gz",
                    unsupported="data/bed/cleavage-sites/utrome.unsupported.e3.t200.gc39.pas3.bed.gz"),
        wildcards=list(epsilon="3", threshold="200", version="39", tpm="3"),
        params=list(radius="20"),
        threads=1
        )
}

################################################################################
                                        # Load Data
################################################################################

## polyASite
df_pas <- read_tsv(snakemake@input$atlas)

## filter by tpm
tpm_pas <- df_pas %>%
    select(starts_with(c("GSM", "SAMEA", "SRX"))) %>%
    as.matrix %>% Matrix(sparse=TRUE) %>%
    `rownames<-`(df_pas$name)

idx_expressed <- which(rowSums(tpm_pas >= as.numeric(snakemake@wildcards$tpm)) > 0)

gr_pas <- df_pas[idx_expressed,] %>%
    select(chrom, chromStart, chromEnd, strand) %>%
    filter(str_detect(chrom, '^([0-9]+|X|Y|M|MT)$')) %>%
    mutate(chrom=str_c("chr", chrom)) %>%
    makeGRangesFromDataFrame(start.field="chromStart", end.field="chromEnd",
                             keep.extra.columns=TRUE)

## Load unvalidated
gr_unvalidated <- import(snakemake@input$unvalidate)

################################################################################
                                        # Filter Sites
################################################################################

gr_supported <- (gr_unvalidated + as.integer(snakemake@params$radius)) %>%
    find_overlaps_directed(gr_pas) %>%
    unique()

################################################################################
                                        # Export Sites
################################################################################

gr_unvalidated %>%
    filter(name %in% gr_supported$name) %>%
    write_bed(file=snakemake@output$supported)

gr_unvalidated %>%
    filter(!(name %in% gr_supported$name)) %>%
    write_bed(file=snakemake@output$unsupported)
