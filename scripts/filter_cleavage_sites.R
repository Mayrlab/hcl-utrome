#!/usr/bin/env Rscript

library(rtracklayer)
library(tidyverse)
library(magrittr)

################################################################################
                                        # Fake Arguments
################################################################################

if (interactive()) {
    Snakemake <- setClass("Snakemake", slots=c(input='list', output='list',
                                               wildcards='list', threads='numeric'))
    snakemake <- Snakemake(
        input=list(bed="data/bed/cleavage-sites/utrome.cleavage.e3.t200.bed.gz",
                   probs="data/cleavage-sites/utrome.classification.e3.t200.tsv.gz"),
        output=list(bed="data/bed/utrome.cleavage.e3.t200.f0.999.bed.gz"),
        wildcards=list(epsilon="3", threshold="200", likelihood="0.9999"),
        threads=1
        )
}

################################################################################
                                        # Load Data
################################################################################

gr_peaks <- import(snakemake@input$bed) %>%
    `seqlevelsStyle<-`("UCSC") %>%
    keepStandardChromosomes(pruning.mode="coarse")
names(gr_peaks) <- elementMetadata(gr_peaks)$name

df_probs <- read_tsv(snakemake@input$probs, col_types='cddi')

################################################################################
                                        # Filter and Export
################################################################################

idx_passing <- df_probs %>%
    filter(prob_true_pA >= as.numeric(snakemake@wildcards$likelihood)) %$%
    peak_name

gr_filtered <- gr_peaks[idx_passing]

export.bed(object=gr_filtered, con=snakemake@output$bed)
