#!/usr/bin/env Rscript

library(GenomicFeatures)
library(plyranges)
library(dplyr)
library(readr)

################################################################################
                                        # Fake Arguments
################################################################################

if (interactive()) {
    Snakemake <- setClass("Snakemake", slots=c(input='list', output='list',
                                               wildcards='list',
                                               params='list', threads='numeric'))
    snakemake <- Snakemake(
        input=list(utrome="data/gff/utrome.e20.t5.gc39.pas3.f0.9999.w500.gtf.gz",
                   gencode="data/gff/gencode.v39.mRNA_ends_found.gff3.gz"),
        output=list(tsv="data/gff/utrome.e20.t5.gc39.pas3.f0.9999.w500.ipa.tsv"),
        wildcards=list(epsilon="20", threshold="5", version="39", tpm="3", likelihood="0.9999"),
        params=list(),
        threads=1
        )
}

################################################################################
                                        # Load Data
################################################################################

gr_introns <- makeTxDbFromGFF(snakemake@input$gencode) %>%
    intronicParts()

gr_cleavage_sites <- read_gff(snakemake@input$utrome,
                              c('type', 'transcript_id', 'gene_id')) %>%
    filter(type == 'transcript') %>%
    select(-c("type")) %>%
    anchor_3p() %>%
    mutate(width=0)

################################################################################
                                        # Intersect cleavage sites
################################################################################

df_ipa <- join_overlap_intersect_directed(gr_cleavage_sites, gr_introns) %>%
    filter(mapply(function (x,y) {x %in% y}, x=gene_id.x, y=gene_id.y)) %>%
    mcols() %>% as_tibble() %>%
    select(transcript_id, gene_id.x) %>%
    rename(gene_id=gene_id.x)

################################################################################
                                        # Export
################################################################################

write_tsv(df_ipa, snakemake@output$tsv)
