#!/usr/bin/env Rscript

library(plyranges)
library(stringr)
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
        input=list(ipa="data/gff/utrome.e20.t5.gc39.pas3.f0.9999.w500.tsv",
                   gtf="data/gff/utrome.e20.t5.gc39.pas3.f0.9999.w500.gtf.gz"),
        output=list(gr="data/granges/utrome_gr_txs.e20.t5.gc39.pas3.f0.9999.w500.Rds"),
        wildcards=list(epsilon="20", threshold="5", version="39", tpm="3",
                       likelihood="0.9999", width="500"),
        params=list(),
        threads=1
        )
}

################################################################################
                                        # Load Data
################################################################################

df_ipa <- read_tsv(snakemake@input$ipa, show_col_types=FALSE)

gr_cleavage_sites <- read_gff(snakemake@input$gtf, genome_info="hg38") %>%
    ## transcripts only
    filter(type == 'transcript') %>%
    select(-c("type")) %>%
    ## annotate intronic & novel isoforms; fix rank
    mutate(is_ipa=transcript_id %in% df_ipa$transcript_id,
           is_novel=str_detect(transcript_id, "UTR"),
           utr_rank=as.integer(utr_rank)) %>%
    ## annotate proximal, distal
    group_by(gene_id) %>%
    mutate(is_proximal=utr_rank == min(utr_rank[!is_ipa]),
           is_distal=utr_rank == max(utr_rank[!is_ipa]),
           utr_count=plyranges::n()) %>%
    ungroup() %>%
    ## set utr type
    mutate(utr_type=ifelse(utr_count > 1, "multi", "single"))

################################################################################
                                        # Export GRanges
################################################################################

saveRDS(gr_cleavage_sites, snakemake@output$gr)
