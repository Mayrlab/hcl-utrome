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
        input=list(ipa="data/gff/utrome.e30.t5.gc39.pas3.f0.9999.w500.ipa.tsv",
                   gtf="data/gff/utrome.e30.t5.gc39.pas3.f0.9999.w500.gtf.gz"),
        output=list(txs="data/gff/df_utrome_txs.e30.t5.gc39.pas3.f0.9999.w500.Rds",
                    genes="data/gff/df_utrome_genes.e30.t5.gc39.pas3.f0.9999.w500.Rds"),
        wildcards=list(epsilon="30", threshold="5", version="39", tpm="3",
                       likelihood="0.9999", width="500"),
        params=list(),
        threads=1
        )
}

################################################################################
                                        # Load Data
################################################################################

df_ipa <- read_tsv(snakemake@input$ipa, show_col_types=FALSE)

df_txs <- read_gff(snakemake@input$gtf, genome_info="hg38") %>%
    ## transcripts only
    filter(type == 'transcript') %>%
    as_tibble() %>%
    select(c("transcript_id", "transcript_name",
             "utr_name", "utr_rank",
             "gene_id", "gene_name")) %>%
    ## annotate intronic & novel isoforms; fix rank
    mutate(is_ipa=transcript_id %in% df_ipa$transcript_id,
           is_novel=str_detect(transcript_id, "UTR"),
           utr_rank=as.integer(utr_rank)) %>%
    ## annotate proximal, distal
    group_by(gene_id) %>%
    mutate(is_proximal=utr_rank == min(utr_rank[!is_ipa]),
           is_distal=utr_rank == max(utr_rank[!is_ipa]),
           utr_count_raw=dplyr::n(),
           utr_count_no_ipa=sum(!is_ipa)) %>%
    ungroup() %>%
    ## set utr type
    mutate(utr_type_raw=ifelse(utr_count_raw > 1, "multi", "single"),
           utr_type_no_ipa=ifelse(utr_count_no_ipa > 1, "multi", "single"))

df_genes <- df_txs %>%
    group_by(gene_id, gene_name) %>%
    summarize(utr_count_raw=dplyr::first(utr_count_raw),
              utr_count_no_ipa=dplyr::first(utr_count_no_ipa),
              utr_type_raw=dplyr::first(utr_type_raw),
              utr_type_no_ipa=dplyr::first(utr_type_no_ipa),
              has_ipa=any(is_ipa))
              
################################################################################
                                        # Export DataFrames
################################################################################

df_txs %>%
    DataFrame(row.names=.$transcript_id) %>%
    saveRDS(snakemake@output$txs)

df_genes %>%
    DataFrame(row.names=.$gene_id) %>%
    saveRDS(snakemake@output$genes)
