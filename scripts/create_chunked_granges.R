#!/usr/bin/env Rscript

library(rtracklayer)
library(plyranges)

################################################################################
                                        # Fake Arguments
################################################################################

if (interactive()) {
    Snakemake <- setClass("Snakemake", slots=c(input='list', output='list',
                                               wildcards='list',
                                               params='list', threads='numeric'))
    snakemake <- Snakemake(
        input=list(utr3="data/gff/txs.utr3.e3.t200.gc39.pas3.f0.9999.gff3.gz",
                   extutr3="data/gff/txs.extutr3.e3.t200.gc39.pas3.f0.9999.gff3.gz",
                   gencode="data/gff/gencode.v39.mRNA_ends_found.gff3.gz"),
        output=list(positive="/fscratch/fanslerm/augmented.positive.chunked.e3.t200.gc39.pas3.f0.9999.Rds",
                    negative="/fscratch/fanslerm/augmented.negative.chunked.e3.t200.gc39.pas3.f0.9999.Rds",),
        wildcards=list(epsilon="3", threshold="200", version="39", tpm="3", likelihood="0.9999"),
        params=list(),
        threads=1
        )
}

################################################################################
                                        # Load Data
################################################################################

## parameters
CHUNK_SIZE=100

message("Loading GENCODE...")
gr_gencode <- import(snakemake@input$gencode, genome='hg38') %>%
    keepStandardChromosomes(pruning.mode="coarse") %>%
    filter(gene_type == 'protein_coding')

message("Loading upstream transcripts...")
gr_upstream <- import(snakemake@input$utr3, genome='hg38') %>%
    keepStandardChromosomes(pruning.mode="coarse") %>%
    filter(gene_type == 'protein_coding')

message("Loading extended transcripts...")
gr_downstream <- import(snakemake@input$extutr3, genome='hg38') %>%
    keepStandardChromosomes(pruning.mode="coarse") %>%
    filter(gene_type == 'protein_coding')

message("Merging positive strand gene models...")
gr_positive <- bind_ranges(gr_gencode, gr_upstream, gr_downstream) %>%
    filter(strand == "+") %>%
    filter(type == 'gene' | (type %in% c('transcript', 'exon') & transcript_type == 'protein_coding'))

message("Chunking positive strand gene models...")
grl_positive <- gr_positive$gene_id %>%
    unique %>%
    { split(., ceiling(seq_along(.)/CHUNK_SIZE)) } %>%
    lapply(function (gene_ids) { gr_positive[gr_positive$gene_id %in% gene_ids] }) %>%
    as("GRangesList")

message("Exporting chunked positive strand GRangesList...")
saveRDS(grl_positive, snakemake@output$positive)

message("Merging negative strand gene models...")
gr_negative <- bind_ranges(gr_gencode, gr_upstream, gr_downstream) %>%
    filter(strand == "-") %>%
    filter(type == 'gene' | (type %in% c('transcript', 'exon') & transcript_type == 'protein_coding'))

message("Chunking negative strand gene models...")
grl_negative <- gr_negative$gene_id %>%
    unique %>%
    { split(., ceiling(seq_along(.)/CHUNK_SIZE)) } %>%
    lapply(function (gene_ids) { gr_negative[gr_negative$gene_id %in% gene_ids] }) %>%
    as("GRangesList")

message("Exporting chunked negative strand GRangesList...")
saveRDS(grl_negative, snakemake@output$negative)

message("Done.")
