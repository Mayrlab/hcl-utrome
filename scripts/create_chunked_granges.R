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
        output=list(plus="/fscratch/fanslerm/augmented.plus.chunked.e3.t200.gc39.pas3.f0.9999.Rds",
                    minus="/fscratch/fanslerm/augmented.minus.chunked.e3.t200.gc39.pas3.f0.9999.Rds",),
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

cat("Loading GENCODE...\n")
gr_gencode <- import(snakemake@input$gencode, genome='hg38') %>%
    keepStandardChromosomes(pruning.mode="coarse") %>%
    filter(gene_type == 'protein_coding')

cat("Loading upstream transcripts...\n")
gr_upstream <- import(snakemake@input$utr3, genome='hg38') %>%
    keepStandardChromosomes(pruning.mode="coarse")

cat("Loading extended transcripts...\n")
gr_downstream <- import(snakemake@input$extutr3, genome='hg38') %>%
    keepStandardChromosomes(pruning.mode="coarse")

cat("Merging plus strand gene models...\n")
gr_plus <- bind_ranges(gr_gencode, gr_upstream, gr_downstream) %>%
    filter(strand == "+") %>%
    filter(type == 'gene' | (type %in% c('transcript', 'exon') & transcript_type == 'protein_coding'))

cat("Chunking plus strand gene models...\n")
grl_plus <- gr_plus$gene_id %>%
    unique %>%
    { split(., ceiling(seq_along(.)/CHUNK_SIZE)) } %>%
    lapply(function (gene_ids) { gr_plus[gr_plus$gene_id %in% gene_ids] }) %>%
    as("GRangesList")

cat("Exporting chunked plus strand GRangesList...\n")
saveRDS(grl_plus, snakemake@output$plus)

cat("Merging minus strand gene models...\n")
gr_minus <- bind_ranges(gr_gencode, gr_upstream, gr_downstream) %>%
    filter(strand == "-") %>%
    filter(type == 'gene' | (type %in% c('transcript', 'exon') & transcript_type == 'protein_coding'))

cat("Chunking minus strand gene models...\n")
grl_minus <- gr_minus$gene_id %>%
    unique %>%
    { split(., ceiling(seq_along(.)/CHUNK_SIZE)) } %>%
    lapply(function (gene_ids) { gr_minus[gr_minus$gene_id %in% gene_ids] }) %>%
    as("GRangesList")

cat("Exporting chunked minus strand GRangesList...\n")
saveRDS(grl_minus, snakemake@output$minus)

cat("Done.\n")
