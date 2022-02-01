#!/usr/bin/env Rscript

library(magrittr)
library(dplyr)
library(plyranges)

################################################################################
                                        # Fake Arguments
################################################################################

if (interactive()) {
    Snakemake <- setClass("Snakemake", slots=c(input='list', output='list',
                                               wildcards='list',
                                               params='list', threads='numeric'))
    snakemake <- Snakemake(
        input=list(supported="data/bed/cleavage-sites/utrome.supported.e3.t200.gc39.pas3.bed.gz",
                   likely="data/bed/cleavage-sites/utrome.likely.e3.t200.gc39.pas3.f0.9999.bed.gz",
                   gencode="data/gff/gencode.v39.mRNA_ends_found.gff3.gz"),
        output=list(utr3="data/bed/cleavage-sites/utrome.utr3.e3.t200.gc39.pas3.f0.9999.bed.gz",
                    extutr3="data/bed/cleavage-sites/utrome.extutr3.e3.t200.gc39.pas3.f0.9999.bed.gz"),
        wildcards=list(epsilon="3", threshold="200", likelihood="0.9999", tpm="3", version="39"),
        params=list(ext_utr3="5000", ext_utr5="1000"),
        threads=3
        )
}

################################################################################
                                        # Load Data
################################################################################

## load annotation
gr_gencode <- read_gff3(snakemake@input$gencode) %>%
    keepStandardChromosomes(pruning.mode="coarse") %>%
    filter(type != 'gene')

## load candidate sites
gr_sites <- c(read_bed(snakemake@input$supported),
              read_bed(snakemake@input$likely))

################################################################################
                                        # Annotate Sites with GENCODE
################################################################################

#' NB: The left join operation expands the GRanges giving an entry for each
#' overlap that is found with the feature set. We can use `unique` to collapse
#' back down to only the original set. This has the effect that only the first
#' match to a feature is retained, since the subsequent ones are marked as
#' `duplicated` and left out.
#'
#' It would be preferably to collapse the column of interest to a CharacterList
#' or simply character with delimiters, but neither plyranges nor GRanges
#' appear to support this functionality without converting out to a data.frame-
#' like object and converting back.

N_SITES = length(gr_sites)
validate_site_count <- function (gr) {
    stopifnot(length(gr) == N_SITES)
    gr
}

## 3' UTR 
gr_sites <- gr_gencode %>%
    ## we only want to extend protein coding genes
    filter(type == "three_prime_UTR", transcript_type == "protein_coding") %>%
    mutate(utr3=gene_id) %>% select(utr3) %>% unique() %>%
    join_overlap_left_directed(x=gr_sites) %>%
    unique() %>%
    validate_site_count()

## 5' UTR 
gr_sites <- gr_gencode %>%
    filter(type == "five_prime_UTR") %>%
    mutate(utr5=gene_id) %>% select(utr5) %>% unique() %>%
    join_overlap_left_directed(x=gr_sites) %>%
    unique() %>%
    validate_site_count()

## exon
gr_sites <- gr_gencode %>%
    filter(type == "exon") %>%
    mutate(exon=gene_id) %>% select(exon) %>% unique() %>%
    join_overlap_left_directed(x=gr_sites) %>%
    unique() %>%
    validate_site_count()

## transcript
gr_sites <- gr_gencode %>%
    filter(type == "transcript") %>%
    mutate(transcript=gene_id) %>% select(transcript) %>% unique() %>%
    join_overlap_left_directed(x=gr_sites) %>%
    unique() %>%
    validate_site_count()

## extended 5' UTR 
gr_sites <- gr_gencode %>%
    filter(type == "five_prime_UTR") %>%
    flank_upstream(as.integer(snakemake@params$ext_utr5)) %>%
    mutate(extutr5=gene_id) %>% select(extutr5) %>% unique() %>%
    join_overlap_left_directed(x=gr_sites) %>%
    unique() %>%
    validate_site_count()

## extended 3' UTR 
gr_sites <- gr_gencode %>%
    filter(type == "three_prime_UTR", transcript_type == "protein_coding") %>%
    flank_downstream(as.integer(snakemake@params$ext_utr3)) %>%
    mutate(extutr3=gene_id) %>% select(extutr3) %>% unique() %>%
    join_overlap_left_directed(x=gr_sites) %>%
    unique() %>%
    validate_site_count()

################################################################################
                                        # Finalize Annotation
################################################################################

gr_sites$site_type <- mcols(gr_sites) %>%
    as_tibble() %$%
    case_when(
        !is.na(utr3) ~ "utr3",
        !is.na(utr5) ~ "utr5",
        !is.na(exon) ~ "exon",
        !is.na(transcript) ~ "intron",
        !is.na(extutr5) ~ "extutr5",
        !is.na(extutr3) ~ "extutr3",
        TRUE ~ "intergenic"
    )

################################################################################
                                        # Export Annotated Sites
################################################################################

## Export terminal exon sites
gr_sites %>%
    filter(site_type == "utr3") %>%
    write_bed(snakemake@output$utr3)

## Export extended 3' UTR sites
gr_sites %>%
    filter(site_type == "extutr3") %>%
    write_bed(snakemake@output$extutr3)
