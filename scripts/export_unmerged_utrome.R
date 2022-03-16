#!/usr/bin/env Rscript

library(GenomicFeatures)
library(BSgenome.Hsapiens.UCSC.hg38)
library(rtracklayer)
library(plyranges)
library(stringr)
library(dplyr)
library(magrittr)
library(BiocParallel)

################################################################################
                                        # Fake Arguments
################################################################################

if (interactive()) {
    Snakemake <- setClass("Snakemake", slots=c(input='list', output='list',
                                               wildcards='list',
                                               params='list', threads='numeric'))
    snakemake <- Snakemake(
        input=list(positive="data/granges/utrome.raw.positive.e3.t200.gc39.pas3.f0.9999.w500.Rds",
                   negative="data/granges/utrome.raw.negative.e3.t200.gc39.pas3.f0.9999.w500.Rds"),
        output=list(gtf="/fscratch/fanslerm/utrome.e3.t200.gc39.pas3.f0.9999.w500.gtf",
                    fa="/fscratch/fanslerm/utrome.e3.t200.gc39.pas3.f0.9999.w500.fa.gz"),
        wildcards=list(epsilon="3", threshold="200", version="39", tpm="3",
                       likelihood="0.9999", width="500"),
        params=list(),
        threads=3
        )
}

################################################################################
                                        # Load Data
################################################################################

## parameters
MAX_TX_LENGTH=as.integer(snakemake@wildcards$width)
CHUNK_SIZE=100

## set parallelization defaults
register(MulticoreParam(as.integer(snakemake@threads)))

## Load References
hg38 <- BSgenome.Hsapiens.UCSC.hg38

## Load data
message("Loading truncated transcripts...")
message("  strand: +")
gr_positive <- readRDS(snakemake@input$positive) %>%
    `names<-`(NULL)

message("  strand: -")
gr_negative <- readRDS(snakemake@input$negative) %>%
    `names<-`(NULL)

gr_genes <- filter(c(gr_positive, gr_negative), type == 'gene') %>%
    `names<-`(.$ID)

gr_txs <- filter(c(gr_positive, gr_negative), type == 'transcript') %>%
    `names<-`(.$ID)

grl_exons <- filter(c(gr_positive, gr_negative), type == 'exon') %>%
    expand_ranges(Parent) %>%
    split(.$Parent) %>%
    as("GRangesList")

################################################################################
                                        # Remove exact duplicates
################################################################################

df_intersect <- findOverlaps(grl_exons, minoverlap=MAX_TX_LENGTH,
                             ignore.strand=FALSE, drop.self=TRUE) %>%
    ## convert to tibble
    as.data.frame() %>% tibble::as_tibble() %>%
    ## consider each pair only once
    filter(queryHits < subjectHits) %>%
    ## determine tx id
    mutate(tx_query=names(grl_exons)[queryHits],
           tx_subject=names(grl_exons)[subjectHits]) %>%
    ## check if adjusted and get source
    mutate(is_adjusted_query=str_detect(tx_query, '(-|\\+)[0-9]+$'),
           is_adjusted_subject=str_detect(tx_subject, '(-|\\+)[0-9]+$'),
           source_query=mcols(gr_txs)[tx_query, "source"],
           source_subject=mcols(gr_txs)[tx_subject, "source"]) %>%
    ## compute adjustment
    mutate(adjustment_query=ifelse(is_adjusted_query, as.integer(str_extract(tx_query, '(-|\\+)[0-9]+$')), 0L),
           adjustment_subject=ifelse(is_adjusted_subject, as.integer(str_extract(tx_subject, '(-|\\+)[0-9]+$')), 0L))

## identify duplicates to remove
idx_tx_remove <- df_intersect %$%
    case_when(
        ## prefer non-adjusted
        (!is_adjusted_query) & ( is_adjusted_subject) ~ tx_subject,
        ( is_adjusted_query) & (!is_adjusted_subject) ~ tx_query,
        ## subsequent tests assume is_augmented_* matches
        ## prefer least adjustment
        abs(adjustment_query) < abs(adjustment_subject) ~ tx_subject,
        abs(adjustment_query) > abs(adjustment_subject) ~ tx_query,
        ## subseqeunt tests assume adjustment_* matches
        ## prefer "HAVANA" (manually curated) txs
        source_query == "HAVANA" & source_subject != "HAVANA" ~ tx_subject,
        source_query != "HAVANA" & source_subject == "HAVANA" ~ tx_query,
        ## prefer earlier txs
        tx_query < tx_subject ~ tx_subject,
        tx_query > tx_subject ~ tx_query,
        ## default
        TRUE ~ sprintf("Unmatched case: %s,%s", tx_query, tx_subject)
    )

## remove duplicates
gr_txs_dedup <- gr_txs[!(names(gr_txs) %in% idx_tx_remove)]
grl_exons_dedup <- grl_exons[!(names(grl_exons) %in% idx_tx_remove)]

################################################################################
                                        # Generate UTR names by position
################################################################################

## positive strand
grl_txs_positive <- gr_txs_dedup %>% 
    filter(strand == "+") %>%
    { unique(.$gene_id) } %>%
    { split(., ceiling(seq_along(.)/CHUNK_SIZE)) } %>%
    lapply(function (gene_ids) { gr_txs_dedup[gr_txs_dedup$gene_id %in% gene_ids] }) %>%
    as("GRangesList")

gr_txs_positive <- bplapply(grl_txs_positive, function (gr) {
    gr %>%
        group_by(gene_id) %>%
        mutate(utr_rank=row_number(end)) %>%
        mutate(utr_name=str_c(gene_name, utr_rank, sep=".")) %>%
        ungroup()
}) %>%
    as("GRangesList") %>% unlist() %>%
    `names<-`(NULL)

## negative strand
grl_txs_negative <- gr_txs_dedup %>% 
    filter(strand == "-") %>%
    { unique(.$gene_id) } %>%
    { split(., ceiling(seq_along(.)/CHUNK_SIZE)) } %>%
    lapply(function (gene_ids) { gr_txs_dedup[gr_txs_dedup$gene_id %in% gene_ids] }) %>%
    as("GRangesList")

gr_txs_negative <- bplapply(grl_txs_negative, function (gr) {
    gr %>%
        group_by(gene_id) %>%
        mutate(utr_rank=row_number(-start)) %>%
        mutate(utr_name=str_c(gene_name, utr_rank, sep=".")) %>%
        ungroup()
}) %>%
    as("GRangesList") %>% unlist() %>%
    `names<-`(NULL)

################################################################################
                                        # Prepare features
################################################################################

message("Preparing transcripts...")
gr_utrome_txs <- bind_ranges(gr_txs_positive, gr_txs_negative) %>%
    `names<-`(NULL) %>%
    select(type, ID, Parent,
           gene_id, gene_name,
           transcript_id, transcript_name,
           utr_name, utr_rank,
           exon_id)

message("Preparing genes...")
gr_utrome_genes <- gr_genes %>%
    filter(ID %in% gr_utrome_txs$gene_id) %>%
    mutate(utr_name=NA_character_, utr_rank=NA_integer_) %>%
    select(type, ID, Parent,
           gene_id, gene_name,
           transcript_id, transcript_name,
           utr_name, utr_rank,
           exon_id)

message("Preparing exons...")
utr_names <- set_names(gr_utrome_txs$utr_name, nm=gr_utrome_txs$transcript_id)
utr_ranks <- set_names(gr_utrome_txs$utr_rank, nm=gr_utrome_txs$transcript_id)

gr_utrome_exons <- grl_exons[gr_utrome_txs$transcript_id] %>%
    unlist() %>% `names<-`(NULL) %>%
    mutate(utr_name=utr_names[transcript_id], utr_rank=utr_ranks[transcript_id]) %>%
    select(type, ID, Parent,
           gene_id, gene_name,
           transcript_id, transcript_name,
           utr_name, utr_rank,
           exon_id)

################################################################################
                                        # Export Annotation
################################################################################

message("Exporting UTRome GTF...")
utrome_name <- str_c("utrome-hcl-gencode_v", snakemake@wildcards$version)
gr_utrome_genes %>% sort() %>%
    export(snakemake@output$gtf, format="GTF", source=utrome_name)
gr_utrome_txs %>% sort() %>%
    export(snakemake@output$gtf, format="GTF", append=TRUE, source=utrome_name)
gr_utrome_exons %>% filter(strand == '+') %>% sort() %>%
    export(snakemake@output$gtf, format="GTF", append=TRUE, source=utrome_name)
gr_utrome_exons %>% filter(strand == '-') %>% sort(decreasing=TRUE) %>%
    export(snakemake@output$gtf, format="GTF", append=TRUE, source=utrome_name)

################################################################################
                                        # Generate and Export Sequences
################################################################################

message("Creating TxDb object from UTRome...")
## Combine hierarchy of elements
gr_utrome <- bind_ranges(gr_utrome_genes, gr_utrome_txs, gr_utrome_exons)
txdb_utrome <- makeTxDbFromGRanges(gr_utrome, taxonomyId=9606,
                                   metadata=data.frame(name=c("Genome"), value=c("hg38")))

message("Exporting UTRome FASTA...\n")
seqs_utrome <- extractTranscriptSeqs(hg38, txdb_utrome, use.names=TRUE)
writeXStringSet(seqs_utrome, snakemake@output$fa, format="fasta", compress=TRUE)

message("All operations complete.")
