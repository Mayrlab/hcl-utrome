#!/usr/bin/env Rscript

library(BSgenome.Hsapiens.UCSC.hg38)
library(rtracklayer)
library(plyranges)
library(stringr)

################################################################################
                                        # Fake Arguments
################################################################################

if (interactive()) {
    Snakemake <- setClass("Snakemake", slots=c(input='list', output='list',
                                               wildcards='list',
                                               params='list', threads='numeric'))
    snakemake <- Snakemake(
        input=list(extutr3="data/bed/cleavage-sites/utrome.extutr3.e3.t200.gc39.pas3.f0.9999.bed.gz",
                   gencode="data/gff/gencode.v39.mRNA_ends_found.gff3.gz"),
        output=list(gff="data/gff/txs.extutr3.e3.t200.gc39.pas3.f0.9999.gff3.gz",
                    gtf="data/gff/txs.extutr3.e3.t200.gc39.pas3.f0.9999.gtf.gz"),
        wildcards=list(epsilon="3", threshold="200", version="39", tpm="3", likelihood="0.9999"),
        params=list(ext_utr3="5000"),
        threads=3
        )
}

################################################################################
                                        # Load Data
################################################################################

## Load References
hg38 <- BSgenome.Hsapiens.UCSC.hg38

cat("Loading GENCODE...\n")
gr_gencode <- import(snakemake@input$gencode, genome='hg38') %>%
    keepStandardChromosomes(pruning.mode="coarse")

cat("Loading downstream cleavage sites...\n")
gr_downstream <- import(snakemake@input$extutr3, genome='hg38') %>%
    keepStandardChromosomes(pruning.mode="coarse") %>%
    mutate(cleavage_site=end(.))

################################################################################
                                        # Intersect and Truncate
################################################################################

cat("Intersecting cleavage sites with extended 3' UTRs...\n")
gr_overlaps <- gr_gencode %>%
    filter(type == "three_prime_UTR", transcript_type == "protein_coding") %>%
    flank_downstream(as.integer(snakemake@params$ext_utr3)) %>%
    find_overlaps_directed(gr_downstream, suffix=c("_gencode", ""))

cat("Extracting intersected transcripts...\n")
gr_txs <- gr_gencode %>%
    filter(type == "transcript", transcript_id %in% gr_overlaps$transcript_id) %>%
    `names<-`(.$transcript_id)

## expand existing txs
cat("Generating new transcripts...\n")
gr_txs_new <- gr_txs[gr_overlaps$transcript_id,] %>%
    ## store previous id
    mutate(transcript_id_old=transcript_id) %>%
    ## copy new transcript ends
    mutate(cleavage_site=gr_overlaps$cleavage_site,
           name=gr_overlaps$name) %>%
    ## compute cleavage site offset
    mutate(offset=ifelse(strand == '+', cleavage_site - end, start - cleavage_site)) %>%
    ## filter to most parsimonious
    group_by(name) %>% filter(offset == min(offset)) %>% ungroup() %>%
    ## update tx id
    mutate(transcript_id=str_c(transcript_id, "-UTR+", offset),
           transcript_name=str_c(transcript_name, "-UTR+", offset),
           ID=str_c(ID, "-UTR+", offset)) %>%
    ## update names
    `names<-`(.$transcript_id) %>%
    ## save old 3' end
    mutate(old_3p_end=ifelse(strand == "+", end, start)) %>%
    ## update 3' ends
    mutate(start=ifelse(strand == "+", start, cleavage_site),
           end  =ifelse(strand == "+", cleavage_site, end)) %>%
    unique()

## extracting corresponding exons
cat("Extracting exons for transcripts...\n")
gr_exons <- gr_gencode %>%
    ## include copy of exon for every transcript
    expand(colnames=c("Parent")) %>%
    ## only consider exons from new parents
    filter(Parent %in% gr_txs_new$transcript_id_old) 

## discard unused objects
rm(list=c("gr_gencode", "gr_txs", "gr_downstream", "gr_overlaps"))
gc()

## update exons
cat("Generating new exons...\n")
gr_exons_new <- gr_txs_new %>%
    ## only need subset of columns
    select(transcript_id, transcript_name, transcript_id_old,
           cleavage_site, name, offset, old_3p_end) %>%
    ## intersect
    join_overlap_intersect_directed(x=gr_exons, suffix=c(".exon", "")) %>%
    ## only keep exact matches
    filter(transcript_id.exon == transcript_id_old) %>%
    ## update ID and Parent
    mutate(ID=str_replace(ID, coll(transcript_id_old), transcript_id),
           Parent=str_replace(Parent, coll(transcript_id_old), transcript_id)) %>%
    ## truncate
    mutate(start=ifelse(old_3p_end == start, cleavage_site, start),
           end  =ifelse(old_3p_end == end, cleavage_site, end)) %>%
    ## remove old columns
    select(-c("transcript_id.exon", "transcript_name.exon", "old_3p_end"))

gr_txs_new <- select(gr_txs_new, -c("old_3p_end"))

################################################################################
                                        # Export New Transcripts
################################################################################

cat("Exporting transcriptome for downstream sites...\n")
cat("  GFF3...\n")
f_gff <- gzfile(snakemake@output$gff, "w")
export(c(gr_txs_new, gr_exons_new), f_gff, format="gff3")

cat("  GTF...\n")
f_gtf <- gzfile(snakemake@output$gtf, "w")
export(c(gr_txs_new, gr_exons_new), f_gtf, format="gtf")
closeAllConnections()

cat("Done.\n")
