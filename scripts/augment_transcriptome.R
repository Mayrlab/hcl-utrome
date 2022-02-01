#!/usr/bin/env Rscript

library(BSgenome.Hsapiens.UCSC.hg38)
library(rtracklayer)
library(plyranges)
library(tibble)
library(dplyr)
library(stringr)
library(BiocParallel)

################################################################################
                                        # Fake Arguments
################################################################################

if (interactive()) {
    Snakemake <- setClass("Snakemake", slots=c(input='list', output='list',
                                               wildcards='list',
                                               params='list', threads='numeric'))
    snakemake <- Snakemake(
        input=list(utr3="data/bed/cleavage-sites/utrome.utr3.e3.t200.gc39.pas3.f0.9999.bed.gz",
                   extutr3="data/bed/cleavage-sites/utrome.extutr3.e3.t200.gc39.pas3.f0.9999.bed.gz",
                   gencode="data/gff/gencode.v39.mRNA_ends_found.gff3.gz"),
        output=list(gff_utr3="data/gff/txs.utr3.e3.t200.gc39.pas3.f0.9999.gff3.gz",
                    gtf_utr3="data/gff/txs.utr3.e3.t200.gc39.pas3.f0.9999.gtf.gz",
                    gff_extutr3="data/gff/txs.extutr3.e3.t200.gc39.pas3.f0.9999.gff3.gz",
                    gtf_extutr3="data/gff/txs.extutr3.e3.t200.gc39.pas3.f0.9999.gtf.gz"),
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

cat("Loading upstream cleavage sites...\n")
gr_upstream <- import(snakemake@input$utr3, genome='hg38') %>%
    keepStandardChromosomes(pruning.mode="coarse") %>%
    mutate(cleavage_site=end(.))

################################################################################
                                        # Intersect and Truncate
################################################################################

cat("Intersecting cleavage sites with 3' UTRs...\n")
gr_overlaps <- gr_gencode %>%
    filter(type == "three_prime_UTR", transcript_type == "protein_coding") %>%
    find_overlaps_directed(gr_upstream, suffix=c("_gencode", ""))

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
    group_by(name) %>% filter(offset == max(offset)) %>% ungroup() %>%
    ## update tx id
    mutate(transcript_id=str_c(transcript_id, "-UTR", offset),
           transcript_name=str_c(transcript_name, "-UTR", offset),
           ID=str_c(ID, "-UTR", offset)) %>%
    ## update names
    `names<-`(.$transcript_id) %>%
    ## update 3' ends
    mutate(start=ifelse(strand == "+", start, cleavage_site),
           end  =ifelse(strand == "+", cleavage_site, end))

## extracting corresponding exons
cat("Extracting exons for transcripts...\n")
gr_exons <- gr_gencode %>%
    ## include copy of exon for every transcript
    expand(colnames=c("Parent")) %>%
    ## only consider exons from new parents
    filter(Parent %in% gr_txs_new$transcript_id_old) 

## update exons
cat("Generating new exons...\n")
gr_exons_new <- gr_txs_new %>%
    ## only need subset of columns
    select(transcript_id, transcript_name, cleavage_site, offset, transcript_id_old) %>%
    ## intersect
    join_overlap_intersect_directed(x=gr_exons, suffix=c(".exon", "")) %>%
    ## only keep exact matches
    filter(transcript_id.exon == transcript_id_old) %>%
    ## update ID and Parent
    mutate(ID=str_replace(ID, coll(transcript_id_old), transcript_id),
           Parent=str_replace(Parent, coll(transcript_id_old), transcript_id)) %>%
    ## truncate
    mutate(start=ifelse(strand == "+", start, ifelse(start > cleavage_site, start, cleavage_site)),
           end  =ifelse(strand == "-",   end, ifelse(end   < cleavage_site,   end, cleavage_site))) %>%
    ## remove old columns
    select(-c("transcript_id.exon", "transcript_name.exon"))


cat("Exporting transcriptome for upstream sites...\n")
cat("  GFF3...\n")
write_gff3(c(gr_txs_new, gr_exons_new), snakemake@output$gff3)

cat("  GTF...\n")
f <- gzfile(snakemake@output$gtf, "w")
export(c(gr_txs_new, gr_exons_new), f, format="gtf")

cat("Done.\n")


##adjust.transcriptEnd <- function (tx.id, end.new, tx.strand) {
##    stopifnot(tx.strand %in% c('-', '+'))

##    tx.gr <- annotation.gr %Q% (!is.na(transcript_id) & transcript_id == tx.id)
##    tx <- tx.gr[tx.gr$type == 'transcript',]

    ## LABELS
##    offset <- if (tx.strand == '+') {
##                  end.new - end(gr.end(tx, ignore.strand = FALSE))
##              } else {
##                  end(gr.end(tx, ignore.strand = FALSE)) - end.new
##              }
    ## update transcript_id
##    tx.id.new <- paste0(tx$transcript_id, "-UTR", if (offset < 0) "" else "+", offset)
##    elementMetadata(tx.gr)[,"transcript_id"] <- tx.id.new
    ## update transcript_name
##    tx.name.new <- paste0(tx$transcript_name, "-UTR", if (offset < 0) "" else "+", offset)
##    elementMetadata(tx.gr)[,"transcript_name"] <- tx.name.new
    ## update ID
##    elementMetadata(tx.gr)[, "ID"] <- str_replace(tx.gr$ID, tx.id, tx.id.new)
    ## update Parent
##    elementMetadata(tx.gr)[, "Parent"] <- str_replace(tx.gr$Parent, tx.id, tx.id.new)

    ## ADJUST TX END
##    end.cur <- end(gr.end(tx, ignore.strand = FALSE))
    ##cat(sprintf("Adjusting %s on %s strand from %d to %d\n", tx.id, tx.strand, end.cur, end.new))

    if (offset > 0) {
        if (tx.strand == '+') {
            ## Find all elements contiguous w/ TX end and extend range end
            end(tx.gr[end(gr.end(tx.gr, ignore.strand = FALSE)) == end.cur,]) <- end.new
        } else {
            ## Find all elements contiguous w/ TX end and extend range end
            start(tx.gr[end(gr.end(tx.gr, ignore.strand = FALSE)) == end.cur,]) <- end.new
        }
##    } else {
##        if (tx.strand == "+") {
            ## Delete all elements that might get removed (e.g., if multiple exons in 3' UTR)
##            tx.gr <- tx.gr[end(gr.start(tx.gr, ignore.strand = FALSE)) < end.new,]
            ## For remaining ranges, clip at new site
##            end(tx.gr[end(gr.end(tx.gr, ignore.strand = FALSE)) > end.new,]) <- end.new
##        } else {
            ## Delete all elements that might get removed (e.g., if multiple exons in 3' UTR)
##            tx.gr <- tx.gr[end(gr.start(tx.gr, ignore.strand = FALSE)) > end.new,]
            ## For remaining ranges, clip at new site
##            start(tx.gr[end(gr.end(tx.gr, ignore.strand = FALSE)) < end.new,]) <- end.new
##        }
    }

    tx.gr
}


cat("Extracting transcripts for upstream sites...\n")
upstream.annotation.gr <- grbind(
    mcmapply(adjust.transcriptEnd,
             tx.id = upstream.txs.gr$transcript_id,
             end.new = end(gr.end(upstream.txs.gr, ignore.strand = FALSE)),
             tx.strand = as.character(strand(upstream.txs.gr)),
             mc.cores = arg.cores, SIMPLIFY = FALSE)
)

cat("Exporting transcriptome for upstream sites...\n")
upstream.annotation.gff3.file <- sprintf("data/gff/%s.txs.utr3.%s.gff3", arg.outPrefix, arg.outSuffix)
upstream.annotation.gtf.file <- sprintf("data/gff/%s.txs.utr3.%s.gtf", arg.outPrefix, arg.outSuffix)
export.gff3(upstream.annotation.gr, upstream.annotation.gff3.file)
export(upstream.annotation.gr, upstream.annotation.gtf.file, format = "gtf")

write_gff(
cat("Clearing cached files...\n")
rm(list=c("upstream.annotation.gr", "upstream.txs.gr", "upstreamUTR.sites.gr"))

    
## Process Extended UTR Sites
cat("Intersecting annotation with downstream sites...\n")
extended.txs.gr <- gr.findoverlaps(
    extendedUTR.sites.gr,
    flank(gr.end(annotation.gr %Q%
                 (type == 'three_prime_UTR' &
                  transcript_type == 'protein_coding'),
                 ignore.strand = FALSE),
          width = arg.extDownstream, start = FALSE),
    scol = c("transcript_id"), ignore.strand = FALSE)

if (all(1:length(extendedUTR.sites.gr) %in% extended.txs.gr$query.id)) {
    cat("All cleavage sites intersected within 5Kb downstream of known 3' UTRs\n")
} else {
    warning("Some cleavage sites in extendedUTR.sites.gr did not intersect known sites!")
}

cat("Extracting transcripts for extended transcriptome...\n")
extended.annotation.gr <- grbind(
    mcmapply(adjust.transcriptEnd,
             tx.id = extended.txs.gr$transcript_id,
             end.new = end(gr.end(extended.txs.gr, ignore.strand = FALSE)),
             tx.strand = as.character(strand(extended.txs.gr)),
             mc.cores = arg.cores, SIMPLIFY = FALSE)
)

cat("Exporting transcripts for downstream sites...\n")
extended.annotation.gff3.file <- sprintf("data/gff/%s.txs.extutr3.%s.gff3", arg.outPrefix, arg.outSuffix)
extended.annotation.gtf.file <- sprintf("data/gff/%s.txs.extutr3.%s.gtf", arg.outPrefix, arg.outSuffix)
export.gff3(extended.annotation.gr, extended.annotation.gff3.file)
export(extended.annotation.gr, extended.annotation.gtf.file, format = "gtf")

cat("All operations complete.")
