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
        input=list(plus="data/granges/utrome.raw.plus.e3.t200.gc39.pas3.f0.9999.w500.Rds",
                   minus="data/granges/utrome.raw.minus.e3.t200.gc39.pas3.f0.9999.w500.Rds"),
        output=list(gtf="data/gff/utrome.e3.t200.gc39.pas3.f0.9999.w500.gtf",
                    fa="data/gff/utrome.e3.t200.gc39.pas3.f0.9999.w500.fa.gz"),
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

## set parallelization defaults
register(MulticoreParam(as.integer(snakemake@threads)))

## Load References
hg38 <- BSgenome.Hsapiens.UCSC.hg38

## Load data
message("Loading truncated transcripts...")
message("  strand: +")
gr_plus <- readRDS(snakemake@input$plus) %>%
    `names<-`(NULL)

message("  strand: -")
gr_minus <- readRDS(snakemake@input$minus) %>%
    `names<-`(NULL)

gr_txs <- filter(c(gr_plus, gr_minus), type == 'transcript') %>%
    `names<-`(.$ID)

grl_exons <- filter(c(gr_plus, gr_minus), type == 'exon') %>%
    expand_ranges(Parent) %>%
    split(.$Parent) %>%
    as("GRangesList")

## remove exact duplicates
df_intersect <- findOverlaps(grl_exons, minoverlap=MAX_TX_LENGTH, ignore.strand=FALSE, drop.self=TRUE) %>%
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


                                            
cat("Deduplicating transcripts...\n")
utrs.gr <- utrs.gr[!duplicated(utrs.gr),]

if (arg.txLength > 0) {
    ## Remove transcripts with less than 50 nt difference
    utrs.intersect <- gr.findoverlaps(utrs.gr, utrs.gr, minoverlap = arg.txLength - 50, by = 'gene_id', ignore.strand = FALSE) %Q% (query.id < subject.id)
    utrs.intersect$nearby.start <- abs(start(utrs.gr[utrs.intersect$query.id]) - start(utrs.gr[utrs.intersect$subject.id])) < 50 
    utrs.intersect$nearby.end <- abs(end(utrs.gr[utrs.intersect$query.id]) - end(utrs.gr[utrs.intersect$subject.id])) < 50 
    utrs.intersect <- utrs.intersect[utrs.intersect$nearby.start & utrs.intersect$nearby.end]

    is.augmented <- function (gr) {
        str_detect(gr$transcript_id, '(-|\\+)[0-9]+')
    }

    ## First, remove all de novo sites that are near GENCODE sites
    utrs.intersect$augmented.subject <- is.augmented(utrs.gr[utrs.intersect$subject.id])
    utrs.intersect$augmented.query <- is.augmented(utrs.gr[utrs.intersect$query.id])
    utrs.removable <- utrs.intersect %Q% xor(augmented.subject, augmented.query)
    
    ids.toRemove <- c(utrs.removable$subject.id[utrs.removable$augmented.subject],
                      utrs.removable$query.id[utrs.removable$augmented.query])
    
    utrs.intersect <- utrs.intersect[!(utrs.intersect$subject.id %in% ids.toRemove | utrs.intersect$query.id %in% ids.toRemove)]
    
    ## Second, for de novo sites, favor the one with smaller augmentation size
    utrs.denovo <- utrs.intersect %Q% (augmented.subject & augmented.query)
    utrs.denovo$pos.subject <- as.numeric(str_replace(utrs.gr$transcript_id[utrs.denovo$subject.id], ".*(-|\\+)(\\d+)$", "\\2"))
    utrs.denovo$pos.query <- as.numeric(str_replace(utrs.gr$transcript_id[utrs.denovo$query.id], ".*(-|\\+)(\\d+)$", "\\2"))
    
    ids.toRemove <- unique(c(ids.toRemove,
                             utrs.denovo$subject.id[utrs.denovo$pos.subject > utrs.denovo$pos.query],
                             utrs.denovo$query.id[utrs.denovo$pos.subject < utrs.denovo$pos.query]))
    
    utrs.intersect <- utrs.intersect[!(utrs.intersect$subject.id %in% ids.toRemove | utrs.intersect$query.id %in% ids.toRemove)]
    
    ## Third, for GENCODE-GENCODE, favor HAVANA (manually annotated) over ENSEMBL (automated pipeline)
    utrs.sources <- data.frame(query.id = utrs.intersect$query.id, subject.id = utrs.intersect$subject.id,
                               query.source = utrs.gr$source[utrs.intersect$query.id], subject.source = utrs.gr$source[utrs.intersect$subject.id])
    ids.toRemove <- unique(c(ids.toRemove,
                             utrs.sources$query.id[utrs.sources$subject.source == "HAVANA" & utrs.sources$query.source == "ENSEMBL"],
                             utrs.sources$subject.id[utrs.sources$query.source == "HAVANA" & utrs.sources$subject.source == "ENSEMBL"]))
    
    utrs.intersect <- utrs.intersect[!(utrs.intersect$subject.id %in% ids.toRemove | utrs.intersect$query.id %in% ids.toRemove)]
    
    ## Fourth, for GENCODE-GENCODE, favor more downstream end sites
    utrs.intersect$match.end <- end(gr.end(utrs.gr[utrs.intersect$subject.id], ignore.strand=FALSE)) == end(gr.end(utrs.gr[utrs.intersect$query.id], ignore.strand=FALSE))
    utrs.intersect.unmatched.pos <- utrs.intersect %Q% (!match.end & strand == '+')
    utrs.intersect.unmatched.neg <- utrs.intersect %Q% (!match.end & strand == '-')
    
    ids.toRemove <- unique(c(ids.toRemove,
                             ifelse(end(utrs.gr[utrs.intersect.unmatched.pos$query.id]) > end(utrs.gr[utrs.intersect.unmatched.pos$subject.id]),
                                    utrs.intersect.unmatched.pos$subject.id, utrs.intersect.unmatched.pos$query.id),
                             ifelse(start(utrs.gr[utrs.intersect.unmatched.neg$query.id]) < start(utrs.gr[utrs.intersect.unmatched.neg$subject.id]),
                                    utrs.intersect.unmatched.neg$subject.id, utrs.intersect.unmatched.neg$query.id)))
    
    utrs.intersect <- utrs.intersect[!(utrs.intersect$subject.id %in% ids.toRemove | utrs.intersect$query.id %in% ids.toRemove)]
    
    ## Finally, for GENCODE-GENCODE with identical spans, keep the transcripts that were discovered earlier
    ids.toRemove <- unique(c(ids.toRemove,
                             ifelse(utrs.gr$transcript_name[utrs.intersect$subject.id] < utrs.gr$transcript_name[utrs.intersect$query.id], utrs.intersect$query.id, utrs.intersect$subject.id)))
    
    stopifnot(length(utrs.intersect[!(utrs.intersect$subject.id %in% ids.toRemove | utrs.intersect$query.id %in% ids.toRemove)]) == 0)
    
    ## Remove all high overlapping transcripts from utrs.gr
    utrs.gr <- utrs.gr[-ids.toRemove]
}

## Generate UTR names by UTR position in gene
utrs.dt <- data.table(start = start(utrs.gr), end = end(utrs.gr), strand = as.character(strand(utrs.gr)), gene = utrs.gr$gene_name, tx_id = utrs.gr$transcript_id)
utrs.dt[, utr.name := ifelse(strand == "+", paste0(gene, '.', order(order(end, decreasing=FALSE))), paste0(gene, '.', order(order(start, decreasing=TRUE)))), by=gene]
setkey(utrs.dt, tx_id)

## Update child entries
cat("Intersecting truncated transcripts with child exons...\n")
children.utrs.dt <- gr.findoverlaps(
    utrs.gr, exons.gr,
    qcol = c("ID"), scol = c("Parent"), ignore.strand = FALSE,
    return.type = 'data.table', mc.cores = arg.cores)
children.utrs.dt <- children.utrs.dt[ID == Parent]

clipChild <- function (tx.start, tx.end, subject.id) {
    cat("\tExtracting exons...", '\n')
    child.gr <- exons.gr[subject.id, ]
    cat("\tUpdating start positions...", '\n')
    start(child.gr) <- pmax.int(start(child.gr), tx.start)
    cat("\tUpdating end positions...", '\n')
    end(child.gr) <- pmin.int(end(child.gr), tx.end)
    child.gr
}

cat("Adjusting exons to truncation bounds...\n")
exons.truncated.gr <- children.utrs.dt[, clipChild(start, end, subject.id)]

## Update genes based on new transcripts
cat("Intersecting genes with truncated transcripts...\n")
genes.txs.dt <- gr.findoverlaps(
    utrs.gr, genes.gr,
    scol = c("ID"), qcol = c("Parent"), ignore.strand = FALSE,
    return.type = 'data.table', mc.cores = arg.cores)
genes.txs.dt <- genes.txs.dt[ID == Parent]

cat("Adjusting gene bounds to match transcript ranges...\n")
## TODO: Improve speed of this method
gene.positions <- genes.txs.dt[, .(start = min(start), end = max(end)), by = subject.id]

genes.truncated.gr <- genes.gr[gene.positions$subject.id, ]
start(genes.truncated.gr) <- gene.positions$start
end(genes.truncated.gr) <- gene.positions$end

## Combine hierarchy of elements
utrome.gr <- grbind(genes.truncated.gr, utrs.gr, exons.truncated.gr)

cat("Creating TxDb object from UTRome annotation...\n")
## TODO: Figure out why some transcripts get dropped!
## TODO: Add more metadata
utrome.txdb <- makeTxDbFromGRanges(utrome.gr, taxonomyId = 10090,
                                   metadata = data.frame(name = c("Genome"), value = c("mm10")))

## Extract exons
cat("Extracting exons from TxDb...\n")
utrome.exons <- exons(utrome.txdb, columns = c('gene_id', 'tx_name', 'exon_name'))
mcols(utrome.exons)$type <- 'exon'
mcols(utrome.exons)$source <- 'GENCODE.vM21_MouseCellAtlas'
mcols(utrome.exons)$transcript_id <- as(utrs.dt[as.character(utrome.exons$tx_name), utr.name], "CharacterList")
mcols(utrome.exons)$transcript_name <- mcols(utrome.exons)$tx_name
mcols(utrome.exons)$exon_id <- as(mcols(utrome.exons)$exon_name, "CharacterList")
mcols(utrome.exons)$tx_id <- NULL
mcols(utrome.exons)$tx_name <- NULL

## Extract transcripts
cat("Extracting transcripts from TxDb...\n")
utrome.txs <- transcripts(utrome.txdb, columns = c('gene_id', 'tx_id', 'tx_name'))
mcols(utrome.txs)$type <- 'transcript'
mcols(utrome.txs)$source <- 'GENCODE.vM21_MouseCellAtlas'
mcols(utrome.txs)$transcript_id <- as(utrs.dt[as.character(utrome.txs$tx_name), utr.name], "CharacterList")
mcols(utrome.txs)$transcript_name <- mcols(utrome.txs)$tx_name
mcols(utrome.txs)$tx_id <- NULL
mcols(utrome.txs)$tx_name <- NULL

## Extract genes
cat("Extracting genes from TxDb...\n")
utrome.genes <- genes(utrome.txdb, columns = c('gene_id'))
mcols(utrome.genes)$type <- 'gene'
mcols(utrome.genes)$source <- 'GENCODE.vM21_MouseCellAtlas'


#########################
##  Export Annotation  ##
#########################
cat("Exporting UTRome GTF...\n")
utrome.gtf.file <- sprintf("data/gff/%s.utrome.%s.gtf", arg.outPrefix, arg.outSuffix)
export(utrome.genes, utrome.gtf.file, format = "GTF")
export(utrome.txs, utrome.gtf.file, format = "GTF", append = TRUE)
export(utrome.exons %Q% (strand == '+'), utrome.gtf.file, format = "GTF", append = TRUE)
export(sort(utrome.exons %Q% (strand == '-'), decreasing = TRUE), utrome.gtf.file, format = "GTF", append = TRUE)


########################
##  Export Sequences  ##
########################
cat("Exporting UTRome FASTA...\n")
utrome.seq.file <- sprintf("data/gff/%s.utrome.%s.fasta", arg.outPrefix, arg.outSuffix)
utrome.seq <- extractTranscriptSeqs(mm10, utrome.txdb, use.names=TRUE)
names(utrome.seq) <- utrs.dt[names(utrome.seq), utr.name]
writeXStringSet(utrome.seq, utrome.seq.file, format = "fasta")

cat("All operations complete.")
