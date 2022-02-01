#!/usr/bin/env Rscript

library(plyranges)
library(BiocParallel)

################################################################################
                                        # Fake Arguments
################################################################################

if (interactive()) {
    Snakemake <- setClass("Snakemake", slots=c(input='list', output='list',
                                               wildcards='list',
                                               params='list', threads='numeric'))
    snakemake <- Snakemake(
        input=list(granges="data/granges/augmented.plus.chunked.e3.t200.gc39.pas3.f0.9999.Rds"),
        output=list(granges="/fscratch/fanslerm/utrome.raw.plus.Rds"),
        wildcards=list(epsilon="3", threshold="200", version="39", tpm="3",
                       likelihood="0.9999", width=500),
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

## Load GRangesList
message("Loading GRangesList...")
grl_plus <- readRDS(snakemake@input$granges)

################################################################################
                                        # Adjust Gene Sizes - Positive Strand
################################################################################

message("Adjusting transcript sizes...")
gr_plus_adjusted <- bplapply(grl_plus, function (gr) {
    # split by feature types
    gr_genes <- filter(gr, type == 'gene')
    gr_txs   <- filter(gr, type == 'transcript')
    gr_exons <- filter(gr, type == 'exon')

    ## Intersect exons with transcripts
    gr_children <- gr_txs %>%
        plyranges::select(ID) %>%
        join_overlap_intersect_directed(x=gr_exons, suffix=c("", ".tx")) %>%
        expand_ranges(Parent) %>%
        filter(Parent == ID.tx)

    ## clipping and removing exons
    gr_exons_adjusted <- gr_children %>%
        ## order from 3' ends
        sort(decreasing=TRUE) %>%
        ## group by tx
        group_by(Parent) %>%
        ## compute tx-relative position of exon
        mutate(tx_length_end=cumsum(width),
               tx_length_start=tx_length_end - width) %>%
        ## exclude exons beyond max tx length
        filter(tx_length_start < MAX_TX_LENGTH) %>%
        ## adjust all exons traversing cutoff
        mutate(start=ifelse(tx_length_end > MAX_TX_LENGTH, start + tx_length_end - MAX_TX_LENGTH, start)) %>%
        ungroup()

    ## check everything is as expected
    gr_exons_adjusted %>%
        group_by(Parent) %>%
        summarize(tx_length=sum(width)) %>%
        { stopifnot(all(.$tx_length <= MAX_TX_LENGTH)) }
    
    ## adjust transcript bounds
    df_txs_bounds <- gr_exons_adjusted %>%
        group_by(Parent) %>%
        summarize(start=min(start), end=max(end))
    
    gr_txs_adjusted <- gr_txs %>% `names<-`(.$ID)
    
    start(gr_txs_adjusted[df_txs_bounds$Parent]) <- df_txs_bounds$start
    end(gr_txs_adjusted[df_txs_bounds$Parent]) <- df_txs_bounds$end
    
    ## filter any missing IDs
    gr_txs_adjusted <- gr_txs_adjusted[gr_txs_adjusted$gene_id %in% gr_genes$ID,]
    gr_exons_adjusted <- gr_exons_adjusted[gr_exons_adjusted$gene_id %in% gr_genes$ID,]

    ## adjust gene bounds
    df_genes_bounds <- gr_txs_adjusted %>%
        expand_ranges(Parent) %>%
        group_by(Parent) %>%
        summarize(start=min(start), end=max(end))

    gr_genes_adjusted <- gr_genes %>%
        `names<-`(.$ID)

    start(gr_genes_adjusted[df_genes_bounds$Parent]) <- df_genes_bounds$start
    end(gr_genes_adjusted[df_genes_bounds$Parent]) <- df_genes_bounds$end
    
    ## merge elements
    gr_chunk <- bind_ranges(gr_genes_adjusted, gr_txs_adjusted, gr_exons_adjusted) %>%
        select(-c("transcript_id_old", "cleavage_site", "name", "offset", "ID.tx", "tx_length_end", "tx_length_start")) %>%
        `names<-`(NULL)

    gr_chunk
}) %>% as("GRangesList") %>% unlist()

message("Exporting clipped GRanges...")
saveRDS(gr_plus_adjusted, snakemake@output$granges)
