#!/usr/bin/env Rscript

library(BSgenome.Hsapiens.UCSC.hg38)
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
        input=list(bed_all="data/bed/cleavage-sites/utrome.cleavage.e3.t200.bed.gz",
                   gencode="data/gff/gencode.v39.mRNA_ends_found.gff3.gz"),
        output=list(validated="data/bed/cleavage-sites/utrome.validated.e3.t200.gc39.bed.gz",
                    unvalidated="data/bed/cleavage-sites/utrome.unvalidated.e3.t200.gc39.bed.gz"),
        wildcards=list(epsilon="3", threshold="200", version="39"),
        params=list(radius="20"),
        threads=1
        )
}

################################################################################
                                        # Load Data
################################################################################

## Load References
hg38 <- BSgenome.Hsapiens.UCSC.hg38

## import transcript ends
gr_gencode <- import(snakemake@input$gencode, genome='hg38',
                     feature.type="transcript", colnames=character()) %>%
    keepStandardChromosomes(pruning.mode="coarse") %>%
    anchor_3p() %>% mutate(width=0)

## Load data
gr_sites_all <- import(snakemake@input$bed_all) %>%
    `seqlevelsStyle<-`("UCSC") %>%
    keepStandardChromosomes(pruning.mode="coarse")
seqlevels(gr_sites_all) <- seqlevels(hg38)
seqinfo(gr_sites_all) <- seqinfo(hg38)
names(gr_sites_all) <- gr_sites_all$name

################################################################################
                                        # Filter Sites
################################################################################

gr_validated <- (gr_sites_all + as.integer(snakemake@params$radius)) %>%
    find_overlaps_directed(gr_gencode) %>%
    unique()

################################################################################
                                        # Export Sites
################################################################################

gr_sites_all %>%
    filter(name %in% gr_validated$name) %>%
    write_bed(file=snakemake@output$validated)

gr_sites_all %>%
    filter(!(name %in% gr_validated$name)) %>%
    write_bed(file=snakemake@output$unvalidated)
