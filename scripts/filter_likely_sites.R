#!/usr/bin/env Rscript

library(plyranges)

################################################################################
                                        # Fake Arguments
################################################################################

if (interactive()) {
    Snakemake <- setClass("Snakemake", slots=c(input='list', output='list',
                                               wildcards='list',
                                               params='list', threads='numeric'))
    snakemake <- Snakemake(
        input=list(unsupported="data/bed/cleavage-sites/utrome.unsupported.e3.t200.gc39.pas3.bed.gz",
                   passing="data/bed/cleavage-sites/utrome.cleavage.e3.t200.f0.9999.bed.gz"),
        output=list(supported="data/bed/cleavage-sites/utrome.supported.e3.t200.gc39.pas3.bed.gz",
                    unsupported="data/bed/cleavage-sites/utrome.unsupported.e3.t200.gc39.pas3.bed.gz"),
        wildcards=list(epsilon="3", threshold="200", version="39",
                       tpm="3", likelihood="0.9999"),
        params=list(),
        threads=1
        )
}

################################################################################
                                        # Load Data
################################################################################

gr_unsupported <- read_bed(snakemake@input$unsupported)
gr_passing <- read_bed(snakemake@input$passing)    

################################################################################
                                        # Export Sites
################################################################################

gr_unsupported %>%
    filter(name %in% gr_passing$name) %>%
    write_bed(file=snakemake@output$likely)

gr_unsupported %>%
    filter(!(name %in% gr_passing$name)) %>%
    write_bed(file=snakemake@output$unlikely)
