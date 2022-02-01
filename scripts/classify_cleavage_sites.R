#!/usr/bin/env Rscript

library(BSgenome.Hsapiens.UCSC.hg38)
library(cleanUpdTSeq)
library(rtracklayer)
library(magrittr)
library(BiocParallel)

################################################################################
                                        # Fake Arguments
################################################################################

if (interactive()) {
    Snakemake <- setClass("Snakemake", slots=c(input='list', output='list', threads='numeric'))
    snakemake <- Snakemake(
        input=list(bed="data/bed/cleavage-sites/utrome.cleavage.e3.t200.bed.gz"),
        output=list(probs="data/bed/cleavage-sites/utrome.classification.e3.t200.tsv.gz"),
        threads=3
        )
}

################################################################################
                                        # Load Data
################################################################################

gr_peaks <- import(snakemake@input$bed) %>%
    `seqlevelsStyle<-`("UCSC") %>%
    keepStandardChromosomes(pruning.mode="coarse")
seqlevels(gr_peaks) <- seqlevels(BSgenome.Hsapiens.UCSC.hg38)
seqinfo(gr_peaks) <- seqinfo(BSgenome.Hsapiens.UCSC.hg38)

## Need row names or else cleanUpdTSeq fails
names(gr_peaks) <- elementMetadata(gr_peaks)$name

################################################################################
                                        # Filter Edge Cases
################################################################################
#' Peaks too close to the telomeres are not able to extract features needed to
#' run the cleanUpdTSeq procedure.

idx_edge_peaks <- gr_peaks %>%
    { start(.) < 40 |
          end(.) > (seqlengths(.)[as.character(seqnames(.))] - 40) } %>%
    set_names(names(gr_peaks)) %>%
    which

## precreate NA entries to be appended to results
res_na <- data.frame(peak_name=names(idx_edge_peaks),
                     prob_fake_pA=NA,
                     prob_true_pA=NA,
                     pred_class=NA)

print("The following peaks could not be classified:")
print(gr_peaks[idx_edge_peaks])

################################################################################
                                        # Classify Peaks
################################################################################

classify_peaks <- function (gr) {
    vec_features <- buildFeatureVector(gr, genome=Hsapiens, sampleType='unknown',
                                       upstream=40L, downstream=30L, wordSize=6L, method="NaiveBayes",
                                       alphabet=c("ACGT"), replaceNAdistance=30L, fetchSeq=TRUE)

    ## load default classifier
    data(classifier)

    predictTestSet(testSet.NaiveBayes=vec_features, outputFile=NULL,
                   classifier=classifier, assignmentCutoff=0.5)
}

res <- gr_peaks[-idx_edge_peaks] %>%
    ## split
    { split(., ceiling(seq_along(.)/1000)) } %>%

    ## apply
    bplapply(classify_peaks, BPPARAM=MulticoreParam(snakemake@threads)) %>%

    ## combine
    do.call(what=rbind) %>%

    ## include NAs
    rbind(res_na)

file_out <- gzfile(snakemake@output$probs, "w")
write.table(res, file_out, sep='\t', row.names=FALSE, quote=FALSE)
