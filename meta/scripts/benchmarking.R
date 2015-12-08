library(ChIPpeakAnno)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
devtools::load_all()

################################################################################
# ~164K rows, ~14s
cpa_small <- function(){

    annoData <- annoGR(TxDb.Hsapiens.UCSC.hg19.knownGene, feature="gene")
    file = '../data/IDH2mut_v_NBM_test.txt'
    myPeakList = toGRanges(file, format="BED", header=FALSE)

    annotatedPeak <- annotatePeakInBatch(myPeakList,
                                         AnnotationData=annoData,
                                         output="overlapping")
}

################################################################################
# ~4M rows, ~5:15m

cpa_medium <- function(){

    annoData <- annoGR(TxDb.Hsapiens.UCSC.hg19.knownGene, feature="gene")
    file = '../data/2607_mc_hmc_perc_meth_test2.txt'
    myPeakList = toGRanges(file, format="BED", header=FALSE)
    annotatedPeak <- annotatePeakInBatch(myPeakList,
                                         AnnotationData=annoData,
                                         output="overlapping")
}
################################################################################
# ~25.4M rows ...
cpa_large <- function(){
    library(TxDb.Hsapiens.UCSC.hg19.knownGene)
    annoData <- annoGR(TxDb.Hsapiens.UCSC.hg19.knownGene, feature="gene")
    file = '../data/IDH2mut_v_NBM_class_comp_trim.bed'
    myPeakList = toGRanges(file, format="BED", header=FALSE)
    annotatedPeak <- annotatePeakInBatch(myPeakList,
                                         AnnotationData=annoData,
                                         output="overlapping")
}

dplyr_small <- function(){
    file = '../data/IDH2mut_v_NBM_test.txt'
    d = read_bed(filename = file, genome = 'hg19', stranded = FALSE, use.score = TRUE)
    annotations = c('cpgs', 'detailed_genes')
    i = intersect_annotations(
        regions = d,
        annotations = annotations,
        genome = 'hg19',
        ignore.strand = T)
    t = annotate_intersections(
        regions = d,
        intersections = i,
        use.score = T)
}

dplyr_medium <- function(){
    file = '../data/2607_mc_hmc_perc_meth_test2.txt'
    d = read_bed(filename = file, genome = 'hg19', stranded = TRUE, use.score = TRUE)
    annotations = c('cpgs', 'detailed_genes')
    i = intersect_annotations(
        regions = d,
        annotations = annotations,
        genome = 'hg19',
        ignore.strand = F)
    t = annotate_intersections(
        regions = d,
        intersections = i,
        use.score = T)

}

test_small <- microbenchmark(cpa_small(), dplyr_small(), times = 10)
test_med <- microbenchmark(cpa_medium(), dplyr_medium(), times = 10)

write.table(test_small, file='../analysis/microbench_small_results.txt', sep='\t')
write.table(test_med, file='../analysis/microbench_medium_results.txt', sep='\t')


