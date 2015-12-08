library(ChIPpeakAnno)

################################################################################
# ~164K rows, ~14s
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
annoData <- annoGR(TxDb.Hsapiens.UCSC.hg19.knownGene, feature="gene")
file = '../data/IDH2mut_v_NBM_test.txt'
myPeakList = toGRanges(file, format="BED", header=FALSE)

annotatedPeak <- annotatePeakInBatch(myPeakList,
                                     AnnotationData=annoData,
                                     output="overlapping")


################################################################################
# ~4M rows, ~5:15m
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
annoData <- annoGR(TxDb.Hsapiens.UCSC.hg19.knownGene, feature="gene")
file = '../data/2607_mc_hmc_perc_meth_test2.txt'
myPeakList = toGRanges(file, format="BED", header=FALSE)
annotatedPeak <- annotatePeakInBatch(myPeakList,
                                     AnnotationData=annoData,
                                     output="overlapping")

################################################################################
# ~25.4M rows ...
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
annoData <- annoGR(TxDb.Hsapiens.UCSC.hg19.knownGene, feature="gene")
file = '../data/IDH2mut_v_NBM_class_comp_trim.bed'
myPeakList = toGRanges(file, format="BED", header=FALSE)
annotatedPeak <- annotatePeakInBatch(myPeakList,
                                     AnnotationData=annoData,
                                     output="overlapping")
