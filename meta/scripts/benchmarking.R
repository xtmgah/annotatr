library(ChIPpeakAnno)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
devtools::load_all()

################################################################################
# ChIPpeakAnno

  ################################
  # ~365K rows
  cpa_small <- function(){
      annoData <- annoGR(TxDb.Hsapiens.UCSC.hg19.knownGene, feature="gene")
      file = '../data/IDH2mut_v_NBM_DM_test.txt'
      myPeakList = toGRanges(file, format="BED", header=FALSE)

      annotatedPeak <- annotatePeakInBatch(myPeakList,
                                           AnnotationData=annoData,
                                           output="overlapping")
  }

  ################################
  # ~4M rows
  cpa_medium <- function(){
      annoData <- annoGR(TxDb.Hsapiens.UCSC.hg19.knownGene, feature="gene")
      file = '../data/2607_mc_hmc_perc_meth_test2.txt'
      myPeakList = toGRanges(file, format="BED", header=FALSE)
      annotatedPeak <- annotatePeakInBatch(myPeakList,
                                           AnnotationData=annoData,
                                           output="overlapping")
  }
  ################################
  # ~25.4M rows
  cpa_large <- function(){
      annoData <- annoGR(TxDb.Hsapiens.UCSC.hg19.knownGene, feature="gene")
      file = '../data/IDH2mut_v_NBM_class_comp_trim.bed'
      myPeakList = toGRanges(file, format="BED", header=FALSE)
      annotatedPeak <- annotatePeakInBatch(myPeakList,
                                           AnnotationData=annoData,
                                           output="overlapping")
  }

################################################################################
# annotatr

  ################################
  # ~365K rows
  annotatr_small <- function(){
      file = '../data/IDH2mut_v_NBM_DM_test.txt'
      d = read_bed(file = file, col_names = FALSE, genome = 'hg19', stranded = FALSE, use.score = TRUE)
      annotations = c('cpgs', 'detailed_genes')
      t = annotate_regions(
          regions = d,
          annotations = annotations,
          genome = 'hg19',
          ignore.strand = T,
          use.score = T)
  }

  ################################
  # ~4M rows
  annotatr_medium <- function(){
      file = '../data/2607_mc_hmc_perc_meth_test2.txt'
      d = read_bed(file = file, col_names = FALSE, genome = 'hg19', stranded = TRUE, use.score = TRUE)
      annotations = c('cpgs', 'detailed_genes')
      t = annotate_regions(
          regions = d,
          annotations = annotations,
          genome = 'hg19',
          ignore.strand = T,
          use.score = T)
  }

  ################################
  # ~25M rows
  annotatr_large() <- function(){
      file = '../data/IDH2mut_v_NBM_class_comp_trim.bed'
      d = read_bed(file = file, col_names = FALSE, genome = 'hg19', stranded = TRUE, use.score = TRUE)
      annotations = c('cpgs', 'detailed_genes')
      t = annotate_regions(
          regions = d,
          annotations = annotations,
          genome = 'hg19',
          ignore.strand = T,
          use.score = T)
  }

test_small = microbenchmark(cpa_small(), annotatr_small(), times = 10)
test_med = microbenchmark(cpa_medium(), annotatr_medium(), times = 10)
test_large = microbenchmark(cpa_large(), annotatr_large(), times = 10)

write.table(test_small, file='../analysis/microbench_small_results.txt', sep='\t', quote=F)
write.table(test_med, file='../analysis/microbench_medium_results.txt', sep='\t', quote=F)
write.table(test_large, file='../analysis/microbench_large_results.txt', sep='\t', quote=F)
