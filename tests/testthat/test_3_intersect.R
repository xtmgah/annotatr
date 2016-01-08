context('Test intersect/annotate module')

################################################################################
# Test errors

test_that('Test error thrown for non-GRanges regions object in annotate_regions()',{
  bed = system.file('extdata', 'test_intersect.bed', package = 'annotatr')
  annotations = c('hg19_cpg_islands')
  d = read_bed(file = bed, genome = 'hg19', stranded = FALSE)

  expect_error(
      annotate_regions(
      regions = bed,
      annotations = annotations,
      ignore.strand = TRUE,
      use.score = FALSE),
    "GRanges")
})

################################################################################
# Test annotate_regions()

test_that('Test a la carte annotations in annotate_regions()',{
  bed = system.file('extdata', 'test_intersect.bed', package = 'annotatr')
  annotations = c('hg19_cpg_islands', 'hg19_knownGenes_promoters')
  d = read_bed(file = bed, genome = 'hg19', stranded = FALSE)

  i = annotate_regions(
    regions = d,
    annotations = annotations,
    ignore.strand = TRUE,
    use.score = FALSE)

  expect_true( dplyr::setequal(unique(i[['annot_type']]), annotations) )
})

test_that('Test a la carte and shortcut annotations in annotate_regions()',{
  bed = system.file('extdata', 'Gm12878_Ezh2_sorted_scores.narrowPeak.gz', package = 'annotatr')
  annotations = c('hg19_cpg_islands', 'hg19_basicgenes')
  d = read_bed(file = bed, genome = 'hg19', stranded = FALSE)

  i = annotate_regions(
    regions = d,
    annotations = annotations,
    ignore.strand = TRUE,
    use.score = FALSE)

  expect_true( dplyr::setequal(unique(i[['annot_type']]), c('hg19_cpg_islands', 'hg19_knownGenes_1to5kb', 'hg19_knownGenes_promoters', 'hg19_knownGenes_5UTRs', 'hg19_knownGenes_exons', 'hg19_knownGenes_introns', 'hg19_knownGenes_3UTRs')) )
})

test_that('Test dual annotation shortcut in annotate_regions()',{
  bed = system.file('extdata', 'Gm12878_Ezh2_sorted_scores.narrowPeak.gz', package = 'annotatr')
  annotations = c('hg19_basicgenes','hg19_cpgs')

  d = read_bed(file = bed, genome = 'hg19', stranded = FALSE)

  i = annotate_regions(
    regions = d,
    annotations = annotations,
    ignore.strand = TRUE,
    use.score = FALSE)

  expect_equal( nrow(i) , expected = 7216)
})

test_that('Custom annotations work in annotate_regions()', {
  a_file = system.file('extdata', 'test_annotations_3.bed', package='annotatr')
  r_file = system.file('extdata', 'test_read_multiple_data_head.bed', package='annotatr')

  r = read_bed(file=r_file, col.names=TRUE, genome='hg19', stranded=FALSE, use.score=FALSE)
  hg19_custom_TFBS = read_annotations(file=a_file, genome='hg19', annotation_name='TFBS')

  annotations = c('hg19_custom_TFBS', 'hg19_cpgs')

  a = annotate_regions(
    regions = r,
    annotations = annotations,
    ignore.strand = TRUE,
    use.score = TRUE
  )

  expect_equal(nrow(a) == 10, expected = TRUE)
})
