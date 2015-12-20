context('Test intersect/annotate module')

################################################################################
# Test errors

test_that('Test error thrown for non-GRanges regions object in annotate_regions()',{
  bed = system.file('extdata', 'test_intersect.bed', package = 'annotatr')
  annotations = c('hg19_cpg_islands')
  d = read_bed(file = bed, genome = 'hg19', stranded = F)

  expect_error(
      annotate_regions(
      regions = bed,
      annotations = annotations,
      ignore.strand = T,
      use.score = F),
    "GRanges")
})

################################################################################
# Test annotate_regions()

test_that('Test a la carte annotations in annotate_regions()',{
  bed = system.file('extdata', 'test_intersect.bed', package = 'annotatr')
  annotations = c('hg19_cpg_islands', 'hg19_knownGenes_promoters')
  d = read_bed(file = bed, genome = 'hg19', stranded = F)

  i = annotate_regions(
    regions = d,
    annotations = annotations,
    ignore.strand = T,
    use.score = F)

  expect_true( dplyr::setequal(unique(i[['annot_type']]), annotations) )
})

test_that('Test a la carte and shortcut annotations in annotate_regions()',{
  bed = system.file('extdata', 'Gm12878_Ezh2_sorted_scores.narrowPeak.gz', package = 'annotatr')
  annotations = c('hg19_cpg_islands', 'hg19_basicgenes')
  d = read_bed(file = bed, genome = 'hg19', stranded = F)

  i = annotate_regions(
    regions = d,
    annotations = annotations,
    ignore.strand = T,
    use.score = F)

  expect_true( dplyr::setequal(unique(i[['annot_type']]), c('hg19_cpg_islands', 'hg19_knownGenes_1to5kb', 'hg19_knownGenes_promoters', 'hg19_knownGenes_5UTRs', 'hg19_knownGenes_exons', 'hg19_knownGenes_introns', 'hg19_knownGenes_3UTRs')) )
})

test_that('Test dual annotation shortcut in annotate_regions()',{
  bed = system.file('extdata', 'Gm12878_Ezh2_sorted_scores.narrowPeak.gz', package = 'annotatr')
  annotations = c('hg19_basicgenes','hg19_cpgs')

  d = read_bed(file = bed, genome = 'hg19', stranded = F)

  i = annotate_regions(
    regions = d,
    annotations = annotations,
    ignore.strand = TRUE,
    use.score = F)

  expect_equal( nrow(i) , expected = 7216)
})
