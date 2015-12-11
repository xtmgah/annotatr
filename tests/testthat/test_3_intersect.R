context('Test intersection and annotation of intersection')

test_that('Test error thrown for non-GRanges regions object in annotate_regions()',{
  bed = system.file('extdata', 'test_intersect.bed', package = 'annotatr')
  annotations = c('hg19_cpg_islands')
  d = read_bed(filename = bed, genome = 'hg19', stranded = F)

  expect_error(
      annotate_regions(
      regions = bed,
      annotations = annotations,
      genome = 'hg20',
      ignore.strand = T,
      use.score = F),
    "GRanges")
})

test_that('Test error thrown for unsupported genome in annotate_regions()',{
  bed = system.file('extdata', 'test_intersect.bed', package = 'annotatr')
  annotations = c('hg19_cpg_islands')
  d = read_bed(filename = bed, genome = 'hg19', stranded = F)

  expect_error(
    annotate_regions(
      regions = d,
      annotations = annotations,
      genome = 'hg20',
      ignore.strand = T,
      use.score = F),
    "unsupported genome")
})

test_that('Test error thrown for unsupported annotation in annotate_regions()',{
  bed = system.file('extdata', 'test_intersect.bed', package = 'annotatr')
  annotations = c('hg19_cpg_islands','big_willy_style')
  d = read_bed(filename = bed, genome = 'hg19', stranded = F)

  expect_error(
    annotate_regions(
      regions = d,
      annotations = annotations,
      genome = 'hg19',
      ignore.strand = T,
      use.score = F),
    "big_willy_style")
})

test_that('Test error thrown for basic_genes and detailed_genes in annotate_regions()',{
  bed = system.file('extdata', 'test_intersect.bed', package = 'annotatr')
  annotations = c('basic_genes', 'detailed_genes')
  d = read_bed(filename = bed, genome = 'hg19', stranded = F)

  expect_error(
    annotate_regions(
      regions = d,
      annotations = annotations,
      genome = 'hg19',
      ignore.strand = T,
      use.score = F),
    "please choose between")
})

test_that('Test a la carte annotations in annotate_regions()',{
  bed = system.file('extdata', 'test_intersect.bed', package = 'annotatr')
  annotations = c('hg19_cpg_islands', 'hg19_knownGenes_promoters')
  d = read_bed(filename = bed, genome = 'hg19', stranded = F)

  i = annotate_regions(
    regions = d,
    annotations = annotations,
    genome = 'hg19',
    ignore.strand = T,
    use.score = F)

  expect_true( dplyr::setequal(unique(i[['annot_type']]), annotations) )
})

test_that('Test a la carte and shortcut annotations in annotate_regions()',{
  bed = system.file('extdata', 'Gm12878_Ezh2_sorted_scores.narrowPeak.gz', package = 'annotatr')
  annotations = c('hg19_cpg_islands', 'basic_genes')
  d = read_bed(filename = bed, genome = 'hg19', stranded = F)

  i = annotate_regions(
    regions = d,
    annotations = annotations,
    genome = 'hg19',
    ignore.strand = T,
    use.score = F)

  expect_true( dplyr::setequal(unique(i[['annot_type']]), c('hg19_cpg_islands', 'hg19_knownGenes_1to5kb', 'hg19_knownGenes_promoters', 'hg19_knownGenes_5UTRs', 'hg19_knownGenes_exons', 'hg19_knownGenes_introns', 'hg19_knownGenes_3UTRs')) )
})

test_that('Test dual annotation shortcut in annotate_regions()',{
  bed = system.file('extdata', 'Gm12878_Ezh2_sorted_scores.narrowPeak.gz', package = 'annotatr')
  annotations = c('basic_genes','cpgs')

  d = read_bed(filename = bed, genome = 'hg19', stranded = F)

  i = annotate_regions(
    regions = d,
    annotations = annotations,
    genome = 'hg19',
    ignore.strand = TRUE,
    use.score = F)

  expect_equal( nrow(i) , expected = 7216)
})
