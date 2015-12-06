context('Test intersection and annotation of intersection')

test_that('Test error thrown for non-GRanges regions object in intersect_annotations()',{
  bed = system.file('extdata', 'test_intersect.bed', package = 'annotatr')
  annotations = c('hg19_cpg_islands')
  d = read_bed(filename = bed, genome = 'hg19', stranded = F)

  expect_error(
    intersect_annotations(
    regions = bed,
    annotations = annotations,
    genome = 'hg20',
    ignore.strand = T,
    "GRanges")
  )
})

test_that('Test error thrown for unsupported genome in intersect_annotations()',{
  bed = system.file('extdata', 'test_intersect.bed', package = 'annotatr')
  annotations = c('hg19_cpg_islands')
  d = read_bed(filename = bed, genome = 'hg19', stranded = F)

  expect_error(
    intersect_annotations(
    regions = d,
    annotations = annotations,
    genome = 'hg20',
    ignore.strand = T,
    "unsupported genome")
  )
})

test_that('Test error thrown for unsupported annotation in intersect_annotations()',{
  bed = system.file('extdata', 'test_intersect.bed', package = 'annotatr')
  annotations = c('hg19_cpg_islands','big_willy_style')
  d = read_bed(filename = bed, genome = 'hg19', stranded = F)

  expect_error(
    intersect_annotations(
    regions = d,
    annotations = annotations,
    genome = 'hg19',
    ignore.strand = T,
    "big_willy_style")
  )
})

test_that('Test error thrown for basic_genes and detailed_genes in intersect_annotations()',{
  bed = system.file('extdata', 'test_intersect.bed', package = 'annotatr')
  annotations = c('basic_genes', 'detailed_genes')
  d = read_bed(filename = bed, genome = 'hg19', stranded = F)

  expect_error(
    intersect_annotations(
    regions = d,
    annotations = annotations,
    genome = 'hg19',
    ignore.strand = T,
    "please choose between")
  )
})

test_that('Test dual annotation shortcut in intersect_annotations()',{
  bed = system.file('extdata', 'Gm12878_Pol2.narrowPeak.gz', package = 'annotatr')
  annotations = c('basic_genes','cpgs')

  d = read_bed(filename = bed, genome = 'hg19', stranded = F)

  i = intersect_annotations(
    regions = d,
    annotations = annotations,
    genome = 'hg19',
    ignore.strand = T)

  expect_equal( length(i) , expected = 10)
})


test_that('Test correct intersections in intersect_annotations()',{
  bed = system.file('extdata', 'test_intersect.bed', package = 'annotatr')
  annotations = c('hg19_cpg_islands','hg19_cpg_shores','hg19_knownGenes_promoters')

  d = read_bed(filename = bed, genome = 'hg19', stranded = F)

  i = intersect_annotations(
    regions = d,
    annotations = annotations,
    genome = 'hg19',
    ignore.strand = T)

  query_hits = sapply(i, function(a){a@queryHits})

  expect_equal( all(query_hits == c(3,2,1)) , expected = TRUE)
})

test_that('Test annotate_intersections() works',{
  bed = system.file('extdata', 'Gm12878_Pol2.narrowPeak.gz', package = 'annotatr')
  annotations = c('basic_genes','cpgs')

  d = read_bed(filename = bed, genome = 'hg19', stranded = F, use.score = FALSE)

  i = intersect_annotations(
    regions = d,
    annotations = annotations,
    genome = 'hg19',
    ignore.strand = T)

  t = annotate_intersections(
    regions = d,
    intersections = i,
    use.score = FALSE)

  # Get number of annotations from i and t and make sure they're the same
  i_lengths = sort(sapply(i, length))
  t_lengths = sort(table(t$annot_type))

  expect_equal( all(i_lengths == t_lengths) , expected = TRUE)
})
