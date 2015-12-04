context('Test annotation functions')

test_that('Test tabulate_intersections().',{
  bed = system.file('extdata', 'Gm12878_Pol2.narrowPeak.gz', package = 'annotatr')
  annotations = c('basic_genes','cpgs')

  d = read_bed(filename = bed, genome = 'hg19', stranded = F)

  i = intersect_annotations(
    regions = d,
    annotations = annotations,
    genome = 'hg19',
    ignore.strand = T)

  t = tabulate_intersections(
    regions = d,
    intersections = i)

  # Get number of annotations from i and t and make sure they're the same
  i_lengths = sort(sapply(i, length))
  t_lengths = sort(table(t$annot_type))

  expect_equal( all(i_lengths == t_lengths) , expected = TRUE)
})
