context('Test summarize module')

test_that('Test error on incorrect input class in summarize_score()', {
  bed = system.file('extdata', 'Gm12878_Ezh2_sorted_scores.narrowPeak.gz', package = 'annotatr')

  expect_error(summarize_score(bed), 'must have class tbl_df')
})

test_that('Test error on incorrect input class in summarize_name()', {
  bed = system.file('extdata', 'Gm12878_Ezh2_sorted_scores.narrowPeak.gz', package = 'annotatr')

  expect_error(summarize_name(bed), 'must have class tbl_df')
})

test_that('Test summarize_annotations()', {
  bed = system.file('extdata', 'Gm12878_Ezh2_sorted_scores.narrowPeak.gz', package = 'annotatr')
  annotations = c('cpgs')

  d = read_bed(filename = bed, genome = 'hg19', stranded = F, use.score = TRUE)

  i = annotate_regions(
    regions = d,
    annotations = annotations,
    genome = 'hg19',
    ignore.strand = T,
    use.score = T)

  s = summarize_annotation(i)

  expect_equal( sum(s[['n']]), expected = 2843)
})

test_that('Test summarize_score()', {
  bed = system.file('extdata', 'Gm12878_Ezh2_sorted_scores.narrowPeak.gz', package = 'annotatr')
  annotations = c('basic_genes','cpgs')

  d = read_bed(filename = bed, genome = 'hg19', stranded = F, use.score = TRUE)

  i = annotate_regions(
    regions = d,
    annotations = annotations,
    genome = 'hg19',
    ignore.strand = T,
    use.score = T)

  s = summarize_score(i)

  expect_equal( mean(s[['mean']]), expected = 25.19482, tolerance = 0.01)
})

test_that('Test summarize_name()', {
  bed = system.file('extdata', 'IDH2mut_v_NBM_names_scores_chr9.txt.gz', package = 'annotatr')
  annotations = c('basic_genes','cpgs')

  d = read_bed(filename = bed, genome = 'hg19', stranded = F, use.score = TRUE)

  i = annotate_regions(
    regions = d,
    annotations = annotations,
    genome = 'hg19',
    ignore.strand = T,
    use.score = TRUE)

  s = summarize_name(i)

  expect_equal( sum(s[['n']]), expected = 19984)
})
