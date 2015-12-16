context('Test summarize module')

################################################################################
# Test errors

bed = system.file('extdata', 'Gm12878_Ezh2_sorted_scores.narrowPeak.gz', package = 'annotatr')

test_that('Test error on incorrect input class in summarize functions', {
  expect_error(summarize_annotation(annotated_regions = bed), 'must have class tbl_df')
  expect_error(summarize_score(annotated_regions = bed), 'must have class tbl_df')
  expect_error(summarize_name(annotated_regions = bed), 'must have class tbl_df')
})

################################################################################
# Test summarize functions

test_that('Test summarize_annotations()', {
  bed = system.file('extdata', 'Gm12878_Ezh2_sorted_scores.narrowPeak.gz', package = 'annotatr')
  annotations = c('cpgs')

  d = read_bed(file = bed, genome = 'hg19', stranded = F, use.score = TRUE)

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

  d = read_bed(file = bed, genome = 'hg19', stranded = F, use.score = TRUE)

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

  d = read_bed(file = bed, genome = 'hg19', stranded = F, use.score = TRUE)

  i = annotate_regions(
    regions = d,
    annotations = annotations,
    genome = 'hg19',
    ignore.strand = T,
    use.score = TRUE)

  s = summarize_name(i)

  expect_equal( sum(s[['n']]), expected = 19984)
})
