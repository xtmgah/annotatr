context('Test summarize module')

################################################################################
# Test errors

bed = system.file('extdata', 'Gm12878_Ezh2_sorted_scores.narrowPeak.gz', package = 'annotatr')

test_that('Test error on incorrect input class in summarize functions', {
  expect_error(summarize_annotations(annotated_regions = bed), 'must have class tbl_df')
  expect_error(summarize_numerical(annotated_regions = bed), 'must have class tbl_df')
  expect_error(summarize_categorical(annotated_regions = bed), 'must have class tbl_df')
})

################################################################################
# Test summarize functions

test_that('Test summarize_annotations()', {
  bed = system.file('extdata', 'Gm12878_Ezh2_sorted_scores.narrowPeak.gz', package = 'annotatr')
  annotations = c('hg19_cpgs')

  d = read_bed(file = bed, genome = 'hg19', stranded = F, use.score = TRUE)

  i = annotate_regions(
    regions = d,
    annotations = annotations,
    ignore.strand = T,
    use.score = T)

  s = summarize_annotations(i)

  expect_equal( sum(s[['n']]), expected = 2843)
})

test_that('Test summarize_numerical()', {
  bed = system.file('extdata', 'Gm12878_Ezh2_sorted_scores.narrowPeak.gz', package = 'annotatr')
  annotations = c('hg19_basicgenes','hg19_cpgs')

  d = read_bed(file = bed, genome = 'hg19', stranded = F, use.score = TRUE)

  i = annotate_regions(
    regions = d,
    annotations = annotations,
    ignore.strand = T,
    use.score = T)

  s = summarize_numerical(
    annotated_regions = i,
    by = c('annot_type', 'annot_id'))

  expect_equal( mean(s[['mean']]), expected = 25.19482, tolerance = 0.01)
})

test_that('Test summarize_numerical() and summarize_categorical() over small data', {
  bed = system.file('extdata', 'test_read_multiple_data_head.bed', package = 'annotatr')
  annotations = c('hg19_cpg_islands', 'hg19_knownGenes_promoters')

  d = read_bed(
    file = bed,
    col.names=TRUE,
    genome = 'hg19',
    stranded = FALSE,
    use.score = TRUE)

  i = annotate_regions(
    regions = d,
    annotations = annotations,
    ignore.strand = T,
    use.score = TRUE)

  # Testing summarize_numerical()
    sn1 = summarize_numerical(
      annotated_regions = i,
      by = c('annot_type', 'annot_id'))
    sn2 = summarize_numerical(
      annotated_regions = i,
      by = c('name'),
      over = c('score', 'mu1', 'mu0'))

  # Testing summarize_categorical()
    sc1 = summarize_categorical(
      annotated_regions = i,
      by = c('annot_type', 'name'))
    sc2 = summarize_categorical(
      annotated_regions = i,
      by = c('annot_type', 'diff_exp'))
    sc3 = summarize_categorical(
      annotated_regions = i,
      by = c('name','diff_exp'))

  expect_equal( sn1[[which(sn1[,'annot_id'] == 'uc010nxq.1'), 'mean']], expected = 66)
  expect_equal( sn1[[which(sn1[,'annot_id'] == 'island:1'), 'mean']], expected = 48)
  expect_equal( sn2[[which(sn2[,'name'] == 'A'), 'mu0_mean']], expected = 25)
  expect_equal( sn2[[which(sn2[,'name'] == 'B'), 'mu1_mean']], expected = 95)

  expect_equal( sc1[[which(sc1[,'annot_type'] == 'hg19_knownGenes_promoters' & sc1[,'name'] == 'A'), 'n']], expected = 2)
  expect_equal( sc2[[which(sc2[,'annot_type'] == 'hg19_cpg_islands' & sc2[,'diff_exp'] == 'Y'), 'n']], expected = 2)
  expect_equal( sc3[[which(sc3[,'name'] == 'A' & sc3[,'diff_exp'] == 'N'), 'n']], expected = 1)
})
