context('Test visualization of summaries')

test_that('Test error on incorrect input class in visualize_score()', {
  bed = system.file('extdata', 'Gm12878_Ezh2_sorted_scores.narrowPeak.gz', package = 'annotatr')

  expect_error(visualize_score(bed), 'summarized_scores must have class grouped_df')
})

test_that('Test visualize_score() on signalValue from ChIP-seq data', {
  bed = system.file('extdata', 'Gm12878_Ezh2_sorted_scores.narrowPeak.gz', package = 'annotatr')
  annotations = c('basic_genes','cpgs')

  d = read_bed(filename = bed, genome = 'hg19', stranded = F, use.score = TRUE)

  i = intersect_annotations(
    regions = d,
    annotations = annotations,
    genome = 'hg19',
    ignore.strand = T)

  t = annotate_intersections(
    regions = d,
    intersections = i,
    use.score = TRUE)

  s = summarize_score(t)

  cpgs_order = c(
    'hg19_cpg_islands',
    'hg19_cpg_shores',
    'hg19_cpg_shelves',
    'hg19_cpg_inter')
  v_cpgs = visualize_score(s, cpgs_order)

  genes_order = c(
    'hg19_knownGenes_1to5kb',
    'hg19_knownGenes_promoters',
    'hg19_knownGenes_5UTRs',
    'hg19_knownGenes_exons',
    'hg19_knownGenes_introns',
    'hg19_knownGenes_3UTRs')
  v_genes = visualize_score(s, genes_order)

  cpgs_error_order = c(
    'hg19_cpg_islands',
    'hg19_cpg_shores',
    'hg19_cpg_shelves',
    'hg19_cpg_inter',
    'walla_walla')

  expect_error( visualize_score(s, cpgs_error_order), 'annotations in annotation_order that are not present')
  expect_equal( dplyr::setequal(class(v_cpgs), c('gg','ggplot')), expected = TRUE)
  expect_equal( dplyr::setequal(class(v_genes), c('gg','ggplot')), expected = TRUE)
})

test_that('Test visualize_score() on methylation data', {
  bed = system.file('extdata', '2607_mc_hmc_perc_meth_chr21.txt.gz', package = 'annotatr')
  annotations = c('basic_genes','cpgs')

  d = read_bed(filename = bed, genome = 'hg19', stranded = F, use.score = TRUE)

  i = intersect_annotations(
    regions = d,
    annotations = annotations,
    genome = 'hg19',
    ignore.strand = T)

  t = annotate_intersections(
    regions = d,
    intersections = i,
    use.score = TRUE)

  s = summarize_score(t)

  cpgs_order = c(
    'hg19_cpg_islands',
    'hg19_cpg_shores',
    'hg19_cpg_shelves',
    'hg19_cpg_inter')
  v_cpgs = visualize_score(s, cpgs_order)

  genes_order = c(
    'hg19_knownGenes_1to5kb',
    'hg19_knownGenes_promoters',
    'hg19_knownGenes_5UTRs',
    'hg19_knownGenes_exons',
    'hg19_knownGenes_introns',
    'hg19_knownGenes_3UTRs')
  v_genes = visualize_score(s, genes_order)

  expect_equal( dplyr::setequal(class(v_cpgs), c('gg','ggplot')), expected = TRUE)
  expect_equal( dplyr::setequal(class(v_genes), c('gg','ggplot')), expected = TRUE)
})

test_that('Test error on incorrect input class in visualize_name()', {
  bed = system.file('extdata', 'Gm12878_Ezh2_sorted_scores.narrowPeak.gz', package = 'annotatr')

  expect_error(visualize_name(bed), 'summarized_names must have class grouped_df')
})

test_that('Test visualize_name() on classification data', {
  bed = system.file('extdata', 'IDH2mut_v_NBM_names_scores_chr9.txt.gz', package = 'annotatr')
  annotations = c('basic_genes','cpgs')

  d = read_bed(filename = bed, genome = 'hg19', stranded = F, use.score = TRUE)

  i = intersect_annotations(
    regions = d,
    annotations = annotations,
    genome = 'hg19',
    ignore.strand = T)

  t = annotate_intersections(
    regions = d,
    intersections = i,
    use.score = TRUE)

  s = summarize_name(t)

  cpgs_order = c(
    'hg19_cpg_islands',
    'hg19_cpg_shores',
    'hg19_cpg_shelves',
    'hg19_cpg_inter')
  data_order = c(
    'DMup',
    'DMdown')
  v_cpgs_counts = visualize_name(s, cpgs_order, data_order, fill=FALSE)
  v_cpgs_proportions = visualize_name(s, cpgs_order, data_order, fill=TRUE)

  genes_order = c(
    'hg19_knownGenes_1to5kb',
    'hg19_knownGenes_promoters',
    'hg19_knownGenes_5UTRs',
    'hg19_knownGenes_exons',
    'hg19_knownGenes_introns',
    'hg19_knownGenes_3UTRs')
  data_order = c(
    'DMup',
    'DMdown')
  v_genes_counts = visualize_name(s, genes_order, data_order, fill=FALSE)
  v_genes_proportions = visualize_name(s, genes_order, data_order, fill=TRUE)

  cpgs_error_order = c(
    'hg19_cpg_islands',
    'hg19_cpg_shores',
    'hg19_cpg_shelves',
    'hg19_cpg_inter',
    'walla_walla')
  data_error_order = c(
    'DMup',
    'DMdown',
    'margaret_thatcher')

  expect_error( visualize_name(s, cpgs_error_order, data_order, fill=FALSE), 'annotations in annotation_order')
  expect_error( visualize_name(s, cpgs_order, data_error_order, fill=FALSE), 'data_names in data_order that are')
  expect_equal( dplyr::setequal(class(v_cpgs_counts), c('gg','ggplot')), expected = TRUE)
  expect_equal( dplyr::setequal(class(v_cpgs_proportions), c('gg','ggplot')), expected = TRUE)
  expect_equal( dplyr::setequal(class(v_genes_counts), c('gg','ggplot')), expected = TRUE)
  expect_equal( dplyr::setequal(class(v_genes_proportions), c('gg','ggplot')), expected = TRUE)
})
