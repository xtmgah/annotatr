context('Test visualization of summaries')

test_that('Test error on incorrect input class in visualize functions', {
  bed = system.file('extdata', 'Gm12878_Ezh2_sorted_scores.narrowPeak.gz', package = 'annotatr')

  expect_error(visualize_annotation(bed), 'summarized_annotations must have class tbl_df')
  expect_error(visualize_score(bed), 'summarized_scores must have class grouped_df')
  expect_error(visualize_name(bed), 'summarized_names must have class grouped_df')
})

test_that('Test visualize_annotation() on signalValue from ChIP-seq data', {
  bed = system.file('extdata', 'Gm12878_Ezh2_sorted_scores.narrowPeak.gz', package = 'annotatr')
  annotations = c('basic_genes','cpgs')

  d = read_bed(filename = bed, genome = 'hg19', stranded = F, use.score = TRUE)

  i = annotate_regions(
    regions = d,
    annotations = annotations,
    genome = 'hg19',
    ignore.strand = T,
    use.score = TRUE)

  s = summarize_annotation(i)

  annots_order = c(
    'hg19_cpg_islands',
    'hg19_cpg_shores',
    'hg19_cpg_shelves',
    'hg19_cpg_inter',
    'hg19_knownGenes_1to5kb',
    'hg19_knownGenes_promoters',
    'hg19_knownGenes_5UTRs',
    'hg19_knownGenes_exons',
    'hg19_knownGenes_introns',
    'hg19_knownGenes_3UTRs')

  v_annots = visualize_annotation(s, annots_order)

  error_order = c('bingbong')

  expect_error( visualize_annotation(s, error_order), 'annotations in annotation_order that are not present')
  expect_equal( dplyr::setequal(class(v_annots), c('gg','ggplot')), expected = TRUE)
})

test_that('Test visualize_score() on signalValue from ChIP-seq data', {
  bed = system.file('extdata', 'Gm12878_Ezh2_sorted_scores.narrowPeak.gz', package = 'annotatr')
  annotations = c('basic_genes','cpgs')

  d = read_bed(filename = bed, genome = 'hg19', stranded = F, use.score = TRUE)

  i = annotate_regions(
    regions = d,
    annotations = annotations,
    genome = 'hg19',
    ignore.strand = T,
    use.score = TRUE)

  s = summarize_score(i)

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

  i = annotate_regions(
    regions = d,
    annotations = annotations,
    genome = 'hg19',
    ignore.strand = T,
    use.score = TRUE)

  s = summarize_score(i)

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

  i = annotate_regions(
    regions = d,
    annotations = annotations,
    genome = 'hg19',
    ignore.strand = T,
    use.score = TRUE)

  s = summarize_name(i)

  ################################################
  # Test one order on CpGs and all positions
  fill_order = c(
    'hg19_cpg_islands',
    'hg19_cpg_shores',
    'hg19_cpg_shelves',
    'hg19_cpg_inter')
  x_order = c(
    'DMup',
    'DMdown')
  v_cpgs_counts_data_annot = visualize_name(summarized_names=s, x='data_name', fill='annot_type',
    x_order = x_order, fill_order = fill_order, position='stack')
  v_cpgs_proportions_data_annot = visualize_name(summarized_names=s, x='data_name', fill='annot_type',
    x_order = x_order, fill_order = fill_order, position='fill')
  v_cpgs_nofill_data = visualize_name(summarized_names=s, x='data_name', fill=NULL,
    x_order = x_order, fill_order = fill_order, position='stack')

  ################################################
  # Test the other order on CpGs and all positions
  x_order = c(
    'hg19_cpg_islands',
    'hg19_cpg_shores',
    'hg19_cpg_shelves',
    'hg19_cpg_inter')
  fill_order = c(
    'DMup',
    'DMdown',
    'noDM')
  v_cpgs_counts_annot_data = visualize_name(summarized_names=s, x='annot_type', fill='data_name',
    x_order = x_order, fill_order = fill_order, position='stack')
  v_cpgs_proportions_annot_data = visualize_name(summarized_names=s, x='annot_type', fill='data_name',
    x_order = x_order, fill_order = fill_order, position='fill')
  v_cpgs_nofill_annot = visualize_name(summarized_names=s, x='annot_type', fill=NULL,
    x_order = x_order, fill_order = fill_order, position='stack')

  ################################################
  # Test on knownGenes and with error order vectors
  fill_order = c(
    'hg19_knownGenes_1to5kb',
    'hg19_knownGenes_promoters',
    'hg19_knownGenes_5UTRs',
    'hg19_knownGenes_exons',
    'hg19_knownGenes_introns',
    'hg19_knownGenes_3UTRs')
  x_order = c(
    'DMup',
    'DMdown')
  v_genes_counts = visualize_name(s, x='data_name', fill='annot_type',
    x_order = x_order, fill_order = fill_order, position='stack')
  v_genes_proportions = visualize_name(s, x='data_name', fill='annot_type',
    x_order = x_order, fill_order = fill_order, position='stack')

  x_order = c(
    'hg19_cpg_islands',
    'hg19_cpg_shores',
    'hg19_cpg_shelves',
    'hg19_cpg_inter')
  fill_order = c(
    'DMup',
    'DMdown',
    'noDM')
  x_error_order = c(
    'hg19_cpg_islands',
    'hg19_cpg_shores',
    'hg19_cpg_shelves',
    'hg19_cpg_inter',
    'walla_walla')
  fill_error_order = c(
    'DMup',
    'DMdown',
    'margaret_thatcher')

  expect_error( visualize_name(s, x='annot_type', fill='data_name', x_order=x_error_order, fill_order=fill_order, position='stack'), 'elements in x_order that are not present')
  expect_error( visualize_name(s, x='annot_type', fill='data_name', x_order=x_order, fill_order=fill_error_order, position='stack'), 'elements in fill_order that are not present')
  expect_equal( dplyr::setequal(class(v_cpgs_counts_data_annot), c('gg','ggplot')), expected = TRUE)
  expect_equal( dplyr::setequal(class(v_cpgs_proportions_annot_data), c('gg','ggplot')), expected = TRUE)
  expect_equal( dplyr::setequal(class(v_cpgs_nofill_data), c('gg','ggplot')), expected = TRUE)
  expect_equal( dplyr::setequal(class(v_cpgs_counts_annot_data), c('gg','ggplot')), expected = TRUE)
  expect_equal( dplyr::setequal(class(v_cpgs_proportions_annot_data), c('gg','ggplot')), expected = TRUE)
  expect_equal( dplyr::setequal(class(v_cpgs_nofill_annot), c('gg','ggplot')), expected = TRUE)
  expect_equal( dplyr::setequal(class(v_genes_counts), c('gg','ggplot')), expected = TRUE)
  expect_equal( dplyr::setequal(class(v_genes_proportions), c('gg','ggplot')), expected = TRUE)
})
