context('Test visualization module')

################################################################################
# Setup objects for visualize_annotation()

  chip = system.file('extdata', 'Gm12878_Ezh2_sorted_scores.narrowPeak.gz', package = 'annotatr')
  annotations = c('basic_genes','cpgs')

  chip_d = read_bed(file = chip, genome = 'hg19', stranded = F, use.score = TRUE)

  chip_r = annotate_regions(
    regions = chip_d,
    annotations = annotations,
    genome = 'hg19',
    ignore.strand = T,
    use.score = TRUE)

  chip_sa = summarize_annotation(chip_r)

################################################################################
# Setup objects for visualize_score()

  # Use the same input file as above, but summarize_score instead
  chip_ss = summarize_score(chip_r)

################################################################################
# Setup objects for visualize_name()

  dm = system.file('extdata', 'IDH2mut_v_NBM_names_scores_chr9.txt.gz', package = 'annotatr')
  annotations = c('basic_genes','cpgs')

  dm_d = read_bed(file = dm, genome = 'hg19', stranded = F, use.score = TRUE)

  dm_r = annotate_regions(
    regions = dm_d,
    annotations = annotations,
    genome = 'hg19',
    ignore.strand = T,
    use.score = TRUE)

  dm_sn = summarize_name(dm_r)

################################################################################
# Setup order vectors and visualizations that will work

  dm_order = c(
    'hyper',
    'hypo',
    'none')
  all_order = c(
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
  cpgs_order = c(
    'hg19_cpg_islands',
    'hg19_cpg_shores',
    'hg19_cpg_shelves',
    'hg19_cpg_inter')
  genes_order = c(
    'hg19_knownGenes_1to5kb',
    'hg19_knownGenes_promoters',
    'hg19_knownGenes_5UTRs',
    'hg19_knownGenes_exons',
    'hg19_knownGenes_introns',
    'hg19_knownGenes_3UTRs')

################################################################################
# Test errors common to all functions

  test_that('Test error on incorrect input class in all visualize functions', {
    expect_error(visualize_annotation(chip), 'summarized_annotations must have class tbl_df')
    expect_error(visualize_score(chip), 'summarized_scores must have class grouped_df')
    expect_error(visualize_name(chip), 'summarized_names must have class grouped_df')
  })

################################################################################
# Test visualize_annotation()

  test_that('Test visualize_annotation() errors', {
    expect_error( visualize_annotation(chip_sa, dm_order), 'annotations in annotation_order that are not present')
  })

  test_that('Test visualize_annotation() success', {
    chip_va_min = visualize_annotation(summarized_annotations = chip_sa)

    chip_va = visualize_annotation(
      summarized_annotations = chip_sa,
      annotation_order = all_order,
      plot_title = 'Testing plot title')

    expect_equal( dplyr::setequal(class(chip_va_min), c('gg','ggplot')), expected = TRUE)
    expect_equal( dplyr::setequal(class(chip_va), c('gg','ggplot')), expected = TRUE)
  })

################################################################################
# Test visualize_score()

  test_that('Test visualize_score() errors', {
    expect_error( visualize_score(chip_ss, dm_order), 'annotations in annotation_order that are not present')
  })

  test_that('Test visualize_score() success', {
    chip_vs_min = visualize_score(
      summarized_scores = chip_ss,
      bin_width = 30)

    chip_vs = visualize_score(
      summarized_scores = chip_ss,
      annotation_order = genes_order,
      x_label = 'Test x-axis label')

    expect_equal( dplyr::setequal(class(chip_vs_min), c('gg','ggplot')), expected = TRUE)
    expect_equal( dplyr::setequal(class(chip_vs), c('gg','ggplot')), expected = TRUE)
  })

################################################################################
# Test visualize_name()

  test_that('Test visualize_name() errors', {
    expect_error(
      visualize_name(summarized_names = dm_sn),
      'argument "x" is missing')

    expect_error(
      visualize_name(
        summarized_names = dm_sn,
        x = 'testing'),
      'column name used for x does not exist in summarized_names')

    expect_error(
      visualize_name(
        summarized_names = dm_sn,
        x = 'name',
        fill = 'testing'),
      'column name used for fill does not exist in summarized_names')

    expect_error(
      visualize_name(
        summarized_names = dm_sn,
        x = 'name',
        fill = 'name'),
      'x cannot equal fill')

    expect_error(
      visualize_name(
        summarized_names = dm_sn,
        x = 'name',
        fill = 'annot_type',
        position = 'no'),
      'position must be one of "stack", "fill"')

    expect_error(
      visualize_name(
        summarized_names = dm_sn,
        x = 'name',
        fill = 'annot_type',
        x_order = cpgs_order),
      'are elements in x_order that are not present')

    expect_error(
      visualize_name(
        summarized_names = dm_sn,
        x = 'name',
        fill = 'annot_type',
        fill_order = dm_order),
      'elements in fill_order that are not present')
  })

test_that('Test visualize_name() success', {
  dm_vn_min = visualize_name(
    summarized_names = dm_sn,
    x = 'name')

  dm_vn = visualize_name(
    summarized_names = dm_sn,
    x = 'name',
    fill = 'annot_type',
    x_order = dm_order,
    fill_order = genes_order,
    position = 'fill',
    legend_title = 'knownGene Annotations',
    x_label = 'DM status',
    y_label = 'Proportion')

  expect_equal( dplyr::setequal(class(dm_vn_min), c('gg','ggplot')), expected = TRUE)
  expect_equal( dplyr::setequal(class(dm_vn), c('gg','ggplot')), expected = TRUE)
})
