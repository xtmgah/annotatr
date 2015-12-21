context('Test visualization module')

################################################################################
# Setup objects for visualize_annotation()

  chip = system.file('extdata', 'Gm12878_Ezh2_sorted_scores.narrowPeak.gz', package = 'annotatr')
  annotations = c('hg19_basicgenes','hg19_cpgs')

  chip_d = read_bed(file = chip, genome = 'hg19', stranded = F, use.score = TRUE)

  chip_r = annotate_regions(
    regions = chip_d,
    annotations = annotations,
    ignore.strand = T,
    use.score = TRUE)

  chip_sa = summarize_annotations(chip_r)

################################################################################
# Setup objects for visualize_numerical()

  # Use the same input file as above, but summarize_numerical instead
  chip_ss = summarize_numerical(
    annotated_regions = chip_r,
    by = c('annot_type', 'annot_id'))

################################################################################
# Setup objects for visualize_categorical()

  dm = system.file('extdata', 'IDH2mut_v_NBM_multi_data_chr9.txt.gz', package = 'annotatr')
  annotations = c('hg19_basicgenes','hg19_cpgs')

  dm_d = read_bed(
    file = dm,
    col.names=c('chr','start','end','DM_status','pval','strand','diff_meth','mu1','mu0'),
    genome = 'hg19',
    stranded = F,
    use.score = TRUE)

  dm_r = annotate_regions(
    regions = dm_d,
    annotations = annotations,
    ignore.strand = T,
    use.score = TRUE)

  dm_sn = summarize_numerical(
    annotated_regions = dm_r,
    over = 'diff_meth')

  dm_sc = summarize_categorical(
    annotated_regions = dm_r,
    by = c('annot_type', 'name'))

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
    expect_error(visualize_numerical(chip), 'tbl must have class tbl_df or grouped_df')
    expect_error(visualize_categorical(chip), 'summarized_cats must have class grouped_df')
  })

################################################################################
# Test visualize_annotation()

  test_that('Test visualize_annotation() errors', {
    expect_error( visualize_annotation(chip_sa, dm_order), 'elements in col_order that are not present')
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
# Test visualize_numerical()

  test_that('Test visualize_numerical() errors', {
    expect_error(
      visualize_numerical(
        tbl = chip_ss,
        x = mean,
        facet_order = dm_order),
      'elements in col_order that are not present')
  })

  test_that('Test visualize_numerical() success', {

    chip_vs_min = visualize_numerical(
      tbl = chip_ss,
      x = 'mean',
      bin_width = 30)

    chip_vs = visualize_numerical(
      tbl = chip_ss,
      x = 'mean',
      facet_order = genes_order,
      x_label = 'Test x-axis label')

    dm_vs_regions_annot = visualize_numerical(
      tbl = dm_r,
      x = 'mu0',
      y = 'mu1',
      facet = 'annot_type',
      facet_order = c('hg19_knownGenes_1to5kb','hg19_knownGenes_promoters','hg19_knownGenes_5UTRs','hg19_knownGenes_3UTRs'),
      plot_title = 'Region Methylation: Group 0 vs Group 1',
      x_label = 'Group 0',
      y_label = 'Group 1')

    dm_vs_regions_name = visualize_numerical(
      tbl = dm_r,
      x = 'mu0',
      y = 'mu1',
      facet = 'name',
      facet_order = c('hyper','hypo','none'),
      plot_title = 'Region Methylation: Group 0 vs Group 1',
      x_label = 'Group 0',
      y_label = 'Group 1')

    dm_vs_sumnum = visualize_numerical(
      tbl = dm_sn,
      x = 'mean',
      facet = 'annot_type',
      facet_order = c('hg19_cpg_islands','hg19_cpg_shores','hg19_cpg_shelves','hg19_cpg_inter'),
      bin_width = 5,
      plot_title = 'Mean Meth. Diff. over CpG Annots.',
      x_label = 'Methylation Difference')

    expect_equal( dplyr::setequal(class(chip_vs_min), c('gg','ggplot')), expected = TRUE)
    expect_equal( dplyr::setequal(class(chip_vs), c('gg','ggplot')), expected = TRUE)
    expect_equal( dplyr::setequal(class(dm_vs_regions_annot), c('gg','ggplot')), expected = TRUE)
    expect_equal( dplyr::setequal(class(dm_vs_regions_name), c('gg','ggplot')), expected = TRUE)
    expect_equal( dplyr::setequal(class(dm_vs_sumnum), c('gg','ggplot')), expected = TRUE)
  })

################################################################################
# Test visualize_categorical()

  test_that('Test visualize_categorical() errors', {
    expect_error(
      visualize_categorical(
        summarized_cats = dm_sc),
      'argument "x" is missing')

    expect_error(
      visualize_categorical(
        summarized_cats = dm_sc,
        x = 'testing'),
      'column name used for x does not exist in summarized_cats')

    expect_error(
      visualize_categorical(
        summarized_cats = dm_sc,
        x = 'name',
        fill = 'testing'),
      'column name used for fill does not exist in summarized_cats')

    expect_error(
      visualize_categorical(
        summarized_cats = dm_sc,
        x = 'name',
        fill = 'name'),
      'x cannot equal fill')

    expect_error(
      visualize_categorical(
        summarized_cats = dm_sc,
        x = 'name',
        fill = 'annot_type',
        position = 'no'),
      'position must be one of "stack", "fill"')

    expect_error(
      visualize_categorical(
        summarized_cats = dm_sc,
        x = 'name',
        fill = 'annot_type',
        x_order = cpgs_order),
      'elements in col_order that are not present')

    expect_error(
      visualize_categorical(
        summarized_cats = dm_sc,
        x = 'name',
        fill = 'annot_type',
        fill_order = dm_order),
      'elements in col_order that are not present')
  })

test_that('Test visualize_categorical() success', {
  dm_vn_min = visualize_categorical(
    summarized_cats = dm_sc,
    x = 'annot_type')

  dm_vn = visualize_categorical(
    summarized_cats = dm_sc,
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
