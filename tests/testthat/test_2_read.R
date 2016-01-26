context('Test read module')

################################################################################
# Test errors in read_bed()

  test_that('Test error for file does not exist.',{
    expect_error(
      read_bed("string", genome = 'hg38', stranded = FALSE),
      'not found'
    )
  })

  test_that('Test invalid genome' , {
    file = system.file('extdata', 'K562_Cjun.narrowPeak.gz', package = 'annotatr')

    expect_warning(read_bed(file, genome = 'hg99', stranded = FALSE, use.score = TRUE),
      'not a supported genome. In order to annotate regions, make sure to load custom annotations with read_annotations()')
  })

  test_that('Test chromosome column error' , {
    file = system.file('extdata', 'test_read_error_chr.bed', package = 'annotatr')

    expect_error(read_bed(file, genome = 'hg19', stranded = TRUE, use.score = TRUE),
      'First column does not appear to be chromosome')
  })

  test_that('Test stranded argument throws error' , {
    file = system.file('extdata', 'test_read_error_strand.bed', package = 'annotatr')

    expect_error(read_bed(file, genome = 'hg19', stranded = TRUE, use.score = TRUE),
      'strand column should contain')
  })

################################################################################
# Test column handling in read_bed()

  test_that('Test 3 column BED parameter',{
    file = system.file('extdata', 'test_read_3cols.bed', package='annotatr')
    gr = read_bed(file, col.names = FALSE, genome = 'hg19', stranded = FALSE, use.score = FALSE)

    expect_warning( read_bed(file, col.names = FALSE, genome = 'hg19', stranded = FALSE, use.score = FALSE),
      'Input is not BED6, and no column names were given')
    expect_equal( dplyr::setequal(dim(GenomicRanges::mcols(gr)), c(3,0)), expected = TRUE )
  })

  test_that('Test col.names parameter',{
    file = system.file('extdata', 'test_read_colnames.bed', package='annotatr')
    gr = read_bed(file, col.names = TRUE, genome = 'hg19', stranded = FALSE, use.score = FALSE)

    expect_equal( dplyr::setequal(colnames(GenomicRanges::mcols(gr)), 'name'), expected = TRUE )
  })

  test_that('Test multiple data columns with header',{
    file = system.file('extdata', 'test_read_multiple_data_head.bed', package='annotatr')
    gr = read_bed(file, col.names = TRUE, genome = 'hg19', stranded = FALSE, use.score = TRUE)

    expect_equal(
      dplyr::setequal(colnames(GenomicRanges::mcols(gr)), c('name','score','pval','mu1','mu0','diff_exp')), expected = TRUE )
  })

  test_that('Test multiple data columns without header',{
    file = system.file('extdata', 'test_read_multiple_data_nohead.bed', package='annotatr')
    gr = read_bed(file,
      col.names = c('chr','start','end','name','diff_meth','strand','pval','mu1','mu0'),
      genome = 'hg19', stranded = FALSE, use.score = TRUE)

    expect_equal( dplyr::setequal(colnames(GenomicRanges::mcols(gr)), c('name','diff_meth','pval','mu1','mu0')), expected = TRUE )
})

################################################################################
# Test read_bed() overall

  test_that('Test correct length when no duplicates' , {
    file = system.file('extdata', 'K562_Cjun.narrowPeak.gz', package = 'annotatr')
    gr = read_bed(file, col.names = FALSE, genome = 'hg19', stranded = FALSE, use.score = FALSE)

    expect_equal(length(gr), expected = 9848)
  })

  test_that('Test stranded argument' , {
    file = system.file('extdata', 'test_with_data_strand.bed', package = 'annotatr')
    gr = read_bed(file, col.names = FALSE, genome = 'hg19', stranded = TRUE, use.score = FALSE)

    expect_equal(length(gr), expected = 3)
  })

  test_that('Unique works on ranges only', {
    file = system.file('extdata','test_duplicates_with_data_strand.bed',package='annotatr')
    gr = read_bed(file, col.names = FALSE, genome = 'hg19', stranded=TRUE, use.score=TRUE)

    expect_equal( length(gr), expected = 3)
  })

  test_that('Test unique statement works', {
    file = system.file('extdata', 'Gm12878_Pol2.narrowPeak.gz', package = 'annotatr')
    gr = read_bed(file, col.names = FALSE, genome = 'hg19', stranded = FALSE, use.score = FALSE)

    expect_equal( 27558 - length(gr), expected = 5821)
  })

  test_that('Test multiple column data', {
    file = system.file('extdata', 'IDH2mut_v_NBM_multi_data_chr9.txt.gz', package = 'annotatr')
    gr = read_bed(
      file,
      col.names=c('chr','start','end','name','pval','strand','diff_meth','mu1','mu0'),
      genome = 'hg19',
      stranded = FALSE,
      use.score = TRUE)

    expect_equal( dplyr::setequal(colnames(GenomicRanges::mcols(gr)), c('name','pval','diff_meth','mu1','mu0')), expected = TRUE)
  })

################################################################################
# Test read_df()

  test_that('Test read_df() works', {
    df = data.frame(
      chr = c('chr1','chr1','chr1','chr1'),
      start = c(10800, 11000, 27800, 29000),
      end = c(10900, 11100, 28800, 29300),
      name = c('A','A','B','B'),
      score = c(87, 45, 34, 62),
      strand = c('*','*','*','*'),
      stringsAsFactors=F)

    gr = read_df(
      df = df,
      genome = 'hg19',
      stranded = FALSE,
      use.score=TRUE)

    expect_equal(class(gr)[1], expected = 'GRanges')
    expect_equal(gr[1]$name, expected = 'A')
  })

################################################################################
# Test pre-read errors in read_annotations()

  test_that('Test error for file does not exist',{
    expect_error(
      read_annotations(file = 'hello', genome = 'hg19', annotation_name = 'test'),
      'not found'
      )
  })

  test_that('Test error for NULL genome',{
    file = system.file('extdata', 'test_annotations_3.bed', package='annotatr')

    expect_error(
      read_annotations(file = file, annotation_name = 'test'),
      'argument "genome" is missing'
      )
  })

################################################################################
# Test column handling in read_annotations()

  test_that('Test chrom error',{
    file = system.file('extdata', 'test_annotations_3_error_chr.bed', package='annotatr')

    expect_error(
      read_annotations(file = file, genome = 'hg19', annotation_name = 'custom'),
      'First column of annotation file does not appear'
    )
  })

  test_that('Test duplicate name warning',{
    file = system.file('extdata', 'test_annotations_4_warning_dup.bed', package='annotatr')

    expect_warning(
      read_annotations(file = file, genome = 'hg19', annotation_name = 'custom'),
      'Some annotations have duplicated names'
    )

    test = read_annotations(file = file, genome = 'hg19', annotation_name = 'custom')
    expect_equal(all(test$ID == c('custom:1', 'custom:2', 'custom:3')), expected = TRUE)
  })

  test_that('Test invalid strand ',{
    file = system.file('extdata', 'test_annotations_5_error_strand.bed', package='annotatr')

    expect_error(
      read_annotations(file = file, genome = 'hg19', annotation_name = 'custom'),
      'Strand column should contain only'
    )
  })

  test_that('Test not 3, 4, or 5 columns', {
    file = system.file('extdata', 'test_read_multiple_data_nohead.bed', package='annotatr')

    expect_error(
      read_annotations(file = file, genome = 'hg19', annotation_name = 'custom'),
      'The file used to create custom annotations must'
      )
  })

################################################################################
# Test read_annotations() overall

  test_that('Test 3 columns works', {
    file = system.file('extdata', 'test_annotations_3.bed', package='annotatr')

    hg19_custom_test = read_annotations(file = file, genome='hg19', annotation_name = 'test')

    expect_equal(length(hg19_custom_test), expected = 3)
    expect_equal(all(hg19_custom_test$ID == c('test:1','test:2','test:3')), expected = TRUE)
    expect_equal(as.logical(GenomeInfoDb::genome(hg19_custom_test) == 'hg19'), expected = TRUE)
  })

  test_that('Test 4 columns works', {
    file = system.file('extdata', 'test_annotations_4.bed', package='annotatr')

    expect_warning(
      read_annotations(file = file, genome='mm14', annotation_name = 'test'),
      'genome is not supported, and cannot be used with any built-in annotations'
      )

    mm14_custom_test = read_annotations(file = file, genome='mm14', annotation_name = 'test')

    expect_equal(length(mm14_custom_test), expected = 3)
    expect_equal(all(mm14_custom_test$ID == c('test1','test2','test3')), expected = TRUE)
    expect_equal(as.logical(GenomeInfoDb::genome(mm14_custom_test) == 'mm14'), expected = TRUE)
  })

  test_that('Test 5 columns works', {
    file = system.file('extdata', 'test_annotations_5.bed', package='annotatr')

    hg19_custom_TFBS = read_annotations(file = file, genome='hg19', annotation_name = 'TFBS')

    expect_equal(length(hg19_custom_TFBS), expected = 3)
    expect_equal(all(hg19_custom_TFBS$ID == c('test1','test3','test2')), expected = TRUE)
    expect_equal(as.logical(GenomeInfoDb::genome(hg19_custom_TFBS) == 'hg19'), expected = TRUE)
  })
