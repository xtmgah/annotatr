context('Test read module')

################################################################################
# Test errors

test_that('Test error for file does not exist.',{
    expect_error(
        read_bed("string", genome = 'hg38', stranded = FALSE),
        'not found'
    )
})

test_that('Test chromosome column error' , {
    file = system.file('extdata', 'test_read_error_chr.bed', package = 'annotatr')

    expect_error(read_bed(file, genome = 'hg19', stranded = TRUE, use.score = TRUE),
      'First column of BED file does not appear to be chromsome')
})

test_that('Test integer column errors' , {
    file1 = system.file('extdata', 'test_read_error_int1.bed', package = 'annotatr')
    file2 = system.file('extdata', 'test_read_error_int2.bed', package = 'annotatr')

    expect_error(read_bed(file1, genome = 'hg19', stranded = TRUE, use.score = TRUE),
      'Second column of BED file must be integer valued')
    expect_error(read_bed(file2, genome = 'hg19', stranded = TRUE, use.score = TRUE),
      'Third column of BED file must be integer valued')
})

test_that('Test stranded argument throws error' , {
    file = system.file('extdata', 'test_read_error_strand.bed', package = 'annotatr')

    expect_error(read_bed(file, genome = 'hg19', stranded = TRUE, use.score = TRUE),
      'strand column should contain')
})

################################################################################
# Test read_bed()

test_that('Test correct length when no duplicates' , {
    file = system.file('extdata', 'K562_Cjun.narrowPeak.gz', package = 'annotatr')
    gr = read_bed(file, genome = 'hg19', stranded = FALSE)

    expect_equal(length(gr), expected = 9848)
})

test_that('Test stranded argument' , {
    file = system.file('extdata', 'test_with_data_strand.bed', package = 'annotatr')
    gr = read_bed(file, genome = 'hg19', stranded = TRUE)

    expect_equal(length(gr), expected = 3)
})

test_that('Unique works on ranges only', {
    file = system.file('extdata','test_duplicates_with_data_strand.bed',package='annotatr')
    gr = read_bed(file, genome = 'hg19', stranded=TRUE, use.score=TRUE)

    expect_equal( length(gr), expected = 3)
  })

test_that('Test unique statement works', {
    file = system.file('extdata', 'Gm12878_Pol2.narrowPeak.gz', package = 'annotatr')
    gr = read_bed(file, genome = 'hg19', stranded = FALSE)

    expect_equal( 27558 - length(gr), expected = 5821)
})
