context('Test read module')

test_that('Test error for file does not exist.',{
    expect_error(
        read_bed("Nonexistent File", genome = 'hg38', stranded = FALSE),
        'not found'
    )
})

test_that('Test correct length when no duplicates' , {
    filename1 = system.file('extdata', 'K562_Cjun.narrowPeak.gz', package = 'annotatr')
    gR1 <- read_bed(filename1, 'hg19', stranded = FALSE)

    expect_equal(length(gR1), expected = 9848)
})

test_that('Test stranded argument' , {
    filename1 = system.file('extdata', 'test_with_data_strand.bed', package = 'annotatr')
    gR1 <- read_bed(filename1, 'hg19', stranded = TRUE)

    expect_equal(length(gR1), expected = 3)
})

test_that('Test stranded argument throws error' , {
    filename1 = system.file('extdata', 'test_with_strand_error.bed', package = 'annotatr')

    expect_error(gR1 <- read_bed(filename1, 'hg19', stranded = TRUE, use.score = TRUE),
      'strand column should contain')
})

test_that('Unique works on ranges only', {
    file = system.file('extdata','test_duplicates_with_data_strand.bed',package='annotatr')
    gr = read_bed(file, 'hg19', stranded=TRUE, use.score=TRUE)

    expect_equal( length(gr), expected = 3)
  })

test_that('Test unique statement works', {
    file = system.file('extdata', 'Gm12878_Pol2.narrowPeak.gz', package = 'annotatr')
    gr = read_bed(file, 'hg19', stranded = FALSE)

    expect_equal( 27558 - length(gr), expected = 5821)
})
