context('Test read_bed() function')

test_that('Test error for file does not exist.',{
    expect_error(
        read_bed("Nonexistent File", genome = 'hg38', stranded = F),
        'not found'
    )
})

test_that('Correct Lengths.' , {
    filename1 = system.file('extdata', 'K562_Cjun.narrowPeak.gz', package = 'annotatr')

    ###NOTE TO SELF: NOT SURE IF DATA HAS + or - STRAND COLUMN ###

    gR1 <- read_bed(filename1, 'hg19', stranded = F)

    expect_equal(length(gR1), expected = 9848)
})



