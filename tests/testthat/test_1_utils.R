context('Test utility functions')

# There are 70 total things for ls() including annot and annots
test_that('Test supported_annotations() results can be loaded',{
  annots = supported_annotations()
  for(annot in annots) {
    data(list = annot, package = "annotatr", envir = environment())
  }

  expect_equal( length(ls()) - 2 , expected = length(supported_annotations()))
})

test_that('Test tidy_annotations()', {
  hg19_annots = c('hg19_cpg_islands','hg19_knownGenes_promoters','hg19_cpg_inter')
  mm9_annots = annots = c('mm9_cpg_islands','mm9_knownGenes_exonsCDSs','mm9_cpg_inter')

  hg19_tidy_annots = tidy_annotations(hg19_annots)
  mm9_tidy_annots = tidy_annotations(mm9_annots)

  expect_equal( all(hg19_tidy_annots == c('CpG islands', 'promoters', 'interCGI')), expected = TRUE)
  expect_equal( all(mm9_tidy_annots == c('CpG islands', 'exonsCDSs', 'interCGI')), expected = TRUE)
})
