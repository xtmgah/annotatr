context('Test annotations can be loaded')

# There are 70 total things for ls() including annot and annots
test_that('All datasets successfully loaded',{
  annots = supported_annotations()
  for(annot in annots) {
    data(list = annot, package = "annotatr", envir = environment())
  }

  expect_equal( length(ls()) - 2 , expected = length(supported_annotations()))
})
