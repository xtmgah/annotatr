context('Test that annotations can be loaded')

annots = supported_annotations()
for(annot in annots) {
  message(sprintf('Loading %s',annot))
  data(list = annot, package = "annotatr")
}

# There are 70 total things for ls() including annot and annots
test_that('All datasets successfully loaded',{
    expect_equal(ls()-2, 68)
})
