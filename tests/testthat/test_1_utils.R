context('Test utility functions')

test_that('Test tidy_annotations()', {
  hg19_annots = c('hg19_cpg_islands','hg19_knownGenes_promoters','hg19_cpg_inter')
  mm9_annots = c('mm9_cpg_islands','mm9_knownGenes_exonsCDSs','mm9_cpg_inter')
  rn4_custom_annots = c('rn4_custom_cpgislands','rn4_custom_TFBS')

  hg19_tidy_annots = tidy_annotations(hg19_annots)
  mm9_tidy_annots = tidy_annotations(mm9_annots)
  rn4_tidy_annots = tidy_annotations(rn4_custom_annots)

  expect_equal( all(names(hg19_tidy_annots) == c('CpG islands', 'promoters', 'interCGI')), expected = TRUE)
  expect_equal( all(names(mm9_tidy_annots) == c('CpG islands', 'exonsCDSs', 'interCGI')), expected = TRUE)
  expect_equal( all(names(rn4_tidy_annots) == c('cpgislands', 'TFBS')), expected = TRUE)
})

test_that('Test check_annotations()', {
  annots1 = c('hg17_knownGenes_promoters','hg19_cpgs')
  annots2 = c('hello','hg19_knownGenes_promoters','hg19_cpgs')
  annots3 = c('hg19_knownGenes_promoters', 'mm9_cpg_islands')
  annots4 = c('mm10_basicgenes', 'mm10_detailedgenes')

  expect_error( check_annotations(annots1), 'not supported. See supported_annotations()' )
  expect_error( check_annotations(annots2), 'not supported. See supported_annotations()' )
  expect_error( check_annotations(annots3), 'genome prefix on all annotations must be the same' )
  expect_error( check_annotations(annots4), 'basicgenes and detailedgenes shortcuts may not be used' )
})

test_that('Test expand_annotations()', {
  annots1 = c('hg19_knownGenes_promoters', 'hg19_knownGenes_exons')
  annots2 = c('mm10_basicgenes', 'mm10_cpgs')
  expanded_annots2 = c('mm10_cpg_islands', 'mm10_cpg_shores', 'mm10_cpg_shelves', 'mm10_cpg_inter', 'mm10_knownGenes_1to5kb', 'mm10_knownGenes_promoters', 'mm10_knownGenes_5UTRs', 'mm10_knownGenes_exons', 'mm10_knownGenes_introns', 'mm10_knownGenes_3UTRs')
  annots3 = c('hg38_cpg_shores', 'hg38_detailedgenes', 'hg38_cpg_islands')
  expanded_annots3 = c('hg38_cpg_islands', 'hg38_cpg_shores', 'hg38_knownGenes_1to5kb', 'hg38_knownGenes_promoters', 'hg38_knownGenes_exons5UTRs', 'hg38_knownGenes_introns5UTRs', 'hg38_knownGenes_exonsCDSs', 'hg38_knownGenes_intronsCDSs', 'hg38_knownGenes_exons3UTRs', 'hg38_knownGenes_introns3UTRs')

  expect_equal( dplyr::setequal(expand_annotations(annots1), annots1), expected = TRUE )
  expect_equal( dplyr::setequal(expand_annotations(annots2), expanded_annots2), expected = TRUE )
  expect_equal( dplyr::setequal(expand_annotations(annots3), expanded_annots3), expected = TRUE )
})
