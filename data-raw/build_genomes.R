genomes = c('hg19','hg38','mm9','mm10')

# Create file for roxygen2 documentation of objects in data/
doc_file = '../R/annotatr_genomes_doc.R'
cat('', file= doc_file, append=F)

for(genome in genomes) {
  message(sprintf('Writing chromosome sizes for %s', genome))

  # Get chromosome sizes, and make seqlengths and seqinfo objects for
  # GRanges object creation
  chrom_sizes = read.table(sprintf('%s.chrom.sizes.txt', genome), sep='\t', header=F, stringsAsFactors=F)
  colnames(chrom_sizes) = c('chrom','length')

  seqlengths = as.integer(chrom_sizes$length)
  names(seqlengths) = chrom_sizes$chrom

  object_name = sprintf('%s_chrom_sizes', genome)
  rdata_file = sprintf('../data/%s_chrom_sizes.RData', genome)
  assign(object_name, seqlengths)

  save(list=c(object_name), file=rdata_file)

  # Write roxygen2 documentation for chrom_sizes to ../R/annotatr_genomes_doc.R
  man = c(
    sprintf("#' %s_chrom_sizes", genome),
    "#' ",
    sprintf("#' Chromosome sizes for %s", genome),
    "#' ",
    "#' Chromosome sizes come from UCSC Genome Browser.",
    "#' ",
    "#' @format A numeric vector with names of chromsomes and values of length of chromsomes.",
    sprintf("#' @name %s_chrom_sizes", genome),
    "#' @keywords datasets",
    sprintf("#' @usage data(%s_chrom_sizes)", genome),
    "NULL",
    ""
  )
  cat(man, sep='\n', file=doc_file, append=T)
}
