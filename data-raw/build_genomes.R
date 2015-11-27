genomes = c('hg19','hg38','mm9','mm10')

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
}
