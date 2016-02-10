# This code transforms data in data-raw/ into the .RData in data/
# NOTE: This script is run in data-raw/build_data.sh
# DO NOT RUN ALONE

# Load the annotatr package such as it is
devtools::load_all()

# Create file for roxygen2 documentation of objects in data/
doc_file = '../R/annotatr_data_doc.R'
cat('', file= doc_file, append=F)

################################################################################
# Process Enhancer related files

enh_files = list.files(pattern='renhtmp', full.names=T)

for(file in enh_files) {
  message(sprintf('Processing %s', file))

  # Get genome and create RData filename
    prefix = gsub('renhtmp_', '', unlist(strsplit(basename(file), '[.]'))[1])
    object_name = prefix
    genome = unlist(strsplit(prefix, '_'))[1]
    rdata_file = sprintf('../data/%s.RData', prefix)

  # Get seqlengths from data/
    size_code = sprintf('%s_chrom_sizes', genome)
    data(list = size_code, package = "annotatr")
    seqlengths = get(size_code)

  # Get tmp data
    tmp_data = read.table(file, sep='\t', header=F, stringsAsFactors=F)
    colnames(tmp_data) = c('chrom','start','end','ID')

  # Construct GRanges object
  # NOTE: GRanges objects are 1-based start and end while the input data is
  # 0-based start and 1-based end. Hence the addition of 1 to start.
    gr_data = GenomicRanges::GRanges(
      seqnames = tmp_data$chrom,
      ranges = IRanges::IRanges(start = tmp_data$start, end = tmp_data$end),
      strand = '*',
      ID = tmp_data$ID,
      seqlengths = seqlengths
    )
    GenomeInfoDb::genome(gr_data) = genome

  # Trim just in case
    gr_data_trimmed = GenomicRanges::trim(gr_data)

  # Give gr_data_trimmed the correct name
    assign(object_name, gr_data_trimmed)

  # Write the GRanges object to ../data/
  # Make sure the name of the object is specific to the object like chipenrich
  save(list=c(object_name), file=rdata_file, compress='xz')
  file.remove(file)

  # Write roxygen2 documentation for dataset to ../R/annotatr_data.R
  man = c(
    sprintf("#' %s", prefix),
    "#' ",
    sprintf("#' A GenomicRanges object for %s annotations.", prefix),
    "#' ",
    "#' Permissive enhancers are taken from the FANTOM5 consortium (http://enhancer.binf.ku.dk/presets/). Description of methods that result in the enhancers can be found in the paper corresponding paper (http://www.ncbi.nlm.nih.gov/pubmed/24670763).",
    "#' ",
    "#' @format A GenomicRanges object.",
    sprintf("#' @name %s", prefix),
    "#' @keywords datasets",
    sprintf("#' @usage data(%s)", prefix),
    "NULL",
    ""
  )
  cat(man, sep='\n', file=doc_file, append=T)

}

################################################################################
# Process CpG related files

cpg_files = list.files(pattern='rcpgtmp', full.names=T)

for(file in cpg_files) {
  message(sprintf('Processing %s', file))

  # Get genome and create RData filename
    prefix = gsub('rcpgtmp_', '', unlist(strsplit(basename(file), '[.]'))[1])
    object_name = prefix
    genome = unlist(strsplit(prefix, '_'))[1]
    rdata_file = sprintf('../data/%s.RData', prefix)

  # Get seqlengths from data/
    size_code = sprintf('%s_chrom_sizes', genome)
    data(list = size_code, package = "annotatr")
    seqlengths = get(size_code)

  # Get tmp data
    tmp_data = read.table(file, sep='\t', header=F, stringsAsFactors=F)
    colnames(tmp_data) = c('chrom','start','end','ID')

  # Construct GRanges object
  # NOTE: GRanges objects are 1-based start and end while the input data is
  # 0-based start and 1-based end. Hence the addition of 1 to start.
    gr_data = GenomicRanges::GRanges(
      seqnames = tmp_data$chrom,
      ranges = IRanges::IRanges(start = tmp_data$start, end = tmp_data$end),
      strand = '*',
      ID = tmp_data$ID,
      seqlengths = seqlengths
    )
    GenomeInfoDb::genome(gr_data) = genome

  # Trim just in case
    gr_data_trimmed = GenomicRanges::trim(gr_data)

  # Give gr_data_trimmed the correct name
    assign(object_name, gr_data_trimmed)

  # Write the GRanges object to ../data/
  # Make sure the name of the object is specific to the object like chipenrich
  save(list=c(object_name), file=rdata_file, compress='xz')
  file.remove(file)

  # Write roxygen2 documentation for dataset to ../R/annotatr_data.R
  man = c(
    sprintf("#' %s", prefix),
    "#' ",
    sprintf("#' A GenomicRanges object for %s annotations.", prefix),
    "#' ",
    "#' CpG islands are taken from the eponymous track in the Regulation group from the corresponding organism in the UCSC Table Browser. CpG shores are defined as flanking 2Kb sequence upstream/downstream of the island boundaries. CpG shelves are defined as flanking 2Kb further upstream/downstream from the CpG shores. InterCGI is the remaining genomic sequence.",
    "#' ",
    "#' @format A GenomicRanges object.",
    sprintf("#' @name %s", prefix),
    "#' @keywords datasets",
    sprintf("#' @usage data(%s)", prefix),
    "NULL",
    ""
  )
  cat(man, sep='\n', file=doc_file, append=T)

}

################################################################################
# Process knownGenes related files

kg_files = list.files(pattern='rkgtmp', full.names=T)

for(file in kg_files) {
  message(sprintf('Processing %s', file))

  # Get genome and create RData filename
    prefix = gsub('rkgtmp_', '', unlist(strsplit(basename(file), '[.]'))[1])
    object_name = prefix
    genome = unlist(strsplit(prefix, '_'))[1]
    rdata_file = sprintf('../data/%s.RData', prefix)

  # Get seqlengths from data/
    size_code = sprintf('%s_chrom_sizes', genome)
    data(list = size_code, package = "annotatr")
    seqlengths = get(size_code)

  # Get tmp data
    tmp_data = read.table(file, sep='\t', header=F, stringsAsFactors=F)
    colnames(tmp_data) = c('chrom','start','end','strand','ID')

  # Construct GRanges object
  # NOTE: GRanges objects are 1-based start and end while the input data is
  # 0-based start and 1-based end. Hence the addition of 1 to start.
    gr_data = GenomicRanges::GRanges(
      seqnames = tmp_data$chrom,
      ranges = IRanges::IRanges(start = tmp_data$start+1, end = tmp_data$end),
      strand = tmp_data$strand,
      ID = tmp_data$ID,
      seqlengths = seqlengths
    )
    GenomeInfoDb::genome(gr_data) = genome

  # Trim just in case
    gr_data_trimmed = GenomicRanges::trim(gr_data)

  # Give gr_data_trimmed the correct name
    assign(object_name, gr_data_trimmed)

  # Write the GRanges object to ../data/
  # Make sure the name of the object is specific to the object like chipenrich
  save(list=c(object_name), file=rdata_file, compress='xz')
  file.remove(file)

  # Write roxygen2 documentation for dataset to ../R/annotatr_data.R
  man = c(
    sprintf("#' %s", prefix),
    "#' ",
    sprintf("#' A GenomicRanges object for %s annotations. KnownGenes are taken from the Genes group, UCSC Genes track, and knownGenes table of the UCSC Table Browser.", prefix),
    "#' Promoter sequence is defined as <1Kb upstream of a transcription start site (TSS, txStart). 1-5Kb upstream sequence is exactly that, upstream of the promoter sequence. 5-prime UTR sequence is defined as the txStart to the cdsStart. Coding sequence (CDS) is defined as cdsStart to cdsEnd. Exons are defined as the matching intervals of exonStarts to exonEnds. Introns are defined as filling the space beween the exons. 3-prime UTR sequence is defined as cdsEnd to txEnd.",
    "#' ",
    "#' More fine-grained definitions include exons and introns that occur in 5-prime UTR, CDS, and 3-prime UTR sequence. First exons and introns are exactly that, while respecting strandedness.",
    "#' ",
    "#' Intergenic annotations are the complement of the union of UCSC Genome Browser gap tracks and knownGene annotation tracks.",
    "#' ",
    "#' @format A GenomicRanges object.",
    sprintf("#' @name %s", prefix),
    "#' @keywords datasets",
    sprintf("#' @usage data(%s)", prefix),
    "NULL",
    ""
  )
  cat(man, sep='\n', file=doc_file, append=T)

}
