#' Convert a bed file into a GenomicRanges object
#'
#' \code{read_bed} reads in data from a BED6 file, checks to see if it follows the appropriate BED (Browser Extensible Data) format, and checks validity of genome. If valid, reads file data into a GenomicRanges object. Unlike in the BED specification, the score column can be continuously valued. The names need not be unique.
#'
#' @param file Path to the file with data. File must exist and be in correct format.
#' @param col.names Either \code{TRUE} (column names are given in the first row of the file), \code{FALSE} (the first row does not have column names), or a character vector of column names (the first row is assumed to not be column names, and the values of the character vector will be the resulting column names in the returned \code{GRanges} object). If there are columns in addition to the BED6 columns, it is recommended that a character vector is supplied so that columns can be referred to by downstream functions. Default \code{FALSE}. See documentation for \code{readr::read_tsv()} for more details.
#' @param genome Gives the genome assembly (human or mouse). Must be one of 'hg19', 'hg38', 'mm9' or 'mm10'
#' @param stranded Logical variable. If TRUE, strand attribute of GenomicRanges object is drawn from 6th column of BED file. Default \code{FALSE}
#' @param use.score Logical variable. If TRUE, score attribute of GenomicRanges object is drawn from 5th column of BED file. Default \code{FALSE}
#'
#' @return A GenomicRanges object with ranges limited by genome and BED file. The GenomicRanges object is sorted if it is detected to be unsorted, and the regions are unique. The name column (4th) in the BED file is \code{name} attribute and the score column (5th) in the BED file is the \code{score} attribute in the returned GenomicRanges object.
#'
#' @examples
#' file = system.file('extdata', 'K562_Cjun.narrowPeak.gz', package = 'annotatr')
#' read_bed(file = file, genome = 'hg19', stranded = FALSE, use.score = FALSE)
#'
#' @export
read_bed <- function(file, col.names=FALSE, genome, stranded = FALSE, use.score = FALSE){
  # Error checking pre-read
    if (!file.exists(file)){
      stop(sprintf('Error: File, %s, not found.', file))
    }
    if (! genome %in% c('hg19','hg38','mm9','mm10')){
      stop('Error: Invalid genome.')
    }

  # Read
  bed <- readr::read_tsv(file = file, col_names = col.names)

  # Error checking post-read
    if (!all(grepl("chr", bed[[1]]))){
      stop('Error: First column of BED file does not appear to be chromsome.')
    }
    if (typeof(head(bed[[2]])) != "integer"){
        stop('Error: Second column of BED file must be integer valued.')
    }
    if (typeof(head(bed[[3]]))!= "integer"){
        stop('Error: Third column of BED file must be integer valued.')
    }

  # Retrieve chromosome sizes for the genome
  size_code = sprintf('%s_chrom_sizes', genome)
  data(list = size_code, package = "annotatr")
  seqlengths = get(size_code)

  # Construct the appropriate strand vector
  if(stranded){
    if(length(base::setdiff( unique(bed[[6]]), c("+","-") ) > 0)) {
      stop("Error: When stranded = T, strand column should contain +/- only.")
    }
    strand = bed[[6]]
  } else {
    strand = rep.int('*', nrow(bed))
  }

  # Construct the appropriate mcols data.frame
  if(use.score && ncol(bed) == 6) {
    mcols = data.frame(name = bed[[4]], score = bed[[5]], stringsAsFactors=F)
  } else if (!use.score && ncol(bed) == 6) {
    mcols = data.frame(name = bed[[4]], stringsAsFactors=F)
  } else if (!use.score && ncol(bed) == 3) {
    mcols = NULL
  } else if (use.score && ncol(bed) > 6) {
    mcols = data.frame(name = bed[[4]], bed[c(5,7:ncol(bed))], stringsAsFactors=F)
  }

  if(!is.null(mcols)) {
    gr <- GenomicRanges::GRanges(
        seqnames = bed[[1]],
        ranges = IRanges::IRanges(start = bed[[2]], end = bed[[3]]),
        strand = strand,
        seqlengths = seqlengths)
    GenomicRanges::mcols(gr) = mcols
  } else {
    gr <- GenomicRanges::GRanges(
        seqnames = bed[[1]],
        ranges = IRanges::IRanges(start = bed[[2]], end = bed[[3]]),
        strand = strand,
        seqlengths = seqlengths)
  }

  # Assign the genome metadata
  GenomeInfoDb::genome(gr) = genome

  # Ensure the ranges do not exceed the chromosome lengths for the genome
  gr <- GenomicRanges::trim(gr)

  # Check if gr is sorted and sort it if it isn't
  # Not sure what the story with this warning is when we do is.unsorted()
  # Might have to do with the is.circular column being NA...
  # "In is.na(x) : is.na() applied to non-(list or vector) of type 'S4'"
  if(suppressWarnings(is.unsorted(gr))) {
      # NOTE: The "natural order" for the elements of a GenomicRanges object is to order them (a) first by sequence level, (b) then by strand, (c) then by start, (d) and finally by width. This way, the space of genomic ranges is totally ordered.
      gr <- sort(gr)
  }

  # Enforce uniqueness of regions
  gr <- unique(gr)

  gr
}
