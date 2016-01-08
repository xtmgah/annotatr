#' Convert a BED file into a GenomicRanges object
#'
#' \code{read_bed} reads in data from a BED6 file, checks to see if it follows the appropriate BED (Browser Extensible Data) format, and checks validity of genome. If valid, reads file data into a GenomicRanges object. Unlike in the BED specification, the score column can be continuously valued. The names need not be unique.
#'
#' @param file Path to the file with data. File must exist and be in correct format.
#' @param col.names Either \code{TRUE} (column names are given in the first row of the file), \code{FALSE} (the first row does not have column names), or a character vector of column names (the first row is assumed to not be column names, and the values of the character vector will be the resulting column names in the returned \code{GRanges} object). If \code{FALSE} and the file is not BED6, a warning will indicate that downstream functions may behave unexpectedly.
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
#' file = system.file('extdata', 'test_read_multiple_data_nohead.bed', package = 'annotatr')
#' read_bed(file = file,
#'   col.names = c('chrom','start','end','name','score','strand','pval','mu1','mu0','diff_exp'),
#'   genome = 'hg19', stranded = FALSE, use.score = TRUE)
#'
#' @export
read_bed = function(file, col.names=FALSE, genome, stranded = FALSE, use.score = FALSE) {
  # Error checking pre-read
    if(!file.exists(file)) {
      stop(sprintf('Error: File, %s, not found.', file))
    }
    if(! genome %in% supported_genomes()) {
      warning('Warning: %s is not a supported genome. In order to annotate regions, make sure to load custom annotations with read_annotations()', genome)
    }

  # Read
  bed <- readr::read_tsv(file = file, col_names = col.names)

  # Error checking post-read
    if(!all(grepl("chr", bed[[1]]))) {
      stop('Error: First column of BED file does not appear to be chromosome.')
    }
    if(typeof(head(bed[[2]])) != "integer") {
        stop('Error: Second column of BED file must be integer valued.')
    }
    if(typeof(head(bed[[3]]))!= "integer") {
        stop('Error: Third column of BED file must be integer valued.')
    }

  # Deal with column names
    if(class(col.names) != 'character' && col.names == FALSE) {
      if(ncol(bed) == 6) {
        colnames(bed) = c('chr','start','end','name','score','strand')
      } else {
        warning('Warning: Input file is not BED6, and no column names were given. Downstream functions may behave oddly.')
      }
    }

  # Construct the appropriate strand vector
  if(stranded) {
    if(length(base::setdiff( unique(bed[[6]]), c("+","-") ) > 0)) {
      stop("Error: When stranded = T, strand column should contain +/- only.")
    }
    strand = bed[[6]]
  } else {
    strand = rep.int('*', nrow(bed))
  }

  # Construct the appropriate mcols data.frame
  if(use.score && ncol(bed) == 6) {
    # This is the vanilla use.score case. This overwrites column name information.
    mcols = data.frame(name = bed[[4]], score = bed[[5]], stringsAsFactors=F)
  } else if(!use.score && ncol(bed) == 6) {
    # This is the !use.score case.
    mcols = data.frame(name = bed[[4]], stringsAsFactors=F)
  } else if(!use.score && ncol(bed) < 6) {
    # Shrink to BED3 case
    mcols = NULL
  } else if(use.score && ncol(bed) > 6) {
    # This is the multiple data column case. name column remains so named.
    # We are expecting the user to name their columns, otherwise we
    # won't possibly know what to do in summarize and visualize.
    mcols = data.frame(bed[c(4:5,7:ncol(bed))], stringsAsFactors=F)
  } else if(!use.score && ncol(bed) > 6) {
    # This is the multiple data column case. name column remains so named.
    # It is conceivable that the user may just tack on data columns to an
    # existing BED file, hence skipping over column 5 (score)
    # if it's all 1000s like in a narrowPeak.
    mcols = data.frame(bed[c(4,7:ncol(bed))], stringsAsFactors=F)
  } else {
    # Shrink to BED3 case
    mcols = NULL
  }

  # Construct the GRanges object
  if(!is.null(mcols)) {
    gr <- GenomicRanges::GRanges(
        seqnames = bed[[1]],
        ranges = IRanges::IRanges(start = bed[[2]], end = bed[[3]]),
        strand = strand)
    GenomicRanges::mcols(gr) = mcols
  } else {
    gr <- GenomicRanges::GRanges(
        seqnames = bed[[1]],
        ranges = IRanges::IRanges(start = bed[[2]], end = bed[[3]]),
        strand = strand)
  }

  if(genome %in% supported_genomes()) {
    # Retrieve chromosome sizes for the genome
    size_code = sprintf('%s_chrom_sizes', genome)
    data(list = size_code, package = "annotatr")
    seqlengths = get(size_code)

    GenomeInfoDb::seqlengths(gr) = seqlengths[GenomeInfoDb::seqlevelsInUse(gr)]
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

#' Convert a BED file into annotations
#'
#' Annotations are genomic regions given as 3, 4, or 5 columned tab-delimited files (chrom, start, end, name, strand). If the input is 3-columned, the annotations are assumed to be unstranded and unnamed. In this case, the strand is set to * and the names are set to the numbered sequence \code{annotation_name:n}. If the input is 4-columned, the content of the name column is used, and the strand is set to *. If the input is 5-columned, the name and strand columns are used. If the annotations are unstranded, use '*' for the strand column.
#'
#' If the names for the annotations are not unique, they will be renamed in the form \code{annotation_name:n}.
#'
#' @param file Path to the file with annotations. This file can be either BED3, BED4, or BED5. If BED3, the assumption is that the annotations are not stranded and not named; if BED4,
#' @param genome Gives the genome assembly. The genome may be one of hg19, hg38, mm9, or mm10, in which case the user can also annotate regions with built-in annotations. Otherwise, the user may only annotate regions with their custom annotations.
#' @param annotation_name A description of the types of annotations present in \code{file}.
#'
#' @return A GenomicRanges object with the name \code{[genome]_custom_[annotation_name]} for use in the \code{annotations} parameter in \code{annotate_regions()} along with built-in annotations.
#'
#' @examples
#' file = system.file('extdata', 'test_annotations_5.bed', package='annotatr')
#' hg19_custom_test = read_annotations(file = file, genome = 'hg19', annotation_name = 'test')
#'
#' @export
read_annotations = function(file, genome, annotation_name) {
  # Error checking pre-read
    if(!file.exists(file)) {
      stop(sprintf('Error: File, %s, not found.', file))
    }
    if(is.null(genome)) {
      stop('Error: A genome must be specified.')
    }

  # Read
  rename = FALSE
  bed <- readr::read_tsv(file = file, col_names = FALSE)

  # Error checking post-read
    if(!all(grepl('chr', bed[[1]]))) {
      stop('Error: First column of annotation file does not appear to be chromosome.')
    }
    if(class(bed[[2]]) != 'integer') {
      stop('Error: Second column of annotation file must be integer valued.')
    }
    if(class(bed[[3]]) != 'integer') {
      stop('Error: Third column of annotation file must be integer valued.')
    }
    if(ncol(bed) > 3) {
      if(any(duplicated(bed[[4]]))) {
        rename = TRUE
        warning('Warning: Some annotations have duplicated names, annotations will be renamed to be unique.')
      }
    }
    if(ncol(bed) == 5) {
      if(length(base::setdiff( unique(bed[[5]]), c("+","-","*") ) > 0)) {
        stop('Error: Strand column should contain only "+", "-", or "*"')
      }
    }

  # Construct the GenomicRanges object
  if(ncol(bed) == 3) {
    gr <- GenomicRanges::GRanges(
        seqnames = bed[[1]],
        ranges = IRanges::IRanges(start = bed[[2]], end = bed[[3]]),
        strand = '*',
        ID = sprintf('%s:%s', annotation_name, 1:nrow(bed)))
  } else if(ncol(bed) == 4) {
    if(rename) {
      gr <- GenomicRanges::GRanges(
          seqnames = bed[[1]],
          ranges = IRanges::IRanges(start = bed[[2]], end = bed[[3]]),
          strand = '*',
          ID = sprintf('%s:%s', annotation_name, 1:nrow(bed)))
    } else {
      gr <- GenomicRanges::GRanges(
          seqnames = bed[[1]],
          ranges = IRanges::IRanges(start = bed[[2]], end = bed[[3]]),
          strand = '*',
          ID = bed[[4]])
    }
  } else if(ncol(bed) == 5) {
    if(rename) {
      gr <- GenomicRanges::GRanges(
          seqnames = bed[[1]],
          ranges = IRanges::IRanges(start = bed[[2]], end = bed[[3]]),
          strand = bed[[5]],
          ID = sprintf('%s:%s', annotation_name, 1:nrow(bed)))
    } else {
      gr <- GenomicRanges::GRanges(
          seqnames = bed[[1]],
          ranges = IRanges::IRanges(start = bed[[2]], end = bed[[3]]),
          strand = bed[[5]],
          ID = bed[[4]])
    }
  } else {
    stop('Error: The file used to create custom annotations must have 3, 4, or 5 columns.')
  }

  # Give the GenomicRanges object the genome annotation
  GenomeInfoDb::genome(gr) = genome

  # Add relevant chromosome information if a supported genome is given
  if(genome %in% supported_genomes()) {
    # Retrieve chromosome sizes for the genome
    size_code = sprintf('%s_chrom_sizes', genome)
    data(list = size_code, package = "annotatr")
    seqlengths = get(size_code)

    GenomeInfoDb::seqlengths(gr) = seqlengths[GenomeInfoDb::seqlevelsInUse(gr)]
  } else {
    warning(sprintf('Warning: The %s genome is not supported, and cannot be used with any built-in annotations', genome))
  }

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

  return(gr)
}
