#' Convert a bed file into a GenomicRanges object
#'
#' \code{read_bed} Reads in data table from file, checks to see if it follows the appropriate
#' BED (Browser Extensible Data) format and if genome is valid.
#' If valid, reads file data into a GenomicRanges object.
#'
#' @param filename Path to the file with data. File must exist and be in correct format.
#' @param genome Gives the genome assembly (human or mouse). Must be one of 'hg19', 'hg38', 'mm9' or 'mm10'
#' @param stranded Logical variable. If TRUE, strand attribute of GenomicRanges object is drawn from 6th column of BED file.
#'
#' @return A GenomicRanges object with ranges determined by genome and bed file. The GenomicRanges object is sorted if it is detected to be unsorted, and the regions are unique.
#'
#' @examples
#' file = system.file('extdata', 'K562_Cjun.narrowPeak.gz', package = 'annotatr')
#' read_bed(filename = file, genome = 'hg19', stranded = FALSE)
#'
#' @export
read_bed <- function(filename, genome, stranded = FALSE){
    if (!file.exists(filename)){
        stop(paste("In read_bed(filename, genome, stranded): filename",
                   filename, " not found"))
    }

    if (! genome %in% c('hg19','hg38','mm9','mm10')){
        stop("in read_bed(filename, genome, stranded): Invalid Genome")
    }

    bed <- read.table(filename, sep = "\t", header = FALSE, stringsAsFactors = FALSE)

    if (!all(grepl("chr", bed[,1]))){
        stop("in read_bed(filename, genome, stranded): file not in
             correct format, missing chr numbers")
    }

    if (typeof(head(bed[,2])) != "integer"){
        stop("in read_bed(filename, genome, stranded): file not in
             correct format, second column not numeric")
    }

    if (typeof(head(bed[,3]))!= "integer"){
        stop("in read_bed(filename, genome, stranded): file not in
             correct format, third column not numeric")
    }

    size_code = sprintf('%s_chrom_sizes', genome)
    data(list = size_code, package = "annotatr")
    seqlengths = get(size_code)

    if(stranded){
        if(bed[1,6] != "+" | bed[1,6] != "-") {
            stop("In read_bed(filename, genome, stranded): with stranded = T,
                 strand column should contain +/- only.")
        }

        gR <- GenomicRanges::GRanges(
          seqnames = bed[,1],
          ranges = IRanges::IRanges(start = bed[,2], end = bed[,3]),
          strand = bed[,6],
          type = bed$type,
          seqlengths = seqlengths)
        GenomeInfoDb::genome(gR) = genome
    } else {
        gR <- GenomicRanges::GRanges(
            seqnames = bed[,1],
            ranges = IRanges::IRanges(start = bed[,2], end = bed[,3]),
            strand = '*',
            type = bed$type,
            seqlengths = seqlengths)
        GenomeInfoDb::genome(gR) = genome
    }

    gR <- GenomicRanges::trim(gR)

    # Check if gR is sorted and sort it if it isn't
    # Not sure what the story with this warning is when we do is.unsorted()
    # Might have to do with the is.circular column being NA...
    # "In is.na(x) : is.na() applied to non-(list or vector) of type 'S4'"
    if(suppressWarnings(is.unsorted(gR))) {
      gR <- sort(gR)
    }

    # Enforce uniqueness of regions
    gR <- unique(gR)

    gR
}
