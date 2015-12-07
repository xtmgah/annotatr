#' Convert a bed file into a GenomicRanges object
#'
#' \code{read_bed} reads in data from a BED6 file, checks to see if it follows the appropriate BED (Browser Extensible Data) format, and checks validity of genome. If valid, reads file data into a GenomicRanges object. Unlike in the BED specification, the score column can be continuously valued. The names need not be unique.
#'
#' @param filename Path to the file with data. File must exist and be in correct format.
#' @param genome Gives the genome assembly (human or mouse). Must be one of 'hg19', 'hg38', 'mm9' or 'mm10'
#' @param stranded Logical variable. If TRUE, strand attribute of GenomicRanges object is drawn from 6th column of BED file.
#' @param use.score Logical variable. If TRUE, score attribute of GenomicRanges object is drawn from 5th column of BED file.
#'
#' @return A GenomicRanges object with ranges limited by genome and BED file. The GenomicRanges object is sorted if it is detected to be unsorted, and the regions are unique. The name column (4th) in the BED file is \code{regionName} attribute and the score column (5th) in the BED file is the \code{score} attribute in the returned GenomicRanges object.
#'
#' @examples
#' file = system.file('extdata', 'K562_Cjun.narrowPeak.gz', package = 'annotatr')
#' read_bed(filename = file, genome = 'hg19', stranded = FALSE, use.score = FALSE)
#'
#' @export
read_bed <- function(filename, genome, stranded = FALSE, use.score = FALSE){
    if (!file.exists(filename)){
        stop(paste("In read_bed(filename, genome, stranded): filename",
                   filename, " not found"))
    }

    if (! genome %in% c('hg19','hg38','mm9','mm10')){
        stop("in read_bed(filename, genome, stranded): Invalid Genome")
    }

    bed <- readr::read_tsv(file = filename, col_names = FALSE)

    if (!all(grepl("chr", bed[[1]]))){
        stop("in read_bed(filename, genome, stranded): file not in
             correct format, missing chr numbers")
    }

    if (typeof(head(bed[[2]])) != "integer"){
        stop("in read_bed(filename, genome, stranded): file not in
             correct format, second column not integer")
    }

    if (typeof(head(bed[[3]]))!= "integer"){
        stop("in read_bed(filename, genome, stranded): file not in
             correct format, third column not integer")
    }

    size_code = sprintf('%s_chrom_sizes', genome)
    data(list = size_code, package = "annotatr")
    seqlengths = get(size_code)

    if(stranded){
        if(length(base::setdiff( unique(bed[[6]]), c("+","-") ) > 0)) {
            stop("In read_bed(filename, genome, stranded): with stranded = T,
                 strand column should contain +/- only.")
        }
        if(use.score) {
            gR <- GenomicRanges::GRanges(
                seqnames = bed[[1]],
                ranges = IRanges::IRanges(start = bed[[2]], end = bed[[3]]),
                strand = bed[[6]],
                regionName = bed[[4]],
                score = bed[[5]],
                seqlengths = seqlengths)
        } else{
            gR <- GenomicRanges::GRanges(
                seqnames = bed[[1]],
                ranges = IRanges::IRanges(start = bed[[2]], end = bed[[3]]),
                strand = bed[[6]],
                regionName = bed[[4]],
                seqlengths = seqlengths)
        }
    } else {
        if(use.score) {
            gR <- GenomicRanges::GRanges(
                seqnames = bed[[1]],
                ranges = IRanges::IRanges(start = bed[[2]], end = bed[[3]]),
                strand = '*',
                regionName = bed[[4]],
                score = bed[[5]],
                seqlengths = seqlengths)
        } else{
            gR <- GenomicRanges::GRanges(
                seqnames = bed[[1]],
                ranges = IRanges::IRanges(start = bed[[2]], end = bed[[3]]),
                strand = '*',
                regionName = bed[[4]],
                seqlengths = seqlengths)
        }
    }
    # Assign the genome metadata
    GenomeInfoDb::genome(gR) = genome

    # Ensure the ranges do not exceed the chromosome lengths for the genome
    gR <- GenomicRanges::trim(gR)

    # Check if gR is sorted and sort it if it isn't
    # Not sure what the story with this warning is when we do is.unsorted()
    # Might have to do with the is.circular column being NA...
    # "In is.na(x) : is.na() applied to non-(list or vector) of type 'S4'"
    if(suppressWarnings(is.unsorted(gR))) {
        # NOTE: The "natural order" for the elements of a GenomicRanges object is to order them (a) first by sequence level, (b) then by strand, (c) then by start, (d) and finally by width. This way, the space of genomic ranges is totally ordered.
        gR <- sort(gR)
    }

    # Enforce uniqueness of regions
    gR <- unique(gR)

    gR
}
