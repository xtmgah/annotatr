#' @export
read.bed <- function(filename, genome, stranded = F){
    if (!file.exists(filename)){
        stop(paste("In read.bed(filename, genome, stranded): file:",
            filename, "not found."))
    }

    if (! genome %in% c('hg19','hg38','mm9','mm10')){
        stop("In read.bed(filename, genome, stranded): Invalid genome.
            Valid genomes are: hg19, hg38, mm9, and mm10.")
    }

    # Read the file
    bed <- read.table(filename, sep = "\t", header = F, stringsAsFactors = F)

    if (!all(grepl("chr", bed[,1]))){
        stop("In read.bed(filename, genome, stranded): file not in
            correct format, chromosomes expected to start with 'chr'.")
    }

    if (typeof(head(bed[,2])) != "integer"){
        stop("In read.bed(filename, genome, stranded): file not in
             correct format, start positions must be integers.")
    }

    if (typeof(head(bed[,3]))!= "integer"){
        stop("In read.bed(filename, genome, stranded): file not in
             correct format, end positions must be integers.")
    }

    if (stranded == T & (bed[1,6] != "+" | bed[1,6] != "-")){
        stop("In read.bed(filename, genome, stranded): with stranded = T,
          strand column should contain +/- only.")
    }

    # Get chromosome sizes from data/
    size_code = sprintf('%s_chrom_sizes', genome)
    data(list = size_code, package = "annotatr")
    seqlengths = get(size_code)

    if(stranded){
        gR <- GenomicRanges::GRanges(
            seqnames = bed[,1],
            ranges = IRanges::IRanges(
            start = bed[,2],
            end = bed[,3]),
            strand = bed[,6],
            type = bed$type,
            seqlengths = seqlengths)
    } else {
        gR <- GenomicRanges::GRanges(seqnames = bed[,1],
            ranges = IRanges::IRanges(
            start = bed[,2],
            end = bed[,3]),
            strand = '*',
            type = bed$type,
            seqlengths = seqlengths)
    }

    gR <- GenomicRanges::trim(gR)
    gR
}
