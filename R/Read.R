read.bed <- function(filename, genome, stranded = F){
    if (!file.exists(filename)){
        stop(paste("In read.bed(filename, genome, stranded): filename",
                   filename, " not found"))
    }

    if (! genome in c('hg19','hg38','mm9','mm10')){
        stop("in read.bed(filename, genome, stranded): Invalid Genome")
    }

    bed <- readr::read.table(filename, sep = "\t", header = F, stringsAsFactors = F)

    if (!grepl("chr", bed[,1])){
        stop("in read.bed(filename, genome, stranded): file not in
        correct format, missing chr numbers")
    }

    if (typeof(head(bed[,2])) != "integer"){
        stop("in read.bed(filename, genome, stranded): file not in
             correct format, second column not numeric")
    }

    if (typeof(head(bed[,3]))!= "integer"){
        stop("in read.bed(filename, genome, stranded): file not in
             correct format, third column not numeric")
    }

    if (stranded == T & (bed[1,6] != "+" | bed[1,6] != "-")){
        stop("in read.bed(filename, genome, stranded): strandedness
             column not in correct format")
    }

    size_code = sprintf('%s_chrom_sizes', genome)
    data(list = size_code, package = "annotatr")
    seqlengths = get(size_code)

    if(stranded){
        gR <- GenomicRanges::GRanges(seqnames = bed[,1],
                      ranges = IRanges::IRanges(
                          start = bed[,2],
                            end = bed[,3]),
                            strand = bed[,6],
                            type = bed$type,
                            seqlengths = seqlengths
    }

    else {
        gR <- GenomicRanges::GRanges(seqnames = bed[,1],
                                     ranges = IRanges::IRanges(
                                         start = bed[,2],
                                         end = bed[,3]),
                                        strand = '*',
                                     type = bed$type,
                                     seqlengths = seqlengths
    }

    gR <- genomicRanges::trim(gR)
}




