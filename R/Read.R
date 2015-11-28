library(readr)
library(GenomicRanges)
library(IRanges)

read.bed <- function(filename, genome, stranded = F){
    if (!file.exists(filename)){
        stop(paste("Error in read.bed(filename, genome, stranded): filename",
                   filename, " not found"))
    }

    bed <- read.table(filename, sep = "\t", header = F, stringsAsFactors = F)

    if (!grepl("chr", bed[,1])){
        stop("Error in read.bed(filename, genome, stranded): file not in
        correct format, missing chr numbers")
    }

    if (typeof(head(bed[,2])) != "numeric" | typeof(head(bed[,2])) != "integer"){
        stop("Error in read.bed(filename, genome, stranded): file not in
             correct format, second column not numeric")
    }

    if (typeof(head(bed[,3])) != "numeric" | typeof(head(bed[,3]))!= "integer"){
        stop("Error in read.bed(filename, genome, stranded): file not in
             correct format, third column not numeric")
    }

    if (stranded == T & (bed[1,6] != "+" | bed[1,6] != "-")){
        stop("Error in read.bed(filename, genome, stranded): strandedness
             column not in correct format")
    }

    gR <- GenomicRanges::GRanges(seqnames = bed[,1],
                      ranges = IRanges::IRanges(
                          start = bed[,2],
                            end = bed[,3])
    )

    if (stranded == TRUE){
        gR.strand = bed[,6]
    }

    genome
}




