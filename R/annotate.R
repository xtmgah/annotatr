# Functions relating to using the intersections to annotate regions to annotations

#' A function to create tabular output of regions intersected with annotations
#'
#' Given output from \code{intersect_annotations()}, create a table for output which shows peak information and corresponding annotation information.
#'
#' @param regions A GenomicRanges object output from \code{read_bed()}.
#' @param intersections A list of Hits objects from \code{intersect_annotations()}.
#' @param use.score Logical variable. Include the "score" for each genomic region in the tabulated results. Score can mean a variety of things, e.g. percent methylation of a CpG/region or fold-change of a ChIP-seq peak.
#'
#' @return A table with columns from the GenomicRanges object for the regions and corresponding annotations
#'
#' @examples
#' bed = system.file('extdata', 'Gm12878_Pol2.narrowPeak.gz', package = 'annotatr')
#' annotations = c('basic_genes')
#'
#' d = read_bed(filename = bed, genome = 'hg19', stranded = FALSE, use.score = FALSE)
#'
#' i = intersect_annotations(
#'   regions = d,
#'   annotations = annotations,
#'   genome = 'hg19',
#'   ignore.strand = TRUE)
#'
#' t = tabulate_intersections(
#'   regions = d,
#'   intersections = i,
#'   use.score = FALSE)
#'
#' @export
tabulate_intersections = function(regions, intersections, use.score = FALSE) {
  # Checks before moving forward
  if(class(regions)[1] != "GRanges") {
    stop('Error in tabulate_intersections(...): regions object is not GRanges.')
  }

  if(class(intersections) != "list") {
    stop('Error in tabulate_intersections(...): intersections must be a list.')
  }

  if(unique(sapply(intersections, class)) != "Hits") {
    stop('Error in tabulate_intersections(...): intersections must be a list of Hits objects.')
  }

  # In the future, we could think about changing this to mclapply to speed it up
  tab_list = lapply(names(intersections), function(n){
    message(sprintf('Tabulating %s', n))

    # Get subsets of
    r_sub = regions[intersections[[n]]@queryHits]
    a_sub = get(n)[intersections[[n]]@subjectHits]

    if(!use.score) {
      # Create the data.frame
      df = data.frame(
        data_chrom = as.character(seqnames(r_sub)),
        data_start = start(r_sub),
        data_end = end(r_sub),
        data_strand = as.character(strand(r_sub)),
        annot_chrom = as.character(seqnames(a_sub)),
        annot_start = start(a_sub),
        annot_end = end(a_sub),
        annot_strand = as.character(strand(a_sub)),
        annot_type = n,
        annot_id = a_sub$ID,
        stringsAsFactors=F
      )
    } else {
      # Create the data.frame
      df = data.frame(
        data_chrom = as.character(seqnames(r_sub)),
        data_start = start(r_sub),
        data_end = end(r_sub),
        data_strand = as.character(strand(r_sub)),
        data_score = r_sub$score,
        annot_chrom = as.character(seqnames(a_sub)),
        annot_start = start(a_sub),
        annot_end = end(a_sub),
        annot_strand = as.character(strand(a_sub)),
        annot_type = n,
        annot_id = a_sub$ID,
        stringsAsFactors=F
      )
    }
  })

  # Combine and sort the list of data.frames into a single data.frame
  tab = Reduce(rbind, tab_list)
  tab = tab[order(tab$data_chrom, tab$data_start, tab$data_end, tab$annot_start, decreasing=F),]

  return(tab)
}
