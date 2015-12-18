#' A function to intersect user data with annotation data
#'
#' Given a GRanges object constructed from user supplied data with read_bed()(),
#' and desired annotation overlaps, use GenomicRanges::findOverlaps() to return a
#' Hits object indicating which elements of the user-data (query) intersect which
#' elements of the annotation data (subject).
#'
#' @param regions The GRanges object of the user-data returned by read_bed().
#' @param annotations A character vector of annotations to overlap with user-data. Valid annotation codes can be found with supported_annotations(). The "basic_genes" shortcut annotates regions to the 1-5Kb, promoter, 5UTR, exon, intron, and 3UTR knownGene regions. The "detailed_genes" shortcut annotates regions to the 1-5Kb, promoter, 5UTR exon/intron, CDS exon/intron, and 3UTR exon/intron knownGene regions. The "cpgs" shortcut annotates regions to the CpG islands, shores, shelves, and interCGI regions. NOTE: basic_genes and detailed_genes annotations cannot be done at the same time.
#' @param genome One of the genomes in supported_genomes(). Should match argument used in read_bed()().
#' @param ignore.strand Logical variable indicating whether strandedness should be respected in findOverlaps(). Default FALSE.
#' @param use.score Logical variable. Include the "score" for each genomic region in the tabulated results. Score can mean a variety of things, e.g. percent methylation of a CpG/region or fold-change of a ChIP-seq peak.
#'
#' @return A \code{dplyr::tbl_df} with columns from the GenomicRanges object for the regions and corresponding annotations.
#'
#' @examples
#' # A very simple example with only 3 genomic regions
#' bed = system.file('extdata', 'test_intersect.bed', package = 'annotatr')
#' annotations = c('hg19_cpg_islands','hg19_cpg_shores','hg19_knownGenes_promoters')
#'
#' d = read_bed(file = bed, genome = 'hg19', stranded = FALSE)
#'
#' i = annotate_regions(
#'   regions = d,
#'   annotations = annotations,
#'   genome = 'hg19',
#'   ignore.strand = TRUE,
#'   use.score = FALSE)
#'
#' # A more complicated example using Gm12878 Pol2 ChIP-seq from ENCODE and an annotation shortcut
#' bed = system.file('extdata', 'Gm12878_Pol2.narrowPeak.gz', package = 'annotatr')
#' annotations = c('basic_genes')
#'
#' d = read_bed(file = bed, genome = 'hg19', stranded = FALSE)
#'
#' i = annotate_regions(
#'   regions = d,
#'   annotations = annotations,
#'   genome = 'hg19',
#'   ignore.strand = TRUE,
#'   use.score = FALSE)
#'
#' @export
annotate_regions = function(regions, annotations, genome, ignore.strand = TRUE, use.score) {
  # Checks before moving forward
  if(class(regions)[1] != "GRanges") {
    stop('Error in annotate_regions(...): regions object is not GRanges.')
  }

  if(!(genome %in% supported_genomes())) {
    stop('Error in annotate_regions(...): unsupported genome given. See supported_genomes().')
  }

  # Check that the annotations are supported, tell the user which are unsupported
  if(!all(annotations %in% c(supported_annotations(),'cpgs','basic_genes','detailed_genes'))) {
    unsupported = base::setdiff(annotations, c(supported_annotations(),'cpgs','basic_genes','detailed_genes'))

    stop(sprintf('Error in annotate_regions(...): "%s" is(are) not supported. See supported_annotations().',
      paste(unsupported, collapse=', ')))
  }

  # Do not allow basic_genes and detailed_genes at the same time
  if('basic_genes' %in% annotations && 'detailed_genes' %in% annotations) {
    stop('Error in annotate_regions(...): please choose between basic_genes and detailed_genes annotations.')
  }

  if('basic_genes' %in% annotations || 'detailed_genes' %in% annotations || 'cpgs' %in% annotations) {
    # Check for shortcut annotation accessors 'cpgs', 'basic_genes', or 'detailed_genes'
    # and create the right annotations based on the genome
    new_annotations = c()
    if('cpgs' %in% annotations) {
      new_annotations = paste(genome, 'cpg', c('islands','shores','shelves','inter'), sep='_')
    }
    if ('basic_genes' %in% annotations) {
      new_annotations = c(new_annotations, paste(genome, 'knownGenes', c('1to5kb','promoters','5UTRs','exons','introns','3UTRs'), sep='_'))
    }
    if ('detailed_genes' %in% annotations) {
      new_annotations = c(new_annotations, paste(genome, 'knownGenes',
        c('1to5kb','promoters','exons5UTRs','introns5UTRs','exonsCDSs','intronsCDSs','exons3UTRs','introns3UTRs'), sep='_'))
    }
    annotations = setdiff(c(annotations, new_annotations), c('basic_genes','detailed_genes','cpgs'))
  }

  # Collect the annotation objects into a GRangesList
  data(list = annotations, package = 'annotatr')

  # Perform the intersections in an lapply (consider using mclapply)
  intersections = lapply(annotations, function(annot){
    message(sprintf('Intersecting %s', annot))

    GenomicRanges::findOverlaps(regions, get(annot), ignore.strand = ignore.strand)
  })
  names(intersections) = annotations

  # Pull in data from the regions and the annotations
  # Think about using mclapply instead of lapply to speed up. May cause some
  # memory issues or overhead.
  if(use.score) {
    tab_list = lapply(names(intersections), function(n){
      message(sprintf('Annotating %s with score', n))

      # Get subsets of
      r_sub = regions[intersections[[n]]@queryHits]
      a_sub = get(n)[intersections[[n]]@subjectHits]

      # Create the dplyr::tbl_df with score
      df_ra = dplyr::tbl_df(data.frame(
        annot_start = start(a_sub),
        annot_end = end(a_sub),
        annot_strand = as.character(strand(a_sub)),
        annot_type = n,
        annot_id = a_sub$ID,
        data_chrom = as.character(seqnames(r_sub)),
        data_start = start(r_sub),
        data_end = end(r_sub),
        data_strand = as.character(strand(r_sub)),
        stringsAsFactors = F))

      df_d = dplyr::tbl_df(as.data.frame(GenomicRanges::mcols(r_sub)))

      df = dplyr::bind_cols(df_ra, df_d)
    })
  } else {
    tab_list = lapply(names(intersections), function(n){
      message(sprintf('Annotating %s without score', n))

      # Get subsets of
      r_sub = regions[intersections[[n]]@queryHits]
      a_sub = get(n)[intersections[[n]]@subjectHits]

      # Create the dplyr::tbl_df without score
      df = dplyr::tbl_df(data.frame(
        annot_start = start(a_sub),
        annot_end = end(a_sub),
        annot_strand = as.character(strand(a_sub)),
        annot_type = n,
        annot_id = a_sub$ID,
        data_chrom = as.character(seqnames(r_sub)),
        data_start = start(r_sub),
        data_end = end(r_sub),
        data_strand = as.character(strand(r_sub)),
        name = r_sub$name,
        stringsAsFactors=F
      ))
    })
  }

  # Combine and sort the list of data.frames into a single data.frame
  message('Combining annotations')
  tab = dplyr::bind_rows(tab_list)
  message('Sorting annotations')
  tab = dplyr::arrange(tab, data_chrom, data_start, data_end, annot_start)

  return(tab)
}
