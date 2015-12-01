# Helper functions

#' Function to inform the user which annotations are available.
#'
#' @return A character vector of accessors for use with data() to load annotations.
#'
#' @examples
#' supported_annotations()
#' data(hg19_cpg_islands)
#' data(hg38_knownGenes_promoters)
#'
#' @export
supported_annotations = function() {
  annots = data(package='annotatr')[['results']][,3]
  annots = annots[grepl('cpg|knownGenes', annots)]
  return(annots)
}

#' Function to inform the user which genomes are supported
#'
#' @return A character vector of genomes for which we have annotations
#'
#' @examples
#' supported_genomes()
#'
#' @export
supported_genomes = function() {
  return(c('hg19','hg38','mm9','mm10'))
}
