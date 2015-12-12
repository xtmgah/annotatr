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

#' Function to tidy up annotation accessors for visualization
#'
#' @param annotations A character vector of annotations, in the order they are to appear in the visualization.
#'
#' @return A character vector of tidied annotation names
#'
#' @examples
#' tidy_annotations(annotations = c('hg19_cpg_islands','hg19_knownGenes_promoters'))
tidy_annotations = function(annotations) {
  tidy = sapply(annotations, function(a){
    tokens = unlist(strsplit(a,'_'))
    if(tokens[2] == 'cpg') {
      if(tokens[3] == 'inter') {
        return('interCGI')
      } else {
        return(paste('CpG', tokens[3]))
      }
    } else if (tokens[2] == 'knownGenes') {
      return(tokens[3])
    }
  })

  return(tidy)
}
