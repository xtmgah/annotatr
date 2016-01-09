# Helper functions

#' Function listing which annotations are available.
#'
#' This includes the shortcuts, which cannot be used in the data(...) calls. To call a shortcut, need to use \code{expand_annotations()}.
#'
#' @return A character vector of available annotations.
#'
#' @examples
#' supported_annotations()
#' data(hg19_cpg_islands)
#' data(hg38_knownGenes_promoters)
#'
#' @export
supported_annotations = function() {
  shortcuts = c('basicgenes','detailedgenes','cpgs')

  annots = data(package='annotatr')[['results']][,3]
  annots = annots[grepl('cpg|knownGenes|enhancers', annots)]
  annots = c(annots, apply(expand.grid(supported_genomes(), shortcuts, stringsAsFactors=F), 1, paste, collapse='_'))
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
#' @return A list of mappings from original annotation names to names ready for visualization.
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
      if(tokens[3] == 'firstexons') {
        return('first exons')
      } else if (tokens[3] == 'firstintrons') {
        return('first introns')
      } else {
        return(tokens[3])
      }
    } else if (tokens[2] == 'enhancers') {
      return('enhancers')
    } else if (tokens[2] == 'custom') {
      return(tokens[3])
    }
  })

  flip_tidy = names(tidy)
  names(flip_tidy) = tidy

  return(as.list(flip_tidy))
}

#' Function to check for valid annotations
#'
#' Gives errors if any annotations are not in supported_annotations() (and they are not in the required custom format), basicgenes and detailedgenes are used, or the genome prefixes are not the same for all annotations.
#'
#' @param annotations A character vector of annotations possibly using the shortcuts
check_annotations = function(annotations) {
  # Pull out any custom annotations before checking
  custom_annotations = grep('custom', annotations, value=T)
  annotations = setdiff(annotations, custom_annotations)

  # Check that the annotations are supported, tell the user which are unsupported
  if( !all(annotations %in% supported_annotations()) ) {
    unsupported = base::setdiff(annotations, supported_annotations())

    stop(sprintf('Error: "%s" is(are) not supported. See supported_annotations().',
      paste(unsupported, collapse=', ')))
  }

  # Do not allow basicgenes and detailedgenes at the same time
  if( any(grepl('basicgenes', annotations)) && any(grepl('detailedgenes', annotations)) ) {
    stop('Error: basicgenes and detailedgenes shortcuts may not be used simultaneously.')
  }

  genomes = sapply(annotations, function(a){
    unlist(strsplit(a, '_'))[1]
  }, USE.NAMES = F)

  # Check for same genome on all annotations
  if( length(unique(genomes)) != 1 ){
    stop('Error: genome prefix on all annotations must be the same.')
  }
}

#' Function to expand annotation shortcuts
#'
#' @param annotations A character vector of annotations, possibly using the shortcut accessors
#'
#' @return A vector of data accession-ized names that are ordered from upstream to downstream in the case of knownGenes and islands to interCGI in the case of cpgs.
expand_annotations = function(annotations) {
  are_basicgenes = any(grepl('basicgenes', annotations))
  are_detailedgenes = any(grepl('detailedgenes', annotations))
  are_cpgs = any(grepl('cpgs', annotations))
  which_are_shortcuts = c(which(grepl('basicgenes', annotations)), which(grepl('detailedgenes', annotations)), which(grepl('cpgs', annotations)))

  # expand_shortcuts() will always be run after check_annotations() so we can be
  # sure that the genome prefixes are the same for all annotaitons.
  genome = unique( sapply(annotations, function(a){ unlist(strsplit(a, '_'))[1] }, USE.NAMES=F) )

  if(are_basicgenes || are_detailedgenes || are_cpgs) {

    # Check for shortcut annotation accessors 'cpgs', 'basicgenes', or 'detailedgenes'
    # and create the right annotations based on the genome
    new_annotations = c()
    remove_shortcuts = c()
    if(are_cpgs) {
      new_annotations = paste(genome, 'cpg', c('islands','shores','shelves','inter'), sep='_')
    }
    if(are_basicgenes) {
      new_annotations = c(new_annotations, paste(genome, 'knownGenes', c('1to5kb','promoters','5UTRs','exons','introns','3UTRs'), sep='_'))
    }
    if(are_detailedgenes) {
      new_annotations = c(new_annotations, paste(genome, 'knownGenes',
        c('1to5kb','promoters','exons5UTRs','introns5UTRs','exonsCDSs','intronsCDSs','exons3UTRs','introns3UTRs'), sep='_'))
    }
    annotations = setdiff(c(annotations, new_annotations), annotations[which_are_shortcuts])
  }

  return(annotations)
}

#' Function to subset a tbl_df or grouped_df by a column
#'
#' @param tbl A \code{tbl_df} or \code{grouped_df}.
#' @param col A string indicating which column of of \code{tbl} to subset and order
#' @param col_order A character vector indicating the order of \code{col}.
#'
#' @return A modified version of \code{summary} with \code{col} subsetted by \code{col_order}.
subset_order_tbl = function(tbl, col, col_order) {
  if(!is.null(col)) {
    # Collect all types in the column
    all_col_names = unique(tbl[[col]])

    # Inherit col_order from the order in tbl
    if(is.null(col_order)) {
      col_order = all_col_names
    }

    # Check set equality of col in the summary and the col_order
    if( !dplyr::setequal(all_col_names, col_order) ) {
      if( all(col_order %in% all_col_names) ) {
        tbl = subset(tbl, tbl[[col]] %in% col_order)
      } else {
        stop('There are elements in col_order that are not present in the corresponding column. Check for typos, and check that the elements in col_order do not have 0 tallies in the summarization.')
      }
    }

    # Convert fill to factor with levels in the correct order
    tbl[[col]] = factor(tbl[[col]], levels = col_order)
    # Also convert the levels to tidy names if fill is annotations
    if(col == 'annot_type') {
      levels(tbl[[col]]) = tidy_annotations(col_order)
    }
  }
  return(tbl)
}
