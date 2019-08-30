#' @return \code{NULL}
#'
#' @title print
#' @param x A tagMepprIndex object
#' @param ... further arguments passed to or from other methods.
#' @export
print.tagMepprIndex <- function(x, ...){

  cat(paste0('tagMepprIndex\n\n'))
  cat(paste0('\tProtocol: ', x$ITR,"\n\n"))
  cat(paste0('\tDirectory: ', dirname(x$index),"\n"))
  cat(paste0('\tFasta: ', basename(x$index),"\n\n"))
  cat(paste0('\tTarget insertion site: ', x$targetInsertionSite,"\n"))
  cat(paste0('\tTIS: ', gsub(x = basename(x$index), pattern = ".fa.gz",
                             replacement = ".tis", fixed = T),"\n"))
  cat(paste0('\tNumber of target insertion sites: ',length(x$TIS), "\n\n"))
}

#' @return \code{NULL}
#'
#' @title print
#' @param x A tagMepprSample object
#' @param ... further arguments passed to or from other methods.
#' @importFrom GenomicRanges reduce
#' @export
print.tagMepprSample <- function(x, ...){

  # name
  # protocol
  cat(paste0('tagMepprSample\n\n'))
  cat(paste0('\tName: ', x$name,"\n"))
  cat(paste0('\tProtocol: ', x$protocol,"\n"))

  if(is.na(x$rev5_fwd3)){
    cat(paste0('\tPrimer checked: no (highly recommended!!!)\n'))
  } else if(x$rev5_fwd3){
    cat(paste0('\tPrimer checked: yes\n'))
  } else {
    cat(paste0('\tPrimer checked: yes (flipped)\n'))
  }

  # if aligned: num(aligned), num(informative)
  if(is.null(x$alignedReadsFWD)){
    cat(paste0("\tAligned: FALSE\n"))
  } else {

    cat(paste0('\tAlignment-folder: ', x$alignmentFolder,"\n"))
    cat(paste0('\tInformative FWD-reads: ', length(unique(x$alignedReadsFWD$readName)),"\n"))
    cat(paste0('\tInformative REV-reads: ', length(unique(x$alignedReadsREV$readName)),"\n"))

  }

  # if FI: num(sig sites)
  if(is.null(x$results)){
    cat(paste0("\tAnalysed: FALSE\n"))
  } else {

    cat(paste0('\tUnique TISs covered: ', length(GenomicRanges::reduce(x$results)),"\n"))
    cat(paste0('\t\tp<0.05: ', length(GenomicRanges::reduce(x$results[x$results$padj < 0.05])),"\n"))

  }

  cat(paste0("\n"))
}

#' Obtain the results
#'
#' Running this will generate a data.frame of results
#'
#' @title results
#' @param tagMepprSample Your sample
#' @param alpha A p-value treshold (default: 1)
#' @param countThreshold A p-value treshold (default: 0)
#' @param orientation Select insertions with a specific orientation
#' @param ... other arguments
#' @examples
#' \dontrun{
#' results(tagMepprSample, alpha = 0.05, countThreshold = 100)
#' }
#' @export
results <- function(tagMepprSample, alpha = 1, countThreshold = 0, orientation = "*") UseMethod("results")


#' @export
#' @importFrom GenomicRanges as.data.frame
results.tagMepprSample <- function(tagMepprSample, alpha = 1, countThreshold = 0, orientation = "*"){
  x = NULL
  if(is.null(tagMepprSample$results)){
    cat(paste0("\tFirst run findInsertions() \n"))
  } else {
    x = GenomicRanges::as.data.frame(tagMepprSample$results)

    if(orientation %in% c("-", "+")){
      x = x[x$strand == orientation,]
    }

    x = x[x$padj < alpha &
            x$fwdCount > countThreshold &
            x$revCount > countThreshold,]
    x$width = NULL
    x$fwdD = round(x$fwdD,digits = 2)
    x$revD = round(x$revD,digits = 2)

    x$pval[x$pval < 2e-16] = 2e-16
    x$padj[x$padj < 2e-16] = 2e-16
    rownames(x) = NULL
    x = x[,-c(6,8,10,12,13)][,c(1:5,7,6,8,9)]
  }
  return(x)
}
