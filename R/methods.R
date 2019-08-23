

#' @export
print.tagMepprIndex <- function(x){

  cat(paste0('tagMepprIndex\n\n'))
  cat(paste0('\tProtocol: ', x$ITR,"\n\n"))
  cat(paste0('\tDirectory: ', dirname(x$index),"\n"))
  cat(paste0('\tFasta: ', basename(x$index),"\n\n"))
  cat(paste0('\tTarget insertion site: ', x$targetInsertionSite,"\n"))
  cat(paste0('\tTIS: ', gsub(x = basename(x$index), pattern = ".fa.gz",
                             replacement = ".tis", fixed = T),"\n"))
  cat(paste0('\tNumber of target insertion sites: ',length(x$TIS), "\n\n"))
}

#' @export
print.tagMepprSample <- function(x){

  # name
  # protocol
  cat(paste0('tagMepprSample\n\n'))
  cat(paste0('\tName: ', x$name,"\n"))
  cat(paste0('\tProtocol: ', x$protocol,"\n"))

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


#' results.tagMepprSample <- function(x){
#'
#'   if(is.null(x$results)){
#'     cat(paste0("\tAnalysed: FALSE\n"))
#'   } else {
#'     print(GenomicRanges::as.data.frame(x$results))
#'   }
#'
#'
#' }
#'

#' results <- function(x){
#'   UseMethod("results", object = "tagMepprSample")
#' }




