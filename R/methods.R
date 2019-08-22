print.tagMepprIndex <- function(x){

  cat(paste0('tagMepprIndex\n\n'))
  cat(paste0('\tProtocol: ', x$ITR,"\n\n"))
  cat(paste0('\tDirectory: ', dirname(x$index),"\n"))
  cat(paste0('\tfasta: ', basename(x$index),"\n\n"))
  cat(paste0('\tTarget insertion site: ', x$targetInsertionSite,"\n"))
  cat(paste0('\tTIS: ', gsub(x = basename(x$index), pattern = ".fa.gz",
                             replacement = ".tis", fixed = T),"\n"))
  cat(paste0('\tNumber of target insertion sites: ',length(x$TIS), "\n\n"))
}


print.tagMepprSample <- function(x){

  # name
  # protocol
  # if aligned: num(aligned), num(informative)
  # if FI: num(sig sites)

}
