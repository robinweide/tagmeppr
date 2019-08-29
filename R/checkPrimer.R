

#' checkPrimer
#'
#' This tool checks the assumptions about the primers.
#'
#' @author Robin H. van der Weide, \email{r.vd.weide@nki.nl}
#' @param fwdPrimer A character-string of the forward-primer used.
#' @param revPrimer A character-string of the reverse-primer used.
#' @param exp The tagMeppr-object of a sample: first run \code{\link{align}}.
#' @param ITR Can take PiggyBac (default), SleepingBeauty, or a path to a 1000xN-padded ITR.fasta.
#' @details
#'
#'
#' The expected general layout for the ITR-sequence looks like this:
#'
#' \code{|---ITR---NNN...NNN---ITR---|}
#'
#' The primers are expected to be 5'-end for the reverse and 3' for the forward:
#'
#' \code{<rev}
#' \code{|---ITR---NNN...NNN---ITR---|}
#' \code{                         fwd>}
#'
#' This tool checks these assumptions and sets the rev5_fwd3 flag to TRUE.
#'
#' @examples
#' \dontrun{
#'
#' C9 = newTagMeppr(F1 = 'clone9_FWD_R1.fq.gz',
#'                  F2 = 'clone9_FWD_R2.fq.gz',
#'                  R1 = 'clone9_REV_R1.fq.gz',
#'                  R2 = 'clone9_REV_R2.fq.gz',
#'                  name = "clone9",
#'                  protocol = 'PiggyBac')
#'
#' checkPrimer(fwdPrimer = "CGTCAATTTTACGCAGACTATC",
#'             revPrimer = "GTACGTCACAATATGATTATCTTTCTAG",
#'             exp =  C9,
#'             ITR = 'PiggyBac')
#'
#' }
#' @return The experiment-object will be updated with the rev5_fwd3-flag, which
#' will tell all downstream analyses if our assumptions are correct.
#'
#' @importFrom Biostrings reverseComplement DNAString readDNAStringSet letterFrequency vmatchPattern
#' @export
checkPrimer <- function(fwdPrimer, revPrimer, exp, ITR = 'PiggyBac'){
  rev5_fwd3 = F

  if(exp$protocol != ITR){
    stop('Protocol given in exp (', exp$protocol,
         ') is not the same as given as ITR (',ITR,').')
  }

  ############################################################# get revComplement
  fwdPrimerCompl = Biostrings::reverseComplement(Biostrings::DNAString(fwdPrimer))
  revPrimerCompl = Biostrings::reverseComplement(Biostrings::DNAString(revPrimer))

  ##################################################################### load ITR
  transposonSeq = NULL
  if(ITR == "PiggyBac"){
    transposonSeq = tagMeppr::PiggyBacITRs
  } else if(ITR == "SleepingBeauty"){
    transposonSeq = tagMeppr::SleepingBeautyITRs
  } else if(grepl(ITR, pattern = ".fa")){
    # check if exists
    if(file.exists(ITR)){
      transposonSeq = Biostrings::readDNAStringSet(filepath = ITR, use.names = T)
      # check if N-padded
      N1k = Biostrings::letterFrequency(transposonSeq, letters = "N") == 1000
      if(!N1k){
        stop('The file ', ITR, " has no padding of 1000 N's between the arms.")
      }
    } else {
      stop('The file ', ITR, ' does not exist.')
    }
  } else {
    stop('Please set ITR to either "PiggyBac", "SleepingBeauty", or as a path to a .fasta-file!')
  }

  ##################################################################### get arms
  NpadRange = Biostrings::vmatchPattern(transposonSeq,
                                        pattern = paste0(rep('N', 1e3), collapse = ''))

  ######################################################################### find

  hitF = Biostrings::vmatchPattern(transposonSeq,pattern = fwdPrimer)[1]
  hitR = Biostrings::vmatchPattern(transposonSeq,pattern = revPrimer)[1]

  if(length(hitF[[1]]) == 0){
    hitF = Biostrings::vmatchPattern(transposonSeq,pattern = fwdPrimerCompl)[1]
  }

  if(length(hitR[[1]]) == 0){
    hitR = Biostrings::vmatchPattern(transposonSeq,pattern = revPrimerCompl)[1]
  }


  ############################################################### check if found
  if(length(hitF[[1]]) == 0){
    stop('No match between fwdPrimer (including revComplement) and the sequence.')
  }

  if(length(hitR[[1]]) == 0){
    stop('No match between revPrimer (including revComplement) and the sequence.')
  }

  ############################################# check if they are on unique arms
  belowF = unlist(hitF[[1]] < NpadRange)
  belowR = unlist(hitR[[1]] < NpadRange)

  if(belowR == belowF){
    if(belowF){
      stop('Primers are both found on the first ITR!')
    } else {
      stop('Primers are both found on the second ITR!')
    }
  }

  ####################### check if start(reverse primer) < start(forward primer)
  if(hitR[[1]] < hitF[[1]]){
    rev5_fwd3 = T
  } else {
    # reverse is on second ITR, which is not what I expect
    rev5_fwd3 = F
  }

  exp$rev5_fwd3 = rev5_fwd3

  ##################################################################### assigner
  tmp = exp
  # get arguments
  name <- sapply(match.call(expand.dots=TRUE)[-1], deparse)
  #find argument postion for exp
  AP = which(names(name) == 'exp')

  assign(name[AP], tmp, envir = parent.frame())

  invisible(rev5_fwd3)

}
