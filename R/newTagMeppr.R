#' newTagMeppr
#'
#' Construct a tagMepprExperiment-object from a tagmap-sample.
#'
#' @author Robin H. van der Weide, \email{r.vd.weide@nki.nl}
#' @param F1 The path to forward_R1.fastq.gz.
#' @param F2 The path to forward_R2.fastq.gz.
#' @param R1 The path to reverse_R1.fastq.gz.
#' @param R2 The path to reverse_R2.fastq.gz.
#' @param name The name of the sample.
#' @param protocol Either PiggyBac or SleepingBeauty are supported, but other
#' protocols are allowed.
#'
#' @examples
#' \dontrun{
#' C1 = newTagMeppr(F1 = '/home/A.Dent/analysis42/clone1_FWD_R1.fq.gz',
#'                  F2 = '/home/A.Dent/analysis42/clone1_FWD_R2.fq.gz',
#'                  R1 = '/home/A.Dent/analysis42/clone1_REV_R1.fq.gz',
#'                  R2 = '/home/A.Dent/analysis42/clone1_REC_R2.fq.gz',
#'                  name = "clone1",
#'                  protocol = 'PiggyBac')
#'
#' }
#' @return A tagMepprSample-object with the following slots:
#'
#' \describe{
#' \item{F1}{The path of forward_R1.fastq.gz.}
#' \item{F2}{The path of forward_R2.fastq.gz.}
#' \item{R1}{The path of reverse_R1.fastq.gz.}
#' \item{R2}{The path of reverse_R2.fastq.gz.}
#' \item{name}{The name of the sample.}
#' \item{protocol}{The name of the protocol.}
#'}
#'
#' @export
#'
newTagMeppr <- function(F1, F2, R1, R2, name, protocol = 'PiggyBac'){

###################################################################### check FQs

  if(!file.exists(F1)){
    stop('F1 does not exist!')
  }

  if(!file.exists(F2)){
    stop('F2 does not exist!')
  }

  if(!file.exists(R1)){
    stop('R1 does not exist!')
  }

  if(!file.exists(R2)){
    stop('R2 does not exist!')
  }

################################################################# check protocol
  if(tolower(protocol) %in% c('piggybac', 'pb')){
    # pb
    protocol = 'PiggyBac'
  } else if(tolower(protocol) %in% c('sleepingbeauty', 'sb')){
    # sb
    protocol = 'SleepingBeauty'
  } else {
    # custom
    warning('Custom protocol chosen.
            Only PiggyBac and SleepingBeauty are fully supported.
            Thread carefully!')
  }

  TMobj = structure(list(F1 ,
                         F2 ,
                         R1 ,
                         R2 ,
                         name = name,
                         protocol = protocol),
                    class = c("tagMepprSample", "list"))

  return(TMobj)
}
