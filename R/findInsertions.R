#' findInsertions
#'
#' Identify highly likely insertion-sites in an tagMeppr-object.
#'
#' @author Robin H. van der Weide, \email{r.vd.weide@nki.nl}
#' @param exp The tagMeppr-object of a sample: first run \code{\link{align}}.
#' @param ref A reference object of \code{\link{makeIndex}} or \code{\link{loadIndex}}.
#' @param padding Add padding to the start/end of reads in relation to the TIS.
#' Default is half of the TIS (i.e. 2bp whit Piggybac), meaning that the mapper
#' can report a match Xbp further than the TIS (i.e. TTAA).
#' @examples
#' \dontrun{
#' reference_hg19_PB = makeIndex(indexPath = '/home/A.Dent/analysis42/',
#'                               bsgenome = 'BSgenome.Hsapiens.UCSC.hg19',
#'                               ITR = 'PiggyBac')
#'
#' C9 = newTagMeppr(F1 = 'clone9_FWD_R1.fq.gz',
#'                  F2 = 'clone9_FWD_R2.fq.gz',
#'                  R1 = 'clone9_REV_R1.fq.gz',
#'                  R2 = 'clone9_REV_R2.fq.gz',
#'                  name = "clone9",
#'                  protocol = 'PiggyBac')
#'
#' align(exp = C9, ref = reference_hg19_PB, cores = 30, empericalCentre = T)
#'
#' findInsertions(exp = C9, ref = reference_hg19_PB)
#' }
#' @return The object will be updated with a results-table:
#' \describe{
#' \item{fwdCount}{Forward read-counts at this insertion-site}
#' \item{fwdOrient}{Forward read-orientation at the insertion}
#' \item{D_FWD}{The imbalance of forward reads on both sides of the
#' site: D = 0 means no imbalance}
#' \item{fwdBinom}{The probability of the imbalance}
#' \item{revCount}{Reverse read-counts at this insertion-site}
#' \item{revOrient}{Reverse read-orientation at the insertion}
#' \item{D_REV}{The imbalance of forward reads on both sides of the site: D = 0
#' means no imbalance}
#' \item{revBinom}{The probability of the imbalance}
#' \item{pMarkerMeppR}{The combined probability of the imbalances}
#' \item{padj}{Bonferroni-adjusted probabilties}
#' \item{orientation}{The orientation of the insertion}
#' }
#'
#' @importFrom dplyr bind_rows
#' @importFrom GenomicRanges countOverlaps makeGRangesFromDataFrame findOverlaps GRanges strand
#' @importFrom IRanges subsetByOverlaps
#' @importFrom metap sump
#' @importFrom S4Vectors queryHits subjectHits mcols
#' @importFrom stats p.adjust
#' @export
findInsertions = function(exp, ref, padding = NULL){

  if(is.null(padding)){
    padding = IRanges::width(ref$TIS[1])/2
  }; if(is.na(padding)){padding = 2}

  BAMlist = list('FWD' = exp$alignedReadsFWD,
                 'REV' = exp$alignedReadsREV)

  if(is.null(BAMlist$FWD)){
    stop('No reads found. Did you run align()?')
  }

  TIS = ref$TIS

  readsGRList = list()
  cntTableList = list()
  TISintersect = lapply(BAMlist, function(x){


    TISwithHits = suppressWarnings(IRanges::subsetByOverlaps(TIS, x))
    TISwithHits = TISwithHits[seqnames(TISwithHits) != exp$insertName]

    # make GRs
    GRT = x[x$beforePad == T]
    GRF = x[x$beforePad == F]

    # intersect with piggyback sites
    TIST = TISwithHits;TIST$cnt = 0
    TIST$cnt = suppressWarnings(GenomicRanges::countOverlaps(TIST,GRT))
    TIST$beforePad = T
    # isect 3
    TISF = TISwithHits; TISF$cnt = 0
    TISF$cnt = suppressWarnings(GenomicRanges::countOverlaps(TISF,GRF))
    TISF$beforePad = F

    cntTable = dplyr::bind_rows(as.data.frame(TISF[TISF$cnt > 1]),
                                as.data.frame(TIST[TIST$cnt > 1]))

    # # ! for each TIS, get reads
    cntTable$D = 0
    cntTable$pBinom = 1

    readsGR = x

    cntTableGR = GenomicRanges::makeGRangesFromDataFrame(cntTable)

    hitDF = as.data.frame(suppressWarnings(GenomicRanges::findOverlaps(readsGR, cntTableGR)))
    hitDF = split(hitDF$queryHits, hitDF$subjectHits)
    hitDF = lapply(hitDF, function(x){as.data.frame(readsGR[x])})

    cntTablelapply = lapply(1:nrow(cntTable), function(TISidx){

      thisTIS = cntTableGR[TISidx]

      # hits = subsetByOverlaps(readsGR, thisTIS) # <- put outside
      HHH = hitDF[[TISidx]]

      HHH$placement = NA

      # down is: reads start at the TIS-start
      HHH$placement[HHH$start >= (thisTIS@ranges@start-padding)] = 'down'
      # up is: reads end at the TIS-end
      HHH$placement[HHH$end <= (thisTIS@ranges@start+4+padding) ] = 'up'

      # add Y
      HHH$ID = -base::rank(rowMeans(HHH[,2:3]), ties.method = "first")

      # do a simple binomtest
      stats = c('up' = sum(HHH$placement == 'up', na.rm = T),
                'down' = sum(HHH$placement == 'down', na.rm = T),
                'none' = sum(is.na(HHH$placement)))

      pBinom = 1
      Dscore = unname((stats[1]/sum(stats[1:2]))-.5)
      Dscore = Dscore*2
      if(sum(stats[1:2]) > 0){
        # hits!
        pBinom = stats::binom.test(stats[1], sum(stats[1:2]))$p.value
      }


      # fill cntTable
      cntTable[TISidx,'D'] = Dscore
      cntTable[TISidx,'pBinom'] = pBinom
      cntTable[TISidx,]
    })

    cntTable = dplyr::bind_rows(cntTablelapply)


    cntTable[is.nan(cntTable$D),"D"] = 0

    list("readsGR" = readsGR, "cntTable" = cntTable[,-c(4,5)])
  })


  FWD = GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = T,
                                                TISintersect$FWD$cntTable)
  REV = GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = T,
                                                TISintersect$REV$cntTable)
  FO = GenomicRanges::findOverlaps(FWD, REV)

  FWD = as.data.frame(FWD)[S4Vectors::queryHits(FO),]
  REV = as.data.frame(REV)[S4Vectors::subjectHits(FO),]

  merged = cbind(FWD[,-c(4:7)],REV[,-c(1:7)])
  colnames(merged)[4:7]  = paste0('fwd',colnames(merged)[4:7])
  colnames(merged)[8:11] = paste0('rev',colnames(merged)[8:11])

  # filter out TISs with the same D in FWD and REV
  merged = merged[!(merged$fwdD > 0 & merged$revD > 0),]
  merged = merged[!(merged$fwdD < 0 & merged$revD < 0),]

  # filter out TISs with a D of 0
  merged = merged[merged$fwdD != 0 ,]
  merged = merged[merged$revD != 0 ,]

  merged = merged[merged$fwdbeforePad != merged$revbeforePad,]
  df = merged
  rownames(df) = NULL
  if(nrow(merged) > 0){

    df[df$fwdpBinom == 0, "fwdpBinom"] = min(c(2.2e-16,df[df$fwdpBinom > 0, "fwdpBinom"]))
    df[df$revpBinom == 0, "revpBinom"] = min(c(2.2e-16,df[df$revpBinom > 0, "revpBinom"]))

    df$fwdpBinom[df$fwdpBinom == Inf] = 2.2e-16
    df$revpBinom[df$revpBinom == Inf] = 2.2e-16

    #  do sumlog
    df$pval = apply(df[,c(7,11)], 1, function(x){y = metap::sump(x)$p})
    df$padj = stats::p.adjust(df$pval)


  }

  if(nrow(df) > 0){
    df = GenomicRanges::makeGRangesFromDataFrame(df, keep.extra.columns = T)

    # __R|F__ = +
    ORIENTATION = ifelse(df$revD < 0 & df$fwdD > 0,
                          'RF',
                          'FR')
    # __F|R__ = -
    ORIENTATIONcheck = ifelse(df$fwdD < 0 & df$revD > 0,
                          'FR',
                          'RF')

    # catch both
    ORIENTATION[ORIENTATION != ORIENTATIONcheck] = NA

    # switch
    if(is.na(exp$rev5_fwd3)){
      warning('checkPrimer not used.
              Will not try to flip orientation.')
    } else if(!exp$rev5_fwd3){
      ORIENTATION = ifelse(ORIENTATION == 'RF', 'FR','RF')
      }

    GenomicRanges::strand(df) =  ifelse(ORIENTATION == 'RF', "+", "-")

    colnames(S4Vectors::mcols( df)) = c('fwdCount',
                                        'fwdOrient',
                                        'fwdD',
                                        'fwdBinom',
                                        'revCount',
                                        'revOrient',
                                        'revD',
                                        'revBinom',
                                        'pval',
                                        'padj')
  } else {
    df = GenomicRanges::GRanges()

    S4Vectors::mcols(df) = data.frame(
                           fwdCount = integer(),
                           fwdOrient = numeric(),
                           fwdD = numeric(),
                           fwdBinom = numeric(),
                           revCount = integer(),
                           revOrient = numeric(),
                           revD = numeric(),
                           revBinom = numeric(),
                           pval = numeric(),
                           padj = numeric())
  }

  exp$results = df

  ##################################################################### assigner
  tmp = exp
  # get arguments
  name <- sapply(match.call(expand.dots=TRUE)[-1], deparse)
  #find argument postion for exp
  AP = which(names(name) == 'exp')

  assign(name[AP], tmp, envir = parent.frame())

  invisible(df)

}




