#' findInsertions
#'
#' Identify highly likely insertion-sites in an tagMeppr-object.
#'
#' @author Robin H. van der Weide, \email{r.vd.weide@nki.nl}
#' @param exp The tagMeppr-object of a sample: first run \code{\link{align}}.
#' @param ref A reference object of \code{\link{makeIndex}} or \code{\link{loadIndex}}.
#' @param minP The minimal adjusted p-value.
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
#' findInsertions(exp = C9, ref = reference_hg19_PB, minP = 0.05)
#' }
#' @return The object will be updated with a results-table:
#' \describe{
#' \item{cnt_FWD}{Forward read-counts at this insertion-site}
#' \item{orientation_FWD}{Forward read-orientation at the insertion}
#' \item{D_FWD}{The imbalance of forward reads on both sides of the
#' site: D = 0 means no imbalance}
#' \item{pBinom_FWD}{The probability of the imbalance}
#' \item{cnt_REV}{Reverse read-counts at this insertion-site}
#' \item{orientation_REV}{Reverse read-orientation at the insertion}
#' \item{D_REV}{The imbalance of forward reads on both sides of the site: D = 0
#' means no imbalance}
#' \item{pBinom_REV}{The probability of the imbalance}
#' \item{pMarkerMeppR}{The combined probability of the imbalances}
#' \item{padj}{Bonferroni-adjusted probabilties}
#' \item{orientation}{The orientation of the insertion}
#' }
#'
#' @importFrom dplyr bind_rows
#' @importFrom GenomicRanges countOverlaps makeGRangesFromDataFrame findOverlaps GRanges
#' @importFrom IRanges subsetByOverlaps
#' @importFrom metap sump
#' @importFrom S4Vectors queryHits subjectHits mcols
#' @importFrom stats p.adjust
#' @export
findInsertions = function(exp, ref, minP = 0.05){

  BAMlist = list('FWD' = exp$alignedReadsFWD,
                 'REV' = exp$alignedReadsREV)

  PBbed = ref$PBsites

  readsGRList = list()
  cntTableList = list()
  PBintersect = lapply(BAMlist, function(x){


    PBwithHits = IRanges::subsetByOverlaps(PBbed, x)
    PBwithHits = PBwithHits[seqnames(PBwithHits) != exp$insertName]

    # make GRs
    GR5 = x[x$orientation == 5]
    GR3 = x[x$orientation == 3]

    # intersect with piggyback sites
    PBbed5 = PBwithHits;PBbed5$cnt = 0
    PBbed5$cnt = suppressWarnings(GenomicRanges::countOverlaps(PBbed5,GR5))
    PBbed5$orientation = 5
    # isect 3
    PBbed3 = PBwithHits; PBbed3$cnt = 0
    PBbed3$cnt = suppressWarnings(GenomicRanges::countOverlaps(PBbed3,GR3))
    PBbed3$orientation = 3

    cntTable = dplyr::bind_rows(as.data.frame(PBbed3[PBbed3$cnt > 1]),
                                as.data.frame(PBbed5[PBbed5$cnt > 1]))

    # # ! for each TTAA, get reads
    cntTable$D = 0
    cntTable$pBinom = 1

    readsGR = x
    S4Vectors::mcols(readsGR)[1:4] = NULL

    cntTableGR = GenomicRanges::makeGRangesFromDataFrame(cntTable)

    hitDF = as.data.frame(GenomicRanges::findOverlaps(readsGR, cntTableGR))
    hitDF = split(hitDF$queryHits, hitDF$subjectHits)
    hitDF = lapply(hitDF, function(x){as.data.frame(readsGR[x])})

    cntTablelapply = lapply(1:nrow(cntTable), function(TTAAidx){

      thisTTAA = cntTableGR[TTAAidx]

      # hits = subsetByOverlaps(readsGR, thisTTAA) # <- put outside
      HHH = hitDF[[TTAAidx]]

      # ! ? add padding ? !
      PBpadding = 1
      HHH$placement = NA
      # down is: reads start at the ttaa-start
      HHH$placement[HHH$start >= (thisTTAA@ranges@start-PBpadding)] = 'down'
      # up is: reads end at the ttaa-end
      HHH$placement[HHH$end <= (thisTTAA@ranges@start+4+PBpadding) ] = 'up'

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
      cntTable[TTAAidx,'D'] = Dscore
      cntTable[TTAAidx,'pBinom'] = pBinom
      cntTable[TTAAidx,]
    })

    cntTable = dplyr::bind_rows(cntTablelapply)


    cntTable[is.nan(cntTable$D),"D"] = 0

    list("readsGR" = readsGR, "cntTable" = cntTable[,-c(4,5)])
  })


  FWD = GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = T,
                                                PBintersect$FWD$cntTable)
  REV = GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = T,
                                                PBintersect$REV$cntTable)
  FO = findOverlaps(FWD, REV)

  FWD = as.data.frame(FWD)[S4Vectors::queryHits(FO),]
  REV = as.data.frame(REV)[S4Vectors::subjectHits(FO),]

  merged = cbind(FWD[,-c(4:7)],REV[,-c(1:7)])
  colnames(merged)[4:7]  = paste0(colnames(merged)[4:7],"_FWD")
  colnames(merged)[8:11] = paste0(colnames(merged)[8:11],"_REV")

  # filter out TTAAs with the same D in FWD and REV
  merged = merged[!(merged$D_FWD > 0 & merged$D_REV > 0),]
  merged = merged[!(merged$D_FWD < 0 & merged$D_REV < 0),]

  # filter out TTAAs with a D of 0
  merged = merged[merged$D_FWD != 0 ,]
  merged = merged[merged$D_REV != 0 ,]

  merged = merged[merged$orientation_FWD != merged$orientation_REV,]
  df = merged

  if(nrow(merged) > 0){

    df[df$pBinom_FWD == 0, "pBinom_FWD"] = min(df[df$pBinom_FWD > 0, "pBinom_FWD"])
    df[df$pBinom_REV == 0, "pBinom_REV"] = min(df[df$pBinom_REV > 0, "pBinom_REV"])

    #  do sumlog
    df$pMarkerMeppR = apply(df[,c(7,11)], 1, function(x){y = metap::sump(x)$p})
    df$padj = stats::p.adjust(df$pMarkerMeppR)


  }

  if(nrow(df[df$padj <= minP, ]) > 0){
    df = GenomicRanges::makeGRangesFromDataFrame(df[df$padj <= minP, ], keep.extra.columns = T)
    df$orientation = paste( df$orientation_FWD, df$orientation_REV,  sep = "->")
  } else {
    df = GenomicRanges::GRanges()

    S4Vectors::mcols(df) = data.frame(cnt_FWD = integer(),
                           orientation_FWD = numeric(),
                           D_FWD = numeric(),
                           pBinom_FWD = numeric(),
                           cnt_REV = integer(),
                           orientation_REV = numeric(),
                           D_REV = numeric(),
                           pBinom_REV = numeric(),
                           pMarkerMeppR = numeric(),
                           padj = numeric(),
                           orientation = character())
  }


  exp$results = df
  tmp = exp
  name <- sapply(match.call(expand.dots=TRUE)[-1], deparse)
  assign(name, tmp, envir = parent.frame())


  invisible(df)

}




