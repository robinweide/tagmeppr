#' tagMepprCol
#'
#' A pallette based on CYRUP (github.com/robinweide/CYRUP)
#'
#' @author Robin H. van der Weide, \email{r.vd.weide@nki.nl}
#' @param n The number of colours needed.
#' @export
tagMepprCol = function(n){
  pal = c("#009bef", "#ff5c49", "#949494", "#fed500")
  return(pal[n])
}


#' runIDgen
#'
#' Generate an unique 10-item ID
#'
#' @author Robin H. van der Weide, \email{r.vd.weide@nki.nl}
#' @param n The number of ID's needed.
#' @export
runIDgen <- function(n = 1) {
  a <- do.call(paste0, replicate(5, sample(LETTERS, n, TRUE), FALSE))
  a <- paste0(a, sprintf("%04d", sample(9999, n, TRUE)), sample(LETTERS, n, TRUE))
  return(a)
}


# split multiple [XS]A fields
splitXSA = function(x, field = 'XA'){

  COI = x[,colnames(x) == field]
  NAMES = x$readName

  # split
  CS = reshape2::colsplit(COI, pattern = ';',LETTERS)
  CS[CS == ""] = NA

  # remove na-rows
  notShit = apply(CS, 1, function(x){sum(is.na(x))}) != 26
  CS = CS[notShit, ]
  NAMES = NAMES[notShit]

  cnts = apply(CS, 1, function(x){sum(!is.na(x))}) == 1

  # add name to singles
  singles = data.frame('XA' = CS[cnts, 1], 'RN' = NAMES[cnts])

  # split multiples
  multiples = CS[!cnts, ]

  if(nrow(multiples) == 1){
    multiples = multiples[, as.numeric(apply(as.data.frame(multiples),2,is.na)) < nrow(multiples)]
  } else  if(nrow(multiples) > 0){
    multiples = multiples[, colSums(apply(as.data.frame(multiples),2,is.na)) < nrow(multiples)]
  } else{
    multiples = c()
  }

  out = c()
  if(length(multiples) > 0){
    Mvect = c()
    for(i in 1:ncol(multiples)){
      Mvect = c(Mvect, stats::setNames(multiples[,i], NAMES[!cnts]))
    }
    multiples = data.frame('XA' = Mvect, 'RN' = names(Mvect))
    multiples = multiples[!is.na(multiples[,1]),]
    out = unique(rbind(singles,multiples))
  } else {
    out = unique(singles)
  }

  return(out)
}


empericalTransposonCentre = function(exp, ref){

  transposonSeq = NULL
  if(ref$ITR == "PiggyBac"){
    transposonSeq = tagMeppr::PiggyBacITRs
  } else if(ref$ITR == "SleepingBeauty"){
    transposonSeq = tagMeppr::SleepingBeautyITRs
  } else if(grepl(ref$ITR, pattern = ".fa")){
    # check if exists
    if(file.exists(ref$ITR)){
      transposonSeq = Biostrings::readDNAStringSet(filepath = ref$ITR, use.names = T)
    } else {
      stop('The file ', ref$ITR, ' does not exist.')
    }
  }

  MID = 0
  # 1. get mid of largest N-padded sequence
  if(Biostrings::letterFrequency(transposonSeq, letters = 'N') > 0){

    Nranges = IRanges::reduce(Biostrings::vmatchPattern("N",transposonSeq)[[1]])
    Nranges = Nranges[base::which.max(S4Vectors::width(Nranges))]

    MID = IRanges::mean(Nranges)

  } else {
  # 2. if no Ns, get coverage-distribution
    SL = S4Vectors::width(transposonSeq)
    names(SL) = base::names(transposonSeq)

    gr.windows <- GenomicRanges::tileGenome(SL, tilewidth=10, cut.last.tile.in.chrom=TRUE)

    COF = IRanges::countOverlaps(gr.windows, exp$FWDBAM)
    COR = IRanges::countOverlaps(gr.windows, exp$REVBAM)

    hits = sort(c(gr.windows[which.max(COF)],
                  gr.windows[which.max(COR)]))
    hits = unname(sort(unlist(as.data.frame( IRanges::ranges(hits)))[1:4])[c(1,4)])

    MID = unname(mean(hits))
  }

  return( MID)

}
