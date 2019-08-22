runIDgen <- function(n = 5000) {
  a <- do.call(paste0, replicate(5, sample(LETTERS, n, TRUE), FALSE))
  paste0(a, sprintf("%04d", sample(9999, n, TRUE)), sample(LETTERS, n, TRUE))
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

empericalTransposonCentre = function(exp){

  scaffoldsSeq = Biostrings::readAAStringSet(filepath = exp$transposonPath, use.names = T)
  SL = S4Vectors::width(scaffoldsSeq[names(scaffoldsSeq) == exp$insertName])
  names(SL) = exp$insertName
  gr.windows <- GenomicRanges::tileGenome(SL, tilewidth=10, cut.last.tile.in.chrom=TRUE)

  COF = IRanges::countOverlaps(gr.windows, exp$FWDBAM)
  COR = IRanges::countOverlaps(gr.windows, exp$REVBAM)

  hits = sort(c(gr.windows[which.max(COF)],
                gr.windows[which.max(COR)]))
  hits = unname(sort(unlist(as.data.frame( IRanges::ranges(hits)))[1:4])[c(1,4)])
  return( unname(mean(hits)))

}
