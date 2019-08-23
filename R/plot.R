#' plotSite
#'
#' Plot the read distribution of an insertion.
#'
#' @author Robin H. van der Weide, \email{r.vd.weide@nki.nl}
#' @param exp A list or singular instance of tagMeppr-object(s):
#' first run \code{\link{findInsertions}}.
#' @param site Which insertion should be plotted?
#' @param forceDetail Even with more than 1000reds, still use the high-detail version.
#' @param maxReads Set the maximum number of reads. If forward and/or reverse
#' exceeds this, they will be downsampled to this number.
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
#'
#' plotSite(exp = C9, site = 1)
#'
#' }
#' @return An object of the class ggplot.
#'
#' @importFrom GenomicRanges	findOverlaps
#' @importFrom GenomicAlignments width
#' @importFrom BiocGenerics	end	start
#' @importFrom dplyr	bind_rows
#' @importFrom GenomicRanges	findOverlaps GRanges	makeGRangesFromDataFrame seqnames
#' @importFrom ggplot2	aes_string annotate element_blank element_line element_rect element_text expand_scale geom_hline geom_polygon geom_rect geom_vline ggplot ggtitle labs scale_fill_manual scale_x_continuous scale_y_continuous theme unit
#' @importFrom IRanges IRanges coverage
#' @importFrom S4Vectors	queryHits runValue subjectHits
#' @importFrom stats	setNames approx
#' @importFrom scales extended_breaks
#' @export
#'
plotSite = function(exp, site = 1, forceDetail = F, maxReads = Inf){

  if(length(exp$results) == 0){
    stop('There are no insertions found!')
  } else if(site > length(exp$results)){
    stop('There are only ',  length(exp$results), ' insertions found.')
  }
  flank = max(GenomicAlignments::width(exp$FWDBAM))

  readsGRList = list(exp$alignedReadsFWD,
                     exp$alignedReadsREV)

  FOfwd = GenomicRanges::findOverlaps(readsGRList[[1]],  exp$results)
  FOrev = GenomicRanges::findOverlaps(readsGRList[[2]],  exp$results)


  # get TTAA
  thisPB = as.data.frame( exp$results[site])
  # get reads overlapping FWD
  FWDreads = as.data.frame(readsGRList[[1]][unique(S4Vectors::queryHits(FOfwd[S4Vectors::subjectHits(FOfwd) == site]))])
  FWDreads$primer = 'forward'
  # get reads overlapping REV
  REVreads = as.data.frame(readsGRList[[2]][unique(S4Vectors::queryHits(FOrev[S4Vectors::subjectHits(FOrev) == site]))])
  REVreads$primer = 'reverse'

  if(nrow(FWDreads) > maxReads){
    FWDreads = FWDreads[sample(1:nrow(FWDreads), size = maxReads),]
  }

  if(nrow(REVreads) > maxReads){
    REVreads = REVreads[sample(1:nrow(REVreads), size = maxReads),]
  }

  reads = rbind(REVreads, FWDreads)
  reads$Y = base::rank(reads$start + (reads$width/2), ties.method = "first")
  reads = reads[reads$Y,]

  gg = NULL
  YT = getYticks(reads, thisPB)
  if(nrow(reads) < 1000 | forceDetail){

      gg = ggplot2::ggplot(reads,
                           mapping = ggplot2::aes_string(xmin = 'start',
                                                         xmax = "end",
                                                         ymin = "Y-0.5",
                                                         ymax = "Y+0.5",
                                                         fill = 'primer'))+
        ggplot2::geom_rect()

  } else {

    ROI = GenomicRanges::GRanges(paste0(thisPB$seqnames, ":",
                                        thisPB$start - flank,"-",
                                        thisPB$end + flank))

    tmp = lapply(list('reverse' = REVreads, 'forward' = FWDreads), function(READS){
      GR = GenomicRanges::makeGRangesFromDataFrame(READS)

      RANGE = range(GR)
      RANGE = c(min(start(RANGE)), max(end(RANGE)))
      RANGE = GenomicRanges::GRanges(unique(as.character(seqnames(GR))),
                      ranges = IRanges::IRanges(RANGE[1], RANGE[2])  )

      C = IRanges::coverage(GR)

      x <- C[RANGE][[1]]
      x <- data.frame(x = c(start(x), end(x)), y = c(S4Vectors::runValue(x),
                                                     S4Vectors::runValue(x)))
      x <- x[order(x$x),]
      x$pos = (BiocGenerics::start(RANGE):BiocGenerics::end(RANGE))[x$x]
      x$x = NULL
      x = unique(x)
      x
    })

    tmp = dplyr::bind_rows(tmp, .id = 'primer')
    tmp = rbind(tmp,
                data.frame(primer = 'forward',y = 1, pos = max(tmp[tmp$primer == 'forward',"pos"])),
                data.frame(primer = 'reverse',y = 1, pos = min(tmp[tmp$primer == 'reverse',"pos"])) )

    # the primer with the most downstream pos goes on top
    ontop = names(which.max(sapply(split(tmp, tmp$primer), function(x){max(x$pos)})))
    notontop = unique(tmp$primer)[!unique(tmp$primer) %in% ontop]
    tmp$primer = factor(tmp$primer, levels = c(ontop,notontop))
    tmp[tmp$primer == ontop, 2] = abs((tmp[tmp$primer == ontop, 2] - max(tmp[tmp$primer == ontop, 2]) ) -1 )

    # add a dummy for space
    pad = max(tmp[tmp$primer == notontop,"y"]) *0.05

    # take to one ontop and add y to it
    tmp[tmp$primer == ontop, "y"] = max(tmp[tmp$primer == notontop, "y"]) + tmp[tmp$primer == ontop, "y"]

    gg = ggplot2::ggplot(tmp, mapping = ggplot2::aes_string(x = "pos",
                                                     y = "y",
                                                     fill = "primer",
                                                     label = "LBL")) +
      ggplot2::geom_polygon()

    }

  gg = gg +
    ggplot2::geom_hline(yintercept = YT[[3]])+
    ggplot2::scale_fill_manual(na.value ="#949494",
                      values = stats::setNames(c("#009bef","#ff5c49"),
                                        c('forward','reverse'))) +
    ggplot2::scale_y_continuous(breaks = YT[[2]], labels = YT[[1]],
                       expand = ggplot2::expand_scale(add = c(5,5)))+
    ggplot2::scale_x_continuous(breaks = unique(floor(tmp$pos / 1e2))*1e2,
                       expand = c(0,0),
                       limits = c(thisPB$start - flank, thisPB$end + flank)) +
    ggplot2::geom_vline(xintercept = thisPB$start)+
    ggplot2::geom_vline(xintercept = thisPB$end)+
    ggplot2::labs(x = as.character(thisPB$seqnames),
         y = 'count')+
    ggplot2::theme(legend.title = ggplot2::element_blank()) +
    ggplot2::annotate('text', x = (thisPB$end + (flank * 0.975 )),
                      y = 0, hjust = 1, vjust = 0, cex = 3,
             label =   paste0( "Dfwd = " ,round(thisPB$D_FWD, digits = 2),
                               "\nDrev = " ,round(thisPB$D_REV, digits = 2),
                               "\np = ",format.pval(thisPB$padj, digits = 3),
                               '\nIntegration = ',thisPB$casetteStrand )) +
    ggplot2::ggtitle(label = exp$name,
                    subtitle = paste(paste(thisPB$seqnames,
                                           thisPB$start,
                                           sep = ':'), thisPB$end, sep = "-")) +
    ggplot2::theme(panel.spacing.y = ggplot2::unit(0, units = 'cm'),
                   panel.background = ggplot2::element_rect(fill = NA,
                                                            color = "black",
                                                            size = 0.5),
                   legend.key = ggplot2::element_rect(fill = NA),
                   panel.grid = ggplot2::element_blank(),
                   axis.line = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_line(colour = "black",
                                                      size = 0.5,
                                                      lineend = "square"),
                   text = ggplot2::element_text(color = "black"),
                   aspect.ratio = 1/2.2,
                   axis.text = ggplot2::element_text(color = "black"),
                   strip.background = ggplot2::element_blank(),
                   strip.text = ggplot2::element_text(color = "black"),
                   panel.border =  ggplot2::element_rect(fill = NA, color = 'black'))
  print(gg)
}



getYticks = function(reads, thisPB){
  # find out wich is bottom
  bottom = names(which.min(sapply(split(reads, reads$primer), function(x){min(x$start)})))

  # count reads
  FREQ = as.data.frame(table(reads$primer))

  # get the maxima nd mid labels
  BMT_labels = c(FREQ[FREQ$Var1 == bottom,2], 0, FREQ[FREQ$Var1 != bottom,2])
  BMT_breaks = c(0, BMT_labels[1], nrow(reads))

  # find more labels
  addedLabels = scales::extended_breaks(4)(0:max(BMT_labels))[-1]

  # approximate breaks of added labels
  botExta = approx(BMT_breaks[1:2], y = BMT_labels[1:2], xout = addedLabels)
  botExta$x = c(BMT_breaks[1],botExta$x)
  botExta$y = c(BMT_labels[1],botExta$y)
  botExta$x = botExta$x[!is.na(botExta$y)]
  botExta$y = botExta$y[!is.na(botExta$y)]


  topExta = stats::approx(BMT_labels[2:3], y = BMT_breaks[2:3], xout = addedLabels)
  topExta$y = c(BMT_labels[1], topExta$y,  BMT_breaks[3])
  topExta$x = c(0,             topExta$x,  BMT_labels[3])
  topExta$x = topExta$x[!is.na(topExta$y)]
  topExta$y = topExta$y[!is.na(topExta$y)]

  # check for too close to ends: has to be 25% away
  if(botExta$y[1] - floor(botExta$y[1] *0.25) < botExta$y[2]){
    botExta$y = botExta$y[-2]
    botExta$x = botExta$x[-2]
  }

  if(topExta$x[length(topExta$x)] - floor(topExta$x[length(topExta$x)]*0.25) < topExta$x[length(topExta$x) - 1]){
    topExta$y = topExta$y[-(length(topExta$x)-1)]
    topExta$x = topExta$x[-(length(topExta$x)-1)]
  }

  out  = data.frame(breaks = c(botExta$x, topExta$y),
                    labels = c(botExta$y, topExta$x))

  out = unique(out[complete.cases(out),])

  YT = list(out$labels, out$breaks, BMT_breaks[[2]])
  return(YT)
}








