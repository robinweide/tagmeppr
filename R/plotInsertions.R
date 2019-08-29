#' plotInsertions
#'
#' Identify highly likely insertion-sites in an tagMeppr-object.
#'
#' @author Robin H. van der Weide, \email{r.vd.weide@nki.nl}
#' @param exp A list or singular instance of tagMeppr-object(s):
#' first run \code{\link{findInsertions}}.
#' @param chrom Which chromosomes should be plotted (even when no insertions are found).
#' @param colours Choose colours for individual tagMeppr-objects (if sideBySide = T).
#' @param sideBySide If true, will show per-object ideogram. If false, will show
#' @param showOrientation Split hits on orientation.
#' @param p Set a p-value treshold
#' all (colour-coded) objects inside one ideogram.
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
#' plotInsertions(exp = C9, chrom = c('chr1','chr2'))
#'
#' }
#' @return An object of the class ggplot.
#'
#' @importFrom GenomicRanges GRanges GRangesList
#' @importFrom stats setNames
#' @importFrom grDevices colorRampPalette
#' @importFrom ggrepel geom_text_repel
#' @importFrom scales hue_pal
#' @importFrom ggplot2 ggplot aes_string coord_cartesian element_blank element_text facet_grid geom_segment labs scale_color_manual scale_x_continuous theme theme_minimal xlab ylab
#' @export
plotInsertions = function(exp, chrom = paste0('chr', c(1:25,'X')), colours = NULL,
                          sideBySide = F, p = 0.05, showOrientation = F){

  toPlot = NULL
  SI = NULL
  # check if exp is a list of exps or a simpe one
  if( 'results' %in%  names(exp) ){
    #################################################################### one exp
    toPlot = exp$results
    SI    = exp$seqinfo

    if(length(exp) == 0){
      toPlot = GenomicRanges::GRanges()
    } else {
      toPlot$name = exp$name
    }

    if(is.null(colours)){
      colours = 'black'
    }

  } else {
    ############################################################### multiple exp
    SI = exp[[1]]$seqinfo

    # check that all have done findInsertions!
    YN = sapply(exp, function(x){'results' %in%  names(x)})
    if( unique(YN) != T){
      stop("Item(s) ", paste0(which(!YN), collapse = ', '),
           " haven't been ran through findInsertions().")
    }

    # add samplename
    toPlot = lapply(exp, function(x){
     tmp = x$results
     if(length(tmp) == 0){
       tmp = NULL
     } else {
       tmp$name = x$name
     }
     tmp
     })

    # combine into one df
    toPlot = unlist(GenomicRanges::GRangesList(toPlot[!sapply(toPlot, is.null)]))
  }

  ############################################################ make the ideogram
  gg = plotIdeo(SI, chrom, orientationFacet = T)

  ################################################## filter on chroms and pvalue
  toPlot = toPlot[toPlot$padj < p,]
  if(!is.null(chrom)){

    toPlotChroms = as.character(GenomicRanges::seqnames(toPlot))
    toPlot = toPlot[toPlotChroms %in% chrom, ]

  }

  ###################################################################### make df
  names(toPlot) = NULL
  toPlot = as.data.frame(toPlot)[,c('seqnames','start','end','name','strand')]
  toPlot = toPlot[stats::complete.cases(toPlot),]
  toPlot = stats::setNames(toPlot,
                           nm = c('chrom',"chromStart", 'chromEnd','C','S'))
  toPlot$chrom = factor(toPlot$chrom, levels = chrom)

  ################################################ no insertions, just plot ideo
  if(length(toPlot) == 0){
    gg = gg + ggplot2::facet_grid( chrom ~ ., switch = "y" )
    gg = addGGBeauty(gg)

    return(gg)
  }

  #################################################### there are lines in toPlot
  if(sideBySide){
    ################################################################## SBS-style
    gg = makeSbS(gg, toPlot = toPlot)
  } else {

    pmtp = prepareMultiToPlot( toPlot, colours )
    toPlot = pmtp[[1]];      colours = pmtp[[2]]

    if(showOrientation){
      ################################################################ ORI-style
      gg = gg + ggplot2::facet_grid( chrom ~ S, switch = "y" )
    } else {
      ############################################################## basic-style
      gg = gg + ggplot2::facet_grid( chrom ~ ., switch = "y" )
    }

    gg = gg +
      ggplot2::geom_segment( data = toPlot,
                           mapping = ggplot2::aes_string(x = 'chromStart',
                                                         xend = 'chromEnd',
                                                         col = "factor(label)"),
                           y = -5,
                           yend = 5) +
    ggplot2::scale_color_manual(values = colours) +
      ggplot2::labs(col = '')

  }

  gg = addGGBeauty(gg)

  return(gg)
}


plotIdeo <- function(SI, chroms, orientationFacet = T){

  ########################################################### mangle into BED-df
  SIdf = data.frame(chrom = SI@seqnames,
                    chromStart = 1,
                    chromEnd =SI@seqlengths)
  SIdf = SIdf[SIdf$chrom %in% chroms, ]
  SIdf$chrom = factor(SIdf$chrom, levels = chroms)

  ################################################################ add stranding
  if(orientationFacet){
    SIdf = rbind(SIdf, SIdf)
    SIdf$S = rep(c("+","-"), each = nrow(SIdf)/2)
  }

  ################################################################### set ggplot
  thisGG =  ggplot2::ggplot() +
            ggplot2::geom_segment( data = SIdf,
                                   mapping = ggplot2::aes_string(x = 'chromStart',
                                                                 xend = 'chromEnd'),
                                   y = 0,
                                   yend = 0,
                                   lineend = "round",
                                   color = "#D8D8D8",
                                   size = 5 )

  return(thisGG)
}




addGGBeauty = function(gg){
  thisGG =  gg + ggplot2::theme_minimal() +
    ggplot2::scale_x_continuous( expand = ggplot2::expand_scale(mult = c(0,0.1)),
                                 breaks = seq( 0, 2.5e8, 50e6 ),
                                 labels = c( 1, seq( 50, 250, 50  ) ) ) +
    ggplot2::coord_cartesian(ylim = c(-5,5),)+
    ggplot2::ylab( "" ) +
    ggplot2::xlab( "Genomic positions (Mb)" ) +
    ggplot2::theme(
      legend.position = "top",
      strip.text.y = ggplot2::element_text(angle = 180,
                                           hjust = 1,
                                           vjust = 0.5),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.spacing.y = ggplot2::unit(.2, "lines")
    )

  return(thisGG)
}



makeSbS = function(toPlot, gg){

  gg = gg +
    ggplot2::facet_grid( chrom ~ C, switch = "y" ) +
    ggplot2::geom_segment( data = toPlot,
                           mapping = ggplot2::aes_string(x = 'chromStart',
                                                         xend = 'chromEnd'),
                           y = -5,
                           yend = 5,
                           col = 'black')
  return(gg)
}


prepareMultiToPlot = function(toPlot, colours){

  # get counts for each integration
  toPlot$ID = apply(toPlot[,1:3], 1, paste0, collapse ="/")
  toPlot$ID = gsub(toPlot$ID, pattern = "[ ]*", replacement = '')
  toPlot = merge(toPlot, table(toPlot$ID), by.x = 'ID', by.y = 'Var1')

  # set labels
  toPlot$label = toPlot$C

  #collapse sites with more than one Freq
  toPlot$label[ toPlot$Freq > 1] = 'multiple'
  toPlot = unique(toPlot[,c(2:4,6:8)])
  toPlot$label = as.factor(toPlot$label)

  # set colours
  labs = levels(toPlot$label)[!grepl(pattern = 'multiple',
                                     x = levels(toPlot$label))]
  if(is.null(colours)){

    colours = scales::hue_pal(l = 60)(length(labs))

  }
  colours = grDevices::colorRampPalette(colours)(length(labs))
  names(colours) = labs

  # if multiples, add colour
  if("multiple" %in% levels(toPlot$label)){
    colours = c(colours ,c('multiple' = 'black'))
  }

  # sample-specific insertions should not be annotated
  toPlot$Freq[ toPlot$Freq == 1] = NA


  return(list(toPlot, colours))
}



