#' makeIndex
#'
#' Generate a hybrid reference genome
#'
#' @author Robin H. van der Weide, \email{r.vd.weide@nki.nl}
#' @param indexPath A path to the folder to store the reference-files.
#' @param bsgenome A BSgenome onject (default: BSgenome.Hsapiens.UCSC.hg19)
#' @param ITR Can take PiggyBac (default), SleepingBeauty, or a path to a 1000xN-padded ITR.fasta.
#' @param targetInsertionSite The target insertion site of your system ("TTAA" for piggyBac, "TA" for sleepingBeauty)
#' @param blockSizeMult Set a blocksize multiplication factor for BWA index. If
#' left NULL, will use calculate optimal value. Change only when absolutely
#' needed and when you have read
#' \href{https://github.com/lh3/bwa/issues/104}{this thread}.
#' @param verbose Do you want a more chatty function?
#'
#' @examples
#' \dontrun{
#' reference_hg19_PB = makeIndex(indexPath = '/home/A.Dent/analysis42/',
#'                               bsgenome = 'BSgenome.Hsapiens.UCSC.hg19',
#'                               ITR = 'PiggyBac')
#' }
#' @return Resulting files stored as [indexPath]/[bsgenome]_[ITR]_tagMeppRindex.*.
#' The .fa file will be the hybrid reference genome, while the .tis file stores
#' all possible insertions sites.
#'
#' The function also returns a tagMepprIndex-object with the following slots:
#'
#' \describe{
#' \item{index}{The path to the hybrid reference genome}
#' \item{ITR}{The name of the used transposon-system or ITR.fa}
#' \item{TIS}{A GRanges-object of all target insertion sites.}
#' \item{targetInsertionSite}{The targetInsertionSite}
#' \item{insertName}{The name of the transposon in the hybrid reference}
#' \item{seqinfo}{A seqinfo object of all chromosomes in index.}
#'}
#'
#' @importFrom Biostrings readDNAStringSet vmatchPattern writeXStringSet
#' @importFrom dplyr bind_rows
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom rtracklayer export.bed
#' @importFrom S4Vectors width
#' @importFrom BSgenome getSeq
#' @importFrom BSgenome.Hsapiens.UCSC.hg19 BSgenome.Hsapiens.UCSC.hg19
#' @importFrom utils write.table
#' @export
#'
makeIndex = function(indexPath, bsgenome = NULL, ITR = "PiggyBac", targetInsertionSite = 'TTAA', blockSizeMult = NULL, verbose = F){


  ################################################################# load bsgenome
  if(verbose){message('Loading references')}

  scaffoldsSeq = NULL
  if(is.null(bsgenome)){
    bsgenome = "BSgenome.Hsapiens.UCSC.hg19"
  }
  scaffoldsSeq = BSgenome::getSeq(get(bsgenome))
  scaffFai =  GenomeInfoDb::seqinfo(get(bsgenome))

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
    } else {
      stop('The file ', ITR, ' does not exist.')
    }
  } else {
    stop('Please set ITR to either "PiggyBac", "SleepingBeauty", or as a path to a .fasta-file!')
  }

  transFai =  GenomeInfoDb::seqinfo(transposonSeq)
  GenomeInfoDb::isCircular(transFai) <- F
  GenomeInfoDb::genome(transFai) <- base::names(transposonSeq)

  index = c(transposonSeq,scaffoldsSeq)

  ################################################################## write index

  if(verbose){message('Writing reference')}

  PAD = gsub(paste0(indexPath,'/',bsgenome,"_",ITR,"_tagMeppRindex.fa.gz"),
             pattern = "//",replacement = '/')

  Biostrings::writeXStringSet(x = index, filepath = PAD, compress = T)


  ########################################################################## FAI

  FAI = suppressWarnings(merge(transFai, scaffFai))
  FAIdf = as.data.frame(FAI)
  FAIPAD =  gsub(x = PAD, pattern = ".fa.gz",
                 replacement = ".fa.gz.fai", fixed = T)

  utils::write.table(x = FAIdf,  file = FAIPAD,quote = F, sep = "\t")

  ######################################################################## index

  if(verbose){message(paste0('BWA index ',
                      "(this happens in the background and could take hours)"))}
  stdBlockSize = 1e7

  if(is.null(blockSizeMult)){
    blockSizeMult = ceiling(
      (sum(as.numeric( S4Vectors::width(index))) / 8) / stdBlockSize
    )
  }

  system2(wait = F, command = 'bwa', args = c('index',
                                              '-b',
                                              format(blockSizeMult*stdBlockSize,
                                                     scientific = FALSE), PAD))



  ##################################################################### find TIS
  if(verbose){message('Scanning insertion-sites')}
  TIS = NULL
  if(ITR == "PiggyBac"){
    targetInsertionSite = 'TTAA'
  } else if(ITR == "SleepingBeauty"){
    targetInsertionSite = 'TA'
  }

  TIS = Biostrings::vmatchPattern(targetInsertionSite, subject = index)

  ############################################################## add name to TIS
  TIS = lapply(seq_along(TIS), function(i){
    seqnames = as.character(names(TIS)[i])
    DF = as.data.frame(TIS[[i]])
    if(nrow(DF) == 0){
      NULL
    } else{
      cbind(seqnames, DF[,1:2])
    }
  })

  TIS = dplyr::bind_rows(TIS)
  TIS$name = targetInsertionSite
  TIS$score = 0
  TIS$strand = "."
  TIS = GenomicRanges::GRanges(TIS  )

  ####################################################################### export
  TISpath = gsub(x = PAD, pattern = ".fa.gz", replacement = ".tis", fixed = T)
  rtracklayer::export.bed(TIS, TISpath)

  ######################################################################## class
  thisTagMepprIndex = structure(list(index = PAD,
                         ITR = ITR,
                         TIS = TIS,
                         targetInsertionSite = targetInsertionSite,
                         insertName = names(transposonSeq),
                         seqinfo = FAI),
                    class = c("tagMepprIndex", "list"))

  ####################################################################### return

  return(thisTagMepprIndex)
}

#' loadIndex
#'
#' Load a hybrid reference genome
#'
#' @author Robin H. van der Weide, \email{r.vd.weide@nki.nl}
#' @param indexPath The path to the *_tagMeppRindex.fa.gz made by
#' \code{\link{makeIndex}}.
#'
#' @examples
#' \dontrun{
#' reference_hg19_PB = makeIndex(indexPath = '/home/A.Dent/analysis42/',
#'                               bsgenome = 'BSgenome.Hsapiens.UCSC.hg19',
#'                               ITR = 'PiggyBac')
#'
#' # stuff happens and reference_T42 is no longer loaded, so lets load it:
#' reference_hg19_PB = loadIndex("BSgenome.Hsapiens.UCSC.hg19_PiggyBac_tagMeppRindex.fa.gz")
#' }
#' @return A tagMepprIndex-object with the following slots:
#'
#' \describe{
#' \item{index}{The path to the .fa.gz of hybrid reference genome}
#' \item{ITR}{The name of the used transposon-system or ITR.fa}
#' \item{TIS}{A GRanges-object of all target insertion sites.}
#' \item{targetInsertionSite}{The targetInsertionSite}
#' \item{insertName}{The name of the transposon in the hybrid reference}
#' \item{seqinfo}{A seqinfo object of all chromosomes in index.}
#'}
#'
#' @importFrom GenomicRanges seqnames
#' @importFrom rtracklayer import.bed
#' @export
#'
loadIndex = function(indexPath){

  # check if fasta indeed exists
  if(!file.exists(indexPath)){
    stop(paste0(indexPath, " does not exist.\nPlease run makeIndex() first!"))
  }

  splitIndex = strsplit(basename(indexPath), split = "_")[[1]]
  ITR = splitIndex[length(splitIndex)-1]

  # check if TISpath indeed exists
  insertName = NULL
  TISpath = gsub(x = indexPath, pattern = ".fa.gz", replacement = ".tis", fixed = T)

  if(file.exists(TISpath)){
    TIS = rtracklayer::import.bed(TISpath)
    # first chrom is transposon name
    insertName = as.character(seqnames(TIS[1]))
  } else {
    stop(paste0(TISpath, " does not exist.\nPlease run makeIndex() first!"))
  }

  targetInsertionSite = S4Vectors::mcols(TIS[1])$name


  ###################################################################### seqinfo
  FAIPAD =  gsub(x = indexPath, pattern = ".fa.gz",
                 replacement = ".fa.gz.fai", fixed = T)

  FAIdf = read.delim(FAIPAD)
  FAI = GenomeInfoDb::Seqinfo(rownames(FAIdf), seqlengths=FAIdf$seqlengths,
                              isCircular=FAIdf$isCircular, genome=FAIdf$genome)



  ######################################################################## class
  thisTagMepprIndex = structure(list(index = indexPath,
                                     ITR = ITR,
                                     TIS = TIS,
                                     targetInsertionSite = targetInsertionSite,
                                     insertName =insertName,
                                     seqinfo = FAI),
                                class = c("tagMepprIndex", "list"))

  ####################################################################### return

  return(thisTagMepprIndex)
}
