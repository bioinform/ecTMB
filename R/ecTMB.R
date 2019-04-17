
"_PACKAGE"

#' MutSet class
#'
#' @docType class
#' @import R6
#' @return Object of class \code{MutSet}
#' @format An \code{\link{R6Class}} object.
#' @examples
#' \dontrun{
#' MutSet$new()
#' }
#' @slot mut table from MAF
#' @slot gid list of gene id for analysis
#' @slot exome table contains all possible mutations
#' @slot exomeGene table contains all posiible mutations for each gene.
#' @slot covar table contains covariates for each gene
#' @slot mutContext The context of 96 mutations
#' @slot samples data.frame store sample info. if BED colunm is specified,
#' the background process will be using the regions limited to bed file.
#' @slot incosistantAnno Variants with incosistant gene annotation from
#' GSMuta reference data .
#' \describe{
#'   \item{\code{new()}}{Create a MutSet }
#'   \item{\code{nonsilent()}}{Get nonsilent mutations table}
#'   \item{\code{silent()}}{Get silent mutations table}
#'   \item{\code{nc()}}{Get noncoding mutations table}
#' }
#'
#' @name MutSet
#' @docType class
#' @exportClass MutSet
#'
MutSet = R6::R6Class(
  'MutSet',
  portable = FALSE,
  public = list(
    mut        = NULL,
    gid        = NULL,
    exome      = NULL,
    exomeGene  = NULL,
    covar      = NULL,
    mutContext = NULL,
    incosistantAnno = NULL,
    samples    = NULL,

    nonsilent = function(gid) {
      mut[mut$Sil == 1,]

    },

    silent = function() {
      mut[mut$Sil == 2,]

    },

    nc = function() {
      mut[mut$Sil == 3,]
    },

    get_nonsil_passengers = function(fraction){
      nosam = aggregate(data=mut[mut$Sil == 1,], Tumor_Sample_Barcode~ Ensembl_gene_id, function(x) length(unique(x)))
      mut.gene = nosam[with(nosam,order(Tumor_Sample_Barcode)),1]
      nonsilent_passenger_gene=gid[! gid %in% mut.gene[round(length(mut.gene)* fraction):length(mut.gene)]]
      return(nonsilent_passenger_gene)
    }

  )
)


#' read and load maf file and necessary reference file
#' @param mutf Path to maf file. Detail see \code{mafCleanup}.
#' @param exomef Path to exome file
#' @param covarf Path to covariable file
#' @param mutContextf Path to mutContext file.
#' @param samplef Path to sample file or data.frame. The column name to specify bed file must be BED
#' @param ref Path to reference genome
#' @export
#' @return a MutSet
readData = function(mutf, exomef, covarf, mutContextf, ref, samplef= NULL){
  # load mutation context file
  mutContext            = read.table(mutContextf)
  colnames(mutContext)  = c("triMut", "rev_triMut", "triMut_code", "rev_triMut_code", "tag")

  # load covariant file
  covar                 = read.table(covarf, header=T, row.names=2, sep="\t")

  # load exome file
  ## format of exome file
  ## pos <\t> tri-nucleotide<\t>A <\t> C <\t> G <\t> T <\t> gene <\t> amino_acid_pos/protein_length
  exome                = loadfile(exomef)
  colnames(exome)      = c("Chromosome","pos", "seq_code","A", "C", "G", "T","gid", "aa_pos");


  # get valide gid
  gid                  = intersect(rownames(covar)[covar$Chromosome %in% c(1:22,"X")],exome$gid)


  # read in maf file
  mafTable             = mafCleanup(mutf, gid = gid)
  mafTable             = retrieve_context(mafTable, ref)
  mafTable$tag         = mutContext[match(mafTable[,"Context"],mutContext[,"triMut_code"]),"tag"]

  # remove inconsistant annotation
  inconMut             = Get_incon_mut(mafTable, exome)
  cat(sprintf("Number of inconsistant annotation mutation: %s out of total %s mutation \n",sum(inconMut), length(inconMut)))
  incosistantAnno      = mafTable[inconMut, ]
  mafTable             = mafTable[!inconMut, ]

  # load sample file
  if(is.null(samplef)){
    samples            = data.frame(SampleID = unique(mafTable$Tumor_Sample_Barcode),
                                    stringsAsFactors = FALSE)
  }else if(is.character(samplef)){
    samples            = read.delim(samplef, stringsAsFactors = FALSE)
  }else if("data.frame" %in% class(samplef)){
    samples            = samplef
  }

  ## check all samples contain at least one mutation.
  samples              = checkMutC(samples, mafTable)



  # load exomeGene file
  ## format of exome.gene #################
  ## EnsembleGeneID <\t> tri-nucleotide+change <\t> consequence <\t> count
  # consequence coding
  #  0 - silent
  #  1 - miss-sense
  #  2 - nonsense
  #  3 - nonstop
  #  4 - TSS
  #  5 - splice

  if("BED" %in% colnames(samples)){
    cat("Bed file is specified for samples\n")
    if(length(unique(samples$BED)) == 1){
      cat("\tAll samples have the same bed regions\n")
      exomeGene = get_exomeGene(exome, Bed = unique(samples$BED), mutContext = mutContext)
      gid       = intersect(gid, exomeGene$gid)  ## need to update gid to bed file.
      MutBed    = data.frame(Chromosome = mafTable[, "Chromosome"],
                             Start = mafTable[, "Start_Position"],
                             End = mafTable[,"End_Position"])
      subMut    = OverLap(MutBed,  regions= unique(samples$BED))
      cat(sprintf("Reduce number of mutation from: %s to %s \n",
                  nrow(MutBed), length(subMut)))
      mafTable  = mafTable[subMut,]
      ## check all samples contain at least one mutation.
    }else{
      cat("\tSamples' bed file are different. exomeGene for each sample will be generated\n")
      exomeGene = lapply(samples$BED, function(x){get_exomeGene(exome, Bed = x, mutContext = mutContext)})
      gids      = lapply(exomeGene, function(x){intersect(gid, x$gid)})
      names(exomeGene)  = names(gids) = samples$SampleID
      gid       = gids
      tmp       = lapply(split(mafTable, mafTable$Tumor_Sample_Barcode),
                         function(x){MutBed    = data.frame(Chromosome = x[, "Chromosome"],
                                                        Start = x[, "Start_Position"],
                                                        End = x[,"End_Position"])
                                     subMut    = OverLap(MutBed, regions= unique(samples$BED))
                                     out       = x[subMut,]
                                     return(out)})
      cat(sprintf("Reduce number of mutation from: %s to %s \n", nrow(mafTable), nrow(do.call(rbind, tmp))))
      mafTable  = do.call(rbind, tmp)
    }
  }else{
    cat("Bed file is not specified for samples. exomeGene for whole exome will be generated.\n")
    exomeGene = get_exomeGene(exome,  mutContext = mutContext)
  }

  samples   = checkMutC(samples, mafTable)



  # generate mutSet
  sset                 = MutSet$new()
  sset$mut             = mafTable
  sset$gid             = gid
  sset$exome           = exome
  sset$exomeGene       = exomeGene
  sset$covar           = covar
  sset$mutContext      = mutContext
  sset$samples         = samples
  sset$incosistantAnno = incosistantAnno

  cat(sprintf(paste0("Total gene analyzed: ",length(gid))),"\n")

  return(sset)
}

#'calMut
#' provide observed and predicted mutation rate
#' @param Data MutSet class. Detail see \code{MutSet}.
#' @param params output from fit_model. Detail see \code{fit_model}
#' @param sampleN Sample names
#' @param gids gene IDs
#' @param bed path to bed file.
#' @importFrom  data.table setkey
#' @export
#' @return a summary of predicted and expected mutation rate.
calMut = function(Data, params, sampleN = NULL, gids = NULL, bed = NULL){
  subData                      = Data$clone()
  if(is.null(sampleN)) sampleN = unique(subData$samples$SampleID)
  if(!is.null(bed)){
    exomeGene = get_exomeGene(subData$exome, Bed = bed, mutContext = subData$mutContext)
    gid       = intersect(subData$gid, exomeGene$gid)  ## need to update gid to bed file.
    MutBed    = data.frame(Chromosome = subData$mut[, "Chromosome"],
                           Start = subData$mut[, "Start_Position"],
                           End = subData$mut[,"End_Position"])
    subMut    = OverLap(MutBed,  regions= bed)
    subData$mut  = subData$mut[subMut,]
    subData$exomeGene = subData$exomeGene
    subData$gid  = subData$gid
  }
  if(is.null(gids))    gids = subData$gid

  ### get parameters from modeling
  gene.mu              = params[["geneMu"]]
  gene.phi             = params[["genePhi"]]
  o_scale              = params[["o_scale"]]
  MRtriProb            = params[['MRtriProb']]
  zero_p               = params[["geneZeroP"]]
  if(!is.null(bed)){ ## if bed file provided, parameter need to recalculated.
    MR_p               = MRtriProb$MR_p     # relative mutation frequency for each tri-nucleotide context #
    ### get gene.prob
    geneProb           = getGeneBgProb(subData$exomeGene, MR_p)
    # count number of mutation per patient#
    mutPerP            = count_Mut(subData$mut)
    mutPerP$count      = mutPerP$count/get_glen(subData$exomeGene, selGid = names(subData$gid)) * 1000000
    ### get gene length
    geneLen            = get_glen(subData$exomeGene, byGene= TRUE)
  }else{
    ### get gene.prob
    geneProb           = params[["geneProb"]]
    ### get mutation per pateint
    mutPerP            = params[["mutPerP"]]
    ### get gene length
    geneLen            = params[["geneLen"]]

  }



  if("data.table" %in% class(geneProb)){
    geneProb           = data.table(geneProb)
  }
  setkey(geneProb, gid, consequence)


  ## get mutation for non-silent/silent mutation ##
  if(!is.null(zero_p)){
    ## zero inflated model
    mut_pre_nonsil       = as.matrix((1-zero_p) * gene.mu * ((geneProb[J(gids,1),]$prob + geneLen[gids]*MRtriProb$MR_p["indel"]) %*%
                                                  t(as.matrix(mutPerP[match(sampleN, mutPerP$Tumor_Sample_Barcode),"count"]))))/exp(o_scale)
    mut_pre_sil          = as.matrix((1-zero_p) * gene.mu * ((geneProb[J(gids,0),]$prob) %*%
                                                  t(as.matrix(mutPerP[match(sampleN, mutPerP$Tumor_Sample_Barcode),"count"]))))/exp(o_scale)
  }else{
    mut_pre_nonsil       = as.matrix(gene.mu * (( geneProb[J(gids,1),]$prob + geneLen[gids]*MRtriProb$MR_p["indel"]) %*%
                                                  t(as.matrix(mutPerP[match(sampleN, mutPerP$Tumor_Sample_Barcode),"count"]))))/exp(o_scale)
    mut_pre_sil          = as.matrix(gene.mu * ((geneProb[J(gids,0),]$prob) %*%
                                                  t(as.matrix(mutPerP[match(sampleN, mutPerP$Tumor_Sample_Barcode),"count"]))))/exp(o_scale)

  }

  rownames(mut_pre_nonsil) = rownames(mut_pre_sil) = names(gene.mu)
  colnames(mut_pre_nonsil) = colnames(mut_pre_sil) = sampleN

  ## get observed mutation count for non-silent and silent mutation ##
  mut_obs_sil          = count_Mut(Data$silent(), gid = gids, sampleN = sampleN)
  mut_obs_nonsil       = count_Mut(Data$nonsilent(), gid = gids, sampleN = sampleN)

  out                  = melt(mut_obs_sil)
  colnames(out)        = c("gid", "sample","mut_obs_sil")
  out$mut_obs_nonsil   = melt(mut_obs_nonsil)$value
  out$mut_pre_sil      = melt(mut_pre_sil)$value
  out$mut_pre_nonsil   = melt(mut_pre_nonsil)$value
  out$geneLen          = geneLen[as.character(out$gid)]
  out$geneMu           = gene.mu[as.character(out$gid)]
  out$geneProb0        = geneProb[J(as.character(out$gid),0),]$prob
  out$geneProb1        = geneProb[J(as.character(out$gid),1),]$prob
  return(out)
}

