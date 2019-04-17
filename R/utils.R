
#' loadfile
#' @param x file path
#' @return content of loaded file
#' @export
#' @examples
#' \dontrun{
#' loadfile(x)
#' }
#'

loadfile = function(x){
  temp_space          =  new.env()
  bar                 =  load(x, temp_space)
  object              =  get(bar, temp_space)
  rm(temp_space)
  return(object)
}

#' OverLap
#' @param Bed targeted bed file. data.frame contain colunms "Chromosome" "Start" "End"
#' @param regions path to bed file or data.frame
#' @param base0 A boolen if it is 0 based.
#' @return a list of row number which targetBed overlap with regions.
#' @importFrom GenomicRanges makeGRangesFromDataFrame findOverlaps
#' @export
#' @examples
#' \dontrun{
#'   exomebed        = data.frame(Chromosome = exome[, "Chromosome"],
#'   Start = exome[, "pos"],
#'   End = exome[,"pos"])
#'   OverLap(exomebed, regions)
#' }
#'
OverLap = function(Bed, regions, base0 = FALSE){
  if(class(regions) == "character"){
    if(grepl(".bed$", regions)){
      regions      = makeGRangesFromDataFrame(read.delim(regions, stringsAsFactors = F),
                                              keep.extra.columns=FALSE,
                                              ignore.strand=TRUE,
                                              start.field="Start",
                                              end.field=c("End"),
                                              strand.field="strand",
                                              starts.in.df.are.0based=FALSE)
    }
  }else if(class(regions) != "GRanges"){
    stop("regions should be either GRanges or path to bed/vcf file")
  }

  if(any(grepl("chr", as.character(seqnames(regions)))) & ! any(grepl("chr", Bed$Chromosome))){
    Bed$Chromosome      = paste0("chr",Bed$Chromosome)
  }


  exomeGR               = makeGRangesFromDataFrame(Bed,
                                           keep.extra.columns=FALSE,
                                           ignore.strand=TRUE,
                                           start.field="Start",
                                           end.field=c("End"),
                                           starts.in.df.are.0based=FALSE)

  over          = suppressWarnings(findOverlaps(exomeGR, regions))
  subexome      = unique(queryHits(over))
  return(subexome)
}



#' Get exomeGene
#' @param exome exome data
#' @param Bed path to bed file or data.frame
#' @param mutContext data.frame contain mutation context.
#' @importFrom data.table data.table
#' @return exomeGene file
#' @export
#' @examples
#' \dontrun{
#' get_exomeGene(exome, Bed, mutContext)
#' }
#'
get_exomeGene = function(exome, Bed = NULL, mutContext){
  if(is.null(Bed)){
    subexome         = exome
  }else{
    exomebed         = data.frame(Chromosome = exome[, "Chromosome"],
                                 Start = exome[, "pos"],
                                 End = exome[,"pos"])
    subexome         = exome[OverLap(exomebed, Bed),]
  }
  tmp              = data.table(data.frame(gid = subexome[,"gid"],
                                           triMut_code = paste0(subexome[,"seq_code"],1),
                                           consequence = subexome[,"A"],
                                           pos = subexome[,"pos"],
                                           stringsAsFactors = FALSE))
  a                = data.frame(tmp[, length(pos), by = 'gid,triMut_code,consequence'])


  tmp              = data.table(data.frame(gid = subexome[,"gid"],
                                           triMut_code = paste0(subexome[,"seq_code"],4),
                                           consequence = subexome[,"C"],
                                           pos = subexome[,"pos"],
                                           stringsAsFactors = FALSE))
  b                = data.frame(tmp[, length(pos), by = 'gid,triMut_code,consequence'])


  tmp              = data.table(data.frame(gid = subexome[,"gid"],
                                           triMut_code = paste0(subexome[,"seq_code"],3),
                                           consequence = subexome[,"G"],
                                           pos = subexome[,"pos"],
                                           stringsAsFactors = FALSE))
  c                = data.frame(tmp[, length(pos), by = 'gid,triMut_code,consequence'])



  tmp              = data.table(data.frame(gid = subexome[,"gid"],
                                           triMut_code = paste0(subexome[,"seq_code"],2),
                                           consequence = subexome[,"T"],
                                           pos = subexome[,"pos"],
                                           stringsAsFactors = FALSE))
  d                = data.frame(tmp[, length(pos), by = 'gid,triMut_code,consequence'])

  ## note replace N with 0
  final            = rbind(a,b,c,d)
  final$triMut_code = sub("N", "0", final$triMut_code)
  final$triMut_code=as.numeric(final$triMut_code)
  final            = final[final$triMut_code %in% mutContext[,"triMut_code"],]



  colnames(final)  = c("gid", "triMut_code", "consequence", "count")
  final$tag        = mutContext[match(final$triMut_code, mutContext$triMut_code), "tag"]

  return(final)
}


#' Count mutatin per patient or per mutation tri-nuclear context.
#' @param mut data frame contain mutation info which must contain 'Start_Position' and 'Tumor_Sample_Barcode' or
#' exomeGene - data.frame.
#' @param gid default is NULL. provide gid list
#' @param sampleN default is NULL. Provide list of sample name. This parameter is only useful when byContext is FALSE
#' @param byContext A boolen. If TRUE, the number of mutation per tri-nucleotide context will be reported.
#' If FALSE, the number of mutation per patient will be reported. Default is FALSE.
#' @return mutCount data.frame contain number of mutation for each patient or
#' a vector contains number of mutation per tri-nucleotide context
#' @importFrom data.table data.table
#' @export
#' @examples
#' \dontrun{
#' count_Mut(mut)
#' }
count_Mut = function(mut, gid = NULL, sampleN = NULL, byContext = FALSE){
  if(byContext){
    if(is.null(gid)) gid = unique(mut[,"Ensembl_gene_id"])

    mutCount        = rep(0,96)
    if("Context" %in% colnames(mut)){
      mutab.snp     = mut[ mut[,"Context"] != "INDEL" & mut[,"Ensembl_gene_id"] %in% gid ,]
      x             = table(mutab.snp[,"tag"])
      y             = x[match(1:96,names(x))]
    }else{
      mutab.snp     = data.table(mut[  as.character(mut[,"gid"]) %in% gid ,])
      x             = mutab.snp[, sum(count), by = 'tag']
      y             = x[match(1:96, x$tag)]$V1
    }

    y[is.na(y)]     = 0
    names(y)        = 1:96
    mutCount        = mutCount + y

  }else{
    ## report per patient
    if(is.null(gid) ){
      if(is.null(sampleN)) sampleN = as.character(unique(mut$Tumor_Sample_Barcode))
      mutCount      = data.frame(Tumor_Sample_Barcode = sampleN, count = 0, stringsAsFactors = F)
      rownames(mutCount) = sampleN
      count        = aggregate(Start_Position ~ Tumor_Sample_Barcode, mut, length)
      mutCount[count$Tumor_Sample_Barcode, "count"] = count$Start_Position
    }else if((! is.null(gid))){
      count           = matrix(0, nrow = length(gid), ncol = length(sampleN))
      rownames(count) = gid
      colnames(count) = sampleN
      mutCount = aggregate(Start_Position ~ Ensembl_gene_id + Tumor_Sample_Barcode, mut,length)
      for (i in 1:nrow(mutCount)){
        if(as.character(mutCount[i,1]) %in% gid & as.character(mutCount[i,2]) %in% sampleN){
          count[as.character(mutCount[i,1]),as.character(mutCount[i,2])] = mutCount[i,3]
        }
      }
      mutCount        = count
    }
    if("Start_Position" %in% colnames(mutCount)){ colnames(mutCount)[which(colnames(mutCount) %in% "Start_Position")] = "count"}
  }
  return(mutCount)
}


#' get_glen
#' @description  get gene length
#' @param exomeGene data.frame
#' @param selGid a vector of selected genes.
#' @param byGene a boolen.If TRUE, a vector of length for each gene will be reported
#' If FALSE, a number of total length for all gene will be reported. Default is FALSE.
#' @export
#' @importFrom data.table data.table
#' @return a number or a vector
#'
get_glen = function(exomeGene, selGid = NULL, byGene= FALSE){
  if(byGene){
    if(!is.null(selGid)){
      u        = data.table(exomeGene[exomeGene$gid %in% selGid, ])
    }else{
      u        = data.table(exomeGene)
    }
    x          = u[, sum(count), by = 'gid']
    out        = x$V1/3
    names(out) = x$gid

  }else{
    if(is.null(selGid)){
      out       = sum(exomeGene[ "count"])/3
    }else{
      out       = sum(exomeGene[exomeGene$gid %in% selGid, "count"])/3
    }
  }
  return(out)

}


## ind: default is c(1,2,3,4,5) which are nonsil mutation
#' getGeneBgProb
#' @description  get background gene muation probabilty without regression
#' @param exomeGene data.frame
#' @param prob probability of each tri-nucleotides mutation context
#' @param consequences default is c(1,2,3,4,5) which are nonsil mutation.
#' @param gid Gene ID list.
#' @importFrom data.table data.table setkey
#' @export
#' @return data.frame background gene muation probabilty without regression
#'
getGeneBgProb = function(exomeGene, prob, consequences = c(1,2,3,4,5), gid = NULL){
  if(!is.null(gid)){
    exomeGene        = exomeGene[exomeGene$gid %in% gid, ]
  }
  gids               = unique(exomeGene$gid)
  df                 = data.table(data.frame(gid = exomeGene[,"gid"],
                                             prob = prob[as.numeric(exomeGene[,"tag"])]* as.numeric(exomeGene[,"count"]),
                                             consequence = as.numeric(exomeGene[,"consequence"] %in% consequences)))
  gene_background_p  = df[,sum(prob), by = 'gid,consequence']
  setkey(gene_background_p, gid, consequence)
  colnames(gene_background_p)[3] = "prob"
  out                = data.table(rbind( gene_background_p[J(gids, 0),], gene_background_p[J(gids, 1),]))
  probmin0           = min(out$prob[out$consequence == 0], na.rm = T)
  probmin1           = min(out$prob[out$consequence ==1 ], na.rm = T)
  out$prob[is.na(out$prob) & out$consequence == 0] = probmin0
  out$prob[is.na(out$prob) & out$consequence == 1] = probmin1
  setkey(out, gid, consequence)
  return(out)
}


#' getsubData
#' @description  get subset of orginal MutSet
#' @param Data MutSet
#' @param sampleN A list of sample names.
#' @param gid A list of gene id.
#' @export
#' @return MutSet a subset of orginal MutSet
#
getsubData = function(Data, sampleN = NULL, gid = NULL){
  subData                 = Data$clone()
  if(!is.null(sampleN)){
    subData$mut           = subData$mut[Data$mut$Tumor_Sample_Barcode %in% sampleN,]
    subData$samples       = subData$samples[Data$samples$SampleID %in% sampleN,]
  }

  if(!is.null(gid)){
    subData$mut           = subData$mut[subData$mut$Ensembl_gene_id %in% gid,]
    subData$gid           = subData$gid[subData$gid %in% gid]
    subData$covar         = subData$cover[gid, ]
    subData$exomeGene     = subData$exomeGene[subData$exomeGene$gid %in% gid, ]
    subData$exome         = subData$exome[subData$exome$gid %in% gid, ]
  }

  return(subData)
}



#' check each sample have at least one mutation.
#' @description check each sample have at least one mutation.
#' @param sample data.frame of sample
#' @param mut Mutation table
#' @return restrict sample to only the ones that contain at least one mutation
#
checkMutC = function(sample, mut){
  if(!all(unique(sample$SampleID) %in% mut$Tumor_Sample_Barcode)){
    cat(sprintf("Total %s out of %s samples with at least one mutation detected.\nOnly the ones with at least one mutation will be used.\n",
                     length(unique(mut$Tumor_Sample_Barcode)), length(unique(sample$SampleID))))
    sample                = sample[sample$SampleID %in% mut$Tumor_Sample_Barcode,]
  }
  return(sample)
}


#' getEnsemblID
#' @description  get ensembl ID from gene symbol
#' @param geneID gene symbol
#' @param geneInfof store the info symbol and ensembl ID.
#' @importFrom limma alias2Symbol
#' @export
#' @return a vector of ensemble ID
#
getEnsemblID = function(geneID, geneInfof = "~/Data/GSMuta_data/GSMutaRData/extdata/ensembl_92_exon_pos.hg38.rda"){
  geneinfo         = loadfile(geneInfof)
  out              = geneinfo$ensembl_gene_id[match(geneID, geneinfo$hgnc_symbol)]
  for(i in which(is.na(out))){
    gene           = geneID[i]
    alias          = alias2Symbol(gene)
    if(length(alias) > 1) alias = alias[1]
    if(gene != alias){
      out[i]       = geneinfo$ensembl_gene_id[match(alias, geneinfo$hgnc_symbol)]
    }
  }
  return(out)
}


#' getGeneSymbol
#' @description  get gene symbol from ensembl ID
#' @param ensemblID ensemble ID
#' @param geneInfof store the info symbol and ensembl ID.
#' @export
#' @return a vector of gene symbol
getGeneSymbol = function(ensemblID, geneInfof = "~/Data/GSMuta_data/GSMutaRData/extdata/ensembl_92_exon_pos.hg38.rda"){
  geneinfo         = loadfile(geneInfof)
  out              = geneinfo$hgnc_symbol[match(ensemblID, geneinfo$ensembl_gene_id)]
  return(out)
}


