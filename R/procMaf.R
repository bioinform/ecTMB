#' maf.cleanup
#' @param maf Input maf file. It can be a path to either rda or tsv file or data.frame
#' @param extraCols Default is NULL which no extra column will be included. If "all", all columns
#' will be output. If an array, only specificed columns will be reported.
#' @param keepNoncoding A boolen. If TRUE, the noncoding variants will be keeped. Default is FALSE,
#' which noncoding variants will be removed.
#' @param nonSilTypes To specify what variant classification should be considered as nonsilent variants.
#' Default is NULL which the default list of classification will be used.
#' @param SilTypes To specify what variant classification should be considered as silent variants.
#' Default is NULL which the default list of classification will be used.
#' @param ncTypes To specify what variant classification should be considered as variants in noncoding regions.
#' Default is NULL which the default list of classification will be used.
#' @param save A boolen. If TRUE, the file will be saved, FALSE, the file will not be saved
#' @param fn string, the output file name
#' @param gid Gene ID list.
#' @return cleaned up maf file.
#' @export
#' @examples
#' \dontrun{
#' maf.cleanup(x, output)
#' }
#'
mafCleanup = function(maf, extraCols = NULL, keepNoncoding = FALSE,
                      save = FALSE, fn = "./output.tsv",
                      nonSilTypes = NULL, SilTypes = NULL, ncTypes = NULL,
                      gid = NULL){

  ## variant classification category--------
  if(is.null(nonSilTypes)){
    nonSilTypes = c("Missense_Mutation","Nonsense_Mutation",
                    "Frame_Shift_Del","Frame_Shift_Ins",
                    "In_Frame_Ins","In_Frame_Del")
  }

  if(is.null(SilTypes)){
    SilTypes = c( "Silent")
  }

  if(is.null(ncTypes)){
    ncTypes = c("3'Flank",  "3'UTR", "5'Flank",  "5'UTR", "IGR", "Intron", "Splice_Region", "RNA", "Nonstop_Mutation",
                  "Translation_Start_Site","De_novo_Start_InFrame", "De_novo_Start_OutOfFrame", "Splice_Site")
  }

  ## load file  ---------
  if(length(class(maf)) ==1 && class(maf) == "character"){
    if(grepl(".rda$", maf)){
      mafTable = loadfile(maf)
    }else{
      mafTable = read.delim(maf, header=TRUE, sep="\t", stringsAsFactors=FALSE, skip=1)
    }
  }else{
    mafTable   = maf
  }

  mafTable     = as.matrix(apply(mafTable,2,as.character))
  mafTable     = mafTable[!is.na(mafTable[,"Chromosome"]),] #### remove positions without chromosome info.
  if( sum(mafTable[, "FILTER"] %in% c("PASS", ".")) != nrow(mafTable)) {
	print("Need to remove non-PASS calls")
	print(table(mafTable[, "FILTER"]))
	mafTable     = mafTable[mafTable[, "FILTER"] == "PASS",]  #### Only use pass calls.
  }
  typeInd      = which(colnames(mafTable)=="Variant_Type")
  classInd     = which(colnames(mafTable)=="Variant_Classification")

  ## sanity check -----
  ### Silent variant should not be INDEL
#  if(any( mafTable[,classInd] %in% c("Silent") & mafTable[,typeInd] %in% c("INS", "DEL") ) ){
#    stop("Maf file contain variants which are silent but belong to INDEL in type. Please double check maf file.")
#  }

  ### "De_novo_Start_InFrame", "De_novo_Start_OutOfFrame" are not longer allowed. Just a warning
  if(any( mafTable[,classInd] %in% c("De_novo_Start_InFrame", "De_novo_Start_OutOfFrame")  ) ){
    warning("Variant_Classification field contain De_novo_Start_InFrame or De_novo_Start_OutOfFrame, which is no longer allowed. It will be treat as In_frame and frame_shift respectively")
  }

  ### check extra type
  if(any( !mafTable[,typeInd] %in% c("SNP", "DNP", "TNP", "ONP", "INS", "DEL") ) ){
    warning("There are extra type beside SNP, DNP, TNP, ONP, INS, DEL in variant_type field. It will be convert to either of above categories")
    runConvExtra = TRUE
  }else{
    runConvExtra = FALSE
  }

  ### check if variants are all on 1-22 and X, Y chromosome. Others will be removed
  if(any( !(mafTable[,"Chromosome"] %in% c(1:24,"X","Y") | mafTable[,"Chromosome"] %in% paste0("chr",c(1:24,"X","Y"))) ) ){
    warning("There are variants on genome contig other than 1-22 and X, Y. These variants will be removed.")
    mafTable      = mafTable[(mafTable[,"Chromosome"] %in% c(1:24,"X","Y") | mafTable[,"Chromosome"] %in% paste0("chr",c(1:24,"X","Y"))), ]
  }

  ### check if nonsilent mutation doesn't have gene info
  if(any( mafTable[,classInd] %in% nonSilTypes & is.na(mafTable[,"Gene"]) ) ){
    warning("There are variants are nonsilent variants but don't have gene info. These variants will be removed.")
    mafTable      = mafTable[!(mafTable[,classInd] %in% nonSilTypes & is.na(mafTable[,"Gene"])), ]
  }

  ## simplify maf type for indel -------
  FrameShiftTypes = c("Frame_Shift_Del", "Frame_Shift_Ins","Nonsense_Mutation",
                      "Splice_Site","Nonstop_Mutation","Translation_Start_Site", "De_novo_Start_OutOfFrame")
  InFrameType     = c("In_Frame_Del", "In_Frame_Ins","Missense_Mutation", "De_novo_Start_InFrame")

  FSrows          = mafTable[,typeInd] %in% c("INS", "DEL") & mafTable[, classInd] %in% FrameShiftTypes
  IFrows          = mafTable[,typeInd] %in% c("INS", "DEL") & mafTable[, classInd] %in% InFrameType

  mafTable[FSrows, typeInd] = "Frame_shift"
  mafTable[IFrows, typeInd] = "In_frame"

  ## deal with extra type
  if(runConvExtra){
    mafTable      = ConvertExtraType(mafTable)
  }

  ## format output
  if( !is.null(extraCols) && extraCols == "all"){
    defaultCols   = c("Ensembl_gene_id", "Gene", "Chromosome", "Start_Position","End_Position",
                      "Variant_Type","Reference_Allele", "Tumor_Seq_Allele1",
                      "Tumor_Seq_Allele2", "Tumor_Sample_Barcode")
    extraCols     = colnames(mafTable)[!colnames(mafTable) %in% defaultCols]
  }
  maffinal        = cbind(Ensembl_gene_id=as.character(mafTable[,"Gene"]),
                          mafTable[,c("Chromosome","Start_Position","End_Position",
                                      "Variant_Type","Reference_Allele",
                                      "Tumor_Seq_Allele1", "Tumor_Seq_Allele2",
                                      "Tumor_Sample_Barcode")],
                          Protein_position = mafTable[,"Protein_position"],
                          Variant_Classification = mafTable[,"Variant_Classification"],
                          mafTable[,extraCols])
  maffinal        = as.data.frame(maf_dnp_converter(maffinal), stringsAsFactors = FALSE)
  maffinal$Sil    = ifelse(maffinal[, "Variant_Classification"] %in% nonSilTypes, 1,
                           ifelse(maffinal[,"Variant_Classification"] %in% SilTypes, 2, 3))
  maffinal$Sil[maffinal$Variant_Type %in% c("Frame_shift","In_frame")] = 1

  if(!keepNoncoding){
    maffinal      = maffinal[maffinal$Sil != 3, ]
  }

  if(!is.null(gid)){
    maffinal = maffinal[maffinal[,"Ensembl_gene_id"] %in% gid,]
  }

  maffinal$Start_Position = as.numeric(maffinal$Start_Position)

  if(save) write.table(maffinal ,file = output,row.names=F,col.names=T,quote=F,sep="\t")

  return(maffinal)
}


#' ConvertExtraType
#' @param maf maf matrix
#' @return maf file with converted extra type
#' @export
#' @examples
#' \dontrun{
#' ConvertExtraType(maf)
#' }
#'
ConvertExtraType = function(maf){
  ExtraType       = which(!(maf[,"Variant_Type"] %in% c("SNP", "DNP", "TNP", "ONP", "Frame_shift", "In_frame", "INS", "DEL")))
  for (i in 1:length(ExtraType)){
    alleles = c(maf[ExtraType[i],"Reference_Allele"], maf[ExtraType[i], "Tumor_Seq_Allele1"],maf[ExtraType[i],  "Tumor_Seq_Allele2"])
    maxchar = max(nchar(alleles))
    if ("-" %in% alleles ) {
      maf[ExtraType[i], "Variant_Type"] = "Frame_shift"
      if (maxchar %% 3 == 0) maf[ExtraType[i], "Variant_Type"] = "In_frame"
    } else{
      maf[ExtraType[i], "Variant_Type"] = "SNP"
      if (maxchar > 1) maf[ExtraType[i], "Variant_Type"] = "DNP"
      if (maxchar > 2) maf[ExtraType[i], "Variant_Type"] = "TNP"
      if (maxchar > 3) maf[ExtraType[i], "Variant_Type"] = "ONP"
    }
  }
  return(maf)
}


#' maf_dnp_converter
#' @param mutab maf matrix
#' @return maf file with converted SNP
#' @export
#' @examples
#' \dontrun{
#' maf_dnp_converter(maf)
#' }
#'
maf_dnp_converter = function(mutab){
  mutab = as.matrix(mutab)
  mutab.snp= mutab[ mutab[,"Variant_Type"]=="SNP" ,]

  if(sum(  mutab[,"Variant_Type"]=="In_frame" | mutab[,"Variant_Type"]=="Frame_shift"  )>0 )
    mutab.indel=   mutab[ mutab[,4]=="In_frame" | mutab[,"Variant_Type"]=="Frame_shift",]   else mutab.indel=NULL

  mutab.dnp= mutab.tnp=mutab.onp=NULL
  if(sum(mutab[,"Variant_Type"]=="DNP")>0){
    mutab.dnp=  (as.matrix(mutab[ mutab[,"Variant_Type"]=="DNP" ,]))
    print(head(mutab.dnp))
    if (sum(mutab[,"Variant_Type"]=="DNP")==1) {mutab.dnp=  t(as.matrix(mutab[ mutab[,"Variant_Type"]=="DNP" ,])) }
    #### convert DNP to SNP ##############
    a= matrix(0, nrow= nrow(mutab.dnp)*2,ncol=ncol(mutab.dnp)   )
    colnames(a)=colnames(mutab)
    for(i in c(1,2,5,9:ncol(mutab.dnp))){
      a[,i]=  rep(mutab.dnp[,i] ,each=2)
    }
    # a[,2]=  rep(mutab.dnp[,2] ,each=2)
    # a[,4]=  rep(mutab.dnp[,5] ,each=2)
    # a[,8]=  rep(mutab.dnp[,9] ,each=2)
    # a[,9]=  rep(mutab.dnp[,10] ,each=2)
    # a[,11]=  rep(mutab.dnp[,11] ,each=2)

    a[2*(1:nrow(mutab.dnp))-1,3]=  mutab.dnp[,3]
    a[2*(1:nrow(mutab.dnp)),3]= as.numeric(mutab.dnp[,3]) + 1
    a[,4]= a[,3]
    a[,6]= unlist(  strsplit( mutab.dnp[,6]  ,""))
    a[,7]= unlist(  strsplit( mutab.dnp[,7] ,""))
    a[,8]= unlist(  strsplit( mutab.dnp[,8] ,""))
    a = a[(a[,6] == a[,7] & a[,7] == a[,8]) == FALSE,]
    mutab.dnp=a
  }
  if(sum(mutab[,"Variant_Type"]=="TNP")>0){
    mutab.tnp=  mutab[ mutab[,"Variant_Type"]=="TNP" ,]
    if (sum(mutab[,"Variant_Type"]=="TNP")==1) {mutab.tnp=  t(as.matrix(mutab[ mutab[,"Variant_Type"]=="TNP" ,])) }
    #### convert TNP to SNP ##############
    a=matrix(0, nrow= nrow(mutab.tnp)*3,ncol=ncol(mutab.tnp)   )
    colnames(a)=colnames(mutab)
    for(i in c(1,2,5,9:ncol(mutab.tnp))){
      a[,i]=  rep(mutab.tnp[,i] ,each=3)
    }
    # a[,1]=  rep(mutab.tnp[,1] ,each=3)
    # a[,2]=  rep(mutab.tnp[,2] ,each=3)
    # a[,4]=  rep(mutab.tnp[,4] ,each=3)
    # a[,8]=  rep(mutab.tnp[,8] ,each=3)
    # a[,9]=  rep(mutab.tnp[,9] ,each=3)
    # a[,10]=  rep(mutab.tnp[,10] ,each=3)

    a[3*(1:nrow(mutab.tnp))-2,3]=  mutab.tnp[,3]
    a[3*(1:nrow(mutab.tnp))-1,3]= as.numeric(mutab.tnp[,3]) + 1
    a[3*(1:nrow(mutab.tnp)),3]= as.numeric(mutab.tnp[,3]) + 2
    a[,4]= a[,3]
    a[,6]= unlist(  strsplit( mutab.tnp[,6]  ,""))
    a[,7]= unlist(  strsplit( mutab.tnp[,7] ,""))
    a[,8]= unlist(  strsplit( mutab.tnp[,8] ,""))
    a = a[(a[,6] == a[,7] & a[,7] == a[,8]) == FALSE,]
    mutab.tnp=a
  }
  if(sum(mutab[,"Variant_Type"]=="ONP")>0)
  {
    mutab.onp.all = mutab[ mutab[,"Variant_Type"]=="ONP" ,]
    if (sum(mutab[,"Variant_Type"]=="ONP")==1) {mutab.onp.all=  t(as.matrix(mutab[ mutab[,"Variant_Type"]=="ONP" ,])) }
    #### convert ONP to SNP ##############
    maxchar = apply(cbind(nchar(mutab.onp.all[,"Reference_Allele"]),
                          nchar(mutab.onp.all[,"Tumor_Seq_Allele1"]),
                          nchar(mutab.onp.all[,"Tumor_Seq_Allele2"])),1,max)
    umaxchar = unique(maxchar)
    mutab.onp = mutab.onp.all
    b = NULL
    for (k in 1:length(umaxchar)){
      mutab.onp = mutab.onp.all[maxchar==umaxchar[k],]
      if (length(mutab.onp) == 8){mutab.onp = t(as.matrix(mutab.onp))}
      a= matrix(0, nrow= nrow(mutab.onp)*umaxchar[k],ncol=ncol(mutab.onp)   )
      colnames(a)=colnames(mutab)
      for(i in c(1,2,5,9:ncol(mutab.onp))){
        a[,i]=  rep(mutab.onp[,i] ,each=umaxchar[k])
      }
      # a[,1]=  rep(mutab.onp[,1] ,each=umaxchar[k])
      # a[,2]=  rep(mutab.onp[,2] ,each=umaxchar[k])
      # a[,4]=  rep(mutab.onp[,4] ,each=umaxchar[k])
      # a[,8]=  rep(mutab.onp[,8] ,each=umaxchar[k])
      # a[,9]=  rep(mutab.onp[,9] ,each=umaxchar[k])
      # a[,10]=  rep(mutab.onp[,10] ,each=umaxchar[k])

      for (len in 1:umaxchar[k]){
        a[umaxchar[k]*(1:nrow(mutab.onp))-(umaxchar[k]-len),3]=  as.numeric(mutab.onp[,3]) + (len-1)
      }
      a[,4] = a[,3]
      a[,6]= unlist(  strsplit( mutab.onp[,6]  ,""))
      a[,7]= unlist(  strsplit( mutab.onp[,7] ,""))
      a[,8]= unlist(  strsplit( mutab.onp[,8] ,""))
      a = a[(a[,6] == a[,7] & a[,7] == a[,8]) == FALSE,]
      b = rbind(b,a)
    }
    mutab.onp = b
  }

  mutab=rbind(mutab.snp,mutab.dnp,mutab.tnp,mutab.onp,mutab.indel)
  ref=mutab[,"Reference_Allele"]      ### reference nucleotide
  mut=ref
  temp= (mutab[,"Tumor_Seq_Allele1"]!=ref)
  mut[temp]=mutab[temp,"Tumor_Seq_Allele1"]

  temp= (mutab[,"Tumor_Seq_Allele2"]!=ref)
  mut[temp]=mutab[temp,"Tumor_Seq_Allele2"]   ### mutated nucleotide

  mutab =  cbind(mutab,mut)
  colnames(mutab)[ncol(mutab)] = "Mut"

  return(mutab)
}

#' retrieve_context
#' @param mutab maf matrix
#' @param ref genome reference path.
#' @return Retrieve mutation context code for each variants
#' @export
#' @examples
#' \dontrun{
#' maf_dnp_converter(mutab, ref, codeRetrieve)
#' }
#'
retrieve_context = function(mutab, ref){
  outDir = "./"
  if(!all(grepl("chr", mutab$Chromosome))){mutab$Chromosome = paste0("chr", mutab$Chromosome)}

  code = file.path(system.file( "perl", package="ecTMB"), "Sequence_Retrieve.pl")
  muttmpfile = tempfile(c("mutab"), tmpdir =outDir, fileext = ".tsv" )

  write.table(mutab, file = muttmpfile, quote = F, sep = "\t", row.names = FALSE)
  # cmd = sprintf("export PATH=/sc1/groups/bfx-red/apps/datainsights/bedtools2/bin:$PATH;perl %s %s %s %s",
  #               code, muttmpfile,ref, sub(".tsv", "out.tsv", muttmpfile))
  cmd = sprintf("perl %s %s %s %s",
                code, muttmpfile,ref, sub(".tsv", "out.tsv", muttmpfile))
  system(cmd)
  out = read.delim(sub(".tsv", "out.tsv", muttmpfile), stringsAsFactors = FALSE)
  cmd = sprintf(" rm %s %s", muttmpfile , sub(".tsv", "out.tsv", muttmpfile))
  system(cmd)

  if(!all(grepl("chr", out[,"Chromosome"]))){
    out[,"Chromosome"] = paste0("chr", out[,"Chromosome"])
  }
  return(out)
}


#' Get_incon_mut
#' Get Variants with incosistant gene annotation from
#' GSMuta reference data
#' @param mutab maf matrix
#' @param exome exome file.
#' @return Get location of variants with incosistant gene annotation from
#' GSMuta reference data
#' @importFrom data.table setkey data.table
#' @export
#' @examples
#' \dontrun{
#' Get_incon_mut(mutab,exome)
#' }
#'
Get_incon_mut = function(mutab,exome){
  X        = unique(data.table(exome[,c("pos", "gid", "seq_code")]))
  setkey(X, pos, gid)
  which    = is.na(X[J(as.integer(mutab$Start_Position),mutab$Ensembl_gene_id),]$seq_code)
  return(which)
}
