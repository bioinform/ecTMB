
# esitmate the Gausian mixture #
#' ExGMM
#' @description fit the data distribution to a mixture of Gaussian mixture
#' @param x number of mutation
#' @param msi msi status
#' @param single A boolen. If True, only one class will be identified.
#' @importFrom mixtools normalmixEM
#' @importFrom dplyr %>% group_by summarize
#' @return A vector mutation count for 96 trinucleotide context
#' @export
#' @examples
#' \dontrun{
#' ExGMM(Data)
#' }
#'
ExGMM = function(x, msi = NULL, single = FALSE){
  if(!single){
    if(is.null(msi)){
      df          = data.frame(x = x, cluster = ifelse( x > log(50), 3,
                                                        ifelse(x > log(9), 2, 1)))
    }else{
      df          = data.frame(x = x, cluster = ifelse( msi %in% "MSI-H", 2,
                                                        ifelse(x > mean(x[msi %in% "MSI-H"]), 3, 1)))
    }
    if(sum(df$cluster > 1) < 10){
      df$cluster[df$cluster > 1] = 1
      k         = length(unique(df$cluster))
    }
    if(sum(df$cluster > 2) < 10){
      df$cluster[df$cluster > 2] = 2
      k         = length(unique(df$cluster))
    }
    summary     = df %>%
      group_by(cluster) %>%
      summarize(mu = mean(x), variance = var(x), std = sd(x), size = n())

    if(length(unique(df$cluster)) == 1){
      single = TRUE
    }else{
      mixmdl = normalmixEM(df$x, mu = summary$mu,  sigma = summary$std,  k = k, arbvar = T, arbmean = T, epsilon = 1e-03)
      plot(mixmdl,which=2)
      lines(density(df$x), lty=2, lwd=2)
    }
  }else{
    df          = data.frame(x = x)
  }

  if(single){
    rmoutlierD = remove_outliers(df$x)
    rmoutlierD = rmoutlierD[!is.na(rmoutlierD)]
    mixmdl = list(x = df$x, mu = mean(rmoutlierD),  sigma = sd(rmoutlierD), lambda = 1)
  }
  return(mixmdl)
}


# remove outliers #
#'remove_outliers
#' @description remove outlier based on 1.5 quantile
#' @param x number of mutation
#' @param na.rm how to deal with na
#' @param ... The extra parameters for quantile.
#' @return outlier is labeled as NA
#' @export
#' @examples
#' \dontrun{
#' remove_outliers(x)
#' }
#'
remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
  H <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  y
}

# predict TMB #
#' pred_TMB
#' @description fit the data distribution to a mixture of Gaussian mixture
#' @param x test mutSet object
#' @param params A list contain result from fit
#' @param WES test mutSet object which contain whole exome data which will be used to generate the truth.
#' @param prior A list contain prior parameters. prior = list(mu = prior_bs$mu, sigma = prior_bs$sigma, lambda = prior_bs$lambda)
#' @param mut.nonsil A boolen. If it is TRUE, the nonsilent mutation will be used for prediction
#' @param span The extension of low and high limit. High = high * span and low = low/span
#' @param low The lower limit.
#' @param high The higher limit.
#' @param gid_nonsil_p A vector of gid for passenger genes
#' @param method A string which specify the method to predict the TMB
#' @param cores A number which specifiy the nubmer of core to use for parallel computing.
#' @param bs Either all or nonsil.
#' @return A data.frame contain the result.
#' @export
#' @examples
#' \dontrun{
#' pred_TMB(Data)
#' }
#'
pred_TMB = function(x, params, WES = NULL, prior = NULL,
                    mut.nonsil = FALSE, gid_nonsil_p = NULL, method = "MLE", cores = 10, span = 1, low = 1, high = 10^4,
                    bs = "nonsil"){
  MR_p               = params$MRtriProb$MR_p     # relative mutation frequency for each tri-nucleotide context #
  o_scale            = params$o_scale
  sampleN            = x$samples$SampleID

  ### use panel genes for prediction
  gene_bg_p          = getGeneBgProb(x$exomeGene, MR_p)
  gid                = x$gid
  selGene            = gid
  if(is.null(gid_nonsil_p) & mut.nonsil){
    gid_nonsil_p     = selGene
    cat("Since gid_nonsil_p is NULL, %s percentage of genes' nonsilent mutationw will be used.\n\n")
  }
  if( mut.nonsil ){
    offset_nonsil    =  params$geneMu[selGene] * as.numeric(gene_bg_p[J(selGene,1),]$prob)
    offset_nonsil[!selGene %in% gid_nonsil_p] = 0
    offset_sil       = params$geneMu[selGene] * as.numeric(gene_bg_p[J(selGene,0),]$prob)
    offset           = offset_nonsil + offset_sil
    sil_mut_matrix   = unname(count_Mut(x$silent(), selGene, sampleN))
    nonsil_mut_matrix= unname(count_Mut(x$nonsilent(), selGene, sampleN))
    nonsil_mut_matrix[!selGene %in% gid_nonsil_p, ] = 0
    y                = nonsil_mut_matrix + sil_mut_matrix
  }else{
    offset             = params$geneMu[selGene] * as.numeric(gene_bg_p[J(selGene,0),]$prob)
    sil_mut_matrix     = unname(count_Mut(x$silent(), selGene, sampleN))
    y                  = sil_mut_matrix
  }

  rownames(y)          = selGene
  colnames(y)          = sampleN
  offsetS              = offset/exp(o_scale) ### rescale 0~1



  # cat(sprintf("In panel, %s patients contain 0 silent (or + nonsilent) mutations.\n", sum(colSums(y) == 0)))
  if(method == "MLE"){
    pred_bs_panel      = getTMBs(y, offsetS,  phi = params$geneDisp[selGene] ,
                                 cores, zero_p = params$geneZero_p[selGene],
                                 prior = prior, span = span, low = low, high = high)
  }else if( method == "MCMC"){
    # TODO
    # mcmc_panel         = Instan(y, offsetS, prior = prior, cores = cores,
    #                             zero_p = params$geneZero_p[selGene],
    #                             geneDisp  = params$geneDisp[selGene],
    #                             span = span, low = low, high = high)
    # pred_bs_panel      = mcmc_panel[sampleN, "mean"]
    # mcmc_panel         = mcmc_panel[sampleN, !colnames(mcmc_panel) %in% "mean"]
    # colnames(mcmc_panel) = paste0("panel_", colnames(mcmc_panel))
  }else{
    stop(sprintf("Prediction method should either MLE or MCMC. %s was provided", method))
  }
  # mut_panel_used     = colSums(y)

  if(bs == "all" ){
    count            = CalTMB(x, sampleN = sampleN, type = "all")
  }else if( bs == "nonsil" ){
    count            = CalTMB(x, sampleN = sampleN, type = "nonsil")
  }

  out                = data.frame(sample = sampleN,
                                  ecTMB_panel_TMB = pred_bs_panel,
                                  count_panel_TMB = count$count)



  ## ground truth
  if(!is.null(WES)){
    if(bs == "all" ){
      WES_TMB            = CalTMB(WES, sampleN = sampleN, type = "all")
    }else if( bs == "nonsil" ){
      WES_TMB            = CalTMB(WES, sampleN = sampleN, type = "nonsil")
    }
    out$WES_TMB          = WES_TMB$count
  }
  # if(method == "MCMC"){
  #   out              = cbind(out, mcmc_panel)
  #   out              = cbind(out, mcmc_wes)
  # }
  # cat(sprintf("Columns %s contain NA.\n", paste(colnames(out)[colSums(is.na(out)) > 0], collapse = ",")))
  # out[is.na(out)]    = 0
  return(out)
}


# predict TMB stat function #
#' getTMBs
#' @description predict TMB stat function
#' @param y A matrix for observed mutation count
#' @param offset offset
#' @param phi gene dispersion
#' @param zero_p gene zero fraction
#' @param prior prior parameters. If it is NULL, no prior will be used.
#' @param span The extension of low and high limit. High = high * span and low = low/span
#' @param low The lower limit.
#' @param high The higher limit.
#' @param cores A number which specifiy the nubmer of core to use for parallel computing.
#' @return A data.frame contain the result.
#' @export
#' @examples
#' \dontrun{
#' getTMBs(Data)
#' }
#'
getTMBs = function(y, offset, phi, cores, zero_p = 0, prior = NULL, span = 10, low = 10^(-4), high = 10^4){
  ngene = nrow(y)
  nsam  = ncol(y)
  if (length(phi) == 1) {
    phi = rep(phi, ngene)
  }

  if (length(zero_p) == 1) {
    zero_p = rep(zero_p, ngene)
  }
  getbs = function(x, offset, disp, zero_p, prior, span = 100, low = 10^(-4), high = 10^4){
    #print(1)
    ind0 = x == 0
    obj = function(bs){
      ## Pr(bi)
      p_bs       = 0
      for (i in 1:length(prior$mu)){
        p_bs     = p_bs + (prior$lambda[i] * 1/sqrt(2*pi*prior$sigma[i]) * exp(-(log(bs) - prior$mu[i])^2/(2 * prior$sigma[i]^2) ))
      }
      if(p_bs == 0) {
        lp_bs     = -10^(200)
      }else{
        lp_bs     = log(p_bs)
      }
      if(any(zero_p != 0)){
        # zero-inflated
        if(all(disp == 0)){
          # zero-inflated Poisson
          return(
            # -log(lp_bs) - ## prior section
            -lp_bs -sum( ( log( zero_p + (1 - zero_p) * exp(-bs * offset) ) )[ind0] ) - # zero section disp.g is zero_p
              sum( ( log(1 - zero_p) - (bs * offset)  + x*log(bs * offset) - lgamma(x+1) )[!ind0] )
          )
        }
      }else{
        if(all(disp == 0)){
          ## regular Poisson
          return( lp_bs - sum( x*log(bs * offset)-offset*bs -lgamma(x+1)) )
        }else{
          ## regular NB
          return ( lp_bs -sum( lgamma(x + 1/disp) - lgamma(1/disp) -lgamma(x+1) -
                                 1/disp*log(1 + bs*disp*offset) + x*( log(bs*offset) - log(1/disp+bs*offset) ) ) )
        }
      }
    }

    obj2 = function(bs){     ## no prior
      if(any(zero_p != 0)){
        # zero-inflated
        if(all(disp == 0)){
          # zero-inflated Poisson
          return(
            -sum( ( log( zero_p + ( 1 - zero_p) * exp(-bs * offset) ) )[ind0] ) - # zero section disp.g is zero_p
              sum( ( log(1-zero_p) - (bs * offset)  + x*log(bs * offset) -lgamma(x+1) )[!ind0] )
          )
        }
      }else{
        if(all(disp == 0)){
          ## regular Poisson
          return(  - sum( x*log(bs * offset)-offset*bs -lgamma(x+1)) )
        }else{
          ## regular NB
          return ( -sum( lgamma(x + 1/disp) - lgamma(1/disp) -lgamma(x+1) -
                                 1/disp*log(1+bs*disp*offset) + x*( log(bs*offset) - log(1/disp+bs*offset) ) ) )
        }
      }
    }

    if(is.null(prior)){
      return(optimize(obj2, interval=c(low/span, high * span))$minimum) ## no prior
    }else{
      return(optimize(obj, interval=c(low/span, high * span))$minimum)  ## with prior
    }
  }
  x <- lapply(apply(y, 2, FUN=list), unlist)
  bs <- unlist(mclapply(x, function(x) getbs(x, offset = offset, disp = phi, zero_p = zero_p, prior = prior,span = span, low = low, high = high) ,mc.cores=cores))
  return(bs)
}


# Classify the group of TMB #
#' AssignClass
#' @description predict TMB stat function
#' @param x predicted TMB
#' @param prior prior parameters.
#' @param type If 'low_high", only class 2 and 3 will be grouped to high.
#' If 'exact', exact class will be reported
#' @param add1 A boolen. If prior was defined with log( x + 1), then it should be TRUE.
#' @return A vector of class
#' @export
#' @examples
#' \dontrun{
#' assignClass(Data)
#' }
#'

assignClass = function(x, prior, type = "exact", add1 = TRUE){
  cla       = function(x, prior){
    pr      = c()
    for(i in 1:length(prior$mu)){
      pr    = c(pr, 1/(prior$sigma[i] * sqrt(2 * pi)) * exp(-(log(x) - prior$mu[i])^2/(2*prior$sigma[i])))
      # pr    = c(pr, 1/(prior$sigma[i] * sqrt(2 * pi)) * exp(-(log(x) - prior$mu[i])^2/(2*prior$sigma[i]))* prior$lambda[i])
    }
    ## for sure low
    if(x < prior$mu[1] & pr[1] < pr[2]){
      pr[1] = 1
      pr[2:length(pr)] = 0
    }

    if(sum(pr == max(pr)) > 1){
      MAXs      = which(pr == max(pr))
      class     = MAXs[which(prior$lambda[MAXs] == max(prior$lambda[MAXs]))]
    }else{
      class     = which(pr == max(pr))
    }
    return(list(pred = class, prob = pr))
  }
  out           = lapply(x, cla, prior = prior)
  report        = list(pred = unlist(lapply(out, function(x){x$pred})),
                    prob = do.call(rbind, lapply(out, function(x){x$prob/sum(x$prob)})))
  if(type == "low_high"){
    report$pred = ifelse(report$pred > 1, "high", "low")
    if(ncol(report$prob) > 2){
      report$prob = cbind(report$prob[,1], rowSums(report$prob[, 2:ncol(report$prob)]))
    }
  }else if(type == "exact"){
    report$pred = ifelse(report$pred == 1, "low", ifelse(report$pred == 2, "high", "extreme"))
  }else{
    stop("Parameter type can only be either low_high or exact.\n")
  }
  return(report)
}



# predict number of mutation #
#' CM
#' @description predict number of mutation
#' @param Data mutSet
#' @param mu mu of the gene
#' @param o_scale o_scale
#' @param MRtriProb MRtriProb
#' @param gids gene IDs
#' @param zero_p zero_p
#' @return A list for all the inputs for Instan
#' @export
#' @examples
#' \dontrun{
#' CM(Data)
#' }
#'
CM = function(Data,  mu,  o_scale, MRtriProb, gids, zero_p = 0){
  gene.mu              = mu
  o_scale              = o_scale
  MRtriProb            = MRtriProb
  subData              = Data$clone()
  subData$mut          = Data$mut[Data$mut$Ensembl_gene_id %in% gids, ]
  sampleN              = as.character(Data$samples$SampleID)

  MR_p               = MRtriProb$MR_p     # relative mutation frequency for each tri-nucleotide context #
  ### get gene.prob
  geneProb           = getGeneBgProb(Data$exomeGene, MR_p, gid = gids)
  # count number of mutation per patient#
  mutPerP            = count_Mut(subData$mut)
  mutPerP$count      = mutPerP$count/get_glen(subData$exomeGene, selGid = names(subData$gid)) * 1000000
  ### get gene length
  geneLen            = get_glen(Data$exomeGene, selGid = gids, byGene= TRUE)



  if("data.table" %in% class(geneProb)){
    geneProb           = data.table(geneProb)
  }
  setkey(geneProb, gid, consequence)


  ## get mutation for non-silent/silent mutation ##
  if(any(zero_p != 0)){
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

  rownames(mut_pre_nonsil) = rownames(mut_pre_sil) = gids
  colnames(mut_pre_nonsil) = colnames(mut_pre_sil) = sampleN

  ## get observed mutation count for non-silent and silent mutation ##
  mut_obs_sil          = count_Mut(subData$silent(), gid = gids, sampleN = sampleN)
  mut_obs_nonsil       = count_Mut(subData$nonsilent(), gid = gids, sampleN = sampleN)

  out                  = melt(mut_obs_sil)
  colnames(out)        = c("gid", "sample","mut_obs_sil")
  out$mut_obs_nonsil   = melt(mut_obs_nonsil)$value
  out$mut_pre_sil      = melt(mut_pre_sil)$value
  out$mut_pre_nonsil   = melt(mut_pre_nonsil)$value
  out$geneLen          = geneLen[as.character(out$gid)]
  out$geneProb0        = geneProb[J(as.character(out$gid),0),]$prob
  out$geneProb1        = geneProb[J(as.character(out$gid),1),]$prob
  return(out)
}
