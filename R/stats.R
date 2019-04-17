#' get data mu for each gene #####
#' @param y input data
#' @param offset offset
#' @return muhat
#'
get_muhat <- function(y, offset){
  return(apply(y/exp(offset), 1, mean))
}


#' get sigma ###
#' @param mu mu
#' @param mu_glm glm generated mu
#' @return sigma
getsigma <- function(mu, mu_glm){
  return(sd(log(mu)-log(mu_glm)))
}

#' get zero_p ###
#' @param zero_p_zip zero_p_zip
#' @return zero_p
#' @export
getzero_p <- function(zero_p_zip){
  return(exp(zero_p_zip)/(1+exp(zero_p_zip)))
}


#' MLE to estimate genewise dispersion ###
#' @param y input y
#' @param mu mu
#' @param offset offset
#' @param cores number of threads
#' @param ZI A boolen. Default is FALSE which Negative binomial model will be used.
#' If TRUE, zero-inflated model will be used.
#' @return mu
get_genewise_dispersion_mle_mu <- function(y, mu, offset, cores, ZI = FALSE){
  ngene       = nrow(y)
  nsam        = ncol(y)
  muA         = exp(offset)*mu
  tmp         = cbind(y, muA)
  getdispersionMLE = function(data){
    #print(1)
    n         = length(data)
    y         = data[1:(n/2)]
    mu        = data[(n/2+1):n]
    ind0      = y == 0
    obj       = function(x){
      if(ZI){
        zero_p = x
        # zero inflated
        return(-sum( log(zero_p + (1-zero_p) * exp(-mu[ind0]))) -
          sum(log(1-zero_p) -mu[!ind0] + y[!ind0] * log(mu[!ind0]) - lgamma(y[!ind0] + 1)))
      }else{
        # regular NB
        phi    = x
        return( -sum( lgamma(y + 1/phi) - lgamma(1/phi) -lgamma(y+1) - 1/phi*log(1+mu*phi) + y*( log(mu) - log(1/phi+mu) ) ))  ###optimize the probability mass function)
      }
    }
    if(ZI){
      return(optimize(obj, interval=c(0, 1))$minimum)
    }else{
      return(optimize(obj, interval=c(10^(-20), 10^20))$minimum)
    }
  }
  x           = lapply(apply(tmp, 1, FUN=list), unlist)
  dispersion  = unlist(mclapply(x, getdispersionMLE ,mc.cores=cores))
  #dispersion <- apply(tmp, 1, getdispersionMLE)
  return(dispersion)
}


#' getmu_post_optimize
#' @description  get optimized mu ########
# mu1 = mu_hat
# mu2 = mu_hat_glm
#' @param y The input data
#' @param offset offset of input data
#' @param mu1 mu_hat
#' @param mu2 mu_hat_glm
#' @param sigma dispersion
#' @param span The extention
#' @param phi dispersion parameter when negative binomial model is used. Otherwise,
#' it will be 0.
#' @param cores number of thread
#' @param zero_p probability of zero portion if zero-inflated model is used. Otherwise,
#' it will be 0
#' @return A vector mutation count for 96 trinucleotide context
#' @examples
#' \dontrun{
#' getmu_post_optimize(mut, gid)
#' }
#'
getmu_post_optimize = function(y, offset, mu1, mu2, phi, sigma, span,cores, zero_p = 0){
  ngene <- nrow(y)
  nsam <- ncol(y)
  if (length(phi)==1){
    phi <- rep(phi, ngene)
  }

  if (length(zero_p)==1){
    zero_p <- rep(zero_p, ngene)
  }
  tmp <- cbind(y, offset, mu1, mu2, phi, zero_p)
  getmu <- function(data, sigma, span){
    #print(1)
    n <- length(data)
    zero_p <- data[n]
    disp <- data[n-1]
    mu2 <- data[n-2]
    mu1 <- data[n-3]
    y <- data[1:((n-4)/2)]
    ind0 = y == 0
    offset <- exp(data[((n-4)/2+1):(n-4)])
    obj <- function(mu){
      if(zero_p != 0){
        # zero-inflated
        if(disp == 0){
          # zero-inflated Poisson
          return( (log(mu)-log(mu2))^2/(2 * sigma^2) - ## prior section
                    sum(log(zero_p + ( 1 - zero_p) * exp((-mu * offset)[ind0]))) - # zero section disp.g is zero_p
                    sum(log(1-zero_p) - (mu * offset)[!ind0]  + y[!ind0]*log((mu * offset)[!ind0]) -lgamma(y[!ind0]+1)))

        }
      }else{
        if(disp == 0){
          ## regular Poisson
          return((log(mu)-log(mu2))^2/(2 * sigma^2) - sum( y*log(mu * offset)-offset*mu -lgamma(y+1)) )
        }else{
          ## regular NB
          return ((log(mu)-log(mu2))^2/(2 * sigma^2) -sum( lgamma(y + 1/disp) - lgamma(1/disp) -lgamma(y+1) -
                                                             1/disp*log(1+mu*disp*offset) + y*( log(mu*offset) - log(1/disp+mu*offset) ) ) )
        }
      }
    }
    return(optimize(obj, interval=c(10^(-10), max(mu1, mu2)*span))$minimum)
  }
  x <- lapply(apply(tmp, 1, FUN=list), unlist)
  mu <- unlist(mclapply(x, function(x) getmu(x, sigma,span) ,mc.cores=cores))
  names(mu) <- names(mu2)
  #mu <- apply(tmp, 1, getmu, sigma, span)
  return(mu)
}

#' get_mu_hat_mle
#' @description  from phi estimate mu
#' @param y input y
#' @param disp disp
#' @param offset offset
#' @param ZI A boolen. Default is FALSE which Negative binomial model will be used.
#' If TRUE, zero-inflated model will be used.
#' @param cores number of threads
#' @return mu
get_mu_hat_mle <- function(y, disp, offset,cores, ZI = FALSE){
  ngene <- nrow(y)
  nsam <- ncol(y)
  tmp <- cbind(y,offset,disp)

  getmuMLE <- function(data){
    n <- length(data)-1
    y <- data[1:(n/2)]
    ind0 = y == 0
    offset <- data[(n/2+1):n]
    disp.g = data[(n+1)]

    obj <- function(x){
      mu <- exp(offset)*x
      if(ZI){
        -sum(log(disp.g + ( 1 - disp.g) * exp(-mu[ind0]))) - # zero section disp.g is zero_p
          sum((log(1-disp.g) - mu[!ind0] ) + y[!ind0]*log(mu[!ind0] ) -lgamma(y[!ind0]+1))
      }else{
        -sum( lgamma(y + 1/disp.g) - lgamma(1/disp.g) -lgamma(y+1) - 1/disp.g*log(1+mu*disp.g) + y*( log(mu) - log(1/disp.g+mu) ) )  ###optimize the probability mass function
      }
    }
    return(optimize(obj, interval=c(10^(-20), 10^20))$minimum)
  }

  x <- lapply(apply(tmp, 1, FUN=list), unlist)
  mu <- unlist(mclapply(x, getmuMLE ,mc.cores=cores))

  return(mu)

}

####### apply NB for mu & phi #####
#' get_mu_phi
#' @description apply NB for mu & phi
#' @param y input y
#' @param offset offset
#' @return mu phi
get_mu_phi = function(y, offset){
  nsam <- ncol(y)
  tmp = cbind(y,offset)
  get_nb <- function(tmp){
    data = data.frame(y=tmp[1:nsam], exp=rep(0,nsam), offset=tmp[(nsam+1):length(tmp)])
    fit = glm.nb(y~exp+offset(offset),data=data)
    return(c(exp(fit$coefficients[1]), 1/fit$theta))
  }
  mu_phi <- apply(tmp,1,get_nb)
  return(mu_phi)
}


#' Calculate mu and zero_p by MLE for ZIP model #
#' @param y observed mutation
#' @param offset offset.
#' @param cores number of cores
#' @importFrom parallel mclapply
#' @return a List of parameters
#' @export
#' @examples
#' \dontrun{
#' get_zip_pairparam_mle(Data, offset)
#' }
#'
get_zip_pairparam_mle = function(y, offset,cores){
  data    = y/offset
  data    = lapply(apply(data, 1, FUN=list), unlist)
  pairparams = mclapply(data,  function(x){ zipMLE(x)$param }, mc.cores = cores)
  return(pairparams)
}

#' Calculate mu and zero_p by MLE for ZIP model #
#' @param x observed mutation
#' @param offset offset.
#' @param tol tol of convergency
#' @return a List of parameters
#' @export
#' @examples
#' \dontrun{
#' zipMLE(Data, offset)
#' }
#'
zipMLE = function (x, offset = 1, tol = 1e-09)
{
  if(all(x == 0)){
    param <- c(0, 1)
    names(param) <- c("lambda", "pi")
    list(iters = 0, loglik = 0, param = param)
  }else{
    x  < - x/offset
    no <- sum(x == 0)
    n <- length(x)
    prop <- no/n
    n1 <- n - no
    x1 <- x[x > 0]
    sx <- sum(x1)
    m <- sx/n
    s <- (sum(x1^2) - m * sx)/(n - 1)
    l1 <- s/m + m - 1
    fx <- m - m * exp(-l1) - l1 + prop * l1
    der <- m * exp(-l1) - 1 + prop
    l2 <- l1 - fx/der
    i <- 2
    while (abs(l2 - l1) > tol) {
      i <- i + 1
      l1 <- l2
      fx <- m - m * exp(-l1) - l1 + prop * l1
      der <- m * exp(-l1) - 1 + prop
      l2 <- l1 - fx/der
      if(i > 10000){
        cat("zipMLE not converge \n.")
        break
      }
    }
    if(i > 10000){
      p <- prop
      l2 <- m
      loglik <- no * log(p + (1 - p) * exp(-l2)) + n1 * log(1 -
                                                              p) + sum( x*log(l2)- l2 -lgamma(x+1))
    }else{
      p <- 1 - m/l2
    }
    loglik <- no * log(p + (1 - p) * exp(-l2)) + n1 * log(1 -
                                                            p) + sum( x*log(l2)- l2 -lgamma(x+1))
    param <- c(l2, p)
    names(param) <- c("lambda", "pi")
    list(iters = i, loglik = loglik, param = param)
  }
}

#' Calculate the relative mutation frequency for each tri-nucleotide context #
#' @param Data mutSet object
#' @param fraction fraction of gene used for nonsilent mutation.
#' @return A vector mutation count for 96 trinucleotide context
#' @export
#' @examples
#' \dontrun{
#' get_bg_MRtri(Data)
#' }
#'
getBgMRtri = function(Data, fraction = 0.6){

  if(class(Data$gid) == "list"){
    cat("\tSamples' bed file are different. Common list of gene will be used to predict sampling bias 'r'.\n")
    gid              = Reduce(intersect, gid)
    gid_nonsil_p     = gid [ gid %in% Data$get_nonsil_passengers(fraction)]
    ## need to be tested ###
  }else{
    gid              = Data$gid
    gid_nonsil_p     = Data$get_nonsil_passengers(fraction)
  }

  sil_all            = count_Mut(Data$exomeGene[Data$exomeGene$consequence == 0,], gid = gid, byContext = TRUE)
  nonsil_p           = count_Mut(Data$exomeGene[Data$exomeGene$consequence != 0,], gid = gid_nonsil_p, byContext = TRUE)
  sil_mut            = count_Mut(Data$silent(), gid = gid, byContext = TRUE)
  nonsil_mut_p       = count_Mut(Data$nonsilent(), gid = gid_nonsil_p, byContext = TRUE)

  valid_ind          = which(sil_all >= 1000 & sil_mut > 0)       # select tri-nucleotide context which has enough incidence
  r                  = mean( (nonsil_mut_p[valid_ind]/ nonsil_p[valid_ind]) * (sil_all[valid_ind] /sil_mut[valid_ind]))      ######### relative to silent mutations
  p_all              = (nonsil_mut_p + sil_mut ) / ( sil_all + r * nonsil_p )  # part of phat
  count_all          = (nonsil_mut_p + sil_mut )

  ref                =  p_all[1]  # reference cate is type 1
  p                  = rep(0,96)
  for (m in 1:96) {
    p[m]             = p_all[m]/ref
  }

  ## calculate p_{indel}, which will be used later for calculation of p_{frameshift} and p_{inframe}
  mut                = Data$nonsilent()
  mut                = mut[mut$Ensembl_gene_id %in% gid_nonsil_p,]
  if(sum(  mut[,"Variant_Type"]=="In_frame" |  mut[,"Variant_Type"]=="Frame_shift"  )>0){
    mutab_indel      = mut[ mut$Variant_Type=="In_frame" | mut$Variant_Type=="Frame_shift",]
    pos              = get_glen(Data$exomeGene, selGid = gid_nonsil_p)
    ind              = nrow(mutab_indel)/pos/r/ref  #### p_{indel}
    indel_c          = nrow(mutab_indel)
  }else{
    ind              = 0
    indel_c          = 0
  }
  p                  = c(p,ind)
  names(p)           = c(1:96,"indel")
  count_all          = c(count_all, indel_c)
  names(count_all)   = c(1:96,"indel")
  return(list(MR_p = p, r = r, ref = ref, count = count_all))
}


# calcuate number of mutation per patient  #
#' CalTMB
#' @description calcuate number of mutation per patient
#' @param x mutSet object
#' @param type Type of mutation to be included. It can be 'all' 'nonsil' 'sil' 'indel' 'snp' 'frameshit'
#' @param sampleN Default is NULL. If it is specified, only subset of samples will be reported.
#' @return A data.frame contain the mutation count
#' @export
#' @examples
#' \dontrun{
#' CalTMB(Data)
#' }
#'
CalTMB = function(x, type = "all", sampleN = NULL){
  if(is.null(sampleN)){
    sampleN   = as.character(x$samples$SampleID)
  }
  if(type == "all"){
    out         = count_Mut(x$mut)
  }else if(type == "sil" ){
    out         = count_Mut(x$silent())
  }else if(type == "nonsil"){
    out         = count_Mut(x$nonsilent())
  }else if(type == "indel"){
    mut         = x$nonsilent()
    out         = count_Mut(mut[mut$Variant_Classification %in% c("Frame_Shift_Del", "Frame_Shift_Ins"),])
  }else if(type == "snp"){
    mut         = x$nonsilent()
    out         = count_Mut(mut[mut$Variant_Type %in% c("SNP"),])
  }else if(type == "frameshift"){
    mut         = x$nonsilent()
    out         = count_Mut(mut[mut$Variant_Type %in% c("Frame_shift"),])
  }else{
    stop("type must be either all, sil or nonsil")
  }
  out$count      = out$count/get_glen(x$exomeGene, selGid = names(x$gid)) * 1000000
  out            = out[sampleN, ]
  return(out)
}



# fit background model with negative-binomial regression #
#' fit_model
#' @description fit background model with negative-binomial regression
#' @param Data mutSet object
#' @param MRtriProb output from get_bg_MRtri
#' @param method Either NB or Poisson
#' @param cores Number of parallel core to be used
#' @param bs.type Either 'all' or nonsil. When all is specificed, all mutation per Mb
#' will be used as sample specific background mutation (TMB). When nonsil is specified,
#' non-synonymous mutation per Mb will be used as sample specific background mutation (TMB).
#' @param mut.nonsil A boolen. If TRUE, non-synonymous mutation will be used for background mutation
#' modeling during training.
#' @param nonsil.fraction Default is 0.6. The percentage of genes whose non-synonymous mutation
#' will be used.
#' @importFrom MASS glm.nb
#' @importFrom parallel mclapply
#' @importFrom data.table data.table setkey
#' @return A vector mutation count for 96 trinucleotide context
#' @export
#' @examples
#' \dontrun{
#' fit_model2(Data)
#' }
#'
fit_model = function(Data, MRtriProb, method = "NB", cores = 1,  bs.type = "nonsil", mut.nonsil = FALSE, nonsil.fraction = 0.6){
  MR_p               = MRtriProb$MR_p     # relative mutation frequency for each tri-nucleotide context #
  gene_bg_p          = getGeneBgProb(Data$exomeGene, MR_p)
  sampleN            = Data$samples$SampleID
  geneLen            = get_glen(Data$exomeGene, selGid = Data$gid, byGene = TRUE)


  # count number of mutation per patient for offset calculation
  if(bs.type == "all"){
    mutPerP          = CalTMB(trainset, sampleN = as.character(trainset$samples$SampleID))
  }else if(bs.type == "nonsil"){
    mutPerP          = CalTMB(trainset, sampleN = as.character(trainset$samples$SampleID), type = "nonsil")
  }
  # for all genes with silent mutations
  gid                = Data$gid
  if(method %in% c("NB", "Poisson")){  ## original method
    selGene          = gid[ gid %in%  Data$silent()[,"Ensembl_gene_id"]]
  }else{
    selGene          = gid           ## new method
  }

  ## calculate offset and observed mutation
  if(mut.nonsil){
    gid_nonsil_p     = Data$get_nonsil_passengers(nonsil.fraction)
    offset_nonsil    = as.matrix(as.numeric(gene_bg_p[J(selGene,1),]$prob + geneLen[selGene]*MR_p["indel"])) %*% t(as.matrix(as.numeric(mutPerP[sampleN,2])))
    offset_nonsil[!selGene %in% gid_nonsil_p,] = 0         ## drivers nonsilent offset will be 0
    offset_sil       = as.matrix(as.numeric(gene_bg_p[J(selGene,0),]$prob)) %*% t(as.matrix(as.numeric(mutPerP[sampleN,2])))
    offset           = log(offset_nonsil + offset_sil)
    sil_mut_matrix   = unname(count_Mut(Data$silent(), selGene, sampleN))
    nonsil_mut_matrix= unname(count_Mut(Data$nonsilent(), selGene, sampleN))
    nonsil_mut_matrix[!selGene %in% gid_nonsil_p, ] = 0   ## drivers nonsilent mutation count will be 0
    y                = nonsil_mut_matrix + sil_mut_matrix
    if(any(is.na(offset)| is.infinite(offset))) offset[is.na(offset) | is.infinite(offset)] = min(offset[!is.infinite(offset)], na.rm = T) - log(2)   ## some gene don't have silent/nonsilent mutation. Set it as 2 fold lower than min
  }else{
    offset             = log(as.matrix(as.numeric(gene_bg_p[J(selGene,0),]$prob)) %*% t(as.matrix(as.numeric(mutPerP[sampleN,2]))))
    if(any(is.na(offset)| is.infinite(offset))) offset[is.na(offset) | is.infinite(offset)] = min(offset[!is.infinite(offset)], na.rm = T) - log(2)    ## some gene don't have silent mutation.
    sil_mut_matrix     = unname(count_Mut(Data$silent(), selGene, sampleN))
    y                 = sil_mut_matrix
  }

  # process covar
  covar              = as.matrix(Data$covar,stringsAsFactors=FALSE)
  which              = rownames(covar)[!is.na(covar[,"expr"]) & !is.na(covar[,"reptime"]) & !is.na(covar[,"hic"])]
  covar[is.na(covar[,"expr"]),"expr"] = mean(as.numeric(covar[,"expr"]),na.rm=T)
  covar[is.na(covar[,"reptime"]),"reptime"] = mean(as.numeric(covar[,"reptime"]),na.rm=T)
  covar[is.na(covar[,"hic"]),"hic"] = mean(as.numeric(covar[,"hic"]),na.rm=T)
  covar_m           = prcomp(apply(covar[,c("expr", "reptime", "hic")],2,as.numeric),scale=T)
  covar[,3:5]       = covar_m$x
  m                 = data.frame(epc = as.numeric(gene_bg_p[J(selGene,0),]$prob),
                                 exp = as.numeric(covar[selGene,"expr"]),
                                 rep = as.numeric(covar[selGene,"reptime"]),
                                 hic = as.numeric(covar[selGene,"hic"]),
                                 or = factor(covar[selGene,"or"], levels =1:2))
  if(length(unique(m$or)) == 1){
    cat("There are no olfactory gene. Olfactory variable will be removed in the model.\n")
    design            = model.matrix(~ exp + rep + hic, m) ### design
  }else{
    design            = model.matrix(~ exp + rep + hic + or, m) ### design
  }
  # fit model
  o_scale           = max(offset)
  offsetS           = offset-o_scale ### rescale 0~1
  ysum              = apply(y, 1, sum)
  offset_sum        = rowSums(exp(offsetS))
  cat(sprintf("Gene total: %s; %s genes have 0 mutation detected synonymous/non-synonymous in training cohort.\n", length(ysum), sum(ysum == 0)))

  if(length(unique(m$or)) == 1){
    dataA             = data.frame(ysum=ysum, exp=design[,"exp"], rep=design[,"rep"], hic=design[,"hic"], offset_sum=offset_sum)
  }else{
    dataA             = data.frame(ysum=ysum, exp=design[,"exp"], rep=design[,"rep"], hic=design[,"hic"], or2=design[,"or2"], offset_sum=offset_sum)
  }
  rownames(dataA)   = selGene
  #  selGene           = selGene[selGene %in% which]     ## update selected gene to gene contain covar info
  dataA_sel         = dataA[selGene[selGene %in% which], ]

  ######### (strategy 1 : negative binomial) ### get all_betas from negative binomial model, needs optimization .....
  if (method == "NB"){
    cat(sprintf("Method: %s \n", method))
    if(length(unique(m$or)) == 1){
      fit             = glm.nb(ysum ~ exp + rep + hic + offset(log(offset_sum)), data = dataA_sel,control=glm.control(maxit=100))
    }else{
      fit             = glm.nb(ysum ~ exp + rep + hic + or2 + offset(log(offset_sum)), data = dataA_sel,control=glm.control(maxit=100))
    }
    ###get mu_hat_glm
    mu_glm          = as.numeric(exp(design%*%as.matrix(fit$coefficients)))
    phi_glm         = fit$theta

    ### NB: joint estimation of mu_hat and phi_hat
    cat("NB: joint estimation of mu_hat and phi_hat\n")
    mu_hat          = get_muhat(y=y, offset=offsetS)
    phi_hat_mle     = get_genewise_dispersion_mle_mu(y=y, mu=mu_hat, offset=offsetS, cores=cores)

    mu_hat_pre      = rep(0,length(mu_hat))
    while(sum(abs(mu_hat-mu_hat_pre)) > 1){
      ###get phi, dispersion
      mu_hat_pre    = mu_hat
      mu_hat        = get_mu_hat_mle(y=y, disp=phi_hat_mle, offset=offsetS, cores=cores)
      phi_hat_mle   = get_genewise_dispersion_mle_mu(y=y, mu=mu_hat, offset=offsetS, cores=cores)
      print(sum(abs(mu_hat-mu_hat_pre)))
    }

    ###get sigma
    cat("Estimate sigma and mu_post.\n")
    sigma          = getsigma(mu=mu_hat, mu_glm=mu_glm)

    mu_post        = getmu_post_optimize(y=y, offset=offsetS, mu1=mu_hat, mu2=mu_glm, phi=phi_hat_mle, sigma = sigma, span=100, cores=cores)
    mu_post_pre    =  rep(0, length(mu_hat))
    while(sum(abs(mu_post-mu_post_pre)) > 1){
      mu_post_pre  = mu_post
      phi_hat_mle  = get_genewise_dispersion_mle_mu(y=y, mu=mu_post, offset=offsetS, cores=cores)
      mu_post      = getmu_post_optimize(y=y, offset=offsetS, mu1=mu_hat, mu2=mu_glm, phi=phi_hat_mle, sigma = sigma, span=100, cores=cores)
      print(sum(abs(mu_post-mu_post_pre)))
    }
    phi_hat        = phi_hat_mle
  }


  ########## (strategy 2: poisson with log-linear) #### get all_betas from Poisson model
  if (method == "Poisson"){
    cat(sprintf("Method: %s \n", method))

    if(length(unique(m$or)) == 1){
      fit           = glm(ysum ~ exp + rep + hic  + offset(log(offset_sum)), data = dataA_sel,control=glm.control(maxit=100),family="poisson")
    }else{
      fit           = glm(ysum ~ exp + rep + hic + or2 + offset(log(offset_sum)), data = dataA_sel,control=glm.control(maxit=100),family="poisson")
    }

    ###get mu_hat_glm
    mu_glm          = as.numeric(exp(design%*%as.matrix(fit$coefficients)))

    ### mu_hat for Poisson regression
    mu_hat          = rowSums(y)/rowSums(exp(offsetS))

    ### get sigma
    sigma           = getsigma(mu=mu_hat, mu_glm=mu_glm)

    ### get optimized mu
    mu_post        = getmu_post_optimize(y=y, offset=offsetS, mu1=mu_hat, mu2=mu_glm, phi=0, sigma = sigma, span=100, cores=cores)
    phi_hat        = rep(0, length(selGene))
  }

  ########## (strategy 3: zero-inflated poisson with log-linear) #### get all_betas from Poisson model
  if (method == "ZIP"){
    cat(sprintf("Method: %s \n", method))
    if(length(unique(m$or)) == 1){
      fit           = zeroinfl(ysum ~ exp + rep + hic, offset = log(offset_sum), data = dataA_sel, control=zeroinfl.control(maxit=1000),dist="poisson")
    }else{
      fit           = zeroinfl(ysum ~ exp + rep + hic + or2, offset = log(offset_sum), data = dataA_sel, control=zeroinfl.control(maxit=1000),dist="poisson")
    }
    mu_zip          = as.numeric(exp(design%*%as.matrix(fit$coefficients$count)))
    zero_p_zip      = getzero_p(as.numeric(design%*%as.matrix(fit$coefficients$zero)))
    names(zero_p_zip) = names(mu_zip) =rownames(design)

    ## get mu_hat and zero_p_hat
    cat("ZIP: joint estimation of mu_hat and zero_p_hat\n")
    pairparams     = get_zip_pairparam_mle(y, offset = exp(offsetS), cores = cores)
    mu_hat         = unlist(lapply(pairparams, function(x){x['lambda']}))
    zero_p_hat     = unlist(lapply(pairparams, function(x){x["pi"]}))
    names(mu_hat)  = sub(".lambda", "", names(mu_hat))
    names(zero_p_hat)  = sub(".pi", "", names(zero_p_hat))

    ### post mu ###################
    ###get sigma
    cat("Estimate sigma and mu_post.\n")
    sigma          = getsigma(mu=mu_hat[mu_hat > 0], mu_glm=mu_zip[mu_hat > 0])
    mu_post        = getmu_post_optimize(y=y, offset=offsetS, mu1=mu_hat, mu2=mu_zip, phi=0, sigma = sigma, span=1000, cores=cores, zero_p = zero_p_zip)
  }



  if(method != "Poisson" & method != "NB" & method != "ZIP") stop("method should be either Poisson, NB, ZIP or ZINB.\n")


  ### gene factor calculation ###
  if(!method %in% c("ZIP", "ZINB")){
    if(length(unique(m$or)) == 1){
      base           = exp(fit$coefficients[1] + fit$coefficients[2] * as.numeric(covar[gid,3]) +
                             fit$coefficients[3] * as.numeric(covar[gid,4]) +
                             fit$coefficients[4] * as.numeric(covar[gid,5]))
    }else{
      base           = exp(fit$coefficients[1] + fit$coefficients[2] * as.numeric(covar[gid,3]) +
                             fit$coefficients[3] * as.numeric(covar[gid,4]) +
                             fit$coefficients[4] * as.numeric(covar[gid,5]) +
                             fit$coefficients[5] * (as.numeric(as.character(covar[gid,6])) - 1))
    }

    names(base)    = gid

    ### assign correction factor based on base ############
    missGene       = gid[gid %in% selGene == FALSE]
    or             = order(base[selGene])
    selAdj         = data.table(data.frame(disp = phi_hat[or], p.base = (base[selGene])[or]))
    selAdj$bin     = ceiling(1:nrow(selAdj)/50)
    dispAvg        = selAdj[, {dispAvg = median(disp); list(dispAvg=dispAvg)}, by='bin']
    setkey(dispAvg,bin)

    missData       = base[missGene]
    missData[is.na(missData)] = 0
    range          = findInterval(missData, selAdj$p.base)
    range[range < 1] =1
    range[range > nrow(selAdj)] = nrow(selAdj)
    missDisp       = dispAvg[J(selAdj[range,]$bin),]$dispAvg

    ### adjustment: dispersion, beta ###
    gene.phi       = rep(0,length(gid))
    names(gene.phi) = gid
    gene.mu        = gene.phi

    gene.phi[selGene] = phi_hat
    gene.phi[missGene] = missDisp

    gene.mu[selGene] = mu_post

    offset         = log(as.matrix(as.numeric(gene_bg_p[J(missGene,0),]$prob)) %*% t(as.matrix(as.numeric(mutPerP[sampleN,2]))))
    offsetS        = offset - o_scale ### rescale
    y              = unname(count_Mut(Data$silent(), missGene, sampleN))
    mu_post_miss   = getmu_post_optimize(y=y, offset=offsetS, mu1=rep(0,length(missGene)), mu2=base[missGene], phi=missDisp, sigma = sigma, span=100, cores=cores)
    #gene.mu[miss.gene] = base[miss.gene]
    mu_post_miss[which(is.na(rowSums(offsetS)))] = base[missGene[which(is.na(rowSums(offsetS)))]]
    gene.mu[missGene] = mu_post_miss
    gene.zeroP     = rep(0, length(gene.mu))
    names(gene.zeroP)= names(gene.mu)
  }else{
    gene.mu        = mu_post
    gene.phi       = rep(0, length(gene.mu))
    names(gene.phi)= names(gene.mu)
    gene.zeroP     = zero_p_zip
  }
  return(list( geneMu = gene.mu, geneProb = gene_bg_p,
               mutPerP = mutPerP , geneDisp = gene.phi,
               o_scale = o_scale, MRtriProb = MRtriProb,
               geneLen = geneLen, geneZero_p = gene.zeroP,
               trainData = dataA_sel,
               otherparams = list(mu_glm = ifelse(method == "ZIP", mu_zip, mu_glm),
                                  disp_glm = ifelse(method == "ZIP", zero_p_zip,
                                                    ifelse(method == "NB", rep(phi_glm, length(gene.mu)) ,rep(0, length(gene.mu)))))
               ))
}




