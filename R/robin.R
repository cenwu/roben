#' fit a robust Bayesian variable selection
#'
#' fit a robust Bayesian variable selection model for G×E interactions.
#'
#' @keywords models
#' @param X the matrix of predictors (genetic factors) without intercept. Each row should be an observation vector. A column of 1 will be added to the X matrix
#' as the intercept.
#' @param Y the response variable. The current version of BVCfit only supports continuous response.
#' @param E a matrix of environmental factors for interactions.
#' @param clin a matrix of clinical variables. Clinical variables are not subject to penalty.
#' @param iterations the number of MCMC iterations.
#' @param burn.in the number of iterations for burn-in.
#' @param sparse logical flag. If TRUE, spike-and-slab priors will be used to shrink coefficients of irrelevant covariates to zero exactly. 'sparse' has effect only when VC=TRUE.
#' @param structural logical flag. If TRUE, the coefficient functions with varying effects and constant effects will be penalized separately. 'structural' has effect only when VC=TRUE.
#' @param hyper a named list of hyperparameters.
#' @param debugging logical flag. If TRUE, progress will be output to the console and extra information will be returned.
#'
#' @export

robin <- function(X, Y, E, clin=NULL, iterations=10000, burn.in=NULL, robust=TRUE, sparse=TRUE, structural=TRUE, hyper=NULL, debugging=FALSE)
{

  # this.call = match.call()
  # intercept = TRUE
  # dat = Data.matrix(X, Y, E, clin, intercept)
  #
  # xx = dat$xx
  # y = dat$y
  # CLC = dat$CLC
  # n = dat$n
  # s = dat$s
  # env = dat$env
  # size = dat$size

  x = as.matrix(X); y = cbind(Y)
  n = nrow(x); s = ncol(x)
  noClin = noE = TRUE
  CLC = NULL
  env = nclc = 0
  intercept = TRUE
  this.call = match.call()

  if(nrow(y) != n)  stop("Length of Y does not match the number of rows of X.");

  if(iterations<1) stop("iterations must be a positive integer.")
  if(is.null(burn.in)){
    BI = floor(iterations)/2
  }else if(burn.in>=1){
    BI = as.integer(burn.in)
  }else{
    stop("burn.in must be a positive integer.")
  }
  if(iterations<=BI) stop("iterations must be larger than burn.in.")


  if(!is.null(clin)){
    clin = as.matrix(clin)
    if(nrow(clin) != n)  stop("clin has a different number of rows than X.");
    if(is.null(colnames(clin))){colnames(clin)=paste("clin.", 1:ncol(clin), sep="")}
    CLC = clin
    noClin = FALSE
    Clin.names = colnames(clin)
  }

  if(intercept){ # add intercept
    CLC = cbind(matrix(1,n,1,dimnames=list(NULL, "IC")), CLC)
  }

  if(!is.null(E)){
    E = as.matrix(E);env = ncol(E)
    if(nrow(E) != n)  stop("E has a different number of rows than X.");
    if(is.null(colnames(E))){colnames(E)=paste("E.", 1:env, sep="")}
    CLC = cbind(CLC, E)
    noE = FALSE
  }else if(!debugging){
    stop("E factors must be provided.")
  }

  # if(any(c(nrow(E), nrow(clin), length(y)) != n)){
  #   stop("Input data (X, Y, E or clin) have different numbers of observations.")
  # }


  CLC.names = colnames(CLC)
  nclc = ncol(CLC)

  if(is.null(colnames(x))){
    G.names = paste("G", 1:s, sep="")
  }else{
    G.names = colnames(x)
  }

  # x = cbind(1, x) # add intercept
  if(!noE){
    size = env+1
    xx = as.data.frame(matrix(0, n, s*(env+1)))
    for(j in 1:s){
      last = j*(env+1); first = last-env
      xx[,first:last] = cbind(x[,j], E*x[,j])
      colnames(xx)[first:last] = c(G.names[j], paste(G.names[j], "E", 1:env, sep=""))
    }
    xx = as.matrix(xx)
  }else{
    xx = x
  }

  if(debugging) message("No. of G: ", s, "\tNo. of E: ", env, "\tNo. of G+GxE: ", ncol(xx), "\n")

  clcxx = cbind(CLC, xx)
  lasso.cv = glmnet::cv.glmnet(clcxx,y,alpha=1,nfolds=5)
  lambda.cv = lasso.cv$lambda.min;
  lasso.fit = glmnet::glmnet(clcxx, y, family="gaussian",alpha=1,nlambda=50)
  coeff.array = as.vector(stats::predict(lasso.fit, s=lambda.cv, type="coefficients"))[-2];

  hatAlpha = coeff.array[c(1:nclc)]
  hatBeta = matrix(coeff.array[-c(1:nclc)], ncol = s)

  # if(VC){
  #   hat.m = coeff.array[1:q]      ## coeff for intercept
  #   hat.r0 = coeff.array[(1:s)+q] ## coeff for constant
  #   hat.r.star = utils::head(coeff.array, ncol(xx))[-(1:(s+q))] ## coeff for varying part
  # }

  # coeff.clc = utils::tail(coeff.array, -ncol(xx)) ## E CLC Z EX ZX
  # hat.clc = coeff.clc[1:nclc]              ## E CLC Z
  # hat.zeta = utils::tail(coeff.clc, -nclc)  ## EX ZX
  #
  # if(!VC){
  #   out = BLasso(xx, y, CLC, EX, ZX, s, iterations, coeff.array[1], coeff.array[2:ncol(xx)], hat.clc, hat.zeta, hyper, debugging)
  #   CC = apply(out$posterior$GS.r0[-c(1:BI),,drop=FALSE], 2, stats::median)
  #   LL = apply(out$posterior$GS.rs[-c(1:BI),,drop=FALSE], 2, stats::median)
  #   names(CC) = names(LL) = Var.names
  #   coeff = list(intercept=stats::median(out$posterior$GS.m[-c(1:BI)]), Z=stats::median(out$posterior$GS.Z[-c(1:BI)]), Main=CC, Interaction=LL)
  # }else{
  #   if(structural){
  #     out = BVC_SI(xx, y, CLC, EX, s, q, iterations, hat.m, hat.r0, hat.r.star, hat.clc, hat.zeta, sparse, hyper, debugging)
  #     INT = apply(out$posterior$GS.m[-c(1:BI),,drop=FALSE], 2, stats::median)
  #     CC = apply(out$posterior$GS.r0[-c(1:BI),,drop=FALSE], 2, stats::median)
  #     VV = apply(out$posterior$GS.rs[-c(1:BI),,drop=FALSE], 2, stats::median)
  #     coeff = cbind(INT, rbind(CC, matrix(VV, nrow = q-1)))
  #   }else{
  #     hat.r = c(rbind(hat.r0, matrix(hat.r.star, nrow = (q-1))))
  #     out = BVC_NS(design$Xns, y, CLC, EX, s, q, iterations, hat.m, hat.r, hat.clc, hat.zeta, sparse, hyper, debugging)
  #     INT = apply(out$posterior$GS.m[-c(1:BI),,drop=FALSE], 2, stats::median)
  #     VV = apply(out$posterior$GS.rs[-c(1:BI),,drop=FALSE], 2, stats::median)
  #     coeff = cbind(INT, matrix(VV, nrow = q))
  #   }
  #   colnames(coeff) = c("intercept",Var.names)
  #   rownames(coeff) = paste("basis", 0:(q-1), sep="")
  # }
  #
  #
  # if(noE && noClin){
  #   coeff.clin = NULL
  #   coeff.E = NULL
  #   coeff.zeta = NULL
  # }else if(!noE && !noClin){
  #   coeff.clc = apply(out$posterior$GS.clc[-c(1:BI),,drop=FALSE], 2, stats::median)
  #   coeff.clin = coeff.clc[-1]
  #   coeff.E = coeff.clc[1]
  #   coeff.zeta = apply(out$posterior$GS.zeta[-c(1:BI),,drop=FALSE], 2, stats::median)
  #   names(coeff.clin) = Clin.names
  #   names(coeff.zeta) = paste("G", 1:s, sep="")
  # }else if(noE){
  #   coeff.clin = apply(out$posterior$GS.clc[-c(1:BI),,drop=FALSE], 2, stats::median)
  #   coeff.E = NULL
  #   coeff.zeta = NULL
  #   names(coeff.clin) = Clin.names
  # }else{
  #   coeff.clin = NULL
  #   coeff.E = apply(out$posterior$GS.clc[-c(1:BI),,drop=FALSE], 2, stats::median)
  #   coeff.zeta = apply(out$posterior$GS.zeta[-c(1:BI),,drop=FALSE], 2, stats::median)
  #   names(coeff.zeta) = paste("G", 1:s, sep="")
  # }
  #
  # coefficient = list(E=coeff.E, clin=coeff.clin, EX=coeff.zeta, ZX=coeff)
  # if(noE){
  #   coefficient$E = NULL
  #   coefficient$EX = NULL
  # }
  # if(noClin) coefficient$clin = NULL
  #
  # fit = list(call = this.call, posterior = out$posterior, coefficient=coefficient, burn.in = BI, iterations=iterations)

  fit = list(call = this.call, burn.in = BI, iterations=iterations, design=xx)

  # if(debugging && sparse) fit$debugList = out$debugList;
  # if(VC) fit$basis = basis;
  #
  # class(fit)=c("BVCfit", class(out))
  fit
}
