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

robin <- function(X, Y, E, clin=NULL, iterations=10000, burn.in=NULL, robust=TRUE, sparse=TRUE, structure=c("sparsegroup","group","individual"), hyper=NULL, debugging=FALSE)
{

  structure = match.arg(structure)
  this.call = match.call()
  intercept = TRUE

  if(iterations<1) stop("iterations must be a positive integer.")
  if(is.null(burn.in)){
    BI = floor(iterations/2)
  }else if(burn.in>=1){
    BI = as.integer(burn.in)
  }else{
    stop("burn.in must be a positive integer.")
  }
  if(iterations<=BI) stop("iterations must be larger than burn.in.")


  dat = Data.matrix(X, Y, E, clin, intercept, debugging)

  xx = dat$xx
  y = dat$y
  CLC = dat$CLC
  n = dat$n; s = dat$s
  env = dat$env
  size = dat$size
  G.names = dat$G.names
  E.names = dat$E.names
  clin.names = dat$clin.names

  CLC.names = colnames(CLC)
  nclc = ncol(CLC)

 if(debugging) message("No. of G: ", s, "\tNo. of E: ", env, "\tNo. of G+GxE: ", ncol(xx), "\n")

  clcxx = cbind(CLC, xx)
  lasso.cv = glmnet::cv.glmnet(clcxx,y,alpha=1,nfolds=5)
  lambda.cv = lasso.cv$lambda.min;
  lasso.fit = glmnet::glmnet(clcxx, y, family="gaussian",alpha=1,nlambda=50)
  coeff.array = as.vector(stats::predict(lasso.fit, s=lambda.cv, type="coefficients"))[-2];

  hatAlpha = coeff.array[c(1:nclc)]
  hatBeta = matrix(coeff.array[-c(1:nclc)], ncol = s)

  if(robust){
    out = Robust(xx, y, CLC, s, size, iterations, hatAlpha, hatBeta, sparse, structure, hyper, debugging)
  }else{
    out = NonRobust(xx, y, CLC, s, size, iterations, hatAlpha, hatBeta, sparse, structure, hyper, debugging)
  }

  coeff.main = apply(out$GS.alpha[-(1:BI),,drop=FALSE], 2, median); names(coeff.main) = CLC.names;
  coeff.GE = matrix(apply(out$GS.beta[-(1:BI),], 2, median), size, dimnames=list(c("main",E.names),G.names))

  Int = coeff.main[1]
  coeff.E = tail(coeff.main, env); #names(coeff.E) = E.names;
  coeff.clin = head(coeff.main[-1], -env)
  if(env>0){
    coeff.clin = head(coeff.main[-1], -env)
  }else{
    coeff.clin = coeff.main[-1]
  }

  coefficient = list(Int=Int, clin=coeff.clin, E=coeff.E, GE=coeff.GE)

  fit = list(call = this.call, posterior = out, coefficient=coefficient, burn.in = BI, iterations=iterations, design=list(xx=xx, CLC=CLC))

  # if(debugging && sparse) fit$debugList = out$debugList;
  class(fit)=c("robin", class(out))
  fit
}
