#' fit a robust Bayesian variable selection
#'
#' fit a robust Bayesian variable selection model for G×E interactions.
#'
#' @keywords models
#' @param X the matrix of predictors (genetic factors) without intercept. Each row should be an observation vector. A column of 1 will be added to the X matrix
#' as the intercept.
#' @param Y the response variable. The current version of robin only supports continuous response.
#' @param E a matrix of environmental factors. The interaction terms between X (G factors) and E will be automatically created and included in the model.
#' @param clin a matrix of clinical variables. Clinical variables are not subject to penalty.
#' @param iterations the number of MCMC iterations.
#' @param burn.in the number of iterations for burn-in.
#' @param robust logical flag. If TRUE, robust methods will be used.
#' @param sparse logical flag. If TRUE, spike-and-slab priors will be used to shrink coefficients of irrelevant covariates to zero exactly.
#' @param structure three choices are available. "sparsegroup" for sparse-group selection, which is a bi-level selection on both group-level and individual-level. "group" for selection on group-level only. "individual" for selection on individual-level only.
#' @param hyper a named list of hyperparameters.
#' @param debugging logical flag. If TRUE, progress will be output to the console and extra information will be returned.
#'
#' @details Consider the data model described in "\code{\link{data}}":
#' \deqn{Y_{i} = \alpha_{0} + \sum_{t=1}^{q}\alpha_{t}Clin_{it} + \sum_{m=1}^{k}\theta_{m}E_{im}  + \sum_{j=1}^{p} \big(U_{ij}^\top\beta_{j}\big) +\epsilon_{i},}
#' where the main and interaction effects of the jth genetic variant is corresponding to the coefficient vector \eqn{\beta_{j}=(\beta_{j1}, \beta_{j2},\ldots,\beta_{jL})^\top}.
#'
#' When \emph{structure="sparsegroup"} (default setting), selection will be conducted on both individual and group levels (bi-level selection):
#' \itemize{
#' \item \strong{Group-level selection:} by determining whether \eqn{||\beta_{j}||_{2}=0}, we can know if the jth genetic variant has any effect at all.
#' \item \strong{Individual-level selection:} investigate whether the jth genetic variant has main effect, G\eqn{\times}E interaction or both, by determining which components in \eqn{\beta_{j}} has non-zero values.
#' }
#' If \emph{structure="group"}, only group-level selection will be conducted on \eqn{||\beta_{j}||_{2}}. If \emph{structure="individual"}, only individual-level selection will be conducted on each \eqn{\beta_{jl}}, (\eqn{l=1,\ldots,L}).
#' \cr
#'
#' When \emph{sparse=TRUE} (default), spike--and--slab priors are imposed on individual and/or group levels to identify important main and interaction effects. Otherwise, Laplacian shrinkage will be used.
#'
#' When \emph{robust=TRUE} (default), the distribution of \eqn{\epsilon_{i}} is defined as a Laplace distribution with density
#' \eqn{
#' f(\epsilon_{i}|\nu) = \frac{\nu}{2}\exp\left\{-\nu |\epsilon_{i}|\right\}
#' }, (\eqn{i=1,\dots,n}), which leads to a Bayesian formulation of LAD regression. If \emph{robust=FALSE}, \eqn{\epsilon_{i}} follows a normal distribution.
#'
#' @seealso \code{\link{GxESelection}}
#'
#' @examples
#' data(GxE_small)
#'
#' ## default method
#' iter = 5000
#' fit=robin(X, Y, E, clin, iterations = iter)
#' fit$coefficient
#'
#' ## Ture values of parameters of mian G effects and interactions
#' coeff$GE
#'
#' ## Compute TP and FP
#' sel = GxESelection(fit)
#' pos = which(sel$indicator != 0)
#' tp = length(intersect(which(coeff$GE != 0), pos))
#' fp = length(pos) - tp
#' list(tp=tp, fp=fp)
#'
#' \donttest{
#' ## alternative: robust group selection
#' fit=robin(X, Y, E, clin, iterations = iter, structure="g")
#' fit$coefficient
#'
#' ## alternative: non-robust sparse group selection
#' fit=robin(X, Y, E, clin, iterations = iter, robust=FALSE)
#' fit$coefficient
#' }
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

  coeff.main = apply(out$GS.alpha[-(1:BI),,drop=FALSE], 2, stats::median); names(coeff.main) = CLC.names;
  coeff.GE = matrix(apply(out$GS.beta[-(1:BI),], 2, stats::median), size, dimnames=list(c("main",E.names),G.names))

  Int = coeff.main[1]
  coeff.E = utils::tail(coeff.main, env); #names(coeff.E) = E.names;
  coeff.clin = utils::head(coeff.main[-1], -env)
  if(env>0){
    coeff.clin = utils::head(coeff.main[-1], -env)
  }else{
    coeff.clin = coeff.main[-1]
  }

  coefficient = list(Int=Int, clin=coeff.clin, E=coeff.E, GE=coeff.GE)

  fit = list(call = this.call, posterior = out, coefficient=coefficient, burn.in = BI, iterations=iterations, design=list(xx=xx, CLC=CLC))

  # if(debugging && sparse) fit$debugList = out$debugList;
  class(fit)=c("robin", class(out))
  fit
}
