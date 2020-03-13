Data.matrix <- function(X, Y, E, clin, intercept)
{
  x = as.matrix(X); y = cbind(Y)
  n = nrow(x); s = ncol(x)
  noClin = noE = TRUE
  CLC = NULL
  env = nclc = 0

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


  # CLC.names = colnames(CLC)
  # nclc = ncol(CLC)

  if(is.null(colnames(x))){
    G.names = paste("G", 1:s, sep="")
  }else{
    G.names = colnames(x)
  }

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

  list(xx=xx, y=y, CLC=CLC, n=n, s=s, env=env, size=size)
}
