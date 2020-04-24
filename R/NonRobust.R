
NonRobust <- function(xx, y, CLC, s, L, max.steps, hatAlpha, hatBeta, sparse, structure, hyper, debugging){
  nclc = ncol(CLC)
  n = nrow(xx)

  a = b = r = 1
  hatSigmaSq = 1
  hatInvTauSq = rep(1, s)
  lambdaSq = 1
  hatPi = 0.5
  sh1 = sh0 = 1
  invSigAlpha0 = diag(10^-3, nclc)
  alpha = gamma = 1

  sh0_0 = sh0_1 = sh0
  sh1_0 = sh1_1 = sh1
  hatPi0 = hatPi1 = hatPi
  r1 = r2 = r


  hatGamma = matrix(1, nrow = L, ncol = s)
  hatInvGammaSq = 1/hatGamma

  progress = ifelse(debugging, 10^(floor(log10(max.steps))-1), 0)

  if(sparse){
    fit=switch (structure,
                "sparsegroup" = BSGL_SS(xx, y, CLC, s, L=L, max.steps, hatAlpha, hatBg=hatBeta, hatSigmaSq, hatGamma=hatGamma, invSigAlpha0,
                                        hatSsq=1, hatPi0=hatPi, hatPi1=hatPi, hatT=1, sh0_1, sh0_0, sh1_1, sh1_0, c=1, d=1, 0.05, progress),

                "group" = BGLPointMass(xx, y, CLC, s, L, max.steps, hatAlpha, c(hatBeta), hatInvTauSq, invSigAlpha0, hatPiStar=hatPi,
                                       lambdaSq, hatSigmaSq, a, b, alpha, gamma, sh1, sh0, progress),

                "individual" = BL_SS(xx, y, CLC, max.steps, hatAlpha, c(hatBeta), hatInvTauSq=rep(1,s*L), invSigAlpha0, hatPi,
                                     lambdaSq, hatSigmaSq, a, b, alpha, gamma, sh1, sh0, progress)
    )
  }else{
    fit=switch (structure,
                "sparsegroup" = BSGL(xx, y, CLC, s, L, max.steps, hatAlpha, hatBeta, hatInvTauSq, hatInvGammaSq, invSigAlpha0,
                                     lambdaSq, lambdaSq, hatSigmaSq, s1=1, s2=1, r1, r2, a, b, progress),

                "group" = BGL(xx, y, CLC, s, L, max.steps, c(hatBeta), hatAlpha, hatInvTauSq, invSigAlpha0,
                              lambdaSq, hatSigmaSq, a, b, alpha, gamma, progress),

                "individual" = BLasso(xx, y, CLC, max.steps, c(hatBeta), hatAlpha, hatInvTauSq=rep(1,s*L), invSigAlpha0,
                                      lambdaSq, hatSigmaSq, a, b, alpha, gamma, progress)
    )
  }

  # out = list(posterior = list( GS.alpha = fit$GS.alpha,
  #                              GS.beta = fit$GS.beta))
  out = list( GS.alpha = fit$GS.alpha,
              GS.beta = fit$GS.beta)

  if(sparse){
    class(out)=c("Sparse", "BVS")
  }else{
    class(out)=c("NonSparse", "BVS")
  }

  out

}
