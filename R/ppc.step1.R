ppc.step1 <- function(y.o, model,
                      sample.cov = NULL, sample.mean = NULL, sample.nobs = NULL, group = NULL,
                      n.groups, constraints = "", WLS.V = NULL, NACOV = NULL,
                      nchains=2, nadapt = 1000, nburnin=5000, nsample=5000, dp=NULL,convergence="manual",
                      imp = NULL,n.r,nsim=5000,
                      post, pT, free.i){

  ps1 <- posterior.step1(y.o=y.o, model=model,
                         sample.cov = sample.cov, sample.mean = sample.mean, sample.nobs = sample.nobs, group = group,
                         constraints = constraints, WLS.V = WLS.V, NACOV = NACOV,
                         nchains=nchains, nadapt = nadapt, nburnin=nburnin, nsample=nsample, dp=dp,
                         imp=imp,convergence=convergence)

  ss1 <- sim.step1(n.r=n.r,nsim=nsim,post=ps1$post,pT=ps1$pT,free.i=ps1$free.i,group=group,n.groups=n.groups)

  results <- list(traceplot=ps1$traceplot,pT=ps1$pT,y.s=ss1$y.s)

  return(results)
}
