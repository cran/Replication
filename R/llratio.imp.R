llratio.imp <- function(step2step3,imp,model,effectsize=FALSE, s.i,
                        sample.cov = NULL, sample.mean = NULL, sample.nobs = NULL,
                        group = NULL, cluster = NULL, constraints = "", WLS.V = NULL, NACOV = NULL,
                        bayes=FALSE,dp=NULL,nchains=2){


  mat <- step2step3$`H0 matrices`
  R <- mat$R
  r <- mat$r
  E <- mat$E

  llratio.i <- matrix(NA)
  pT <- list()
  m <- imp$m


  if(bayes==TRUE){
    for (i in 1:m){
      print(i)
      data.r.imp <- complete(imp,i)
      data.r.imp <- data.frame(sapply(data.r.imp , function(x) as.numeric(as.character(x))))
      fit_r <- bsem(model, data=data.r.imp,
                    n.chains=nchains, dp=dp,
                    sample.cov = sample.cov, sample.mean = sample.mean, sample.nobs = sample.nobs,
                    group = group, constraints = "", WLS.V = WLS.V, NACOV = NACOV)
      pT[[i]] <- parameterTable(fit_r)
      pT.r <- parameterTable(fit_r)
      free.i <- which(pT.r$free!=0)
      BKcov.r <- lavInspect(fit_r,"vcov")

      if (sum(duplicated(rownames(BKcov.r)))>0){
        pT.r <- pT.r[free.i,][-which(duplicated(rownames(BKcov.r))),]
        Q.r <- pT.r$est
        BKcov.r <- BKcov.r[!duplicated(rownames(BKcov.r)), !duplicated(colnames(BKcov.r))]
      }else{
        Q.r <- pT.r$est[free.i]
      }

      if(effectsize==TRUE){
        s <- vector()
        for (j in 1:length(r)){s[j] <- pT$est[ pT$id == s.i[j] ]}
        r.e <- r*s
      }

      llratio.i[i] <-tryCatch(llratio.f(BKcov=BKcov.r,Q=Q.r,R=R,r=r,E=E),
                              error=function(e) NA)
    }
  }else{

    for(i in 1:m){
      print(i)
      data.r.imp <- complete(imp,i)
      data.r.imp <- data.frame(sapply(data.r.imp , function(x) as.numeric(as.character(x))))
      fit_r <- sem(model, data=data.r.imp,
                   sample.cov = sample.cov, sample.mean = sample.mean, sample.nobs = sample.nobs,
                   group = group, constraints = "", WLS.V = WLS.V, NACOV = NACOV)

      BKcov.r <- lavInspect(fit_r,"vcov")
      pT[[i]] <- parameterTable(fit_r)
      pT.r <- parameterTable(fit_r)
      free.i <- which(pT.r$free!=0)

      if (sum(duplicated(rownames(BKcov.r)))>0){
        pT.r <- pT.r[free.i,][-which(duplicated(rownames(BKcov.r))),]
        Q.r <- pT.r$est
        BKcov.r <- BKcov.r[!duplicated(rownames(BKcov.r)), !duplicated(colnames(BKcov.r))]
      }else{
        Q.r <- pT.r$est[free.i]
      }

      if(effectsize==TRUE){
        s <- vector()
        for (j in 1:length(r)){s[j] <- pT$est[ pT$id == s.i[j] ]}
        r.e <- r*s

        llratio.i[i] <-tryCatch(llratio.f(BKcov=BKcov.r,Q=Q.r,R=R,r=r.e,E=E),
                                  error=function(e) NA)

      }else{

        llratio.i[i] <-tryCatch(llratio.f(BKcov=BKcov.r,Q=Q.r,R=R,r=r,E=E),
                                  error=function(e) NA)}

    }}

  pT.m.est <- rowMeans(sapply(pT,'[[',"est"), na.rm=TRUE)
  pT.m.se <- rowMeans(sapply(pT,'[[',"se"), na.rm=TRUE)
  pT <- pT[[1]]
  pT$est <- pT.m.est; pT$se <- pT.m.se
  pT <- pT[which(pT$free!=0),]
  pT <- pT[!(duplicated(pT$label))|pT$label=="",c(2:4,11,14:15)]

  llratio.p <- function(x){sum(step2step3$llratio.s>=x)/length(step2step3$llratio.s)}
  pvals <- unlist(lapply(llratio.i,llratio.p))

  results <- list(pT=pT,llratio.i=llratio.i,pvals=pvals)

  return(results)

}
