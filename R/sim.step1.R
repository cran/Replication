sim.step1 <- function(n.r,nsim,post,pT,free.i,group=NULL,n.groups){

  #simulate new data
  fix.i <- which(pT$free==0&pT$op!=":=")
  est.fix<-pT$est[fix.i]

  #define syntax
  lhs <- pT$lhs                 #left hand side
  op <- pT$op                   #operator
  op[which(op=="~1")] <- "~"
  rhs <- pT$rhs                 #right hand side
  rhs[which(rhs=="")] <- "1"
  #separate fixed and free syntax
  lhs.free <- lhs[free.i]; op.free <- op[free.i]; rhs.free <- rhs[free.i]
  lhs.fix <- lhs[fix.i]  ; op.fix <- op[fix.i]  ; rhs.fix <- rhs[fix.i]

  y.s <- list()
  it=0
  Popmodel.line <- list()         #population model part with pop-values from posterior
  Popmodel.line.fix <- list()     #fixed population model part
  post.sample <- round(runif(nsim,1,length(post[[1]])),0)
  print("Creating y.s",quote=FALSE)
  pb = txtProgressBar(min = 0, max = length(post.sample), initial = 0,char = "x",style=3, width = 50)
  for (i in post.sample){
    it=it+1
    setTxtProgressBar(pb,it)
    if(is.null(group)==TRUE){
      for (j in 1:(length(free.i))){
        #left hand side, operator, posterior sample, *, right hand side
        Popmodel.line[j] <- paste(lhs.free[j],op.free[j],post[[j]][i],"*",rhs.free[j])}
      if(length(fix.i)!=0){
        for (j in 1:(length(fix.i))){
          #left hand side, operator, fixed value, *, right hand side
          Popmodel.line.fix[j] <- paste(lhs.fix[j],op.fix[j],est.fix[j],"*",rhs.fix[j])}

        #combine free and fixed lines
        Popmodel.lines <- do.call("paste",c(Popmodel.line,sep="\n "))
        Popmodel.lines.fix <- do.call("paste",c(Popmodel.line.fix,sep="\n "))
        Popmodel <- paste(Popmodel.lines.fix,"\n ",Popmodel.lines)
      }else{
        Popmodel.lines <- do.call("paste",c(Popmodel.line,sep="\n "))
        Popmodel.lines.fix <- do.call("paste",c(Popmodel.line.fix,sep="\n "))
        Popmodel <- paste(Popmodel.lines.fix,"\n ",Popmodel.lines)}
    }

    if(is.null(group)==FALSE){
      unique <- length(free.i)/n.groups
      for (j in 1:unique){
        Popmodel.line[j] <- paste(lhs.free[j],op.free[j],
                                  "c(",post[[j]][i],",",post[[j+unique]][i],")",
                                  "*",rhs.free[j])}
      if(length(fix.i)!=0){
        for (j in 1:(length(fix.i)/n.groups)){
          #left hand side, operator, fixed value, *, right hand side
          Popmodel.line.fix[j] <- paste(lhs.fix[j],op.fix[j],est.fix[j],"*",rhs.fix[j])}

        #combine free and fixed lines
        Popmodel.lines <- do.call("paste",c(Popmodel.line,sep="\n "))
        Popmodel.lines.fix <- do.call("paste",c(Popmodel.line.fix,sep="\n "))
        Popmodel <- paste(Popmodel.lines.fix,"\n ",Popmodel.lines)
      }else{
        Popmodel.lines <- do.call("paste",c(Popmodel.line,sep="\n "))
        Popmodel.lines.fix <- do.call("paste",c(Popmodel.line.fix,sep="\n "))
        Popmodel <- paste(Popmodel.lines.fix,"\n ",Popmodel.lines)}
    }

    y.s[[it]] <- simulateData(model=Popmodel,sample.nobs=n.r) #create new datasets
  }

  cat("\n")
  results <- list(y.s=y.s)
  return(results)

}
