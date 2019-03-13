ppc.plot <- function(llratio.s,llratio.r){

  Mode <- function(x) {
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]}

  h <- hist(llratio.s[-which(llratio.s==Mode(llratio.s))],breaks=20,plot=FALSE)
  if(sum(llratio.s==Mode(llratio.s))/length(llratio.s)>=.10){
    hist(llratio.s[-which(llratio.s==Mode(llratio.s))],ylim=c(0,max(sum(llratio.s==Mode(llratio.s)),h$counts[1])),
         xlim=c(0,max(max(llratio.s),max(llratio.r))),
         xlab=expression(italic(D)),ylab="Frequency",main="",breaks=20)
    segments(x0=Mode(llratio.s),y0=0,x1=Mode(llratio.s),y1=sum(llratio.s==Mode(llratio.s)),col="black",lwd=5)
    abline(v=llratio.r,col="red")
  }else{
    hist(llratio.s,main="",freq=TRUE,breaks=seq(0,max(llratio.s),length.out=40),
         xlim=c(0,max(max(llratio.s),llratio.r)),xlab=expression(italic(D)))
    abline(v=llratio.r,col="red")}
}
