#log likelihood ratio for the parameters ####
llratio.f <- function(BKcov,Q,R,r=NULL,E=0){
  Dmat = 2*ginv(BKcov)
  dvec = 2*(Q%*% ginv(BKcov))
  if(is.null(r)){solveQP = solve.QP(Dmat, dvec = dvec, t(R), meq = E, factorized = FALSE)}else{
    solveQP = solve.QP(Dmat, dvec = dvec, t(R), r, meq = E, factorized = FALSE)}
  tildeQ = solveQP$solution
  Q = solveQP$unconstrained.solution
  #K = length(Q)
  #logHu=as.numeric( ( -K/2*log(2*pi) )-( 0.5*log(det(BKcov) ) ) )
  #logHm=as.numeric( ( -K/2*log(2*pi) )-( 0.5*log(det(BKcov) ) )-( 0.5* t(Q- tildeQ)%*%ginv(BKcov)%*% (Q-tildeQ)))
  #llratio=(logHu-logHm)
  llratio = (as.numeric(0.5* t(Q- tildeQ)%*%ginv(BKcov)%*% (Q-tildeQ)))
  return(llratio)}

