#'Fobj_ALLOC_SAMPLE
#'@title Objective function of the optimal allocation in stratified sampling
#'@description Calculates sum of sample sizes associated  with all strata
#'@references 1. Brito, J. A. M.; Silva, P. L. N. ; Semaan, G. S. ; Maculan, N.(2015). Integer Programming Formulations Applied to Optimal Allocation in Stratified Sampling.
#'Survey Methodology, v. 41, p. 427-442.
#'@examples
#'Nh<-c(212,84,61)
#'Sh2<-c(723.1,2693.4,36231.7)
#'Y<-80548
#'H<-3
#'S<-list(Nh=Nh,Y=Y,cvt=0.05)
#'nh<-c(150,42,30)
#'n<-Fobj_ALLOC_SAMPLE(Sh2,nh,S)
#'@param  Sh2 Variance by stratum
#'@param  nh Number of units in the sample in the hth stratum
#'@param  S  List with Nh, Y and cvt(cv target) (is necessary to define Exa2=S)
#'@return \item{n}{= sum of nh}
#'@import stats
#'@import utils
#'@export

Fobj_ALLOC_SAMPLE<-function(Sh2,nh,S)
{
  Nh<-S$Nh
  Y<-S$Y
  cvt<-S$cvt
  cv<-sqrt(sum(Nh^2*(1/nh-1/Nh)*Sh2))/Y
  n<-sum(nh)
  P<-cv/cvt
  if(P>1) {n<-n+2000^P}
  return(n)
}
