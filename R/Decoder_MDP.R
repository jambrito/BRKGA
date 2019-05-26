#'Decoder_MDP
#'@title Decoder implemented to maximum diversity problem
#'@description Builds a feasible solution associated with a subset M of N
#'@examples
#'data(Data_MDP1)
#'D<-Data_MDP1
#'N<-dim(D)[1]
#'M<-10
#'RDi<-apply(as.matrix(1:N),1,function(i) sum(D[i,]))/sum((D[upper.tri(D)]))
#'u<-popgen(n=N,1)
#'OD<-list(M=M,RDi=RDi)
#'Xm<-Decoder_MDP(u,OD)
#'@param  u    Chromossome
#'@param  OD   List with Ratio distance matrix and cardinality of M (is necessary to define Exa1 = OD)
#'@return \item{Xm}{= elements of M}
#'@import stats
#'@import utils
#'@export

Decoder_MDP<-function(u,OD)
{
  Rdi<-OD$RDi
  m<-OD$M
  Xm<-order(u*Rdi,decreasing = TRUE)[1:m]
  return(Xm)
}
