#'Decoder_MEDOID
#'@title Decoder implemented to k-medoids Problem
#'@description Determines k objects associated with medoids
#'@details k-medoids problem is fully described in Kaufman and Rousseeuw (1989)
#'@references 1. Kaufman L. e Rousseeuw P.J. (1989). Finding Groups in Data. An Introduction to Cluster Analysis.   Wiley-Interscience Publication.
#'@references 2. Zhang, Q. and Couloigner, I. (2005). A New Efficient k-medoid Algorithm for Spatial Clustering. Lecture Notes in Computer Science, v3482, 181-189.
#'@examples
#'N<-nrow(iris)
#'u<-popgen(n=N,p=1)
#'k<-3
#'medoids<-Decoder_MEDOID(u,k)
#'@param  u chromossome
#'@param  k number of clusters (is necessary to define Exa1 = k)
#'@return \item{medoids}{= vector with k-medoids (objects) }
#'@import stats
#'@import utils
#'@export

Decoder_MEDOID<-function(u,k)
{
  It<-as.numeric(cut(u,breaks=k,labels=1:k))
  Gk<-apply(as.matrix(1:k),1,function(ki) (which(It==ki)))
  if (is.list(Gk))
  {z<-lapply(Gk,function(ei) ei[which.min(abs(mean(u[ei])-u[ei]))])
  medoids<-unlist(z)
  }else
  {
    z<-apply(Gk,2,function(ei) ei[which.min(abs(mean(u[ei])-u[ei]))])
    medoids<-z
  }

  return(medoids)
}
