#'Fobj_MEDOID
#'@title Objective function of the k-medoids problem
#'@description Computes the sum of the distances for each object from near
#'object associated with medoid
#'@references 1. Kaufman L. e Rousseeuw P.J. (1989). Finding Groups in Data. An Introduction to Cluster Analysis.   Wiley-Interscience Publication.
#'@references 2. Zhang, Q. and Couloigner, I. (2005). A New Efficient k-medoid Algorithm for Spatial Clustering. Lecture Notes in Computer Science, v3482, 181-189.
#'@examples
#'D<-as.matrix(dist(iris[,1:4]))
#'N<-nrow(iris)
#'k<-3 #clusters
#'x<-sample(N,k)  #selected medoids
#'fobj<-Fobj_MEDOID(D,x,k)
#'@param D    distance matrix
#'@param x    medoids
#'@param k    number of clusters (is necessary to define Exa2=k)
#'@return \item{S}{= sum of distances}
#'@import stats
#'@import utils
#'@export

Fobj_MEDOID<-function(D,x,k)
{ n<-dim(D)[1]
objs<-1:n
clusters_S<-t(apply(as.matrix(objs),1,function(ki) D[x,ki]))
clusters<-apply(clusters_S,1,function(x) which.min(x))
ng<-as.numeric(table(clusters))-1
Soma=sum(apply(as.matrix(1:k),1,function(kc) sum(D[which(clusters==kc),x[kc]])/ng[kc]))
return(Soma)
}
