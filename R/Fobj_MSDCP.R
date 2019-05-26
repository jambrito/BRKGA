#'Fobj_MSDCP
#'@title  Objective function of the minimum sum distance clustering problem
#'@references 1. Friedman, J. H. and Meulman J.J.  (2004). Clustering objects on subsets of attributes.
#'#' Journal Royal Statistics Society B, Part 4,66:815-849.
#'@references 2. Hansen, P. and Jaumard, B. (1997). Cluster Analysis and Mathematical Programming.
#'#'Mathematical  Programming, 79:191-215.
#'@references 3. Rao, M.R. (1971). Cluster Analysis and Mathematical Programming.
#'Journal of American Statistical Association, 66:622-626.
#'@examples
#'k<-3 #Clusters
#'Distance<-as.matrix(dist(iris[,1:4]))
#'N<-nrow(iris)
#'u<-popgen(n=N,p=1)
#'cluster<-Decoder_MSDCP(u,k)
#'SD<-Fobj_MSDCP(Distance,cluster,k)
#'@param  D           distance Matrix
#'@param  cluster     clustering Vector
#'@param  k           number of clusters (is necessary to define Exa2=k)
#'@return \item{SD}{= sum of distances}
#'@import stats
#'@import utils
#'@export

Fobj_MSDCP<-function(D,cluster,k)
{
  dc<-rep(0,k)
  for(i in 1:k)
  {
    wc=which(cluster==i)
    cb<-t(combn(wc,2))
    dc[i]<-sum(D[cb])
  }
  return(sum(dc))
}
