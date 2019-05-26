#'Decoder_MSDCP
#'@title  Decoder implemented to minimum sum distance clustering problem
#'@description Builds a clustering vector
#'@examples
#'N<-nrow(iris)
#'k<-3
#'u<-popgen(n=N,p=1)
#'cluster<-Decoder_MSDCP(u,k)
#'@param  u  chromossome
#'@param  k  number of clusters (is necessary to define Exa1 = k)
#'@return \item{cluster}{ = clustering vector}
#'@import stats
#'@import utils
#'@export

Decoder_MSDCP<-function(u,k)
{
  cluster<-findInterval(u,seq(0,1,1/k))
  return(cluster)
}
