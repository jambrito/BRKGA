#'Decoder_ALLOC_SAMPLE
#'@title Decoder implemented to optimal allocation in stratified sampling
#'@description Builds a solution associated with sample sizes allocated by stratum
#'@examples
#'H<-3
#'u<-popgen(n=H,p=1)
#'Nh<-c(212,84,61)
#'nh<-Decoder_ALLOC_SAMPLE(u,Nh)
#'@param  u  chromossome with H position (number of strata)
#'@param  Nh total Units in hth stratum (is necessary to define Exa1 = Nh)
#'@return \item{nh}{= number of units in the sample in the hth stratum}
#'@import stats
#'@import utils
#'@export

Decoder_ALLOC_SAMPLE<-function(u,Nh)
{
  nh<-2+round(u*(Nh-2))
  return(nh)
}
