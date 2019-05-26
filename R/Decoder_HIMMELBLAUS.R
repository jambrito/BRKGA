#'Decoder_HIMMELBLAUS
#'@title Decoder implemented to Himmelblau's function
#'@description Calculates x e y values to HimmelBlau's Function
#'@examples
#'u<-popgen(n=2,p=1)
#'XY<-Decoder_HIMMELBLAUS(u)
#'@param  u  chromossome  with size = 2 (number of variables)
#'@return \item{XY}{= vector with x e y values}
#'@import stats
#'@import utils
#'@export

Decoder_HIMMELBLAUS<-function(u)
{
  x=-5+u[1]*10
  y=-5+u[2]*10
  return(c(x,y))
}
