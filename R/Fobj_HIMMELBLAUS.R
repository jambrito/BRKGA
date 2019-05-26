#'Fobj_HIMMELBLAUS
#'@title Himmelblau's function
#'@description Calculates value of the Himmelblau's function
#'@references 1. Himmelblau, D. (1972). Applied Nonlinear Programming. McGraw-Hill.
#'@examples
#'F<-Fobj_HIMMELBLAUS(x=c(1.1,1.5))
#'@param  Data  Function terms
#'@param  x     Point (x,y)
#'@return \item{F}{= function Value}
#'@import stats
#'@import utils
#'@export

Fobj_HIMMELBLAUS<-function(Data=NULL,x)
{
  fobj<-(x[1]^2+x[2]-11)^2+(x[1]+x[2]^2-7)^2
  return(fobj)
}
