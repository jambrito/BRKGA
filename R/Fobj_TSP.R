#'Fobj_TSP
#'@title Objective function of the travelling salesman problem
#'@description Computes the sum of the distances based on the route x and
#'D matrix (matrix with the distances between cities)
#'@references 1. Applegate David L. et al. (2006). The travelling salesman problem: a computational Study. Princeton: Princeton University Press.
#'@references 2. Garey, M. and Johnson, D. S (1990). Computer and Intractibility: A guide to the Theory of NPCompleteness, Freeman, San Francisco.
#'@examples
#'data(Data_capitals)
#'D<-as.matrix(dist(Data_capitals))
#'Ncities<-nrow(Data_capitals)
#'x<-sample(Ncities,Ncities)
#'fobj<-Fobj_TSP(D,x)
#'@param D    Distance matrix
#'@param x    Route
#'@return \item{sD}{= sum of distances}
#'@import stats
#'@import utils
#'@export

Fobj_TSP<-function(D,x)
{
  n<-length(x)
  route<-cbind(x[-n],x[-1])
  sumDist<-sum(D[route])
  return(sumDist)
}

