#'Decoder_TSP
#'@title Decoder implemented to tsp problem
#'@description Builds a feasible solution associated a route (order of travel of n cities).
#'@details TSP
#'@references 1. Applegate David L. et al. (2006). The Travelling Salesman Problem: a computational study. Princeton: Princeton University Press.
#'@references 2. Garey, M. and Johnson, D. S (1990). Computer and Intractibility: A guide to the Theory of NPCompleteness, Freeman, San Francisco.
#'@references 3. Kumar, K.(2012). Traveling Salesman Problem (TSP): A Comparative Analysis: Optimization Techniques for Traveling Salesman Problem. LAP LAMBERT Academic Publishing
#'@examples
#'data(Data_capitals)
#'Ncities<-nrow(Data_capitals)
#'u<-popgen(n=Ncities,p=1)
#'x<-Decoder_TSP(u)
#'@param  u chromossome (n positions)
#'@return \item{x}{= Route with sequence of cities}
#'@import stats
#'@import utils
#'@export
Decoder_TSP<-function(u)
{Xr<-order(u);
return(c(Xr,Xr[1]))
}
