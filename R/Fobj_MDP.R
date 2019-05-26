#'Fobj_MDP
#'@title Objective function of the maximum diversity problem
#'@description Calculates the sum of all distances associated with M
#'@references 1. Martí, R., Gallego, M., Duarte, A. e Pardo, E. G. (2013) Heuristics and metaheuristics for the maximum diversity problem. Journal of Heuristics 19(4): 591-615.
#'@references 2. Silva, G.C., Ochi, L.S. e Martins, S.L. (2004).
#'Experimental Comparison of Greedy Randomized Adaptive Search Procedures for the Maximum Diversity Problem.
#'In Experimental and Efficient Algorithms, volume 3059 of Lecture Notes in Computer Science, p. 498b-512.
#'Springer Berlin / Heidelberg.
#'@examples
#'data(DataMDP1)
#'D<-Data_MDP1
#'N<-dim(D)[1]
#'M<-10
#'xm<-sample(N,M)
#'SD<-Fobj_MDP(D,xm)

#'@param  D  distance matrix
#'@param  Xm elements of M
#'@return \item{SD}{= sum all distances associated with M}
#'@import stats
#'@import utils
#'@export

Fobj_MDP<-function(D,Xm)
{
  XC<-t(combn(Xm,2))
  SD<-sum(D[XC])
  return(SD)
}
