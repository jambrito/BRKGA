#'Fobj_KNAPSACK
#'@title Objective function of the knapsack problem
#'@description Computes sum of profits based on itens allocated in knapsack
#'@references 1. Papadimitriou, C. H., and K. Steiglitz (1998). Combinatorial Optimization: Algorithms and Complexity. Dover Publications.
#'@references 2. Martello, S. and Toth, P. (1990). Knapsack Problems - Algorithms and Computer Implementations. Wiley.
#'@examples
#'wi<-c(40,50,30,10,10,40,30)
#'li<-c(40,80,10,10,4,20,60)
#'x<-c(1,0,0,1,1,0,1)
#'C<-100
#'Datalw<-cbind(li,wi)
#'FP<-Fobj_KNAPSACK(Datalw,x,C)
#'@param D    matrix nx2 with weights (wi) and profits (li)
#'@param x    items Knapsack
#'@param C    capacity of knapsack (is necessary to define Exa2=C)
#'@return \item{FP}{= profit associated items of knapscak}
#'@import stats
#'@import utils
#'@export

Fobj_KNAPSACK<-function(D,x,C)
{
  li<-D[,1]
  wi<-D[,2]
  gi<-x%*%wi
  P<-max(0,(gi-C))^3
  FP<-li%*%x-P
  return(FP)
}
