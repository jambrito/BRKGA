#'Decoder_KNAPSACK
#'@title Decoder implemented to knapsack problem
#'@description Determines items that will be loaded in the knapsack.
#'@references 1. Papadimitriou, C. H., and K. Steiglitz (1998). Combinatorial Optimization: Algorithms and Complexity. Dover Publications.
#'@examples
#'nitems=7
#'u<-popgen(n=nitems,p=1)
#'xis<-Decoder_KNAPSACK(u)
#'@param  u chromossome  (with n positions associated with items)
#'@return \item{xis}{= vector with 0 and 1s associated items of the knapsack}
#'@import stats
#'@import utils
#'@export

Decoder_KNAPSACK<-function(u)
{
  xi<-round(u)
  return(xi)
}


