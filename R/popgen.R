#'popgen
#'@title Population Generation
#'@description Generates p vectors of random-keys
#'@references 1. Gonçalves, J.R. and Resende, M.G.C. (2011). Biased random-key genetic algorithms for combinatorial optimization,
# Journal of Heuristics, 17, p 487-525.
#'@examples X<-popgen(n=20,p=10)
#'@param n number of genes in the chromosome associated with a solution
#'@param p number of elements (chromosomes) in the population
#'@return \item{pop}{ = p x n real-valued matrix}
#'@import stats
#'@export
popgen<-function(n,p)
{
  pop<-matrix(runif(n*p),nrow=p,ncol=n)
  return(pop)
}
