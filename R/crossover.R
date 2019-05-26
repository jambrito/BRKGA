#'crossover
#'@title Crossover Operator
#'@description This function performs the uniform crossover between two chromosomes
#'@references 1. Gonçalves, J.R. and Resende, M.G.C. (2011). Biased random-key genetic algorithms for combinatorial optimization,
# Journal of Heuristics, 17, p 487-525.
#'@references 2. Spears, W.M., DeJong, K.A.(1991). On the virtues of parameterized uniform crossover.
#'In: Proceedings of the Fourth International Conference on Genetic Algorithms, pp. 230b.
#'@examples X<-popgen(10,10)
#'Y<-crossover(gelite=X[1:3,],gnelite=X[4:10,],rc=0.7,p=10,pe=3,pm=3,n=10)
#'@param gelite   matrix with  elite chromosomes
#'@param gnelite  matrix with nonelite chromosomes
#'@param rc       crossover probability
#'@param p        number of elements (chromosomes) in the population
#'@param pe       number of elite  chromosomes
#'@param pm       number of mutant chromossomes
#'@param n        number of genes in the chromosomes associated with a solution
#'@return \item{gnew}{= matrix pxn with chromosomes produced by crossover}
#'@import stats
#'@export


crossover<-function(gelite,gnelite,rc,p,pe,pm,n)
{
  uniform<-function(ge,gn,ab,rc,n)
  {
    pr<-runif(n,0,1)
    gp<-which(pr<rc)
    ngp<-setdiff((1:n),gp)
    gn<-rep(0,n)
    gn[gp]<-gelite[ab[1],gp]
    gn[ngp]<-gnelite[ab[2],ngp]
    return(gn)
  }
  pnelite<-p-pe
  psort<-p-pe-pm
  ab=cbind(sample(pe,psort,replace=TRUE),sample(pnelite,psort,replace=TRUE))
  gnew=t(apply(ab,1,function(ab) uniform(gelite,gnelite,ab,rc,n)))
  return(gnew)
}
