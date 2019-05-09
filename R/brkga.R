#'popgen
#'@title Population Generation
#'@description Generates p vectors of random-keys
#'@references 1. Gonçalves, J.R. and Resende, M.G.C. (2011). Biased random-key genetic algorithms for combinatorial optimization,
# Journal of Heuristics, 17, p 487-525.
#'@param n Number of genes in the chromosome associated with a solution
#'@param p Number of elements (chromosomes) in the population
#'@return \item{pop}{The population is stored in the p x n real-valued matrix}
#'@import stats
#'@export
popgen<-function(n,p)
{
pop<-matrix(runif(n*p),nrow=p,ncol=n)
return(pop)
}

#'crossover
#'@title Crossover Operator
#'@description This function performs the uniform crossover between two chromosomes
#'@references 1. Gonçalves, J.R. and Resende, M.G.C. (2011). Biased random-key genetic algorithms for combinatorial optimization,
# Journal of Heuristics, 17, p 487-525.
#'@references 2. Spears, W.M., DeJong, K.A.(1991). On the virtues of parameterized uniform crossover.
#'In: Proceedings of the Fourth International Conference on Genetic Algorithms, pp. 230b.
#'@examples X<-popgen(10,10)
#'Y<-crossover(gelite=X[1:3,],gnelite=X[4:10,],rc=0.7,p=10,pe=3,pm=3,n=10)
#'@param gelite   Matrix with  elite chromosomes
#'@param gnelite  Matrix with nonelite chromosomes
#'@param rc       Crossover probability
#'@param p        Number of elements (chromosomes) in the population
#'@param pe       Number of elite  chromosomes
#'@param pm       Number of mutant chromossomes
#'@param n        Number of genes in the chromosomes associated with a solution
#'@return \item{gnew}{Matrix pxn with chromosomes produced by crossover}
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


#'brkga
#'@title BRKGA Algorithm
#'@description This function applies BRKGA algorithm to a problem considering
#'Objective Function and Decoder defined by user
#'@references 1. Gonçalves, J.F. and Resende, M.G.C. (2011).
#'Biased random-key genetic algorithms for combinatorial optimization,
#'Journal of Heuristics, 17, p 487-525.
#'@references 2. Gonçalves, J.F and Resende M.G.C. (2018). Biased random-key genetic programming.
#'Handbook of Heuristics, Edited by: Martí, R., Pardalos, P.M. Pardalos and Resende M.G.C.
#', pp. 23-37, Springer.
#'@author Jose Brito (jambrito@gmail.com), Gustavo Semaan (gustavosemaan@gmail.com) and Augusto Fadel (augustofadel@gmail.com).
#'@param Data     Vector or Matrix (data to calculate objective function)
#'@param Fo       Objective function defined by user
#'@param Dc       Decoder defined by user
#'@param rc       Crossover probability
#'@param pe       Percentual of elite  chromosomes
#'@param pm       Percentual of mutant chromosomes
#'@param n        Number of genes in the chromosomes associated with a solution
#'@param p        Number of elements (chromosomes) in the population
#'@param ng       Number of generations of the brkga algorithm
#'@param ngw      Number of generations without improvement
#'@param MaxTime  Maximum CPU Time (seconds)
#'@param MAX      Argument that determines maximization (TRUE) or minimization problem (FALSE)
#'@param Exa1     Extra Argument to Decoder
#'@param Exa2     Extra Argument to Objective Function
#'@return \item{fbest}{Best value of Objective Function}
#'@return \item{gbest}{Best Solution}
#'@return \item{cpu_time}{Cpu time in seconds}
#'@import stats
#'@import utils
#'@export
#'
#'@examples
#'#Example1 - TSP - Travelling Salesman Problem
#'data(Data_capitals)
#'D<-as.matrix(dist(Data_capitals))
#'Ncities<-nrow(Data_capitals)
#'s<-brkga(Data=D,Fo=Fobj_TSP,Dc=Decoder_TSP,ng=100,n=Ncities,p=50)
#'
#'#Example2 - Clustering Problem: k-medois
#'Distance<-as.matrix(dist(iris[,1:4]))
#'N<-nrow(iris)
#'k<-3 #Clusters
#\dontrun{s<-brkga(Data=Distance,Fo=Fobj_MEDOID,Dc=Decoder_MEDOID,n=N,p=50,ng=100,MaxTime=2,Exa1=k,Exa2=k)}
#'
#'#Example3 - Knapsack Problem
#'wi<-c(40,50,30,10,10,40,30)
#'li<-c(40,80,10,10,4,20,60)
#'C<-100
#'Datalw<-cbind(li,wi)
#'s<-brkga(Data=Datalw,Fo=Fobj_KNAPSACK,Dc=Decoder_KNAPSACK,n=length(wi),p=10,ng=100,Exa2=C,MAX=TRUE)
#'
#'wi<-c(8,10,39,94,32,88,64,90,20,63,71,7,17,99,77,57,3,26,43,55,95,53,34,62,74)
#'li<-c(23,18,40,87,34,27,44,46,47,65,90,26,6,22,31,12,10,54,71,36,30,63,59,45,79)
#'C<-400
#'Datalw<-cbind(li,wi)
#\dontrun{s<-brkga(Data=Datalw,Fo=Fobj_KNAPSACK,Dc=Decoder_KNAPSACK,n=length(wi),ngw=25,Exa2=C,MAX=TRUE)}
#'
#'#Example4 - MDP Problem - Maximum Diversity Problem
#'data(DataMDP1)
#'D<-Data_MDP1
#'N<-dim(D)[1]
#'M<-10
#'RDi<-apply(as.matrix(1:N),1,function(i) sum(D[i,]))/sum((D[upper.tri(D)]))
#'Exa1<-list(M=M,RDi=RDi)
#\dontrun{s<-brkga(Data=D,Fo=Fobj_MDP,Dc=Decoder_MDP,n=N,p=75,MaxTime=2,Exa1=Exa1,MAX=TRUE)}
#'
#'#Example5 - Optimal Allocation  in Stratified Sampling
#'Nh<-c(212,84,61)
#'Sh2<-c(723.1,2693.4,36231.7)
#'Y<-80548
#'H<-3
#'Exa2<-list(Nh=Nh,Y=Y,cvt=0.05)
#'s<-brkga(Data=Sh2,Fo=Fobj_ALLOC_SAMPLE,Dc=Decoder_ALLOC_SAMPLE,n=H,p=1000,
#'MaxTime=3,Exa1=Nh,Exa2=Exa2,MAX=FALSE)
#'
#'#Example6 - Minimization of Himmelblaus Function
#'nv<-2 #Number of variables
#'s<-brkga(Data=NULL,Fo=Fobj_HIMMELBLAUS,Dc=Decoder_HIMMELBLAUS,n=nv,p=100)
#'
#'#Example7 - MSDCP - Minimum Sum Distance Clustering Problem
#'k<-3 #Clusters
#'Distance<-as.matrix(dist(iris[,1:4]))
#'N<-nrow(iris)
#\dontrun{s<-brkga(Data=Distance,Fo=Fobj_MSDCP,Dc=Decoder_MSDCP,n=N,p=50,MaxTime=7,Exa1=k,Exa2=k)}


brkga<-function(Data,Fo,Dc,rc=0.7,pe=0.2,pm=0.2,n,p=100,ng=2000,ngw=500,MaxTime=3600,MAX=FALSE,Exa1=NULL,Exa2=NULL)
{
  cpu_time<-proc.time()
  time_iter<-0
  f<-popgen(n,p)
  if(is.null(Exa1)==TRUE)
    {g<-t(apply(f,1,function(x) Dc(x)))}
   else {g<-t(apply(f,1,function(x) Dc(x,Exa1)))}
  if(is.null(Exa2)==TRUE)
  { ft<-apply(g,1,function(x) Fo(Data,x))}
  else {ft<-apply(g,1,function(x) Fo(Data,x,Exa2))}

  fbest<-ifelse(MAX==TRUE,-Inf,Inf)
  i<-0
  pelite<-round(pe*p)
  pmutant<-round(pm*p)
  ngwb<-0
  while((i<ng) & (ngwb<=ngw) & (time_iter<MaxTime))
     {i<-i+1
      ngwb<-ngwb+1
      pq<-order(ft,decreasing = MAX)
      f<-f[pq,] #Sorting by Fitness
      g<-g[pq,]
      fmin<-ft[pq[1]]
      ft<-ft[pq]
      if (MAX==FALSE)
       {
        if (fmin<fbest)
          {fbest<-fmin
           gbest<-g[1,]
           solution_best<-c(fbest,gbest)
           cat("Best Solution Generation ",i," = ",fbest,"\n")
           flush.console()
           ibest<-i
           ngwb<-0
          }
      } else
        {
          if (fmin>fbest)
          {fbest<-fmin
          gbest<-g[1,]
          solution_best<-c(fbest,gbest)
          cat("Best Solution Generation ",i," = ",fbest,"\n")
          flush.console()
          ibest<-i
          }

        }
      felite<-f[1:pelite,]
      fnonelite<-f[(pelite+1):p,] #Non-Elite
      fmutant<-popgen(n,pmutant)
      fnovos<-crossover(felite,fnonelite,rc,p,pelite,pmutant,n)
      fnew<-rbind(fmutant,fnovos)
      if(is.null(Exa1)==TRUE)
      { gnew<-t(apply(fnew,1,function(x) Dc(x)))}
      else {gnew<-t(apply(fnew,1,function(x) Dc(x,Exa1)))}

      g<-rbind(g[1:pelite,],gnew)
      f<-rbind(felite,fnew)
      glk<-g[(pelite+1):p,]
      if (is.null(Exa2)==TRUE)
      {ftk<-apply(glk,1,function(x) Fo(Data,x))}
      else
      {ftk<-apply(glk,1,function(x) Fo(Data,x,Exa2))}

      ft<-c(ft[1:pelite],ftk)
      time_iter<-(proc.time()-cpu_time)[3]
     }
  cpu_time<-(proc.time()-cpu_time)[3]
  return(list(fbest=fbest,gbest=gbest,cpu_time=cpu_time))
}


#'Fobj_TSP
#'@title Objective Function of the Travelling Salesman Problem
#'@description Computes the sum of the distances based on the route x and
#'D matrix (matrix with the distances between cities)
#'@references 1. Applegate David L. et al. (2006). The travelling salesman problem: a computational Study. Princeton: Princeton University Press.
#'@references 2. Garey, M. and Johnson, D. S (1990). Computer and Intractibility: A guide to the Theory of NPCompleteness, Freeman, San Francisco.
#'@param D    Distance matrix
#'@param x    Route
#'@return \item{sD}{Sum of distances}
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


#'Decoder_TSP
#'@title Decoder implemented to TSP Problem
#'@description Builds a feasible solution associated a route. Order of travel of n cities.
#'@details TSP Problem is most important problems literatura (see, Kumar, 2012)
#'@references 1. Applegate David L. et al. (2006). The Travelling Salesman Problem: a computational study. Princeton: Princeton University Press.
#'@references 2. Garey, M. and Johnson, D. S (1990). Computer and Intractibility: A guide to the Theory of NPCompleteness, Freeman, San Francisco.
#'@references 3. Kumar, K.(2012). Traveling Salesman Problem (TSP): A Comparative Analysis: Optimization Techniques for Traveling Salesman Problem. LAP LAMBERT Academic Publishing
#'@param  u chromossome
#'@return \item{x}{Route}
#'@import stats
#'@import utils
#'@export
Decoder_TSP<-function(u)
{Xr<-order(u);
 return(c(Xr,Xr[1]))
}


#'Fobj_MEDOID
#'@title Objective Function of the k-medoids Problem
#'@description Computes the sum of the distances for each object from near
#'object associated with medoid
#'@references 1. Kaufman L. e Rousseeuw P.J. (1989). Finding Groups in Data. An Introduction to Cluster Analysis.   Wiley-Interscience Publication.
#'@references 2. Zhang, Q. and Couloigner, I. (2005). A New Efficient k-medoid Algorithm for Spatial Clustering. Lecture Notes in Computer Science, v3482, 181-189.
#'@param D    Distance matrix
#'@param x    Medoids
#'@param k    Number of Clusters
#'@return \item{Soma}{Sum of distances total}
#'@import stats
#'@import utils
#'@export

Fobj_MEDOID<-function(D,x,k)
{ n<-dim(D)[1]
  objs<-1:n
  clusters_S<-t(apply(as.matrix(objs),1,function(ki) D[x,ki]))
  clusters<-apply(clusters_S,1,function(x) which.min(x))
  ng<-as.numeric(table(clusters))-1
  Soma=sum(apply(as.matrix(1:k),1,function(kc) sum(D[which(clusters==kc),x[kc]])/ng[kc]))
  return(Soma)
}

#'Decoder_MEDOID
#'@title Decoder implemented to k-medoids Problem
#'@description Determines k objects associated with medoids
#'@details k-medoids Problem is fully described in Kaufman and Rousseeuw (1989)
#'@references 1. Kaufman L. e Rousseeuw P.J. (1989). Finding Groups in Data. An Introduction to Cluster Analysis.   Wiley-Interscience Publication.
#'@references 2. Zhang, Q. and Couloigner, I. (2005). A New Efficient k-medoid Algorithm for Spatial Clustering. Lecture Notes in Computer Science, v3482, 181-189.
#'@param  u chromossome
#'@param  k number of clusters (is necessary to define Exa1 = k)
#'@return \item{medoids}{Vector with k-medoids}
#'@import stats
#'@import utils
#'@export

Decoder_MEDOID<-function(u,k)
{
It<-as.numeric(cut(u,breaks=k,labels=1:k))
Gk<-apply(as.matrix(1:k),1,function(ki) (which(It==ki)))
if (is.list(Gk))
 {z<-lapply(Gk,function(ei) ei[which.min(abs(mean(u[ei])-u[ei]))])
 medoids<-unlist(z)
}else
  {
  z<-apply(Gk,2,function(ei) ei[which.min(abs(mean(u[ei])-u[ei]))])
  medoids<-z
  }

return(medoids)
}


#'Fobj_KNAPSACK
#'@title Objective Function of the Knapsack Problem
#'@description Computes Sum of Profits based on itens allocated in knapsack
#'@references 1. Papadimitriou, C. H., and K. Steiglitz (1998). Combinatorial Optimization: Algorithms and Complexity. Dover Publications.
#'@references 2. Martello, S. and Toth, P. (1990). Knapsack Problems - Algorithms and Computer Implementations. Wiley.
#'@param D    Distance matrix
#'@param x    Items Knapsack
#'@param C    Capacity of knapsack
#'@return \item{FP}{Profit associated items of knapscak}
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

#'Decoder_Knapsack
#'@title Decoder implemented to Knapsack Problem
#'@description Determines items that will be loaded in the knapsack.
#'@references 1. Papadimitriou, C. H., and K. Steiglitz (1998). Combinatorial Optimization: Algorithms and Complexity. Dover Publications.
#'@param  u chromossome
#'@return \item{xis}{Vector with 0 and 1s associated items of the knapsack}
#'@import stats
#'@import utils
#'@export

Decoder_KNAPSACK<-function(u)
{
xi<-round(u)
return(xi)
}

#'Fobj_MDP
#'@title Objective Function of the Maximum Diversity Problem
#'@description Calculates the sum of all distances associated with M
#'@references 1. Martí, R., Gallego, M., Duarte, A. e Pardo, E. G. (2013) Heuristics and metaheuristics for the maximum diversity problem. Journal of Heuristics 19(4): 591-615.
#'@references 2. Silva, G.C., Ochi, L.S. e Martins, S.L. (2004).
#'Experimental Comparison of Greedy Randomized Adaptive Search Procedures for the Maximum Diversity Problem.
#'In Experimental and Efficient Algorithms, volume 3059 of Lecture Notes in Computer Science, p. 498b-512.
#'Springer Berlin / Heidelberg.

#'@param  D  distance matrix
#'@param  Xm Elements of M
#'@return \item{SD}{Sum all distances associated with M}
#'@import stats
#'@import utils
#'@export

Fobj_MDP<-function(D,Xm)
{
  XC<-t(combn(Xm,2))
  SD<-sum(D[XC])
  return(SD)
}



#'Decoder_MDP
#'@title Decoder implemented to Maximum Diversity Problem
#'@description Builds a feasible solution associated with a subset M of N
#'@param  u    Chromossome
#'@param  OD   List with Ratio distance matrix and cardinality of M (is necessary to define Exa1 = OD)
#'@return \item{Xm}{Elements of M}
#'@import stats
#'@import utils
#'@export

Decoder_MDP<-function(u,OD)
{
  Rdi<-OD$RDi
  m<-OD$M
  Xm<-order(u*Rdi,decreasing = TRUE)[1:m]
  return(Xm)
}


#'Fobj_ALLOC_SAMPLE
#'@title Objective Function of the Optimal Allocation in Stratified Sampling
#'@description Calculates sum of sample sizes associated  with all strata
#'@references 1. Brito, J. A. M.; Silva, P. L. N. ; Semaan, G. S. ; Maculan, N.(2015). Integer Programming Formulations Applied to Optimal Allocation in Stratified Sampling.
#'Survey Methodology, v. 41, p. 427-442.
#'@param  Sh2 Variance by stratum
#'@param  nh Number of units in the sample in the hth stratum
#'@param  S  List with Nh, Y and cvt (cv target)
#'@return \item{n}{Sum of nh}
#'@import stats
#'@import utils
#'@export

Fobj_ALLOC_SAMPLE<-function(Sh2,nh,S)
{
  Nh<-S$Nh
  Y<-S$Y
  cvt<-S$cvt
  cv<-sqrt(sum(Nh^2*(1/nh-1/Nh)*Sh2))/Y
  n<-sum(nh)
  P<-cv/cvt
  if(P>1) {n<-n+2000^P}
  return(n)
}


#'Decoder_ALLOC_SAMPLE
#'@title Decoder implemented to Optimal Allocation in Stratified Sampling
#'@description Builds a solution associated with sample sizes allocated by stratum
#'@param  u  chromossome
#'@param  Nh Total Units in hth stratum (is necessary to define Exa1 = Nh)
#'@return \item{nh}{Number of units in the sample in the hth stratum}
#'@import stats
#'@import utils
#'@export

Decoder_ALLOC_SAMPLE<-function(u,Nh)
{
  nh<-2+round(u*(Nh-2))
  return(nh)
}


#'Fobj_HIMMELBLAUS
#'@title Himmelblau's function
#'@description Calculates value of the Himmelblau's function
#'@references 1. Himmelblau, D. (1972). Applied Nonlinear Programming. McGraw-Hill.
#'@param  Data  Function terms
#'@param  x     Point (x,y)
#'@return \item{F}{Function Value}
#'@import stats
#'@import utils
#'@export

Fobj_HIMMELBLAUS<-function(Data=NULL,x)
{
fobj<-(x[1]^2+x[2]-11)^2+(x[1]+x[2]^2-7)^2
return(fobj)
}


#'Decoder_HIMMELBLAUS
#'@title Decoder implemented to Himmelblau's function
#'@description Calculates x e y values to HimmelBlau's Function
#'@param  u  chromossome
#'@return \item{XY}{Vector with x e y values}
#'@import stats
#'@import utils
#'@export

Decoder_HIMMELBLAUS<-function(u)
{
x=-5+u[1]*10
y=-5+u[2]*10
return(c(x,y))
}


#'Fobj_MSDCP
#'@title  Objective Function of the Minimum Sum Distance Clustering Problem
#'@references 1. Friedman, J. H. and Meulman J.J.  (2004). Clustering objects on subsets of attributes.
#'#' Journal Royal Statistics Society B, Part 4,66:815-849.
#'@references 2. Hansen, P. and Jaumard, B. (1997). Cluster Analysis and Mathematical Programming.
#'#'Mathematical  Programming, 79:191-215.
#'@references 3. Rao, M.R. (1971). Cluster Analysis and Mathematical Programming.
#'Journal of American Statistical Association, 66:622-626.
#'@param  D           Distance Matrix
#'@param  cluster     Clustering Vector
#'@param  k           Number of clusters
#'@return \item{SD}{Sum Distances}
#'@import stats
#'@import utils
#'@export

Fobj_MSDCP<-function(D,cluster,k)
{
dc<-rep(0,k)
for(i in 1:k)
 {
  wc=which(cluster==i)
  cb<-t(combn(wc,2))
  dc[i]<-sum(D[cb])
 }
return(sum(dc))
}


#'Decoder_MSDCP
#'@title  Decoder implemented to Minimum Sum Distance Clustering Problem
#'@description Builds a clustering vector
#'@param  u  chromossome
#'@param  k  number of clusters (is necessary to define Exa1 = k)
#'@return \item{cluster}{clustering vector}
#'@import stats
#'@import utils
#'@export

Decoder_MSDCP<-function(u,k)
{
cluster<-findInterval(u,seq(0,1,1/k))
return(cluster)
}






