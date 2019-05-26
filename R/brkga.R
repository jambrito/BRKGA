#'brkga
#'@title BRKGA Algorithm
#'@description This function applies brkga algorithm to a problem considering
#'objective function and decoder defined by user
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
#'s<-brkga(Data=D,Fo=Fobj_TSP,Dc=Decoder_TSP,ng=20,n=Ncities,p=50)
#'
#'#Example2 - Clustering Problem: k-medois
#'Distance<-as.matrix(dist(iris[,1:4]))
#'N<-nrow(iris)
#'k<-3 #Clusters
#'s<-brkga(Data=Distance,Fo=Fobj_MEDOID,Dc=Decoder_MEDOID,n=N,p=25,ng=20,MaxTime=2,Exa1=k,Exa2=k)
#'
#'#Example3 - Knapsack Problem
#'wi<-c(40,50,30,10,10,40,30)
#'li<-c(40,80,10,10,4,20,60)
#'C<-100
#'Datalw<-cbind(li,wi)
#'s<-brkga(Data=Datalw,Fo=Fobj_KNAPSACK,Dc=Decoder_KNAPSACK,n=length(wi),p=50,ng=20,Exa2=C,MAX=TRUE)
#'
#'wi<-c(8,10,39,94,32,88,64,90,20,63,71,7,17,99,77,57,3,26,43,55,95,53,34,62,74)
#'li<-c(23,18,40,87,34,27,44,46,47,65,90,26,6,22,31,12,10,54,71,36,30,63,59,45,79)
#'C<-400
#'Datalw<-cbind(li,wi)
#'s<-brkga(Data=Datalw,Fo=Fobj_KNAPSACK,Dc=Decoder_KNAPSACK,n=length(wi),ng=50,ngw=25,Exa2=C,MAX=TRUE)
#'
#'#Example4 - MDP Problem - Maximum Diversity Problem
#'data(Data_MDP1)
#'D<-Data_MDP1
#'N<-dim(D)[1]
#'M<-10
#'RDi<-apply(as.matrix(1:N),1,function(i) sum(D[i,]))/sum((D[upper.tri(D)]))
#'Exa1<-list(M=M,RDi=RDi)
#'s<-brkga(Data=D,Fo=Fobj_MDP,Dc=Decoder_MDP,n=N,ng=10,p=50,MaxTime=2,Exa1=Exa1,MAX=TRUE)
#'
#'#Example5 - Optimal Allocation  in Stratified Sampling
#'Nh<-c(212,84,61)
#'Sh2<-c(723.1,2693.4,36231.7)
#'Y<-80548
#'H<-3
#'Exa2<-list(Nh=Nh,Y=Y,cvt=0.05)
#'s<-brkga(Data=Sh2,Fo=Fobj_ALLOC_SAMPLE,Dc=Decoder_ALLOC_SAMPLE,n=H,p=1000,ng=50,
#'MaxTime=2,Exa1=Nh,Exa2=Exa2,MAX=FALSE)
#'
#'#Example6 - Minimization of Himmelblaus Function
#'nv<-2 #Number of variables
#'s<-brkga(Data=NULL,Fo=Fobj_HIMMELBLAUS,Dc=Decoder_HIMMELBLAUS,n=nv,p=100,ng=200)
#'
#'#Example7 - MSDCP - Minimum Sum Distance Clustering Problem
#'k<-3 #Clusters
#'Distance<-as.matrix(dist(iris[,1:4]))
#'N<-nrow(iris)
#'s<-brkga(Data=Distance,Fo=Fobj_MSDCP,Dc=Decoder_MSDCP,n=N,p=50,ng=10,MaxTime=2,Exa1=k,Exa2=k)


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












