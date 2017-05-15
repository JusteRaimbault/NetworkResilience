
setwd(paste0(Sys.getenv('MONITORAT'),'/L2AnalyseSpatiale/Partiel/NetworkResilience'))

source('functions.R')

realnetworks = c("idf","lacourtine","londonM25","lyon","paris","randstad")
measures = c(gamma,normalizedBetweenness,shortestPathMeasures)#,
#clustCoef,louvainModularity)

library(doParallel)
cl <- makeCluster(24,outfile='logrealcorrs')
registerDoParallel(cl)

startTime = proc.time()[3]

res <- foreach(i=1:(4*length(realnetworks))) %dopar% {
  source('functions.R')
  vals <- bootstrapCorrelation("real",n=0,measures=measures,nbootstrap=2,realname=realnetworks[(i%%6)+1])
  res=list()
  res[[realnetworks[(i%%6)+1]]]<-vals
  return(res)
}

stopCluster(cl)

show(paste0("Ellapsed Time : ",proc.time()[3]-startTime))

save(res,file='realCorrelations_test.RData')


