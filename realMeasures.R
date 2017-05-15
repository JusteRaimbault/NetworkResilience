
# compute measures on real networks

setwd(paste0(Sys.getenv('MONITORAT'),'/L2AnalyseSpatiale/Partiel/NetworkResilience'))

source('functions.R')

realnetworks = c("idf","lacourtine","londonM25","lyon","paris","randstad")
measures = c(gamma,normalizedBetweenness)#,shortestPathMeasures,
             #clustCoef,louvainModularity)

library(doParallel)
cl <- makeCluster(6,outfile='log')
registerDoParallel(cl)

startTime = proc.time()[3]

res <- foreach(i=1:length(realnetworks)) %dopar% {
  source('functions.R')
  g=generateNetwork("real",realname=realnetworks[i])
  vals <- computeDeterministic(g,measures)
  res=list()
  res[[realnetworks[i]]]<-vals
  return(res)
}

stopCluster(cl)

show(paste0("Ellapsed Time : ",proc.time()[3]-startTime))

save(res,file='real.RData')


