# bootstrap synthetic networks

setwd(paste0(Sys.getenv('MONITORAT'),'/L2AnalyseSpatiale/Partiel/NetworkResilience'))

source('functions.R')

networktypes = c("lattice","pa-age","random","tree")
measures = c(gamma,normalizedBetweenness,shortestPathMeasures,clustCoef,louvainModularity)

library(doParallel)
cl <- makeCluster(50,outfile='logsynth')
registerDoParallel(cl)

startTime = proc.time()[3]

# 50 bootstraps in //
res <- foreach(i=1:50) %dopar% {
  source('functions.R')
  
  res=list()
  for(type in networktypes){
    vals = bootstrapMeasures(type,1000,measures,10000)
    res[[type]]<-vals
  }
  
  return(res)
}

stopCluster(cl)

show(paste0("Ellapsed Time : ",proc.time()[3]-startTime))

save(res,file='synthetic.RData')

