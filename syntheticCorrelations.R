
setwd(paste0(Sys.getenv('MONITORAT'),'/L2AnalyseSpatiale/Partiel/NetworkResilience'))

source('functions.R')

networktypes = c("lattice","pa-age","random","tree")
measures = c(gamma,normalizedBetweenness,shortestPathMeasures)

library(doParallel)
cl <- makeCluster(32,outfile='logsynthcorrs')
registerDoParallel(cl)

startTime = proc.time()[3]

res <- foreach(i=1:(8*length(networktypes))) %dopar% {
  source('functions.R')
  vals <- bootstrapCorrelation(type=networktypes[(i%%4)+1],n=1000,measures=measures,nbootstrap=100)
  res=list()
  res[[networktypes[(i%%4)+1]]]<-vals
  return(res)
}

stopCluster(cl)

show(paste0("Ellapsed Time : ",proc.time()[3]-startTime))

save(res,file='synthCorrelations.RData')
