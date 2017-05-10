
# construct networks

setwd(paste0(Sys.getenv('MONITORAT'),'/L2AnalyseSpatiale/Partiel/NetworkResilience'))

source('functions.R')

extentfiles = c("idf","lacourtine","londonM25","lyon","paris","randstad")

for(file in extentfiles){
  g = getRoadNetwork(paste0("data/",file))
  save(g,file=paste0('data/',file,'.RData'))
}


