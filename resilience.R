

##
# Complex networks : measures and vulnerability

library(igraph)

setwd(paste0(Sys.getenv('MONITORAT'),'/L2AnalyseSpatiale/Partiel/NetworkResilience'))

source('functions.R')

n = 100
# real road network ; optimized network ?
networktypes = c("lattice","pa-age","real")
realnetworks = c("idf","lacourtine","londonM25","lyon","paris","randstad")
indicators = c("stats","gamma","betweenness","closeness","transitivity","efficiency","modularity")


## tests
g=generateNetwork("lattice",n)
bw = betweenness(g)
bw*2/(vcount(g)*(vcount(g)-1))

g=generateNetwork("pa-age",n)
transitivity(g)

alphabws = c();bws=c()
for(b in 1:10000){
  g=generateNetwork("pa-age",n);bw=normalizedBetweenness(g)
  alphabws=append(alphabws,bw$alphaBetweenness);bws=append(bws,bw$meanBetweenness)
}

g=generateNetwork("real",realname="idf")

deltaMeasure(g,c(0.1,0.2,0.3),normalizedBetweenness)
deltaMeasure(g,c(0.1,0.2,0.3),efficiency)


# Q : how do ∆measure correlate with ∆efficiency (= measure of resilience ?)
removals = seq(from=0.05,to=0.5,by=0.05)
measure=normalizedBetweenness
nbootstrap = 1000
deltaEff=c();deltaM=c()
for(b in 1:nbootstrap){
  #g=generateNetwork("pa-age",n)
  g=generateNetwork("lattice",n)
  deltaEff=append(deltaEff,deltaMeasure(g,removals,efficiency))
  deltaM=append(deltaM,deltaMeasure(g,removals,measure))
}

cor.test(deltaEff,deltaM)






