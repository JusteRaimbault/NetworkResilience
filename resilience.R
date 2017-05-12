

##
# Complex networks : measures and vulnerability

library(igraph)

setwd(paste0(Sys.getenv('MONITORAT'),'/L2AnalyseSpatiale/Partiel/NetworkResilience'))

source('functions.R')

n = 1000
# real road network ; optimized network ?
networktypes = c("lattice","pa-age","random","tree","real")
realnetworks = c("idf","lacourtine","londonM25","lyon","paris","randstad")
indicators = c("stats","gamma","betweenness","diameter","closeness","transitivity","efficiency","modularity")
measures = c(gamma,normalizedBetweenness,shortestPathMeasures,clustCoef,louvainModularity)

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
paris=generateNetwork("real",realname="paris")
randstad=generateNetwork("real",realname="randstad")
lacourtine=generateNetwork("real",realname="lacourtine")

res <- computeDeterministic(lacourtine,measures)

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


# plots
g=generateNetwork("real",realname="paris")
png(filename = 'figures/paris.png',width = 30,height=30,units = 'cm',res=600)
plot(g,vertex.size=0.05,vertex.label=NA,edge.color='black',edge.width=1.0)
dev.off()

g=generateNetwork("real",realname="idf")
png(filename = 'figures/idf.png',width = 30,height=30,units = 'cm',res=600)
plot(g,vertex.size=0.05,vertex.label=NA,edge.color='black',edge.width=1.0)
dev.off()

g=generateNetwork("real",realname="lacourtine")
png(filename = 'figures/lacourtine.png',width = 30,height=30,units = 'cm',res=600)
plot(g,vertex.size=0.05,vertex.label=NA,edge.color='black',edge.width=1.0)
dev.off()

g=generateNetwork("real",realname="londonM25")
png(filename = 'figures/londonM25.png',width = 30,height=30,units = 'cm',res=600)
plot(g,vertex.size=0.05,vertex.label=NA,edge.color='black',edge.width=1.0)
dev.off()

g=generateNetwork("real",realname="lyon")
png(filename = 'figures/lyon.png',width = 30,height=30,units = 'cm',res=600)
plot(g,vertex.size=0.05,vertex.label=NA,edge.color='black',edge.width=1.0)
dev.off()

g=generateNetwork("real",realname="randstad")
png(filename = 'figures/randstad.png',width = 30,height=30,units = 'cm',res=600)
plot(g,vertex.size=0.05,vertex.label=NA,edge.color='black',edge.width=1.0)
dev.off()

g=generateNetwork("random",n=1000)
png(filename = 'figures/random.png',width = 30,height=30,units = 'cm',res=600)
plot(g,vertex.size=0.05,vertex.label=NA,edge.color='black',edge.width=1.0,
     layout=layout_with_fr,margin=c(-0.2,-0.2,-0.2,-0.2)
     )
dev.off()

g=generateNetwork("pa-age",n=1000)
png(filename = 'figures/pa-age.png',width = 30,height=30,units = 'cm',res=600)
plot(g,vertex.size=0.05,vertex.label=NA,edge.color='black',edge.width=1.0
     ,layout=layout_with_fr#,margin=c(-0.2,-0.2,-0.2,-0.2)
     )
dev.off()

g=generateNetwork("lattice",n=1000)
png(filename = 'figures/lattice.png',width = 30,height=30,units = 'cm',res=600)
plot(g,vertex.size=0.05,vertex.label=NA,edge.color='black',edge.width=1.0
     #,layout=layout_on_grid
     )
dev.off()

g=generateNetwork("tree",n=1000)
png(filename = 'figures/tree.png',width = 30,height=30,units = 'cm',res=600)
plot(g,vertex.size=0.05,vertex.label=NA,edge.color='black',edge.width=1.0,
     layout=layout_as_tree
     )
dev.off()



####
#  test bootstrap

test <- bootstrapMeasures("random",1000,measures,1000)





