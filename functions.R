
#library(RPostgreSQL)
#library(rgeos)
#library(rgdal)
#library(raster)
library(igraph)


# functions

defaultParams<-function(){
  return(list(
    erdosrenyiProba = 0.005,
    latticeEdgesProportion = 0.65,
    pa.m=5,pa.exp=1,aging.exp=-2,aging.bin=100,
    tree.children = 3
  ))
}


#'
#' @description random network of the given type and number of nodes
generateNetwork<-function(type,n=0,params=defaultParams(),realname=""){
  
  if(type=="random"){
    return(erdos.renyi.game(n=n,params$erdosrenyiProba , type ="gnp"))
  }
  
  if(type=="lattice"){
    # already perturbate the lattice by sampling edges
    g=make_lattice(length=floor(sqrt(n)),dim=2,directed = FALSE)
    # get the layout before removing edges
    layout = layout_on_grid(g);V(g)$x=layout[,1];V(g)$y=layout[,2]
    return(subgraph.edges(g,eids = sample.int(n = vcount(g),size=floor(vcount(g)*params$latticeEdgesProportion),replace = F)))
  }
  
  if(type=="pa-age"){
    return(simplify(sample_pa_age(n,m=params$pa.m, params$pa.exp, aging.exp=params$aging.exp, aging.bin=params$aging.bin,directed = FALSE)))
  }
  
  if(type=="tree"){
    return(make_tree(n=n,children=params$tree.children,mode="undirected"))
  }
  
  if(type=="real"){
    load(paste0('data/',realname,'.RData'))
    return(g)
  }
  
}

gamma<-function(g){
  return(list(
    vcount=vcount(g),
    ecount=ecount(g),
    gamma = 2*ecount(g)/(vcount(g)*(vcount(g)-1)),
    meanDegree = mean(degree(g))
  )
  )
}

normalizedBetweenness<-function(g){
  bw = edge_betweenness(g)*2/(vcount(g)*(vcount(g)-1))
  y=sort(log(bw),decreasing=T)
  reg = lm(data=data.frame(x=log(1:length(which(is.finite(y)))),y=y[is.finite(y)]),formula = y~x)
  return(
    list(
      meanBetweenness = mean(bw),
      alphaBetweenness = reg$coefficients[2]
    )
  )
}


#'
#' @description includes closeness, efficiency, and diameter
#'      note : distances are not weighted for comparison purposes between synthetic and real
#'      , meaning that we consider only topological distance.
shortestPathMeasures<-function(g){
  distmat = distances(g)
  distmatfinite = distmat
  distmatfinite[!is.finite(distmatfinite)]=0
  # get diameter
  diameter = max(distmatfinite)
  #show(diameter)
  # get closeness
  closenesses = (vcount(g)-1) / rowSums(distmatfinite[rowSums(distmatfinite)>0,])
  #show(closenesses)
  y=sort(log(closenesses),decreasing=T)
  reg = lm(data=data.frame(x=log(1:length(which(is.finite(y)))),y=y[is.finite(y)]),formula = y~x)
  # compute efficiency
  diag(distmat)<-Inf
  efficiency=mean(1/distmat)
  return(list(
    diameter=diameter,
    efficiency=efficiency,
    meanCloseness=mean(closenesses),
    alphaCloseness=reg$coefficients[2]
  ))
}

#'
#'
clustCoef<-function(g){
  return(list(transitivity=transitivity(g)))
}

#'
#'
louvainModularity<-function(g){
  com=cluster_louvain(g)
  return(list(
    modularity = max(com$modularity)
  ))
}



#'
#' @description each function in measure should return a named list
#'    in bootstrap cases, returns 
bootstrapMeasures <- function(type,n,measures,nbootstrap,allValues = F){
    vals = list()
    for(b in 1:nbootstrap){
      if(b%%10==0){show(b)}
      g = generateNetwork(type,n)
      currentvals = computeDeterministic(g,measures)
      for(measure in names(currentvals)){
        if(!measure%in%names(vals)){vals[[measure]]=c(currentvals[[measure]])}
        else{vals[[measure]]=append(vals[[measure]],currentvals[[measure]])}   
      }
    }
    show(names(vals))
    if(allValues==T){return(vals)}
    else{
      # compute mean and sd
      res=list()
      for(measure in names(vals)){
        res[[measure]]=mean(vals[[measure]])
        res[[paste0(measure,"Sd")]]=sd(vals[[measure]])
      }
      return(res)
    }
}


#'
#' @description computes mesures on a fixed network
computeDeterministic<-function(g,measures){
  res=list()
  for(measure in measures){
    val = measure(g)
    for(resname in names(val)){
      res[[resname]]=val[[resname]]
    }
  }
  return(res)
}


#deltaMeasure <- function(g,removed_prop,func){
#  res=c()
#  for(alpha in removed_prop){
#    subgraph = subgraph.edges(g,sample.int(n = ecount(g),size = (1-alpha)*ecount(g),replace = FALSE),delete.vertices = F)
#    res=append(res,func(subgraph)-func(g))
#  }
#  return(res)
#}

deltaMeasures<-function(v1,v2){
  res = list()
  for(measure in names(v1)){
    res[[measure]]=v2[[measure]]-v1[[measure]]
  }
  return(res)
}

#'
#' @description estimate correlation between measures variations
bootstrapCorrelation <- function(type,n,measures,nbootstrap){
  vals = list()
  for(b in 1:nbootstrap){
    if(b%%10==0){show(b)}
    g = generateNetwork(type,n)
    baselinevals = computeDeterministic(g,measures)
    for(alpha in seq(from=0.05,to=0.5,by=0.05)){
      subgraph = subgraph.edges(g,sample.int(n = ecount(g),size = (1-alpha)*ecount(g),replace = FALSE),delete.vertices = F)
      currentvals = computeDeterministic(subgraph,measures)
      deltavals = deltaMeasures(baselinevals,currentvals)
      for(measure in names(deltavals)){
        if(!measure%in%names(vals)){vals[[measure]]=c(deltavals[[measure]])}
        else{vals[[measure]]=append(vals[[measure]],deltavals[[measure]])}  
      }
    }
  }
  # estimate correlation
  res=list()
  for(measure in names(vals)){
    if(measure!="efficiency"){
      corr = cor.test(vals[[measure]],vals$efficiency)
      res[[measure]] = corr$estimate
      res[[paste0(measure,"_inf")]] = corr$conf.int[1]
      res[[paste0(measure,"_sup")]] = corr$conf.int[2]
    }
  }
  return(res)
}





#########
## functions to retrieve real nw from postgis database

#'
#' global.dbport=5433;global.dbuser="juste";global.dbhost="";global.nwdb='nw_simpl_4'
graphEdgesFromBase<-function(lonmin,latmin,lonmax,latmax,dbname='nw_simpl_4',dbport=5433,dbuser="juste",dbhost=""){
  pgsqlcon = dbConnect(dbDriver("PostgreSQL"), dbname=dbname,user=dbuser,port=dbport,host=dbhost)
  q = paste0(
    "SELECT origin,destination,length,speed,roadtype FROM links",
    " WHERE ST_Intersects(ST_MakeEnvelope(",lonmin,",",latmin,",",lonmax,",",latmax,",4326),","geography);")
  show(q)
  query = dbSendQuery(pgsqlcon,q)
  data = fetch(query,n=-1)
  dbDisconnect(pgsqlcon)
  if(length(data)==0){return(list())}
  res=list(edgelist=data.frame(from=data$origin,to=data$destination),speed=data$speed,type=data$speed,length=data$length)
  return(res)
}

#'
#' called as 
#' g = graphFromEdges(
#'    graphEdgesFromBase(lonmin,latmin,lonmax,latmax,dbname=global.nwdb),
#'    densraster,from_query = FALSE)
#' densraster <- raster(paste0(Sys.getenv("CN_HOME"),"/Data/PopulationDensity/raw/density_wgs84.tif"))
graphFromEdges<-function(edgelist,densraster,from_query=TRUE){
  if(is.null(edgelist$edgelist)){return(make_empty_graph())}
  if(from_query==TRUE){edgesmat=matrix(data=as.character(unlist(edgelist$edgelist)),ncol=2,byrow=TRUE);}
  else{edgesmat=edgelist$edgelist}
  g = graph_from_data_frame(data.frame(edgesmat,speed=edgelist$speed,type=edgelist$type),directed=FALSE)
  gcoords = xyFromCell(densraster,as.numeric(V(g)$name))
  V(g)$x=gcoords[,1];V(g)$y=gcoords[,2]
  E(g)$length=edgelist$length
  gg=simplify(g,edge.attr.comb="min")
  return(gg)
}


#'
#' @description get graph object given a spatial extent
getRoadNetwork <- function(extentfile){
  densraster <- raster(paste0(Sys.getenv("CN_HOME"),"/Data/PopulationDensity/raw/density_wgs84.tif"))
  path = strsplit(x = extentfile,split="/")[[1]]
  extent <- readOGR(paste0(path[1:(length(path)-1)],collapse = "/"),path[length(path)])
  bbox = bbox(extent)
  lonmin = bbox[1,1];latmin=bbox[2,1];lonmax=bbox[1,2];latmax=bbox[2,2]
  g = graphFromEdges(
         graphEdgesFromBase(lonmin,latmin,lonmax,latmax),
         densraster,from_query = FALSE)
  return(g)
}







