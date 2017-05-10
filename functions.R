
library(RPostgreSQL)
library(rgeos)
library(rgdal)
library(raster)
library(igraph)


# functions

defaultParams<-function(){
  return(list(
    
  ))
}


#'
#' @description random network of the given type and number of nodes
generateNetwork<-function(type,n,params=defaultParams()){
  
  if(type=="lattice"){
    return(make_lattice(length=floor(sqrt(n)),dim=2,directed = FALSE))
  }
  
  if(type=="pa-age"){
    return(sample_pa_age(n,m=10, pa.exp=1, aging.exp=-3, aging.bin=1000,directed = FALSE))
  }
  
}

normalizedBetweenness<-function(g){
  return(mean(betweenness(g)*2/(vcount(g)*(vcount(g)-1))))
}

efficiency<-function(g){
  distmat = distances(g)
  diag(distmat)<-Inf
  return(mean(1/distmat))
}


deltaMeasure <- function(g,removed_prop,func){
  res=c()
  for(alpha in removed_prop){
    subgraph = subgraph.edges(g,sample.int(n = ecount(g),size = (1-alpha)*ecount(g),replace = FALSE),delete.vertices = F)
    res=append(res,func(subgraph)-func(g))
  }
  return(res)
}


#'
#' @description each function in measure should return a named list
estimateRandomMeasures <- function(type,n,measures,nbootstrap){
  vals = list()
  for(b in 1:nbootstrap){
    if(b%%100==0){show(b)}
    g = generateNetwork(type,n)
    for(measure in measures){
      val = measure(g)
      
    }
  }
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
         graphEdgesFromBase(lonmin,latmin,lonmax,latmax,dbname=global.nwdb),
         densraster,from_query = FALSE)
  return(g)
}







