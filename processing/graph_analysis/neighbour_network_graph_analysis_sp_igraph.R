# Author : Louis Manière
# CNRS HiSoMA UMR 5189
# Desert Networks ERC project
# Date : 16/07/2021

# network analysis on antic camel neighbour networks to extract node centrality measures and clustering

# sessionInfo and packages in neighbour_network_graph_analysis_sp_igraph_sessionInfo_packages.txt

# installation packages
# if(!"remotes" %in% installed.packages()) {
#   install.packages("remotes")
# }
# 
# cran_pkgs = c(
#   "sp",
#   "rgdal",
#   "shp2graph",
#   "igraph",
#   "sf"
# )
# 
# remotes::install_cran(cran_pkgs)

# Packages
library(sp)
library(rgdal)
library(shp2graph)
library(igraph)
library(sf)

# path to files

# network file name
network_name <- 'roman_neighbour_network_40km_nonbuilt'
site_name <- 'roman_sites_nonbuilt'
extension <- '.shp'

# working folder
wd_folder <- file.path(dirname(dirname(getwd())))

# input
# spatial network path
path_spatial_network <- paste(wd_folder, '/output/neighbour_network_40km/',network_name ,extension, sep='')
path_sites <- paste(wd_folder, '/input/sites/',site_name ,extension, sep='')

# output
# graph nodes shapefile
path_graph_nodes <- paste(wd_folder,'/output/neighbour_network_analysis_sp_igraph/' ,network_name ,'_node_sp_igraph' ,extension, sep='')

# read input files
network <- readOGR(dsn = path_spatial_network)
sites <- st_read(path_sites)

# check network connection
connect <- nt.connect(network)

# read network with shp2graph package (see documentation for output value)
network_read <- readshpnw (network, ELComputed=TRUE, longlat=FALSE, Detailed=FALSE,ea.prop=NULL)

# transform shp2graph network to igraph graph
graph_network <- nel2igraph(network_read[[2]], network_read[[3]], weight = NULL, eadf = NULL, Directed = FALSE)

# plot
plot(graph_network, 
     vertex.size=4,
     vertex.label.dist=0.5,
     vertex.color="green",
     edge.arrow.size=0.5)

# clustering coefficient
# local coefficient
local_cluster_coef =transitivity(graph_network,type = "local")

# global coefficient
average_cluster_coef = transitivity(graph_network,type = "average")

# centrality indicators
# degre
centrality_degree<- degree(graph = graph_network, 
                           v = V(graph_network), 
                           mode = "all", 
                           loops = TRUE, 
                           normalized = FALSE)
# closeness
centrality_closeness <- closeness(graph = graph_network, 
                                  vids = V(graph_network), 
                                  mode = c("all"), 
                                  weights = NULL, 
                                  normalized = FALSE)
# betweenness 
centrality_betweenness <- betweenness(graph = graph_network, 
                                      v = V(graph_network), 
                                      directed = FALSE, 
                                      weights = NULL, 
                                      nobigint = TRUE,
                                      normalized = FALSE)
# eigen value
centrality_eigen <- eigen_centrality(graph = graph_network, 
                                     directed = FALSE, 
                                     scale = TRUE,
                                     weights = NULL, 
                                     options = arpack_defaults)

# edge betweenness community structure
edge_betweenness <- cluster_edge_betweenness(graph_network, weights = E(graph_network)$weight,
                                               directed = FALSE, edge.betweenness = TRUE, merges = TRUE,
                                               bridges = TRUE, modularity = TRUE, membership = TRUE)
plot(edge_betweenness,graph_network)

V(graph_network)$betweenness_cluster <- membership(edge_betweenness)

# Louvain community structure
louvain <- cluster_louvain(graph_network, weights = NULL)
plot(louvain,graph_network)
V(graph_network)$louvain_cluster <- membership(louvain)

# Clauset and Newman aggregative greedy community structure
greedy <- cluster_fast_greedy(graph_network, merges = TRUE, modularity = TRUE,
                              membership = TRUE, weights = E(graph_network)$weight)
plot(greedy,graph_network)
V(graph_network)$greedy_cluster <- membership(greedy)

# walktrap community structure
walktrap <- cluster_walktrap(graph_network, weights = E(graph_network)$weight, steps = 4,
                 merges = TRUE, modularity = TRUE, membership = TRUE)
plot(walktrap,graph_network)
V(graph_network)$walktrap_cluster <- membership(walktrap)


# vertex/node data in dataframe
df_graph_network <- data.frame(V = as.vector(V(graph_network)), 
                               x=vertex_attr(graph = graph_network, name = "x", index = V(graph_network)),
                               y=vertex_attr(graph = graph_network, name = "y", index = V(graph_network)),
                               degree = centrality_degree,
                               clsnes = centrality_closeness,
                               btwnes = centrality_betweenness,
                               eigen = centrality_eigen$vector,
                               btwnes_c = V(graph_network)$betweenness_cluster,
                               louv_c = V(graph_network)$louvain_cluster,
                               greed_c = V(graph_network)$greedy_cluster,
                               walk_c = V(graph_network)$walktrap_cluster
                               )

# extract xy coordinate in separate dataframe
xy <- data.frame(x=df_graph_network$x, y=df_graph_network$y)
# convert dataframe to sp SpatialPointsDataFrame
spdf_graph_nodes <- SpatialPointsDataFrame(coords = xy, 
                                            data = df_graph_network, 
                                            coords.nrs = numeric(0), 
                                            proj4string = CRS("+proj=utm +zone=36 +datum=WGS84 +units=m +no_defs "), 
                                            match.ID = TRUE, 
                                            bbox = NULL)

# sf conversion
nodes_sf <- st_as_sf(spdf_graph_nodes)

# spatial join to site layer
nodes_sf <- st_join(
  nodes_sf,
  sites,
  join = st_is_within_distance,
  dist = 0.001
)

# write
# nodes
st_write(nodes_sf, path_graph_nodes, delete_layer = TRUE)
