# Author : Louis Manière
# CNRS HiSoMA UMR 5189
# Desert Networks ERC project
# Date : 16/07/2021

# network analysis on antic camel neighbour networks to extract node centrality measures and clustering with sfnetworks

# documentation : https://luukvdmeer.github.io/sfnetworks/index.html
# sessionInfo and packages in neighbour_network_graph_analysis_sfnetworks_sessionInfo_packages.txt

# installation packages
# if(!"remotes" %in% installed.packages()) {
#   install.packages("remotes")
# }
# 
# cran_pkgs = c(
#   "sfnetworks",
#   "sf",
#   "tidygraph",
#   "tidyverse",
#   "igraph"
# )
# 
# remotes::install_cran(cran_pkgs)

# packages
library(sfnetworks)
library(sf)
library(tidygraph)
library(tidyverse)
library(igraph)

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
path_graph_nodes <- paste(wd_folder,'/output/neighbour_network_analysis_sfnetworks/' ,network_name ,'_node_sfnetworks' ,extension, sep='')
# graph edges shapefile
path_graph_edges <- paste(wd_folder,'/output/neighbour_network_analysis_sfnetworks/' ,network_name ,'_edge_sfnetworks' ,extension, sep='')

# read input files
network <- st_read(path_spatial_network)
sites <- st_read(path_sites)
ggplot(data = network) + geom_sf()

# Give each edge a unique index
edges <- network %>%
  mutate(edgeID = c(1:n()))

# Create nodes at the start and end point of each edge
nodes <- edges %>%
  st_coordinates() %>%
  as_tibble() %>%
  rename(edgeID = L1) %>%
  group_by(edgeID) %>%
  slice(c(1, n())) %>%
  ungroup() %>%
  mutate(start_end = rep(c('start', 'end'), times = n()/2))

# Give each node a unique index
nodes <- nodes %>%
  mutate(xy = paste(.$X, .$Y)) %>% 
  mutate(nodeID = group_indices(., factor(xy, levels = unique(xy)))) %>%
  select(-xy)

# Combine the node indices with the edges
source_nodes <- nodes %>%
  filter(start_end == 'start') %>%
  pull(nodeID)

target_nodes <- nodes %>%
  filter(start_end == 'end') %>%
  pull(nodeID)

edges = edges %>%
  mutate(from = source_nodes, to = target_nodes)

# Remove duplicate nodes
nodes <- nodes %>%
  distinct(nodeID, .keep_all = TRUE) %>%
  select(-c(edgeID, start_end)) %>%
  st_as_sf(coords = c('X', 'Y')) %>%
  st_set_crs(st_crs(edges))

# create sfnetwork object
net <- sfnetwork(nodes, edges,node_key = "nodeID", directed = FALSE)

# length and centrality measure
net <- net %>%
  activate("edges") %>%
  mutate(weight = edge_length()) %>%
  mutate(btwnes = centrality_edge_betweenness(weights = weight)) %>%
  activate("nodes") %>%
  mutate(degree = centrality_degree()) %>%
  mutate(btwnes = centrality_betweenness(weights = weight, directed = FALSE))

# cluster
# edges
net = net %>%
  morph(to_linegraph) %>%
  mutate(btwnes_c = group_edge_betweenness(weights = NULL, directed = FALSE)) %>%
  unmorph()

# nodes
net <- net %>%
  activate ("nodes") %>%
  mutate(btwnes_c = group_edge_betweenness(weights = NULL, directed = FALSE)) %>%
  mutate(louv_c = group_louvain(weights = NULL)) %>%
  mutate(greed_c = group_fast_greedy(weights = NULL))

# plot
autoplot(net)

# sf extraction
# nodes
nodes_sf <- net %>%
  activate("nodes") %>%
  st_as_sf()
# edges
edges_sf <- st_as_sf(net, "edges")

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
# edges
st_write(edges_sf, path_graph_edges, delete_layer = TRUE)
