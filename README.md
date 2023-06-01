# network_analysis

16/07/2021
Louis Manière
ERC Desert Networks - CNRS
--------------

Metadata file :
DicoGIS_network_analysis_2021-08-20.xls

--------------
neighbour network analysis

Folder structure : 	

+---input

|   +---networks = camel networks for roman and ptolemaic period with and without watering places (nowp : without watering place, nonbuilt = non built settlement)

|   \---sites = archaeological sites and watering places used to build the camel networks (nowp : without watering place, nonbuilt = non built settlement)

+---output

|   +---neighbour_network_40km = = ptolemaic and roman neighbour networks output from neigbour_network processing (nowp : without watering place, nonbuilt = non built settlement)

|   +---neighbour_network_analysis_sfnetworks = spatial network nodes and edges with centrality measures and clustering compute wit sp and igraph package from processing/graph_analysis/neighbour_network_graph_analysis_sfnetworks.R

|   \---neighbour_network_analysis_sp_igraph = network nodes with centrality measures and clustering compute wit sp and igraph package from processing/graph_analysis/neighbour_network_graph_analysis_sp_igraph.R

\---processing

    +---graph_analysis = compute neighbour network by distance threshold (pyQGIS 3 python script)
    
    \---neighbour_network = extract nodes from shapefile network and compute node indicators (R scripts)

Fields name description :
Centrality measures
documentation : https://dshizuka.github.io/networkanalysis/04_measuring.html
degree = degree centrality
clsnes = closeness centrality
btwnes = betweenness centrality
eigen = eigenvector centrality

community detection
Documentation = https://dshizuka.github.io/networkanalysis/05_community.html
btwnes_c = betweenness clustering
louv_c = Louvain clustering
greed_c = fastgreedy clustering
walk_c = walktrap clustering
