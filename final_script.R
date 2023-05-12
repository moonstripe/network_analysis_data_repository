setwd("/Users/kojinglick/Documents/School/Middlebury/2023 Spring/DPPG 8565/final")

rm(list=ls())
gc()
library(statnet)

ac_data <- read.csv("data/edgelist_ac.csv", header = TRUE)

ac_network <- network(ac_data, directed=TRUE, loops=TRUE, multiple=TRUE)

# Read the mapping of vertex IDs to names from a CSV file
id_to_name_df <- read.csv("data/name_map.csv", header = TRUE)

# Convert the data frame to a named vector for easy mapping
id_to_name_map <- setNames(id_to_name_df$title, id_to_name_df$id)

name_to_id_map <- setNames(id_to_name_df$id, id_to_name_df$title)

network.vertex.names(ac_network) <- id_to_name_df$title

ac_network %v% "vertex.names"

ac_network


??sna::centralization
??sna::betweenness

network.vertex.names(ac_network)[124]

sort(sna::betweenness(ac_network, cmode="undirected", tmaxdev=F), decreasing=TRUE)[2]

sort(sna::betweenness(ac_network, cmode="undirected", tmaxdev=F), decreasing=TRUE)[2]

which(sna::betweenness(ac_network, cmode="undirected", tmaxdev=FALSE) == sort(sna::betweenness(ac_network, cmode="undirected", tmaxdev=FALSE), decreasing=TRUE)[1])
# centralization in statnet
max(sna::betweenness(ac_network, cmode="undirected", tmaxdev=FALSE))
sna::centralization(ac_network, sna::betweenness, normalize=FALSE, mode="graph", diag=TRUE)

components_weak <- component.dist(ac_network, connected="weak")
components_weak

# cliques
comemb <- clique.census(ac_network, clique.comembership="sum")$clique.comemb # Co-membership matrix

CoNet <- network(comemb, matrix.type="adjacency",
                 directed=FALSE,
                 ignore.eval=FALSE,   # Keep edge weights
                 names.eval="weight") # names.eval stores edge weights
#   as edge attribute named "weight".

CoWeights <- as.sociomatrix(CoNet,"weight")


# k-cores
kco <- kcores(ac_network)

gplot(ac_network,
      vertex.col=kco,
      usearrows=FALSE,
      displaylabels=FALSE)

ac_network[kco>8, kco>8]

gplot(ac_network[kco>8, kco>8],
      usearrows=FALSE,
      displaylabels=TRUE)

ac_network[kco>8, kco>8]

?cug.test

cugBetSize <- cug.test(ac_network,
                       sna::centralization,
                       FUN.arg=list(FUN=sna::betweenness),
                       reps=10,
                       mode="graph", 
                       cmode="size",
                       diag = FALSE)

cugBetSize
cugBetEdges <- cug.test(ac_network,
                        sna::centralization,
                        FUN.arg=list(FUN=sna::betweenness, mode="graph", diag=FALSE), 
                        mode="graph", 
                        cmode="edges",
                        diag = FALSE)

# # Girvan-Newman

library(intergraph)
library(igraph)

ac_igraph <- asIgraph(ac_network)

ac_undirected <- simplify(
  as.undirected(ac_igraph, mode="collapse"),
  remove.loops = TRUE
)


# ac_articulation
ac_art <- ac_undirected

V(ac_art)$names <- V(ac_art)$vertex.names

ac_art

get_node_name <- function(graph, node_id) {
  node_name <- as.list(V(graph)$vertex.name)[node_id]
  return(node_name)
}

ac_as.list <- as.list(articulation_points(ac_art))

articulation_points(ac_art)

articulation_name_list <- vector()

for (i in articulation_points(ac_art)) {
  node_name <- as.list(V(ac_art)$names)[i]  
  articulation_name_list[length(articulation_name_list)+1] <- node_name
}

articulation_df <- as.data.frame(cbind(articulation_name_list))
articulation_df
articulation_df <- apply(articulation_df,2,as.character)
write.csv(articulation_df, "data/articulation_name_list.csv")

# plot
V(ac_art)$color = ifelse(V(ac_art) %in%
                                  articulation_points(ac_art),   
                                "salmon", "lightblue")

V(ac_art)$size = ifelse(V(ac_art) %in%
                                  articulation_points(ac_art),   
                                2, 0.1)

V(ac_art)$label = ifelse(V(ac_art) %in% articulation_points(ac_art),
                                V(ac_art)$names, "")


# Idk about bridges, not super useful

E(ac_art)$size = ifelse(E(ac_art) %in% bridges(ac_art),
                               0.1, 5) 

E(ac_art)$color = ifelse(E(ac_art) %in% bridges(ac_art),
                                "lightblue", "salmon") 

plot(ac_art, layout=layout_with_kk(ac_art))

# Plot largest Bi-Component: Only one bi-component, not really helpful

ac_bi <- ac_undirected

biconnected_components(ac_bi)

largest_component <- lapply(biconnected_components(ac_bi)$components, 
                            length) %>% which.max()

for (i in biconnected_components(ac_bi)$components) {
 if (length(i) > 2) {
   print(i)
 }
}

largest_component

V(ac_bi)$color <- ifelse(V(ac_bi) %in%
                       biconnected_components(ac_bi)$components[[largest_component]],
                     "salmon","lightblue")

plot(ac_bi)

# reach centrality: largest reach = AntifashGordon Lawsuit

ac_reach <- ac_undirected

reach2<-function(x){
  r=vector(length=vcount(x))
  for (i in 1:vcount(x)){
    n=neighborhood(x,2,nodes=i)
    ni=unlist(n)
    l=length(ni)
    r[i]=(l)}
  r}

# Function for 3-step reach
reach3<-function(x){
  r=vector(length=vcount(x))
  for (i in 1:vcount(x)){
    n=neighborhood(x,3,nodes=i)
    ni=unlist(n)
    l=length(ni)
    r[i]=(l)}
  r}


Reach_2 <- reach2(ac_reach)    # Note the differences between the object
Reach_3 <- reach3(ac_reach)    #   names and the function names!

for (i in Reach_2) {
  if (i>1100) {
    print(which(Reach_2 == i))
    index <- which(Reach_2 == i)
    print(as.list(V(ac_reach)$vertex.names)[index])
  }
}

# Edge Betweenness Highest Edge Betweenness edges involve AntiFash Gordon Lawsuit
ac_eb <- ac_undirected

edge_betweenness(ac_eb)

E(ac_eb)$width <- edge_betweenness(ac_eb)

# ugh too much data
plot.igraph(ac_eb,
            edge.width = igraph::edge.betweenness(ac_eb)+1,              # The "+1" was added to make edgewidths non-zero.
            edge.color = heat.colors(igraph::edge.betweenness(ac_eb)+1), # The "+1" was added to make edgewidths non-zero.
            vertex.shape="sphere",  # Here, we are using sphere because it looks cool.
            vertex.size=20,
            vertex.label.font=2,    # Here, we are using bold font.
            vertex.color="lightgreen")

Edge.Betweenness <- as.data.frame(edge_betweenness(ac_eb, directed = FALSE))
el <- get.edgelist(ac_eb)
el <- as.data.frame(el)
el$EB <- Edge.Betweenness

V1Names <- list()
V2Names <- list()

for (i in el$V1) {
  V1Names[length(V1Names)+1] <- as.list(V(ac_eb)$vertex.names)[i]
}

for (i in el$V2) {
  V2Names[length(V2Names)+1] <- as.list(V(ac_eb)$vertex.names)[i]
}

el$V1Names <- V1Names
el$V2Names <- V2Names
el
class(el)

df <- apply(el,2,as.character)
write.csv(df, "data/edge_betweenness.csv") 

# Burt's Constraint

ac_constraint <- ac_undirected

const <- constraint(ac_constraint)
invConstraint <- 1.125 - const  # (Inverse constraint = brokerage potential)
invConstraint

name <- list()
loop_index = 1
for (i in invConstraint) {
  
  name[length(name)+1] <- as.list(V(ac_constraint)$vertex.names)[loop_index]
  
  loop_index <- loop_index + 1
}

const_df <- as.data.frame(cbind(name, invConstraint))
const_df <- apply(const_df,2,as.character)

write.csv(const_df, "data/inverse_constraint.csv")

# centralization (Betweenness, Closeness)

ac_centr_betw <- ac_undirected

bw_centr <- centr_betw(ac_centr_betw, directed=F, normalized=F) 

bw_centr

bw_centr_tmax <- centr_betw_tmax(ac_undirected, directed=F)

bw_centr_tmax

tmax <- (1151-1)*(1151-2)/2

bw_centr$centralization/tmax



max(betweenness(ac_undirected))

cl_centr <- centr_clo(ac_centr_betw)

cl_centr

bw

bw_centr$centralization

centr_bw_list <- list()
loop_index = 1
for (i in bw_centr$res) {
  
  centr_bw_list[length(centr_bw_list)+1] <- as.list(bw_centr$res)[loop_index]
  
  loop_index <- loop_index + 1
}

centr_cl_list <- list()
loop_index = 1
for (i in cl_centr$res) {
  
  centr_cl_list[length(centr_cl_list)+1] <- as.list(cl_centr$res)[loop_index]
  
  loop_index <- loop_index + 1
}
#closeness csv
centr_cl_df <- as.data.frame(cbind(name, centr_cl_list))

centr_cl_df <-  apply(centr_cl_df, 2, as.character)

write.csv(centr_cl_df, "data/closeness_centralization_df.csv")


#btw csv
centr_bw_df <- as.data.frame(cbind(name, centr_bw_list))

centr_bw_df <-  apply(centr_bw_df, 2, as.character)

write.csv(centr_bw_df, "data/betweenness_centralization_df.csv")

# sd centralization (Betweenness)

sd(betweenness(ac_centr_betw))

sd(eigen_centrality(ac_centr_betw)$vector)

# Brokerage roles (help me god)

# split the network into clusters

ac_cliques <-  network(ac_data, directed=FALSE, loops=TRUE, multiple=TRUE)

# Read the mapping of vertex IDs to names from a CSV file
id_to_name_df <- read.csv("data/name_map.csv", header = TRUE)

# Convert the data frame to a named vector for easy mapping
id_to_name_map <- setNames(id_to_name_df$title, id_to_name_df$id)

name_to_id_map <- setNames(id_to_name_df$id, id_to_name_df$title)

network.vertex.names(ac_cliques) <- id_to_name_df$title

comemb <- clique.census(ac_cliques, clique.comembership="sum")$clique.comemb # Co-membership matrix

coords <- cmdscale(1/(1+comemb))        # Perform an MDS  (multidimensional scaling) 
# the more overlap, 
# the "closer" the vertices

CoNet <- network(comemb, matrix.type="adjacency", 
                 directed=FALSE,
                 ignore.eval=FALSE,   # Keep edge weights
                 names.eval="weight") # names.eval stores edge weights
#   as edge attribute named "weight".

CoNet

CoWeights <- as.sociomatrix(CoNet,"weight") 

gplot(CoNet,coord=coords, 
      usearrows=FALSE, 
      displaylabels=TRUE,
      edge.lwd=CoNet%e%'weight'*2) # Use the edge attribute, "weight" (multiplied by 2

# cluster dendrogram

hc <- hclust(as.dist(1/(CoWeights+1)))      # Hierarchal clustering by co-membership

plot(hc)                                    # Plot the dendrogram

# Optional: 
rect.hclust(hc, h=0.6)                      # Plot a cutoff point

