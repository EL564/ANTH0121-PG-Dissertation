################################################################################
#
#
#(5) PT_network_transformations
#session includes network descptives and transformation (combining edge-level and node-level data)
#author: Ella Lipscombe 
#
#
################################################################################

#start session


#set working directory
setwd("D:/MSc HEB 2023-2025/Modules 2024-2025/ANTH0121 Dissertation/Methodology/R script")

#shortcuts:
PT_edge_list <- read.csv ("PT_edge_list.csv", header = TRUE) #code to import master data set
PT_ego_list <- read.csv ("PT_ego_list.csv", header = TRUE) #code to import master data set

load("net_Q_NBT.RData") #not clean Q_NBT combined network (i.e., combined edge and node attributes)
load("net_Q_NBT_clean.RData") #this is the cleaned up version post-checks (outliers and missing data removed)

#load packages 
install.packages(c("network", "sna", "ergm", "statnet", "dplyr"))
library(network)
library(sna)
library(ergm)
library(statnet)  
library(dplyr)
library(ggplot2)

################################################################################

#Q network transformation  ####

#treat the network as binary, regardless of multiple ties from the same ego (i.e., strength)
Q_support <- PT_edge_list %>%
  filter(Support_type == "Q") %>%
  filter(Ego != Alter) %>%          # remove loops here
  distinct(Ego, Alter, .keep_all = TRUE)

#network construction
net_Q_support <- network(Q_support, directed = TRUE, matrix.type = "edgelist", loops = TRUE)

#make sure the order matches. We assume 'Ego' in PT_ego_list corresponds to node labels in the network
for(attr_name in names(PT_ego_list)[names(PT_ego_list) != "Ego"]) {
  attr_values <- setNames(PT_ego_list[[attr_name]], PT_ego_list$Ego)
  # Reorder according to node names
  ordered_attr <- attr_values[network.vertex.names(net_Q_support)]
  set.vertex.attribute(net_Q_support, attr_name, ordered_attr)
}

#review network (ego and edge attributes)
list.vertex.attributes(net_Q_support)
get.vertex.attribute(net_Q_support, "Age")
list.edge.attributes(net_Q_support)
library(network)

#get vector of vertex attribute "Age"
ages <- get.vertex.attribute(net_Q_support, "Age")

#find vertices with NA in Age
vertices_to_remove <- which(is.na(ages))

#remove those vertices from the network
net_Q_support_clean <- delete.vertices(net_Q_support, vertices_to_remove)
net_Q_support <- net_Q_support_clean

#save economic network
save(net_Q_support, file = "net_Q_support.RData")
load("net_Q_support.RData")


#K network transformation  ####

#treat the network as binary, regardless of multiple ties from the same ego (i.e., strength)
K_support <- PT_edge_list %>%
  filter(Support_type == "K") %>%
  distinct(Ego, Alter, .keep_all = TRUE)

#network construction
net_K_support <- network(K_support, directed = TRUE, matrix.type = "edgelist", loops = TRUE)

#make sure the order matches. We assume 'Ego' in PT_ego_list corresponds to node labels in the network
for(attr_name in names(PT_ego_list)[names(PT_ego_list) != "Ego"]) {
  attr_values <- setNames(PT_ego_list[[attr_name]], PT_ego_list$Ego)
  # Reorder according to node names
  ordered_attr <- attr_values[network.vertex.names(net_K_support)]
  set.vertex.attribute(net_K_support, attr_name, ordered_attr)
}

#review network (ego and edge attributes)
list.vertex.attributes(net_K_support)
list.edge.attributes(net_K_support)

#save economic network
save(net_K_support, file = "net_K_support.RData")
load("net_K_support.RData")

#Q: NBT transformation ####

#treat the network as binary, regardless of multiple ties from the same ego (i.e., strength)
Q_NBT_edges_clean <- PT_edge_list %>%
  filter(Support_type == "Q", NBTs == "YES") %>%
  distinct(Ego, Alter, .keep_all = TRUE)

#network construction
net_Q_NBT <- network(Q_NBT_edges_clean, directed = TRUE, matrix.type = "edgelist")

#now list edge attributes
list.edge.attributes(net_Q_NBT)
get.edge.attribute(net_Q_NBT, "Dyadic_distance")       

#add node (vertex) attributes 
head(PT_ego_list)
str(PT_ego_list)
network.vertex.names(net_Q_NBT)

#make sure the order matches. We assume 'Ego' in PT_ego_list corresponds to node labels in the network
for(attr_name in names(PT_ego_list)[names(PT_ego_list) != "Ego"]) {
  attr_values <- setNames(PT_ego_list[[attr_name]], PT_ego_list$Ego)
  # Reorder according to node names
  ordered_attr <- attr_values[network.vertex.names(net_Q_NBT)]
  set.vertex.attribute(net_Q_NBT, attr_name, ordered_attr)
}

#review network (ego and edge attributes)
list.vertex.attributes(net_Q_NBT)
get.vertex.attribute(net_Q_NBT, "Kin_density_QNBT")
list.edge.attributes(net_Q_NBT)

#save NBT network
save(net_Q_NBT, file = "net_Q_NBT.RData")
load("net_Q_NBT.RData")

#descriptives 

#network density 
193/1095*100

#in-degree 
Q_NBT_in_degree <- Q_NBTs[which(Q_NBTs$Alter %in% Q_NBTs$Ego), ] 
Q_NBT_in_degree <- as.data.frame(table(Q_NBTs$Alter))
colnames(Q_NBT_in_degree) <- c("Ego_ID", "in_degree")
mean(Q_NBT_in_degree$in_degree) #mean in-degree tie counts among those who received at least one in-degree nomination

hist(Q_NBT_in_degree$in_degree)
shapiro.test(Q_NBT_in_degree$in_degree)
median(Q_NBT_in_degree$in_degree)
IQR(Q_NBT_in_degree$in_degree)

#out-degree
Q_NBT_out_degree <- as.data.frame(table(Q_NBTs$Ego))  # or use Q_NBT_out_degree$Ego if you filtered above
colnames(Q_NBT_out_degree) <- c("Ego_ID", "out_degree")
mean(Q_NBT_out_degree$out_degree) # mean of out-degrees where at least one tie has been sent 
sd(Q_NBT_out_degree$out_degree)

median(Q_NBT_out_degree$out_degree)
IQR(Q_NBT_out_degree$out_degree)

#density
gden(net_Q_NBT) #very sparce representation of the whole dyadic network. Less than 1.8% meaning only 1.8% of all possible directed ties are present in this network

#number of nodes and edges
network.size(net_Q_NBT)
network.edgecount(net_Q_NBT)
sum(degree(net_Q_NBT, cmode = "outdegree") > 0) #number of unique nodes 

#unique egos and alter
unique_egos_NBT <- length(unique(Q_NBT_edges_unique$Ego))
unique_alters_NBT <- length(unique(Q_NBT_edges_unique$Alter))

cat("Q_NBT_edges_unique:\n")
cat("Unique Egos:", unique_egos_NBT, "\n")
cat("Unique Alters:", unique_alters_NBT, "\n\n")
57+81

#component analysis (if there are disconnected parts)
components(net_Q_NBT) #62 disconnected subgraphs (i.e., weakly connected components) suggesting the network of Q-based NBTs are highly fragmented, with many small groups

#kinship 
list.vertex.attributes(net_Q_NBT)

kin_density_values <- get.vertex.attribute(net_Q_NBT, "Kin_density_QNBT") #extract kin density attribute from the network
kin_density_values <- get.vertex.attribute(net_Q_NBT, "Kin_density_QNBT")

mean_kin_density <- mean(kin_density_values, na.rm = TRUE) #calculate the mean kin density (excluding any NAs)
mean_kin_density
sd(kin_density_values, na.rm =TRUE)
kin_density_values

#calculate the mean kin density (excluding any NAs)
mean_kin_density <- mean(kin_density_values, na.rm = TRUE)
mean_kin_density
sd(kin_density_values, na.rm=TRUE)

#distances
list.edge.attributes(net_Q_NBT)
dyadic_distances <- get.edge.attribute(net_Q_NBT, "Dyadic_distance") #extract the Dyadic_distance edge attribute
mean_dyadic_distance <- mean(dyadic_distances, na.rm = TRUE)
mean_dyadic_distance
sd(dyadic_distances, na.rm = TRUE)

#simple and final NBT network 
edges_NBT_clean <- as.data.frame.network(net_Q_NBT_clean, matrix.type = "edgelist")

#create network object from edges_NBT dyad list
net_Q_NBT_clean <- network(edges_NBT_clean[, c(".tail", ".head")], 
                           directed = TRUE,  # Change to FALSE if undirected
                           matrix.type = "edgelist")

##revised
library(network)
library(intergraph)  # for conversion functions
library(igraph)
library(scales)      # for rescale()

#convert to igraph
g <- asIgraph(net_Q_NBT_clean)

#extract attributes
age_vec    <- V(g)$Age
gender_vec <- V(g)$Gender

#scale node size based on age
node_sizes <- rescale(age_vec, to = c(5, 20))  # adjust size range

#shape by gender
node_shapes <- ifelse(gender_vec == "M", "circle", "square")

#identify isolates
isolate_nodes <- degree(g, mode = "all") == 0

#colors for isolates and non-isolates
vertex_fill <- ifelse(
  isolate_nodes,
  adjustcolor("olivedrab1", alpha.f = 0.3),  # alpha 0.3 = ~70% transparent
  "olivedrab2"
)

#borders: isolates get olivedrab1, others olivedrab2
vertex_border <- ifelse(isolate_nodes, "olivedrab1", "olivedrab2")

#use Fruchterman–Reingold layout
layout_fr <- layout_with_fr(g)

#more spaced-out Fruchterman–Reingold layout
layout_fr <- layout_with_fr(
  g,
  niter = 2000,        # more iterations for better separation
  repulserad = vcount(g)^3  # bigger = more space
)

#scale coordinates to stretch layout
layout_fr_default <- layout_with_fr(g)

#plot with more space
jpeg("net_Q_NBT_clean_plot_isolates.jpeg", width = 8, height = 8, units = "in", res = 300)
plot(
  g,
  layout             = layout_fr_default,
  vertex.size        = node_sizes,
  vertex.color       = vertex_fill,
  vertex.frame.color = "black",
  vertex.shape       = node_shapes,
  edge.color         = "red",
  edge.arrow.size    = 0.4,
  edge.width         = 1.5,
  vertex.label       = NA,
  margin             = 0
)
dev.off()

# Reset plot settings
par(mfrow = c(1, 1))


################################################################################

#end session 

################################################################################
