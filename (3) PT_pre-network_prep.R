################################################################################
#
#
#(3) PT_pre-network_prep
#session includes: subsets of edge data for economic (Q), behavioural (K), and needs-based netwokr (Q: NBT)
#author: Ella Lipscombe 
#
#
################################################################################


#start session 

#packages 
library(igraph)
library(dplyr)
library(tidyr)
library(sf)
library(ggplot2)
library(ggraph)
library(tidygraph)
library(RColorBrewer)
library(ggtext)  
library(grid)    

#shortcut:
setwd("D:/MSc HEB 2023-2025/Modules 2024-2025/ANTH0121 Dissertation/Methodology/R script")

SU_combined  <- read.csv ("SU_combined.csv", header = TRUE)
PT_dyads_kin_pedigree <- read.csv ("PT_dyads_kin_pedigree.csv", header = TRUE) #code to import master data set
PT_indiv <- read.csv("PT_Indiv.csv", header=TRUE) #raw data set containing individual-level data and other ENDOW variables
SU_distances <- read.csv ("su_distances.csv", header = TRUE) #raw dataset of distances between sharing unit
PT_full_dyads<- read.csv ("PT_full_dyads.csv", header = TRUE) #code to import master data set

#final dyadic dataset with kin and distances merged
PT_dyads_kin_distance <- read.csv ("PT_dyads_kin_distance.csv", header = TRUE) #code to import master data set

#social network subset shortcuts 
PT_dyads_Q <- read.csv ("PT_dyads_Q.csv", header = TRUE) 
PT_dyads_K <- read.csv ("PT_dyads_K.csv", header = TRUE) 
Q_NBTs <- read.csv ("Q_NBTs.csv", header = TRUE) 

################################################################################

#data prep

#merging edge data distance to the full dyadic data set ####


# 1. standardise column names and types in SU_distances
SU_distances <- SU_distances %>%
  rename(Ego_SU = sui, Alter_SU = suj) %>%
  mutate(across(c(Ego_SU, Alter_SU), as.character))

# 2. standardise case in both datasets
PT_dyads_kin_pedigree <- PT_dyads_kin_pedigree %>%
  mutate(across(c(Ego_SU, Alter_SU), ~ toupper(trimws(.))))

SU_distances <- SU_distances %>%
  mutate(across(c(Ego_SU, Alter_SU), ~ toupper(trimws(.))))

# 3. define helper function to generate undirected pair keys
make_pair_key <- function(x, y) {
  paste0(pmin(x, y), "_", pmax(x, y))
}

# 4. create pair_key in both datasets
PT_dyads_kin_pedigree <- PT_dyads_kin_pedigree %>%
  mutate(pair_key = make_pair_key(Ego_SU, Alter_SU))

SU_distances <- SU_distances %>%
  mutate(pair_key = make_pair_key(Ego_SU, Alter_SU))

# 5. join distance to dyads
SU_distances <- SU_distances %>%
  distinct(pair_key, .keep_all = TRUE)

PT_dyads_kin_pedigree <- PT_dyads_kin_pedigree %>%
  left_join(SU_distances %>% select(pair_key, distance), by = "pair_key")

# 6. drop pair_key if no longer needed
PT_dyads_kin_pedigree <- PT_dyads_kin_pedigree %>% select(-pair_key)

# 7. change column name to Distance 
PT_dyads_kin_pedigree <- PT_dyads_kin_pedigree %>%
  rename(Dyadic_distance = distance)

# 8. write as csv
write.csv(PT_dyads_kin_pedigree, "PT_dyads_kin_distance.csv", row.names = FALSE) #save cleaned master dataset
PT_dyads_kin_distance <- read.csv ("PT_dyads_kin_distance.csv", header = TRUE) #code to import master data set

#merge ego community status to the full dyadic data set ####

#rename variables and values  
names(SU_combined)[names(SU_combined) == "SharingUnitID"] <- "Ego_SU"
names(SU_combined)[names(SU_combined) == "DisplacementStatus"] <- "Community"

SU_combined$Community[SU_combined$Community == 1] <- 2
SU_combined$Community[SU_combined$Community == 0] <- 1

#join SU_combined to PT_dyads_kin_classified by Ego
PT_full_dyads <- merge(PT_dyads_kin_classified, SU_combined,
                   by = "Ego",
                   all.x = TRUE)

#create a new dataset by removing the last 5 columns
PT_full_dyads <- PT_full_dyads[, -( (ncol(PT_full_dyads)-4) : ncol(PT_full_dyads) )]

write.csv(PT_full_dyads, "PT_full_dyads.csv", row.names = FALSE) #save cleaned master dataset
PT_full_dyads<- read.csv ("PT_full_dyads.csv", header = TRUE) #code to import master data set

#subsets of edge-level data frames ####


#Q: in-degree ties = how many times the head of household appears in the alter column

class(PT_full_dyads$Community)
PT_full_dyads$Community <- as.factor(PT_full_dyads$Community) #change community variable to factor

#step 1. subset data into Q questions only (1-5)
PT_dyads_Q <- PT_dyads_kin_distance[which(PT_dyads_kin_distance$Question %in% c("1", "2", "3", "5")), ] #create new data object of subset 

write.csv(PT_dyads_Q, "PT_dyads_Q.csv", row.names = FALSE) #save cleaned master dataset
PT_dyads_Q <- read.csv ("PT_dyads_Q.csv", header = TRUE) #code to import master data set

#step 2. (sample-level) subset data into those that contain an ego match from SU_combined in the alter column 
Q_in_degree <- PT_dyads_Q[which(PT_dyads_Q$Alter %in% SU_combined$Ego), ]

write.csv(Q_in_degree, "Q_in_degree.csv", row.names = FALSE) #save cleaned master dataset
Q_in_degree <- read.csv ("Q_in_degree.csv", header = TRUE) #code to import master data set
#K: in-degree ties = how many times the head of household appears in the alter column


#step 1. subset data into Q questions only (6-9) 
PT_dyads_K <- PT_dyads_kin_distance[which(PT_dyads_kin_distance$Question %in% c("6", "7", "8", "9")), ] #create new data object of subset 

write.csv(PT_dyads_K, "PT_dyads_K.csv", row.names = FALSE) #save cleaned master dataset
PT_dyads_K <- read.csv ("PT_dyads_K.csv", header = TRUE) #code to import master data set

#step 2. (sample-level) subset data into those that contain an ego match from SU_combined in the alter column 
K_in_degree <- PT_dyads_K[which(PT_dyads_K$Alter %in% SU_combined$Ego), ]

write.csv(K_in_degree, "K_in_degree.csv", row.names = FALSE) #save cleaned master dataset
K_in_degree <- read.csv ("K_in_degree.csv", header = TRUE) #code to import master data set


#(Q): economic out-degree = the count of nominations from ego to all alters (not used in thesis)

#step 1. (sample-level) rename PT_dyads)Q as dyads show outward ties
Q_out_degree <- PT_dyads_Q

#step 2. 
write.csv(Q_out_degree, "Q_out_degree.csv", row.names = FALSE) #save cleaned master dataset
Q_out_degree <- read.csv ("Q_out_degree.csv", header = TRUE) #code to import master data set


#(K): caring out-degree = the count of nominations from ego to alter (not used in thesis)

#step 1. (sample-level) rename PT_dyads_K as dyads show outward ties 
K_out_degree <- PT_dyads_K

#step 2. 
write.csv(K_out_degree, "K_out_degree.csv", row.names = FALSE) #save cleaned master dataset
K_out_degree <- read.csv ("Q_out_degree.csv", header = TRUE) #code to import master data set

#Needs-based transfer subset
Q_NBTs <- PT_dyads_Q[which(PT_dyads_Q$Question %in% c("1", "3")), ] #create new data object of subset 
Q_NBTs <- Q_NBTs %>% distinct(Ego, Alter, .keep_all = TRUE) #unique (duplicates removed)

write.csv(Q_NBTs, "Q_NBTs.csv", row.names = FALSE) #save
Q_NBTs <- read.csv ("Q_NBTs.csv", header = TRUE) #import


################################################################################

#drawings

#Q drawings; K drawings ####

#Q (economic)

#change community to factor 
class(Q_in_degree$Community)
Q_in_degree$Community <- as.factor(Q_in_degree$Community)

#no. nodes:
num_nodes_Q <- vcount(Draw_Q_in_degree)
cat("Number of nodes in Draw_Q_in_degree:", num_nodes_Q, "\n")

#no. edges 
num_edges_Q <- ecount(Draw_Q_in_degree)
cat("Number of edges in Draw_Q_in_degree:", num_edges_Q, "\n")

#simple network drawing 
Draw_Q<- graph_from_data_frame(net_Q_support, directed = TRUE)

plot(Draw_Q,
     vertex.label = NA,
     vertex.size = 5,
     edge.arrow.size = 0.3,
     main = "Q-Type In-Degree Support Network")

#more complex social network drawing
q_tbl_graph <- as_tbl_graph(Draw_Q)

ggraph(q_tbl_graph, layout = "fr") +  # Or try "kk" for Kamada-Kawai layout
  geom_edge_link() +
  geom_node_point() +
  theme_void() +
  ggtitle("Q-Type Support Network")

#spatial social networks drawing
#1.match vertex names (Ego IDs) to coordinates in PT_full_dyads
coords_df <- PT_full_dyads[match(V(Draw_Q)$name, PT_full_dyads$Ego), c("Longitude", "Latitude")]

#2.convert to matrix for spatial layout
layout_coords <- as.matrix(coords_df)

#3.plot the network using spatial coordinates
plot(Draw_Q,
     layout = layout_coords,
     vertex.label = NA,
     vertex.size = 5,
     edge.arrow.size = 0.3,
     main = "Spatial Q-Type In-Degree Network")

#4.colour code

#ensure both are character for matching
V(Draw_Q)$name <- as.character(V(Draw_Q)$name)
PT_full_dyads$Ego <- as.character(PT_full_dyads$Ego)

#match and assign Community to vertices
V(Draw_Q)$Community <- PT_full_dyads$Community[match(V(Draw_Q)$name, PT_full_dyads$Ego)]

#assign colors: 1 = blue, 2 = orange (customisable)
vertex_colors <- ifelse(V(Draw_Q)$Community == 1, "mediumblue", "darkorange1")

pdf("SQ_plot.pdf", width = 7, height = 7)

plot(Draw_Q,
     layout = layout_coords,
     vertex.label = NA,
     vertex.size = 3,
     vertex.color = vertex_colors,
     edge.arrow.size = 0.5)

dev.off()

#5.convert igraph to tidygraph
q_tbl_graph <- as_tbl_graph(Draw_Q)

q_tbl_graph <- q_tbl_graph %>%
  mutate(Community = as.factor(Community))

#plot with community-based color
Q_in_SN <- ggraph(q_tbl_graph) +
  geom_edge_link(color = "gray70", alpha = 0.6) +
  geom_node_point(aes(color = Community), size = 4) +
  scale_color_manual(values = c("1" = "mediumblue", "2" = "darkorange1")) +
  theme_void() +
  ggtitle("Spatial Q-Type Support Network by Community") +
  theme(legend.position = "right")

Q_in_SN #name of social network 

#modifications. . .

#adding weights 
#call target nodes (alters) of each edge
to_ids <- as.character(ends(Draw_Q_in_degree, es = E(Draw_Q_in_degree), names = TRUE)[, 2])

#count nominations per alter
nomination_counts <- table(to_ids)

#convert to named numeric vector
nomination_vector <- as.numeric(nomination_counts)
names(nomination_vector) <- names(nomination_counts)

#assign weights to each edge based on target node
E(Draw_Q_in_degree)$weight <- nomination_vector[to_ids]

q_tbl_graph <- as_tbl_graph(Draw_Q_in_degree)

q_tbl_graph <- q_tbl_graph %>%
  activate(nodes) %>%
  mutate(Community = as.factor(Community)) #change community to factor

Q_in_SN <- ggraph(q_tbl_graph) +
  geom_edge_link(aes(width = weight), color = "gray70", alpha = 0.5, show.legend = FALSE) +  # Edge weights
  geom_node_point(aes(color = Community), size = 4, alpha = 0.8) +      # Node transparency
  scale_color_manual(
    values = c("1" = "mediumblue", "2" = "darkorange1"),
    name = "Community",
    labels = c("Community 1", "Community 2")
  ) +
  scale_edge_width(range = c(0.5, 2.5)) +  # Adjust edge width scaling
  theme_void() +
  ggtitle("Spatial Q-Type Support Network by Community") +
  theme(legend.position = "right")


Q_in_SN

#add node sized by age

#match age to each node by Ego ID
V(Draw_Q_in_degree)$Age <- PT_full_dyads$Age[match(V(Draw_Q_in_degree)$name, PT_full_dyads$Ego)]

q_tbl_graph <- as_tbl_graph(Draw_Q_in_degree)

q_tbl_graph <- q_tbl_graph %>%
  activate(nodes) %>%
  mutate(Community = as.factor(Community)) #change community to factor

Q_in_SN <- ggraph(q_tbl_graph) +
  geom_edge_link(aes(width = weight), color = "gray70", alpha = 0.5, show.legend = FALSE) +
  geom_node_point(aes(color = Community, size = Age), alpha = 0.6, show.legend = FALSE) +
  scale_color_manual(
    values = c("1" = "mediumblue", "2" = "darkorange1"),
    name = "Community",
    labels = c("Community 1", "Community 2")
  ) +
  scale_edge_width(range = c(0.5, 2.5)) +
  scale_size(range = c(2, 8), name = "Age") +
  theme_void() +
  ggtitle("Spatial Q-Type Support Network") +
  theme(legend.position = "right")

Q_in_SN

#add arrows 

Q_in_SN <- ggraph(q_tbl_graph)+
  geom_edge_link(
    aes(width = weight),
    color = "gray70",
    alpha = 0.4,
    show.legend = FALSE,
    arrow = arrow(length = unit(2, "mm"), type = "closed"),
    end_cap = circle(4, 'mm')
  ) +
  geom_node_point(aes(color = Community, size = Age), alpha = 0.8) +  # show legend for color
  scale_color_manual(
    values = c("1" = "mediumblue", "2" = "darkorange1"),
  ) +
  scale_edge_width(range = c(0.2, 1.2)) +
  scale_size(range = c(3, 10), guide = "none") +  # hide size (age) legend
  theme_void() +
  theme(legend.position = "none")

Q_in_SN

ggsave("Q_In-Degree_SN_plot.pdf", plot = Q_SN, width = 8, height = 6, units = "in")
ggsave("Q_SN_plot_landscape.pdf", plot = Q_SN, width = 11, height = 8, units = "in")

#K IN-DEGREE (SOCIAL NETWORK NO. 2)

#change community to factor 
class(K_in_degree$Community)
K_in_degree$Community <- as.factor(K_in_degree$Community)

#simple network drawing 
Draw_K_in_degree <- graph_from_data_frame(K_in_degree, directed = TRUE)

#no. nodes:
num_nodes_K <- vcount(Draw_K_in_degree)
cat("Number of nodes in Draw_K_in_degree:", num_nodes_Q, "\n")

#no. edges 
num_edges_Q <- ecount(Draw_K_in_degree)
cat("Number of edges in Draw_K_in_degree:", num_edges_Q, "\n")

plot(Draw_K_in_degree,
     vertex.label = NA,
     vertex.size = 5,
     edge.arrow.size = 0.3,
     main = "K-Type In-Degree Support Network")

#more complex social network drawing
q_tbl_graph <- as_tbl_graph(Draw_K_in_degree)

ggraph(q_tbl_graph, layout = "fr") +  # Or try "kk" for Kamada-Kawai layout
  geom_edge_link() +
  geom_node_point() +
  theme_void() +
  ggtitle("K-Type Support Network")

#spatial social networks drawing
#1.match vertex names (Ego IDs) to coordinates in PT_full_dyads
coords_df <- PT_full_dyads[match(V(Draw_K_in_degree)$name, PT_full_dyads$Ego), c("Longitude", "Latitude")]

#2.convert to matrix for spatial layout
layout_coords <- as.matrix(coords_df)

#3.plot the network using spatial coordinates
plot(Draw_K_in_degree,
     layout = layout_coords,
     vertex.label = NA,
     vertex.size = 5,
     edge.arrow.size = 0.3,
     main = "Spatial K-Type In-Degree Network")

#4.colour code

#ensure both are character for matching
V(Draw_K_in_degree)$name <- as.character(V(Draw_K_in_degree)$name)
PT_full_dyads$Ego <- as.character(PT_full_dyads$Ego)

#match and assign Community to vertices
V(Draw_K_in_degree)$Community <- PT_full_dyads$Community[match(V(Draw_K_in_degree)$name, PT_full_dyads$Ego)]

#assign colors: 1 = blue, 2 = orange (customisable)
vertex_colors <- ifelse(V(Draw_K_in_degree)$Community == 1, "mediumblue", "darkorange1")

pdf("SQ_plot.pdf", width = 7, height = 7)

plot(Draw_Q_in_degree,
     layout = layout_coords,
     vertex.label = NA,
     vertex.size = 6,
     vertex.color = vertex_colors,
     edge.arrow.size = 0.5)

dev.off()

#ensure both are character for matching
V(Draw_K_in_degree)$name <- as.character(V(Draw_K_in_degree)$name)
PT_full_dyads$Ego <- as.character(PT_full_dyads$Ego)

#match and assign Community to vertices
V(Draw_K_in_degree)$Community <- PT_full_dyads$Community[match(V(Draw_K_in_degree)$name, PT_full_dyads$Ego)]

#assign colors
vertex_colors <- ifelse(V(Draw_K_in_degree)$Community == 1, "mediumblue", "darkorange1")

#plot
pdf("SQ_K_plot.pdf", width = 7, height = 7)

plot(Draw_K_in_degree,
     layout = layout_coords,
     rescale = TRUE,  
     margin = 0.2,
     vertex.label = NA,
     vertex.size = 3,
     vertex.color = vertex_colors,
     edge.arrow.size = 0.3)

dev.off()


#5.convert igraph to tidygraph
q_tbl_graph <- as_tbl_graph(Draw_K_in_degree)

q_tbl_graph <- q_tbl_graph %>%
  mutate(Community = as.factor(Community))

#plot with community-based color
K_in_SN <- ggraph(q_tbl_graph) +
  geom_edge_link(color = "gray70", alpha = 0.6) +
  geom_node_point(aes(color = Community), size = 4) +
  scale_color_manual(values = c("1" = "mediumblue", "2" = "darkorange1")) +
  theme_void() +
  ggtitle("Spatial K-Type Support Network by Community") +
  theme(legend.position = "right")

K_in_SN #name of social network 

#modifications. . .

#adding weights 
# call target nodes (alters) of each edge
to_ids <- as.character(ends(Draw_K_in_degree, es = E(Draw_K_in_degree), names = TRUE)[, 2])

#count nominations per alter
nomination_counts <- table(to_ids)

#convert to named numeric vector
nomination_vector <- as.numeric(nomination_counts)
names(nomination_vector) <- names(nomination_counts)

#assign weights to each edge based on target node
E(Draw_K_in_degree)$weight <- nomination_vector[to_ids]

q_tbl_graph <- as_tbl_graph(Draw_K_in_degree)

q_tbl_graph <- q_tbl_graph %>%
  activate(nodes) %>%
  mutate(Community = as.factor(Community)) #change community to factor

K_in_SN <- ggraph(q_tbl_graph) +
  geom_edge_link(aes(width = weight), color = "gray70", alpha = 0.5, show.legend = FALSE) +  # Edge weights
  geom_node_point(aes(color = Community), size = 4, alpha = 0.8) +      # Node transparency
  scale_color_manual(
    values = c("1" = "mediumblue", "2" = "darkorange1"),
    name = "Community",
    labels = c("Community 1", "Community 2")
  ) +
  scale_edge_width(range = c(0.5, 2.5)) +  # Adjust edge width scaling
  theme_void() +
  ggtitle("Spatial K-Type Support Network by Community") +
  theme(legend.position = "right")

K_in_SN

#add node sized by age

#match age to each node by Ego ID
V(Draw_K_in_degree)$Age <- PT_full_dyads$Age[match(V(Draw_K_in_degree)$name, PT_full_dyads$Ego)]

q_tbl_graph <- as_tbl_graph(Draw_K_in_degree)

q_tbl_graph <- q_tbl_graph %>%
  activate(nodes) %>%
  mutate(Community = as.factor(Community)) #change community to factor


K_in_SN <- ggraph(q_tbl_graph) +
  geom_edge_link(aes(width = weight), color = "gray70", alpha = 0.5, show.legend = FALSE) +
  geom_node_point(aes(color = Community, size = Age), alpha = 0.6, show.legend = FALSE) +
  scale_color_manual(
    values = c("1" = "mediumblue", "2" = "darkorange1"),
    name = "Community",
    labels = c("Community 1", "Community 2")
  ) +
  scale_edge_width(range = c(0.5, 2.5)) +
  scale_size(range = c(2, 8), name = "Age") +
  theme_void() +
  ggtitle("Spatial K-Type Support Network") +
  theme(legend.position = "right")

K_in_SN

#add arrows 

K_in_SN <- ggraph(q_tbl_graph)+
  geom_edge_link(
    aes(width = weight),
    color = "gray70",
    alpha = 0.4,
    show.legend = FALSE,
    arrow = arrow(length = unit(2, "mm"), type = "closed"),
    end_cap = circle(4, 'mm')
  ) +
  geom_node_point(aes(color = Community, size = Age), alpha = 0.8) +  # show legend for color
  scale_color_manual(
    values = c("1" = "mediumblue", "2" = "darkorange1"),
  ) +
  scale_edge_width(range = c(0.2, 1.2)) +
  scale_size(range = c(3, 10), guide = "none") +  # hide size (age) legend
  theme_void() +
  theme(legend.position = "none")

K_in_SN

ggsave("K_In-Degree_SN_plot.pdf", plot = K_in_SN, width = 8, height = 6, units = "in")
ggsave("K_SN_plot_landscape.pdf", plot = K_in_SN, width = 11, height = 8, units = "in")

#Fruchterman-Reingold force-directed layout (layout = "fr")

#OUT-DEGREE at the individual-level (not used in thesis)

#Q OUT-DEGREE 

#change community to factor 
class(PT_dyads_Q$Community)
PT_dyads_Q$Community <- as.factor(PT_dyads_Q$Community)

#simple network drawing 
Draw_Q_out_degree <- graph_from_data_frame(PT_dyads_Q, directed = TRUE)

plot(Draw_Q_out_degree,
     vertex.label = NA,
     vertex.size = 5,
     edge.arrow.size = 0.3,
     main = "Q-Type Out-Degree Support Network")


#more complex social network drawing
q_tbl_graph <- as_tbl_graph(Draw_Q_out_degree)

ggraph(q_tbl_graph, layout = "fr") +  # Or try "kk" for Kamada-Kawai layout
  geom_edge_link() +
  geom_node_point() +
  theme_void() +
  ggtitle("K-Type Support Network")


#spatial social networks drawing
#1.match vertex names (Ego IDs) to coordinates in PT_full_dyads
coords_df <- PT_full_dyads[match(V(Draw_Q_out_degree)$name, PT_full_dyads$Ego), c("Longitude", "Latitude")]

#2.convert to matrix for spatial layout
layout_coords <- as.matrix(coords_df)

#3.plot the network using spatial coordinates
plot(Draw_Q_out_degree,
     layout = layout_coords,
     vertex.label = NA,
     vertex.size = 5,
     edge.arrow.size = 0.3,
     main = "Spatial Q-Type Out-Degree Network")


#4.colour code

#ensure both are character for matching
V(Draw_Q_out_degree)$name <- as.character(V(Draw_Q_out_degree)$name)
PT_full_dyads$Ego <- as.character(PT_full_dyads$Ego)

#match and assign Community to vertices
V(Draw_Q_out_degree)$Community <- PT_full_dyads$Community[match(V(Draw_Q_out_degree)$name, PT_full_dyads$Ego)]

#assign colors: 1 = blue, 2 = orange (customisable)
vertex_colors <- ifelse(V(Draw_Q_out_degree)$Community == 1, "mediumblue", "darkorange1")

pdf("SQ_plot2.pdf", width = 7, height = 7)

plot(Draw_Q_out_degree,
     layout = layout_coords,
     vertex.label = NA,
     vertex.size = 3,
     vertex.color = vertex_colors,
     edge.arrow.size = 0.5)

dev.off()

#5.convert igraph to tidygraph
q_tbl_graph <- as_tbl_graph(Draw_Q_out_degree)

q_tbl_graph <- q_tbl_graph %>%
  mutate(Community = as.factor(Community))

#plot with community-based color
Q_out_SN <- ggraph(q_tbl_graph) +
  geom_edge_link(color = "gray70", alpha = 0.6) +
  geom_node_point(aes(color = Community), size = 4) +
  scale_color_manual(values = c("1" = "mediumblue", "2" = "darkorange1")) +
  theme_void() +
  ggtitle("Spatial Q-Type out-degree Support Network by Community") +
  theme(legend.position = "right")

Q_out_SN #name of social network 

#modifications. . .

#adding weights 
#call target nodes (alters) of each edge
to_ids <- as.character(ends(Draw_Q_out_degree, es = E(Draw_Q_out_degree), names = TRUE)[, 2])

#count nominations per alter
nomination_counts <- table(to_ids)

#convert to named numeric vector
nomination_vector <- as.numeric(nomination_counts)
names(nomination_vector) <- names(nomination_counts)

#assign weights to each edge based on target node
E(Draw_Q_out_degree)$weight <- nomination_vector[to_ids]

q_tbl_graph <- as_tbl_graph(Draw_Q_out_degree)

q_tbl_graph <- q_tbl_graph %>%
  activate(nodes) %>%
  mutate(Community = as.factor(Community)) #change community to factor

Q_out_SN <- ggraph(q_tbl_graph) +
  geom_edge_link(aes(width = weight), color = "gray70", alpha = 0.5, show.legend = FALSE) +  # Edge weights
  geom_node_point(aes(color = Community), size = 4, alpha = 0.8) +      # Node transparency
  scale_color_manual(
    values = c("1" = "mediumblue", "2" = "darkorange1"),
    name = "Community",
    labels = c("Community 1", "Community 2")
  ) +
  scale_edge_width(range = c(0.5, 2.5)) +  # Adjust edge width scaling
  theme_void() +
  ggtitle("Spatial Q-Type Support Network by Community") +
  theme(legend.position = "right")


Q_out_SN <- ggraph(q_tbl_graph) +
  geom_edge_link(
    aes(width = weight),
    color = "gray70",
    alpha = 0.5,
    show.legend = FALSE,
    arrow = arrow(length = unit(3, "mm"), type = "closed"),
    end_cap = circle(3, "mm")
  ) +
  geom_node_point(aes(color = Community), size = 4, alpha = 0.8) +
  scale_color_manual(
    values = c("1" = "mediumblue", "2" = "darkorange1"),
    name = "Community",
    labels = c("Community 1", "Community 2")
  ) +
  scale_edge_width(range = c(0.5, 2.5)) +
  theme_void() +
  ggtitle("Spatial Q-Type Support Network by Community") +
  theme(legend.position = "right")


#add NAs as other 

q_tbl_graph <- q_tbl_graph %>%
  activate(nodes) %>%
  mutate(
    Community = as.character(Community),            
    Community = ifelse(is.na(Community), "Other", Community),  # replace NAs
    Community = factor(Community, levels = c("1", "2", "Other"))  #set factor levels explicitly
  )

q_tbl_graph %>% activate(nodes) %>% as_tibble() %>% 
  count(Community)

library(ggraph)
library(grid)

Q_out_SN <- ggraph(q_tbl_graph) +
  geom_edge_link(
    aes(width = weight),
    color = "gray70",
    alpha = 0.4,
    show.legend = FALSE,
    arrow = arrow(length = unit(2, "mm"), type = "closed"),
    end_cap = circle(4, "mm")
  ) +
  geom_node_point(aes(color = Community), size = 4, alpha = 0.8) +  # Removed size = Age
  scale_color_manual(
    values = c("1" = "mediumblue", "2" = "darkorange1", "Other" = "gray70")
  ) +
  scale_edge_width(range = c(0.2, 1.2)) +
  theme_void() +
  theme(legend.position = "none")

Q_out_SN

ggsave("Q_out_SN_landscape.pdf", plot =Q_out_SN, width = 11, height = 8, units = "in")


#K out-degree 

#change community to factor 
class(PT_dyads_K$Community)
PT_dyads_K$Community <- as.factor(PT_dyads_K$Community)

#simple network drawing 
Draw_K_out_degree <- graph_from_data_frame(PT_dyads_K, directed = TRUE)

plot(Draw_K_out_degree,
     vertex.label = NA,
     vertex.size = 5,
     edge.arrow.size = 0.3,
     main = "K-Type Out-Degree Support Network")


#more complex social network drawing
q_tbl_graph <- as_tbl_graph(Draw_K_out_degree)

ggraph(q_tbl_graph, layout = "fr") +  # Or try "kk" for Kamada-Kawai layout
  geom_edge_link() +
  geom_node_point() +
  theme_void() +
  ggtitle("K-Type Support Network")


#spatial social networks drawing
#1.match vertex names (Ego IDs) to coordinates in PT_full_dyads
coords_df <- PT_full_dyads[match(V(Draw_K_out_degree)$name, PT_full_dyads$Ego), c("Longitude", "Latitude")]

#2.convert to matrix for spatial layout
layout_coords <- as.matrix(coords_df)

#3.plot the network using spatial coordinates
plot(Draw_K_out_degree,
     layout = layout_coords,
     vertex.label = NA,
     vertex.size = 5,
     edge.arrow.size = 0.3,
     main = "Spatial K-Type Out-Degree Network")


#4.colour code

#ensure both are character for matching
V(Draw_K_out_degree)$name <- as.character(V(Draw_K_out_degree)$name)
PT_full_dyads$Ego <- as.character(PT_full_dyads$Ego)

#match and assign Community to vertices
V(Draw_K_out_degree)$Community <- PT_full_dyads$Community[match(V(Draw_K_out_degree)$name, PT_full_dyads$Ego)]

#assign colors: 1 = blue, 2 = orange (customisable)
vertex_colors <- ifelse(V(Draw_K_out_degree)$Community == 1, "mediumblue", "darkorange1")

pdf("SQ_plot1.pdf", width = 7, height = 7)

plot(Draw_K_out_degree,
     layout = layout_coords,
     vertex.label = NA,
     vertex.size = 3,
     vertex.color = vertex_colors,
     edge.arrow.size = 0.5)

dev.off()

#5.convert igraph to tidygraph
q_tbl_graph <- as_tbl_graph(Draw_K_out_degree)

q_tbl_graph <- q_tbl_graph %>%
  mutate(Community = as.factor(Community))

#plot with community-based color
K_out_SN <- ggraph(q_tbl_graph) +
  geom_edge_link(color = "gray70", alpha = 0.6) +
  geom_node_point(aes(color = Community), size = 4) +
  scale_color_manual(values = c("1" = "mediumblue", "2" = "darkorange1")) +
  theme_void() +
  ggtitle("Spatial K-Type out-degree Support Network by Community") +
  theme(legend.position = "right")

K_out_SN #name of social network 

#modifications. . .

#adding weights 
#call target nodes (alters) of each edge
to_ids <- as.character(ends(Draw_K_out_degree, es = E(Draw_K_out_degree), names = TRUE)[, 2])

#count nominations per alter
nomination_counts <- table(to_ids)

#convert to named numeric vector
nomination_vector <- as.numeric(nomination_counts)
names(nomination_vector) <- names(nomination_counts)

#assign weights to each edge based on target node
E(Draw_K_out_degree)$weight <- nomination_vector[to_ids]

q_tbl_graph <- as_tbl_graph(Draw_K_out_degree)

q_tbl_graph <- q_tbl_graph %>%
  activate(nodes) %>%
  mutate(Community = as.factor(Community)) #change community to factor

K_out_SN <- ggraph(q_tbl_graph) +
  geom_edge_link(aes(width = weight), color = "gray70", alpha = 0.5, show.legend = FALSE) +  # Edge weights
  geom_node_point(aes(color = Community), size = 4, alpha = 0.8) +      # Node transparency
  scale_color_manual(
    values = c("1" = "mediumblue", "2" = "darkorange1"),
    name = "Community",
    labels = c("Community 1", "Community 2")
  ) +
  scale_edge_width(range = c(0.5, 2.5)) +  # Adjust edge width scaling
  theme_void() +
  ggtitle("Spatial Q-Type Support Network by Community") +
  theme(legend.position = "right")


K_out_SN <- ggraph(q_tbl_graph) +
  geom_edge_link(
    aes(width = weight),
    color = "gray70",
    alpha = 0.5,
    show.legend = FALSE,
    arrow = arrow(length = unit(3, "mm"), type = "closed"),
    end_cap = circle(3, "mm")
  ) +
  geom_node_point(aes(color = Community), size = 4, alpha = 0.8) +
  scale_color_manual(
    values = c("1" = "mediumblue", "2" = "darkorange1"),
    name = "Community",
    labels = c("Community 1", "Community 2")
  ) +
  scale_edge_width(range = c(0.5, 2.5)) +
  theme_void() +
  ggtitle("Spatial K-Type Support Network by Community") +
  theme(legend.position = "right")


#add NAs as other 

q_tbl_graph <- q_tbl_graph %>%
  activate(nodes) %>%
  mutate(
    Community = as.character(Community),              # convert to char to modify
    Community = ifelse(is.na(Community), "Other", Community),  # replace NAs
    Community = factor(Community, levels = c("1", "2", "Other"))  # set factor levels explicitly
  )

q_tbl_graph %>% activate(nodes) %>% as_tibble() %>% 
  count(Community)

library(ggraph)
library(grid)

K_out_SN <- ggraph(q_tbl_graph) +
  geom_edge_link(
    aes(width = weight),
    color = "gray70",
    alpha = 0.4,
    show.legend = FALSE,
    arrow = arrow(length = unit(2, "mm"), type = "closed"),
    end_cap = circle(4, "mm")
  ) +
  geom_node_point(aes(color = Community), size = 4, alpha = 0.8) +  # Removed size = Age
  scale_color_manual(
    values = c("1" = "mediumblue", "2" = "darkorange1", "Other" = "gray70")
  ) +
  scale_edge_width(range = c(0.2, 1.2)) +
  theme_void() +
  theme(legend.position = "none")

K_out_SN

ggsave("K_out_SN_landscape.pdf", plot =K_out_SN, width = 11, height = 8, units = "in")

#NBTs ####

#Q NBTs

length(unique(edges_NBT$.tail))
names(edges_NBT)[names(edges_NBT) == ".tail"] <- "Ego"
names(edges_NBT)[names(edges_NBT) == ".head"] <- "Alter"
community2_edges <- edges_NBT[edges_NBT$Community == 2, ]
community1_edges <- edges_NBT[edges_NBT$Community == 1, ]

unique_egos_2 <- unique(community2_edges$Ego) #16 community 2
unique_egos_1 <- unique(community1_edges$Ego) #27 community 1

#count unique egos in Community 2
ego_count_community2 <- length(unique_egos_community2)
#check 'Community' is a factor in PT_ego_list
names(PT_ego_list)[names(PT_ego_list) == "DisplacementStatus"] <- "Community"

PT_ego_list$Community[PT_ego_list$Community == 1] <- 2
PT_ego_list$Community[PT_ego_list$Community == 0] <- 1

PT_ego_list$Community <- as.factor(PT_ego_list$Community)
class(PT_ego_list$Community)

#create igraph object from edge list
Draw_Q_NBTs <- graph_from_data_frame(edges_NBT, directed = TRUE)

#convert to tbl_graph and attach node attributes
q_tbl_graph <- as_tbl_graph(Draw_Q_NBTs) %>%
  activate(nodes) %>%
  left_join(PT_ego_list, by = c("name" = "Ego"))  #replace 'Ego' if your ID column is named differently

#draw network plot

Q_NBTs_plot <- ggraph(q_tbl_graph, layout = "fr") +
  geom_edge_link(
    width = 0.5,            # fixed edge width instead of mapping to weight
    color = "gray70",
    alpha = 0.4,
    show.legend = FALSE,
    arrow = arrow(length = unit(2, "mm"), type = "closed"),
    end_cap = circle(4, 'mm')
  ) +
  geom_node_point(aes(color = Community, size = Age), alpha = 0.8) +
  scale_color_manual(
    values = c("1" = "mediumblue", "2" = "darkorange1")
  ) +
  scale_size(range = c(3, 10), guide = "none") +
  theme_void() +
  theme(legend.position = "none")


#display the plot
print(Q_NBTs_plot)

q_tbl_graph %>%
  as_tibble() %>%
  count(Community)


#add weights

#extract target nodes (alters) of each edge in Draw_Q_NBTs
to_ids <- as.character(ends(Draw_Q_NBTs, es = E(Draw_Q_NBTs), names = TRUE)[, 2])

#count how many nominations (edges) each alter receives
nomination_counts <- table(to_ids)

#convert to named numeric vector
nomination_vector <- as.numeric(nomination_counts)
names(nomination_vector) <- names(nomination_counts)

#assign weights to edges based on the nomination count of the target node
E(Draw_Q_NBTs)$weight <- nomination_vector[to_ids]

#convert to tidygraph and join node attributes
q_tbl_graph <- as_tbl_graph(Draw_Q_NBTs) %>%
  activate(nodes) %>%
  left_join(PT_ego_list, by = c("name" = "Ego"))  #adjust if necessary

#plot with edge widths scaled by weight
Q_NBTs_plot <- ggraph(q_tbl_graph, layout = "fr") +
  geom_edge_link(
    aes(width = weight),
    color = "gray70",
    alpha = 0.4,
    show.legend = FALSE,
    arrow = arrow(length = unit(2, "mm"), type = "closed"),
    end_cap = circle(4, 'mm')
  ) +
  geom_node_point(aes(color = Community, size = Age), alpha = 0.8) +
  scale_color_manual(values = c("1" = "mediumblue", "2" = "darkorange1")) +
  scale_edge_width(range = c(0.2, 1.2)) +
  scale_size(range = c(3, 10), guide = "none") +
  theme_void() +
  theme(legend.position = "none")

#display the plot
print(Q_NBTs_plot) 

ggsave("Q_NBTs_plot.pdf", plot = Q_NBTs_plot, width = 10, height = 8, units = "in")


################################################################################

#end session

################################################################################