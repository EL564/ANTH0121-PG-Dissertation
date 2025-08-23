################################################################################
#
#
#(8) PT_ERGM_NBT (Q: NBT)_network
#includes: pre-diagnostic tests, ERGM analysis, and post-diagnostic tests (GOOF) 
#author: Ella Lipscombe 
#note: in august choose whether to keep or remove relatedness
#
################################################################################

#start session 

#working directory
setwd("D:/MSc HEB 2023-2025/Modules 2024-2025/ANTH0121 Dissertation/Methodology/R script")

PT_ego_list <- read.csv ("PT_ego_list.csv", header = TRUE) #code to import master data set
PT_edge_list <- read.csv ("PT_edge_list.csv", header = TRUE) #code to import master data set
pairwise_distances_full <- read.csv ("pairwise_distances_full.csv", header = TRUE) #code to import master data set
pedigree_dyads_full <- read.csv ("pedigree_dyads_full.csv", header = TRUE)

load("net_Q_NBT.RData") #Q_NBT combined network (i.e., combined edge and node attrbutes)
load("net_Q_NBT_clean.RData") #USE THIS cleaned up version post-checks (outliers and missing data removed)

#load packages
library(network)
library(sna)
library(ergm)
library(statnet)  
library(dplyr)
library(car)

#[1] pre-model diagnostic checks [inspect the network variables] ####
summary(net_Q_NBT)
list.vertex.attributes(net_Q_NBT)
list.edge.attributes(net_Q_NBT)
get.vertex.attribute(net_Q_NBT, "Kin_density_QNBT")


#eliminate missing kin_density_NBT values from network 
missing_nodes <- which(is.na(get.vertex.attribute(net_Q_NBT, "Kin_density_QNBT")))
network.vertex.names(net_Q_NBT)[missing_nodes]

net_Q_NBT_trimmed <- net_Q_NBT
delete.vertices(net_Q_NBT_trimmed, missing_nodes)
summary(net_Q_NBT_trimmed) #note for discussion: 193 dyads in original NBT network but only 90 dyads after the trim (because kin density was not attainable across certain individuals, affecting 93 dyads)
get.vertex.attribute(net_Q_NBT_trimmed, "Kin_density_QNBT")


#trim incomplete wealth data
complete_wealth <- !is.na(net_Q_NBT_trimmed %v% "MaterialWealth")

#create new network with complete cases
net_Q_NBT_complete <- net_Q_NBT_trimmed
net_Q_NBT_complete <- delete.vertices(net_Q_NBT_complete, which(!complete_wealth))

#trim dyadic distances 
summary(net_Q_NBT_complete %e% "Dyadic_distance")
sum(is.na(net_Q_NBT_complete %e% "Dyadic_distance"))
complete_dist <- !is.na(net_Q_NBT_complete %e% "Dyadic_distance")

#create new network without NA distance edges
net_Q_NBT_final <- delete.edges(net_Q_NBT_complete, which(!complete_dist))

rm(net_Q_NBT_trimmed, net_Q_NBT_complete) #remove older subsets to clear up environment 

#look at actual values
head(net_Q_NBT_final %e% "Dyadic_distance", 100) #89 observations for distance (because there are 88 dyads)

#*pre-model diagnostic checks* 
#extract variables of interest relevant to H1
kin_density <- get.vertex.attribute(net_Q_NBT_final, "Kin_density_QNBT")
material_wealth <- get.vertex.attribute(net_Q_NBT_final, "MaterialWealth")
age <- get.vertex.attribute(net_Q_NBT_final, "Age")
gender <- get.vertex.attribute(net_Q_NBT_final, "Gender")
dyadic_distance <- get.edge.attribute(net_Q_NBT_final, "Dyadic_distance")

#(1) outliers [inspect data] 
detect_outliers <- function(x, var_name) {
  Q1 <- quantile(x, 0.25, na.rm = TRUE)
  Q3 <- quantile(x, 0.75, na.rm = TRUE)
  IQR <- Q3 - Q1
  lower_bound <- Q1 - 1.5 * IQR
  upper_bound <- Q3 + 1.5 * IQR
  
  outliers <- which(x < lower_bound | x > upper_bound)
  
  cat("\n", var_name, ":\n")
  cat("  Lower bound:", lower_bound, "\n")
  cat("  Upper bound:", upper_bound, "\n")
  cat("  Number of outliers:", length(outliers), "\n")
  cat("  Outlier indices:", outliers, "\n")
  cat("  Outlier values:", x[outliers], "\n")
  
  return(outliers)
}

#check each variable for outliers
outliers_kin <- detect_outliers(kin_density, "Kin Density")
outliers_wealth <- detect_outliers(material_wealth, "Material Wealth")
outliers_age <- detect_outliers(age, "Age")
outlier_distance <- detect_outliers(dyadic_distance, "Dyadic_Distance")

#boxplots for visual inspection
par(mfrow = c(2, 2))
boxplot(kin_density, main = "Kin Density", ylab = "Density")
boxplot(material_wealth, main = "Material Wealth", ylab = "Wealth")
boxplot(age, main = "Age", ylab = "Years")
boxplot(dyadic_distance, main = "Dyadic Distance", ylab = "Dyadic Distance")

#export boxplots for appendix 
png("boxplots_outlier_detection.png", width = 10, height = 8, units = "in", res = 300)
par(mfrow = c(2, 2))
boxplot(kin_density, main = "Kin Density", ylab = "Density")
boxplot(material_wealth, main = "Material Wealth", ylab = "Wealth")
boxplot(age, main = "Age", ylab = "Years")
boxplot(dyadic_distance, main = "Dyadic Distance", ylab = "Distance")
dev.off()
cat("Boxplots exported to: boxplots_outlier_detection.png\n")

#summary statistics
cat("\n")
cat("SUMMARY STATISTICS:\n")
cat("-------------------\n")
summary_stats <- rbind(
  "Kin Density" = summary(kin_density),
  "Material Wealth" = summary(material_wealth),
  "Age" = summary(age),
  "Dyadic Distance" = summary(dyadic_distance)
)

print(summary_stats)

#export summary statistics for appendix
write.csv(summary_stats, "summary_statistics_variables.csv", row.names = TRUE)
cat("Summary statistics exported to: summary_statistics_variables.csv\n")

#notes for results: after careful inspection of outliers, I am choosing to remove the outliers in kin density (as there are 3 that reflect self nomination) and also the extreme values in dyadic distance as they are skewing the data
#create a copy first
net_Q_NBT_clean <- net_Q_NBT_final

#(2) normal distribution 
kin_density <- net_Q_NBT_clean %v% "Kin_density_QNBT"
material_wealth <- net_Q_NBT_clean %v% "MaterialWealth"
age <- net_Q_NBT_clean %v% "Age"
dyadicdistance <- net_Q_NBT_clean %e% "Dyadic_distance"

get.edge.attribute(net_Q_NBT_clean, "Dyadic_distance")

#remove distance outliers
net_Q_NBT_clean <- delete.edges(net_Q_NBT_clean, c(52, 68, 78))

get.edge.attribute(net_Q_NBT_clean, "Dyadic_distance") #outliers have been removed! so 52 dyads remaining

par(mfrow = c(2, 2))

hist(kin_density, main = "Kin Density", xlab = "Kin Density", col = "grey", breaks = 20)
qqnorm(kin_density); qqline(kin_density)

hist(material_wealth, main = "Material Wealth", xlab = "Material Wealth", col = "grey", breaks = 20)
qqnorm(material_wealth); qqline(material_wealth)

hist(age, main = "Age", xlab = "Age", col = "grey", breaks = 20)
qqnorm(age); qqline(age)

hist(dyadicdistance, main = "Dyadic Distance", xlab = "Dyadic Distance", col = "grey", breaks = 20)
qqnorm(dyadicdistance); qqline(dyadicdistance)

shapiro.test(kin_density)      #kin_density
shapiro.test(material_wealth)  #material_wealth
shapiro.test(age)              #age
shapiro.test(dyadicdistance)   #dyadic distance

summary(kin_density)
summary(material_wealth)
summary(age)
summary(dyadicdistance)

#(3) multicollinarity [to avoid overfitting]

#check your network structure
class(net_Q_NBT_clean)
network.size(net_Q_NBT_clean) #49 egos in network (now outliers and missing data from key variables are removed)

#list all vertex attributes
list.vertex.attributes(net_Q_NBT_clean)

#extract node attributes if interest 
#based on your ERGM, extract these node attributes:
node_data <- data.frame(
  Kin_density_QNBT = net_Q_NBT_clean %v% "Kin_density_QNBT",
  MaterialWealth = net_Q_NBT_clean %v% "MaterialWealth",
  Age = net_Q_NBT_clean %v% "Age",
  Gender = net_Q_NBT_clean %v% "Gender"
)

#check the structure and lengths
str(node_data)
sapply(node_data, length)
sapply(node_data, class)

kin_density <- net_Q_NBT_clean %v% "Kin_density_QNBT"
material_wealth <- net_Q_NBT_clean %v% "MaterialWealth"
age <- net_Q_NBT_clean %v% "Age"
gender <- net_Q_NBT_clean %v% "Gender"
 
class(kin_density)
class(material_wealth)
class(age)
class(gender) 

node_data <- data.frame(
  Kin_density_QNBT = as.numeric(kin_density),
  MaterialWealth = as.numeric(material_wealth),
  Age = as.numeric(age),
  Gender = as.factor(gender))

str(node_data)
dim(node_data) #45 egos; 4 variables

library(car)
lm_model <- lm(I(rep(1, length(kin_density))) ~ kin_density + material_wealth + age)
vif(lm_model)

#save clean NBT network
save(net_Q_NBT_clean, file = "net_Q_NBT_clean.RData")
load("net_Q_NBT_clean.RData")

#[2] matrix formation####

#inspect attributes 
class(net_Q_NBT_clean)
as.edgelist(net_Q_NBT_clean) #dyadic attributes 
edges_NBT <- as.data.frame.network(net_Q_NBT_clean, matrix.type = "edgelist") #88 dyadic observations
summary(net_Q_NBT_clean) #summary of 87 dyadic observations

list.vertex.attributes(net_Q_NBT_clean) #list all dyadic attributes 
list.edge.attributes(net_Q_NBT_clean) #list all edge attributes 

net_Q_NBT_clean
summary(net_Q_NBT_clean)
length(unique(edges_NBT$.tail)) #40 egos 
length(unique(edges_NBT$.head)) #33 alters 
network.size(net_Q_NBT_clean) #49 nodes in total
network.edgecount(net_Q_NBT_clean) #87 edges in total
network.density(net_Q_NBT_clean) #network density = 3.69898%
summary(net_Q_NBT_clean ~ idegree(0:10) + odegree(0:10)) #indegree and outdegree descriptives

summary(net_Q_NBT_clean ~ idegree(0:10) + odegree(0:10)) #indegree and outdegree descriptives
deg_summary <- summary(net_Q_NBT_clean ~ idegree(0:10) + odegree(0:10))
deg_summary

library(sna)
recip_value <- reciprocity(net_Q_NBT_clean)
print(recip_value)

barplot(deg_summary,
        ylab = "Number of Nodes",
        las = 2,
        col = c("yellow", "pink"))


#covariate descriptives
#--> kin density 
summary(get.vertex.attribute(net_Q_NBT_clean, "Kin_density_QNBT"))
median(get.vertex.attribute(net_Q_NBT_clean, "Kin_density_QNBT"))
IQR(get.vertex.attribute(net_Q_NBT_clean, "Kin_density_QNBT"))

#--> material wealth 
summary(get.vertex.attribute(net_Q_NBT_clean, "MaterialWealth"))
IQR(get.vertex.attribute(net_Q_NBT_clean, "MaterialWealth"))

#--> age 
summary(get.vertex.attribute(net_Q_NBT_clean, "Age"))
sd(get.vertex.attribute(net_Q_NBT_clean, "Age"))

#--> dyadic distance
dyadic_distance <- get.edge.attribute(net_Q_NBT_clean, "Dyadic_distance")
summary(dyadic_distance)
IQR(get.edge.attribute(net_Q_NBT_clean, "Dyadic_distance"))

#scale covariates (converting the values into z-scores) -- maybe delete
#note: kin density and relatedness is already a value between 0-1
net_Q_NBT_clean %v% "Age" <- as.vector(scale(net_Q_NBT_clean %v% "Age"))
net_Q_NBT_clean %v% "Material Wealth" <- as.vector(scale(net_Q_NBT_clean %v% "MaterialWealth"))
net_Q_NBT_clean %e% "Dyadic_distance" <- as.vector(scale(net_Q_NBT_clean %e% "Dyadic_distance"))

---------------
#distance matrix
pairwise_distances <- read.csv ("su_distances.csv", header = TRUE) #code to import master data set

#add prefix "PTSU" to make sure the sharing unit IDs match PT_ego_list
library(dplyr)
library(stringr)

pairwise_distances <- pairwise_distances %>%
  mutate(
    sui = paste0("PTSU", str_remove_all(sui, regex("(?i)ptsu[_]*"))),
    suj = paste0("PTSU", str_remove_all(suj, regex("(?i)ptsu[_]*")))
  )

#rename columns
pairwise_distances <- pairwise_distances %>%
  rename(
    `SUID(ego)` = sui,
    `SUID(alter)` = suj
  )

#subset ego list 
ego_su <- PT_ego_list %>%
  select(Ego, SharingUnitID)

#rename egos
ego_su <- ego_su %>%
  rename(`SUID(ego)` = SharingUnitID)

#subset alter list 
alter_su <- PT_ego_list %>%
  select(Ego, SharingUnitID)

#rename altersu and alter
alter_su <- alter_su %>%
  rename(`SUID(alter)` = SharingUnitID)

alter_su <- alter_su %>%
  rename(`Alter` = Ego)

#combine ego and alter to su in pairwise dyads
#join Ego names
pairwise_distances <- pairwise_distances %>%
  left_join(ego_su, by = "SUID(ego)")

# Join Alter names (need to rename for second join)
pairwise_distances <- pairwise_distances %>%
  left_join(
    alter_su, by = "SUID(alter)"
  )

#save as pairwise_distances_full
write.csv(pairwise_distances, "pairwise_distances_full.csv", row.names = FALSE)

#--> continue from here 
#subset the list based on my network
network_verts <- network.vertex.names(net_Q_NBT_clean)

pairwise_distances_QNBT <- pairwise_distances_full %>%
  filter(Ego %in% network_verts & Alter %in% network_verts)

#checks
unique_egos <- unique(c(pairwise_distances_QNBT$Ego, pairwise_distances_QNBT$Alter))
setdiff(network_verts, unique_egos)  # Should return character(0)

#matrix formation 
library(tidyr)
library(dplyr)
library(tibble)

#pivot to wide format: rows = Ego, columns = Alter, values = distance
distance_QNBT_matrix <- pairwise_distances_QNBT%>%
  select(Ego, Alter, distance) %>%
  pivot_wider(names_from = Alter, values_from = distance) %>%
  column_to_rownames("Ego") %>%
  as.matrix()

#match vertex order
vertex_order <- network.vertex.names(net_Q_NBT_clean)
distance_QNBT_matrix <- distance_QNBT_matrix[vertex_order, vertex_order]

---------------
  #kinship matrix 
  PT_dyads_kin_pedigree <- read.csv ("PT_dyads_kin_pedigree.csv", header = TRUE) #code to import master data set

#now. .  subset the list based on my network
network_verts <- network.vertex.names(net_K_REC_clean)

pairwise_kinship_KREC <- PT_dyads_kin_pedigree %>%
  filter(Ego %in% network_verts & Alter %in% network_verts)

#checks
unique_egos <- unique(c(pairwise_kinship_KREC$Ego, pairwise_kinship_KREC$Alter))
setdiff(network_verts, unique_egos)  # Should return character(0)

#remove duplicates 
pairwise_kinship_KREC$dyad_id <- apply(pairwise_kinship_KREC[c("Ego", "Alter")], 1, function(x) paste(sort(x), collapse = "_"))

#keep only the first instance of each unique dyad
pairwise_kinship_KREC <- pairwise_kinship_KREC[!duplicated(pairwise_kinship_KREC$dyad_id), ]

#drop the helper column
pairwise_kinship_KREC$dyad_id <- NULL

#matrix formation 
library(tidyr)
library(dplyr)
library(tibble)

#pivot to wide format: rows = Ego, columns = Alter, values = distance
kinship_KREC_matrix <- pairwise_kinship_KREC %>%
  select(Ego, Alter, kinship_coeff) %>%
  pivot_wider(names_from = Alter, values_from = kinship_coeff) %>%
  column_to_rownames("Ego") %>%
  as.matrix()

#match vertex order
vertex_order <- network.vertex.names(net_K_REC_clean)
distance_KREC_matrix <- distance_KREC_matrix[vertex_order, vertex_order]

#kinship matrix 
#import kinship metadata
kinship_meta <- read.csv ("PT_Kinship.csv", header = TRUE) 

#subset full kinship metadata to head of households only 
kinship_meta_su <- kinship_meta[kinship_meta$Ego %in% PT_ego_list$Ego, ]

#load kinship2
library(kinship2)

# Get list of all individual IDs
ego_ids <- kinship_meta_su$Ego

# Identify dad and mom IDs not in ego IDs
missing_dads <- unique(kinship_meta_su$Father[!kinship_meta_su$Father %in% ego_ids & kinship_meta_su$Father != 0])
missing_moms <- unique(kinship_meta_su$Mother[!kinship_meta_su$Mother %in% ego_ids & kinship_meta_su$Mother != 0])

missing_parents <- unique(c(missing_dads, missing_moms))
missing_parents

# Create placeholder rows
placeholder_parents <- data.frame(
  Ego = missing_parents,
  Father = 0,
  Mother = 0,
  Sex = F
)

# Combine with your original kinship_meta_su
kinship_meta_su <- rbind(kinship_meta_su, placeholder_parents)

kinship_meta_su$Sex <- 1  # assume all male to start
kinship_meta_su$Sex[kinship_meta_su$Ego %in% kinship_meta_su$Mother] <- 2  # overwrite with female if listed as mother

# Convert "0" strings to NA
kinship_meta_su$Father[kinship_meta_su$Father == "0"] <- NA
kinship_meta_su$Mother[kinship_meta_su$Mother == "0"] <- NA

ped_su <- with(kinship_meta_su, pedigree(id = Ego, dadid = Father, momid = Mother, sex = Sex))

#generate the kinship coefficient matrix
kin_matrix <- kinship(ped_su)

#inspect the matrix
print(dim(kin_matrix))
print(head(kin_matrix))

#save kin_matrix
#convert kinship matrix to long (tidy) format
kin_df <- as.data.frame(as.table(kin_matrix))
colnames(kin_df) <- c("Ego", "Alter", "Kinship_Coefficient")

#save to CSV
write.csv(kin_df, "pedigree_dyads_full.csv", row.names = FALSE)

#--> continue from here 
#load data (see shortcut)
library(dplyr)
library(tidyr)
library(network)

#get vertex IDs from your network
vertex_ids <- network.vertex.names(net_Q_NBT_clean)

#subset your pedigree dyads to only those where both Ego and Alter are in the network
pedigree_dyads_QNBT <- pedigree_dyads_full %>%
  filter(Ego %in% vertex_ids & Alter %in% vertex_ids)

#create a 'dyad_id' for undirected duplicate removal (sorted concatenation)
pedigree_dyads_QNBT <- pedigree_dyads_QNBT %>%
  mutate(dyad_id = ifelse(Ego < Alter, paste(Ego, Alter, sep = "-"), paste(Alter, Ego, sep = "-")))

#remove duplicates based on dyad_id
pedigree_dyads_QNBT <- pedigree_dyads_QNBT %>%
  distinct(dyad_id, .keep_all = TRUE) %>%
  select(-dyad_id)

#pivot to wide matrix format
kinship_QNBT_matrix <- pedigree_dyads_QNBT %>%
  pivot_wider(names_from = Alter, values_from = Kinship_Coefficient) %>%
  column_to_rownames("Ego")

#reorder matrix rows and columns to match network vertex order
vertex_order <- vertex_ids
kinship_QNBT_matrix <- kinship_QNBT_matrix[vertex_order, vertex_order]

#replace any NAs (optional; replace with 0 or mean)
kinship_QNBT_matrix[is.na(kinship_QNBT_matrix)] <- 0
kinship_QNBT_matrix <- as.matrix(kinship_QNBT_matrix)

#finding missing vertex

n <- network.size(net_Q_NBT_clean)
vertex_names <- network.vertex.names(net_Q_NBT_clean)
print(n)              # Should be 49
print(vertex_names)   

setdiff(vertex_names, rownames(kinship_QNBT_matrix))
setdiff(rownames(kinship_QNBT_matrix), vertex_names)

vertex_names <- network.vertex.names(net_Q_NBT_clean)

#reorder rows and columns to network vertex order
kinship_QNBT_matrix <- kinship_QNBT_matrix[vertex_names, vertex_names]

#rescale 
kinship_QNBT_matrix_rescaled <- kinship_QNBT_matrix / 0.1

#[3] ERGM ####
#MODELS -- propensity to form a NBT (outcome) ~ covariates 
#in order of increasing complexity . . .

model_QNBT_1 <- ergm(net_Q_NBT_clean ~ edges, 
                     control = control.ergm(MCMLE.maxit = 100, MCMC.samplesize = 200,MCMC.burnin = 200, MCMC.interval = 100,seed = 123))

model_QNBT_2 <- ergm(net_Q_NBT_clean ~ edges 
                      + mutual + gwesp(decay = 0.5, fixed = TRUE) + gwdsp(decay = 0.5),
                       control = control.ergm(MCMLE.maxit = 100, MCMC.samplesize = 200,MCMC.burnin = 200, MCMC.interval = 100,seed = 123))

model_QNBT_3 <- ergm(net_Q_NBT_clean ~ edges 
                     + mutual + gwesp(decay = 0.5, fixed = TRUE) + gwdsp(decay = 0.5)
                     + nodecov("Age") + nodematch("Gender") + nodefactor ("Gender"), 
                      control = control.ergm(MCMLE.maxit = 100, MCMC.samplesize = 200,MCMC.burnin = 200, MCMC.interval = 100,seed = 123))

model_QNBT_4 <- ergm(net_Q_NBT_clean ~ edges 
                     + mutual + gwesp(decay = 0.5, fixed = TRUE) + gwdsp(decay = 0.5)
                     + nodecov("Age") + nodematch("Gender") + nodefactor ("Gender") 
                     + nodecov("MaterialWealth") + nodeocov("MaterialWealth") + absdiff("MaterialWealth"),
                     control = control.ergm(MCMLE.maxit = 100, MCMC.samplesize = 200,MCMC.burnin = 200, MCMC.interval = 100,seed = 123))

model_QNBT_5 <- ergm(net_Q_NBT_clean ~ edges 
                     + mutual + gwdsp(decay = 0.5, fixed = TRUE) + gwesp(decay = 0.5, fixed = TRUE) 
                     + nodecov("Age") + nodematch("Gender") + nodefactor ("Gender") 
                     + nodecov("MaterialWealth") + nodeocov("MaterialWealth") + absdiff("MaterialWealth")
                     + nodecov("Kin_density_QNBT") +  edgecov (kinship_QNBT_matrix_rescaled) +  edgecov(distance_QNBT_matrix),
                     control = control.ergm(MCMLE.maxit = 100, MCMC.samplesize = 200,MCMC.burnin = 200, MCMC.interval = 100,seed = 123))


save(model_QNBT_1, model_QNBT_2, model_QNBT_3, model_QNBT_4, model_QNBT_5, file = "models_QNBT_all.RData")
load("models_QNBT_all.RData")

summary(model_QNBT_1) # baseline (edges only)
summary(model_QNBT_2) # structural (mutual, GWESP, GDESP)
summary(model_QNBT_3) # control (age and gender)
summary(model_QNBT_4) # covariates (wealth)
summary(model_QNBT_5) # key predictors

AIC(model_QNBT_1, model_QNBT_2, model_QNBT_3, model_QNBT_4, model_QNBT_5)
logLik(model_QNBT_1)
logLik(model_QNBT_2)
logLik(model_QNBT_3)
logLik(model_QNBT_4)
logLik(model_QNBT_5)

length(coef(model_QNBT_1)) 
length(coef(model_QNBT_2)) 
length(coef(model_QNBT_3)) 
length(coef(model_QNBT_4)) 
length(coef(model_QNBT_5)) 


#model_QNBT_5 (best fitting model) report logistic regression results 
#regression coefficients 
a <- -4.353 #edges
b1 <- 1.424 #mutual
b2 <- 6.605e-01  #kinship matrix
b3 <- 5.658e-05 #ego material wealth
b4 <- -5.011e-05 #alter material wealth
b5 <- -1.181e-05 #wealth homophily
b6 <-  9.099e-03 #age 
b7 <-  3.016e-01 #gender homophily
b8 <- 2.249e-01 #gender (m)
b9 <-  -4.738e-05 #distance matrix 
b10 <-  1.094e+00 #edge-wise
b11 <-  -2.458e-01 #dyad-wise


#baseline odds of NBTs 
exp(a) #0.01286815

#OR associated with covariates vs baseline
exp(b1)  #4.153702
exp(b2)  #27.16692
exp(b3)  #1.000057
exp(b4)  #0.9999499
exp(b5)  #0.9999882
exp(b6)  #1.009141
exp(b7)  #1.35202
exp(b8)  #1.252197
exp(b9)  #0.9999526
exp(b10) #2.986195
exp(b11) #0.7820786
 exp(-0.3511)

exp(confint(model_QNBT_5)) #CI intervals
exp(a)/(1+exp(a)) #predicted proportion for the baseline group = 0.0158993

#table of results 
summary_model <- summary(model_QNBT_6)

coefs <- summary_model$coefficients
estimates <- coefs[, "Estimate"]
se <- coefs[, "Std. Error"]
p_values <- coefs[, "Pr(>|z|)"]

#OR and CIs
OR <- exp(estimates)
CI_lower <- exp(estimates - 1.96 * se)
CI_upper <- exp(estimates + 1.96 * se)

variables <- c("Baseline", "Reciprocity", "Kin Relatedness Matrix", "Ego Wealth", "Alter Wealth","Wealth Homophily", "Ego Age", "Gender Homophily", "Gender M", "Distance Matrix", "Edge-Wise", "Dyad-Wise")

table_results <- data.frame(
  Variable = variables,
  Coefficient = round(estimates, 3),
  OR = round(OR, 3),
  `95% CI` = paste0("(", round(CI_lower, 2), ", ", round(CI_upper, 2), ")"),
  `p-value` = round(p_values, 3),
  check.names = FALSE  # So "95% CI" doesn't become "X95..CI"
)

print(table_results, row.names = FALSE)

#[4] post-model checks ####
# goodness of fit model 

#GOOF
gof_step1 <- gof(model_QNBT_6)
print(gof_step1)
plot(gof_step1)
par(mfrow = c(2, 3))
plot(gof_step1)

library(latticeExtra)

mcmc.diagnostics(model_QNBT_6)
par(mfrow = c(3, 3))  # or c(4, 3) depending on number of parameters




################################################################################

#end session 

################################################################################