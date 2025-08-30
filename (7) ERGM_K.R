################################################################################
#
#
#(7) PT_ERGM_behavioural(K)_network
#includes: pre-diagnostic tests, ERGM analysis, and post-diagnostic tests (GOOF) 
#author: Ella Lipscombe 
#note: in august choose whether to keep or remove relatedness
#
################################################################################

#start session

#working directory
setwd("D:/MSc HEB 2023-2025/Modules 2024-2025/ANTH0121 Dissertation/Methodology/R script")

#shortcuts 
PT_ego_list <- read.csv ("PT_ego_list.csv", header = TRUE) #code to import master data set
PT_edge_list <- read.csv ("PT_edge_list.csv", header = TRUE) #code to import master data set
pairwise_distances_full <- read.csv ("pairwise_distances_full.csv", header = TRUE) #code to import master data set
pedigree_dyads_full <- read.csv ("pedigree_dyads_full.csv", header = TRUE)

load("net_K_support.RData")
load("net_K_support_clean.RData") #use this one 

#load packages
library(network)
library(sna)
library(ergm)
library(statnet)  
library(dplyr)
library(car)

#[0] correlation and linear regression ####
edges <- as.data.frame(as.edgelist(net_K_support))

edges$kin_relatedness <- get.edge.attribute(net_K_support, "kinship_coeff")
edges$distance <- get.edge.attribute(net_K_support, "Dyadic_distance")

shapiro.test(edges$kin_relatedness) #not normal
shapiro.test(edges$distance) #not normal 

cor.test(edges$kin_relatedness, edges$distance, method = "spearman")

#[1] pre-model diagnostic checks [inspect the network variables] ####
summary(net_K_support)
list.vertex.attributes(net_K_support)
list.edge.attributes(net_K_support)

#trim incomplete wealth data
complete_wealth <- !is.na(net_K_support %v% "MaterialWealth")

#create new network with complete cases
net_K_clean <- net_K_support
net_K_clean <- delete.vertices(net_K_clean, which(!complete_wealth))

#trim dyadic distances 
summary(net_K_clean %e% "Dyadic_distance")
sum(is.na(net_K_clean %e% "Dyadic_distance"))
complete_dist <- !is.na(net_K_clean %e% "Dyadic_distance")

#create new network without NA distance edges
net_K_support <- delete.edges(net_K_clean, which(!complete_dist))

#look at actual values
head(net_K_support %e% "Dyadic_distance", 100) #89 observations for distance (because there are 88 dyads)

#*pre-model diagnostic checks* 
#extract variables of interest relevant to H1
material_wealth <- get.vertex.attribute(net_K_support, "MaterialWealth")
age <- get.vertex.attribute(net_K_support, "Age")
gender <- get.vertex.attribute(net_K_support, "Gender")
dyadic_distance <- get.edge.attribute(net_K_support, "Dyadic_distance")

#(1) outliers [inspect data] 
detect_outliers <- function(x, var_name) {
  K1 <- quantile(x, 0.25, na.rm = TRUE)
  K3 <- quantile(x, 0.75, na.rm = TRUE)
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
png("KSUPPORT_boxplots_outlier_detection.png", width = 10, height = 8, units = "in", res = 300)
par(mfrow = c(2, 2))
boxplot(kin_density, main = "Kin Density", ylab = "Density")
boxplot(material_wealth, main = "Material Wealth", ylab = "Wealth")
boxplot(age, main = "Age", ylab = "Years")
boxplot(dyadic_distance, main = "Dyadic Distance", ylab = "Distance")
dev.off()
cat("Boxplots exported to: KSUPPORTboxplots_outlier_detection.png\n")

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

#(2) normal distribution 
material_wealth <- net_K_support %v% "MaterialWealth"
age <- net_K_support %v% "Age"
dyadicdistance <- net_K_support %e% "Dyadic_distance"

get.edge.attribute(net_K_support, "Dyadic_distance")

#remove distance outliers
net_K_support <- delete.edges(net_K_support , c(69, 99, 122, 142))

get.edge.attribute(net_K_support , "Dyadic_distance") #outliers have been removed! so 52 dyads remaining

par(mfrow = c(2, 2))

hist(kin_density, main = "Kin Density", xlab = "Kin Density", col = "grey", breaks = 20)
KKnorm(kin_density); KKline(kin_density)

hist(material_wealth, main = "Material Wealth", xlab = "Material Wealth", col = "grey", breaks = 20)
qqnorm(material_wealth); KKline(material_wealth)

hist(age, main = "Age", xlab = "Age", col = "grey", breaks = 20)
qqnorm(age); KKline(age)

hist(dyadicdistance, main = "Dyadic Distance", xlab = "Dyadic Distance", col = "grey", breaks = 20)
qqnorm(dyadicdistance); KKline(dyadicdistance)

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
class(net_K_support)
network.size(net_K_support) #54 egos in network (now outliers and missing data from key variables are removed)

#list all vertex attributes
list.vertex.attributes(net_K_support)

#extract node attributes if interest 
#based on your ERGM, extract these node attributes:
node_data <- data.frame(
  Kin_density_Ksupport = net_K_support_clean %v% "Kin_density_Ksupport",
  MaterialWealth = net_K_support_clean %v% "MaterialWealth",
  Age = net_K_support_clean %v% "Age",
  Gender = net_K_support_clean %v% "Gender"
)

#check the structure and lengths
str(node_data)
sapply(node_data, length)
sapply(node_data, class)

kin_density <- net_K_support_clean %v% "Kin_density_Ksupport"
material_wealth <- net_K_support_clean %v% "MaterialWealth"
age <- net_K_support_clean %v% "Age"
gender <- net_K_support_clean %v% "Gender"

class(kin_density)
class(material_wealth)
class(age)
class(gender) 

node_data <- data.frame(
  Kin_density_Ksupport = as.numeric(kin_density),
  MaterialWealth = as.numeric(material_wealth),
  Age = as.numeric(age),
  Gender = as.factor(gender))

str(node_data)
dim(node_data) #45 egos; 4 variables

library(car)
lm_model <- lm(I(rep(1, length(kin_density))) ~ kin_density + material_wealth + age)
vif(lm_model)

#save clean support network
save(net_K_support, file = "net_K_support.RData")
load("net_K_support.RData")

#[2] matrix formation####

#inspect attributes 
class(net_K_support)
as.edgelist(net_K_support) #dyadic attributes 
edges_support <- as.data.frame.network(net_K_support, matrix.type = "edgelist") #88 dyadic observations
summary(net_K_support) #summary of 127 dyadic observations

list.vertex.attributes(net_K_support) #list all dyadic attributes 
list.edge.attributes(net_K_support) #list all edge attributes 

net_K_support
plot(net_K_clean)
summary(net_K_support)
length(unique(edges_support$.tail)) #54 egos 
length(unique(edges_support$.head)) #91 alters 
network.size(net_K_support) #105 nodes in total
network.edgecount(net_K_support) #250 edges in total
network.density(net_K_support) #network density = 0.0531
summary(net_K_support ~ idegree(0:10) + odegree(0:10)) #indegree and outdegree descriptives
deg_summary <- summary(net_K_support ~ idegree(0:10) + odegree(0:10))
deg_summary

library(sna)
recip_value <- reciprocity(net_K_support)
print(recip_value)

K_degree_plot <- barplot(deg_summary,
        ylab = "Number of Nodes",
        las = 2,
        col = c("darkgreen"))

#covariate descriptives

#--> material wealth 
summary(get.vertex.attribute(net_K_clean, "MaterialWealth"))
IQR(get.vertex.attribute(net_K_clean, "MaterialWealth"), na.rm = TRUE)

#--> age 
summary(get.vertex.attribute(net_K_clean, "Age"))
sd(get.vertex.attribute(net_K_clean, "Age"), na.rm = TRUE)

#--> dyadic distance
dyadic_distance <- get.edge.attribute(net_K_support, "Dyadic_distance")
summary(dyadic_distance)
IQR(get.edge.attribute(net_K_support, "Dyadic_distance"))

#--> kin coef
kin <- get.edge.attribute(net_K_support, "kinship_coeff")
summary(kin)
sd(get.edge.attribute(net_K_support, "kinship_coeff"))

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
network_verts <- network.vertex.names(net_K_support)

pairwise_distances_K_support <- pairwise_distances_full %>%
  filter(Ego %in% network_verts & Alter %in% network_verts)

#checks
unique_egos <- unique(c(pairwise_distances_K_support$Ego, pairwise_distances_K_support$Alter))
setdiff(network_verts, unique_egos)  # Should return character(0)

#matrix formation 
library(tidyr)
library(dplyr)
library(tibble)

#pivot to wide format: rows = Ego, columns = Alter, values = distance
distance_K_support_matrix <- pairwise_distances_K_support%>%
  select(Ego, Alter, distance) %>%
  pivot_wider(names_from = Alter, values_from = distance) %>%
  column_to_rownames("Ego") %>%
  as.matrix()

#match vertex order
vertex_order <- network.vertex.names(net_K_support)

---------------
  #kinship matrix 
  PT_dyads_kin_pedigree <- read.csv ("PT_dyads_kin_pedigree.csv", header = TRUE) #code to import master data set

#now. .  subset the list based on my network
network_verts <- network.vertex.names(net_K_REC_clean)

pairwise_kinship_KREC <- PT_dyads_kin_pedigree %>%
  filter(Ego %in% network_verts & Alter %in% network_verts)

#checks
unique_egos <- uniKue(c(pairwise_kinship_KREC$Ego, pairwise_kinship_KREC$Alter))
setdiff(network_verts, uniKue_egos)  # Should return character(0)

#remove duplicates 
pairwise_kinship_KREC$dyad_id <- apply(pairwise_kinship_KREC[c("Ego", "Alter")], 1, function(x) paste(sort(x), collapse = "_"))

#keep only the first instance of each uniKue dyad
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
missing_dads <- uniKue(kinship_meta_su$Father[!kinship_meta_su$Father %in% ego_ids & kinship_meta_su$Father != 0])
missing_moms <- uniKue(kinship_meta_su$Mother[!kinship_meta_su$Mother %in% ego_ids & kinship_meta_su$Mother != 0])

missing_parents <- uniKue(c(missing_dads, missing_moms))
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
vertex_ids <- network.vertex.names(net_K_support)

#subset your pedigree dyads to only those where both Ego and Alter are in the network
pedigree_dyads_K_support <- pedigree_dyads_full %>%
  filter(Ego %in% vertex_ids & Alter %in% vertex_ids)

#create a 'dyad_id' for undirected duplicate removal (sorted concatenation)
pedigree_dyads_K_support <- pedigree_dyads_K_support %>%
  mutate(dyad_id = ifelse(Ego < Alter, paste(Ego, Alter, sep = "-"), paste(Alter, Ego, sep = "-")))

#remove duplicates based on dyad_id
pedigree_dyads_K_support <- pedigree_dyads_K_support %>%
  distinct(dyad_id, .keep_all = TRUE) %>%
  select(-dyad_id)

#pivot to wide matrix format
kinship_K_support_matrix <- pedigree_dyads_K_support %>%
  pivot_wider(names_from = Alter, values_from = Kinship_Coefficient) %>%
  column_to_rownames("Ego")

#reorder matrix rows and columns to match network vertex order
vertex_order <- vertex_ids
kinship_K_support_matrix <- kinship_K_support_matrix[vertex_order, vertex_order]

#replace any NAs (optional; replace with 0 or mean)
kinship_K_support_matrix[is.na(kinship_K_support_matrix)] <- 0
kinship_K_support_matrix <- as.matrix(kinship_K_support_matrix)

#finding missing vertex

n <- network.size(net_K_support)
vertex_names <- network.vertex.names(net_K_support)
print(n)              # Should be 54
print(vertex_names)   

setdiff(vertex_names, rownames(kinship_K_support_matrix))
setdiff(rownames(kinship_K_support_matrix), vertex_names)

vertex_names <- network.vertex.names(net_K_support)

#reorder rows and columns to network vertex order
kinship_K_support_matrix <- kinship_K_support_matrix[vertex_names, vertex_names]

#rescale 
kinship_K_support_matrix_rescaled <- kinship_K_support_matrix / 0.1


#correlation between node attributes
vertex_data <- data.frame(
  Age = get.vertex.attribute(net_K_support, "Age"),
  Gender = as.numeric(as.factor(get.vertex.attribute(net_K_support, "Gender"))),
  MaterialWealth = get.vertex.attribute(net_K_support, "MaterialWealth")
)

cor(vertex_data, use = "pairwise.complete.obs")

install.packages("Rglpk")
library(Rglpk)


#[3] ERGM ####
#MODELS -- propensity to form a support (outcome) ~ covariates 
#in order of increasing complexity . . .
model_K_support_1 <- ergm(net_K_support ~ edges, 
                     control = control.ergm(MCMLE.maxit = 100, MCMC.samplesize = 200,MCMC.burnin = 200, MCMC.interval = 100,seed = 123))

model_K_support_2 <- ergm(net_K_support ~ edges 
                     + mutual, control = control.ergm(MCMLE.maxit = 100, MCMC.samplesize = 200,MCMC.burnin = 200, MCMC.interval = 100,seed = 123))
                    
model_K_support_3 <- ergm(net_K_support ~ edges 
                     + mutual 
                     + nodecov("Age") + nodematch("Gender") + nodefactor ("Gender"), 
                     control = control.ergm(MCMLE.maxit = 100, MCMC.samplesize = 200,MCMC.burnin = 200, MCMC.interval = 100,seed = 123))

model_K_support_4 <- ergm(net_K_support ~ edges 
                     + mutual 
                     + nodecov("Age") + nodematch("Gender") + nodefactor ("Gender") 
                     + nodecov("MaterialWealth") + nodeocov("MaterialWealth") + absdiff("MaterialWealth"),
                     control = control.ergm(MCMLE.maxit = 100, MCMC.samplesize = 200,MCMC.burnin = 200, MCMC.interval = 100,seed = 123))

model_K_support_5 <- ergm (net_K_support ~ edges + mutual + nodecov("Age") 
                           + nodematch("Gender") + nodefactor("Gender") + edgecov(kinship_K_support_matrix_rescaled)
                           + edgecov(distance_K_support_matrix), control = control.ergm(MCMLE.maxit = 100, MCMC.samplesize = 200, MCMC.burnin = 200, MCMC.interval = 100, seed = 123))

save(model_K_support_1, model_K_support_2, model_K_support_3, model_K_support_4, model_K_support_5, file = "models_K_support_all.RData")
load("models_K_support_all.RData")

summary(model_K_support_1) #baseline (edges only)
summary(model_K_support_2) #+ control variables (gender and age)
summary(model_K_support_3) #+ kin 
summary(model_K_support_4) #+ reciprocity 
summary(model_K_support_5) #+ ego wealth + alter wealth + wealth homophily + structural variables (+ control variables)


AIC(model_K_support_1, model_K_support_2, model_K_support_3, model_K_support_4, model_K_support_5)
logLik(model_K_support_1)
logLik(model_K_support_2)
logLik(model_K_support_3)
logLik(model_K_support_4)
logLik(model_K_support_5)

length(coef(model_K_support_1)) 
length(coef(model_K_support_2)) 
length(coef(model_K_support_3)) 
length(coef(model_K_support_4)) 
length(coef(model_K_support_5)) 
length(coef(model_K_support_6)) 
948-785

#model_K_support_5 (best fitting model) report logistic regression results 
#regression coefficients 
a <-  -5.039e+00 #edges
b1 <-   2.755e+00 #mutual
b2 <- 1.838e+00#kinship matrix
b3 <- 1.784e-02  #age 
b4 <-  7.625e-01 #gender homophily
b5 <- -6.448e-02 #gender (m)
b6 <-    -1.183e-04#distance matrix 


#baseline odds of supports 
exp(a) 
#OR associated with covariates vs baseline
exp(b1)  15.72104
exp(b2)  6.283958
exp(b3)  1.018
exp(b4)  2.143629
exp(b5)  0.9375549
exp(b6)  0.9998817


exp(confint(model_K_support_5)) #CI intervals
exp(a)/(1+exp(a)) #predicted proportion for the baseline group = 0.0158993

#table of results 
summary_model <- summary(model_K_support_5)

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
gof_step1 <- gof(model_K_support_5)
print(gof_step1)
plot(gof_step1)
par(mfrow = c(2, 3))
plot(gof_step1)

library(latticeExtra)

mcmc.diagnostics(model_Ksupport_6)
par(mfrow = c(3, 3))  # or c(4, 3) depending on number of parameters




################################################################################

#end session 

################################################################################