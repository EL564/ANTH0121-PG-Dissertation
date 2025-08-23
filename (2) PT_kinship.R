################################################################################
#
#
#(2) PT_kinship 
#Consanguineous kinship (i.e., blood-related kin, excluding affinal kin such as in-laws and spouses)
#author: Ella Lipscombe 
#
################################################################################

#start session 

#shortcut (saved csvs from this session)
setwd("D:/MSc HEB 2023-2025/Modules 2024-2025/ANTH0121 Dissertation/Methodology/R script")
PT_dyads_kin_clean <- read.csv ("PT_dyads_kin_clean.csv", header = TRUE) #code to import master data set
PT_dyads_kin_pedigree <- read.csv ("PT_dyads_kin_pedigree.csv", header = TRUE) #code to import master data set

#packages
library(igraph)
library(dplyr)
library(tidyr)
library(sf)
library(ggplot2)
library(kinship2)

#Full pedigree ####
#extract unique individuals from both ego and alter columns
ego_data <- data.frame (
  id = PT_dyads_kin_clean$Ego,
  dadid = PT_dyads_kin_clean$Ego_Father,
  momid = PT_dyads_kin_clean$Ego_Mother,
  su = PT_dyads_kin_clean$Ego_SU,
  rel_to_hh = PT_dyads_kin_clean$Ego_relationshiptoHH
)

alter_data <- data.frame (
  id = PT_dyads_kin_clean$Alter,
  dadid = PT_dyads_kin_clean$Alter_Father,
  momid = PT_dyads_kin_clean$Alter_Mother,
  su = PT_dyads_kin_clean$Alter_SU,
  rel_to_hh = PT_dyads_kin_clean$Alter_relationshiptoHH
)

#combine and get unique individuals
individuals <- rbind(ego_data, alter_data)
individuals <- individuals[!duplicated(individuals$id), ]

#create sex variable based on parent roles
#start with NA for everyone
individuals$sex <- NA

#anyone who appears as a father gets "M" 
individuals$sex[individuals$id %in% individuals$dadid] <- "M"

#anyone who appears as a mother gets "F"
individuals$sex[individuals$id %in% individuals$momid] <- "F"

#for anyone still NA (not a parent in the data), assign a default
individuals$sex[is.na(individuals$sex)] <- "M"

#get all missing parent IDs (both fathers and mothers)
all_parent_ids <- unique(c(individuals$dadid, individuals$momid))
all_parent_ids <- all_parent_ids[!is.na(all_parent_ids) & all_parent_ids != 0]

missing_parents <- setdiff(all_parent_ids, individuals$id)

print("Missing parents that need to be added:")
print(missing_parents)

#create records for these missing parents
missing_parent_records <- data.frame(
  id = missing_parents,
  dadid = 0,  #parent unknown
  momid = 0,  #parent unknown
  su = NA,    #unknown sharing unit
  rel_to_hh = NA,  #unknown relationship to household
  sex = ifelse(missing_parents %in% individuals$dadid, "M", "F")
)

#add these missing parents to the individuals dataset
individuals <- rbind(individuals, missing_parent_records)

#convert 0s to NAs for missing parents
individuals$dadid[individuals$dadid == 0 | is.na(individuals$dadid)] <- NA
individuals$momid[individuals$momid == 0 | is.na(individuals$momid)] <- NA

#create the pedigree
ped <- pedigree(id = individuals$id,
                dadid = individuals$dadid, 
                momid = individuals$momid,
                sex = individuals$sex)


#calculate kinship matrix
kinship_matrix <- kinship(ped)

#basic summary of all kinship coefficients
summary(as.vector(kinship_matrix))

#mean kinship coefficient
mean(kinship_matrix, na.rm = TRUE)

#standard deviation
sd(as.vector(kinship_matrix), na.rm = TRUE)

#range
range(kinship_matrix, na.rm = TRUE)

#get frequency table of kinship values
kinship_values <- as.vector(kinship_matrix)
table(kinship_values)

#or for more detailed breakdown
library(dplyr)
kinship_df <- data.frame(kinship = as.vector(kinship_matrix))
kinship_df %>% 
  count(kinship, sort = TRUE) %>%
  mutate(percentage = n/sum(n)*100)

#extract kinship coefficients for your original dyadic data
PT_dyads_kin_clean$kinship_coeff <- mapply(function(ego, alter) {
  kinship_matrix[as.character(ego), as.character(alter)]
}, PT_dyads_kin_clean$Ego, PT_dyads_kin_clean$Alter)

#summary of kinship for your actual ties
summary(PT_dyads_kin_clean$kinship_coeff)
table(PT_dyads_kin_clean$kinship_coeff)

#categorize your dyadic relationships
PT_dyads_kin_clean$relationship_type <- categorize_kinship(PT_dyads_kin_clean$kinship_coeff)
table(PT_dyads_kin_clean$relationship_type)

#make sure you have kinship coefficients for your dyadic data
PT_dyads_kin_clean$kinship_coeff <- mapply(function(ego, alter) {
  kinship_matrix[as.character(ego), as.character(alter)]
}, PT_dyads_kin_clean$Ego, PT_dyads_kin_clean$Alter)

#categorise
categorize_kinship <- function(k) {
  case_when(
    k == 0.5 ~ "Parent-Child/Full Siblings",
    k == 0.25 ~ "Grandparent-Grandchild/Half Siblings/Uncle-Nephew",
    k == 0.125 ~ "First Cousins/Great Grandparent-Great Grandchild",
    k == 0.0625 ~ "First Cousins Once Removed",
    k == 0 ~ "Unrelated",
    k > 0 & k < 0.0625 ~ "Second Cousins",
    TRUE ~ "Other"
  )
}

#histogram of kinship coefficients
hist(as.vector(kinship_matrix), 
     main = "Distribution of Kinship Coefficients",
     xlab = "Kinship Coefficient",
     breaks = 20)

#for your dyadic data specifically
hist(PT_dyads_kin_clean$kinship_coeff,
     xlab = "Kinship Coefficient")

write.csv(PT_dyads_kin_clean, "PT_dyads_kin_pedigree.csv", row.names = FALSE)

PT_dyads_kin_pedigree$is_kin <- PT_dyads_kin_pedigree$kinship_coeff > 0

################################################################################

#end session

################################################################################
