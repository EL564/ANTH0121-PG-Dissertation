################################################################################
#
#
#(4) PT_kin_density
#session includes: kin density calculation
#author: Ella Lipscombe 
#
################################################################################

#start session 

#set working directory
setwd("D:/MSc HEB 2023-2025/Modules 2024-2025/ANTH0121 Dissertation/Methodology/R script")

#shortcut:
PT_master <- read.csv ("PT_master.csv", header = TRUE) #code to import master data set
PT_full_dyads <- read.csv ("PT_full_dyads.csv", header = TRUE) #code to import master data set
SU_combined <- read.csv ("SU_combined.csv", header = TRUE) #import 

PT_edge_list <- read.csv ("PT_edge_list.csv", header = TRUE) #code to import master data set
PT_ego_list <- read.csv ("PT_ego_list.csv", header = TRUE) #code to import master data set

Q_NBTs <- read.csv ("Q_NBTs.csv", header = TRUE) 
Q_NBTs_kin_density  <- read.csv ("Q_NBTs_kin_density .csv", header = TRUE) #import

################################################################################

#Q: NBT kin density

#load packages
library(dplyr)
library(tidyr)
library(purrr)

#create an alter-alter data frame, grouped by ego
alter_pairs_df <- Q_NBTs %>%
  group_by(Ego) %>%
  summarise(
    alters = list(unique(Alter)),
    .groups = "drop"
  ) %>%
  mutate(
    alter_pairs = map(alters, ~ {
      alters_vec <- .x
      n_alters <- length(alters_vec)
      
      if (n_alters < 2) {
        #include ego with one alter: create one row with Alter1 = alter, Alter2 = NA
        tibble(
          Alter1 = alters_vec,
          Alter2 = NA_character_
        )
      } else {
        #create all combinations of alters (2 at a time)
        pairs <- combn(alters_vec, 2)
        tibble(
          Alter1 = pairs[1, ],
          Alter2 = pairs[2, ]
        )
      }
    })
  ) %>%
  select(Ego, alter_pairs) %>%
  unnest(cols = c(alter_pairs))

#is alter kin to ego? 

#select the relevant ego-alter kinship info from Q_reciprocity
ego_alter_kin_df <- Q_NBTs %>%
  select(Ego, Alter, is_kin)

#join to get Alter1 kin info
alter_pairs_df <- alter_pairs_df %>%
  left_join(
    ego_alter_kin_df,
    by = c("Ego", "Alter1" = "Alter")
  ) %>%
  rename(Alter1_is_kin_to_ego = is_kin)

#join again to get Alter2 kin info
alter_pairs_df <- alter_pairs_df %>%
  left_join(
    ego_alter_kin_df,
    by = c("Ego", "Alter2" = "Alter")
  ) %>%
  rename(Alter2_is_kin_to_ego = is_kin)

#view result
alter_pairs_df

#count the alter-alter kin ties 
alter_pairs_df <- alter_pairs_df %>%
  mutate(
    alter_alter_kin = ifelse(
      Alter1_is_kin_to_ego & Alter2_is_kin_to_ego,
      "YES",
      "NO"
    )
  )

alter_pairs_df <- alter_pairs_df %>%
  group_by(Ego) %>%
  mutate(
    alter_alter_kin_count = sum(alter_alter_kin == "YES", na.rm = TRUE)
  ) %>%
  ungroup()

#aggregate by ego

alter_alter_kin_summary <- alter_pairs_df %>%
  group_by(Ego) %>%
  summarise(
    total_alter_alter_kin = sum(alter_alter_kin == "YES", na.rm = TRUE),
    total_alter_alter_ties = n(),  # number of alter–alter pairs per ego
    .groups = "drop"
  ) %>%
  #add ego–alter ties count from edges_NBT
  left_join(
    Q_NBTs %>%
      group_by(Ego) %>%
      summarise(total_ego_alter_ties = n(), .groups = "drop"),
    by = "Ego"
  ) %>%
  #calculate total ties = alter–alter ties + ego–alter ties
  mutate(
    total_ties = total_alter_alter_ties + total_ego_alter_ties
  )


#subset the Q NBT list to include ego, alter, is_kin, and dyadic distance only 

Q_NBTs_subset <- Q_NBTs %>%
  group_by(Ego) %>%
  summarise(
    ego_alter_kin_count = sum(is_kin == TRUE, na.rm = TRUE),
    .groups = "drop"
  )

#join data under a new data frame called kin density 

kin_density_NBT <- Q_NBTs_subset %>%
  left_join(alter_alter_kin_summary, by = "Ego")


#calculate kin density as a new variable

#variables required:
#ego_kin_count: columns Ego, ego_alter_kin_count
#alter_alter_kin_summary: columns Ego, total_alter_alter_kin, total_alter_alter_ties, total_ego_alter_ties, total_ties

Q_NBTs_kin_density <- kin_density_NBT %>%
  mutate(
    kin_density_NBT = ifelse(
      total_ties > 0,
      (ego_alter_kin_count + total_alter_alter_kin) / total_ties,
      NA_real_
    )
  )

#save 
write.csv(Q_NBTs_kin_density , "Q_NBTs_kin_density .csv", row.names = FALSE) #save
Q_NBTs_kin_density  <- read.csv ("Q_NBTs_kin_density .csv", header = TRUE) #import


################################################################################

#edge-list preparation for ERGM
#to include the following variables: ego, alter, question no., support type, in-degree (Y/N), out-degree (Y/N), reciprocity (Y/N), NBTs (Y/N), is kin?, kin relatedness, kin density, dyadic distance, displacement status, community, distance from PA 

#vector of column 
PT_master <- PT_dyads_kin_distance

#create new, empty columns in master dataset 
# Add new empty variables (set to NA)
PT_master$Support_type <- NA
PT_master$In_degree <- NA
PT_master$Out_degree <- NA
PT_master$Reciprocity <- NA
PT_master$NBTs <- NA
PT_master$Kin_density_QNBTs <- NA

#move support type next to questions
after_col <- "Question"

cols <- names(PT_master) #get current column names

#remove support_type from current position
PT_master$support_type -> temp_support
PT_master$support_type <- NULL

pos <- match(after_col, cols) #find position of 'questions'

new_order <- c(cols[1:pos], "Support_type", cols[(pos+1):length(cols)]) #rebuild column order: everything up to pos, then support_type, then rest

PT_master <- PT_master[, new_order] #reorder dataframe columns

PT_master$Support_type.1 <- NULL


#1.fill in/merge variables 

PT_master$Support_type <- ifelse(PT_master$Question >= 1 & PT_master$Question <= 5, "Q", "K")

#2.fill in in-degree 
PT_master$In_degree <- apply(PT_master, 1, function(row) {
  ego <- row["Ego"]
  alter <- row["Alter"]
  #look for a row where Ego and Alter are reversed
  any(PT_master$Ego == alter & PT_master$Alter == ego)
})

#convert logical TRUE/FALSE to "YES"/"NO"
PT_master$In_degree <- ifelse(PT_master$In_degree, "YES", "NO")
sum(PT_master$In_degree == "YES", na.rm = TRUE)
sum(PT_master$In_degree == "NO", na.rm = TRUE)

#fill in out-degrees
PT_master$Out_degree <- "YES"

#fill in reciprocity 
PT_master$Reciprocity <- ifelse(
  PT_master$Ego != PT_master$Alter & 
    paste(PT_master$Alter, PT_master$Ego) %in% paste(PT_master$Ego, PT_master$Alter),
  "YES", "NO"
)

#fill in needs-based transfers 
PT_master$NBTs <- ifelse(PT_master$Question %in% c(1, 3) & PT_master$Support_type == "Q", "YES", "NO")

#remove the last three variables in PT_master 
PT_master <- PT_master[, -((ncol(PT_master)-2):ncol(PT_master))]

#change name of PT_master to edge_list
PT_edge_list <- PT_master

write.csv(PT_edge_list, "PT_edge_list.csv", row.names = FALSE) #save cleaned master dataset
PT_edge_list <- read.csv ("PT_edge_list.csv", header = TRUE) #code to import master data set

################################################################################

#ego list for ERGM

#change name of SU_combined
PT_ego_list <- SU_combined

#merge all network-specific kin density 
#ensure all variable names are consistent before merging

names(Q_NBTs_kin_density)[names(Q_NBTs_kin_density) == "kin_density"] <- "Kin_density_QNBTs"

PT_ego_list <- merge(PT_ego_list, Q_rec_kin_density[, c("Ego", "Kin_density_QREC")], by = "Ego", all.x = TRUE)
PT_ego_list <- merge(PT_ego_list, K_rec_kin_density[, c("Ego", "Kin_density_KREC")], by = "Ego", all.x = TRUE)
PT_ego_list <- merge(PT_ego_list, Q_NBTs_kin_density[, c("Ego", "Kin_density_QNBTs")], by = "Ego", all.x = TRUE)

PT_ego_list <- PT_ego_list %>%
  left_join(Q_NBTs_kin_density %>% select(Ego, kin_density_NBT), by = "Ego") #new NBT kin density
PT_ego_list$Kin_density_QNBTs <- NULL #remove old kin density NBT list 
names(PT_ego_list)[names(PT_ego_list) == "kin_density_NBT"] <- "Kin_density_QNBT" #rename new kin density

write.csv(PT_ego_list, "PT_ego_list.csv", row.names = FALSE) #save cleaned master dataset
PT_ego_list <- read.csv ("PT_ego_list.csv", header = TRUE) #code to import master data set

###############################################################################

#end session

################################################################################