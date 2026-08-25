#Identify unique TAD boundaries (+ conserved TAD boundaries) for Meridionalis
library(readr)
library(stringr)
library(dplyr)

merid_total_boundaries <- read.table("meridionalis.boundary.bed",header=T,sep=",")
merid_total_boundaries1 <- merid_total_boundaries[,1]
all_cons_boundaries <- read.table("TAD_groups_copy.csv",header=T,sep=",")
merid_cons_boundaries <- all_cons_boundaries %>%
  pull(2) %>%
  na.omit() %>%
  paste(collapse = ";") %>%         # Combine all rows into one string
  str_split(";") %>%                # Split on semicolons
  unlist() %>%
  str_trim()                        # Remove leading/trailing whitespace
merid_unique_boundaries <- setdiff(merid_total_boundaries1,merid_cons_boundaries) #identify number of unique boundaries

full_info <- merid_total_boundaries[merid_total_boundaries$tad_name %in% unique_boundaries,]
write.csv(full_info, "unique_merid_boundaries.csv", row.names = FALSE)

#Code below is to format it for "comparative analysis boundaries 5kb" script
unique_boundaries <- full_info[,2:4] #removes boundary names
unique_boundaries[,1] <- gsub("meridionalis_","",unique_boundaries[,1]) #removed meridionalis prefix from chr column
colnames(unique_boundaries) <- NULL #removes header
write.table(unique_boundaries, file = "merid_unique_boundaries_ANON.bed", 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE) #ANON means annonymous as we removed the names



