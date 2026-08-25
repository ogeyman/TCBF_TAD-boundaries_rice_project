#Identify boundary groups with no Meridionalis presence
library(readr)
library(stringr)
library(dplyr)

all_cons_boundaries <- read.table("TAD_groups_copy.csv",header=T,sep=",")
no_merid_groups <- all_cons_boundaries %>%
  filter(meridionalis == "" &
    apply(.[,3:6],1,function(row) all(row != "" & !is.na(row)))) %>% select(3:6)

#A list of all species boundaries in the extracted groups
azu_boundaries_NM <- no_merid_groups %>% #NM means no meridionalis
  pull(1) %>%
  na.omit() %>%
  .[. != ""] %>%
  paste(collapse = ";") %>%         # Combine all rows into one string
  str_split(";") %>%                # Split on semicolons
  unlist() %>%
  str_trim()                        # Remove leading/trailing whitespace
rufi_boundaries_NM <- no_merid_groups %>% #NM means no meridionalis
  pull(2) %>%
  na.omit() %>%
  .[. != ""] %>%
  paste(collapse = ";") %>%         # Combine all rows into one string
  str_split(";") %>%                # Split on semicolons
  unlist() %>%
  str_trim()                        # Remove leading/trailing whitespace
nipp_boundaries_NM <- no_merid_groups %>% #NM means no meridionalis
  pull(3) %>%
  na.omit() %>%
  .[. != ""] %>%
  paste(collapse = ";") %>%         # Combine all rows into one string
  str_split(";") %>%                # Split on semicolons
  unlist() %>%
  str_trim()                        # Remove leading/trailing whitespace
IR64_boundaries_NM <- no_merid_groups %>% #NM means no meridionalis
  pull(4) %>%
  na.omit() %>%
  .[. != ""] %>%
  paste(collapse = ";") %>%         # Combine all rows into one string
  str_split(";") %>%                # Split on semicolons
  unlist() %>%
  str_trim()                        # Remove leading/trailing whitespace

#The full information of all species' total boundaries
nipp_all_boundaries <- read.table("nipponbare.boundary.bed",header=F,sep=",")
azu_all_boundaries <- read.table("azucena.boundary.bed",header=F,sep=",")
IR64_all_boundaries <- read.table("IR64.boundary.bed",header=F,sep=",")
rufi_all_boundaries <- read.table("rufipogon.boundary.bed",header=F,sep=",")

#Matching selected boundaries with corresponding data - CONJOINT 
full_info_N <- nipp_all_boundaries[nipp_all_boundaries[,1] %in% nipp_boundaries_NM,]
full_info_A <- azu_all_boundaries[azu_all_boundaries[,1] %in% azu_boundaries_NM,]
full_info_I <- IR64_all_boundaries[IR64_all_boundaries[,1] %in% IR64_boundaries_NM,]
full_info_R <- rufi_all_boundaries[rufi_all_boundaries[,1] %in% rufi_boundaries_NM,]

#Write full file to computer
write.csv(full_info_N, "nipponbare_boundaries_NM.csv", row.names = FALSE)
write.csv(full_info_A, "azucena_boundaries_NM.csv", row.names = FALSE)
write.csv(full_info_I, "IR64_boundaries_NM.csv", row.names = FALSE)
write.csv(full_info_R, "rufipogon_boundaries_NM.csv", row.names = FALSE)

#Changing format for Jupiter Notebook analysis scripts
N_boundaries <- full_info_N[,2:4]
N_boundaries[,1] <- gsub("nipponbare_","",N_boundaries[,1]) #removed nipponbare prefix from chr column
N_boundaries <- N_boundaries[-1,]
colnames(N_boundaries) <- NULL

A_boundaries <- full_info_A[,2:4]
A_boundaries[,1] <- gsub("azucena_","",A_boundaries[,1]) #removed azucena prefix from chr column
A_boundaries <- A_boundaries[-1,]
colnames(A_boundaries) <- NULL

I_boundaries <- full_info_I[,2:4]
I_boundaries[,1] <- gsub("IR64_","",I_boundaries[,1]) #removed IR64 prefix from chr column
I_boundaries <- I_boundaries[-1,]
colnames(I_boundaries) <- NULL

R_boundaries <- full_info_R[,2:4]
R_boundaries[,1] <- gsub("rufipogon_","",R_boundaries[,1]) #removed rufipogon prefix from chr column
R_boundaries <- R_boundaries[-1,]
colnames(R_boundaries) <- NULL

#Write finalized, annonymized files to computer
write.table(N_boundaries, file = "nipponbare_boundaries_NM_ANON.bed", 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE) #ANON means annonymous as we removed the names
write.table(A_boundaries, file = "azucena_boundaries_NM_ANON.bed", 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE) #ANON means annonymous as we removed the names
write.table(I_boundaries, file = "IR64_boundaries_NM_ANON.bed", 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE) #ANON means annonymous as we removed the names
write.table(R_boundaries, file = "rufipogon_boundaries_NM_ANON.bed", 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE) #ANON means annonymous as we removed the names