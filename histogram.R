#calculate # of groups overall with a certain boundary total (across all species)
library(googlesheets4)
TAD_groups <- read_sheet("1nYCwVnE4Fm0XRWNQuYxQAS7or4EjZD0o6VKaBO47dUc",sheet="Boundaries_Groups")
TAD_groups <- TAD_groups[,-1]

TAD_groups_with_sum <- transform(TAD_groups, sum = rowSums(TAD_groups))