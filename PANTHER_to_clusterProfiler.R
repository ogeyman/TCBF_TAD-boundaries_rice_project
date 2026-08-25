#plotting PANTHER GO terms on clusterProfiler
library(googlesheets4)
library(dplyr)
library(stringr)
library(clusterProfiler)
library(enrichplot)

#FOR HC
panther <- read_sheet("1VRpPU8Q6aKESwGkFafnihkuFiGxa8Te8blWwtSvddvU",sheet='BP - HC')

panther <- panther[c(1:8,10),]

names(panther) <- names(panther) %>% 
  str_replace_all("\\s+", "_") %>% 
  str_replace_all("[^A-Za-z0-9_]", "")

# Extract GO ID & description
panther <- panther %>% 
  mutate(
    ID = str_extract(Term, "GO:\\d+"),
    Description = str_trim(str_remove(Term, "\\(GO:\\d+\\)"))
  )

# Identify correct column names
names(panther)[names(panther) == "_in_REF"] <- "ref_col"
names(panther)[names(panther) == "_in_QUERY"] <- "query_col"

# Compute denominators
total_query <- max(panther$query_col, na.rm = TRUE)
total_ref   <- max(panther$ref_col, na.rm = TRUE)

# Add required clusterProfiler columns
panther <- panther %>%
  mutate(
    GeneRatio = paste0(.data$query_col, "/", total_query),
    BgRatio   = paste0(.data$ref_col, "/", total_ref),
    pvalue    = Raw_Pvalue,
    p.adjust  = FDR,
    qvalue    = FDR,
    Count     = .data$query_col,
    geneID    = NA
  )

# Reorder into the format enrichResult expects
cp_df <- panther %>%
  select(ID, Description, GeneRatio, BgRatio,
         pvalue, p.adjust, qvalue, geneID, Count)

# Create an enrichResult object
ego_panther <- new(
  "enrichResult",
  result = cp_df,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  qvalueCutoff = 1,
  gene = panther$ID,    # or any non-empty vector
  ontology = "BP"
)

# Plotting
dotplot(ego_panther,title="Enriched Genes in HC TAD Boundaries")


#FOR LC
panther <- read_sheet("1VRpPU8Q6aKESwGkFafnihkuFiGxa8Te8blWwtSvddvU",sheet='BP - LC')

names(panther) <- names(panther) %>% 
  str_replace_all("\\s+", "_") %>% 
  str_replace_all("[^A-Za-z0-9_]", "")

# Extract GO ID & description
panther <- panther %>% 
  mutate(
    ID = str_extract(Term, "GO:\\d+"),
    Description = str_trim(str_remove(Term, "\\(GO:\\d+\\)"))
  )

# Identify correct column names
names(panther)[names(panther) == "_in_REF"] <- "ref_col"
names(panther)[names(panther) == "_in_QUERY"] <- "query_col"

# Compute denominators
total_query <- max(panther$query_col, na.rm = TRUE)
total_ref   <- max(panther$ref_col, na.rm = TRUE)

# Add required clusterProfiler columns
panther <- panther %>%
  mutate(
    GeneRatio = paste0(.data$query_col, "/", total_query),
    BgRatio   = paste0(.data$ref_col, "/", total_ref),
    pvalue    = Raw_Pvalue,
    p.adjust  = FDR,
    qvalue    = FDR,
    Count     = .data$query_col,
    geneID    = NA
  )

# Reorder into the format enrichResult expects
cp_df <- panther %>%
  select(ID, Description, GeneRatio, BgRatio,
         pvalue, p.adjust, qvalue, geneID, Count)

# Create an enrichResult object
ego_panther <- new(
  "enrichResult",
  result = cp_df,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  qvalueCutoff = 1,
  gene = panther$ID,    # or any non-empty vector
  ontology = "BP"
)

# Plotting
dotplot(ego_panther,title="Enriched Genes in LC TAD Boundaries")
