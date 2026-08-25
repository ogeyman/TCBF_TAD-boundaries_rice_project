#Convert PANTHER results to PRISM graphing

##HC
library(googlesheets4)
panther <- read_sheet("1VRpPU8Q6aKESwGkFafnihkuFiGxa8Te8blWwtSvddvU",sheet='BP - HC')
panther <- panther[c(1:8,10),]
panther_prism <- panther[,c("Term","Fold enrichment","FDR")]
panther_prism$Term <- sub(" \\(GO:.*\\)$", "", panther_prism$Term) #keeps only GO descriptions, not ID
write.csv(panther_prism,"panther_data_for_prism.csv",row.names=FALSE)