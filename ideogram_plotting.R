#Ideogram of conserved nipponbare TAD boundaries

HC_boundaries_colored <- read.table("HC_boundaries.bed")
C_boundaries_colored <- read.table("C_boundaries.bed")
MC_boundaries_colored <- read.table("MC_boundaries.bed")
LC_boundaries_colored <- read.table("LC_boundaries.bed")
unique_boundaries_colored <- read.csv("unique_boundaries.csv",header=TRUE)
unique_boundaries_col <- unique_boundaries_colored[,2:4]

HC_boundaries_colored$color <- "acen"
C_boundaries_colored$color <- "gpos50"
MC_boundaries_colored$color <- "gpos75"
LC_boundaries_colored$color <- "stalk"
unique_boundaries_col$color <- "stalk"

colnames(HC_boundaries_colored) <- c("chr","start","end","gieStain")
colnames(C_boundaries_colored) <- c("chr","start","end","gieStain")
colnames(MC_boundaries_colored) <- c("chr","start","end","gieStain")
colnames(LC_boundaries_colored) <- c("chr","start","end","gieStain")
colnames(unique_boundaries_col) <- c("chr","start","end","gieStain")

HC_boundaries_col <- HC_boundaries_colored[-1,]
C_boundaries_col <- C_boundaries_colored[-1,]
MC_boundaries_col <- MC_boundaries_colored[-1,]
LC_boundaries_col <- LC_boundaries_colored[-1,]

#HC_boundaries_col$chr <- gsub("chr", "", HC_boundaries_col$chr)
#C_boundaries_col$chr <- gsub("chr", "", C_boundaries_col$chr)
#MC_boundaries_col$chr <- gsub("chr", "", MC_boundaries_col$chr)
#LC_boundaries_col$chr <- gsub("chr", "", LC_boundaries_col$chr)

HC_boundaries_col$name <- paste0("HC", seq_len(nrow(HC_boundaries_col)))
C_boundaries_col$name <- paste0("C", seq_len(nrow(C_boundaries_col)))
MC_boundaries_col$name <- paste0("MC", seq_len(nrow(MC_boundaries_col)))
LC_boundaries_col$name <- paste0("LC", seq_len(nrow(LC_boundaries_col)))
unique_boundaries_col$name <- paste0("UNIQ", seq_len(nrow(unique_boundaries_col)))

HC_boundaries_col <- HC_boundaries_col[,c("chr","start","end","name","gieStain")]
C_boundaries_col <- C_boundaries_col[,c("chr","start","end","name","gieStain")]
MC_boundaries_col <- MC_boundaries_col[,c("chr","start","end","name","gieStain")]
LC_boundaries_col <- LC_boundaries_col[,c("chr","start","end","name","gieStain")]
unique_boundaries_col <- unique_boundaries_col[,c("chr","start","end","name","gieStain")]

all_boundaries <- rbind(HC_boundaries_col,C_boundaries_col,MC_boundaries_col,LC_boundaries_col,unique_boundaries_col)

nipp_chr_coords <- read.table("nipp_chr_coords.txt",header=TRUE)
#nipp_chr_coords$chr <- gsub("chr", "", nipp_chr_coords$chr)

centromeres <- data.frame(
  start = c("16,765,400.0","13,617,389.0","19,954,695.0","9,955,623 ","12,501,000.0",
             "15,612,461.0","12,210,944.0","12,961,751.0","6,771,892.0","8,320,734.",
             "12,889,381.0","12,038,689.0"),
  end = c("17,553,481.0","14,214,518.0","20,491,681.0","10,112,504","12,621,066.0",
           "16,338,095.0","12,436,258.0","13,064,074.0","7,241,524.0","8,690,470.0",
           "14,055,729.0","12,499,575.0"),
  stringsAsFactors=FALSE)
centromeres$start <- as.numeric(gsub("[,\\s]", "", centromeres$start))
centromeres$end <- as.numeric(gsub("[,\\s]", "", centromeres$end))
centromeres$name <- paste0("CTM", seq_len(nrow(centromeres)))
centromeres$gieStain <- "acen"
centromeres$chr <- nipp_chr_coords$chr
centromeres <- centromeres[,c("chr","start","end","name","gieStain")]

library(karyoploteR)
nipp.genome <- toGRanges(nipp_chr_coords)
nipp.cytobands <- toGRanges(all_cons_boundaries)
HC_nipp.cytobands <- toGRanges(HC_boundaries_col)
centromeres.cytobands <- toGRanges(centromeres)

kp1 <- plotKaryotype(genome = nipp.genome, cytobands = nipp.cytobands)
kp2 <- plotKaryotype(genome=nipp.genome,cytobands=HC_nipp.cytobands)

kpA <- plotKaryotype(genome = nipp.genome, cytobands = centromeres.cytobands)
kpB <- plotKaryotype(genome = nipp.genome)
