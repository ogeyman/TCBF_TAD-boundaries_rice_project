######PHYLOGENY

###2D Phylogenetic Tree
library(ape)
library(reshape2)

#mash_distances.tab file is obtained from performing the mash pipeline in Terminal
tab <- read.table("mash_distances.tab",header=FALSE)
df <- tab[,c(1,2,3)]
df$V1 <- gsub(".fasta", "", df$V1)
df$V2 <- gsub(".fasta", "", df$V2)
df$V1 <- gsub("_updated", "", df$V1)
df$V2 <- gsub("_updated", "", df$V2)

df$V1 <- gsub("azucena", "O. sativa japonica azucena", df$V1)
df$V2 <- gsub("azucena", "O. sativa japonica azucena", df$V2)
df$V1 <- gsub("nipp", "O. sativa japonica nipponbare", df$V1)
df$V2 <- gsub("nipp", "O. sativa japonica nipponbare", df$V2)
df$V1 <- gsub("mer", "O. meridionalis", df$V1)
df$V2 <- gsub("mer", "O. meridionalis", df$V2)
df$V1 <- gsub("IR64", "O. sativa indica", df$V1)
df$V2 <- gsub("IR64", "O. sativa indica", df$V2)
df$V1 <- gsub("rufi", "O. rufipogon", df$V1)
df$V2 <- gsub("rufi", "O. rufipogon", df$V2)
df$V1 <- gsub("perr", "L. perrieri", df$V1)
df$V2 <- gsub("perr", "L. perrieri", df$V2)

#rearranging
mat <- acast(df, V1 ~ V2, value.var="V3")
mat[lower.tri(mat)] <- t(mat)[lower.tri(mat)]

dist_2d <- as.dist(mat)
tree_2d <- nj(dist_2d)

tree_2d_rooted <- root(tree_2d,outgroup="L. perrieri",resolve.root=TRUE)
plot(tree_2d_rooted,cex=1.2)

###3D-Genomic Phylogenetic Tree
library(vegan)
library(ape)
library(phytools)
library(googlesheets4)
library(ggtree)
library(ggplot2)
library(dplyr)

TADb_data <- read_sheet("1nYCwVnE4Fm0XRWNQuYxQAS7or4EjZD0o6VKaBO47dUc",sheet="Boundaries_Groups")
binary_data <- ifelse(TADb_data >= 1, 1, 0)
binary_data <- binary_data[,-1]

jac_dist <- vegdist(t(binary_data), method="jaccard",binary=TRUE)
jac_matrix <- as.matrix(jac_dist)
colnames(jac_matrix) <- c("O. meridionalis","O. sativa japonica azucena","O. rufipogon","O. sativa japonica nipponbare","O. sativa indica")
row.names(jac_matrix) <- c("O. meridionalis","O. sativa japonica azucena","O. rufipogon","O. sativa japonica nipponbare","O. sativa indica")

#unrooted
tree_3D <- nj(jac_matrix)
plot(tree_3D)

#rooted using midpoint rooting method
tree_3D_mid <- midpoint.root(tree_3D)
plot(tree_3D_mid)

#bootstrapping 
# boot_vals <- boot.phylo(
#   tree_nj,
#   binary_data,
#   FUN = function(x) nj(vegdist(t(x),method="jaccard",binary=TRUE)),
#   B=1000
# )
# boot_percent_clean <- round(boot_percent, 1)
# 
# plot(tree_nj)
# nodelabels(boot_percent_clean,cex=0.8)

#plotting using ggtree (rooted version)
gt <- ggtree(tree_3D_mid)
boot_df <- data.frame(
  node = c(6, 7, 8)
)
gt$data <- dplyr::left_join(gt$data, boot_df, by = "node")

final_plot <- gt +
  geom_tiplab(size = 5, fontface = "italic")
scaled_plot <- final_plot + xlim(0,0.5)
scaled_plot

###MANTEL TEST (2d vs 3d distance matrices)
library(vegan)
phylo <- mantel(dist_2d,jac_dist,method="pearson",permutations=119,na.rm=FALSE)
phylo

###MANTEL DENSITY PLOT
# null distribution of permuted R values
null_vals <- phylo$perm

# observed R
obs_R <- phylo$statistic

# make density plot
df <- data.frame(null = null_vals)

ggplot(df, aes(x = null)) +
  geom_density(fill = "#C7E9F1", alpha = 1, color = "#2C3E50") +
  geom_vline(xintercept = obs_R, color = "red", linewidth = 1.2) +
  annotate("text",
           x = obs_R, y = max(density(null_vals)$y),
           label = paste("Observed R =", round(obs_R, 4)),
           hjust = -0.1, vjust = -1, color = "red", fontface="bold") +
  labs(
    title = "Mantel Test: Null Distribution of Pearson Correlations",
    x = "Mantel",
    y = "density"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black")
  )

ggplot(df, aes(x = null)) +
  geom_density(linewidth = 1, color = "black") +
  geom_vline(xintercept = obs_R, color = "red", linewidth = 1.2) +
  labs(
    x = "Mantel",
    y = "density",
    title="Mantel Test: Null Distribution of Pearson Correlations"
  ) +
  theme_bw(base_size = 16) +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    plot.title = element_text(size = 17)
  )

###MANTEL SCATTER PLOT
library(ggplot2)

dist2D <- as.vector(dist_2d)
dist3D <- as.vector(jac_dist)

# Extract pairwise distances
vec2D <- dist2D[lower.tri(dist2D)]
vec3D <- dist3D[lower.tri(dist3D)]

df_scatter <- data.frame(dist2D = vec2D, dist3D = vec3D)

ggplot(df_scatter, aes(x = dist2D, y = dist3D)) +
  geom_point(size=4, color="#1f78b4") +
  geom_smooth(method="lm", se=FALSE, color="red", linewidth=1.2) +
  labs(
    x = "2D Evolutionary Distance",
    y = "3D Structural Distance",
    title = "Correlation Between 2D and 3D Genomic Distances"
  ) +
  theme_bw(base_size = 16) +
  theme(
    axis.title = element_text(face = "bold",size=15),
    axis.text = element_text(color = "black"),
    plot.title = element_text(size = 15)
  )
