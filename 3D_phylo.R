#3D-Genomic Phylogenetic Tree

library(vegan)
library(ape)
library(googlesheets4)
library(ggtree)
library(ggplot2)
TADb_data <- read_sheet("1nYCwVnE4Fm0XRWNQuYxQAS7or4EjZD0o6VKaBO47dUc",sheet="Boundaries_Groups")
binary_data <- ifelse(TADb_data >= 1, 1, 0)
binary_data <- binary_data[,-1]

jac_dist <- vegdist(t(binary_data), method="jaccard",binary=TRUE)
jac_matrix <- as.matrix(jac_dist)

#unrooted
tree_nj <- nj(jac_matrix)
plot(tree_nj)

#rooted with O. meridionalis
tree_rooted <- root(tree_nj,outgroup="meridionalis",resolve.root=TRUE)
plot(tree_rooted)

#bootstrapping
boot_vals <- boot.phylo(
  tree_nj,
  binary_data,
  FUN = function(x) nj(vegdist(t(x),method="jaccard",binary=TRUE)),
  B=1000
)
boot_percent_clean <- round(boot_percent, 1)

plot(tree_nj)
nodelabels(boot_percent_clean,cex=0.8)

#plotting using ggtree (rooted version)
gt <- ggtree(tree_rooted)
boot_df <- data.frame(
  node = c(6, 7, 8),
  bootstrap = c(NA, 70.1, 100)
)
gt$data <- dplyr::left_join(gt$data, boot_df, by = "node")

final_plot <- gt +
  geom_tiplab(size = 5, fontface = "italic") +
  geom_label2(
    aes(label = bootstrap),
    fill = "#C7E9F1",
    color = "black",
    size = 4,
    fontface = "bold",
    label.padding = unit(0.2, "lines"),
    label.size = 0.3,
    na.rm = TRUE
  ) +
  ggtitle("3D-Genomic Phylogeny of Oryza Species") +
  theme_tree2() +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0),
    text = element_text(family = "Helvetica")
  )
scaled_plot <- final_plot +
  coord_cartesian(xlim = c(0, 0.45))  # adjust 0.35 as needed
scaled_plot

