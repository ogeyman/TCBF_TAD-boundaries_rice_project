#PCA plot of conserved TAD boundary data

library(googlesheets4)
gs4_auth()
all_cboundary_bitable <- read_sheet("1nYCwVnE4Fm0XRWNQuYxQAS7or4EjZD0o6VKaBO47dUc",sheet="Boundaries_Groups")

transposed <- t(all_cboundary_bitable[, -1]) %>% as.data.frame()
colnames(transposed) <- all_cboundary_bitable$Groups
rownames(transposed) <- colnames(all_cboundary_bitable)[-1]

transposed <- transposed[, apply(transposed, 2, function(x) var(as.numeric(x)) != 0)] #removes constant or 0 values
pca <- prcomp(transposed,scale=TRUE)

library(ggplot2)
pca_df <- as.data.frame(pca$x)
pca_df$species <- rownames(pca_df)

# Plot only species points
ggplot(pca_df, aes(x = PC1, y = PC2, label = species)) +
  geom_point(size = 5, color = "#0072B2") +
  geom_text(vjust = -1, size = 4) +
  labs(title = "PCA: Species Based on Conserved TAD Groups",
       x = paste0("PC1 (", round(summary(pca)$importance[2, 1] * 100, 1), "%)"),
       y = paste0("PC2 (", round(summary(pca)$importance[2, 4] * 100, 1), "%)")) +
  theme_minimal()

#3D PCA
library(plotly)

pca_df <- as.data.frame(pca$x)
pca_df$species <- rownames(pca_df)

plot_ly(pca_df, 
        x = ~PC1, y = ~PC2, z = ~PC3,
        type = "scatter3d", mode = "markers+text",
        text = ~species, textposition = "top center",
        marker = list(size = 6, color = "#1f77b4")) %>%
  layout(title = "3D PCA of Species Based on Conserved TAD Groups",
         scene = list(xaxis = list(title = paste0("PC1 (", round(summary(pca)$importance[2, 1] * 100, 1), "%)")),
                      yaxis = list(title = paste0("PC2 (", round(summary(pca)$importance[2, 2] * 100, 1), "%)")),
                      zaxis = list(title = paste0("PC3 (", round(summary(pca)$importance[2, 3] * 100, 1), "%)"))))

