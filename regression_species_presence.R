#regression model on species presence in conservation groups
library(ggplot2)
library(dplyr)
library(FSA)

df <- data.frame(
  Conservation = factor(rep(c("HC", "C", "MC", "LC"), each = 5), levels = c("LC", "MC", "C", "HC")),
  Species = rep(c("Azucena", "Nipponbare", "IR64", "Meridionalis", "Rufipogon"), times = 4),
  Presence = c(
    100.00, 100.00, 100.00, 100.00, 100.00,      # HC
    86.36, 91.76, 83.81, 53.69, 84.38,            # C
    66.99, 72.50, 67.98, 36.54, 55.99,            # MC
    42.15, 46.22, 47.28, 29.00, 35.35             # LC
  )
)
df$ConservationNum <- as.numeric(df$Conservation)

ggplot(df, aes(x = ConservationNum, y = Presence, color = Species)) +
  geom_point(size = 3) +
  geom_smooth(method = "loess", se = FALSE, linewidth = 1) +
  scale_x_continuous(breaks = 1:4, labels = c("LC", "MC", "C", "HC")) +
  labs(
    title = "Species Presence Across TAD Boundary Conservation Levels",
    x = "Conservation Category",
    y = "% of Boundary Groups Containing Species"
  ) +
  theme_minimal(base_size = 14)

###statistical analysis
#R^2 and whole p-value (differences in species)
lm_model <- lm(Presence ~ Conservation + Species, data = df)
summary(lm_model)

#p-value (differences in conservation categories)
kruskal_result <- kruskal.test(Percentage ~ as.factor(Conservation), data = df)

#r
cor(df$Conservation,df$Percentage,method="pearson")