library(readr)
library(tidyverse)
library(ggplot2)
library(sampling)

### MU284 ####
data(MU284)
?MU284
# Y = Population en 1985
# X = Population en 1975
# Création de strates : petites / Moyennes / Grandes villes
# Plus grande croissance que prévue, au contraire, fort déclin de la ville 
# => Changement de strates
ggplot(MU284, aes(x = P85)) +
  geom_histogram(bins = 30, fill = "steelblue", color = "white") +
  labs(title = paste("Histogramme de la population en 1985"),
       x = "Population", y = "Fréquence") +
  theme_minimal()

ggplot(MU284, aes(x = P75)) +
  geom_histogram(bins = 30, fill = "steelblue", color = "white") +
  labs(title = paste("Histogramme de la population en 1975"),
       x = "Population", y = "Fréquence") +
  theme_minimal()

ggplot(df, aes(x = P75, y = P85)) +
  geom_point(color = "steelblue", alpha = 0.7) +
  labs(
    title = "Scatter plot entre P75 et P85",
    x = "Population en 1975 (P75)",
    y = "Population en 1985 (P85)"
  ) +
  theme_minimal()

boxplot(MU284$P85, MU284$P75)

### Identification des outliers 
Q1  <- quantile(MU284$P85, 0.25)
Q3  <- quantile(MU284$P85, 0.75)
IQR <- Q3 - Q1

lower <- Q1 - 1.5 * IQR
upper <- Q3 + 1.5 * IQR

outliers <- MU284$LABEL[MU284$P85 < lower | MU284$P85 > upper]
MU284$outlier <- ifelse(seq_len(nrow(MU284)) %in% outliers, 1, 0)

### Création des strates
# Par quartiles

MU284$strate <- cut(MU284$P75,
                 breaks = quantile(MU284$P75, probs = c(0, 1/3, 2/3, 1)),
                 include.lowest = TRUE,
                 labels = c("basse", "moyenne", "haute"))

### Estimation du total (dans la population)
# P85
P85_pop <- sum(MU284$P85)
P75_pop <- sum(MU284$P75)
