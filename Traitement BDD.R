library(readr)
library(tidyverse)
library(ggplot2)

### Import of the dataset
housing <- read_csv("housing.csv")

### First Analysis of the dataset
columns <- setdiff(names(housing), c("longitude", "latitude", "ocean_proximity"))
for (col in columns) {
  print(col)
    p <- ggplot(housing, aes(x = .data[[col]])) +
      geom_histogram(bins = 30, fill = "steelblue", color = "white") +
      labs(title = paste("Histogramme de", col),
           x = col, y = "Fréquence") +
      theme_minimal()
    print(p)
    }


### Creation of the outliers

