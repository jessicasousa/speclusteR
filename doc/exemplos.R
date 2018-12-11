## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ------------------------------------------------------------------------
library(ggplot2)
set.seed(2018)
n <- 150
r <- rnorm(n, 5, .25)
theta <- runif(n, 0, 2 * pi)
c1 <- data.frame(x = rnorm(n), y = rnorm(n))
c2 <- data.frame(x = r * cos(theta), y = r * sin(theta))
points <- rbind(c1, c2)
km <- stats::kmeans(points, centers = 2)
cluster <- factor(km$cluster)
ggplot(points, aes(x = x, y = y, color = cluster)) +
geom_point()

## ------------------------------------------------------------------------
pointsm <- data.matrix(points)
specc <- speclusteR::apply_spectral_clustering(pointsm, k = 2,  type = 3)
cluster <- factor(specc$cluster)
#plot
ggplot(points, aes(x = x, y = y, color = cluster)) +
geom_point()

## ---- fig.align='center', fig.show='hold'--------------------------------
require(gridExtra)

dados <- speclusteR::spec_data
espirais <- dados[['espirais']]

df <- as.data.frame(espirais)
names(df) <- c('x','y')

specc <- speclusteR::apply_spectral_clustering(espirais, k = 3,  type = 3)
km <- stats::kmeans(df, centers = 3)
cluster_spec <-  factor(specc$cluster)
cluster_kmeans <- factor(km$cluster)
#plot
plot1 <- ggplot(df, aes(x = x, y = y, color = cluster_spec)) +
         geom_point()
plot2 <- ggplot(df, aes(x = x, y = y, color = cluster_kmeans)) +
         geom_point()

grid.arrange(plot1, plot2, nrow=2)

