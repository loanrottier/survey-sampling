library(readr)
library(tidyverse)
library(ggplot2)
library(sampling)
library(survey)

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

ggplot(MU284, aes(x = P85)) +
  geom_density() +
  labs(title = paste("Densité de la population en 1985"),
       x = "Population", y = "Fréquence") +
  theme_minimal()

ggplot(MU284, aes(x = P75)) +
  geom_histogram(bins = 30, fill = "steelblue", color = "white") +
  labs(title = paste("Histogramme de la population en 1975"),
       x = "Population", y = "Fréquence") +
  theme_minimal()

ggplot(MU284, aes(x = P75, y = P85)) +
  geom_point(color = "steelblue", alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "red") +
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

ggplot(MU284, aes(x = P75, y = P85, color = factor(outlier))) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, color = "lightgreen") +
  labs(
    title = "Scatter plot entre P75 et P85",
    x = "Population en 1975 (P75)",
    y = "Population en 1985 (P85)"
  ) +
  theme_minimal()

### Création des strates
# Par quartiles

MU284$strate <- cut(MU284$P75,
                 breaks = quantile(MU284$P75, probs = c(0, 1/2, 28.5/30, 1)),
                 include.lowest = TRUE,
                 labels = c("basse", "moyenne", "haute"))

MU284 <- MU284 %>%
  mutate(strate_id = case_when(
    strate == "basse"   ~ 1,
    strate == "moyenne" ~ 2,
    strate == "haute"   ~ 3
  ))

ggplot(MU284, aes(x = P75, y = P85, color = factor(strate))) +
  geom_point() +
  labs(
    title = "Scatter plot entre P75 et P85",
    x = "Population en 1975 (P75)",
    y = "Population en 1985 (P85)"
  ) +
  theme_minimal()

### Estimation du total (dans la population)
# Population
P85_pop <- sum(MU284$P85)
P75_pop <- sum(MU284$P75)

#Sample
#Sampling design: all of the big cities
#some of the middle and small size
N_H <- table(MU284$strate)
# calcul n_h
n <- 70 #choix de n
table(MU284$strate)
std <- c()
for (i in 1:3){
  std[i] <- sqrt(var(MU284$P85[MU284$strate_id == i]))
}
tab <- table(MU284$strate)*std
sum_NL <- sum(tab)
n_H <- n * tab / sum_NL
round(n_H) #pas ok car n_3 = 26 > N_3 = 15

n_3 <- table(MU284$strate)['haute']
n_12 <- n - n_3

tab <- table(MU284$strate)[1:2]*std[1:2]
sum_NL <- sum(tab)
n_H <- n_12 * tab / sum_NL
n_H[3] <- n_3 
n_H <- round(n_H)

# sample
stsrswor = strata(MU284,"strate",size=c(n_H[2],n_H[1],n_H[3]),method="srswor")
stsrswor_data = getdata(MU284,stsrswor)

#y_ht
poids_stsi=c(rep(rep(N_H[2]/n_H[2],n_H[2]),N_H[1]/n_H[1],n_H[1]),rep(N_H[3]/n_H[3],n_H[3]))
ech.stsi=svydesign(id=~LABEL,strata=~strate,weights=poids_stsi,
                   fpc=c(rep(N_H[1], n_H[1]), rep(N_H[2], n_H[2]), rep(N_H[3], n_H[3])),
                   data=stsrswor_data)
res=svytotal(~P85,ech.stsi) #Warning mais je sais pas d'où ça vient
res

#Estimation from the sample
#Winsorization
#sum(min(y_i,R))
# R aléatoire
R = 200
stsrswor_data$P85_R <- min(stsrswor_data$P_85, R)

R_quantile = quantile(stsrswor_data$P85, probs = 0.90)
stsrswor_data$P85_R_quantile <- min(stsrswor_data$P_85, R_quantile) 

ech.stsi=svydesign(id=~LABEL,strata=~strate,weights=poids_stsi,
                   fpc=c(rep(N_H[1], n_H[1]), rep(N_H[2], n_H[2]), rep(N_H[3], n_H[3])),
                   data=stsrswor_data)

res_R <- svytotal(~P85_R, ech.stsi)
res_R_mean <- svymean(~P85_R, ech.stsi)
res_R_quantile_mean <- svymean(~P85_R_quantile, ech.stsi)

#R - µ/ n-1 = E(max(Y-R,0))
# On estime l'espérance puis on résoud l'équation : H(R) = R - µ/ n-1 - E(max(Y-R,0)) = 0
mu_pop <- mean(MU284$P85)
n <- length(stsrswor_data) 
H_R <- function(R, mu, n, sample){
  terme1 <- (R - mu) / (n - 1)
  terme2 <- sum(pmax(sample - R, 0))/n
  return(terme1 - terme2)
}
R_esp_max <- uniroot(f = H_R, interval = c(min(MU284$P85), max(MU284$P85)),
                     mu = mu_pop, n = n, sample = stsrswor_data$P85)$root

stsrswor_data$P85_R_esp_max <- min(stsrswor_data$P85, R_esp_max)

ech.stsi=svydesign(id=~LABEL,strata=~strate,weights=poids_stsi,
                   fpc=c(rep(N_H[1], n_H[1]), rep(N_H[2], n_H[2]), rep(N_H[3], n_H[3])),
                   data=stsrswor_data)

res_R_esp_max <- svytotal(~P85_R_esp_max, ech.stsi)
#pas optimal (résultat de l'article)
res_R_esp_max_mean <- svymean(~P85_R_esp_max, ech.stsi)

#R = 2ème plus rand y_i
#dans l'échantillon
P85_ordered <- sort(unique(stsrswor_data$P85),decreasing = TRUE)
P85_n <- P85_ordered[1]
P85_n1 <- P85_ordered[2]
y_i1 <- mean(stsrswor_data$P85) - (P85_n - P85_n1) / n
