### Survey sampling ######
library(tidyverse) 
library(ggplot2) 
library(sampling) 
library(survey)
library(robsurvey)

data(MU284)

Q1 <- quantile(MU284$P85, 0.25)
Q3  <- quantile(MU284$P85, 0.75) 
IQR <- Q3 - Q1
lower <- Q1 - 1.5 * IQR 
upper <- Q3 + 1.5 * IQR
outliers <- MU284$LABEL[MU284$P85 < lower | MU284$P85 > upper]
MU284$outlier <- ifelse(seq_len(nrow(MU284)) %in% outliers, 1, 0)

N <- nrow(MU284)
n <- 70

# Creation of the strata
MU284$strate <- cut(MU284$P75, breaks = quantile(MU284$P75, probs = c(0, 1/2, 28.5/30, 1)), include.lowest = TRUE, labels = c("basse", "moyenne", "haute"))
MU284 <- MU284 %>% mutate(strate_id = case_when(
  strate == "basse" ~ 1,
  strate == "moyenne" ~ 2,
  strate == "haute" ~ 3 ))

ggplot(MU284, aes(x = P75, y = P85, color = factor(strate))) + 
  geom_point() + 
  labs( title = "Scatter plot entre P75 et P85", x = "Population en 1975 (P75)", y = "Population en 1985 (P85)" ) + 
  theme_minimal()

# Neyman allocation 
N_H <- table(MU284$strate)
n <- 70 
table(MU284$strate) 
std <- c() 
for (i in 1:3){ 
  std[i] <- sqrt(var(MU284$P85[MU284$strate_id == i])) 
} 
tab <- table(MU284$strate)*std 
sum_NL <- sum(tab) 
n_H <- n * tab / sum_NL 
round(n_H) #pas ok car n_3 = 34 > N_3 = 15

n_3 <- table(MU284$strate)['haute'] 
n_12 <- n - n_3
tab <- table(MU284$strate)[1:2] * std[1:2] 
sum_NL <- sum(tab) 
n_H <- n_12 * tab / sum_NL 
n_H['haute'] <- n_3 
n_H <- round(n_H)
n_H

#### Mean - Winsorization ######
set.seed(1)

ech_index <- srswor(n, N)
echantillon_srswor <- MU284[ech_index == 1, ]
poids_srswor <- rep(N/n, n)
ech.srswor <- svydesign(id=~LABEL, weights=poids_srswor, fpc=rep(N, n), data=echantillon_srswor)

P85_mean_pop <- mean(MU284$P85)
#Horvitz thomson estimator for the mean (with standart deviation)
estimate_mean <- svymean(~P85, ech.srswor)

### Once-winsorized mean #####
## manually
P85_sample <- echantillon_srswor$P85
y_bar <- mean(P85_sample)
P85_sorted <- sort(P85_sample, decreasing = TRUE)
y_n <- P85_sorted[1]       # y_(n)
y_n_moins_1 <- P85_sorted[2] # y_(n-1)

y_bar_1_manuel <- y_bar - (y_n - y_n_moins_1) / n
y_bar_1_manuel
## robsurvey package
svymean_1_winsorized <- svymean_k_winsorized(~P85, ech.srswor, k = 1)

### Winsorization #####
## Try different k
means_k_winso <- numeric(15)
for (k in 1:15){
  means_k_winso[k] <- coef(svymean_k_winsorized(~P85, ech.srswor, k = k))
}
plot(means_k_winso, ylim = c(20,30), xlab = "Value of k", ylab = "mean")
abline(h = P85_mean_pop, col = "red", lwd = 2)

## Optimal R (theory)
H_R <- function(R, mu, n, sample){
  terme1 <- (R - mu) / (n - 1)
  terme2 <- sum(pmax(sample - R, 0))/n
  return(terme1 - terme2)
}
R_esp_max <- uniroot(f = H_R, interval = c(min(MU284$P85), max(MU284$P85)),
                     mu = P85_mean_pop, n = n, sample = echantillon_srswor$P85)$root
MU284 %>% arrange(desc(P85)) %>% select(P85) %>% head(10) #compare R with the value
# k = 2

## Searl approach
quantiles <- seq(0.9, 0.99, by = 0.01)
mse_values <- numeric(length(quantiles))
for (i in seq_along(quantiles)) {
  UB <- quantiles[i]
  winsorized_mean <- svymean_winsorized(~P85, ech.srswor, LB = 0, UB = UB)
  mse_values[i] <- mse(winsorized_mean)
}
mse_values
best_quantile <- quantiles[which.min(mse_values)]
# best quantile : 0.98

### With Monte Carlo #####
nsimu <- 1000
estimations_HT <- numeric(nsimu)
estimations_k1 <- numeric(nsimu) 
estimations_UB_fixed <- numeric(nsimu)

for (i in 1:nsimu) {
  
  ech_index <- srswor(n, N)
  sample_si <- MU284[ech_index == 1, ]
  
  weight_si <- rep(N/n, n)
  ech.si <- svydesign(id=~LABEL, weights=weight_si, fpc=rep(N, n), data=sample_si)
  #for comparison with H-T estimator
  mean_HT <- svymean(~P85, ech.si)
  estimations_HT[i] <- coef(mean_HT)
  
  # winsorized mean with k = 1
  mean_k1 <- svymean_k_winsorized(~P85, ech.si, k = 1)
  estimations_k1[i] <- coef(mean_k1) 
  
  # Compute the mean with searl's fixed UB 
  mean_UB_fixed <- svymean_winsorized(~P85, ech.si, LB = 0, UB = 0.97)
  estimations_UB_fixed[i] <- coef(mean_UB_fixed) 
  
}

cat("with rivest method (k=1)")
bias_k1 <- mean(estimations_k1) - P85_mean_pop
var_k1 <- var(estimations_k1)
mse_k1 <- bias_k1^2 + var_k1

cat("with searl method (fixed UB)")
bias_UB_fixed <- mean(estimations_UB_fixed) - P85_mean_pop
var_UB_fixed <- var(estimations_UB_fixed)
mse_UB_fixed <- bias_UB_fixed^2 + var_UB_fixed

cat("comparison with Horvitz-Thomson estimator")
bias_HT <- mean(estimations_HT) - P85_mean_pop
var_HT <- var(estimations_HT)
mse_HT <- bias_HT^2 + var_HT

### Total - Winsorization ######
true_total <- sum(MU284$P85)

results_ht <- numeric(nsimu)
results_k1 <- numeric(nsimu)
results_ub <- numeric(nsimu)

for (i in 1:nsimu) {
  
  s <- srswor(n, N)
  sample_data <- MU284[s == 1, ]
  
  sample_data$weights <- N / n
  sample_data$fpc <- N
  
  design_obj <- svydesign(
    id = ~1,
    weights = ~weights,
    fpc = ~fpc,
    data = sample_data
  )
  
  results_ht[i] <- coef(svytotal(~P85, design_obj))
  results_k1[i] <- coef(svytotal_k_winsorized(~P85, design_obj, k = 1))
  results_ub[i] <- coef(svytotal_winsorized(~P85, design_obj, LB = 0, UB = 0.97))
}

mse_ht <- mean((results_ht - true_total)^2)
mse_k1 <- mean((results_k1 - true_total)^2)
mse_ub <- mean((results_ub - true_total)^2)

var_ht <- var(results_ht)
var_k1 <- var(results_k1)
var_ub <- var(results_ub)

bias_ht <- mean(results_ht) - true_total
bias_k1 <- mean(results_k1) - true_total
bias_ub <- mean(results_ub) - true_total

final_results <- data.frame(
  Estimator = c("Horvitz-Thompson", "Winsorized (k=1)", "Winsorized (UB=0.95)"),
  MSE = c(mse_ht, mse_k1, mse_ub),
  Variance = c(var_ht, var_k1, var_ub),
  Bias = c(bias_ht, bias_k1, bias_ub)
)

print(final_results)

## Total with STSRSWOR
estimations_HT_strat <- numeric(nsimu)
estimations_k1_strat <- numeric(nsimu)
estimations_UB_fixed_strat <- numeric(nsimu)
for (i in 1:nsimu) {
  sample_si_list <- list()
  
  data_basse <- MU284[MU284$strate == "basse", ]
  N_basse <- nrow(data_basse)
  n_basse <- n_H["basse"]
  ech_index_basse <- srswor(n_basse, N_basse)
  sample_basse <- data_basse[ech_index_basse == 1, ]
  
  data_moyenne <- MU284[MU284$strate == "moyenne", ]
  N_moyenne <- nrow(data_moyenne)
  n_moyenne <- n_H["moyenne"]
  ech_index_moyenne <- srswor(n_moyenne, N_moyenne)
  sample_moyenne <- data_moyenne[ech_index_moyenne == 1, ]
  
  #take all values in strat haute
  data_haute <- MU284[MU284$strate == "haute", ]
  
  
  sample_si <-rbind(sample_basse, sample_moyenne, data_haute)
  weight_si <- c(rep(N_basse/n_basse, n_basse), rep(N_moyenne/n_moyenne, n_moyenne), rep(1, n_H["haute"]))
  ech.si <- svydesign(id=~LABEL, strata=~strate, weights=weight_si, fpc=c(rep(N_basse, n_basse), rep(N_moyenne, n_moyenne), rep(N_H["haute"], n_H["haute"])), data=sample_si)
  
  #for comparison with H-T estimator
  res_ht <- svyby(~P85, ~strate, ech.si, svytotal)
  estimations_HT_strat[i] <- sum(res_ht$P85)
  
  
  # winsorized total with k = 1 in each stratum
  res_k1 <- svyby(~P85, ~strate, ech.si, svytotal_k_winsorized, k = 5)
  total_k1_mixte <- res_k1$P85[res_k1$strate == "basse"] + 
    res_k1$P85[res_k1$strate == "moyenne"] + 
    res_ht$P85[res_ht$strate == "haute"]
  
  estimations_k1_strat[i] <- total_k1_mixte
  
  # Compute the total with searl's fixed UB
  res_ub <- svyby(~P85, ~strate, ech.si, svytotal_winsorized, LB = 0, UB = 0.8)
  total_ub_mixte <- res_ub$P85[res_ub$strate == "basse"] + 
    res_ub$P85[res_ub$strate == "moyenne"] + 
    res_ht$P85[res_ht$strate == "haute"]
  estimations_UB_fixed_strat[i] <- total_ub_mixte
}

Y_total <- sum(MU284$P85)
Y_HT <- mean(estimations_HT_strat)
Y_UB_fixed_total <- mean(estimations_UB_fixed_strat)
Y_searl_total <- mean(estimations_k1_strat)

cat("compare mse")
mse_ht <- mean((estimations_HT_strat - Y_total)^2)
mse_k1 <- mean((estimations_k1_strat - Y_total)^2)
mse_ub <- mean((estimations_UB_fixed_strat - Y_total)^2)

cat("compare bias")
bias_ht <- mean(estimations_HT_strat) - Y_total
bias_k1 <- mean(estimations_k1_strat) - Y_total
bias_ub <- mean(estimations_UB_fixed_strat) - Y_total

cat("compare variance")
var_ht <- var(estimations_HT_strat)
var_k1 <- var(estimations_k1_strat)
var_ub <- var(estimations_UB_fixed_strat)

# variance and mse are higher with the winsorization methods than with the H-T estimator.
# due to the fact that we dont have outlier in our sample

## With outliers in each stratum
# Creation of the outliers dataset
MU284_OUTLIERS <- MU284

indices_basse <- which(MU284_OUTLIERS$strate == "basse")

#change the first 2 indices in strat basse and multiply it by 4 and 8
MU284_OUTLIERS$P85[indices_basse[1]] <- MU284_OUTLIERS$P85[indices_basse[1]] * 4
MU284_OUTLIERS$P85[indices_basse[2]] <- MU284_OUTLIERS$P85[indices_basse[2]] * 8

indices_moyenne <- which(MU284_OUTLIERS$strate == "moyenne")
#change the first 2 indices in strat moyenne and multiply it by 4 and 8
MU284_OUTLIERS$P85[indices_moyenne[1]] <- MU284_OUTLIERS$P85[indices_moyenne[1]] * 5
MU284_OUTLIERS$P85[indices_moyenne[2]] <- MU284_OUTLIERS$P85[indices_moyenne[2]] * 4

print("see the data with outliers in comparaison without outliers")
by(MU284_OUTLIERS$P85, MU284_OUTLIERS$strate, summary)

# Winsorization
set.seed(2)
nsimu <- 1000 
estimations_HT_strat <- numeric(nsimu)
estimations_k1_strat <- numeric(nsimu)
estimations_UB_fixed_strat <- numeric(nsimu)
for (i in 1:nsimu) {
  sample_si_list <- list()
  
  data_basse <- MU284_OUTLIERS[MU284_OUTLIERS$strate == "basse", ]
  N_basse <- nrow(data_basse)
  n_basse <- n_H["basse"]
  ech_index_basse <- srswor(n_basse, N_basse)
  sample_basse <- data_basse[ech_index_basse == 1, ]
  
  data_moyenne <- MU284_OUTLIERS[MU284_OUTLIERS$strate == "moyenne", ]
  N_moyenne <- nrow(data_moyenne)
  n_moyenne <- n_H["moyenne"]
  ech_index_moyenne <- srswor(n_moyenne, N_moyenne)
  sample_moyenne <- data_moyenne[ech_index_moyenne == 1, ]
  
  #take all values in strat haute
  data_haute <- MU284_OUTLIERS[MU284_OUTLIERS$strate == "haute", ]
  
  
  sample_si <-rbind(sample_basse, sample_moyenne, data_haute)
  weight_si <- c(rep(N_basse/n_basse, n_basse), rep(N_moyenne/n_moyenne, n_moyenne), rep(1, n_H["haute"]))
  ech.si <- svydesign(id=~LABEL, strata=~strate, weights=weight_si, fpc=c(rep(N_basse, n_basse), rep(N_moyenne, n_moyenne), rep(N_H["haute"], n_H["haute"])), data=sample_si)
  
  #for comparison with H-T estimator
  res_ht <- svyby(~P85, ~strate, ech.si, svytotal)
  estimations_HT_strat[i] <- sum(res_ht$P85)
  
  
  # winsorized total with k = 1 in each stratum
  res_k1 <- svyby(~P85, ~strate, ech.si, svytotal_k_winsorized, k = 1)
  total_k1_mixte <- res_k1$P85[res_k1$strate == "basse"] + 
    res_k1$P85[res_k1$strate == "moyenne"] + 
    res_ht$P85[res_ht$strate == "haute"]
  
  estimations_k1_strat[i] <- total_k1_mixte
  
  # Compute the total with searl's fixed UB
  res_ub <- svyby(~P85, ~strate, ech.si, svytotal_winsorized, LB = 0, UB = 0.9)
  total_ub_mixte <- res_ub$P85[res_ub$strate == "basse"] + 
    res_ub$P85[res_ub$strate == "moyenne"] + 
    res_ht$P85[res_ht$strate == "haute"]
  estimations_UB_fixed_strat[i] <- total_ub_mixte
}

Y_total <- sum(MU284_OUTLIERS$P85)
Y_HT <- mean(estimations_HT_strat)
Y_UB_fixed_total <- mean(estimations_UB_fixed_strat)
Y_searl_total <- mean(estimations_k1_strat)

cat("compare bias")
bias_ht <- mean(estimations_HT_strat) - Y_total
bias_k1 <- mean(estimations_k1_strat) - Y_total
bias_ub <- mean(estimations_UB_fixed_strat) - Y_total

cat("compare variance")
var_ht <- var(estimations_HT_strat)
var_k1 <- var(estimations_k1_strat)
var_ub <- var(estimations_UB_fixed_strat)

cat("compare mse")
mse_ht <- mean((estimations_HT_strat - Y_total)^2)
mse_k1 <- mean((estimations_k1_strat - Y_total)^2)
mse_ub <- mean((estimations_UB_fixed_strat - Y_total)^2)

### Auxiliary information ######
set.seed(12)
sum_x <- numeric(nsimu)
calibrated_x <- numeric(nsimu)
estimation_y <- numeric(nsimu)
estimation_y_weighted <- numeric(nsimu)

for (i in 1:nsimu) {
  
  ech_index <- srswor(n, N)
  sample_si <- MU284[ech_index == 1, ]
  
  sample_si$weight_si <- rep(N/n, n)
  ech.si <- svydesign(id=~LABEL, weights=~weight_si, fpc=rep(N, n), data=sample_si)
  
  
  #calibrated estimator on x
  sum_x[i] <- sum(sample_si$P75) * N / n
  calibrated_x[i]<- svytotal(~P75, ech.si)
  #estimation of y
  estimation_y[i] <- svytotal(~P85, ech.si)
  estimation_y_weighted[i] <- weighted_total(sample_si$P85, sample_si$weight_si)
  
}

print(sum(MU284$P75))
print(mean(sum_x))
print(mean(calibrated_x)) #iso to the previous line
#the estimation is ok then we can use the same weight to compute the T_y
print(sum(MU284$P85))
print(mean(estimation_y))
print(mean(estimation_y_weighted))

# with weight corrected by the ratio on the auxiliary variable
estimation_total_corrected_weigths <- numeric(nsimu)
mean <- mean(sum_x)
true_x <- sum(MU284$P75)

for(i in 1:nsimu){
  ech_index <- srswor(n, N)
  sample_si <- MU284[ech_index == 1, ]
  
  sample_si$weight_si <- rep(N/n, n) /  mean * true_x
  ech.si <- svydesign(id=~LABEL, weights=~weight_si, fpc=rep(N, n), data=sample_si)
  
  estimation_total_corrected_weigths[i] <- weighted_total(sample_si$P85, sample_si$weight_si)
}

print(mean(estimation_y_cor_weight))

# model
model <- lm(P85 ~ P75, data = sample_si)
beta <- coef(model)['P75']
plot(model$residuals)
sample_si$residuals <- model$residuals
sample_si$weight_model <- sample_si$weight_si
sample_si[sample_si$residuals > 10,]$weight_model <- sample_si[sample_si$residuals > 10,]$weight_model * 0.9
weighted_total(sample_si$P85, sample_si$weight_si)
weighted_total(sample_si$P85, sample_si$weight_model)

### Stratum jumpers ######

# Creation of the dataset (adding stratum jumpers)
nb_stratum_jumpers <- 1
MU284_jump <- MU284 
index <- sample(MU284_jump[MU284_jump$strate == "haute",]$LABEL, nb_stratum_jumpers)
MU284_jump[MU284$LABEL %in% index,]$strate <- "basse" 
MU284_jump[MU284$LABEL %in% index,]$strate_id <- "1" 
MU284_jump[MU284$LABEL %in% index,]$P75 <- 10 #mean of P75 in the "basse" strata 
MU284_jump$stratum_jumper <- ifelse(MU284_jump$LABEL %in% index, 1, 0)

ggplot(MU284_jump, aes(x = P75, y = P85, color = factor(strate))) + 
  geom_point() + 
  labs( title = "Scatter plot entre P75 et P85", x = "Population en 1975 (P75)", y = "Population en 1985 (P85)" ) + 
  theme_minimal()

# Compute the allocation
N_H <- table(MU284_jump$strate)
# calcul n_h
n <- 70 #choix de n
table(MU284_jump$strate) 
std <- c() 
for (i in 1:3){ 
  std[i] <- sqrt(var(MU284_jump$P85[MU284_jump$strate_id == i])) 
} 
tab <- table(MU284_jump$strate)*std 
sum_NL <- sum(tab) 
n_H <- n * tab / sum_NL 
round(n_H) #not ok car n_3 > N_3

# compute again the allocation by fixing n_3
n_3 <- table(MU284_jump$strate)['haute'] 
n_12 <- n - n_3
tab <- table(MU284_jump$strate)[1:2] * std[1:2] 
sum_NL <- sum(tab) 
n_H <- n_12 * tab / sum_NL 
n_H['haute'] <- n_3 
n_H <- round(n_H)
n_H

# draw the sample
set.seed(12)
ech1 <- sample(1:nrow(MU284_jump[MU284_jump$strate=="basse",]), size = n_H["basse"])
ech2 <- sample(1:nrow(MU284_jump[MU284_jump$strate=="moyenne",]), size = n_H["moyenne"])
ech3 <- sample(1:nrow(MU284_jump[MU284_jump$strate=="haute",]), size = n_H["haute"])
sample1 <- MU284_jump[MU284_jump$strate=="basse",]
sample1 <- sample1[ech1,]
sample1 <- sample1[1:(nrow(sample1)-nb_stratum_jumpers),]
sample1 <- rbind(sample1, MU284_jump[MU284_jump$stratum_jumper == 1,])
sample2 <- MU284_jump[MU284_jump$strate=="moyenne",]
sample2 <- sample2[ech2,]
sample3 <- MU284_jump[MU284_jump$strate=="haute",]
sample3 <- sample3[ech3,]
stsrswor_data_jump <- rbind(sample1, sample2, sample3) 

w1 <- nrow(MU284_jump[MU284_jump$strate == "basse",]) / nrow(stsrswor_data_jump[stsrswor_data_jump$strate == "basse",])
w2 <- nrow(MU284_jump[MU284_jump$strate == "moyenne",]) / nrow(stsrswor_data_jump[stsrswor_data_jump$strate == "moyenne",])
w3 <- nrow(MU284_jump[MU284_jump$strate == "haute",]) / nrow(stsrswor_data_jump[stsrswor_data_jump$strate == "haute",])

poids_stsi=c(rep(w1,n_H[1]),rep(w2,n_H[2]),rep(w3,n_H[3]))
stsrswor_data_jump$poids <- poids_stsi
ech.stsi=svydesign(id=~LABEL,strata=~strate,weights=poids_stsi, fpc=c(rep(N_H[1], n_H[1]), rep(N_H[2], n_H[2]), rep(N_H[3], n_H[3])), data=stsrswor_data_jump) 

svytotal(~P85, ech.stsi) # SJ have higher rate than it should

# changing weights
stsrswor_data_jump$poids_corrige <- poids_stsi
stsrswor_data_jump[stsrswor_data_jump$stratum_jumper == 1,]$poids_corrige <- w3

stsrswor_data_jump$collection_stratum <- stsrswor_data_jump$strate
stsrswor_data_jump[stsrswor_data_jump$stratum_jumper == 1,]$collection_stratum <- "haute"
moyennes <- stsrswor_data_jump %>%
  group_by(collection_stratum) %>%
  summarise(poids_moyen = mean(poids)) %>%
  ungroup()
stsrswor_data_jump <- stsrswor_data_jump %>%
  left_join(moyennes, by = "collection_stratum")

stsrswor_data_jump$poids_moyen_SJ <- stsrswor_data_jump$poids
stsrswor_data_jump[stsrswor_data_jump$stratum_jumper == 1,]$poids_moyen_SJ <- stsrswor_data_jump[stsrswor_data_jump$stratum_jumper == 1,]$poids_moyen
ratio <- sum(stsrswor_data_jump$poids_moyen) / sum(stsrswor_data_jump$poids_moyen_SJ)
stsrswor_data_jump$poids_ratio <- stsrswor_data_jump$poids_moyen_SJ * ratio  

weighted_total(stsrswor_data_jump$P85, stsrswor_data_jump$poids)
weighted_total(stsrswor_data_jump$P85, stsrswor_data_jump$poids_corrige)
weighted_total(stsrswor_data_jump$P85, stsrswor_data_jump$poids_moyen)
weighted_total(stsrswor_data_jump$P85, stsrswor_data_jump$poids_ratio)
