##############################################################################################################################################
#### Spatial and Racial/ethnic Disparities in Cardiovascular Mortality Attributable to PM2.5 Components in the Contiguous United States   ####
#### Cleaned R code for the Sensitivity analyses: non-linear relationship                                                                 ####
#### Ying Hu, Lingzhi Chu, Stefano Renzetti, Emma Zang, Ijeoma Opara, Yuan Lu, Erica S. Spatz, Harlan M. Krumholz, Kai Chen               ####
#### Yale University, New Haven, CT                                                                                                       ####
#### Feb. 20, 2025                                                                                                                        ####
##############################################################################################################################################

rm(list = ls())
graphics.off()
cat("\014")

library(gWQS)
library(fixest)
library(dplyr)
library(splines)
library(multcomp)
library(parallel)
library(future)
library(future.apply)

lag <- 11
start_year <- 2001
end_year <- 2020

# gWQS_feglm settings
usecores <- 30
rh <- 100
b <- 10
nq <- 10
validation <- 0.5
lambda_value <- 10
seed <- 1000
mix_name <- c('BC', 'NH4', 'NIT', 'OM', 'SO4', 'SOIL')
covar_name <- c("temp1", "temp2", "temp3", "temp4", "temp5")
fix_name <- c("GEOID", "year_month")
cluster_name <- "GEOID"

# file path
healthdata.path <- "file path"
expdata.path <- "file path"
code.path <- "file path"
save.path <- "file path"

# Load data and function
exposure.data <- readRDS(paste0(expdata.path, "exposure_data.rds"))
health.data <- readRDS(paste0(healthdata.path, "health_data.rds"))
all.data <- merge(health.data, exposure.data, c("GEOID", "year_month"), all = TRUE)
source(paste0(code.path, paste0("gWQS_feglm.R")))

# Calculate lagged exposure
for (pollutant in c ("BC", "NH4", "NIT", "OM", "SO4", "SOIL")) {
  pollutant.sub <- all.data[, c(pollutant, paste0(pollutant, "_lag", 1:lag))]
  assign(pollutant, base::rowMeans(pollutant.sub))
}
temp.mean <- all.data[, c("tmean")]
death.rate <- all.data$monthly.rate.adj

# Clean data
model.data <- data.frame(y = death.rate,
                         BC = BC,
                         NH4 = NH4,
                         NIT = NIT,
                         OM = OM,
                         SO4 = SO4,
                         SOIL = SOIL,
                         temp1 = ns(temp.mean,5)[,1],
                         temp2 = ns(temp.mean,5)[,2],
                         temp3 = ns(temp.mean,5)[,3],
                         temp4 = ns(temp.mean,5)[,4],
                         temp5 = ns(temp.mean,5)[,5],
                         GEOID = paste("id", all.data$GEOID),
                         year = all.data$Year,
                         year_month = all.data$year_month)
model.data <- na.omit(model.data)
model.data <- model.data %>%
  filter(year >= start_year & year <= end_year)

weight_vector <- as.vector(as.numeric(t(weight_vector)))
data.fit <- model.data[ ,fix_name]
data.fit[,covar_name] <- model.data[, covar_name]
data.fit$y <- model.data$y

q_f <- gwqs_rank(model.data, mix_name, nq)
bQ <- q_f$Q
cr.data <- readRDS(paste0(save.path, "result_file_name.rds"))
weights <- as.matrix(cr.data[["Output"]][["p"]]["Mean", 1:length(mix_name)])
weight_vector <- as.vector(as.numeric(t(weights)))
data.fit[, "wqs"] <- as.numeric(as.matrix(bQ) %*% weight_vector)
df_value <- 5
feglm_model <- feglm(y ~ ns(wqs, df = df_value) + temp1 + temp2 + temp3 + temp4 + temp5 | year_month + GEOID, 
                     family = "gaussian",
                     cluster = cluster_name,
                     data = data.fit)

# Extract fixed effects from feglm model
fix_ef <- fixef(feglm_model)

# Copy valid_data and adjust for fixed effects
model.data.nfix <- data.fit
for (fe_var in c("year_month", "GEOID")) {
  fe_levels <- unique(model.data.nfix[[fe_var]])
  for (level in fe_levels) {
    group_indices <- which(model.data.nfix[[fe_var]] == level)
    model.data.nfix[group_indices, "y"] <- model.data.nfix[group_indices, "y"] - fix_ef[[fe_var]][[level]]
  }
}

# Fit a standard GLM without fixed effects
formula2 <- as.formula(paste("y ~ ns(wqs, df = ", df_value, ") +", paste(c("temp1", "temp2", "temp3", "temp4", "temp5"), collapse = " + ")))
result.glm <- glm(formula2, data = model.data.nfix, family = "gaussian")

# Extract coefficients and standard errors from both models
glm.coef <- coef(result.glm)
feglm.coef <- coef(feglm_model)
coef_names <- intersect(names(glm.coef), names(feglm.coef))
glm.se <- sqrt(diag(vcov(result.glm)))[coef_names]
feglm.se <- sqrt(diag(vcov(feglm_model)))[coef_names]
scale_factors <- feglm.se / glm.se

# Prepare data for visualization
term <- paste0("ns(wqs, df = ", df_value, ")")
plot_data <- termplot(result.glm, terms = term, partial.resid = TRUE, se = TRUE, plot = FALSE)
plot_data[[1]]$y <- (plot_data[[1]]$y - min(plot_data[[1]]$y)) * 1e6
plot_data[[1]]$se <- plot_data[[1]]$se * 1e6 

# Interpolation of scaling factors
x_values <- seq(min(plot_data[[1]]$x, na.rm = TRUE), max(plot_data[[1]]$x, na.rm = TRUE), length.out = df_value)
interpolated_factors <- approx(x = x_values, y = scale_factors[1:df_value], xout = plot_data[[1]]$x)$y
adjusted_se <- plot_data[[1]]$se * interpolated_factors
plot_data[[1]]$se_adj <- adjusted_se

# Compute confidence intervals
RR.L <- plot_data[[1]]$y - 1.96 * adjusted_se
RR.U <- plot_data[[1]]$y + 1.96 * adjusted_se


plot(plot_data[[1]]$x, plot_data[[1]]$y, type = "l", col = "black", lwd = 2,
     xlab = "WQS", ylab = "Death Rate")
lines(plot_data[[1]]$x, RR.U, col = "darkgrey", lty = 2)
lines(plot_data[[1]]$x, RR.L, col = "darkgrey", lty = 2)