##############################################################################################################################################
#### Spatial and Racial/ethnic Disparities in Cardiovascular Mortality Attributable to PM2.5 Components in the Contiguous United States   ####
#### Cleaned R code for the Sensitivity analyses: Adjust for NO2, O3, and dewT (2001 - 2016)                                              ####
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
end_year <- 2016

# gWQS_feglm settings
usecores <- 30
rh <- 100
b <- 10
nq <- 10
validation <- 0.5
lambda_value <- 10
seed <- 1000
mix_name <- c('BC', 'NH4', 'NIT', 'OM', 'SO4', 'SOIL')
covar_name <- c("temp1", "temp2", "temp3", "temp4", "temp5", "no2", "o3", "td")
fix_name <- c("GEOID", "year_month")
cluster_name <- "GEOID"

# file path
healthdata.path <- "file path"
expdata.path <- "file path"
code.path <- "file path"
save.path <- "file path"

# Load data and function
exposure.data <- readRDS(paste0(expdata.path, "exposure_data.rds"))
no2 <- readRDS(paste0(expdata.path, "no2.rds"))
o3 <- readRDS(paste0(expdata.path, "o3.rds"))
dewT <- readRDS(paste0(expdata.path, "dewT.rds"))
exposure.data <- merge(exposure.data, no2, c("GEOID", "Year", "Month"), all = TRUE)
exposure.data <- merge(exposure.data, o3, c("GEOID", "Year", "Month"), all = TRUE)
exposure.data <- merge(exposure.data, dewT, c("GEOID", "Year", "Month", "year_month"), all = TRUE)
health.data <- readRDS(paste0(healthdata.path, "health_data.rds"))
all.data <- merge(health.data, exposure.data, c("GEOID", "year_month"), all = TRUE)
source(paste0(code.path, paste0("gWQS_feglm.R")))

# Calculate lagged exposure
for (pollutant in c ("BC", "NH4", "NIT", "OM", "SO4", "SOIL", "no2", "o3")) {
  pollutant.sub <- all.data[, c(pollutant, paste0(pollutant, "_lag", 1:lag))]
  assign(pollutant, base::rowMeans(pollutant.sub))
}
temp.mean <- all.data[, c("tmean")]
td.mean <- all.data[, c("tdmean")]
death.rate <- all.data$monthly.rate.adj

# Clean data
model.data <- data.frame(y = death.rate,
                         no2= no2,
                         o3 = o3,
                         td = td.mean,
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

time.start <- Sys.time()
result <- data.frame()
result <- gWQS_feglm(model.data = model.data,
                     usecores = usecores, rh = rh, b = b, nq = nq, validation = validation, lambda_value = lambda_value, seed = seed,
                     mix_name = mix_name, covar_name = covar_name, fix_name = fix_name, cluster_name = cluster_name)
print(result$Output$p)
print(Sys.time()-time.start)

saveRDS(result, file = paste0(save.path, "result_sa_no2o3td.rds"))