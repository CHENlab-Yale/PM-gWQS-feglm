##############################################################################################################################################
#### Spatial and Racial/ethnic Disparities in Cardiovascular Mortality Attributable to PM2.5 Components in the Contiguous United States   ####
#### Cleaned R code for the burden of disease estimation                                                                                  ####
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
library(ggplot2)
library(sf)
library(tigris)
library(patchwork)

lag <- 11
start_year <- 2001
end_year <- 2020
trend_years <- c(2001,2020)
mix_name <- c('BC', 'NH4', 'NIT', 'OM', 'SO4', 'SOIL')

shp.path <- "file_path"
healthdata.path <- "file_path"
expdata.path <- "file_path"
save.path <- "file_path"

shp.data <- st_read(shp.path)
states <- states(cb = TRUE)
nation <- nation()
states <- st_transform(states, crs = 4269)
nation <- st_transform(nation, crs = 4269)
bbox <- st_bbox(c(xmin = -125, ymin = 20, xmax = -66, ymax = 50), crs = 4269)
contiguous_states <- st_crop(states, bbox)
contiguous_nation <- st_crop(nation, bbox)

# Load data
exposure.data <- readRDS(paste0(expdata.path, "exposure_data.rds"))
pop.data <- readRDS(paste0(expdata.path, "health_data.rds"))[,c("GEOID", "year_month","pop.Total")]
all.data <- merge(pop.data, exposure.data, c("GEOID", "year_month"))

for (pollutant in c ("BC", "NH4", "NIT", "OM", "SO4", "SOIL", "PM25")) {
  pollutant.sub <- all.data[, c(pollutant, paste0(pollutant, "_lag", 1:lag))]
  assign(pollutant, base::rowMeans(pollutant.sub))
}
model.data <- data.frame(PM25 = PM25,
                         BC = BC,
                         NH4 = NH4,
                         NIT = NIT,
                         OM = OM,
                         SO4 = SO4,
                         SOIL = SOIL,
                         GEOID = all.data$GEOID,
                         state = as.numeric(substr(all.data$GEOID,1,2)),
                         year = all.data$Year,
                         month = all.data$Month,
                         Pop = all.data$pop.Total,
                         year_month = all.data$year_month)
model.data <- na.omit(model.data)
model.data <- model.data %>%
  filter(year >= start_year & year <= end_year)


cr.data <- readRDS(paste0(save.path, "result_file_name.rds"))
weights <- as.matrix(cr.data[["Output"]][["p"]]["Mean", 1:length(mix_name)])
colnames(weights) <- mix_name
coef.mean <- cr.data[["Output"]][["p"]]["Mean", "Estimate"]
coef.LCI <- cr.data[["Output"]][["p"]]["LCI", "Estimate"]
coef.UCI <- cr.data[["Output"]][["p"]]["UCI", "Estimate"]
nq <- cr.data[["parameter.setting"]][["nq"]]
q_f <- gwqs_rank(model.data, mix_name, nq)
bQ <- q_f$Q
model.data$wqs <- as.numeric(as.matrix(bQ) %*% as.vector(weights))
model.data$deaths <- coef.mean*model.data$wqs*model.data$Pop
model.data$deaths.LCI <- coef.LCI*model.data$wqs*model.data$Pop
model.data$deaths.UCI <- coef.UCI*model.data$wqs*model.data$Pop


result.county <- aggregate(deaths ~ GEOID, data = model.data, mean, na.rm =TRUE) #per month
result.county$Pop <- aggregate(Pop ~ GEOID, data = model.data, mean, na.rm =TRUE)$Pop
result.county$death.rate <- result.county$deaths/aggregate(Pop ~ GEOID, data = model.data, mean, na.rm =TRUE)$Pop*1e6 #per month per 1e6
for (yeari in trend_years) {
  model.data.year <- model.data %>%
    filter(year == yeari)
  result.county[,paste0("death.rate.", yeari)] <- aggregate(deaths ~ GEOID, data = model.data.year, mean, na.rm =TRUE)$deaths/aggregate(Pop ~ GEOID, data = model.data.year, mean, na.rm =TRUE)$Pop*1e6
}
result.year <- aggregate(deaths ~ year, data = model.data, sum)
result.year$deaths.LCI <- aggregate(deaths.LCI ~ year, data = model.data, sum)$deaths.LCI
result.year$deaths.UCI <- aggregate(deaths.UCI ~ year, data = model.data, sum)$deaths.UCI

for (AP in mix_name){
  model.data[, c(paste0(AP, "_contri"))] <- weights[,AP]*bQ[, AP]/model.data$wqs
  model.data$temp.j <- weights[, AP]*bQ[, AP]*model.data$Pop
  model.data$temp.mix <- model.data$wqs*model.data$Pop
  result.county[, AP] <- aggregate(temp.j ~ GEOID, data = model.data, sum)$temp.j/aggregate(temp.mix ~ GEOID, data = model.data, sum)$temp.mix
  result.year[, AP] <- aggregate(temp.j ~ year, data = model.data, sum)$temp.j/aggregate(temp.mix ~ year, data = model.data, sum)$temp.mix
  for (yeari in trend_years) {
    model.data.year <- model.data %>%
      filter(year == yeari)
    result.county[, paste0(AP, ".", yeari)] <- aggregate(temp.j ~ GEOID, data = model.data.year, sum, na.rm =TRUE)$temp.j/aggregate(temp.mix ~ GEOID, data = model.data.year, sum, na.rm =TRUE)$temp.mix
  }
}
write.csv(result.year, file = paste0(save.path, "burden_year.csv"))

result.county$MAX.com <- apply(result.county[,mix_name], 1, function(x){
  names(x)[which.max(x)]
})
for (yeari in trend_years) {
  result.county[, paste0("MAX.com.", yeari)] <- apply(
    result.county[, paste0(mix_name, ".", yeari)], 
    1, 
    function(x) {
      if (is.null(names(x)) || all(is.na(x))) {
        return("No Data")
      }
      names(x)[which.max(x)]
    }
  )
  result.county[,paste0("MAX.com.", yeari)] <- sub(paste0("\\.", yeari), "", result.county[,paste0("MAX.com.", yeari)])
}


result.map <- merge(shp.data, result.county, c("GEOID"))
saveRDS(result.map, file = paste0(save.path, "burden_mapdata.rds"))

for (yeari in trend_years) {
  # Figure 3C and D
  name.year <- paste0("MAX.com.", yeari)
  max.map <- ggplot() +
    geom_sf(data = result.map, aes(fill = .data[[name.year]]), color = NA) +
    scale_fill_manual(values = c("BC" = "#67000D",
                                 "NH4" = "#EF3B2C",
                                 "NIT" = "#FEE0D2",
                                 "OM" = "#C6DBEF",
                                 "SO4" = "#6BAED6",
                                 "SOIL" = "#2171B5"),
                      name = "Components") +
    geom_sf(data = contiguous_states, fill = NA, color = "white", size = 1) +
    geom_sf(data = contiguous_nation, fill = NA, color = "black", size = 2) +
    theme_minimal() +
    theme(plot.title = element_text(size = 8, face = "bold"),
          legend.position = "none",
          legend.title = element_text(size = 6),
          legend.text = element_text(size = 6),
          legend.position.inside = c(0.1, 0.1),
          legend.justification = c(0, 0),
          legend.background = element_rect(fill = "white", color = NA),
          legend.key = element_rect(color = "black", fill = NA),
          legend.key.width = unit(0.2, "inches"),
          legend.key.height = unit(0.2, "inches"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank())
  ggsave(paste0(save.path, "max_com_", yeari, "_map.pdf"), plot = max.map, width = 3.5, height = 2.5)
  
  # Figure 3A and B
  print(max(result.map[[paste0("death.rate.", yeari)]] ,na.ram=TRUE))
  LIM <- 15
  result.map[[paste0("death.rate.", yeari)]] <- ifelse( result.map[[paste0("death.rate.", yeari)]] > LIM, LIM,  result.map[[paste0("death.rate.", yeari)]])
  
  death.rate.map <- ggplot() +
    geom_sf(data = result.map, aes(fill = .data[[paste0("death.rate.", yeari)]]), color = NA) +
    scale_fill_gradientn(colors = c("white", "#EF3B2C", "#67000D"),
                         values = c(0, 0.5, 1),
                         limits = c(0, LIM),
                         breaks = seq(0, LIM, by = 5),
                         name = "C") +
    geom_sf(data = contiguous_states, fill = NA, color = "white", size = 1) +
    geom_sf(data = contiguous_nation, fill = NA, color = "black", size = 2) +
    theme_minimal() +
    theme(plot.title = element_text(size = 8, face = "bold"),
          legend.title = element_text(size = 6),
          legend.text = element_text(size = 6),
          legend.position = "bottom",
          legend.direction = "horizontal",
          legend.background = element_rect(fill = "white", color = NA),
          legend.key = element_rect(color = "black", fill = NA),
          legend.key.width = unit(0.25, "inches"),
          legend.key.height = unit(0.1, "inches"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank())
  ggsave(paste0(save.path, "death_rate_", yeari, "_map.pdf"), plot = death.rate.map, width = 3.5, height = 2.5)
}
