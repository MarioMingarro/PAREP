source("Dependencies/Fun.R")

gc(reset = T)


tic()


dir_present_climate_data <- "C:/A_TRABAJO/A_GABRIEL/REPRESENTATIVIDAD/TEST_SIMPLE/P/"
dir_future_climate_data <- "C:/A_TRABAJO/A_GABRIEL/REPRESENTATIVIDAD/TEST_SIMPLE/F/"
dir_result <- "C:/A_TRABAJO/A_GABRIEL/REPRESENTATIVIDAD/TEST_SIMPLE/"
#OLD/TEST_PNAC/"
#Peninsula_Iberica_89.shp"
#national_parks.shp"
polygon <- read_sf("C:/A_TRABAJO/A_GABRIEL/REPRESENTATIVIDAD/TEST_SIMPLE/AES.shp")


# Create name object

year <- "999"
model <- "TEST"

names <- polygon$Id


# Crear las subcarpetas 'presente' y 'futuro' dentro de 'dir_result'
dir_present <- paste0(dir_result, "Present/")
dir_fut <- paste0(dir_result, "Future/")
dir_future <- paste0(dir_fut, "TEST/")
dir_charts <- paste0(dir_result, "Charts/")

# Crear las carpetas si no existen
if (!dir.exists(dir_result)) {
  dir.create(dir_result)
}
if (!dir.exists(dir_present)) {
  dir.create(dir_present)
}

if (!dir.exists(dir_fut)) {
  dir.create(dir_fut)
}

if (!dir.exists(dir_futu)) {
  dir.create(dir_futu)
}

if (!dir.exists(dir_future)) {
  dir.create(dir_future)
}
if (!dir.exists(dir_charts)) {
  dir.create(dir_charts)
}

# CLIMATE ----
## Load data ----
present_climatic_variables <- terra::rast(list.files(dir_present_climate_data, "\\.tif$", full.names = T))

future_climatic_variables <- terra::rast(list.files(dir_future_climate_data, "\\.tif$", full.names = T))


names(future_climatic_variables) <- names(present_climatic_variables)

# Reference system ----
terra::crs(present_climatic_variables)
reference_system <-"EPSG:4326" 

polygon <- st_transform(polygon, crs(reference_system))

# Crop raster to study area
# Raster to df
data_present_climatic_variables <- terra::as.data.frame(present_climatic_variables, xy = TRUE)
data_future_climatic_variables <- terra::as.data.frame(future_climatic_variables, xy = TRUE)
# Delete NA
data_present_climatic_variables<-na.omit(data_present_climatic_variables)
data_future_climatic_variables <-na.omit(data_future_climatic_variables)

# Pre VIF
jpeg(paste0(dir_result, "All variables correlation.jpeg"), quality = 75, width = 1200, height = 600)
par(mfrow = c(1, 2))
corrplot(cor(data_present_climatic_variables[3:length(data_present_climatic_variables)]),
         method = "number", type = "upper", title = "Present variables correlation")
corrplot(cor(data_future_climatic_variables[3:length(data_future_climatic_variables)]),
         method = "number", type = "upper", title = "Future variables correlation")
dev.off()

# VIF
vif_r <- vif_filter(present_climatic_variables, th = 10)
drop <- vif_r$excluded

data_present_climatic_variables <- data_present_climatic_variables[!names(data_present_climatic_variables) %in% drop]
data_future_climatic_variables <- data_future_climatic_variables[!names(data_future_climatic_variables) %in% drop]

present_climatic_variables <- terra::subset(present_climatic_variables, !names(present_climatic_variables) %in% drop)
future_climatic_variables <- terra::subset(future_climatic_variables, !names(future_climatic_variables) %in% drop)

# Post VIF
jpeg(paste0(dir_result, "VIF selected variables correlation.jpeg"), quality = 75, width = 1200, height = 600)
par(mfrow = c(1, 2))
corrplot(cor(data_present_climatic_variables[3:length(data_present_climatic_variables)]),
         method = "number", type = "upper", title = "Present variables correlation")
corrplot(cor(data_future_climatic_variables[3:length(data_future_climatic_variables)]),
         method = "number", type = "upper", title = "Future variables correlation")
dev.off()

# Add field period 
data_present_climatic_variables <- dplyr::mutate(data_present_climatic_variables, Period = c("Present"),  .after = "y")
data_future_climatic_variables  <- dplyr::mutate(data_future_climatic_variables, Period = c("Future"),  .after = "y")

# Join two dataset
colnames(data_future_climatic_variables) <- colnames(data_present_climatic_variables)
#data <- rbind(data_present_climatic_variables, data_future_climatic_variables)




tic()
for(j in 1:length(names)){
  pa_mh_present_future(j, th= .95)
}
toc()