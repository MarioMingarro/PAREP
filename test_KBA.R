
source("Dependencies/Fun.R")

gc(reset = T)


tic()

dir_present_climate_data <- "C:/A_TRABAJO/A_GABRIEL/REPRESENTATIVIDAD/CLIMA/PRESENT/"
dir_future_climate_data <- "C:/A_TRABAJO/A_GABRIEL/REPRESENTATIVIDAD/CLIMA/FUTURE/GFDL/"
dir_result <- "C:/A_TRABAJO/A_GABRIEL/REPRESENTATIVIDAD/FAST_TEST/"

study_area <- read_sf("C:/A_TRABAJO/A_GABRIEL/REPRESENTATIVIDAD/FAST_TEST/MURCIA.shp")
polygon <- read_sf("C:/A_TRABAJO/A_GABRIEL/REPRESENTATIVIDAD//FAST_TEST/PAS_murcia.shp")



# Create name object

year <- "2070"
model <- "GFDL"

names <- polygon$NatName


# Crear las subcarpetas 'presente' y 'futuro' dentro de 'dir_result'
dir_present <- paste0(dir_result, "Present/")
dir_fut <- paste0(dir_result, "Future/")
dir_futu <- paste0(dir_fut, year,"/")
dir_future <- paste0(dir_futu, model,"/")
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
# Crear un vector con los nombres de las variables a excluir
exclude_vars <- c("bio8", "bio9", "bio18", "bio19")

# Crear un patrón de expresión regular para excluir estas variables
exclude_pattern <- paste0("bio(", paste(gsub("bio", "", exclude_vars), collapse = "|"), ")")

# Seleccionar las variables deseadas
present_climatic_variables <- subset(present_climatic_variables, grep(exclude_pattern, names(present_climatic_variables), invert = TRUE, value = TRUE))

future_climatic_variables <- terra::rast(list.files(dir_future_climate_data, "\\.tif$", full.names = T))
future_climatic_variables <- subset(future_climatic_variables, grep(exclude_pattern, names(future_climatic_variables), invert = TRUE, value = TRUE))


names(present_climatic_variables) <- c("CHELSA_bio1","CHELSA_bio10","CHELSA_bio11","CHELSA_bio12","CHELSA_bio13","CHELSA_bio14",
                                       "CHELSA_bio15","CHELSA_bio16","CHELSA_bio17","CHELSA_bio2",
                                       "CHELSA_bio3","CHELSA_bio4","CHELSA_bio5","CHELSA_bio6","CHELSA_bio7")

names(future_climatic_variables) <- names(present_climatic_variables)

# Reference system ----
terra::crs(present_climatic_variables)
reference_system <-"EPSG:4326" 

study_area <- st_transform(study_area, crs(reference_system))
polygon <- st_transform(polygon, crs(reference_system))


polygon <- st_intersection(st_crop(polygon, st_bbox(study_area)), study_area)

plot(study_area$geometry)
plot(polygon$geometry, add = T, col = "red")

# Crop raster to study area
present_climatic_variables <-  terra::mask (crop(present_climatic_variables, study_area), study_area)
future_climatic_variables  <-  terra::mask(crop(future_climatic_variables,  study_area), study_area)
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
data <- rbind(data_present_climatic_variables, data_future_climatic_variables)

tic()
for(j in 1:length(names)){
  pa_mh_present_future2(j, th= .95)
}
toc()
