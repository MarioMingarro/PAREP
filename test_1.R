# library(RStoolbox)
# library(stringr)
source("Dependencies/Fun.R")

gc(reset = T)
dir_present_climate_data <- "C:/A_TRABAJO/A_GABRIEL/REPRESENTATIVIDAD/CLIMA/PRESENT/"
dir_future_climate_data <- "C:/A_TRABAJO/A_GABRIEL/REPRESENTATIVIDAD/CLIMA/FUTURO/IPSL/"
dir_result <- "C:/A_TRABAJO/A_GABRIEL/REPRESENTATIVIDAD/FAST_TEST/"

study_area <- read_sf("C:/A_TRABAJO/A_GABRIEL/REPRESENTATIVIDAD/FAST_TEST/MURCIA.shp")
polygon <- read_sf("C:/A_TRABAJO/A_GABRIEL/REPRESENTATIVIDAD/FAST_TEST/PA.shp")
plot(study_area$geometry)
plot(polygon$geometry, add = T)

# Create name object
names <- polygon$NatName
year <- "2070"
model <- "GFDL"


# Crear las subcarpetas 'presente' y 'futuro' dentro de 'dir_result'
dir_present <- paste0(dir_result, "Present/")
dir_fut <- paste0(dir_result, "Future/")
dir_futu <- paste0(dir_fut, year,"/")
dir_future <- paste0(dir_futu, model,"/")

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

# CLIMATE ----
## Load data ----
present_climatic_variables <- terra::rast(list.files(dir_present_climate_data, "\\.tif$", full.names = T))
present_climatic_variables <- present_climatic_variables[[1:17]]

future_climatic_variables <- terra::rast(list.files(dir_future_climate_data, "\\.tif$", full.names = T))
future_climatic_variables <- future_climatic_variables[[1:17]]

names(present_climatic_variables) <- c("CHELSA_bio1","CHELSA_bio10","CHELSA_bio11","CHELSA_bio12","CHELSA_bio13","CHELSA_bio14",
                                       "CHELSA_bio15","CHELSA_bio16","CHELSA_bio17","CHELSA_bio18","CHELSA_bio19","CHELSA_bio2",
                                       "CHELSA_bio3","CHELSA_bio4","CHELSA_bio5","CHELSA_bio6","CHELSA_bio7")

names(future_climatic_variables) <- names(present_climatic_variables)

# Reference system ----
terra::crs(present_climatic_variables)
reference_system <-"EPSG:4326" 

study_area <- st_transform(study_area, crs(reference_system))
polygon <- st_transform(polygon, crs(reference_system))


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


# Create name object
names <- polygon$NatName


tic()
for(j in 1:length(names)){
  pa_mh_present_future(j)
}
toc()



tic()
for(j in 1:length(names)){
  pa_mh_present(j)
}
for(j in 1:length(names)){
  pa_mh_future(j)
}
toc()



pa_mh_present_future

puntos_todos_p <- terra::as.points(mh_raster_p)
puntos_todos_p <- sf::st_as_sf(puntos_todos_p)
colnames(puntos_todos_p) <- c("mh", "geometry")
puntos_dentro <- sf::st_intersection(puntos_todos_p, pol)
puntos_todos_p$inv_mh <- 1 / puntos_todos_p$mh

th_mh_p <- quantile(na.omit(puntos_dentro$mh), probs = th)

coords <- st_coordinates(puntos_todos_p)
puntos_todos_p <- cbind(puntos_todos_p, coords)

# Distancia euclidea ----

# Calcular la distancia mínima desde cada punto del raster a los puntos dentro del polígono
dist <- sf::st_distance(puntos_todos_p, puntos_dentro)

# Para obtener la distancia mínima por cada punto, tomamos el valor mínimo de cada fila
dist <- apply(dist, 1, min)

# Agregar las distancias mínimas como un atributo a los puntos del raster
puntos_todos$dist <- dist

# Convertir las distancias a kilómetros y redondear
puntos_todos$dist <- round(puntos_todos$dist / 1000, 0)
puntos_todos$dist[puntos_todos$dist == 0] <- NA


# Autocorrelacion ----

nb <- spdep::dnearneigh(coords, 0, 0.015, use_s2 = TRUE, dwithin = TRUE)
lw <- spdep::nb2listw(nb, style = "W", zero.policy = TRUE)

LM <- spdep::localmoran(puntos_todos_p$inv_mh, lw)
LM_df <- as.data.frame(LM)

puntos_todos_p$SA <- LM_df$Ii
puntos_todos_p$SA_sig <- LM_df$`Pr(z != E(Ii))`

puntos_todos_p$th <- case_when(
  puntos_todos_p$mh > th_mh_p ~ 0,
  puntos_todos_p$mh <= th_mh_p  ~ 1
)

puntos_todos_p$th <- as.numeric(puntos_todos_p$th)

res <- res(mh_raster_p)
bbox <- ext(mh_raster_p)

nrows <- round((bbox[4] - bbox[3]) / res[2])
ncols <- round((bbox[2] - bbox[1]) / res[1])
raster_template <- rast(ext = bbox, nrows = nrows, ncols = ncols)
puntos_vect_p <- vect(puntos_todos_p)

raster <- terra::rasterize(puntos_vect_p, raster_template, field = "th")
crs(raster) <- crs(mh_raster_p)
plot(raster)
writeRaster(raster, paste0(dir_result, "Pre_th_", names[j], ".tif"), overwrite = TRUE)

raster <- terra::rasterize(puntos_vect_p, raster_template, field = "SA")
crs(raster) <- crs(mh_raster_p)
writeRaster(raster, paste0(dir_result, "Pre_SA_", names[j], ".tif"), overwrite = TRUE)



rm(list = setdiff(ls(), c("data_present_climatic_variables", "data_future_climatic_variables", "dir_future_climate_data",
                          "dir_present_climate_data", "dir_result", "future_climatic_variables", "polygon", 
                          "present_climatic_variables", "reference_system", "study_area", "names", "model", "year"
)))  


# Cargar las bibliotecas necesarias
library(sf)       # Para trabajar con datos espaciales
library(spdep)    # Para análisis de autocorrelación espacial (Moran's I)
library(dplyr)    # Para manipulación de datos
library(geosphere) # Para cálculo de distancias geográficas

# Paso 1: Crear el conjunto de puntos dentro del polígono y los puntos en toda el área
# Supongamos que ya tienes los objetos `puntos_todos_p` (todos los puntos) y `pol` (el polígono)

# Intersección para obtener los puntos dentro del polígono
puntos_dentro <- st_intersection(puntos_todos_p, pol)

# Distancia euclidea ----

# Calcular la distancia mínima desde cada punto del raster a los puntos dentro del polígono
dist <- sf::st_distance(puntos_todos, puntos_dentro)

# Para obtener la distancia mínima por cada punto, tomamos el valor mínimo de cada fila
dist <- apply(dist, 1, min)

# Agregar las distancias mínimas como un atributo a los puntos del raster
puntos_todos$dist <- dist

# Convertir las distancias a kilómetros y redondear
puntos_todos$dist <- round(puntos_todos$dist / 1000, 0)
puntos_todos$dist[puntos_todos$dist == 0] <- NA

# Paso 2: Calcular la distancia de cada punto en `puntos_todos_p` al borde del polígono
# Convertimos el polígono a una geometría de líneas para calcular la distancia al borde
borde_poligono <- st_cast(pol, "MULTILINESTRING")

# Calcula la distancia de cada punto al borde del polígono para todos los puntos
puntos_todos_p <- puntos_todos_p %>%
  mutate(distancia_al_borde = st_distance(geometry, borde_poligono))



# Paso 3: Crear una matriz de pesos espaciales en función de la distancia entre todos los puntos
# Primero transformamos a coordenadas proyectadas si los datos están en lat/lon para mejor manejo de distancias
puntos_todos_proj <- st_transform(puntos_todos_p, crs = 3857)  # CRS métrico (Web Mercator)

# Extraer coordenadas para crear la matriz de pesos de distancia
coords <- st_coordinates(puntos_todos_proj)
# Calcular la matriz de distancia

library(spdep)

# Paso 1: Definir el umbral de distancia en metros (por ejemplo, 500 metros)
umbral_distancia <- 500000

# Paso 2: Crear una lista de vecinos basados en distancias menores al umbral
lista_vecinos <- dnearneigh(coords, 0, umbral_distancia)

# Paso 3: Crear la lista de pesos espaciales con el umbral de distancia
matriz_pesos <- nb2listw(lista_vecinos, style = "W", zero.policy = TRUE)

# Calcular el índice de Moran Local
moran_local <- localmoran(puntos_todos_p$variable_de_interes, matriz_pesos)


# Paso 5: Guardar los valores de autocorrelación y distancia en una nueva tabla de resultados
# Extraemos solo los puntos dentro del polígono para los resultados finales
resultados <- puntos_todos_p %>%
  mutate(moran_I = moran_local[,1],     # Ii de Moran (autocorrelación local)
         distancia_al_borde = as.numeric(distancia_al_borde))

# Paso 6: Visualizar los resultados en un gráfico de autocorrelación vs. distancia al borde
plot(resultados$distancia_al_borde, resultados$moran_I,
     xlab = "Distancia al borde del polígono",
     ylab = "Autocorrelación Local de Moran (Ii)",
     main = "Variación de la autocorrelación local de Moran con la distancia al borde")
abline(lm(moran_I ~ distancia_al_borde, data = resultados), col = "blue")

