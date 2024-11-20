library(caret)

library(RStoolbox)
library(stringr)
library(corrplot)

library(terra)
library(sf)
library(tidyverse)
library(caret)
library(tictoc)
library(spdep)


# Load packages ----
packages.to.use <- c("corrplot", "terra", "tictoc", "tidyverse","sf")

packages.to.use <- unique(packages.to.use)

for(package in packages.to.use) {
  print(package)
  if( ! package %in% rownames(installed.packages()) ) { install.packages(package ) }
  if( ! package %in% rownames(installed.packages()) ) { stop("Error on package instalation") }
  suppressWarnings( library(package, character.only = TRUE) )
}


# VIF----

vif_filter <- function(x, th = 10) {
  
  # Función para calcular el VIF de todas las variables
  calc_vif <- function(df) {
    vif_values <- sapply(1:ncol(df), function(i) {
      formula <- as.formula(paste(names(df)[i], "~ ."))  # Formula de lm con todas las demás variables
      model <- lm(formula, data = df)
      return(1 / (1 - summary(model)$r.squared))
    })
    names(vif_values) <- colnames(df)
    return(vif_values)
  }
  
  # Convertir el objeto Raster a un data frame
  x <- as.data.frame(x, na.rm = TRUE)
  
  
  # Eliminar variables multicolineales iterativamente
  exc <- character(0)  # Lista de variables excluidas
  while (TRUE) {
    v <- calc_vif(x)  # Calcular el VIF
    if (max(v) < th) break  # Terminar si no hay VIF por encima del umbral
    ex <- names(v)[which.max(v)]  # Variable con mayor VIF
    exc <- c(exc, ex)  # Agregar a la lista de excluidos
    x <- x[, !(colnames(x) %in% ex)]  # Remover variable
  }
  
  # Crear lista de resultados
  result <- list(
    variables = colnames(x),
    excluded = exc,
    corMatrix = cor(x, method = "pearson"),
    results = data.frame(Variables = names(v), VIF = v)
  )
  
  
  print(result)
  return(result)
}

# Mahalanobis distance ----
## PRESENT ----
# Función para procesar los datos de la serie presente
pa_mh_present <- function(j, th = .95) {
  data_p <- data_present_climatic_variables
  mh_p <- data.frame(matrix(1,    
                            nrow = nrow(data_p),
                            ncol = length(names)))
  
  names(mh_p) <- names
  
  for (i in 1:nrow(polygon)){
    pol <- polygon[i,]
    raster_polygon <- terra::mask(terra::crop(present_climatic_variables, pol), pol)
    data_polygon <- terra::as.data.frame(raster_polygon, xy = TRUE)
    data_polygon <- na.omit(data_polygon)
    
    mh <- mahalanobis(data_p[,4:length(data_p)], 
                      colMeans(data_polygon[,3:length(data_polygon)]), 
                      cov(data_p[,4:length(data_p)]), 
                      inverted = F)
    
    mh_p[,i] <- mh
  }
  
  mh_p <- cbind(data_p[,c(1:3)], mh_p)
  mh_raster_p <- terra::rast(mh_p[, c(1:2, j+3)], crs = reference_system)
  names(mh_raster_p) <- colnames(mh_p[j+3])
  plot(mh_raster_p)
  writeRaster(mh_raster_p, paste0(dir_present, "MH_PRESENT_", names[j], ".tif"), overwrite = TRUE)
  assign(paste0("MH_PRESENT_", names[j]), mh_raster_p, envir = .GlobalEnv)
  
  
  puntos_todos_p <- terra::as.points(mh_raster_p)
  puntos_todos_p <- sf::st_as_sf(puntos_todos_p)
  colnames(puntos_todos_p) <- c("mh", "geometry")
  puntos_dentro_p <- sf::st_intersection(puntos_todos_p, pol)
  
  th_mh_p <- quantile(na.omit(puntos_dentro_p$mh), probs = th)
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
  writeRaster(raster, paste0(dir_present, "TH_MH_PRESENT_", names[j], ".tif"), overwrite = TRUE)
  assign(paste0("TH_MH_PRESENT_", names[j]), mh_raster_p, envir = .GlobalEnv)

}



## FUTURE----
pa_mh_future <- function(j, th = .95) {
  data_f <- data_future_climatic_variables
  mh_f <- data.frame(matrix(1,    
                            nrow = nrow(data_f),
                            ncol = length(names)))
  names(mh_f) <- names
  
  for (i in 1:nrow(polygon)){
    pol <- polygon[i,]
    raster_polygon <- terra::mask(terra::crop(future_climatic_variables, pol), pol)
    data_polygon <- terra::as.data.frame(raster_polygon, xy = TRUE)
    data_polygon <- na.omit(data_polygon)
    
    mh <- mahalanobis(data_f[,4:length(data_f)], 
                      colMeans(data_polygon[,3:length(data_polygon)]), 
                      cov(data_f[,4:length(data_f)]), 
                      inverted = F)
    
    mh_f[,i] <- mh
  }
  
  mh_f <- cbind(data_f[,c(1:3)], mh_f)
  mh_raster_f <- terra::rast(mh_f[, c(1:2, j+3)], crs = reference_system)
  names(mh_raster_f) <- colnames(mh_f[j+3])
  plot(mh_raster_f)
  writeRaster(mh_raster_f, paste0(dir_future, "MH_", model, "_", year, "_", names[j], ".tif"), overwrite = TRUE)
  assign(paste0("MH_", model, "_", year, "_", names[j]), mh_raster_f, envir = .GlobalEnv)
  
  puntos_todos_f <- terra::as.points(mh_raster_f)
  puntos_todos_f <- sf::st_as_sf(puntos_todos_f)
  colnames(puntos_todos_f) <- c("mh", "geometry")
  puntos_dentro_f <- sf::st_intersection(puntos_todos_f, pol)
  
  th_mh_f <- quantile(na.omit(puntos_dentro_f$mh), probs = th)
  puntos_todos_f$th <- case_when(
    puntos_todos_f$mh > th_mh_f ~ 0,
    puntos_todos_f$mh <= th_mh_f  ~ 1
  )
  puntos_todos_f$th <- as.numeric(puntos_todos_f$th)
  
  res <- res(mh_raster_f)
  bbox <- ext(mh_raster_f)
  
  nrows <- round((bbox[4] - bbox[3]) / res[2])
  ncols <- round((bbox[2] - bbox[1]) / res[1])
  raster_template <- rast(ext = bbox, nrows = nrows, ncols = ncols)
  puntos_vect_f <- vect(puntos_todos_f)
  
  raster <- terra::rasterize(puntos_vect_f, raster_template, field = "th")
  crs(raster) <- crs(mh_raster_f)
  plot(raster)
  writeRaster(raster,
              paste0(dir_future, "TH_MH_", model, "_", year, "_", names[j], ".tif"),
              overwrite = TRUE)
  
  assign(paste0("TH_MH_", model, "_", year, "_", names[j]), mh_raster_p, envir = .GlobalEnv)
}

## PRESENT-FUTURE----
# Función para procesar los datos de la serie presente
pa_mh_present_future <- function(j, th = .95) {
  
  # Presente
  data_p <- data_present_climatic_variables
  
  mh_p <- data.frame(matrix(1,    
                            nrow = nrow(data_p),
                            ncol = length(names)))
  
  names(mh_p) <- names
  
  for (i in 1:nrow(polygon)){
    pol <- polygon[i,]
    raster_polygon <- terra::mask(terra::crop(present_climatic_variables, pol), pol)
    data_polygon <- terra::as.data.frame(raster_polygon, xy = TRUE)
    data_polygon <- na.omit(data_polygon)
    
    mh <- mahalanobis(data_p[,4:length(data_p)], 
                      colMeans(data_polygon[,3:length(data_polygon)]), 
                      cov(data_p[,4:length(data_p)]), 
                      inverted = F)
    
    mh_p[,i] <- mh
  }
  
  mh_p <- cbind(data_p[,c(1:3)], mh_p)
  mh_raster_p <- terra::rast(mh_p[, c(1:2, j+3)], crs = reference_system)
  names(mh_raster_p) <- colnames(mh_p[j+3])
  plot(mh_raster_p)
  writeRaster(mh_raster_p, paste0(dir_present, "MH_PRESENT_", names[j], ".tif"), overwrite = TRUE)
  assign(paste0("MH_PRESENT_", names[j]), mh_raster_p, envir = .GlobalEnv)
  
  puntos_todos_p <- terra::as.points(mh_raster_p)
  puntos_todos_p <- sf::st_as_sf(puntos_todos_p)
  colnames(puntos_todos_p) <- c("mh", "geometry")
  puntos_dentro_p <- sf::st_intersection(puntos_todos_p, pol)
  
  th_mh_p <- quantile(na.omit(puntos_dentro_p$mh), probs = th)
  
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
  writeRaster(raster, paste0(dir_present, "TH_MH_PRESENT_", names[j], ".tif"), overwrite = TRUE)
  assign(paste0("TH_MH_PRESENT_", names[j]), mh_raster_p, envir = .GlobalEnv)
  
  # Futuro
  data_p_f <- rbind(data_present_climatic_variables, data_future_climatic_variables)
  
  mh_p_f <- data.frame(matrix(1,    
                            nrow = nrow(data_p_f),
                            ncol = length(names)))
  
  names(mh_p_f) <- names
  
  for (i in 1:nrow(polygon)) {
    pol <- polygon[i, ]
    raster_polygon <- terra::mask(terra::crop(present_climatic_variables, pol), pol)
    data_polygon <- terra::as.data.frame(raster_polygon, xy = TRUE)
    data_polygon <- na.omit(data_polygon)
    
    mh <- mahalanobis(data_p_f[, 4:length(data_p_f)], colMeans(data_polygon[, 3:length(data_polygon)]), cov(data_p_f[, 4:length(data_p_f)]), inverted = F)
    
    mh_p_f[, i] <- mh
  }
  
  mh_p_f <- cbind(data_p_f[, c(1:3)], mh_p_f)
  
  mh_raster_p_f <- dplyr::filter(mh_p_f, Period == "Future")
  mh_raster_p_f <- terra::rast(mh_raster_p_f[, c(1:2, j+3)], crs = reference_system)
  names(mh_raster_p_f) <- colnames(mh_p_f[j+3])
  plot(mh_raster_p_f)
  writeRaster(mh_raster_p_f,
              paste0(dir_future, "MH_", model, "_", year, "_", names[j], ".tif"),
              overwrite = TRUE)
  
  puntos_todos_f <- terra::as.points(mh_raster_p_f)
  puntos_todos_f <- sf::st_as_sf(puntos_todos_f)
  colnames(puntos_todos_f) <- c("mh", "geometry")
  puntos_dentro_f <- sf::st_intersection(puntos_todos_f, pol)
  
  puntos_todos_f$th <- case_when(
    puntos_todos_f$mh > th_mh_p ~ 0,
    puntos_todos_f$mh <= th_mh_p  ~ 1
  )
  puntos_todos_f$th <- as.numeric(puntos_todos_f$th)
  
  res <- res(mh_raster_p_f)
  bbox <- ext(mh_raster_p_f)
  
  nrows <- round((bbox[4] - bbox[3]) / res[2])
  ncols <- round((bbox[2] - bbox[1]) / res[1])
  raster_template <- rast(ext = bbox, nrows = nrows, ncols = ncols)
  puntos_vect_f <- vect(puntos_todos_f)
  
  raster <- terra::rasterize(puntos_vect_f, raster_template, field = "th")
  crs(raster) <- crs(mh_raster_p_f)
  plot(raster)
  writeRaster(raster,
              paste0(dir_future, "TH_MH_", model, "_", year, "_", names[j], ".tif"),
              overwrite = TRUE)
  
  assign(paste0("TH_MH_", model, "_", year, "_", names[j]), mh_raster_p, envir = .GlobalEnv)
  
}



## PRESENT-FUTURE 2----
# Función para procesar los datos de la serie presente
pa_mh_present_future <- function(j, th = .95) {
  
  # Presente
  data_p <- data_present_climatic_variables
  
  mh_p <- data.frame(matrix(1,    
                            nrow = nrow(data_p),
                            ncol = length(names)))
  
  names(mh_p) <- names
  
  for (i in 1:nrow(polygon)){
    pol <- polygon[i,]
    raster_polygon <- terra::mask(terra::crop(present_climatic_variables, pol), pol)
    data_polygon <- terra::as.data.frame(raster_polygon, xy = TRUE)
    data_polygon <- na.omit(data_polygon)
    
    mh <- mahalanobis(data_p[,4:length(data_p)], 
                      colMeans(data_polygon[,3:length(data_polygon)]), 
                      cov(data_p[,4:length(data_p)]), 
                      inverted = F)
    
    mh_p[,i] <- mh
  }
  
  mh_p <- cbind(data_p[,c(1:3)], mh_p)
  mh_raster_p <- terra::rast(mh_p[, c(1:2, j+3)], crs = reference_system)
  names(mh_raster_p) <- colnames(mh_p[j+3])
  plot(mh_raster_p)
  writeRaster(mh_raster_p, paste0(dir_present, "MH_PRESENT_", names[j], ".tif"), overwrite = TRUE)
  assign(paste0("MH_PRESENT_", names[j]), mh_raster_p, envir = .GlobalEnv)
  
  puntos_todos_p <- terra::as.points(mh_raster_p)
  puntos_todos_p <- sf::st_as_sf(puntos_todos_p)
  colnames(puntos_todos_p) <- c("mh", "geometry")
  puntos_dentro_p <- sf::st_intersection(puntos_todos_p, pol)
  
  th_mh_p <- quantile(na.omit(puntos_dentro_p$mh), probs = th)
  
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
  writeRaster(raster, paste0(dir_present, "TH_MH_PRESENT_", names[j], ".tif"), overwrite = TRUE)
  assign(paste0("TH_MH_PRESENT_", names[j]), mh_raster_p, envir = .GlobalEnv)
  
  # Futuro
  data_p_f <- rbind(data_present_climatic_variables, data_future_climatic_variables)
  
  mh_p_f <- data.frame(matrix(1,    
                              nrow = nrow(data_p_f),
                              ncol = length(names)))
  
  names(mh_p_f) <- names
  
  for (i in 1:nrow(polygon)) {
    pol <- polygon[i, ]
    raster_polygon <- terra::mask(terra::crop(present_climatic_variables, pol), pol)
    data_polygon <- terra::as.data.frame(raster_polygon, xy = TRUE)
    data_polygon <- na.omit(data_polygon)
    
    mh <- mahalanobis(data_p_f[, 4:length(data_p_f)], colMeans(data_polygon[, 3:length(data_polygon)]), cov(data_p_f[, 4:length(data_p_f)]), inverted = F)
    
    mh_p_f[, i] <- mh
  }
  
  mh_p_f <- cbind(data_p_f[, c(1:3)], mh_p_f)
  
  mh_raster_p_f <- dplyr::filter(mh_p_f, Period == "Future")
  mh_raster_p_f <- terra::rast(mh_raster_p_f[, c(1:2, j+3)], crs = reference_system)
  names(mh_raster_p_f) <- colnames(mh_p_f[j+3])
  plot(mh_raster_p_f)
  writeRaster(mh_raster_p_f,
              paste0(dir_future, "MH_", model, "_", year, "_", names[j], ".tif"),
              overwrite = TRUE)
  
  puntos_todos_f <- terra::as.points(mh_raster_p_f)
  puntos_todos_f <- sf::st_as_sf(puntos_todos_f)
  colnames(puntos_todos_f) <- c("mh", "geometry")
  puntos_dentro_f <- sf::st_intersection(puntos_todos_f, pol)
  
  puntos_todos_f$th <- case_when(
    puntos_todos_f$mh > th_mh_p ~ 0,
    puntos_todos_f$mh <= th_mh_p  ~ 1
  )
  puntos_todos_f$th <- as.numeric(puntos_todos_f$th)
  
  res <- res(mh_raster_p_f)
  bbox <- ext(mh_raster_p_f)
  
  nrows <- round((bbox[4] - bbox[3]) / res[2])
  ncols <- round((bbox[2] - bbox[1]) / res[1])
  raster_template <- rast(ext = bbox, nrows = nrows, ncols = ncols)
  puntos_vect_f <- vect(puntos_todos_f)
  
  raster <- terra::rasterize(puntos_vect_f, raster_template, field = "th")
  crs(raster) <- crs(mh_raster_p_f)
  plot(raster)
  writeRaster(raster,
              paste0(dir_future, "TH_MH_", model, "_", year, "_", names[j], ".tif"),
              overwrite = TRUE)
  
  assign(paste0("TH_MH_", model, "_", year, "_", names[j]), mh_raster_p, envir = .GlobalEnv)
  
  
  dist <- sf::st_distance(puntos_todos_p, puntos_dentro_p)
  
  # Para obtener la distancia mínima por cada punto, tomamos el valor mínimo de cada fila
  dist <- apply(dist, 1, min)
  
  # Agregar las distancias mínimas como un atributo a los puntos del raster
  puntos_todos_p$dist <- round(dist/ 1000, 0)
  puntos_todos_f$dist <- round(dist/ 1000, 0)
  
  
  puntos_todos_p <- cbind(puntos_todos_p, st_coordinates(puntos_todos_p))
  puntos_todos_f <- cbind(puntos_todos_f, st_coordinates(puntos_todos_f))
  #write.csv2(puntos_todos_f, "C:/A_TRABAJO/A_GABRIEL/REPRESENTATIVIDAD/FAST_TEST/GEODA/puntos_todos_p.csv" )
  rango_breaks <- seq(0, max(puntos_todos_p$dist, na.rm = TRUE), by = 10)
  
  puntos_todos_p <- puntos_todos_p %>%
    mutate(rango_distancia = cut(dist, breaks = rango_breaks, include.lowest = TRUE))
  puntos_todos_f <- puntos_todos_f %>%
    mutate(rango_distancia = cut(dist, breaks = rango_breaks, include.lowest = TRUE))
  
  resultados_p <- puntos_todos_p %>%
    filter(!is.na(rango_distancia)) %>%  # Filtrar filas con rango de distancia válido
    group_by(rango_distancia) %>%  # Agrupar por rango de distancia
    summarise(
      n_total = n(),  # Número total de observaciones en el rango
      n_ceros = sum(th == 0, na.rm = TRUE),  # Número de ceros en 'th'
      n_unos = sum(th == 1, na.rm = TRUE)  # Número de unos en 'th'
    ) %>%
    mutate(
      porcentaje_ceros = (n_ceros / n_total) * 100,  # Porcentaje de ceros
      porcentaje_unos = (n_unos / n_total) * 100  # Porcentaje de unos
    )
  
  resultados_f <- puntos_todos_f %>%
    filter(!is.na(rango_distancia)) %>%  # Filtrar filas con rango de distancia válido
    group_by(rango_distancia) %>%  # Agrupar por rango de distancia
    summarise(
      n_total = n(),  # Número total de observaciones en el rango
      n_ceros = sum(th == 0, na.rm = TRUE),  # Número de ceros en 'th'
      n_unos = sum(th == 1, na.rm = TRUE)  # Número de unos en 'th'
    ) %>%
    mutate(
      porcentaje_ceros = (n_ceros / n_total) * 100,  # Porcentaje de ceros
      porcentaje_unos = (n_unos / n_total) * 100  # Porcentaje de unos
    )
  
  
  ggplot() +
    geom_point(data =resultados_p, aes(x = rango_distancia, y = porcentaje_unos, color = "Present"), size = 3) +  # Puntos para los 'Unos'
    geom_line(data =resultados_p, aes(x = rango_distancia, y = porcentaje_unos, color = "Present"), size = 1, group = 1) +  # Línea para conectar los puntos
    geom_point(data =resultados_f, aes(x = rango_distancia, y = porcentaje_unos, color = "Future"), size = 3) +  # Puntos para los 'Unos'
    geom_line(data =resultados_f, aes(x = rango_distancia, y = porcentaje_unos, color = "Future"), size = 1, group = 1) +
    labs(
      title = "Distribución de ceros y unos según el rango de distancia",
      x = "Distancia (km)",
      y = "Porcentaje receptoras",
      fill = "Valores de mh",
      color = "Valores de mh"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

kk <- select(puntos_todos_p, c(X, Y, th, mh, dist))
pp <- select(puntos_todos_f, c(X, Y, th, mh, dist))

ss <- cbind(kk, pp)
colnames(ss)
modelo <- lm(ss$mh ~ ss$mh.1)
summary(modelo)
ggplot(ss, aes(x = mh.1, y = mh)) +
  geom_point() +
  geom_smooth(method = "lm", color = "blue", se = FALSE) +
  labs(
    title = "Relación entre mh y mh.1",
    x = "mh.1 (independiente)",
    y = "mh (dependiente)"
  ) +
  theme_minimal()
## PRESENT ----
# Función para procesar los datos de la serie presente
pa_mh_present <- function(j, th = .95) {
  data_p <- data_present_climatic_variables
  mh_p <- data.frame(matrix(1,    
                            nrow = nrow(data_p),
                            ncol = length(names)))
  
  names(mh_p) <- names
  
  for (i in 1:nrow(polygon)){
    pol <- polygon[i,]
    raster_polygon <- terra::mask(terra::crop(present_climatic_variables, pol), pol)
    data_polygon <- terra::as.data.frame(raster_polygon, xy = TRUE)
    data_polygon <- na.omit(data_polygon)
    
    mh <- mahalanobis(data_p[,4:length(data_p)], 
                      colMeans(data_polygon[,3:length(data_polygon)]), 
                      cov(data_p[,4:length(data_p)]), 
                      inverted = F)
    
    mh_p[,i] <- mh
  }
  
  mh_p <- cbind(data_p[,c(1:3)], mh_p)
  mh_raster_p <- terra::rast(mh_p[, c(1:2, j+3)], crs = reference_system)
  names(mh_raster_p) <- colnames(mh_p[j+3])
  plot(mh_raster_p)
  writeRaster(mh_raster_p, paste0(dir_present, "MH_PRESENT_", names[j], ".tif"), overwrite = TRUE)
  assign(paste0("MH_PRESENT_", names[j]), mh_raster_p, envir = .GlobalEnv)
  
  
  puntos_todos_p <- terra::as.points(mh_raster_p)
  puntos_todos_p <- sf::st_as_sf(puntos_todos_p)
  colnames(puntos_todos_p) <- c("mh", "geometry")
  puntos_dentro_p <- sf::st_intersection(puntos_todos_p, pol)
  
  th_mh_p <- quantile(na.omit(puntos_dentro_p$mh), probs = th)
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
  writeRaster(raster, paste0(dir_present, "TH_MH_present_", names[j], ".tif"), overwrite = TRUE)
  assign(paste0("TH_MH_present_", names[j]), mh_raster_p, envir = .GlobalEnv)
  
  
  dist <- sf::st_distance(puntos_todos_p, puntos_dentro_p)
  
  # Para obtener la distancia mínima por cada punto, tomamos el valor mínimo de cada fila
  dist <- apply(dist, 1, min)
  
  # Agregar las distancias mínimas como un atributo a los puntos del raster
  puntos_todos_p$dist <- dist
  
  # Convertir las distancias a kilómetros y redondear
  puntos_todos_p$dist <- round(puntos_todos_p$dist / 1000, 0)
  puntos_todos_p$dist_na <- puntos_todos_p$dist 
  puntos_todos_p$dist_na[puntos_todos_p$dist_na == 0] <- NA
  puntos_todos_p <- cbind(puntos_todos_p, st_coordinates(puntos_todos_p))
  write.csv2(puntos_todos_p, "C:/A_TRABAJO/A_GABRIEL/REPRESENTATIVIDAD/FAST_TEST/GEODA/puntos_todos_p.csv" )
  rango_breaks <- seq(0, max(puntos_todos_p$dist, na.rm = TRUE), by = 10)
  puntos_todos_p <- puntos_todos_p %>%
    mutate(rango_distancia = cut(dist, breaks = rango_breaks, include.lowest = TRUE))
  
  resultados <- puntos_todos_p %>%
    filter(!is.na(rango_distancia)) %>%  # Filtrar filas con rango de distancia válido
    group_by(rango_distancia) %>%  # Agrupar por rango de distancia
    summarise(
      n_total = n(),  # Número total de observaciones en el rango
      n_ceros = sum(th == 0, na.rm = TRUE),  # Número de ceros en 'th'
      n_unos = sum(th == 1, na.rm = TRUE)  # Número de unos en 'th'
    ) %>%
    mutate(
      porcentaje_ceros = (n_ceros / n_total) * 100,  # Porcentaje de ceros
      porcentaje_unos = (n_unos / n_total) * 100  # Porcentaje de unos
    )
  
  
  ggplot(resultados, aes(x = rango_distancia)) +
    geom_point(aes(y = porcentaje_unos, color = "Unos"), size = 3) +  # Puntos para los 'Unos'
    geom_line(aes(y = porcentaje_unos, color = "Unos"), size = 1, group = 1) +  # Línea para conectar los puntos
    labs(
      title = "Distribución de ceros y unos según el rango de distancia",
      x = "Rango de distancia",
      y = "Porcentaje",
      fill = "Valores de mh",
      color = "Valores de mh"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }
  ###############################################################################################################################################################
  
  
pa_mh_present(j)
  
  
  ## SORENSEN
  
  # Instalar y cargar los paquetes necesarios
  install.packages("raster")
  install.packages("ggplot2")
  library(raster)
  library(ggplot2)
  
  
  
  
  # Paso 3: Definir dos subconjuntos de celdas (por ejemplo, puedes tomar dos áreas distintas del raster)
  # Convertir el raster a puntos para extraer sus coordenadas
  puntos_raster <- terra::as.points(raster)
  puntos_raster <- sf::st_as_sf(puntos_raster)
  colnames(puntos_raster) <- c("mh", "geometry")
  
  # Seleccionar las celdas que están dentro del polígono
  celdas_A <- puntos_dentro_p$mh
  
  # Paso 4: Calcular el índice de Sorensen iterativamente por distancia
  # Definir una función para calcular el índice de Sorensen
  calcular_sorensen <- function(celdas_A, celdas_B) {
    interseccion <- length(intersect(rownames(celdas_A), rownames(celdas_B)))  # Número de celdas comunes
    indice_sorensen <- (2 * interseccion) / (nrow(celdas_A) + nrow(celdas_B))  # Índice de Sorensen
    return(indice_sorensen)
  }
  
  # Paso 5: Calcular las distancias entre las celdas de A y otras celdas
  # Calcular distancias entre las celdas de A y todas las demás celdas
  distancias <- puntos_todos_p$dist
  
  # Paso 6: Iterar sobre las distancias y calcular el índice de Sorensen para cada distancia
  # Aquí agrupamos las distancias en intervalos y calculamos el índice de Sorensen para cada intervalo
  
  # Crear un dataframe para guardar los resultados
  resultados <- data.frame(distancia = numeric(), indice_sorensen = numeric())
  
  # Establecer los rangos de distancia (por ejemplo, de 0 a 100 metros, con un intervalo de 10 metros)
  rangos_distancia <- seq(0, max(distancias), by = 10)
  
  # Iterar sobre cada rango de distancia
  for (rango in rangos_distancia) {
    # Seleccionar las celdas dentro del rango de distancia
    celdas_en_rango <- puntos_raster[distancias <= rango, ]
    
    # Calcular el índice de Sorensen entre las celdas de A y las celdas en el rango
    indice <- calcular_sorensen(celdas_A, celdas_en_rango)
    
    # Guardar el resultado
    resultados <- rbind(resultados, data.frame(distancia = rango, indice_sorensen = indice))
  }
  
  # Paso 7: Graficar la variación del índice de Sorensen con la distancia
  ggplot(resultados, aes(x = distancia, y = indice_sorensen)) +
    geom_line() +
    labs(title = "Variación del índice de Sorensen con la distancia",
         x = "Distancia Euclidiana", y = "Índice de Sorensen")
  
  ###############################################################################################################################################################
  library(ncf)
  coords <- st_coordinates(puntos_todos_p)  # Extrae coordenadas (X, Y)
  variable <- puntos_todos_p$tu_variable    # Variable numérica de interés
  
  # 2. Calcular el correlograma espacial
  correlog <- correlog(coords[, 1], coords[, 2], variable, increment = 0.089, resamp = 1)
  
  # 3. Visualizar el correlograma
  plot(correlog, main = "Correlograma Espacial", xlab = "Distancia", ylab = "Autocorrelación")
  
  # Definir los rangos de distancia
  distances <- seq(0, 50, by = 10) # Puedes ajustar el rango y el incremento
  
  # Crear un dataframe vacío para almacenar los resultados
  results <- data.frame(
    Distance = numeric(),
    MoranI = numeric(),
    PValue = numeric()
  )
 
  # Iterar sobre cada distancia
  for (i in 2:length(distances)) {
    # Filtrar los puntos dentro de la distancia actual
    subset_points <- filter(puntos_todos_p, dist > distances[i-1] & dist <= distances[i] )
    
    # Verificar que haya suficientes puntos para calcular Moran's I
    if (nrow(subset_points) > 1) {
      # Obtener las coordenadas
      coords <- st_coordinates(subset_points)
      
      # Crear vecinos espaciales dentro de la distancia actual
      neighbors <- dnearneigh(coords, distances[i-1], distances[i], use_s2 = TRUE, dwithin = TRUE)
      weights <- nb2listw(neighbors, style = "W", zero.policy = TRUE)
      
      # Calcular Moran's I
      moran_result <- moran.mc(subset_points$mh, weights, nsim = 499, zero.policy = TRUE)
      
      # Almacenar los resultados
      results <- rbind(results, data.frame(
        Distance = distances[i],
        MoranI = moran_result$statistic,
        PValue = min(moran_result$p.value, 1 - moran_result$p.value)
      ))
    }
  }
  
  # Graficar los resultados
  ggplot(results, aes(x = Distance, y = MoranI)) +
    geom_line(color = "blue") +
    geom_point(aes(color = PValue < 0.01), size = 3) +
    scale_color_manual(values = c("grey" = "grey", "red" = "red"), guide = "none") +
    geom_text(aes(label = round(PValue, 3)), vjust = -0.5, size = 3) +
    labs(
      x = "Distance (meters)",
      y = "Moran's I",
      title = "Moran's I vs Distance",
      subtitle = "Significant points in red (p < 0.01)"
    ) +
    theme_minimal() +
    theme(
      panel.grid.major = element_line(color = "#EEEEEE"),
      panel.grid.minor = element_blank()
    )
 
  
  
  s.center <- st_point_on_surface(polygon)
  s.coord <- st_coordinates(puntos_todos_p)
  start <- 0 # Starting distance in meters (the From)
  end <- 2000 # Ending distance in meters (the To)
  incr <- 1000 # Distance increment (which also defines the annulus width)
  incr.v <- seq(start, end, incr)
  
  
  # Define empty vector elements to store the I and p-values
  morI.mc <- vector()
  sign.mc <- vector()
  
  # Loop through each distance band
  for (i in (2:length(incr.v))) {
    s.dist <- dnearneigh(s.coord, incr.v[i - 1], incr.v[i])
    s.lw <- nb2listw(s.dist, style = "W", zero.policy=T)
    s.mor <- moran.mc(puntos_todos_p$mh, s.lw, nsim=5, zero.policy = TRUE)
    sign.mc[i] <- s.mor$p.value
    morI.mc[i] <- s.mor$statistic
  }
  
  # Modify p-value to reflect extremes at each end
  sign.mc <- ifelse(sign.mc > 0.5, 1 - sign.mc, sign.mc)
  # First, generate an empty plot
  plot(morI.mc ~ eval(incr.v - incr * 0.5), type = "n", ann = FALSE, axes = FALSE)
  
  # Set the background plot to grey then add white grids
  u <- par("usr") # Get plot are coordinates
  rect(u[1], u[3], u[2], u[4], col = "#EEEEEE", border = NA)
  axis(1, lab = ((incr.v) / 1000), at = (incr.v), tck = 1, col = "#FFFFFF", lty = 1)
  axis(2, tck = 1, col = "#FFFFFF", lty = 1, labels = FALSE)
  
  # Add the theoretical "no autocorelation" line
  abline(h = -1 / (length(s$Income)), col = "grey20")
  
  # Add the plot to the canvas
  par(new = TRUE)
  plot(morI.mc ~ eval(incr.v - incr * 0.5),
       type = "b", xaxt = "n", las = 1,
       xlab = "Distance (km)", ylab = "Moran's I")
  points(morI.mc ~ eval(incr.v - incr * 0.5), 
         col = ifelse(sign.mc < 0.01, "red", "grey"), 
         pch = 16, cex = 2.0)
  
  # Add numeric values to points
  text(eval(incr.v - incr * 0.5), morI.mc, round(sign.mc,3), pos = 3, cex = 0.5)
  
  
  
}

# calculate great cirlce distances with fields package
library(fields)
# popdists matrix
popdists <- as.matrix(rdist.earth(cbind(puntos_todos_p$lon, puntos_todos_p$lat), miles = F, R = NULL))
diag(popdists) <- 0

# autocorr function
autocorr <- function(w,x,dist=500){
  aa <- ceiling(max(w)/dist)
  dists <- seq(0,aa*dist,dist)
  cors <- NULL
  for(i in 1:aa){
    w1 <- ifelse(w > dists[i] & w <= dists[i+1], 1, 0) 
    w2 <- w1
    for(j in 1:dim(w1)[1]){
      nu <- sum(w1[j,])
      if(nu>0){
        w2[j,] <- w1[j,]/nu
      }  
    }
    lag <- w2 %*% x
    cors <- c(cors,cor(x,lag))
  }
  return(cors)
}


library(ggplot2)
library(dplyr)

# Calcular el punto central de los polígonos
s.center <- st_point_on_surface(polygon)

# Extraer las coordenadas
s.coord <- st_coordinates(puntos_todos_p)

# Definir el rango de distancias y el incremento
start <- 0 # Distancia inicial (en metros)
end <- 100 # Distancia final (en metros)
incr <- 50 # Incremento de distancia (ancho del anillo)
incr.v <- seq(start, end, by = incr)

# Crear un dataframe para almacenar los resultados
results <- data.frame(
  Distance = numeric(),
  MoranI = numeric(),
  PValue = numeric()
)

# Iterar sobre cada banda de distancia
for (i in 2:length(incr.v)) {
  # Crear vecinos por distancia
  s.dist <- dnearneigh(s.coord, d1 = incr.v[i - 1], d2 = incr.v[i])
  
  # Convertir a lista de pesos
  s.lw <- nb2listw(s.dist, style = "W", zero.policy = TRUE)
  
  # Calcular Moran's I con simulación Monte Carlo
  s.mor <- moran.test(puntos_todos_p$mh, s.lw,  zero.policy = TRUE)
  
  # Guardar resultados
  results <- rbind(results, data.frame(
    Distance = (incr.v[i - 1] + incr.v[i]) / 2,
    MoranI = s.mor$statistic,
    PValue = min(s.mor$p.value, 1 - s.mor$p.value)
  ))
}

# Graficar los resultados con ggplot2
ggplot(results, aes(x = Distance / 1000, y = MoranI)) +
  geom_hline(yintercept = -1 / length(puntos_todos_p$mh), linetype = "dashed", color = "grey20") +
  geom_line(color = "blue") +
  geom_point(aes(color = PValue < 0.01), size = 3) +
  scale_color_manual(values = c("grey" = "grey", "red" = "red"), guide = "none") +
  geom_text(aes(label = round(PValue, 3)), vjust = -0.5, size = 3) +
  labs(
    x = "Distance (km)",
    y = "Moran's I",
    title = "Moran's I vs Distance",
    subtitle = "Significant points in red (p < 0.01)"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_line(color = "#EEEEEE"),
    panel.grid.minor = element_blank()
  )
write.csv2(puntos_todos_p, "C:/A_TRABAJO/A_GABRIEL/REPRESENTATIVIDAD/FAST_TEST/GEODA/PUNTOS_TODOS.csv" )

