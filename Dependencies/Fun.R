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
  writeRaster(mh_raster_p, paste0(dir_present, "MH_present_", names[j], ".tif"), overwrite = TRUE)
  assign(paste0("MH_present_", names[j]), mh_raster_p, envir = .GlobalEnv)
  
  
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
pa_mh_present_future <- function(j) {
  
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
  writeRaster(mh_raster_p, paste0(dir_present, "MH_present_", names[j], ".tif"), overwrite = TRUE)
  assign(paste0("MH_present_", names[j]), mh_raster_p, envir = .GlobalEnv)
  
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
  writeRaster(raster, paste0(dir_present, "TH_MH_present", names[j], ".tif"), overwrite = TRUE)
  assign(paste0("TH_MH_present_", names[j]), mh_raster_p, envir = .GlobalEnv)
  
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
