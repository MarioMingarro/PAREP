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
  
  
  pol <- polygon[j, ]
  raster_polygon <- terra::mask(terra::crop(present_climatic_variables, pol), pol)
  data_polygon <- terra::as.data.frame(raster_polygon, xy = TRUE)
  data_polygon <- na.omit(data_polygon)
  
  mh <- mahalanobis(data_p[, 4:length(data_p)], colMeans(data_polygon[, 3:length(data_polygon)]), cov(data_p[, 4:length(data_p)]), inverted = F)
  
  mh_p <- cbind(data_p[,c(1:3)], mh)
  mh_raster_p <- terra::rast(mh_p[, c(1:2, j+3)], crs = reference_system)
  names(mh_raster_p) <- colnames(mh_p[j+3])
  plot(mh_raster_p, main = paste0( "Present mh distance of ", names[j]))
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
  plot(raster, main = paste0("Present mh distance threshold of ", names[j]))
  writeRaster(raster, paste0(dir_present, "TH_MH_PRESENT_", names[j], ".tif"), overwrite = TRUE)
  assign(paste0("TH_MH_PRESENT_", names[j]), mh_raster_p, envir = .GlobalEnv)
  
  # Futuro
  data_p_f <- rbind(data_present_climatic_variables, data_future_climatic_variables)
  
  mh_p_f <- data.frame(matrix(1,    
                            nrow = nrow(data_p_f),
                            ncol = length(names)))
  
  names(mh_p_f) <- names
  
  raster_polygon <- terra::mask(terra::crop(present_climatic_variables, pol), pol)
  data_polygon <- terra::as.data.frame(raster_polygon, xy = TRUE)
  data_polygon <- na.omit(data_polygon)
  
  mh <- mahalanobis(data_p_f[, 4:length(data_p_f)], colMeans(data_polygon[, 3:length(data_polygon)]), cov(data_p_f[, 4:length(data_p_f)]), inverted = F)
  
  
  mh_p_f <- cbind(data_p_f[, c(1:3)], mh)
  
  mh_raster_p_f <- dplyr::filter(mh_p_f, Period == "Future")
  mh_raster_p_f <- terra::rast(mh_raster_p_f[, c(1:2, j+3)], crs = reference_system)
  names(mh_raster_p_f) <- colnames(mh_p_f[j+3])
  plot(mh_raster_p_f, main = paste0("Present mh distance of ", names[j]))
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
pa_mh_present_future2 <- function(j, th = .95) {
  
  # Presente
  data_p <- data_present_climatic_variables
  
  mh_p <- data.frame(matrix(1,    
                            nrow = nrow(data_p),
                            ncol = length(names)))
  
  names(mh_p) <- names
  
  
  pol <- polygon[j, ]
  raster_polygon <- terra::mask(terra::crop(present_climatic_variables, pol), pol)
  data_polygon <- terra::as.data.frame(raster_polygon, xy = TRUE)
  data_polygon <- na.omit(data_polygon)
  
  mh <- mahalanobis(data_p[, 4:length(data_p)], colMeans(data_polygon[, 3:length(data_polygon)]), cov(data_p[, 4:length(data_p)]), inverted = F)
  
  mh_p <- cbind(data_p[,c(1:3)], mh)
  mh_raster_p <- terra::rast(mh_p[, c(1:2, 4)], crs = reference_system)
  names(mh_raster_p) <- colnames(mh_p[4])
  plot(mh_raster_p, main = paste0( "Present mh distance of ", names[j]))
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
  plot(raster, main = paste0("Present mh distance threshold of ", names[j]))
  writeRaster(raster, paste0(dir_present, "TH_MH_PRESENT_", names[j], ".tif"), overwrite = TRUE)
  assign(paste0("TH_MH_PRESENT_", names[j]), mh_raster_p, envir = .GlobalEnv)
  
  # Futuro
  data_p_f <- rbind(data_present_climatic_variables, data_future_climatic_variables)
  
  mh_p_f <- data.frame(matrix(1,    
                              nrow = nrow(data_p_f),
                              ncol = length(names)))
  
  names(mh_p_f) <- names
  
  raster_polygon <- terra::mask(terra::crop(present_climatic_variables, pol), pol)
  data_polygon <- terra::as.data.frame(raster_polygon, xy = TRUE)
  data_polygon <- na.omit(data_polygon)
  
  mh <- mahalanobis(data_p_f[, 4:length(data_p_f)], colMeans(data_polygon[, 3:length(data_polygon)]), cov(data_p_f[, 4:length(data_p_f)]), inverted = F)
  
  
  mh_p_f <- cbind(data_p_f[, c(1:3)], mh)
  
  mh_raster_p_f <- dplyr::filter(mh_p_f, Period == "Future")
  mh_raster_p_f <- terra::rast(mh_raster_p_f[, c(1:2, 4)], crs = reference_system)
  names(mh_raster_p_f) <- colnames(mh_p_f[4])
  plot(mh_raster_p_f, main = paste0(year," ",model," Future mh distance of ", names[j]))
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
  plot(raster, main = paste0(year," ",model," mh distance threshold of ", names[j]))
  writeRaster(raster,
              paste0(dir_future, "TH_MH_", model, "_", year, "_", names[j], ".tif"),
              overwrite = TRUE)
  
  assign(paste0("TH_MH_", model, "_", year, "_", names[j]), mh_raster_p, envir = .GlobalEnv)
  
  
  dist <- sf::st_distance(puntos_todos_p, pol)
  
  # Para obtener la distancia mínima por cada punto, tomamos el valor mínimo de cada fila
  dist <- apply(dist, 1, min)
  
  # Agregar las distancias mínimas como un atributo a los puntos del raster
  puntos_todos_p$dist <- round(dist/ 1000, 0)
  puntos_todos_f$dist <- round(dist/ 1000, 0)
  
  
  puntos_todos_p <- cbind(puntos_todos_p, st_coordinates(puntos_todos_p))
  puntos_todos_f <- cbind(puntos_todos_f, st_coordinates(puntos_todos_f))
  #write.csv2(puntos_todos_f, "C:/A_TRABAJO/A_GABRIEL/REPRESENTATIVIDAD/FAST_TEST/GEODA/puntos_todos_p.csv" )
  # rango_breaks <- seq(0, max(puntos_todos_p$dist, na.rm = TRUE), by = 10)
  # 
  # puntos_todos_p <- puntos_todos_p %>%
  #   mutate(rango_distancia = cut(dist, breaks = rango_breaks, include.lowest = TRUE))
  # puntos_todos_f <- puntos_todos_f %>%
  #   mutate(rango_distancia = cut(dist, breaks = rango_breaks, include.lowest = TRUE))
  
  resultados_p <- puntos_todos_p %>%
    filter(!is.na(dist)) %>%  # Filtrar filas con rango de distancia válido
    group_by(dist) %>%  # Agrupar por rango de distancia
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
    filter(!is.na(dist)) %>%  # Filtrar filas con rango de distancia válido
    group_by(dist) %>%  # Agrupar por rango de distancia
    summarise(
      n_total = n(),  # Número total de observaciones en el rango
      n_ceros = sum(th == 0, na.rm = TRUE),  # Número de ceros en 'th'
      n_unos = sum(th == 1, na.rm = TRUE)  # Número de unos en 'th'
    ) %>%
    mutate(
      porcentaje_ceros = (n_ceros / n_total) * 100,  # Porcentaje de ceros
      porcentaje_unos = (n_unos / n_total) * 100  # Porcentaje de unos
    )
  
  
  # Calcular el área bajo la curva acumulada usando la regla del trapecio
  calc_auc_acumulado <- function(x, y) {
    auc_acum <- cumsum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
    return(c(0, auc_acum))  # Añadimos un 0 al inicio para comenzar desde 0
  }
  
  resultados_p_df <- st_drop_geometry(resultados_p)  # Elimina la geometría
  resultados_f_df <- st_drop_geometry(resultados_f)
  
  # Unir los datos y asegurarnos que las columnas sean correctas
  resultados_comb <- resultados_p_df %>%
    rename(presente = porcentaje_unos) %>%
    inner_join(
      resultados_f_df %>% rename(futuro = porcentaje_unos),
      by = "dist"
    )
  
  # Asegurarnos de que no haya valores NA en las columnas clave
  resultados_comb <- resultados_comb %>%
    filter(!is.na(presente) & !is.na(futuro))
  
  # Calcular la diferencia acumulada de AUC
  resultados_comb <- resultados_comb %>%
    mutate(
      auc_presente_acum = calc_auc_acumulado(dist, presente),
      auc_futuro_acum = calc_auc_acumulado(dist, futuro),
      auc_dif_acum = abs(auc_presente_acum - auc_futuro_acum)
    )
  
  # Calculamos el factor de escala para auc_dif_acum
  max_auc <- max(resultados_comb$auc_dif_acum, na.rm = TRUE)
  max_porcentaje <- max(resultados_comb$presente, na.rm = TRUE)
  
  # El factor de escala para el eje derecho (densidad)
  factor_escala <- max_auc / max_porcentaje
  
  # Graficar con los dos ejes
  p <- ggplot() +
    # Área para los valores "Present"
    geom_area(data = resultados_comb, aes(x = dist, y = presente, fill = "Present"), alpha = 0.5, size = 3) +
    # Área para los valores "Future"
    geom_area(data = resultados_comb, aes(x = dist, y = futuro, fill = "Future"), alpha = 0.5, size = 3) +
    # Línea para la diferencia acumulada de densidad (AUC), escalada
    geom_line(data = resultados_comb, aes(x = dist, y = auc_dif_acum / factor_escala, group = 1), color = "blue", size = 1) + 
    # Títulos y etiquetas
    labs(
      title = "% de areas receptoras acorde a distancia y diferencia acumulada entre escenarios",
      x = "Distancia (km)",
      y = "Porcentaje receptoras",
      fill = "Escenario",
      color = "Escenario"
    ) +
    # Estilo minimalista
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_manual(values = c("Present" = "green", "Future" = "red")) +
    # Segunda escala en el eje y derecho para la diferencia acumulada (AUC)
    scale_y_continuous(
      name = "Porcentaje receptoras",  # Eje izquierdo
      limits = c(0, 100),  # Aseguramos que el rango en Y izquierdo esté de 0 a 100
      sec.axis = sec_axis(
        trans = ~ . * factor_escala,  # Escala de transformación ajustada para el eje derecho
        name = "Diferencia Acumulada",  # Eje derecho
        labels = scales::comma
      )  # Eje derecho con la transformación aplicada
    )
  print(p)
  
}






