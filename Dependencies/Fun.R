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

library(tidyterra)
# Load packages ----
packages.to.use <- c("corrplot", "terra", "tictoc", "tidyverse","sf", "pracma", "tidyterra")

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
pa_mh_present_future <- function(j, th = .9) {
  color_map <- c("0" = "coral3", "1" = "aquamarine3")
  
  
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
  plot(raster, main = paste0("Present mh distance threshold of ", names[j]), col = color_map[as.character(sort(unique(values(raster))))])
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
  
  # UMBRAL presente futuro
  mh_raster_p_f_p <- dplyr::filter(mh_p_f, Period == "Present")
  mh_raster_p_f_p <- terra::rast(mh_raster_p_f_p[, c(1:2, 4)], crs = reference_system)
  names(mh_raster_p_f_p) <- colnames(mh_p_f[4])
  
  puntos_todos_f_p <- terra::as.points(mh_raster_p_f_p)
  puntos_todos_f_p <- sf::st_as_sf(puntos_todos_f_p)
  colnames(puntos_todos_f_p) <- c("mh", "geometry")
  puntos_dentro_f_p <- sf::st_intersection(puntos_todos_f_p, pol)
  
  th_mh_p_f <- quantile(na.omit(puntos_dentro_f_p$mh), probs = th)
  
  ####
  
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
    puntos_todos_f$mh > th_mh_p_f ~ 0,
    puntos_todos_f$mh <= th_mh_p_f  ~ 1
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
  plot(raster, main = paste0(year," ",model," mh distance threshold of ", names[j]), col = color_map[as.character(sort(unique(values(raster))))])
  writeRaster(raster,
              paste0(dir_future, "TH_MH_", model, "_", year, "_", names[j], ".tif"),
              overwrite = TRUE)
  
  assign(paste0("TH_MH_", model, "_", year, "_", names[j]), mh_raster_p, envir = .GlobalEnv)
  
  
  dist <- sf::st_distance(puntos_todos_p, pol)
  
  # Para obtener la distancia mínima por cada punto, tomamos el valor mínimo de cada fila
  dist <- apply(dist, 1, min)
  
  # Agregar las distancias mínimas como un atributo a los puntos del raster
  puntos_todos_p$dist <- round(dist/ 1000, 0)
  dist <- sf::st_distance(puntos_todos_f, pol)
  
  # Para obtener la distancia mínima por cada punto, tomamos el valor mínimo de cada fila
  dist <- apply(dist, 1, min)
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
      n_total_p = n(),  # Número total de observaciones en el rango
      n_ceros_p = sum(th == 0, na.rm = TRUE),  # Número de ceros en 'th'
      n_unos_p = sum(th == 1, na.rm = TRUE)  # Número de unos en 'th'
    ) %>%
    mutate(
      porcentaje_ceros_p = (n_ceros_p / n_total_p) * 100,  # Porcentaje de ceros
      porcentaje_unos_p = (n_unos_p / n_total_p) * 100  # Porcentaje de unos
    )
  
  resultados_f <- puntos_todos_f %>%
    filter(!is.na(dist)) %>%  # Filtrar filas con rango de distancia válido
    group_by(dist) %>%  # Agrupar por rango de distancia
    summarise(
      n_total_f = n(),  # Número total de observaciones en el rango
      n_ceros_f = sum(th == 0, na.rm = TRUE),  # Número de ceros en 'th'
      n_unos_f = sum(th == 1, na.rm = TRUE)  # Número de unos en 'th'
    ) %>%
    mutate(
      porcentaje_ceros_f = (n_ceros_f / n_total_f) * 100,  # Porcentaje de ceros
      porcentaje_unos_f = (n_unos_f / n_total_f) * 100  # Porcentaje de unos
    )
  
  
  resultados_p_df <- st_drop_geometry(resultados_p)  # Elimina la geometría
  resultados_f_df <- st_drop_geometry(resultados_f)
  
  # Unir los datos y asegurarnos que las columnas sean correctas
  resultados_comb <- resultados_p_df %>%
    inner_join(
      resultados_f_df,
      by = "dist"
    )
  resultados_comb <- resultados_comb %>%
    mutate(
      n_unos_p_acum = cumsum(n_unos_p),  # Acumulación de n_unos_p
      n_unos_f_acum = cumsum(n_unos_f)   # Acumulación de n_unos_f
    )
  max_n_unos_p_acum <- max(resultados_comb$n_unos_p_acum, na.rm = TRUE)
  factor_escala <- max_n_unos_p_acum / 100
  
  area_presente <- round(trapz(resultados_comb$dist, resultados_comb$porcentaje_unos_p/100), 2)
  area_futuro <- round(trapz(resultados_comb$dist, resultados_comb$porcentaje_unos_f/100),2)
  area_compartida <- round(trapz(resultados_comb$dist, pmin(resultados_comb$porcentaje_unos_p/100, resultados_comb$porcentaje_unos_f/100)),2)
  
    # Graficar con los dos ejes
  p <- ggplot() +
    # Área para los valores "Present"
    geom_area(data = resultados_comb, aes(x = dist, y = porcentaje_unos_p, fill = "Present"), alpha = 0.5, size = 3) +
    # Área para los valores "Future"
    geom_area(data = resultados_comb, aes(x = dist, y = porcentaje_unos_f, fill = "Future"), alpha = 0.5, size = 3) +
    # Línea para n_unos_p acumulado
    geom_line(data = resultados_comb, aes(x = dist, y = n_unos_p_acum / factor_escala, group = 1), color = "aquamarine4", size = 1) + 
    # Línea para n_unos_f acumulado
    geom_line(data = resultados_comb, aes(x = dist, y = n_unos_f_acum / factor_escala, group = 2), color = "coral3", size = 1) + 
    # Línea horizontal escalada al eje derecho
    geom_segment(
      aes(
        x = min(resultados_comb$dist),  # Inicio del segmento (mínima distancia)
        xend = max(resultados_comb$dist), # Fin del segmento (máxima distancia)
        y = nrow(puntos_dentro_p) / factor_escala, # Escalado para eje derecho
        yend = nrow(puntos_dentro_p) / factor_escala # Escalado para eje derecho
      ),
      color = "black", 
      linetype = "dashed",
      size = 1
    ) +
    # Títulos y etiquetas
    labs(
      title = names[j],
      subtitle = paste0(year, "_", model),
      x = "Distancia (km)",
      y = "Porcentaje receptoras",
      fill = "Escenario",
      color = "Escenario",
      caption = paste0("AUC Presente = ", area_presente, 
                       "   AUC Futuro = ", area_futuro, 
                       "   AUC Compartida = ", area_compartida)
    ) +
    # Estilo minimalista
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_manual(values = c("Present" = "aquamarine3", "Future" = "coral3")) +
    # Segunda escala en el eje y derecho para la diferencia acumulada (AUC)
    scale_y_continuous(
      name = "Porcentaje receptoras",  # Eje izquierdo
      limits = c(0, 100),  # Aseguramos que el rango en Y izquierdo esté de 0 a 100
      sec.axis = sec_axis(
        trans = ~ . * factor_escala,  # Escala de transformación ajustada para el eje derecho
        name = "Celdas receptoras acumuladas",  # Eje derecho
        labels = scales::comma
      )  # Eje derecho con la transformación aplicada
    )
  ggsave(filename = paste0(dir_charts, names[j], "_scenarios_difference.jpeg"), plot = p, width = 10, height = 8)
  print(p)
  
}  


## PRESENT-FUTURE 2----
# Función para procesar los datos de la serie presente
pa_mh_present_future2 <- function(j, th = .9) {
  color_map <- c("0" = "coral3", "1" = "aquamarine3")
  
  
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
  lines(pol, add = T, col = "red")
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
  
  raster_th_p <- terra::rasterize(puntos_vect_p, raster_template, field = "th")
  crs(raster_th_p) <- crs(mh_raster_p)
  plot(raster_th_p, main = paste0("Present mh distance threshold of ", names[j]), col = color_map[as.character(sort(unique(values(raster_th_p))))])
  lines(pol, add = T, col = "black")
  writeRaster(raster_th_p, paste0(dir_present, "TH_MH_PRESENT_", names[j], ".tif"), overwrite = TRUE)
  assign(paste0("TH_MH_PRESENT_", names[j]), raster_th_p, envir = .GlobalEnv)
  
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
  
  # UMBRAL presente futuro
  mh_raster_p_f_p <- dplyr::filter(mh_p_f, Period == "Present")
  mh_raster_p_f_p <- terra::rast(mh_raster_p_f_p[, c(1:2, 4)], crs = reference_system)
  names(mh_raster_p_f_p) <- colnames(mh_p_f[4])
  
  puntos_todos_f_p <- terra::as.points(mh_raster_p_f_p)
  puntos_todos_f_p <- sf::st_as_sf(puntos_todos_f_p)
  colnames(puntos_todos_f_p) <- c("mh", "geometry")
  puntos_dentro_f_p <- sf::st_intersection(puntos_todos_f_p, pol)
  
  th_mh_p_f <- quantile(na.omit(puntos_dentro_f_p$mh), probs = th)
  
  ####
  
  mh_raster_p_f <- dplyr::filter(mh_p_f, Period == "Future")
  mh_raster_p_f <- terra::rast(mh_raster_p_f[, c(1:2, 4)], crs = reference_system)
  names(mh_raster_p_f) <- colnames(mh_p_f[4])
  plot(mh_raster_p_f, main = paste0(year," ",model," Future mh distance of ", names[j]))
  lines(pol, add = T, col = "red")
  writeRaster(mh_raster_p_f,
              paste0(dir_future, "MH_", model, "_", year, "_", names[j], ".tif"),
              overwrite = TRUE)
  
  puntos_todos_f <- terra::as.points(mh_raster_p_f)
  puntos_todos_f <- sf::st_as_sf(puntos_todos_f)
  colnames(puntos_todos_f) <- c("mh", "geometry")
  puntos_dentro_f <- sf::st_intersection(puntos_todos_f, pol)
  
  
  puntos_todos_f$th <- case_when(
    puntos_todos_f$mh > th_mh_p_f ~ 0,
    puntos_todos_f$mh <= th_mh_p_f  ~ 1
  )
  puntos_todos_f$th <- as.numeric(puntos_todos_f$th)
  
  res <- res(mh_raster_p_f)
  bbox <- ext(mh_raster_p_f)
  
  nrows <- round((bbox[4] - bbox[3]) / res[2])
  ncols <- round((bbox[2] - bbox[1]) / res[1])
  raster_template <- rast(ext = bbox, nrows = nrows, ncols = ncols)
  puntos_vect_f <- vect(puntos_todos_f)
  
  raster_th_f <- terra::rasterize(puntos_vect_f, raster_template, field = "th")
  crs(raster_th_f) <- crs(mh_raster_p_f)
  plot(raster_th_f, main = paste0(year," ",model," mh distance threshold of ", names[j]), col = color_map[as.character(sort(unique(values(raster_th_f))))])
  lines(pol, add = T, col = "black")
  writeRaster(raster_th_f,
              paste0(dir_future, "TH_MH_", model, "_", year, "_", names[j], ".tif"),
              overwrite = TRUE)
  
  assign(paste0("TH_MH_", model, "_", year, "_", names[j]), raster_th_f, envir = .GlobalEnv)
  
  
  # 1. Superficie compartida (1 en ambos rasters)
  compartida <- raster_th_p * raster_th_f
  
  # 2. Superficie en raster1 pero no en raster2
  solo_presente <- (raster_th_p == 1) * (raster_th_f == 0)
  
  # 3. Superficie en raster2 pero no en raster1
  solo_futura <- (raster_th_p == 0) * (raster_th_f == 1)
  
  # Crear un raster final con todas las categorías
  raster_th_p_f <- compartida + (solo_presente * 2) + (solo_futura * 3)
  
  resultado_factor <- as.factor(raster_th_p_f)
  
  l <- ggplot() +
    geom_spatraster(data = resultado_factor) +
    geom_sf(data = pol, color = "black", fill = NA)+
    scale_fill_manual(
      values = c("0" = "grey90", 
                 "1" = "gold", 
                 "2" = "aquamarine3", 
                 "3" = "coral3"),
      labels = c("0" = "Zonas sin representatividad",
                 "1" = "Representatividad estable",
                 "2" = "Representatividad presente",
                 "3" = "Representatividad futura"),
      na.value = "transparent"
    ) +
    coord_sf() +
    theme_minimal() +
    labs(title = paste0(names[j]),
         fill = "Valor")
  
  ggsave(filename = paste0(dir_charts, names[j], "_rep_shared.jpeg"), plot = l, width = 10, height = 8)
  print(l)
  
  crs(raster_th_p_f) <- crs(mh_raster_p_f)
  writeRaster(raster_th_p_f,
              paste0(dir_future, "TH_MH_PRESENT_", model, "_", year, "_", names[j], ".tif"),
              overwrite = TRUE)
  
  assign(paste0("TH_MH_", model, "_", year, "_", names[j]), raster_th_f, envir = .GlobalEnv)  
  dist <- sf::st_distance(puntos_todos_p, pol)
  
  # Para obtener la distancia mínima por cada punto, tomamos el valor mínimo de cada fila
  dist <- apply(dist, 1, min)
  
  # Agregar las distancias mínimas como un atributo a los puntos del raster
  puntos_todos_p$dist <- round(dist/ 1000, 0)
  dist <- sf::st_distance(puntos_todos_f, pol)
  
  # Para obtener la distancia mínima por cada punto, tomamos el valor mínimo de cada fila
  dist <- apply(dist, 1, min)
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
      n_total_p = n(),  # Número total de observaciones en el rango
      n_ceros_p = sum(th == 0, na.rm = TRUE),  # Número de ceros en 'th'
      n_unos_p = sum(th == 1, na.rm = TRUE)  # Número de unos en 'th'
    ) %>%
    mutate(
      porcentaje_ceros_p = (n_ceros_p / n_total_p) * 100,  # Porcentaje de ceros
      porcentaje_unos_p = (n_unos_p / n_total_p) * 100  # Porcentaje de unos
    )
  
  resultados_f <- puntos_todos_f %>%
    filter(!is.na(dist)) %>%  # Filtrar filas con rango de distancia válido
    group_by(dist) %>%  # Agrupar por rango de distancia
    summarise(
      n_total_f = n(),  # Número total de observaciones en el rango
      n_ceros_f = sum(th == 0, na.rm = TRUE),  # Número de ceros en 'th'
      n_unos_f = sum(th == 1, na.rm = TRUE)  # Número de unos en 'th'
    ) %>%
    mutate(
      porcentaje_ceros_f = (n_ceros_f / n_total_f) * 100,  # Porcentaje de ceros
      porcentaje_unos_f = (n_unos_f / n_total_f) * 100  # Porcentaje de unos
    )
  
  
  resultados_p_df <- st_drop_geometry(resultados_p)  # Elimina la geometría
  resultados_f_df <- st_drop_geometry(resultados_f)
  
  # Unir los datos y asegurarnos que las columnas sean correctas
  resultados_comb <- resultados_p_df %>%
    inner_join(
      resultados_f_df,
      by = "dist"
    )
  resultados_comb <- resultados_comb %>%
    mutate(
      n_unos_p_acum = cumsum(n_unos_p),  # Acumulación de n_unos_p
      n_unos_f_acum = cumsum(n_unos_f)   # Acumulación de n_unos_f
    )
  max_n_unos_p_acum <- max(resultados_comb$n_unos_p_acum, na.rm = TRUE)
  factor_escala <- max_n_unos_p_acum / 100
  
  area_presente <- round(trapz(resultados_comb$dist, resultados_comb$porcentaje_unos_p/100), 2)
  area_futuro <- round(trapz(resultados_comb$dist, resultados_comb$porcentaje_unos_f/100),2)
  area_compartida <- round(trapz(resultados_comb$dist, pmin(resultados_comb$porcentaje_unos_p/100, resultados_comb$porcentaje_unos_f/100)),2)
  
  # Graficar con los dos ejes
  p <- ggplot() +
    # Área para los valores "Present"
    geom_area(data = resultados_comb, aes(x = dist, y = porcentaje_unos_p, fill = "Present"), alpha = 0.5, size = 3) +
    # Área para los valores "Future"
    geom_area(data = resultados_comb, aes(x = dist, y = porcentaje_unos_f, fill = "Future"), alpha = 0.5, size = 3) +
    # Línea para n_unos_p acumulado
    geom_line(data = resultados_comb, aes(x = dist, y = n_unos_p_acum / factor_escala, group = 1), color = "aquamarine4", size = 1) + 
    # Línea para n_unos_f acumulado
    geom_line(data = resultados_comb, aes(x = dist, y = n_unos_f_acum / factor_escala, group = 2), color = "coral3", size = 1) + 
    # Línea horizontal escalada al eje derecho
    geom_segment(
      aes(
        x = min(resultados_comb$dist),  # Inicio del segmento (mínima distancia)
        xend = max(resultados_comb$dist), # Fin del segmento (máxima distancia)
        y = nrow(puntos_dentro_p) / factor_escala, # Escalado para eje derecho
        yend = nrow(puntos_dentro_p) / factor_escala # Escalado para eje derecho
      ),
      color = "black", 
      linetype = "dashed",
      size = 1
    ) +
    # Títulos y etiquetas
    labs(
      title = names[j],
      subtitle = paste0(year, "_", model),
      x = "Distancia (km)",
      y = "Porcentaje receptoras",
      fill = "Escenario",
      color = "Escenario",
      caption = paste0("AUC Presente = ", area_presente, 
                       "   AUC Futuro = ", area_futuro, 
                       "   AUC Compartida = ", area_compartida)
    ) +
    # Estilo minimalista
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_manual(values = c("Present" = "aquamarine3", "Future" = "coral3")) +
    # Segunda escala en el eje y derecho para la diferencia acumulada (AUC)
    scale_y_continuous(
      name = "Porcentaje receptoras",  # Eje izquierdo
      limits = c(0, 100),  # Aseguramos que el rango en Y izquierdo esté de 0 a 100
      sec.axis = sec_axis(
        trans = ~ . * factor_escala,  # Escala de transformación ajustada para el eje derecho
        name = "Celdas receptoras acumuladas",  # Eje derecho
        labels = scales::comma
      )  # Eje derecho con la transformación aplicada
    )
  ggsave(filename = paste0(dir_charts, names[j], "_scenarios_difference.jpeg"), plot = p, width = 10, height = 8)
  print(p)
  
}  
