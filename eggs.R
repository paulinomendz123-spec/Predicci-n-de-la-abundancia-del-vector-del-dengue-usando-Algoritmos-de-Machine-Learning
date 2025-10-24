# Predicción de la abundancia del vector del dengue usando algoritmos de Machine Learning

# Script para la Integración de Covariables Geoespaciales y Socio-Demográficas
# Autor: Benjamín Paulino Mendoza Contreras y Mariana Ramírez Hernández 
# Fecha: 23/10/2025

## 1. Configuración Inicial y Carga de Librerías

# Instala paquetes si es necesario: install.packages(c("sf", "terra", "geodata", "readr", "dplyr", "stringr", "lubridate", "worldmet"))
library(sf)
library(terra)
library(geodata)
library(readr)
library(dplyr)
library(stringr)
library(lubridate)
library(worldmet)

# Crear carpeta para guardar datos (se ignorará si ya existe)
dir.create("data_covariables", showWarnings = FALSE)

## 2. Carga y Preparación de Datos de Ovitrampas
eggs_data <- read_csv("eggs_data.csv") # Asegúrate de que este archivo esté en el directorio correcto.
eggs_data <- eggs_data %>%
  rename(lon = x, lat = y) %>%
  mutate(id = 1:n())

# Convertir a objeto espacial sf
eggs_sf <- st_as_sf(eggs_data, coords = c("lon", "lat"), crs = 4326)

# Definir los límites de recorte (Bbox de Yucatán, México y alrededores)
yuc_bounds <- ext(-91.5, -87.5, 18.0, 22.5)

## 3. Covariables Estáticas Globales (Elevación, Cobertura del Suelo, Temp. Media Histórica)

# Función de extracción generalizada para evitar repetición
extract_and_bind <- function(var_name, eggs_sf, eggs_data, bounds, new_name, path = "data_covariables/") {
  rast <- geodata::worldclim_global(var = var_name, res = 0.5, path = path)
  rast_crop <- crop(rast, bounds)
  extract_val <- extract(rast_crop, eggs_sf, ID = FALSE)
  names(extract_val) <- new_name
  eggs_data <- cbind(eggs_data, extract_val)
  return(eggs_data)
}

# 3.1. Elevación (SRTM)
eggs_data <- extract_and_bind("elev", eggs_sf, eggs_data, yuc_bounds, "elev_srtm")

# 3.2. Cobertura del Suelo (Built-Up Fraction)
eggs_data <- extract_and_bind("built", eggs_sf, eggs_data, yuc_bounds, "built_frac")

# 3.3. Temperatura Media Histórica (TAVG)
eggs_data <- extract_and_bind("tavg", eggs_sf, eggs_data, yuc_bounds, "temp_media_hist")

# Limpieza y selección de variables estáticas
eggs_data_limpio <- eggs_data %>%
  select(
    id, year, week, lon, lat, eggs,
    elev_srtm, built_frac, temp_media_hist
  )

write_csv(eggs_data_limpio, "data_covariables/eggs_data_parcial_covariables.csv")
cat("✅ Covariables estáticas (elevación, built, tavg) añadidas y exportadas.\n")


## 4. Integración de Covariables Sociodemográficas (INEGI - Por Localidad)

# Cargar la base de datos limpia de la Sección 3
eggs_data_limpio <- read_csv("data_covariables/eggs_data_parcial_covariables.csv", show_col_types = FALSE)
eggs_sf_limpio <- st_as_sf(eggs_data_limpio, coords = c("lon", "lat"), crs = 4326)
ruta_ageb_shp <- "ruta/a/tu/shapefile_ageb.shp" # ¡ACTUALIZA ESTA RUTA!
ruta_censo_csv <- "ruta/a/tu/ITER_CENSO_2020.csv" # ¡ACTUALIZA ESTA RUTA!

# 4.A. Preparación de Datos INEGI
ageb_sf <- st_read(ruta_ageb_shp) %>%
  st_transform(crs = 4326)

censo_datos <- read_csv(ruta_censo_csv, show_col_types = FALSE)

censo_datos_vars <- censo_datos %>%
  mutate(
    CVE_LOC = paste0(
      str_pad(MUN, 3, pad = "0"),
      str_pad(LOC, 4, pad = "0")
    )
  ) %>%
  select(CVE_LOC, POBTOT, TVIVHAB, GRAPROES) %>%
  rename(POB_TOTAL = POBTOT) %>%
  mutate(Densidad_Pob_LOC = POB_TOTAL / TVIVHAB) %>%
  select(CVE_LOC, Densidad_Pob_LOC, GRAPROES)


# 4.B. Unión de Datos (Spatial Join)
ageb_con_datos <- ageb_sf %>%
  mutate(
    CVE_LOC = paste0(
      str_pad(CVE_MUN, 3, pad = "0"),
      str_pad(CVE_LOC, 4, pad = "0")
    ),
    CVE_LOC = as.character(CVE_LOC)
  ) %>%
  left_join(censo_datos_vars, by = "CVE_LOC") %>%
  filter(!is.na(Densidad_Pob_LOC))

eggs_cov_inegi <- st_join(
  eggs_sf_limpio,
  ageb_con_datos %>% select(Densidad_Pob_LOC, GRAPROES),
  join = st_intersects
)

eggs_data_semi_completa <- eggs_cov_inegi %>%
  st_drop_geometry() %>%
  as_tibble()

# 4.C. Exportar Resultado Intermedio
write_csv(eggs_data_semi_completa, "data_covariables/eggs_data_semi_completa_localidad.csv")

## 5. Covariables Climáticas Dinámicas (WorldClim - Mensuales)

# Cargar la base de datos con covariables sociodemográficas
eggs_data_final_estatica <- read_csv("data_covariables/eggs_data_semi_completa_localidad.csv", show_col_types = FALSE)
eggs_sf_clima <- st_as_sf(eggs_data_final_estatica, coords = c("lon", "lat"), crs = 4326)

# 5.A. Descarga y Extracción de Datos Mensuales

cat("Descargando datos mensuales de Precipitación (Prec) y Temperatura Máx (Tmax)...\n")
eggs_data_final_clima <- eggs_data_final_estatica # Inicializar

# 1. DESCARGA Y EXTRACCIÓN DE TEMPERATURA MÁXIMA (TMAX)
tmax_rast <- geodata::worldclim_global(var = "tmax", res = 0.5, path = "data_covariables/")
tmax_yucatan <- crop(tmax_rast, yuc_bounds)
tmax_extract <- extract(tmax_yucatan, eggs_sf_clima, ID = FALSE)
names(tmax_extract) <- paste0("tmax_", str_pad(1:12, 2, pad = "0"))
eggs_data_final_clima <- cbind(eggs_data_final_clima, tmax_extract)

# 2. DESCARGA Y EXTRACCIÓN DE PRECIPITACIÓN (PREC)
prcp_rast <- geodata::worldclim_global(var = "prec", res = 0.5, path = "data_covariables/")
prcp_yucatan <- crop(prcp_rast, yuc_bounds)
prcp_extract <- extract(prcp_yucatan, eggs_sf_clima, ID = FALSE)
names(prcp_extract) <- paste0("prcp_", str_pad(1:12, 2, pad = "0"))
eggs_data_final_clima <- cbind(eggs_data_final_clima, prcp_extract)

# 5.B. Resultado Final
write_csv(eggs_data_final_clima, "data_covariables/eggs_data_base_COMPLETA_FINAL_MENSUAL.csv")
