
## ============================================================
## 07_UAV – SCRIPT ÚNICO
## ORTOMOSAICOS, CORRECCIÓN DE DESFASE, GLI, GRVI Y EXTRACCIÓN
## ============================================================
##
## Carpeta de trabajo:
##   D:/UDCA_2017_CLASE/46_ARTICULOS_2025/
##   ARTICULO_CEBADA_DOCTORADO/AJUSTES_ARTÍCULO/07_UAV
##
## Entradas esperadas:
##   - Dos ortomosaicos .tif o .tiff:
##       * uno correspondiente al 09 de noviembre;
##       * uno correspondiente al 14 de noviembre.
##
## Entradas auxiliares:
##   - Se intenta recuperar automáticamente:
##       a) los 56 nodos de cultivo;
##       b) el polígono del lote;
##     desde los objetos ya generados en las secciones 04 o 06.
##
## Salidas:
##   07_UAV/
##   ├── 01_AUDITORIA
##   ├── 02_ALINEACION
##   ├── 03_INDICES
##   ├── 04_CAMBIOS
##   ├── 05_EXTRACCION_NODOS
##   ├── 06_FIGURAS
##   ├── 07_TABLAS
##   └── 08_OBJETOS_R
##
## ÍNDICES:
##
##   GLI = (2G - R - B) / (2G + R + B)
##
##   GRVI = (G - R) / (G + R)
##
## PROPÓSITO:
##   1. Auditar CRS, resolución, extensión y bandas.
##   2. Identificar automáticamente los ortomosaicos del 09 y 14.
##   3. Reproyectar ambos a una malla común.
##   4. Estimar y corregir un posible desplazamiento horizontal
##      del ortomosaico del 09 frente al del 14.
##   5. Calcular GLI y GRVI para ambas fechas.
##   6. Calcular cambios 14 – 09.
##   7. Extraer valores puntuales y promedios en buffers alrededor
##      de los 56 nodos de cultivo.
##   8. Exportar TIFF, CSV, XLSX, PNG y objetos R reproducibles.
##
## ADVERTENCIAS:
##   - El ajuste automático corrige únicamente traslación X–Y.
##   - No corrige rotación, deformación local ni errores de relieve.
##   - La alineación debe verificarse visualmente antes de interpretar
##     diferencias temporales.
##   - El orden de bandas se asume R = 1, G = 2, B = 3.
##     Modifique los parámetros si los TIFF usan otro orden.
## ============================================================


## 0. Paquetes -------------------------------------------------

paquetes <- c(
  "terra",
  "sf",
  "dplyr",
  "tidyr",
  "purrr",
  "stringr",
  "readr",
  "writexl",
  "ggplot2",
  "patchwork",
  "tibble"
)

instalar <- paquetes[
  !paquetes %in% rownames(installed.packages())
]

if (length(instalar) > 0) {
  install.packages(instalar)
}

library(terra)
library(sf)
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(readr)
library(writexl)
library(ggplot2)
library(patchwork)
library(tibble)


## 1. Parámetros editables -------------------------------------

ruta_uav <- paste0(
  "D:/UDCA_2017_CLASE/46_ARTICULOS_2025/",
  "ARTICULO_CEBADA_DOCTORADO/AJUSTES_ARTÍCULO/07_UAV"
)

# Orden de bandas del ortomosaico.
banda_R <- 1
banda_G <- 2
banda_B <- 3

# Fecha que se utilizará como referencia geométrica.
# I2 y la geometría del lote se trabajaron sobre el vuelo del 14.
fecha_referencia <- "14"

# Búsqueda de desplazamiento máximo en metros.
# Ajuste este valor si el desfase visual supera este rango.
desplazamiento_max_m <- 5

# Resolución usada en la búsqueda automática de desplazamiento.
# Se agrega la imagen para acelerar el cálculo.
resolucion_busqueda_m <- 0.25

# Paso de búsqueda del desplazamiento en metros.
paso_desplazamiento_m <- 0.25

# Número máximo de celdas válidas usadas para evaluar correlación.
muestra_max_correlacion <- 150000

# Buffers para extracción alrededor de cada nodo.
buffers_m <- c(
  0,
  0.50,
  1.00
)

# CRS métrico esperado del proyecto.
crs_objetivo <- "EPSG:32618"

# Método de remuestreo para imágenes continuas.
metodo_remuestreo <- "bilinear"

# Rango permitido para índices.
rango_min_indice <- -1
rango_max_indice <- 1

# Presentación cartográfica de GLI y GRVI.
# TRUE: clases o bandas comunes para ambas fechas.
# FALSE: escala continua independiente por índice.
usar_bandas_clasificadas <- TRUE

# Bandas fijas de interpretación gráfica para GLI y GRVI.
cortes_indices <- c(
  -1.00,
  -0.50,
  -0.25,
   0.00,
   0.25,
   0.50,
   0.75,
   1.00
)

etiquetas_bandas_indices <- c(
  "-1.00 a -0.50",
  "-0.50 a -0.25",
  "-0.25 a 0.00",
  "0.00 a 0.25",
  "0.25 a 0.50",
  "0.50 a 0.75",
  "0.75 a 1.00"
)

# Semilla reproducible.
set.seed(2026)


## 2. Estructura de salida -------------------------------------

ruta_auditoria <- file.path(
  ruta_uav,
  "01_AUDITORIA"
)

ruta_alineacion <- file.path(
  ruta_uav,
  "02_ALINEACION"
)

ruta_indices <- file.path(
  ruta_uav,
  "03_INDICES"
)

ruta_cambios <- file.path(
  ruta_uav,
  "04_CAMBIOS"
)

ruta_extraccion <- file.path(
  ruta_uav,
  "05_EXTRACCION_NODOS"
)

ruta_figuras <- file.path(
  ruta_uav,
  "06_FIGURAS"
)

ruta_tablas <- file.path(
  ruta_uav,
  "07_TABLAS"
)

ruta_objetos <- file.path(
  ruta_uav,
  "08_OBJETOS_R"
)

rutas_salida <- c(
  ruta_auditoria,
  ruta_alineacion,
  ruta_indices,
  ruta_cambios,
  ruta_extraccion,
  ruta_figuras,
  ruta_tablas,
  ruta_objetos
)

walk(
  rutas_salida,
  function(ruta) {
    dir.create(
      ruta,
      recursive = TRUE,
      showWarnings = FALSE
    )
  }
)


## 3. Localizar ortomosaicos -----------------------------------

archivos_tif <- list.files(
  ruta_uav,
  pattern = "\\.(tif|tiff)$",
  full.names = TRUE,
  ignore.case = TRUE,
  recursive = FALSE
)

if (length(archivos_tif) < 2) {
  stop(
    "Se requieren al menos dos archivos TIFF en:\n",
    ruta_uav
  )
}

identificar_fecha <- function(ruta) {

  nombre <- basename(
    ruta
  ) |>
    str_to_lower()

  case_when(
    str_detect(
      nombre,
      "(^|[^0-9])0?9([^0-9]|$)|09[_-]?11|nov[^a-z]*0?9"
    ) ~
      "09",

    str_detect(
      nombre,
      "(^|[^0-9])14([^0-9]|$)|14[_-]?11|nov[^a-z]*14"
    ) ~
      "14",

    TRUE ~
      NA_character_
  )
}

tabla_archivos <- tibble(
  Ruta = archivos_tif,
  Archivo = basename(
    archivos_tif
  ),
  Fecha_detectada = map_chr(
    archivos_tif,
    identificar_fecha
  )
)

archivo_09 <- tabla_archivos |>
  filter(
    Fecha_detectada == "09"
  ) |>
  slice_head(
    n = 1
  ) |>
  pull(
    Ruta
  )

archivo_14 <- tabla_archivos |>
  filter(
    Fecha_detectada == "14"
  ) |>
  slice_head(
    n = 1
  ) |>
  pull(
    Ruta
  )

if (
  length(archivo_09) == 0 ||
    length(archivo_14) == 0
) {

  stop(
    paste0(
      "No fue posible identificar automáticamente los TIFF del 09 ",
      "y 14 de noviembre.\n\nArchivos encontrados:\n",
      paste(
        basename(
          archivos_tif
        ),
        collapse = "\n"
      ),
      "\n\nRenombre los archivos incluyendo 09 y 14."
    )
  )
}


## 4. Cargar ortomosaicos --------------------------------------

orto_09_original <- rast(
  archivo_09
)

orto_14_original <- rast(
  archivo_14
)

if (
  nlyr(orto_09_original) <
    max(
      c(
        banda_R,
        banda_G,
        banda_B
      )
    )
) {
  stop(
    "El ortomosaico del 09 no contiene las bandas RGB indicadas."
  )
}

if (
  nlyr(orto_14_original) <
    max(
      c(
        banda_R,
        banda_G,
        banda_B
      )
    )
) {
  stop(
    "El ortomosaico del 14 no contiene las bandas RGB indicadas."
  )
}


## 5. Auditoría básica -----------------------------------------

auditar_raster <- function(
  raster_obj,
  fecha,
  ruta_archivo
) {

  tibble(
    Fecha = fecha,
    Archivo = basename(
      ruta_archivo
    ),
    Capas = nlyr(
      raster_obj
    ),
    Filas = nrow(
      raster_obj
    ),
    Columnas = ncol(
      raster_obj
    ),
    Resolucion_X = res(
      raster_obj
    )[1],
    Resolucion_Y = res(
      raster_obj
    )[2],
    X_min = xmin(
      raster_obj
    ),
    X_max = xmax(
      raster_obj
    ),
    Y_min = ymin(
      raster_obj
    ),
    Y_max = ymax(
      raster_obj
    ),
    CRS = crs(
      raster_obj,
      proj = TRUE
    ),
    Tiene_CRS = crs(
      raster_obj
    ) != ""
  )
}

tabla_auditoria_raster <- bind_rows(
  auditar_raster(
    orto_09_original,
    "09",
    archivo_09
  ),
  auditar_raster(
    orto_14_original,
    "14",
    archivo_14
  )
)

if (
  any(
    !tabla_auditoria_raster$Tiene_CRS
  )
) {
  stop(
    "Uno o ambos ortomosaicos no tienen CRS definido."
  )
}


## 6. Recuperar nodos y lote -----------------------------------

ruta_obj_06A <- paste0(
  "D:/UDCA_2017_CLASE/46_ARTICULOS_2025/",
  "ARTICULO_CEBADA_DOCTORADO/AJUSTES_ARTÍCULO/",
  "06_SUELO_CULTIVO/01_BASE_INTEGRADA_AUDITADA/",
  "objetos_R/objetos_base_integrada_auditada_06A.rds"
)

ruta_obj_04 <- paste0(
  "D:/UDCA_2017_CLASE/46_ARTICULOS_2025/",
  "ARTICULO_CEBADA_DOCTORADO/AJUSTES_ARTÍCULO/",
  "04_TRATAMIENTOS/02_INTEGRACION_SUELO_CULTIVO/",
  "objetos_R/objetos_integracion_suelo_cultivo_tratamientos.rds"
)

obj_auxiliar <- NULL
fuente_objeto <- NA_character_

if (file.exists(ruta_obj_06A)) {

  obj_auxiliar <- readRDS(
    ruta_obj_06A
  )

  fuente_objeto <- ruta_obj_06A

} else if (file.exists(ruta_obj_04)) {

  obj_auxiliar <- readRDS(
    ruta_obj_04
  )

  fuente_objeto <- ruta_obj_04
}

nodos_sf <- NULL
lote_sf <- NULL
base_nodos <- NULL

if (!is.null(obj_auxiliar)) {

  if (!is.null(obj_auxiliar$base_auditada)) {

    base_nodos <- as_tibble(
      obj_auxiliar$base_auditada
    )

  } else if (!is.null(obj_auxiliar$base_integrada)) {

    base_nodos <- as_tibble(
      obj_auxiliar$base_integrada
    )
  }

  if (!is.null(obj_auxiliar$cultivo_sf)) {

    nodos_sf <- st_as_sf(
      obj_auxiliar$cultivo_sf
    )
  }

  if (!is.null(obj_auxiliar$lote)) {

    lote_sf <- st_as_sf(
      obj_auxiliar$lote
    )
  }
}

if (
  is.null(nodos_sf) &&
    !is.null(base_nodos) &&
    all(
      c(
        "X",
        "Y"
      ) %in% names(
        base_nodos
      )
    )
) {

  nodos_sf <- st_as_sf(
    base_nodos,
    coords = c(
      "X",
      "Y"
    ),
    crs = crs_objetivo,
    remove = FALSE
  )
}

if (!is.null(nodos_sf)) {

  if (is.na(st_crs(nodos_sf))) {
    st_crs(nodos_sf) <- crs_objetivo
  }

  nodos_sf <- st_transform(
    nodos_sf,
    crs_objetivo
  )
}

if (!is.null(lote_sf)) {

  if (is.na(st_crs(lote_sf))) {
    st_crs(lote_sf) <- crs_objetivo
  }

  lote_sf <- st_transform(
    lote_sf,
    crs_objetivo
  ) |>
    st_make_valid()
}


## 7. Reproyección a CRS común ---------------------------------

orto_09_utm <- if (
  crs(
    orto_09_original,
    proj = TRUE
  ) != crs_objetivo
) {
  project(
    orto_09_original,
    crs_objetivo,
    method = metodo_remuestreo
  )
} else {
  orto_09_original
}

orto_14_utm <- if (
  crs(
    orto_14_original,
    proj = TRUE
  ) != crs_objetivo
) {
  project(
    orto_14_original,
    crs_objetivo,
    method = metodo_remuestreo
  )
} else {
  orto_14_original
}


## 8. Seleccionar bandas RGB -----------------------------------

seleccionar_rgb <- function(
  orto
) {

  rgb <- orto[[c(
    banda_R,
    banda_G,
    banda_B
  )]]

  names(
    rgb
  ) <- c(
    "R",
    "G",
    "B"
  )

  rgb
}

rgb_09 <- seleccionar_rgb(
  orto_09_utm
)

rgb_14 <- seleccionar_rgb(
  orto_14_utm
)


## 9. Normalizar radiometría para alineación -------------------

normalizar_banda <- function(x) {

  cuantiles <- global(
    x,
    fun = function(v) {
      quantile(
        v,
        probs = c(
          0.02,
          0.98
        ),
        na.rm = TRUE,
        names = FALSE
      )
    }
  )

  q02 <- as.numeric(
    cuantiles[1, 1]
  )

  q98 <- as.numeric(
    cuantiles[1, 2]
  )

  salida <- (
    x -
      q02
  ) /
    (
      q98 -
        q02
    )

  clamp(
    salida,
    lower = 0,
    upper = 1,
    values = TRUE
  )
}

crear_luminancia <- function(rgb) {

  r <- normalizar_banda(
    rgb[["R"]]
  )

  g <- normalizar_banda(
    rgb[["G"]]
  )

  b <- normalizar_banda(
    rgb[["B"]]
  )

  0.299 *
    r +
    0.587 *
      g +
    0.114 *
      b
}

lum_09 <- crear_luminancia(
  rgb_09
)

lum_14 <- crear_luminancia(
  rgb_14
)


## 10. Preparar referencia común -------------------------------

# Se toma como referencia geométrica el ortomosaico del 14.
# El ortomosaico del 09 será desplazado y remuestreado sobre esta malla.
referencia_geom <- lum_14

# Agregar a resolución de búsqueda.
factor_agregacion <- max(
  1,
  round(
    resolucion_busqueda_m /
      mean(
        res(
          referencia_geom
        )
      )
  )
)

lum_14_busqueda <- aggregate(
  lum_14,
  fact = factor_agregacion,
  fun = mean,
  na.rm = TRUE
)

lum_09_busqueda_base <- resample(
  lum_09,
  lum_14_busqueda,
  method = metodo_remuestreo
)


## 11. Máscara espacial de búsqueda ----------------------------

mascara_busqueda <- !is.na(
  lum_14_busqueda
) &
  !is.na(
    lum_09_busqueda_base
)

if (!is.null(lote_sf)) {

  lote_vect <- vect(
    lote_sf
  )

  mascara_lote <- rasterize(
    lote_vect,
    lum_14_busqueda,
    field = 1,
    background = NA
  )

  mascara_busqueda <- mascara_busqueda &
    !is.na(
      mascara_lote
    )
}


## 12. Función de evaluación de desplazamiento -----------------

evaluar_desplazamiento <- function(
  dx_m,
  dy_m
) {

  # Se desplaza el ortomosaico móvil del 09 sobre la referencia del 14.
  desplazado <- shift(
    lum_09_busqueda_base,
    dx = dx_m,
    dy = dy_m
  )

  desplazado <- resample(
    desplazado,
    lum_14_busqueda,
    method = metodo_remuestreo
  )

  valores <- c(
    lum_14_busqueda,
    desplazado,
    mascara_busqueda
  ) |>
    as.data.frame(
      na.rm = FALSE
    )

  names(
    valores
  ) <- c(
    "Referencia",
    "Movil",
    "Mascara"
  )

  valores <- valores |>
    filter(
      !is.na(
        Mascara
      ),
      is.finite(
        Referencia
      ),
      is.finite(
        Movil
      )
    )

  if (nrow(valores) < 100) {

    return(
      tibble(
        DX_m = dx_m,
        DY_m = dy_m,
        N_celdas = nrow(
          valores
        ),
        Correlacion = NA_real_,
        RMSE = NA_real_
      )
    )
  }

  if (nrow(valores) > muestra_max_correlacion) {

    valores <- valores |>
      slice_sample(
        n = muestra_max_correlacion
      )
  }

  tibble(
    DX_m = dx_m,
    DY_m = dy_m,
    N_celdas = nrow(
      valores
    ),
    Correlacion = suppressWarnings(
      cor(
        valores$Referencia,
        valores$Movil,
        method = "pearson",
        use = "complete.obs"
      )
    ),
    RMSE = sqrt(
      mean(
        (
          valores$Referencia -
            valores$Movil
        )^2,
        na.rm = TRUE
      )
    )
  )
}


## 13. Búsqueda automática del desplazamiento ------------------

secuencia_desplazamiento <- seq(
  -desplazamiento_max_m,
  desplazamiento_max_m,
  by = paso_desplazamiento_m
)

grilla_desplazamientos <- tidyr::crossing(
  DX_m = secuencia_desplazamiento,
  DY_m = secuencia_desplazamiento
)

message(
  "Evaluando ",
  nrow(
    grilla_desplazamientos
  ),
  " combinaciones de desplazamiento..."
)

tabla_desplazamientos <- pmap_dfr(
  grilla_desplazamientos,
  function(
    DX_m,
    DY_m
  ) {
    evaluar_desplazamiento(
      dx_m = DX_m,
      dy_m = DY_m
    )
  }
)

mejor_desplazamiento <- tabla_desplazamientos |>
  filter(
    is.finite(
      Correlacion
    ),
    is.finite(
      RMSE
    )
  ) |>
  arrange(
    desc(
      Correlacion
    ),
    RMSE,
    abs(
      DX_m
    ) +
      abs(
        DY_m
      )
  ) |>
  slice_head(
    n = 1
  )

if (nrow(mejor_desplazamiento) == 0) {
  stop(
    "No fue posible estimar un desplazamiento válido."
  )
}

dx_optimo <- mejor_desplazamiento$DX_m[1]
dy_optimo <- mejor_desplazamiento$DY_m[1]


## 14. Aplicar desplazamiento al ortomosaico del 09 ------------

rgb_09_desplazado <- shift(
  rgb_09,
  dx = dx_optimo,
  dy = dy_optimo
)

# El 14 conserva su geometría y malla original.
rgb_14_alineado <- rgb_14

# El 09 se remuestrea sobre la malla de referencia del 14.
rgb_09_alineado <- resample(
  rgb_09_desplazado,
  rgb_14_alineado,
  method = metodo_remuestreo
)

# Recortar a extensión común efectiva.
extension_comun <- intersect(
  ext(
    rgb_14_alineado
  ),
  ext(
    rgb_09_alineado
  )
)

rgb_14_alineado <- crop(
  rgb_14_alineado,
  extension_comun
)

rgb_09_alineado <- crop(
  rgb_09_alineado,
  extension_comun
)

# Enmascarar al lote cuando esté disponible.
if (!is.null(lote_sf)) {

  lote_vect <- vect(
    lote_sf
  )

  rgb_14_alineado <- mask(
    crop(
      rgb_14_alineado,
      lote_vect
    ),
    lote_vect
  )

  rgb_09_alineado <- mask(
    crop(
      rgb_09_alineado,
      lote_vect
    ),
    lote_vect
  )
}


## 15. Guardar RGB alineados -----------------------------------

writeRaster(
  rgb_09_alineado,
  file.path(
    ruta_alineacion,
    "Ortomosaico_RGB_09_alineado.tif"
  ),
  overwrite = TRUE,
  gdal = c(
    "COMPRESS=LZW",
    "TILED=YES"
  )
)

writeRaster(
  rgb_14_alineado,
  file.path(
    ruta_alineacion,
    "Ortomosaico_RGB_14_alineado.tif"
  ),
  overwrite = TRUE,
  gdal = c(
    "COMPRESS=LZW",
    "TILED=YES"
  )
)


## 16. Calcular GLI y GRVI -------------------------------------

calcular_indices <- function(
  rgb,
  fecha
) {

  r <- as.numeric(
    1
  )

  R <- rgb[["R"]]
  G <- rgb[["G"]]
  B <- rgb[["B"]]

  denominador_gli <- 2 *
    G +
    R +
    B

  denominador_grvi <- G +
    R

  GLI <- ifel(
    abs(
      denominador_gli
    ) >
      .Machine$double.eps,
    (
      2 *
        G -
        R -
        B
    ) /
      denominador_gli,
    NA
  )

  GRVI <- ifel(
    abs(
      denominador_grvi
    ) >
      .Machine$double.eps,
    (
      G -
        R
    ) /
      denominador_grvi,
    NA
  )

  GLI <- clamp(
    GLI,
    lower = rango_min_indice,
    upper = rango_max_indice,
    values = TRUE
  )

  GRVI <- clamp(
    GRVI,
    lower = rango_min_indice,
    upper = rango_max_indice,
    values = TRUE
  )

  names(
    GLI
  ) <- paste0(
    "GLI_",
    fecha
  )

  names(
    GRVI
  ) <- paste0(
    "GRVI_",
    fecha
  )

  c(
    GLI,
    GRVI
  )
}

indices_09 <- calcular_indices(
  rgb_09_alineado,
  "09"
)

indices_14 <- calcular_indices(
  rgb_14_alineado,
  "14"
)

writeRaster(
  indices_09[["GLI_09"]],
  file.path(
    ruta_indices,
    "GLI_09_noviembre.tif"
  ),
  overwrite = TRUE,
  gdal = c(
    "COMPRESS=LZW",
    "TILED=YES"
  )
)

writeRaster(
  indices_09[["GRVI_09"]],
  file.path(
    ruta_indices,
    "GRVI_09_noviembre.tif"
  ),
  overwrite = TRUE,
  gdal = c(
    "COMPRESS=LZW",
    "TILED=YES"
  )
)

writeRaster(
  indices_14[["GLI_14"]],
  file.path(
    ruta_indices,
    "GLI_14_noviembre.tif"
  ),
  overwrite = TRUE,
  gdal = c(
    "COMPRESS=LZW",
    "TILED=YES"
  )
)

writeRaster(
  indices_14[["GRVI_14"]],
  file.path(
    ruta_indices,
    "GRVI_14_noviembre.tif"
  ),
  overwrite = TRUE,
  gdal = c(
    "COMPRESS=LZW",
    "TILED=YES"
  )
)


## 17. Cambios temporales 14 – 09 ------------------------------

delta_gli <- indices_14[["GLI_14"]] -
  indices_09[["GLI_09"]]

delta_grvi <- indices_14[["GRVI_14"]] -
  indices_09[["GRVI_09"]]

names(
  delta_gli
) <- "Delta_GLI_14_menos_09"

names(
  delta_grvi
) <- "Delta_GRVI_14_menos_09"

writeRaster(
  delta_gli,
  file.path(
    ruta_cambios,
    "Delta_GLI_14_menos_09.tif"
  ),
  overwrite = TRUE,
  gdal = c(
    "COMPRESS=LZW",
    "TILED=YES"
  )
)

writeRaster(
  delta_grvi,
  file.path(
    ruta_cambios,
    "Delta_GRVI_14_menos_09.tif"
  ),
  overwrite = TRUE,
  gdal = c(
    "COMPRESS=LZW",
    "TILED=YES"
  )
)


## 18. Auditoría estadística de índices ------------------------

resumir_raster <- function(
  raster_obj,
  variable
) {

  valores <- values(
    raster_obj,
    mat = FALSE
  )

  valores <- valores[
    is.finite(
      valores
    )
  ]

  if (length(valores) == 0) {

    return(
      tibble(
        Variable = variable,
        n_celdas = 0,
        Minimo = NA_real_,
        Q05 = NA_real_,
        Mediana = NA_real_,
        Media = NA_real_,
        Q95 = NA_real_,
        Maximo = NA_real_,
        DE = NA_real_
      )
    )
  }

  tibble(
    Variable = variable,
    n_celdas = length(
      valores
    ),
    Minimo = min(
      valores
    ),
    Q05 = quantile(
      valores,
      0.05,
      names = FALSE
    ),
    Mediana = median(
      valores
    ),
    Media = mean(
      valores
    ),
    Q95 = quantile(
      valores,
      0.95,
      names = FALSE
    ),
    Maximo = max(
      valores
    ),
    DE = sd(
      valores
    )
  )
}

tabla_resumen_indices <- bind_rows(
  resumir_raster(
    indices_09[["GLI_09"]],
    "GLI_09"
  ),
  resumir_raster(
    indices_09[["GRVI_09"]],
    "GRVI_09"
  ),
  resumir_raster(
    indices_14[["GLI_14"]],
    "GLI_14"
  ),
  resumir_raster(
    indices_14[["GRVI_14"]],
    "GRVI_14"
  ),
  resumir_raster(
    delta_gli,
    "Delta_GLI_14_menos_09"
  ),
  resumir_raster(
    delta_grvi,
    "Delta_GRVI_14_menos_09"
  )
)


## 19. Extracción en nodos -------------------------------------

extraer_en_nodos <- function(
  raster_obj,
  nodos,
  buffer_m,
  nombre_variable
) {

  if (is.null(nodos)) {
    return(
      tibble()
    )
  }

  nodos_vect <- vect(
    nodos
  )

  if (buffer_m == 0) {

    extraido <- terra::extract(
      raster_obj,
      nodos_vect,
      method = "bilinear"
    )

  } else {

    extraido <- terra::extract(
      raster_obj,
      nodos_vect,
      buffer = buffer_m,
      fun = mean,
      na.rm = TRUE
    )
  }

  valor_columna <- setdiff(
    names(
      extraido
    ),
    "ID"
  )[1]

  nodo_id <- if (
    "Nodo" %in% names(
      nodos
    )
  ) {
    nodos$Nodo
  } else {
    seq_len(
      nrow(
        nodos
      )
    )
  }

  tibble(
    Nodo = nodo_id,
    Variable = nombre_variable,
    Buffer_m = buffer_m,
    Valor = extraido[[valor_columna]]
  )
}

lista_rasters_extraccion <- list(
  GLI_09 = indices_09[["GLI_09"]],
  GRVI_09 = indices_09[["GRVI_09"]],
  GLI_14 = indices_14[["GLI_14"]],
  GRVI_14 = indices_14[["GRVI_14"]],
  Delta_GLI = delta_gli,
  Delta_GRVI = delta_grvi
)

tabla_extraccion_larga <- tibble()

if (!is.null(nodos_sf)) {

  resultados_extraccion <- list()
  contador <- 1

  for (nombre_raster in names(lista_rasters_extraccion)) {

    for (buffer_actual in buffers_m) {

      resultados_extraccion[[contador]] <- extraer_en_nodos(
        raster_obj = lista_rasters_extraccion[[nombre_raster]],
        nodos = nodos_sf,
        buffer_m = buffer_actual,
        nombre_variable = nombre_raster
      )

      contador <- contador + 1
    }
  }

  tabla_extraccion_larga <- bind_rows(
    resultados_extraccion
  )
}


## 20. Base ancha por nodo -------------------------------------

tabla_extraccion_ancha <- tibble()

if (nrow(tabla_extraccion_larga) > 0) {

  tabla_extraccion_ancha <- tabla_extraccion_larga |>
    mutate(
      Nombre_columna = paste0(
        Variable,
        "_buffer_",
        str_replace_all(
          format(
            Buffer_m,
            trim = TRUE,
            scientific = FALSE
          ),
          "\\.",
          "_"
        ),
        "m"
      )
    ) |>
    select(
      Nodo,
      Nombre_columna,
      Valor
    ) |>
    pivot_wider(
      names_from = Nombre_columna,
      values_from = Valor
    )

  if (!is.null(base_nodos)) {

    tabla_extraccion_ancha <- base_nodos |>
      select(
        any_of(
          c(
            "Nodo",
            "Experimental_Unit",
            "Treatment_ID",
            "Treatment_base",
            "Dose_percent",
            "Omission_description",
            "Replication",
            "Macollas",
            "Altura",
            "Diametro",
            "NDVI",
            "Nitratos",
            "Potasio_savia",
            "Calcio_savia",
            "X",
            "Y"
          )
        )
      ) |>
      left_join(
        tabla_extraccion_ancha,
        by = "Nodo"
      )
  }
}


## 21. Correlaciones exploratorias UAV–cultivo -----------------

tabla_correlaciones_uav_cultivo <- tibble()

if (
  nrow(tabla_extraccion_ancha) > 0 &&
    !is.null(base_nodos)
) {

  variables_uav <- names(
    tabla_extraccion_ancha
  )[
    str_detect(
      names(
        tabla_extraccion_ancha
      ),
      "^(GLI|GRVI|Delta_)"
    )
  ]

  variables_cultivo <- intersect(
    c(
      "Macollas",
      "Altura",
      "Diametro",
      "NDVI",
      "Nitratos",
      "Potasio_savia",
      "Calcio_savia"
    ),
    names(
      tabla_extraccion_ancha
    )
  )

  pares <- tidyr::crossing(
    Variable_UAV = variables_uav,
    Variable_cultivo = variables_cultivo
  )

  tabla_correlaciones_uav_cultivo <- pmap_dfr(
    pares,
    function(
      Variable_UAV,
      Variable_cultivo
    ) {

      x <- tabla_extraccion_ancha[[Variable_UAV]]
      y <- tabla_extraccion_ancha[[Variable_cultivo]]

      validos <- is.finite(
        x
      ) &
        is.finite(
          y
        )

      if (sum(validos) < 4) {

        return(
          tibble(
            Variable_UAV = Variable_UAV,
            Variable_cultivo = Variable_cultivo,
            n = sum(
              validos
            ),
            Pearson = NA_real_,
            Spearman = NA_real_
          )
        )
      }

      tibble(
        Variable_UAV = Variable_UAV,
        Variable_cultivo = Variable_cultivo,
        n = sum(
          validos
        ),
        Pearson = suppressWarnings(
          cor(
            x[validos],
            y[validos],
            method = "pearson"
          )
        ),
        Spearman = suppressWarnings(
          cor(
            x[validos],
            y[validos],
            method = "spearman"
          )
        )
      )
    }
  )
}


## 22. Figuras de alineación -----------------------------------

tabla_superficie_desplazamiento <- tabla_desplazamientos |>
  filter(
    is.finite(
      Correlacion
    )
  )

g_desplazamiento <- ggplot(
  tabla_superficie_desplazamiento,
  aes(
    x = DX_m,
    y = DY_m,
    fill = Correlacion
  )
) +
  geom_raster() +
  geom_point(
    data = mejor_desplazamiento,
    aes(
      x = DX_m,
      y = DY_m
    ),
    inherit.aes = FALSE,
    shape = 4,
    size = 4,
    stroke = 1.2
  ) +
  coord_equal() +
  labs(
    title = "Búsqueda automática de desplazamiento",
    subtitle = paste0(
      "Óptimo: DX = ",
      round(
        dx_optimo,
        2
      ),
      " m; DY = ",
      round(
        dy_optimo,
        2
      ),
      " m"
    ),
    x = "Desplazamiento X (m)",
    y = "Desplazamiento Y (m)",
    fill = "Correlación"
  ) +
  theme_minimal()

ggsave(
  filename = file.path(
    ruta_figuras,
    "01_superficie_busqueda_desplazamiento.png"
  ),
  plot = g_desplazamiento,
  width = 9,
  height = 7,
  dpi = 400
)


## 23. Figuras de índices --------------------------------------

raster_a_df <- function(
  raster_obj,
  nombre,
  indice
) {

  df <- as.data.frame(
    raster_obj,
    xy = TRUE,
    na.rm = TRUE
  )

  names(
    df
  )[3] <- "Valor"

  df |>
    mutate(
      Variable = nombre,
      Indice = indice
    )
}

df_gli <- bind_rows(
  raster_a_df(
    indices_09[["GLI_09"]],
    "09 Nov",
    "GLI"
  ),
  raster_a_df(
    indices_14[["GLI_14"]],
    "14 Nov",
    "GLI"
  )
)

df_grvi <- bind_rows(
  raster_a_df(
    indices_09[["GRVI_09"]],
    "09 Nov",
    "GRVI"
  ),
  raster_a_df(
    indices_14[["GRVI_14"]],
    "14 Nov",
    "GRVI"
  )
)

clasificar_bandas <- function(df) {

  df |>
    mutate(
      Banda = cut(
        Valor,
        breaks = cortes_indices,
        labels = etiquetas_bandas_indices,
        include.lowest = TRUE,
        right = TRUE,
        ordered_result = TRUE
      )
    )
}

if (usar_bandas_clasificadas) {

  df_gli_plot <- clasificar_bandas(
    df_gli
  )

  df_grvi_plot <- clasificar_bandas(
    df_grvi
  )

  g_gli <- ggplot(
    df_gli_plot,
    aes(
      x = x,
      y = y,
      fill = Banda
    )
  ) +
    geom_raster() +
    facet_wrap(
      ~ Variable,
      ncol = 2
    ) +
    coord_equal() +
    labs(
      title = "GLI por fecha",
      subtitle = "Bandas comunes para facilitar la comparación temporal",
      x = "X (m)",
      y = "Y (m)",
      fill = "Banda GLI"
    ) +
    theme_minimal()

  g_grvi <- ggplot(
    df_grvi_plot,
    aes(
      x = x,
      y = y,
      fill = Banda
    )
  ) +
    geom_raster() +
    facet_wrap(
      ~ Variable,
      ncol = 2
    ) +
    coord_equal() +
    labs(
      title = "GRVI por fecha",
      subtitle = "Bandas comunes para facilitar la comparación temporal",
      x = "X (m)",
      y = "Y (m)",
      fill = "Banda GRVI"
    ) +
    theme_minimal()

} else {

  limites_gli <- quantile(
    df_gli$Valor,
    probs = c(
      0.02,
      0.98
    ),
    na.rm = TRUE,
    names = FALSE
  )

  limites_grvi <- quantile(
    df_grvi$Valor,
    probs = c(
      0.02,
      0.98
    ),
    na.rm = TRUE,
    names = FALSE
  )

  g_gli <- ggplot(
    df_gli,
    aes(
      x = x,
      y = y,
      fill = Valor
    )
  ) +
    geom_raster() +
    facet_wrap(
      ~ Variable,
      ncol = 2
    ) +
    coord_equal() +
    scale_fill_viridis_c(
      limits = limites_gli,
      oob = scales::squish
    ) +
    labs(
      title = "GLI por fecha",
      subtitle = "Escala continua común para ambas fechas",
      x = "X (m)",
      y = "Y (m)",
      fill = "GLI"
    ) +
    theme_minimal()

  g_grvi <- ggplot(
    df_grvi,
    aes(
      x = x,
      y = y,
      fill = Valor
    )
  ) +
    geom_raster() +
    facet_wrap(
      ~ Variable,
      ncol = 2
    ) +
    coord_equal() +
    scale_fill_viridis_c(
      limits = limites_grvi,
      oob = scales::squish
    ) +
    labs(
      title = "GRVI por fecha",
      subtitle = "Escala continua común para ambas fechas",
      x = "X (m)",
      y = "Y (m)",
      fill = "GRVI"
    ) +
    theme_minimal()
}

panel_indices <- g_gli /
  g_grvi +
  plot_annotation(
    title = paste0(
      "Índices UAV alineados sobre la referencia del 14 de noviembre"
    )
  )

ggsave(
  filename = file.path(
    ruta_figuras,
    "02_panel_GLI_GRVI_fechas.png"
  ),
  plot = panel_indices,
  width = 14,
  height = 13,
  dpi = 400
)


df_cambios <- bind_rows(
  raster_a_df(
    delta_gli,
    "Δ GLI: 14 – 09",
    "Delta_GLI"
  ),
  raster_a_df(
    delta_grvi,
    "Δ GRVI: 14 – 09",
    "Delta_GRVI"
  )
)

limite_delta <- quantile(
  abs(
    df_cambios$Valor
  ),
  probs = 0.98,
  na.rm = TRUE,
  names = FALSE
)

if (
  !is.finite(
    limite_delta
  ) ||
    limite_delta <= 0
) {
  limite_delta <- max(
    abs(
      df_cambios$Valor
    ),
    na.rm = TRUE
  )
}

g_cambios <- ggplot(
  df_cambios,
  aes(
    x = x,
    y = y,
    fill = Valor
  )
) +
  geom_raster() +
  facet_wrap(
    ~ Variable,
    ncol = 2
  ) +
  coord_equal() +
  scale_fill_gradient2(
    midpoint = 0,
    limits = c(
      -limite_delta,
      limite_delta
    ),
    oob = scales::squish
  ) +
  labs(
    title = "Cambio temporal de índices UAV",
    subtitle = paste0(
      "Valores positivos indican aumento entre el 09 y 14; ",
      "malla de referencia: 14 de noviembre"
    ),
    x = "X (m)",
    y = "Y (m)",
    fill = "Cambio"
  ) +
  theme_minimal()

ggsave(
  filename = file.path(
    ruta_figuras,
    "03_panel_cambios_GLI_GRVI.png"
  ),
  plot = g_cambios,
  width = 13,
  height = 6.5,
  dpi = 400
)


## 24. Figura de extracción en nodos ---------------------------

if (nrow(tabla_extraccion_larga) > 0) {

  nodos_coords <- nodos_sf |>
    mutate(
      X = st_coordinates(
        nodos_sf
      )[, 1],
      Y = st_coordinates(
        nodos_sf
      )[, 2]
    ) |>
    st_drop_geometry() |>
    select(
      Nodo,
      X,
      Y
    )

  df_nodos_indices <- tabla_extraccion_larga |>
    filter(
      Buffer_m == 0,
      Variable %in% c(
        "GLI_09",
        "GLI_14",
        "GRVI_09",
        "GRVI_14"
      )
    ) |>
    left_join(
      nodos_coords,
      by = "Nodo"
    )

  g_nodos <- ggplot(
    df_nodos_indices,
    aes(
      x = X,
      y = Y,
      fill = Valor
    )
  ) +
    geom_point(
      shape = 21,
      size = 3
    ) +
    facet_wrap(
      ~ Variable,
      ncol = 2
    ) +
    coord_equal() +
    labs(
      title = "GLI y GRVI extraídos en los nodos de cultivo",
      x = "X (m)",
      y = "Y (m)",
      fill = "Valor"
    ) +
    theme_minimal()

  ggsave(
    filename = file.path(
      ruta_figuras,
      "04_indices_extraidos_nodos.png"
    ),
    plot = g_nodos,
    width = 12,
    height = 9,
    dpi = 400
  )
}


## 25. Exportar tablas -----------------------------------------

tabla_control_alineacion <- tibble(
  Parametro = c(
    "Archivo referencia",
    "Archivo móvil",
    "DX óptimo (m)",
    "DY óptimo (m)",
    "Correlación óptima",
    "RMSE óptimo",
    "Rango máximo buscado (m)",
    "Paso de búsqueda (m)",
    "Resolución de búsqueda (m)",
    "CRS objetivo"
  ),
  Valor = c(
    basename(
      archivo_14
    ),
    basename(
      archivo_09
    ),
    dx_optimo,
    dy_optimo,
    mejor_desplazamiento$Correlacion[1],
    mejor_desplazamiento$RMSE[1],
    desplazamiento_max_m,
    paso_desplazamiento_m,
    resolucion_busqueda_m,
    crs_objetivo
  )
)

write_csv(
  tabla_auditoria_raster,
  file.path(
    ruta_auditoria,
    "Auditoria_ortomosaicos.csv"
  )
)

write_csv(
  tabla_desplazamientos,
  file.path(
    ruta_alineacion,
    "Evaluacion_desplazamientos.csv"
  )
)

write_csv(
  tabla_control_alineacion,
  file.path(
    ruta_alineacion,
    "Control_alineacion.csv"
  )
)

write_csv(
  tabla_resumen_indices,
  file.path(
    ruta_indices,
    "Resumen_estadistico_indices.csv"
  )
)

if (nrow(tabla_extraccion_larga) > 0) {

  write_csv(
    tabla_extraccion_larga,
    file.path(
      ruta_extraccion,
      "Extraccion_indices_nodos_formato_largo.csv"
    )
  )

  write_csv(
    tabla_extraccion_ancha,
    file.path(
      ruta_extraccion,
      "Base_nodos_indices_UAV.csv"
    )
  )
}

if (nrow(tabla_correlaciones_uav_cultivo) > 0) {

  write_csv(
    tabla_correlaciones_uav_cultivo,
    file.path(
      ruta_extraccion,
      "Correlaciones_exploratorias_UAV_cultivo.csv"
    )
  )
}


## 26. Libro Excel ---------------------------------------------

hojas_excel <- list(
  Auditoria_ortomosaicos =
    tabla_auditoria_raster,

  Control_alineacion =
    tabla_control_alineacion,

  Desplazamientos =
    tabla_desplazamientos,

  Resumen_indices =
    tabla_resumen_indices
)

if (nrow(tabla_extraccion_larga) > 0) {

  hojas_excel$Extraccion_larga <-
    tabla_extraccion_larga

  hojas_excel$Base_nodos_UAV <-
    tabla_extraccion_ancha
}

if (nrow(tabla_correlaciones_uav_cultivo) > 0) {

  hojas_excel$Correlaciones_UAV_cultivo <-
    tabla_correlaciones_uav_cultivo
}

write_xlsx(
  hojas_excel,
  path = file.path(
    ruta_tablas,
    "Resultados_UAV_GLI_GRVI.xlsx"
  )
)


## 27. Guardar objeto R ----------------------------------------

saveRDS(
  list(
    archivo_09 = archivo_09,
    archivo_14 = archivo_14,
    auditoria = tabla_auditoria_raster,
    desplazamientos = tabla_desplazamientos,
    mejor_desplazamiento = mejor_desplazamiento,
    rgb_09_alineado = rgb_09_alineado,
    rgb_14_alineado = rgb_14_alineado,
    indices_09 = indices_09,
    indices_14 = indices_14,
    delta_gli = delta_gli,
    delta_grvi = delta_grvi,
    nodos_sf = nodos_sf,
    lote_sf = lote_sf,
    extraccion_larga = tabla_extraccion_larga,
    extraccion_ancha = tabla_extraccion_ancha,
    correlaciones_uav_cultivo = tabla_correlaciones_uav_cultivo
  ),
  file = file.path(
    ruta_objetos,
    "objetos_UAV_GLI_GRVI.rds"
  )
)


## 28. Registro de sesión --------------------------------------

capture.output(
  sessionInfo(),
  file = file.path(
    ruta_uav,
    "sessionInfo_07_UAV.txt"
  )
)


## 29. Resumen en consola --------------------------------------

cat(
  "\n============================================\n"
)

cat(
  "SCRIPT 07_UAV FINALIZADO\n\n"
)

cat(
  "Ortomosaico 09:\n",
  archivo_09,
  "\n\n"
)

cat(
  "Ortomosaico 14:\n",
  archivo_14,
  "\n\n"
)

cat(
  "Desplazamiento óptimo aplicado al 09:\n"
)

print(
  mejor_desplazamiento
)

cat(
  "\nResumen de índices:\n"
)

print(
  tabla_resumen_indices
)

if (nrow(tabla_extraccion_larga) > 0) {

  cat(
    "\nNodos con extracción UAV:\n",
    n_distinct(
      tabla_extraccion_larga$Nodo
    ),
    "\n"
  )
}

cat(
  "\nResultados almacenados en:\n",
  ruta_uav,
  "\n"
)

cat(
  "============================================\n"
)
