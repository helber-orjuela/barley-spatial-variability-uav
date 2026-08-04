
## ============================================================
## 04_TRATAMIENTOS – SCRIPT 02
## INTEGRAR TRATAMIENTOS, CULTIVO Y SUPERFICIES DE SUELO
## ============================================================
##
## Entradas:
##   1. Tratamientos y cultivo:
##      04_TRATAMIENTOS/01_MATRIZ_MAESTRA/
##      nodos_tratamientos_lote_2/objetos_R/
##      objetos_tratamientos_56_nodos.rds
##
##   2. Modelos geoestadísticos finales de suelo I2:
##      03_ANISOTROPIA_I2/objetos_R/
##      objetos_anisotropia_I2.rds
##
##   3. Resultados de Kriging I2:
##      04_KRIGING_I2/objetos_R/
##      objetos_kriging_I2.rds
##
## Salida:
##   04_TRATAMIENTOS/02_INTEGRACION_SUELO_CULTIVO/
##
## Propósito:
##   1. Predecir directamente las variables de suelo en las
##      coordenadas de los 56 nodos de cultivo.
##   2. Asociar predicción, varianza y error estándar de Kriging
##      con tratamiento, dosis, omisión y repetición.
##   3. Calcular distancia al nodo de suelo más cercano.
##   4. Crear una base integrada para los modelos de tratamientos.
##
## CRITERIOS:
##   - La unidad de análisis es cada nodo/unidad observada.
##   - No se usan los valores ráster por grid_id.
##   - El Kriging se evalúa directamente en las coordenadas
##     de cultivo, evitando desajustes entre mallas.
##   - N, P y K se conservan como índices instrumentales nominales.
##   - pH y Hv se marcan como resultados exploratorios.
##   - La incertidumbre de calibración OLS de CE no está incluida
##     en la varianza de Kriging.
## ============================================================


## 0. Paquetes -------------------------------------------------

paquetes <- c(
  "sf",
  "gstat",
  "dplyr",
  "tidyr",
  "purrr",
  "ggplot2",
  "writexl",
  "readr",
  "patchwork"
)

instalar <- paquetes[
  !paquetes %in% rownames(installed.packages())
]

if (length(instalar) > 0) {
  install.packages(instalar)
}

library(sf)
library(gstat)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(writexl)
library(readr)
library(patchwork)


## 1. Rutas ----------------------------------------------------

ruta_base <- paste0(
  "D:/UDCA_2017_CLASE/46_ARTICULOS_2025/",
  "ARTICULO_CEBADA_DOCTORADO/AJUSTES_ARTÍCULO"
)

ruta_tratamientos <- file.path(
  ruta_base,
  "04_TRATAMIENTOS",
  "01_MATRIZ_MAESTRA",
  "nodos_tratamientos_lote_2",
  "objetos_R",
  "objetos_tratamientos_56_nodos.rds"
)

ruta_modelos_suelo <- file.path(
  ruta_base,
  "03_ANISOTROPIA_I2",
  "objetos_R",
  "objetos_anisotropia_I2.rds"
)

ruta_kriging <- file.path(
  ruta_base,
  "04_KRIGING_I2",
  "objetos_R",
  "objetos_kriging_I2.rds"
)

ruta_salida <- file.path(
  ruta_base,
  "04_TRATAMIENTOS",
  "02_INTEGRACION_SUELO_CULTIVO"
)

ruta_tablas <- file.path(
  ruta_salida,
  "tablas"
)

ruta_datos <- file.path(
  ruta_salida,
  "datos_procesados"
)

ruta_figuras <- file.path(
  ruta_salida,
  "figuras"
)

ruta_objetos <- file.path(
  ruta_salida,
  "objetos_R"
)

archivos_necesarios <- c(
  ruta_tratamientos,
  ruta_modelos_suelo,
  ruta_kriging
)

faltantes <- archivos_necesarios[
  !file.exists(
    archivos_necesarios
  )
]

if (length(faltantes) > 0) {
  stop(
    "Faltan los siguientes archivos:\n",
    paste(
      faltantes,
      collapse = "\n"
    )
  )
}

dir.create(
  ruta_tablas,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  ruta_datos,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  ruta_figuras,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  ruta_objetos,
  recursive = TRUE,
  showWarnings = FALSE
)


## 2. Cargar objetos -------------------------------------------

obj_trat <- readRDS(
  ruta_tratamientos
)

obj_modelos <- readRDS(
  ruta_modelos_suelo
)

obj_kriging <- readRDS(
  ruta_kriging
)

matriz_cultivo <- obj_trat$matriz_observada
cultivo_sf <- obj_trat$nodos_sf

suelo_sf <- obj_kriging$suelo_sf
lote <- obj_kriging$lote

modelos_finales <- obj_modelos$modelos_finales
tabla_modelos_finales <- obj_modelos$tabla_modelos_finales

variables_suelo <- intersect(
  names(
    modelos_finales
  ),
  names(
    suelo_sf
  )
)

if (length(variables_suelo) == 0) {
  stop(
    "No existen variables de suelo comunes entre modelos y datos."
  )
}


## 3. Homologar CRS y validar geometría ------------------------

if (is.na(st_crs(cultivo_sf))) {
  stop(
    "Los nodos de cultivo no tienen CRS."
  )
}

if (is.na(st_crs(suelo_sf))) {
  stop(
    "Los nodos de suelo no tienen CRS."
  )
}

crs_objetivo <- st_crs(
  suelo_sf
)

cultivo_sf <- st_transform(
  cultivo_sf,
  crs_objetivo
)

lote <- st_transform(
  lote,
  crs_objetivo
)

lote_union <- st_union(
  st_make_valid(
    lote
  )
)

dentro_lote <- lengths(
  st_intersects(
    cultivo_sf,
    lote_union
  )
) > 0

distancia_lote <- as.numeric(
  st_distance(
    cultivo_sf,
    lote_union
  )
)

tabla_coherencia <- matriz_cultivo |>
  transmute(
    Nodo,
    Experimental_Unit,
    Treatment_ID,
    X,
    Y,
    Dentro_lote_suelo_I2 = dentro_lote,
    Distancia_al_lote_m = distancia_lote
  )

resumen_coherencia <- tibble(
  Indicador = c(
    "Nodos de cultivo",
    "Nodos dentro o sobre el lote I2",
    "Nodos fuera del lote I2",
    "Porcentaje dentro del lote",
    "Distancia máxima externa (m)"
  ),
  Valor = c(
    nrow(
      cultivo_sf
    ),
    sum(
      dentro_lote
    ),
    sum(
      !dentro_lote
    ),
    round(
      100 *
        mean(
          dentro_lote
        ),
      2
    ),
    ifelse(
      any(
        !dentro_lote
      ),
      max(
        distancia_lote[
          !dentro_lote
        ],
        na.rm = TRUE
      ),
      0
    )
  )
)

print(
  resumen_coherencia
)


## 4. Distancia al nodo de suelo más cercano -------------------

matriz_distancias <- st_distance(
  cultivo_sf,
  suelo_sf
)

indice_vecino <- apply(
  matriz_distancias,
  1,
  which.min
)

distancia_vecino <- apply(
  matriz_distancias,
  1,
  min
)

nodo_suelo_id <- if (
  "Nodo" %in% names(
    suelo_sf
  )
) {
  suelo_sf$Nodo[
    indice_vecino
  ]
} else {
  indice_vecino
}

tabla_vecino_suelo <- tibble(
  Nodo = matriz_cultivo$Nodo,
  Nodo_suelo_mas_cercano = nodo_suelo_id,
  Distancia_suelo_mas_cercano_m = as.numeric(
    distancia_vecino
  )
)


## 5. Función de predicción puntual ----------------------------

predecir_variable <- function(
  variable,
  modelo,
  datos_suelo,
  puntos_cultivo
) {

  if (is.null(modelo)) {
    return(NULL)
  }

  formula_variable <- as.formula(
    paste(
      variable,
      "~ 1"
    )
  )

  pred <- try(
    krige(
      formula = formula_variable,
      locations = datos_suelo,
      newdata = puntos_cultivo,
      model = modelo,
      debug.level = 0
    ),
    silent = TRUE
  )

  if (inherits(pred, "try-error")) {
    warning(
      "Falló la extracción de ",
      variable,
      ": ",
      as.character(
        pred
      )
    )
    return(NULL)
  }

  tibble(
    Nodo = matriz_cultivo$Nodo,
    Variable_suelo = variable,
    Prediccion = pred$var1.pred,
    Varianza_kriging = pred$var1.var,
    Error_estandar_kriging = sqrt(
      pmax(
        pred$var1.var,
        0
      )
    )
  )
}


## 6. Predecir todas las variables -----------------------------

predicciones_lista <- map(
  variables_suelo,
  function(variable) {

    predecir_variable(
      variable = variable,
      modelo = modelos_finales[[variable]],
      datos_suelo = suelo_sf,
      puntos_cultivo = cultivo_sf
    )
  }
)

names(
  predicciones_lista
) <- variables_suelo

predicciones_lista <- predicciones_lista[
  !map_lgl(
    predicciones_lista,
    is.null
  )
]

if (length(predicciones_lista) == 0) {
  stop(
    "No se pudo predecir ninguna variable de suelo."
  )
}

predicciones_largas <- bind_rows(
  predicciones_lista
)


## 7. Base ancha de predicciones -------------------------------

pred_valores <- predicciones_largas |>
  select(
    Nodo,
    Variable_suelo,
    Prediccion
  ) |>
  pivot_wider(
    names_from = Variable_suelo,
    values_from = Prediccion,
    names_glue = "Suelo_{Variable_suelo}"
  )

pred_varianzas <- predicciones_largas |>
  select(
    Nodo,
    Variable_suelo,
    Varianza_kriging
  ) |>
  pivot_wider(
    names_from = Variable_suelo,
    values_from = Varianza_kriging,
    names_glue = "VarK_{Variable_suelo}"
  )

pred_errores <- predicciones_largas |>
  select(
    Nodo,
    Variable_suelo,
    Error_estandar_kriging
  ) |>
  pivot_wider(
    names_from = Variable_suelo,
    values_from = Error_estandar_kriging,
    names_glue = "SEK_{Variable_suelo}"
  )


## 8. Construir base integrada ---------------------------------

base_integrada <- matriz_cultivo |>
  left_join(
    tabla_vecino_suelo,
    by = "Nodo"
  ) |>
  left_join(
    tabla_coherencia |>
      select(
        Nodo,
        Dentro_lote_suelo_I2,
        Distancia_al_lote_m
      ),
    by = "Nodo"
  ) |>
  left_join(
    pred_valores,
    by = "Nodo"
  ) |>
  left_join(
    pred_varianzas,
    by = "Nodo"
  ) |>
  left_join(
    pred_errores,
    by = "Nodo"
  ) |>
  mutate(
    Flag_pH = "Exploratorio_sensor_no_validado",
    Flag_Hv = "Exploratorio_estructura_espacial_limitada",
    Flag_NPK = "Indices_instrumentales_nominales",
    Flag_CE_OLS = "Calibracion_exploratoria_sin_incertidumbre_adicional",
    Diferencia_dias_suelo_cultivo = 5
  )


## 9. Resumen por tratamiento ----------------------------------

variables_cultivo <- c(
  "Macollas",
  "Altura",
  "Diametro",
  "NDVI",
  "Nitratos",
  "Potasio_savia",
  "Calcio_savia"
)

columnas_suelo_pred <- names(
  base_integrada
)[
  str_detect(
    names(
      base_integrada
    ),
    "^Suelo_"
  )
]

resumen_tratamiento <- base_integrada |>
  group_by(
    Treatment_ID,
    Treatment_base,
    Dose_code,
    Dose_percent,
    Omission_description
  ) |>
  summarise(
    n = n(),
    Repeticiones = paste(
      sort(
        unique(
          Replication
        )
      ),
      collapse = ", "
    ),
    across(
      all_of(
        c(
          variables_cultivo,
          columnas_suelo_pred
        )
      ),
      list(
        Media = ~ mean(
          .x,
          na.rm = TRUE
        ),
        DE = ~ sd(
          .x,
          na.rm = TRUE
        )
      ),
      .names = "{.col}_{.fn}"
    ),
    Distancia_suelo_media_m = mean(
      Distancia_suelo_mas_cercano_m,
      na.rm = TRUE
    ),
    .groups = "drop"
  )


## 10. Tabla de incertidumbre ----------------------------------

tabla_incertidumbre <- predicciones_largas |>
  group_by(
    Variable_suelo
  ) |>
  summarise(
    n = n(),
    Pred_min = min(
      Prediccion,
      na.rm = TRUE
    ),
    Pred_media = mean(
      Prediccion,
      na.rm = TRUE
    ),
    Pred_max = max(
      Prediccion,
      na.rm = TRUE
    ),
    Varianza_media = mean(
      Varianza_kriging,
      na.rm = TRUE
    ),
    Error_estandar_medio = mean(
      Error_estandar_kriging,
      na.rm = TRUE
    ),
    Error_estandar_max = max(
      Error_estandar_kriging,
      na.rm = TRUE
    ),
    .groups = "drop"
  )


## 11. Figuras -------------------------------------------------

lote_borde <- st_boundary(
  lote_union
)

g_coherencia <- ggplot() +
  geom_sf(
    data = lote_borde,
    fill = NA,
    linewidth = 0.55
  ) +
  geom_sf(
    data = cultivo_sf |>
      mutate(
        Estado = ifelse(
          dentro_lote,
          "Dentro",
          "Fuera"
        )
      ),
    aes(
      shape = Estado
    ),
    size = 2.2
  ) +
  coord_sf(
    expand = TRUE
  ) +
  labs(
    title =
      "Coherencia espacial de los nodos de cultivo",
    subtitle = paste0(
      sum(
        dentro_lote
      ),
      " de ",
      nrow(
        cultivo_sf
      ),
      " nodos dentro o sobre el lote de suelo I2"
    ),
    x = "X (m)",
    y = "Y (m)",
    shape = "Estado"
  ) +
  theme_minimal()

ggsave(
  filename = file.path(
    ruta_figuras,
    "Coherencia_nodos_cultivo_lote_suelo.png"
  ),
  plot = g_coherencia,
  width = 8,
  height = 7,
  dpi = 400
)

datos_figura <- base_integrada |>
  select(
    Nodo,
    X,
    Y,
    Treatment_ID,
    Dose_percent,
    all_of(
      columnas_suelo_pred
    )
  ) |>
  pivot_longer(
    cols = all_of(
      columnas_suelo_pred
    ),
    names_to = "Variable",
    values_to = "Prediccion"
  ) |>
  mutate(
    Variable = str_remove(
      Variable,
      "^Suelo_"
    )
  )

g_superficies_puntos <- ggplot(
  datos_figura,
  aes(
    x = X,
    y = Y,
    fill = Prediccion
  )
) +
  geom_point(
    shape = 21,
    size = 2.7
  ) +
  facet_wrap(
    ~ Variable,
    ncol = 3
  ) +
  coord_equal() +
  labs(
    title =
      "Variables de suelo predichas en los 56 nodos de cultivo",
    subtitle =
      "Predicción puntual mediante los modelos geoestadísticos de I2",
    x = "X (m)",
    y = "Y (m)",
    fill = "Predicción"
  ) +
  theme_minimal()

ggsave(
  filename = file.path(
    ruta_figuras,
    "Suelo_predicho_en_nodos_cultivo.png"
  ),
  plot = g_superficies_puntos,
  width = 13,
  height = 10,
  dpi = 400
)

g_distancias <- ggplot(
  base_integrada,
  aes(
    x = Distancia_suelo_mas_cercano_m
  )
) +
  geom_histogram(
    bins = 10,
    color = "white"
  ) +
  labs(
    title =
      "Distancia de los nodos de cultivo al nodo de suelo más cercano",
    x = "Distancia (m)",
    y = "Frecuencia"
  ) +
  theme_minimal()

ggsave(
  filename = file.path(
    ruta_figuras,
    "Distancia_cultivo_nodo_suelo.png"
  ),
  plot = g_distancias,
  width = 8,
  height = 6,
  dpi = 400
)


## 12. Exportar -------------------------------------------------

write_csv(
  base_integrada,
  file.path(
    ruta_datos,
    "Base_integrada_tratamientos_suelo_cultivo_56_nodos.csv"
  )
)

write_csv(
  predicciones_largas,
  file.path(
    ruta_datos,
    "Predicciones_suelo_en_nodos_cultivo_formato_largo.csv"
  )
)

write_csv(
  resumen_tratamiento,
  file.path(
    ruta_datos,
    "Resumen_integrado_por_tratamiento.csv"
  )
)

write_xlsx(
  list(
    Resumen_coherencia =
      resumen_coherencia,

    Coherencia_nodos =
      tabla_coherencia,

    Distancia_nodo_suelo =
      tabla_vecino_suelo,

    Incertidumbre_prediccion =
      tabla_incertidumbre,

    Modelos_suelo =
      tabla_modelos_finales,

    Base_integrada =
      base_integrada,

    Predicciones_largas =
      predicciones_largas,

    Resumen_tratamiento =
      resumen_tratamiento
  ),
  path = file.path(
    ruta_tablas,
    "Resultados_integracion_suelo_cultivo_tratamientos.xlsx"
  )
)


## 13. Guardar objeto ------------------------------------------

saveRDS(
  list(
    base_integrada =
      base_integrada,

    cultivo_sf =
      cultivo_sf,

    suelo_sf =
      suelo_sf,

    lote =
      lote_union,

    predicciones_largas =
      predicciones_largas,

    resumen_tratamiento =
      resumen_tratamiento,

    tabla_incertidumbre =
      tabla_incertidumbre,

    tabla_coherencia =
      tabla_coherencia
  ),
  file = file.path(
    ruta_objetos,
    "objetos_integracion_suelo_cultivo_tratamientos.rds"
  )
)


## 14. Registro de sesión --------------------------------------

capture.output(
  sessionInfo(),
  file = file.path(
    ruta_salida,
    "sessionInfo_04_02.txt"
  )
)


## 15. Resumen en consola --------------------------------------

cat(
  "\n============================================\n"
)

cat(
  "SCRIPT 04-02 FINALIZADO\n\n"
)

cat(
  "Nodos integrados:\n",
  nrow(
    base_integrada
  ),
  "\n"
)

cat(
  "Variables de suelo extraídas:\n",
  paste(
    variables_suelo,
    collapse = ", "
  ),
  "\n"
)

cat(
  "Distancia media al nodo de suelo más cercano:\n",
  round(
    mean(
      base_integrada$Distancia_suelo_mas_cercano_m,
      na.rm = TRUE
    ),
    3
  ),
  " m\n"
)

cat(
  "Objeto para Script 03:\n",
  file.path(
    ruta_objetos,
    "objetos_integracion_suelo_cultivo_tratamientos.rds"
  ),
  "\n"
)

cat(
  "============================================\n"
)
