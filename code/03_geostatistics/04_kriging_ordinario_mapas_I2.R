
## ============================================================
## SCRIPT 04. KRIGING ORDINARIO Y MAPAS DE INCERTIDUMBRE – I2
## Artículo cebada: revisión mayor
##
## Entrada:
##   AJUSTES_ARTÍCULO/03_ANISOTROPIA_I2/objetos_R/
##   objetos_anisotropia_I2.rds
##
## Salidas:
##   AJUSTES_ARTÍCULO/04_KRIGING_I2/
##
## Objetivos:
##   1. Generar una malla regular recortada al lote.
##   2. Ejecutar Kriging ordinario con el modelo final recomendado.
##   3. Exportar mapas de predicción y varianza de Kriging.
##   4. Comparar CE original y CE corregida OLS.
##   5. Mantener pH como resultado diagnóstico, no como inferencia
##      agronómica principal.
##
## IMPORTANTE:
##   - N, P y K se interpretan como índices instrumentales nominales.
##   - CE corregida OLS es una estimación exploratoria calibrada.
##   - La incertidumbre de calibración de CE no está incluida dentro
##     de la varianza de Kriging; se reporta por separado.
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
  "patchwork",
  "readr"
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
library(patchwork)
library(readr)


## 1. Rutas ----------------------------------------------------

ruta_base <- paste0(
  "D:/UDCA_2017_CLASE/46_ARTICULOS_2025/",
  "ARTICULO_CEBADA_DOCTORADO/AJUSTES_ARTÍCULO"
)

ruta_objeto <- file.path(
  ruta_base,
  "03_ANISOTROPIA_I2",
  "objetos_R",
  "objetos_anisotropia_I2.rds"
)

ruta_salida <- file.path(
  ruta_base,
  "04_KRIGING_I2"
)

ruta_tablas <- file.path(
  ruta_salida,
  "tablas"
)

ruta_figuras <- file.path(
  ruta_salida,
  "figuras"
)

ruta_raster <- file.path(
  ruta_salida,
  "superficies_csv"
)

ruta_objetos <- file.path(
  ruta_salida,
  "objetos_R"
)

if (!file.exists(ruta_objeto)) {
  stop(
    "No existe el objeto del Script 03:\n",
    ruta_objeto,
    "\nEjecuta primero el Script 03."
  )
}

dir.create(
  ruta_tablas,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  ruta_figuras,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  ruta_raster,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  ruta_objetos,
  recursive = TRUE,
  showWarnings = FALSE
)


## 2. Cargar objetos del Script 03 -----------------------------

obj <- readRDS(
  ruta_objeto
)

suelo_sf <- obj$suelo_sf
lote <- obj$lote
variables_modelar <- obj$variables
modelos_finales <- obj$modelos_finales
tabla_modelos_finales <- obj$tabla_modelos_finales
tabla_comparacion_iso_anis <- obj$tabla_comparacion_iso_anis

if (length(modelos_finales) == 0) {
  stop(
    "No existen modelos finales disponibles para Kriging."
  )
}

# Validación básica de CRS.
if (is.na(st_crs(suelo_sf))) {
  stop(
    "Los puntos de suelo no tienen CRS."
  )
}

if (is.na(st_crs(lote))) {
  stop(
    "El lote no tiene CRS."
  )
}

lote <- st_transform(
  lote,
  st_crs(suelo_sf)
)


## 3. Clasificación de variables -------------------------------

clasificacion_variables <- tibble(
  Variable = variables_modelar,
  Categoria = c(
    "Índice instrumental nominal",
    "Índice instrumental nominal",
    "Índice instrumental nominal",
    "Diagnóstico exploratorio",
    "Variable física",
    "Variable física",
    "Variable proximal principal",
    "Estimación exploratoria calibrada"
  ),
  Incluir_panel_principal = c(
    FALSE,
    FALSE,
    FALSE,
    FALSE,
    TRUE,
    TRUE,
    TRUE,
    TRUE
  ),
  Nota = c(
    "Redundante con CE; presentar en suplemento",
    "Redundante con CE; presentar en suplemento",
    "Redundante con CE; presentar en suplemento",
    "Sensor no validado frente a laboratorio",
    "Interpretar con cautela por baja variabilidad",
    "Interpretación principal",
    "Lectura proximal sin calibrar",
    "Transformación OLS específica del conjunto"
  )
)


## 4. Crear malla regular recortada al lote ---------------------

bbox_lote <- st_bbox(
  lote
)

ancho_lote <- as.numeric(
  bbox_lote["xmax"] -
    bbox_lote["xmin"]
)

alto_lote <- as.numeric(
  bbox_lote["ymax"] -
    bbox_lote["ymin"]
)

# Aproximadamente 120 celdas a lo largo de la dimensión mayor.
# Se evita una resolución excesivamente fina para un lote pequeño.
dimension_mayor <- max(
  ancho_lote,
  alto_lote
)

cellsize_m <- dimension_mayor / 120

grid_centros <- st_make_grid(
  lote,
  cellsize = cellsize_m,
  what = "centers",
  square = TRUE
)

grid_sf <- st_sf(
  geometry = grid_centros
)

grid_sf <- grid_sf[
  lengths(
    st_intersects(
      grid_sf,
      lote
    )
  ) > 0,
]

if (nrow(grid_sf) == 0) {
  stop(
    "La malla de predicción quedó vacía."
  )
}

coords_grid <- st_coordinates(
  grid_sf
)

grid_sf$grid_id <- seq_len(
  nrow(grid_sf)
)

tabla_malla <- tibble(
  Indicador = c(
    "Ancho del lote (m)",
    "Alto del lote (m)",
    "Tamaño de celda (m)",
    "Número de centros de predicción"
  ),
  Valor = c(
    ancho_lote,
    alto_lote,
    cellsize_m,
    nrow(grid_sf)
  )
)

print(
  tabla_malla
)


## 5. Función de Kriging ---------------------------------------

ejecutar_kriging <- function(
    variable,
    datos_sf,
    grid_sf,
    modelo
) {
  
  formula_variable <- as.formula(
    paste(
      variable,
      "~ 1"
    )
  )
  
  resultado <- try(
    krige(
      formula = formula_variable,
      locations = datos_sf,
      newdata = grid_sf,
      model = modelo,
      debug.level = 0
    ),
    silent = TRUE
  )
  
  if (inherits(resultado, "try-error")) {
    return(
      list(
        exito = FALSE,
        resultado_sf = NULL,
        resultado_df = NULL,
        error = as.character(resultado)
      )
    )
  }
  
  coords <- st_coordinates(
    resultado
  )
  
  # Verificar que el número de predicciones coincida con la malla
  if (nrow(resultado) != nrow(grid_sf)) {
    stop(
      "El número de predicciones no coincide con el número ",
      "de puntos de la malla para la variable: ",
      variable
    )
  }
  
  resultado_df <- resultado |>
    st_drop_geometry() |>
    mutate(
      grid_id = grid_sf$grid_id,
      X = coords[, 1],
      Y = coords[, 2],
      Variable = variable,
      Prediccion = var1.pred,
      Varianza_kriging = var1.var,
      Error_estandar_kriging = sqrt(
        pmax(
          var1.var,
          0
        )
      ),
      .before = 1
    )
  
  list(
    exito = TRUE,
    resultado_sf = resultado,
    resultado_df = resultado_df,
    error = NA_character_
  )
}

## 6. Ejecutar Kriging para todas las variables ----------------

resultados_kriging <- list()
errores_kriging <- list()

for (variable in variables_modelar) {

  message(
    "Kriging de ",
    variable
  )

  modelo <- modelos_finales[[variable]]

  if (is.null(modelo)) {

    errores_kriging[[variable]] <- tibble(
      Variable = variable,
      Error = "No existe modelo final"
    )

    next
  }

  resultado <- ejecutar_kriging(
    variable = variable,
    datos_sf = suelo_sf,
    grid_sf = grid_sf,
    modelo = modelo
  )

  if (!resultado$exito) {

    errores_kriging[[variable]] <- tibble(
      Variable = variable,
      Error = resultado$error
    )

    next
  }

  resultados_kriging[[variable]] <- resultado

  write_csv(
    resultado$resultado_df,
    file.path(
      ruta_raster,
      paste0(
        "Superficie_",
        variable,
        "_I2.csv"
      )
    )
  )
}


tabla_errores_kriging <- bind_rows(
  errores_kriging
)

if (length(resultados_kriging) == 0) {
  stop(
    "No se pudo interpolar ninguna variable."
  )
}


## 7. Resumen de superficies ----------------------------------

tabla_resumen_superficies <- imap_dfr(
  resultados_kriging,
  function(resultado, variable) {

    df <- resultado$resultado_df

    tibble(
      Variable = variable,
      n_celdas = nrow(df),
      Pred_min = min(
        df$Prediccion,
        na.rm = TRUE
      ),
      Pred_media = mean(
        df$Prediccion,
        na.rm = TRUE
      ),
      Pred_mediana = median(
        df$Prediccion,
        na.rm = TRUE
      ),
      Pred_max = max(
        df$Prediccion,
        na.rm = TRUE
      ),
      Varianza_min = min(
        df$Varianza_kriging,
        na.rm = TRUE
      ),
      Varianza_media = mean(
        df$Varianza_kriging,
        na.rm = TRUE
      ),
      Varianza_max = max(
        df$Varianza_kriging,
        na.rm = TRUE
      ),
      Error_estandar_medio = mean(
        df$Error_estandar_kriging,
        na.rm = TRUE
      )
    )
  }
)

print(
  tabla_resumen_superficies
)


## 8. Preparar borde y puntos observados -----------------------

lote_borde <- st_boundary(
  st_union(
    lote
  )
)

puntos_df <- suelo_sf |>
  st_drop_geometry() |>
  select(
    Nodo,
    X,
    Y,
    all_of(
      variables_modelar
    )
  ) |>
  pivot_longer(
    cols = all_of(
      variables_modelar
    ),
    names_to = "Variable",
    values_to = "Observado"
  )


## 9. Funciones gráficas ---------------------------------------

crear_mapa_prediccion <- function(
  variable,
  resultado_df
) {

  puntos_variable <- puntos_df |>
    filter(
      Variable == variable
    )

  ggplot() +
    geom_raster(
      data = resultado_df,
      aes(
        x = X,
        y = Y,
        fill = Prediccion
      )
    ) +
    geom_contour(
      data = resultado_df,
      aes(
        x = X,
        y = Y,
        z = Prediccion
      ),
      bins = 8,
      linewidth = 0.25,
      alpha = 0.55
    ) +
    geom_sf(
      data = lote_borde,
      fill = NA,
      linewidth = 0.45
    ) +
    geom_point(
      data = puntos_variable,
      aes(
        x = X,
        y = Y
      ),
      size = 0.9
    ) +
    coord_sf(
      expand = FALSE
    ) +
    labs(
      title = variable,
      subtitle = paste0(
        "Modelo ",
        tabla_modelos_finales |>
          filter(
            Variable == variable
          ) |>
          pull(
            Tipo_modelo_final
          )
      ),
      x = "X (m)",
      y = "Y (m)",
      fill = "Predicción"
    ) +
    theme_minimal() +
    theme(
      panel.grid.minor = element_blank()
    )
}


crear_mapa_varianza <- function(
  variable,
  resultado_df
) {

  ggplot() +
    geom_raster(
      data = resultado_df,
      aes(
        x = X,
        y = Y,
        fill = Varianza_kriging
      )
    ) +
    geom_sf(
      data = lote_borde,
      fill = NA,
      linewidth = 0.45
    ) +
    geom_point(
      data = suelo_sf |>
        st_drop_geometry(),
      aes(
        x = X,
        y = Y
      ),
      size = 0.8
    ) +
    coord_sf(
      expand = FALSE
    ) +
    labs(
      title = variable,
      x = "X (m)",
      y = "Y (m)",
      fill = "Varianza"
    ) +
    theme_minimal() +
    theme(
      panel.grid.minor = element_blank()
    )
}


## 10. Generar figuras individuales ----------------------------

mapas_prediccion <- list()
mapas_varianza <- list()

for (variable in names(resultados_kriging)) {

  resultado_df <- resultados_kriging[[variable]]$resultado_df

  mapa_pred <- crear_mapa_prediccion(
    variable,
    resultado_df
  )

  mapa_var <- crear_mapa_varianza(
    variable,
    resultado_df
  )

  mapas_prediccion[[variable]] <- mapa_pred
  mapas_varianza[[variable]] <- mapa_var

  ggsave(
    filename = file.path(
      ruta_figuras,
      paste0(
        "Prediccion_",
        variable,
        "_I2.png"
      )
    ),
    plot = mapa_pred,
    width = 7,
    height = 6,
    dpi = 400
  )

  ggsave(
    filename = file.path(
      ruta_figuras,
      paste0(
        "Varianza_",
        variable,
        "_I2.png"
      )
    ),
    plot = mapa_var,
    width = 7,
    height = 6,
    dpi = 400
  )
}


## 11. Panel principal -----------------------------------------

variables_principales <- clasificacion_variables |>
  filter(
    Incluir_panel_principal
  ) |>
  pull(
    Variable
  )

variables_principales <- intersect(
  variables_principales,
  names(
    mapas_prediccion
  )
)

panel_principal_prediccion <- wrap_plots(
  mapas_prediccion[
    variables_principales
  ],
  ncol = 2
) +
  plot_annotation(
    title =
      "Kriging ordinario de variables principales de suelo en I2"
  )

panel_principal_varianza <- wrap_plots(
  mapas_varianza[
    variables_principales
  ],
  ncol = 2
) +
  plot_annotation(
    title =
      "Varianza de Kriging de variables principales de suelo en I2"
  )

ggsave(
  filename = file.path(
    ruta_figuras,
    "Panel_principal_prediccion_I2.png"
  ),
  plot = panel_principal_prediccion,
  width = 13,
  height = 11,
  dpi = 400
)

ggsave(
  filename = file.path(
    ruta_figuras,
    "Panel_principal_varianza_I2.png"
  ),
  plot = panel_principal_varianza,
  width = 13,
  height = 11,
  dpi = 400
)


## 12. Panel suplementario N, P, K y pH ------------------------

variables_suplementarias <- c(
  "N",
  "P",
  "K",
  "pH"
)

variables_suplementarias <- intersect(
  variables_suplementarias,
  names(
    mapas_prediccion
  )
)

panel_suplementario <- wrap_plots(
  mapas_prediccion[
    variables_suplementarias
  ],
  ncol = 2
) +
  plot_annotation(
    title =
      "Superficies suplementarias de índices instrumentales y pH"
  )

ggsave(
  filename = file.path(
    ruta_figuras,
    "Panel_suplementario_NPK_pH_I2.png"
  ),
  plot = panel_suplementario,
  width = 13,
  height = 11,
  dpi = 400
)


## 13. Comparación CE original y corregida ---------------------

if (
  all(
    c(
      "CE",
      "CE_corregida_OLS"
    ) %in% names(
      resultados_kriging
    )
  )
) {

  ce_original <- resultados_kriging[["CE"]]$resultado_df |>
    select(
      grid_id,
      X,
      Y,
      CE_pred = Prediccion,
      CE_var = Varianza_kriging
    )

  ce_corregida <- resultados_kriging[["CE_corregida_OLS"]]$resultado_df |>
    select(
      grid_id,
      CE_corregida_pred = Prediccion,
      CE_corregida_var = Varianza_kriging
    )

  comparacion_ce <- ce_original |>
    left_join(
      ce_corregida,
      by = "grid_id"
    ) |>
    mutate(
      Diferencia_prediccion =
        CE_corregida_pred -
        CE_pred,
      Relacion_varianzas =
        CE_corregida_var /
        CE_var
    )

  resumen_ce <- tibble(
    Indicador = c(
      "Correlación espacial entre predicciones",
      "Diferencia media corregida-original",
      "Diferencia mínima",
      "Diferencia máxima",
      "Relación media de varianzas"
    ),
    Valor = c(
      cor(
        comparacion_ce$CE_pred,
        comparacion_ce$CE_corregida_pred,
        use = "complete.obs"
      ),
      mean(
        comparacion_ce$Diferencia_prediccion,
        na.rm = TRUE
      ),
      min(
        comparacion_ce$Diferencia_prediccion,
        na.rm = TRUE
      ),
      max(
        comparacion_ce$Diferencia_prediccion,
        na.rm = TRUE
      ),
      mean(
        comparacion_ce$Relacion_varianzas[
          is.finite(
            comparacion_ce$Relacion_varianzas
          )
        ],
        na.rm = TRUE
      )
    )
  )

  mapa_diferencia_ce <- ggplot() +
    geom_raster(
      data = comparacion_ce,
      aes(
        x = X,
        y = Y,
        fill = Diferencia_prediccion
      )
    ) +
    geom_sf(
      data = lote_borde,
      fill = NA,
      linewidth = 0.45
    ) +
    coord_sf(
      expand = FALSE
    ) +
    labs(
      title =
        "Diferencia entre CE corregida y CE original",
      subtitle =
        "La transformación OLS conserva el patrón espacial",
      x = "X (m)",
      y = "Y (m)",
      fill = "Δ CE"
    ) +
    theme_minimal()

  panel_ce <- (
    mapas_prediccion[["CE"]] |
      mapas_prediccion[["CE_corregida_OLS"]]
  ) / mapa_diferencia_ce

  ggsave(
    filename = file.path(
      ruta_figuras,
      "Comparacion_CE_original_corregida_I2.png"
    ),
    plot = panel_ce,
    width = 13,
    height = 12,
    dpi = 400
  )

  write_csv(
    comparacion_ce,
    file.path(
      ruta_raster,
      "Comparacion_CE_original_corregida_I2.csv"
    )
  )

} else {

  resumen_ce <- tibble(
    Indicador =
      "No fue posible comparar CE",
    Valor =
      NA_real_
  )

  comparacion_ce <- NULL
}


## 14. Exportar tabla general ----------------------------------

write_xlsx(
  list(
    Configuracion_malla =
      tabla_malla,

    Clasificacion_variables =
      clasificacion_variables,

    Modelos_finales =
      tabla_modelos_finales,

    Comparacion_iso_anis =
      tabla_comparacion_iso_anis,

    Resumen_superficies =
      tabla_resumen_superficies,

    Comparacion_CE =
      resumen_ce,

    Errores_kriging =
      tabla_errores_kriging
  ),
  path = file.path(
    ruta_tablas,
    "Resultados_kriging_I2.xlsx"
  )
)


## 15. Guardar objetos para Script 05 --------------------------

saveRDS(
  list(
    suelo_sf =
      suelo_sf,

    lote =
      lote,

    grid_sf =
      grid_sf,

    resultados_kriging =
      resultados_kriging,

    tabla_resumen_superficies =
      tabla_resumen_superficies,

    tabla_modelos_finales =
      tabla_modelos_finales,

    clasificacion_variables =
      clasificacion_variables,

    comparacion_ce =
      comparacion_ce
  ),
  file = file.path(
    ruta_objetos,
    "objetos_kriging_I2.rds"
  )
)


## 16. Registro de sesión --------------------------------------

capture.output(
  sessionInfo(),
  file = file.path(
    ruta_salida,
    "sessionInfo_Script04.txt"
  )
)


## 17. Resumen en consola --------------------------------------

cat(
  "\n============================================\n"
)

cat(
  "SCRIPT 04 FINALIZADO\n\n"
)

cat(
  "Variables interpoladas:\n",
  paste(
    names(
      resultados_kriging
    ),
    collapse = ", "
  ),
  "\n\n"
)

cat(
  "Resultados guardados en:\n",
  ruta_salida,
  "\n\n"
)

cat(
  "Objeto para Script 05:\n",
  file.path(
    ruta_objetos,
    "objetos_kriging_I2.rds"
  ),
  "\n"
)

cat(
  "============================================\n"
)
