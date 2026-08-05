
## ============================================================
## I1 – SCRIPT C REFORMULADO
## KRIGING COMPACTO Y COMPARACIÓN TEMPORAL I1–I2
## ============================================================
##
## CAMBIO PRINCIPAL:
##   I1 utiliza exclusivamente el polígono ubicado en:
##
##   AJUSTES_ARTÍCULO/DATOS BASE/LOTE_I1
##
##   El script busca automáticamente el primer archivo .shp
##   disponible dentro de esa carpeta.
##
## OBJETIVOS:
##   1. Verificar coherencia espacial entre nodos I1 y lote I1.
##   2. Generar Kriging I1 recortado al polígono correcto.
##   3. Evitar comparar I1 e I2 mediante grid_id de mallas distintas.
##   4. Construir una malla común únicamente sobre el área de
##      intersección entre los polígonos I1 e I2.
##   5. Recalcular I1 e I2 sobre esa malla común antes de obtener
##      diferencias espaciales.
##
## IMPORTANTE:
##   - No se desplazan ni corrigen coordenadas automáticamente.
##   - Los puntos por fuera del lote se reportan para auditoría.
##   - La comparación I1–I2 es exploratoria, no una serie temporal formal.
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

ruta_i1 <- file.path(
  ruta_base,
  "06_I1_COMPACTO",
  "02_SEMIVARIOGRAMAS",
  "objetos_R",
  "objetos_semivariogramas_I1.rds"
)

ruta_i2_kriging <- file.path(
  ruta_base,
  "04_KRIGING_I2",
  "objetos_R",
  "objetos_kriging_I2.rds"
)

ruta_i2_modelos <- file.path(
  ruta_base,
  "03_ANISOTROPIA_I2",
  "objetos_R",
  "objetos_anisotropia_I2.rds"
)

ruta_carpeta_lote_i1 <- file.path(
  ruta_base,
  "DATOS BASE",
  "LOTE_I1"
)

ruta_salida <- file.path(
  ruta_base,
  "06_I1_COMPACTO",
  "03_KRIGING_TEMPORAL"
)

ruta_figuras <- file.path(
  ruta_salida,
  "figuras"
)

ruta_tablas <- file.path(
  ruta_salida,
  "tablas"
)

ruta_superficies <- file.path(
  ruta_salida,
  "superficies_csv"
)

ruta_objetos <- file.path(
  ruta_salida,
  "objetos_R"
)

ruta_auditoria <- file.path(
  ruta_salida,
  "auditoria_espacial"
)


## 2. Verificar archivos ---------------------------------------

if (!file.exists(ruta_i1)) {
  stop(
    "No existe el objeto I1:\n",
    ruta_i1
  )
}

if (!file.exists(ruta_i2_kriging)) {
  stop(
    "No existe el objeto de Kriging I2:\n",
    ruta_i2_kriging
  )
}

if (!file.exists(ruta_i2_modelos)) {
  stop(
    "No existe el objeto de modelos I2:\n",
    ruta_i2_modelos
  )
}

if (!dir.exists(ruta_carpeta_lote_i1)) {
  stop(
    "No existe la carpeta del lote I1:\n",
    ruta_carpeta_lote_i1
  )
}

archivos_shp_i1 <- list.files(
  ruta_carpeta_lote_i1,
  pattern = "\\.shp$",
  full.names = TRUE,
  ignore.case = TRUE
)

if (length(archivos_shp_i1) == 0) {
  stop(
    "No se encontró ningún archivo .shp dentro de:\n",
    ruta_carpeta_lote_i1
  )
}

if (length(archivos_shp_i1) > 1) {
  warning(
    "Se encontraron varios shapefiles en LOTE_I1. ",
    "Se utilizará el primero:\n",
    archivos_shp_i1[1]
  )
}

ruta_lote_i1 <- archivos_shp_i1[1]


## 3. Crear carpetas de salida ---------------------------------

dir.create(
  ruta_figuras,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  ruta_tablas,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  ruta_superficies,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  ruta_objetos,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  ruta_auditoria,
  recursive = TRUE,
  showWarnings = FALSE
)


## 4. Cargar objetos y lote I1 ---------------------------------

obj_i1 <- readRDS(
  ruta_i1
)

obj_i2_k <- readRDS(
  ruta_i2_kriging
)

obj_i2_m <- readRDS(
  ruta_i2_modelos
)

i1_sf <- obj_i1$i1_sf
variables_i1 <- obj_i1$variables
modelos_i1 <- obj_i1$modelos_finales
tabla_modelos_i1 <- obj_i1$tabla_seleccionados

i2_sf <- obj_i2_k$suelo_sf
lote_i2 <- obj_i2_k$lote
tabla_superficies_i2 <- obj_i2_k$tabla_resumen_superficies
modelos_i2 <- obj_i2_m$modelos_finales

lote_i1 <- st_read(
  ruta_lote_i1,
  quiet = TRUE
)


## 5. Homologar CRS --------------------------------------------

if (is.na(st_crs(i1_sf))) {
  stop(
    "Los nodos I1 no tienen CRS definido."
  )
}

crs_objetivo <- st_crs(
  i1_sf
)

if (is.na(st_crs(lote_i1))) {
  stop(
    "El shapefile del lote I1 no tiene CRS definido. ",
    "Asigna correctamente el CRS antes de continuar."
  )
}

lote_i1 <- st_transform(
  lote_i1,
  crs_objetivo
)

i2_sf <- st_transform(
  i2_sf,
  crs_objetivo
)

lote_i2 <- st_transform(
  lote_i2,
  crs_objetivo
)

# Reparar geometrías si fuera necesario.
lote_i1 <- st_make_valid(
  lote_i1
)

lote_i2 <- st_make_valid(
  lote_i2
)

# Unificar multipartes.
lote_i1_union <- st_union(
  lote_i1
)

lote_i2_union <- st_union(
  lote_i2
)


## 6. Auditoría espacial nodos I1 – lote I1 --------------------

dentro_i1 <- lengths(
  st_intersects(
    i1_sf,
    lote_i1_union
  )
) > 0

distancia_lote_i1 <- as.numeric(
  st_distance(
    i1_sf,
    lote_i1_union
  )
)

tabla_coherencia_i1 <- i1_sf |>
  st_drop_geometry() |>
  transmute(
    Nodo,
    X,
    Y,
    Dentro_lote_I1 = dentro_i1,
    Distancia_al_lote_m = distancia_lote_i1
  )

resumen_coherencia_i1 <- tibble(
  Indicador = c(
    "Número total de nodos I1",
    "Nodos dentro o sobre el borde",
    "Nodos fuera del lote",
    "Porcentaje dentro del lote",
    "Distancia máxima de nodos externos al lote (m)",
    "CRS nodos I1",
    "CRS lote I1",
    "Shapefile utilizado"
  ),
  Valor = c(
    nrow(i1_sf),
    sum(dentro_i1),
    sum(!dentro_i1),
    round(
      100 * mean(dentro_i1),
      2
    ),
    ifelse(
      any(!dentro_i1),
      round(
        max(
          distancia_lote_i1[!dentro_i1],
          na.rm = TRUE
        ),
        4
      ),
      0
    ),
    st_crs(i1_sf)$input,
    st_crs(lote_i1)$input,
    basename(ruta_lote_i1)
  )
)

print(
  resumen_coherencia_i1
)

if (any(!dentro_i1)) {
  warning(
    sum(!dentro_i1),
    " nodos I1 están fuera del lote I1. ",
    "El script no modifica sus coordenadas; revise la auditoría."
  )
}


## 7. Figura de coherencia espacial ----------------------------

g_coherencia <- ggplot() +
  geom_sf(
    data = lote_i1_union,
    fill = NA,
    linewidth = 0.8
  ) +
  geom_sf(
    data = i1_sf |>
      mutate(
        Estado = ifelse(
          dentro_i1,
          "Dentro",
          "Fuera"
        )
      ),
    aes(
      shape = Estado
    ),
    size = 2.6
  ) +
  coord_sf(
    expand = TRUE
  ) +
  labs(
    title =
      "Coherencia espacial entre nodos y polígono de I1",
    subtitle = paste0(
      sum(dentro_i1),
      " de ",
      nrow(i1_sf),
      " nodos dentro o sobre el borde"
    ),
    x = "X (m)",
    y = "Y (m)",
    shape = "Estado"
  ) +
  theme_minimal()

ggsave(
  filename = file.path(
    ruta_auditoria,
    "Coherencia_nodos_lote_I1.png"
  ),
  plot = g_coherencia,
  width = 8,
  height = 7,
  dpi = 400
)


## 8. Función para crear malla ---------------------------------

crear_malla <- function(
  poligono,
  numero_celdas_dimension_mayor = 120
) {

  bbox <- st_bbox(
    poligono
  )

  ancho <- as.numeric(
    bbox["xmax"] -
      bbox["xmin"]
  )

  alto <- as.numeric(
    bbox["ymax"] -
      bbox["ymin"]
  )

  dimension_mayor <- max(
    ancho,
    alto
  )

  cellsize <- dimension_mayor /
    numero_celdas_dimension_mayor

  centros <- st_make_grid(
    poligono,
    cellsize = cellsize,
    what = "centers",
    square = TRUE
  )

  grid <- st_sf(
    geometry = centros
  )

  grid <- grid[
    lengths(
      st_intersects(
        grid,
        poligono
      )
    ) > 0,
  ]

  grid$grid_id <- seq_len(
    nrow(grid)
  )

  list(
    grid = grid,
    cellsize = cellsize,
    ancho = ancho,
    alto = alto
  )
}


## 9. Malla exclusiva del lote I1 ------------------------------

malla_i1 <- crear_malla(
  lote_i1_union
)

grid_i1 <- malla_i1$grid

if (nrow(grid_i1) == 0) {
  stop(
    "La malla I1 quedó vacía."
  )
}


## 10. Función general de Kriging ------------------------------

ejecutar_kriging <- function(
  variable,
  modelo,
  puntos_sf,
  grid_sf
) {

  if (is.null(modelo)) {
    return(NULL)
  }

  if (!variable %in% names(puntos_sf)) {
    return(NULL)
  }

  fml <- as.formula(
    paste(
      variable,
      "~ 1"
    )
  )

  resultado <- try(
    krige(
      formula = fml,
      locations = puntos_sf,
      newdata = grid_sf,
      model = modelo,
      debug.level = 0
    ),
    silent = TRUE
  )

  if (inherits(resultado, "try-error")) {
    warning(
      "Falló el Kriging de ",
      variable,
      ": ",
      as.character(resultado)
    )
    return(NULL)
  }

  coords <- st_coordinates(
    resultado
  )

  if (nrow(resultado) != nrow(grid_sf)) {
    stop(
      "El número de predicciones no coincide con la malla para ",
      variable
    )
  }

  resultado |>
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
}


## 11. Kriging I1 con el polígono correcto ---------------------

superficies_i1 <- map(
  variables_i1,
  function(variable) {

    ejecutar_kriging(
      variable = variable,
      modelo = modelos_i1[[variable]],
      puntos_sf = i1_sf,
      grid_sf = grid_i1
    )
  }
)

names(
  superficies_i1
) <- variables_i1

superficies_i1 <- superficies_i1[
  !map_lgl(
    superficies_i1,
    is.null
  )
]

if (length(superficies_i1) == 0) {
  stop(
    "No se generó ninguna superficie I1."
  )
}


## 12. Exportar superficies I1 ---------------------------------

imap(
  superficies_i1,
  function(df, variable) {

    write_csv(
      df,
      file.path(
        ruta_superficies,
        paste0(
          "Superficie_I1_",
          variable,
          ".csv"
        )
      )
    )
  }
)


## 13. Preparar puntos y borde I1 ------------------------------

lote_i1_borde <- st_boundary(
  lote_i1_union
)

puntos_i1 <- i1_sf |>
  st_drop_geometry() |>
  select(
    Nodo,
    X,
    Y,
    all_of(
      variables_i1
    )
  ) |>
  pivot_longer(
    cols = all_of(
      variables_i1
    ),
    names_to = "Variable",
    values_to = "Observado"
  )


## 14. Función de mapa I1 --------------------------------------

crear_mapa_i1 <- function(
  df,
  variable
) {

  puntos_var <- puntos_i1 |>
    filter(
      .data$Variable == variable
    )

  ggplot() +
    geom_raster(
      data = df,
      aes(
        x = X,
        y = Y,
        fill = Prediccion
      )
    ) +
    geom_contour(
      data = df,
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
      data = lote_i1_borde,
      fill = NA,
      linewidth = 0.55
    ) +
    geom_point(
      data = puntos_var,
      aes(
        x = X,
        y = Y,
        shape = ifelse(
          Dentro_lote_I1,
          "Dentro",
          "Fuera"
        )
      ),
      size = 1.2
    ) +
    coord_sf(
      expand = FALSE
    ) +
    labs(
      title = paste0(
        variable,
        " – I1"
      ),
      x = "X (m)",
      y = "Y (m)",
      fill = "Predicción",
      shape = "Nodo"
    ) +
    theme_minimal()
}


# Incorporar estado espacial a puntos_i1.
puntos_i1 <- puntos_i1 |>
  left_join(
    tabla_coherencia_i1 |>
      select(
        Nodo,
        Dentro_lote_I1
      ),
    by = "Nodo"
  )

mapas_i1 <- imap(
  superficies_i1,
  crear_mapa_i1
)

panel_i1 <- wrap_plots(
  mapas_i1,
  ncol = 2
) +
  plot_annotation(
    title =
      "Kriging ordinario compacto de variables de suelo en I1",
    subtitle =
      "Superficies recortadas al polígono específico de I1"
  )

ggsave(
  filename = file.path(
    ruta_figuras,
    "Panel_kriging_compacto_I1.png"
  ),
  plot = panel_i1,
  width = 14,
  height = 16,
  dpi = 400
)


## 15. Resumen de superficies I1 -------------------------------

resumen_i1 <- imap_dfr(
  superficies_i1,
  function(df, variable) {

    tibble(
      Variable = variable,
      Instante = "I1",
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
      Varianza_media = mean(
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

resumen_i2 <- tabla_superficies_i2 |>
  transmute(
    Variable,
    Instante = "I2",
    Pred_min,
    Pred_media,
    Pred_mediana,
    Pred_max,
    Varianza_media,
    Error_estandar_medio
  )

tabla_temporal <- bind_rows(
  resumen_i1,
  resumen_i2
) |>
  arrange(
    Variable,
    Instante
  )


## 16. Intersección espacial I1–I2 -----------------------------

area_i1 <- as.numeric(
  st_area(
    lote_i1_union
  )
)

area_i2 <- as.numeric(
  st_area(
    lote_i2_union
  )
)

lote_comun <- suppressWarnings(
  st_intersection(
    lote_i1_union,
    lote_i2_union
  )
)

area_comun <- ifelse(
  length(lote_comun) > 0 &&
    !all(
      st_is_empty(
        lote_comun
      )
    ),
  as.numeric(
    st_area(
      lote_comun
    )
  ),
  0
)

porcentaje_i1_cubierto <- ifelse(
  area_i1 > 0,
  100 * area_comun / area_i1,
  NA_real_
)

porcentaje_i2_cubierto <- ifelse(
  area_i2 > 0,
  100 * area_comun / area_i2,
  NA_real_
)

tabla_superposicion_lotes <- tibble(
  Indicador = c(
    "Área lote I1 (m²)",
    "Área lote I2 (m²)",
    "Área común I1–I2 (m²)",
    "Porcentaje de I1 cubierto por área común",
    "Porcentaje de I2 cubierto por área común"
  ),
  Valor = c(
    area_i1,
    area_i2,
    area_comun,
    porcentaje_i1_cubierto,
    porcentaje_i2_cubierto
  )
)

print(
  tabla_superposicion_lotes
)


## 17. Figura de superposición de lotes ------------------------

g_superposicion <- ggplot() +
  geom_sf(
    data = lote_i1_union,
    fill = NA,
    linewidth = 0.8
  ) +
  geom_sf(
    data = lote_i2_union,
    fill = NA,
    linewidth = 0.8,
    linetype = "dashed"
  ) +
  {
    if (
      area_comun > 0
    ) {
      geom_sf(
        data = lote_comun,
        alpha = 0.25
      )
    }
  } +
  coord_sf(
    expand = TRUE
  ) +
  labs(
    title =
      "Superposición espacial de los polígonos I1 e I2",
    subtitle = paste0(
      "Área común: ",
      round(
        area_comun,
        2
      ),
      " m²"
    ),
    x = "X (m)",
    y = "Y (m)"
  ) +
  theme_minimal()

ggsave(
  filename = file.path(
    ruta_auditoria,
    "Superposicion_lotes_I1_I2.png"
  ),
  plot = g_superposicion,
  width = 8,
  height = 7,
  dpi = 400
)


## 18. Comparación temporal sobre malla común ------------------

comparaciones_temporales <- list()
mapas_diferencia <- list()

variables_comunes <- intersect(
  names(
    superficies_i1
  ),
  intersect(
    names(
      modelos_i1
    ),
    names(
      modelos_i2
    )
  )
)

if (area_comun > 0) {

  malla_comun <- crear_malla(
    lote_comun
  )

  grid_comun <- malla_comun$grid

  if (nrow(grid_comun) == 0) {
    warning(
      "El área común existe, pero la malla común quedó vacía."
    )
  } else {

    lote_comun_borde <- st_boundary(
      lote_comun
    )

    for (variable in variables_comunes) {

      superficie_i1_comun <- ejecutar_kriging(
        variable = variable,
        modelo = modelos_i1[[variable]],
        puntos_sf = i1_sf,
        grid_sf = grid_comun
      )

      superficie_i2_comun <- ejecutar_kriging(
        variable = variable,
        modelo = modelos_i2[[variable]],
        puntos_sf = i2_sf,
        grid_sf = grid_comun
      )

      if (
        is.null(
          superficie_i1_comun
        ) ||
          is.null(
            superficie_i2_comun
          )
      ) {
        next
      }

      comparacion <- superficie_i1_comun |>
        select(
          grid_id,
          X,
          Y,
          I1 = Prediccion,
          Var_I1 = Varianza_kriging
        ) |>
        left_join(
          superficie_i2_comun |>
            select(
              grid_id,
              I2 = Prediccion,
              Var_I2 = Varianza_kriging
            ),
          by = "grid_id"
        ) |>
        mutate(
          Diferencia_I2_I1 =
            I2 -
            I1,
          Cambio_relativo_porcentaje = ifelse(
            I1 != 0,
            100 *
              (
                I2 -
                  I1
              ) /
              abs(
                I1
              ),
            NA_real_
          ),
          Variable = variable
        )

      comparaciones_temporales[[variable]] <- comparacion

      mapa <- ggplot() +
        geom_raster(
          data = comparacion,
          aes(
            x = X,
            y = Y,
            fill = Diferencia_I2_I1
          )
        ) +
        geom_sf(
          data = lote_comun_borde,
          fill = NA,
          linewidth = 0.55
        ) +
        coord_sf(
          expand = FALSE
        ) +
        labs(
          title = paste0(
            variable,
            ": I2 − I1"
          ),
          subtitle =
            "Comparación restringida al área espacial común",
          x = "X (m)",
          y = "Y (m)",
          fill = "Diferencia"
        ) +
        theme_minimal()

      mapas_diferencia[[variable]] <- mapa

      write_csv(
        comparacion,
        file.path(
          ruta_superficies,
          paste0(
            "Comparacion_area_comun_I1_I2_",
            variable,
            ".csv"
          )
        )
      )
    }
  }

} else {

  warning(
    "Los polígonos I1 e I2 no presentan intersección espacial. ",
    "No se generarán mapas de diferencia."
  )

  grid_comun <- NULL
}


## 19. Panel de diferencias ------------------------------------

if (length(mapas_diferencia) > 0) {

  panel_diferencias <- wrap_plots(
    mapas_diferencia,
    ncol = 2
  ) +
    plot_annotation(
      title =
        "Cambios espaciales exploratorios entre I1 e I2",
      subtitle =
        "Predicciones recalculadas sobre una malla común"
    )

  ggsave(
    filename = file.path(
      ruta_figuras,
      "Panel_diferencias_exploratorias_I2_menos_I1.png"
    ),
    plot = panel_diferencias,
    width = 14,
    height = 16,
    dpi = 400
  )
}


## 20. Resumen de cambios --------------------------------------

if (length(comparaciones_temporales) > 0) {

  tabla_cambios <- imap_dfr(
    comparaciones_temporales,
    function(df, variable) {

      tibble(
        Variable = variable,
        n_celdas_area_comun = nrow(df),
        Diferencia_media_I2_I1 = mean(
          df$Diferencia_I2_I1,
          na.rm = TRUE
        ),
        Diferencia_mediana_I2_I1 = median(
          df$Diferencia_I2_I1,
          na.rm = TRUE
        ),
        Diferencia_minima = min(
          df$Diferencia_I2_I1,
          na.rm = TRUE
        ),
        Diferencia_maxima = max(
          df$Diferencia_I2_I1,
          na.rm = TRUE
        ),
        Correlacion_espacial_I1_I2 = suppressWarnings(
          cor(
            df$I1,
            df$I2,
            method = "spearman",
            use = "complete.obs"
          )
        )
      )
    }
  )

} else {

  tabla_cambios <- tibble(
    Variable = character(),
    n_celdas_area_comun = integer(),
    Diferencia_media_I2_I1 = numeric(),
    Diferencia_mediana_I2_I1 = numeric(),
    Diferencia_minima = numeric(),
    Diferencia_maxima = numeric(),
    Correlacion_espacial_I1_I2 = numeric()
  )
}


## 21. Exportar tablas -----------------------------------------

write_xlsx(
  list(
    Auditoria_nodos_lote_I1 =
      tabla_coherencia_i1,

    Resumen_coherencia_I1 =
      resumen_coherencia_i1,

    Superposicion_lotes =
      tabla_superposicion_lotes,

    Modelos_I1 =
      tabla_modelos_i1,

    Resumen_I1_I2 =
      tabla_temporal,

    Cambios_area_comun =
      tabla_cambios
  ),
  path = file.path(
    ruta_tablas,
    "Resultados_temporales_compactos_I1_I2.xlsx"
  )
)


## 22. Guardar objeto final ------------------------------------

saveRDS(
  list(
    i1_sf =
      i1_sf,

    i2_sf =
      i2_sf,

    lote_i1 =
      lote_i1_union,

    lote_i2 =
      lote_i2_union,

    lote_comun =
      lote_comun,

    grid_i1 =
      grid_i1,

    grid_comun =
      grid_comun,

    superficies_i1 =
      superficies_i1,

    comparaciones_temporales =
      comparaciones_temporales,

    resumen_temporal =
      tabla_temporal,

    cambios_temporales =
      tabla_cambios,

    auditoria_coherencia =
      tabla_coherencia_i1
  ),
  file = file.path(
    ruta_objetos,
    "objetos_temporales_I1_I2.rds"
  )
)


## 23. Registro de sesión --------------------------------------

capture.output(
  sessionInfo(),
  file = file.path(
    ruta_salida,
    "sessionInfo_I1_C.R.txt"
  )
)


## 24. Resumen en consola --------------------------------------

cat(
  "\n============================================\n"
)

cat(
  "SCRIPT I1-C REFORMULADO FINALIZADO\n\n"
)

cat(
  "Shapefile I1 utilizado:\n",
  ruta_lote_i1,
  "\n\n"
)

cat(
  "Nodos I1 dentro del lote:\n",
  sum(
    dentro_i1
  ),
  " de ",
  nrow(
    i1_sf
  ),
  "\n\n"
)

cat(
  "Área común I1–I2:\n",
  round(
    area_comun,
    2
  ),
  " m²\n\n"
)

cat(
  "Variables comparadas sobre área común:\n",
  paste(
    names(
      comparaciones_temporales
    ),
    collapse = ", "
  ),
  "\n\n"
)

cat(
  "Advertencia: la comparación I1–I2 es exploratoria y no ",
  "constituye una serie temporal formal.\n"
)

cat(
  "============================================\n"
)
