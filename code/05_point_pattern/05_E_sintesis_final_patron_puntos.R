
## ============================================================
## 05_PATRON_PUNTOS – SCRIPT 05-E
## SÍNTESIS FINAL, FIGURAS PARA PUBLICACIÓN Y TABLAS DE SOPORTE
## ============================================================
##
## Entradas:
##   05_PATRON_PUNTOS/01_DIAGNOSTICO_GLOBAL/objetos_R/
##     objetos_diagnostico_global_patron_puntos.rds
##
##   05_PATRON_PUNTOS/02_PATRONES_POR_GRUPO/objetos_R/
##     objetos_patrones_por_grupo.rds
##
##   05_PATRON_PUNTOS/03_ASOCIACION_SUELO_CULTIVO/objetos_R/
##     objetos_asociacion_espacial_suelo_cultivo.rds
##
##   05_PATRON_PUNTOS/04_HOTSPOTS_LISA/objetos_R/
##     objetos_hotspots_LISA_Gi.rds
##
## Salida:
##   05_PATRON_PUNTOS/05_SINTESIS_PUBLICACION/
##
## PROPÓSITO
##   1. Integrar los resultados de los Scripts 05-A a 05-D.
##   2. Generar tablas compactas para manuscrito y suplemento.
##   3. Seleccionar asociaciones suelo–cultivo prioritarias.
##   4. Resumir patrones globales, por dosis y locales.
##   5. Crear paneles finales de alta resolución.
##
## PRINCIPIOS DE INTERPRETACIÓN
##   - Los 56 puntos representan nodos de muestreo, no plantas
##     individuales exhaustivamente georreferenciadas.
##   - La regularidad global se atribuye primero al diseño de
##     muestreo y al arreglo experimental.
##   - Las variables de suelo fueron predichas mediante Kriging.
##   - N, P y K son índices instrumentales nominales.
##   - pH y Hv se mantienen como variables exploratorias.
##   - CE corregida OLS conserva el orden espacial de la CE original.
##   - La interpretación prioriza magnitud, estabilidad, dirección,
##     coherencia espacial y comportamiento residual.
## ============================================================


## 0. Paquetes -------------------------------------------------

paquetes <- c(
  "dplyr",
  "tidyr",
  "purrr",
  "stringr",
  "ggplot2",
  "readr",
  "writexl",
  "sf",
  "patchwork",
  "scales"
)

instalar <- paquetes[
  !paquetes %in% rownames(installed.packages())
]

if (length(instalar) > 0) {
  install.packages(instalar)
}

library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(ggplot2)
library(readr)
library(writexl)
library(sf)
library(patchwork)
library(scales)


## 1. Rutas ----------------------------------------------------

ruta_base <- paste0(
  "D:/UDCA_2017_CLASE/46_ARTICULOS_2025/",
  "ARTICULO_CEBADA_DOCTORADO/AJUSTES_ARTÍCULO"
)

ruta_A <- file.path(
  ruta_base,
  "05_PATRON_PUNTOS",
  "01_DIAGNOSTICO_GLOBAL",
  "objetos_R",
  "objetos_diagnostico_global_patron_puntos.rds"
)

ruta_B <- file.path(
  ruta_base,
  "05_PATRON_PUNTOS",
  "02_PATRONES_POR_GRUPO",
  "objetos_R",
  "objetos_patrones_por_grupo.rds"
)

ruta_C <- file.path(
  ruta_base,
  "05_PATRON_PUNTOS",
  "03_ASOCIACION_SUELO_CULTIVO",
  "objetos_R",
  "objetos_asociacion_espacial_suelo_cultivo.rds"
)

ruta_D <- file.path(
  ruta_base,
  "05_PATRON_PUNTOS",
  "04_HOTSPOTS_LISA",
  "objetos_R",
  "objetos_hotspots_LISA_Gi.rds"
)

rutas_entrada <- c(
  ruta_A,
  ruta_B,
  ruta_C,
  ruta_D
)

faltantes <- rutas_entrada[
  !file.exists(rutas_entrada)
]

if (length(faltantes) > 0) {
  stop(
    "Faltan los siguientes objetos:\n",
    paste(
      faltantes,
      collapse = "\n"
    )
  )
}

ruta_salida <- file.path(
  ruta_base,
  "05_PATRON_PUNTOS",
  "05_SINTESIS_PUBLICACION"
)

ruta_tablas <- file.path(
  ruta_salida,
  "tablas"
)

ruta_figuras <- file.path(
  ruta_salida,
  "figuras"
)

ruta_datos <- file.path(
  ruta_salida,
  "datos_procesados"
)

ruta_objetos <- file.path(
  ruta_salida,
  "objetos_R"
)

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
  ruta_datos,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  ruta_objetos,
  recursive = TRUE,
  showWarnings = FALSE
)


## 2. Cargar objetos -------------------------------------------

obj_A <- readRDS(ruta_A)
obj_B <- readRDS(ruta_B)
obj_C <- readRDS(ruta_C)
obj_D <- readRDS(ruta_D)

datos <- obj_C$datos
cultivo_sf <- obj_C$cultivo_sf
lote <- obj_C$lote

tabla_clark_global <- obj_A$tabla_clark_evans
tabla_envolventes_globales <- obj_A$tabla_resumen_envolventes
interpretacion_global <- obj_A$interpretacion_global

tabla_clark_grupos <- obj_B$tabla_clark
tabla_envolventes_grupos <- obj_B$resumen_envolventes
sintesis_grupos <- obj_B$sintesis_grupos
kernel_dosis <- obj_B$kernel_dosis

tabla_correlaciones <- obj_C$tabla_correlaciones
tabla_moran <- obj_C$tabla_moran
tabla_asociacion <- obj_C$tabla_sintesis
asociacion_distancia <- obj_C$asociacion_distancia
comparacion_ce <- obj_C$comparacion_ce

tabla_lisa <- obj_D$tabla_lisa
tabla_gi <- obj_D$tabla_gi
resumen_lisa <- obj_D$resumen_lisa
resumen_gi <- obj_D$resumen_gi
resumen_estabilidad <- obj_D$resumen_estabilidad
sintesis_local <- obj_D$tabla_sintesis


## 3. Catálogo y etiquetas -------------------------------------

variables_cultivo <- c(
  "Macollas",
  "Altura",
  "Diametro",
  "NDVI",
  "Nitratos",
  "Potasio_savia",
  "Calcio_savia"
)

variables_suelo <- c(
  "Suelo_CE",
  "Suelo_CE_corregida_OLS",
  "Suelo_N",
  "Suelo_P",
  "Suelo_K",
  "Suelo_pH",
  "Suelo_T",
  "Suelo_Hv"
)

etiquetas <- c(
  Macollas = "Tillers",
  Altura = "Plant height",
  Diametro = "Stem diameter",
  NDVI = "NDVI",
  Nitratos = "Sap nitrate",
  Potasio_savia = "Sap potassium",
  Calcio_savia = "Sap calcium",
  Suelo_CE = "Soil EC",
  Suelo_CE_corregida_OLS = "OLS-corrected soil EC",
  Suelo_N = "Soil N index",
  Suelo_P = "Soil P index",
  Suelo_K = "Soil K index",
  Suelo_pH = "Soil pH",
  Suelo_T = "Soil temperature",
  Suelo_Hv = "Soil moisture index"
)

etiquetar <- function(x) {

  salida <- unname(
    etiquetas[
      x
    ]
  )

  salida[
    is.na(
      salida
    )
  ] <- x[
    is.na(
      salida
    )
  ]

  salida
}


## 4. Síntesis del patrón global -------------------------------

tabla_patron_global <- tabla_clark_global |>
  transmute(
    Nivel = "Patrón global",
    Grupo = "56 nodos",
    n,
    Area_m2,
    Intensidad_nodos_m2,
    R_Clark_Evans,
    Media_NN_observada_m,
    Sim_Q025_m,
    Sim_Q975_m,
    Clasificacion
  )

tabla_KL_global <- tabla_envolventes_globales |>
  filter(
    Funcion %in% c(
      "K",
      "L"
    )
  ) |>
  select(
    Funcion,
    Porcentaje_dentro,
    Porcentaje_sobre,
    Porcentaje_bajo,
    r_primera_salida_m,
    Direccion_principal
  )


## 5. Síntesis por dosis y repetición --------------------------

tabla_patron_grupos <- tabla_clark_grupos |>
  left_join(
    tabla_envolventes_grupos |>
      filter(
        Funcion == "L"
      ) |>
      select(
        Tipo_grupo,
        Grupo,
        L_dentro = Porcentaje_dentro,
        L_sobre = Porcentaje_sobre,
        L_bajo = Porcentaje_bajo,
        L_interpretacion = Interpretacion
      ),
    by = c(
      "Tipo_grupo",
      "Grupo"
    )
  ) |>
  arrange(
    Tipo_grupo,
    Grupo
  )


## 6. Selección de asociaciones prioritarias -------------------

tabla_asociaciones_prioritarias <- tabla_asociacion |>
  mutate(
    Suelo_etiqueta = etiquetar(
      Variable_suelo
    ),
    Cultivo_etiqueta = etiquetar(
      Variable_cultivo
    ),
    Magnitud_Spearman_residual = abs(
      Spearman_residual
    ),
    Magnitud_Moran_residual = abs(
      Moran_bivariado_residual
    ),
    Puntuacion_integrada =
      0.60 *
        Magnitud_Spearman_residual +
      0.40 *
        Magnitud_Moran_residual
  ) |>
  arrange(
    desc(
      Puntuacion_integrada
    )
  ) |>
  mutate(
    Rango = row_number()
  )

top_asociaciones <- tabla_asociaciones_prioritarias |>
  slice_head(
    n = min(
      12,
      nrow(
        tabla_asociaciones_prioritarias
      )
    )
  )


## 7. Síntesis de hotspots y LISA ------------------------------

tabla_local_residual <- sintesis_local |>
  mutate(
    Variable_etiqueta = etiquetar(
      Variable
    ),
    Grupo = case_when(
      Variable %in% variables_cultivo ~
        "Cultivo",

      Variable %in% variables_suelo ~
        "Suelo",

      TRUE ~
        "Otro"
    ),
    Nodos_locales_residuales =
      Nodos_LISA_clasificados +
      Nodos_Gi_clasificados
  ) |>
  arrange(
    desc(
      Nodos_locales_residuales
    )
  )

tabla_local_prioritaria <- tabla_local_residual |>
  filter(
    Nodos_locales_residuales > 0
  )


## 8. Tabla integrada para manuscrito --------------------------

tabla_manuscrito <- top_asociaciones |>
  select(
    Rango,
    Variable_suelo,
    Suelo_etiqueta,
    Variable_cultivo,
    Cultivo_etiqueta,
    Spearman_original,
    Spearman_residual,
    Moran_bivariado_residual,
    Clasificacion_Moran,
    Estabilidad_direccion,
    Prioridad_interpretacion
  ) |>
  left_join(
    tabla_local_residual |>
      select(
        Variable,
        Nodos_LISA_clasificados,
        Nodos_Gi_clasificados,
        Evidencia_local_residual
      ),
    by = c(
      "Variable_cultivo" =
        "Variable"
    )
  ) |>
  rename(
    Cultivo_LISA_residual =
      Nodos_LISA_clasificados,
    Cultivo_Gi_residual =
      Nodos_Gi_clasificados,
    Cultivo_evidencia_local =
      Evidencia_local_residual
  )


## 9. Figura 1: patrón global y dosis ---------------------------

lote_sf <- st_as_sf(
  lote
)

lote_sf <- st_transform(
  lote_sf,
  st_crs(
    cultivo_sf
  )
)

lote_borde <- st_boundary(
  st_union(
    lote_sf
  )
)

datos_figura <- datos |>
  mutate(
    Dosis = case_when(
      as.character(
        Dose_percent
      ) %in% c(
        "100",
        "D100"
      ) ~ "100%",

      as.character(
        Dose_percent
      ) %in% c(
        "75",
        "D75"
      ) ~ "75%",

      TRUE ~ as.character(
        Dose_percent
      )
    )
  )

g_patron <- ggplot() +
  geom_sf(
    data = lote_borde,
    inherit.aes = FALSE,
    fill = NA,
    linewidth = 0.55
  ) +
  geom_point(
    data = datos_figura,
    aes(
      x = X,
      y = Y,
      shape = Dosis
    ),
    size = 2
  ) +
  coord_sf(
    expand = FALSE
  ) +
  labs(
    title = "Sampling-node pattern by fertilizer dose",
    subtitle = paste0(
      "Clark–Evans R = ",
      round(
        tabla_clark_global$R_Clark_Evans[1],
        2
      ),
      "; ",
      tabla_clark_global$Clasificacion[1]
    ),
    x = "X (m)",
    y = "Y (m)",
    shape = "Dose"
  ) +
  theme_minimal()

g_kernel <- ggplot(
  kernel_dosis,
  aes(
    x = X,
    y = Y,
    fill = Densidad_nodos_m2
  )
) +
  geom_raster() +
  geom_sf(
    data = lote_borde,
    inherit.aes = FALSE,
    fill = NA,
    linewidth = 0.45
  ) +
  facet_wrap(
    ~ Grupo,
    ncol = 2
  ) +
  coord_sf(
    expand = FALSE
  ) +
  scale_x_continuous(
    breaks = scales::breaks_pretty(n = 4),
    labels = function(x) {
      paste0(
        formatC(
          abs(x),
          format = "f",
          digits = 5
        ),
        "°W"
      )
    },
    guide = guide_axis(
      angle = 45,
      check.overlap = TRUE
    )
  ) +
  labs(
    title = "Kernel density by fertilizer dose",
    x = "X (m)",
    y = "Y (m)",
    fill = "Nodes/m²"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(
      size = 6.5,
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    panel.spacing.x = grid::unit(
      1.2,
      "lines"
    ),
    plot.margin = margin(
      6,
      10,
      16,
      6
    )
  )

panel_patron <- g_patron +
  g_kernel +
  plot_annotation(
    title = "Spatial configuration of the experimental sampling nodes"
  )

ggsave(
  filename = file.path(
    ruta_figuras,
    "Figura_05E_01_patron_global_dosis.png"
  ),
  plot = panel_patron,
  width = 15,
  height = 7.5,
  dpi = 500
)


## 10. Figura 2: correlación y Moran bivariado -----------------

heat_cor <- tabla_correlaciones |>
  mutate(
    Suelo = etiquetar(
      Variable_suelo
    ),
    Cultivo = etiquetar(
      Variable_cultivo
    )
  )

g_spearman <- ggplot(
  heat_cor,
  aes(
    x = Cultivo,
    y = Suelo,
    fill = Spearman_residual
  )
) +
  geom_tile() +
  geom_text(
    aes(
      label = sprintf(
        "%.2f",
        Spearman_residual
      )
    ),
    size = 2.6
  ) +
  labs(
    title = "Residual Spearman association",
    subtitle = "Adjusted for replication, dose, and omission treatment",
    x = "Crop variable",
    y = "Soil variable",
    fill = "Spearman ρ"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

heat_moran <- tabla_moran |>
  filter(
    Escenario == "Residual_ajustado"
  ) |>
  mutate(
    Suelo = etiquetar(
      Variable_suelo
    ),
    Cultivo = etiquetar(
      Variable_cultivo
    )
  )

g_moran <- ggplot(
  heat_moran,
  aes(
    x = Cultivo,
    y = Suelo,
    fill = Moran_bivariado
  )
) +
  geom_tile() +
  geom_text(
    aes(
      label = sprintf(
        "%.2f",
        Moran_bivariado
      )
    ),
    size = 2.6
  ) +
  labs(
    title = "Residual bivariate Moran association",
    subtitle = "Soil value at node i versus spatial lag of crop response",
    x = "Crop variable",
    y = "Soil variable",
    fill = "Bivariate I"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

panel_asociaciones <- g_spearman +
  g_moran +
  plot_annotation(
    title = "Residual spatial relationships between soil and crop variables"
  )

ggsave(
  filename = file.path(
    ruta_figuras,
    "Figura_05E_02_asociaciones_residuales.png"
  ),
  plot = panel_asociaciones,
  width = 16,
  height = 8,
  dpi = 500
)


## 11. Figura 3: asociaciones prioritarias ---------------------

top_plot <- top_asociaciones |>
  mutate(
    Par = paste0(
      Cultivo_etiqueta,
      " ~ ",
      Suelo_etiqueta
    ),
    Par = factor(
      Par,
      levels = rev(
        unique(
          Par
        )
      )
    )
  )

g_top <- ggplot(
  top_plot,
  aes(
    x = Spearman_residual,
    y = Par
  )
) +
  geom_vline(
    xintercept = 0,
    linetype = "dotted"
  ) +
  geom_point(
    aes(
      size = Magnitud_Moran_residual,
      shape = Coherencia_correlacion_espacial
    )
  ) +
  labs(
    title = "Priority residual soil–crop associations",
    subtitle = "Point size represents the magnitude of residual bivariate Moran I",
    x = "Residual Spearman ρ",
    y = NULL,
    size = "|Bivariate I|",
    shape = "Directional coherence"
  ) +
  theme_minimal()

ggsave(
  filename = file.path(
    ruta_figuras,
    "Figura_05E_03_asociaciones_prioritarias.png"
  ),
  plot = g_top,
  width = 11,
  height = 8,
  dpi = 500
)


## 12. Figura 4: estructura local residual ---------------------

local_plot <- tabla_local_residual |>
  mutate(
    Variable_etiqueta = factor(
      Variable_etiqueta,
      levels = rev(
        Variable_etiqueta[
          order(
            Nodos_locales_residuales
          )
        ]
      )
    )
  ) |>
  select(
    Variable_etiqueta,
    Grupo,
    Nodos_LISA_clasificados,
    Nodos_Gi_clasificados
  ) |>
  pivot_longer(
    cols = c(
      Nodos_LISA_clasificados,
      Nodos_Gi_clasificados
    ),
    names_to = "Metodo",
    values_to = "Nodos"
  ) |>
  mutate(
    Metodo = dplyr::recode(
      Metodo,
      "Nodos_LISA_clasificados" = "LISA",
      "Nodos_Gi_clasificados"   = "Getis–Ord Gi*"
    )
  )

g_local <- ggplot(
  local_plot,
  aes(
    x = Nodos,
    y = Variable_etiqueta,
    shape = Metodo
  )
) +
  geom_point(
    size = 2.5,
    position = position_dodge(
      width = 0.45
    )
  ) +
  facet_wrap(
    ~ Grupo,
    scales = "free_y",
    ncol = 2
  ) +
  labs(
    title = "Residual local spatial structure",
    subtitle = "Number of nodes classified outside the permutation reference",
    x = "Classified nodes",
    y = NULL,
    shape = "Method"
  ) +
  theme_minimal()

ggsave(
  filename = file.path(
    ruta_figuras,
    "Figura_05E_04_estructura_local_residual.png"
  ),
  plot = g_local,
  width = 12,
  height = 8,
  dpi = 500
)


## 13. Figura 5: NDVI como caso prioritario --------------------

variables_ndvi <- c(
  "NDVI",
  "Suelo_CE",
  "Suelo_CE_corregida_OLS",
  "Suelo_N",
  "Suelo_P",
  "Suelo_K"
)

lisa_ndvi <- tabla_lisa |>
  filter(
    Escenario == "Residual_ajustado",
    Variable %in% variables_ndvi
  ) |>
  mutate(
    Variable_etiqueta = etiquetar(
      Variable
    )
  )

gi_ndvi <- tabla_gi |>
  filter(
    Escenario == "Residual_ajustado",
    Variable %in% variables_ndvi
  ) |>
  mutate(
    Variable_etiqueta = etiquetar(
      Variable
    )
  )

g_lisa_ndvi <- ggplot() +
  geom_sf(
    data = lote_borde,
    inherit.aes = FALSE,
    fill = NA,
    linewidth = 0.45
  ) +
  geom_point(
    data = lisa_ndvi,
    aes(
      x = X,
      y = Y,
      shape = Clase_LISA
    ),
    size = 1.8
  ) +
  facet_wrap(
    ~ Variable_etiqueta,
    ncol = 3
  ) +
  coord_sf(
    expand = FALSE
  ) +
  scale_x_continuous(
    breaks = scales::breaks_pretty(n = 4),
    guide = guide_axis(
      angle = 45,
      check.overlap = TRUE
    )
  ) +
  labs(
    title = "Residual LISA patterns",
    x = "X (m)",
    y = "Y (m)",
    shape = "LISA class"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(
      size = 6.5,
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    panel.spacing.x = grid::unit(
      1.0,
      "lines"
    ),
    plot.margin = margin(
      6,
      10,
      14,
      6
    )
  )

g_gi_ndvi <- ggplot() +
  geom_sf(
    data = lote_borde,
    inherit.aes = FALSE,
    fill = NA,
    linewidth = 0.45
  ) +
  geom_point(
    data = gi_ndvi,
    aes(
      x = X,
      y = Y,
      shape = Clase_Gi
    ),
    size = 1.8
  ) +
  facet_wrap(
    ~ Variable_etiqueta,
    ncol = 3
  ) +
  coord_sf(
    expand = FALSE
  ) +
  scale_x_continuous(
    breaks = scales::breaks_pretty(n = 4),
    guide = guide_axis(
      angle = 45,
      check.overlap = TRUE
    )
  ) +
  labs(
    title = "Residual hotspots and coldspots",
    x = "X (m)",
    y = "Y (m)",
    shape = "Gi* class"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(
      size = 6.5,
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    panel.spacing.x = grid::unit(
      1.0,
      "lines"
    ),
    plot.margin = margin(
      6,
      10,
      14,
      6
    )
  )

panel_ndvi <- g_lisa_ndvi /
  g_gi_ndvi +
  plot_annotation(
    title = "Local residual patterns of NDVI and associated soil variables"
  )

ggsave(
  filename = file.path(
    ruta_figuras,
    "Figura_05E_05_NDVI_patrones_locales.png"
  ),
  plot = panel_ndvi,
  width = 17,
  height = 14,
  dpi = 500
)


## 14. Figura general de síntesis ------------------------------

panel_general <- (
  panel_patron /
    panel_asociaciones
) /
  (
    g_top +
      g_local
  ) +
  plot_annotation(
    title = "Integrated point-pattern and soil–crop spatial analysis"
  )

ggsave(
  filename = file.path(
    ruta_figuras,
    "Figura_05E_06_panel_general_sintesis.png"
  ),
  plot = panel_general,
  width = 18,
  height = 22,
  dpi = 500
)


## 15. Redacción automática prudente ---------------------------

hallazgo_global <- paste0(
  "The 56 sampling nodes exhibited a regular global pattern ",
  "(Clark–Evans R = ",
  round(
    tabla_clark_global$R_Clark_Evans[1],
    2
  ),
  "), consistent with the experimental sampling arrangement. ",
  "Ripley’s K and Besag’s L supported regularity at part of the ",
  "evaluated distance range."
)

hallazgo_dosis <- sintesis_grupos |>
  filter(
    Tipo_grupo == "Dosis"
  ) |>
  mutate(
    Texto = paste0(
      Grupo,
      ": ",
      Clasificacion,
      "; ",
      L_Interpretacion
    )
  ) |>
  pull(
    Texto
  ) |>
  paste(
    collapse = " | "
  )

hallazgo_asociacion <- top_asociaciones |>
  slice_head(
    n = min(
      5,
      nrow(
        top_asociaciones
      )
    )
  ) |>
  mutate(
    Texto = paste0(
      Cultivo_etiqueta,
      " versus ",
      Suelo_etiqueta,
      ": residual Spearman = ",
      round(
        Spearman_residual,
        2
      ),
      ", residual bivariate Moran I = ",
      round(
        Moran_bivariado_residual,
        2
      )
    )
  ) |>
  pull(
    Texto
  ) |>
  paste(
    collapse = " | "
  )

hallazgo_local <- tabla_local_residual |>
  slice_head(
    n = min(
      5,
      nrow(
        tabla_local_residual
      )
    )
  ) |>
  mutate(
    Texto = paste0(
      Variable_etiqueta,
      ": ",
      Nodos_LISA_clasificados,
      " LISA nodes and ",
      Nodos_Gi_clasificados,
      " Gi* nodes after adjustment"
    )
  ) |>
  pull(
    Texto
  ) |>
  paste(
    collapse = " | "
  )

tabla_textos <- tibble(
  Seccion = c(
    "Global point pattern",
    "Patterns by fertilizer dose",
    "Priority soil–crop associations",
    "Residual local structure",
    "Methodological caution"
  ),
  Texto = c(
    hallazgo_global,
    hallazgo_dosis,
    hallazgo_asociacion,
    hallazgo_local,
    paste0(
      "All point-pattern results refer to the sampled nodes rather ",
      "than an exhaustive census of individual plants. Soil variables ",
      "were predicted by Kriging, and N, P, and K should be interpreted ",
      "as nominal instrumental indices."
    )
  )
)


## 16. Exportar datos procesados -------------------------------

write_csv(
  tabla_patron_global,
  file.path(
    ruta_datos,
    "Sintesis_patron_global.csv"
  )
)

write_csv(
  tabla_patron_grupos,
  file.path(
    ruta_datos,
    "Sintesis_patrones_por_grupo.csv"
  )
)

write_csv(
  tabla_asociaciones_prioritarias,
  file.path(
    ruta_datos,
    "Asociaciones_suelo_cultivo_priorizadas.csv"
  )
)

write_csv(
  tabla_local_residual,
  file.path(
    ruta_datos,
    "Sintesis_estructura_local_residual.csv"
  )
)

write_csv(
  tabla_manuscrito,
  file.path(
    ruta_datos,
    "Tabla_compacta_para_manuscrito.csv"
  )
)

write_csv(
  tabla_textos,
  file.path(
    ruta_datos,
    "Textos_base_resultados_discusion.csv"
  )
)


## 17. Exportar libro de síntesis ------------------------------

write_xlsx(
  list(
    Resumen_ejecutivo =
      tabla_textos,

    Patron_global =
      tabla_patron_global,

    K_L_global =
      tabla_KL_global,

    Patron_por_grupo =
      tabla_patron_grupos,

    Asociaciones_prioritarias =
      tabla_asociaciones_prioritarias,

    Tabla_manuscrito =
      tabla_manuscrito,

    Estructura_local =
      tabla_local_residual,

    LISA_resumen =
      resumen_lisa,

    Gi_resumen =
      resumen_gi,

    Comparacion_CE =
      comparacion_ce,

    Asociacion_distancia =
      asociacion_distancia,

    Advertencias =
      interpretacion_global
  ),
  path = file.path(
    ruta_tablas,
    "Sintesis_final_patron_puntos_publicacion.xlsx"
  )
)


## 18. Guardar objeto final ------------------------------------

saveRDS(
  list(
    datos = datos,
    cultivo_sf = cultivo_sf,
    lote = lote,
    tabla_patron_global = tabla_patron_global,
    tabla_patron_grupos = tabla_patron_grupos,
    tabla_asociaciones_prioritarias = tabla_asociaciones_prioritarias,
    tabla_local_residual = tabla_local_residual,
    tabla_manuscrito = tabla_manuscrito,
    tabla_textos = tabla_textos
  ),
  file = file.path(
    ruta_objetos,
    "objetos_sintesis_final_patron_puntos.rds"
  )
)


## 19. Registro de sesión --------------------------------------

capture.output(
  sessionInfo(),
  file = file.path(
    ruta_salida,
    "sessionInfo_05_E.txt"
  )
)


## 20. Resumen en consola --------------------------------------

cat(
  "\n============================================\n"
)

cat(
  "SCRIPT 05-E FINALIZADO\n\n"
)

cat(
  "Patrón global:\n"
)

print(
  tabla_patron_global |>
    select(
      n,
      R_Clark_Evans,
      Clasificacion
    )
)

cat(
  "\nCinco asociaciones prioritarias:\n"
)

print(
  top_asociaciones |>
    select(
      Variable_suelo,
      Variable_cultivo,
      Spearman_residual,
      Moran_bivariado_residual,
      Prioridad_interpretacion
    ) |>
    slice_head(
      n = min(
        5,
        nrow(
          top_asociaciones
        )
      )
    )
)

cat(
  "\nVariables con mayor estructura local residual:\n"
)

print(
  tabla_local_residual |>
    select(
      Variable,
      Nodos_LISA_clasificados,
      Nodos_Gi_clasificados,
      Evidencia_local_residual
    ) |>
    slice_head(
      n = min(
        8,
        nrow(
          tabla_local_residual
        )
      )
    )
)

cat(
  "\nResultados guardados en:\n",
  ruta_salida,
  "\n"
)

cat(
  "============================================\n"
)
