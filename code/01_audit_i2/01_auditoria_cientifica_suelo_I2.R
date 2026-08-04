
## ============================================================
## SCRIPT 01. AUDITORÍA CIENTÍFICA DE DATOS DE SUELO – I2
## Artículo cebada: revisión mayor
##
## Entrada:
##   AJUSTES_ARTÍCULO/DATOS BASE/Nodos_Lote_2_B.csv
##
## Salidas:
##   AJUSTES_ARTÍCULO/01_AUDITORIA_I2/
##
## Propósitos:
##   1. Verificar estructura, coordenadas, duplicados y valores faltantes.
##   2. Resumir distribución y variabilidad de las variables de suelo.
##   3. Incorporar CE corregida mediante la ecuación OLS completa (n = 28).
##   4. Examinar correlaciones sin usar p-valores como criterio central.
##   5. Preparar la base depurada para el análisis geoestadístico.
## ============================================================


## 0. Paquetes -------------------------------------------------

paquetes <- c(
  "readr",
  "dplyr",
  "tidyr",
  "purrr",
  "ggplot2",
  "e1071",
  "writexl",
  "patchwork"
)

instalar <- paquetes[
  !paquetes %in% rownames(installed.packages())
]

if (length(instalar) > 0) {
  install.packages(instalar)
}

library(readr)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(e1071)
library(writexl)
library(patchwork)


## 1. Rutas ----------------------------------------------------

ruta_base <- paste0(
  "D:/UDCA_2017_CLASE/46_ARTICULOS_2025/",
  "ARTICULO_CEBADA_DOCTORADO/AJUSTES_ARTÍCULO"
)

ruta_datos <- file.path(
  ruta_base,
  "DATOS BASE",
  "Nodos_Lote_2_B.csv"
)

ruta_salida <- file.path(
  ruta_base,
  "01_AUDITORIA_I2"
)

ruta_figuras <- file.path(
  ruta_salida,
  "figuras"
)

ruta_tablas <- file.path(
  ruta_salida,
  "tablas"
)

ruta_datos_procesados <- file.path(
  ruta_salida,
  "datos_procesados"
)

if (!dir.exists(ruta_base)) {
  stop(
    "No existe la carpeta base:\n",
    ruta_base
  )
}

if (!file.exists(ruta_datos)) {
  stop(
    "No existe el archivo:\n",
    ruta_datos,
    "\n\nArchivos disponibles en DATOS BASE:\n",
    paste(
      list.files(
        file.path(ruta_base, "DATOS BASE")
      ),
      collapse = "\n"
    )
  )
}

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
  ruta_datos_procesados,
  recursive = TRUE,
  showWarnings = FALSE
)


## 2. Lectura y normalización de nombres -----------------------

suelo_i2 <- read_csv(
  ruta_datos,
  show_col_types = FALSE,
  trim_ws = TRUE
)

# El archivo original contiene "N " con un espacio final.
names(suelo_i2) <- trimws(names(suelo_i2))

# Normalizar nombres esperados.
nombres_requeridos <- c(
  "Nodo",
  "X",
  "Y",
  "N",
  "P",
  "K",
  "pH",
  "T",
  "Hv",
  "CE"
)

faltantes_columnas <- setdiff(
  nombres_requeridos,
  names(suelo_i2)
)

if (length(faltantes_columnas) > 0) {
  stop(
    "Faltan las siguientes columnas:\n",
    paste(
      faltantes_columnas,
      collapse = ", "
    ),
    "\n\nColumnas encontradas:\n",
    paste(
      names(suelo_i2),
      collapse = ", "
    )
  )
}

suelo_i2 <- suelo_i2 |>
  select(
    all_of(nombres_requeridos)
  ) |>
  mutate(
    Nodo = as.integer(Nodo),
    across(
      c(X, Y, N, P, K, pH, T, Hv, CE),
      as.numeric
    )
  )


## 3. Incorporar CE corregida OLS ------------------------------

# Ecuación seleccionada usando las 28 muestras:
# CE_lab = 266.6932 + 0.9303 × CE_sensor

intercepto_ols <- 266.6932
pendiente_ols <- 0.9303

suelo_i2 <- suelo_i2 |>
  mutate(
    CE_corregida_OLS =
      intercepto_ols +
      pendiente_ols * CE
  )


## 4. Verificaciones estructurales -----------------------------

tabla_estructura <- tibble(
  Indicador = c(
    "Número de filas",
    "Número de columnas",
    "Nodos únicos",
    "Coordenadas XY únicas",
    "Duplicados de Nodo",
    "Duplicados exactos de coordenadas",
    "Valores faltantes totales"
  ),
  Valor = c(
    nrow(suelo_i2),
    ncol(suelo_i2),
    n_distinct(suelo_i2$Nodo),
    n_distinct(
      paste(
        suelo_i2$X,
        suelo_i2$Y,
        sep = "_"
      )
    ),
    sum(duplicated(suelo_i2$Nodo)),
    sum(
      duplicated(
        suelo_i2 |>
          select(X, Y)
      )
    ),
    sum(is.na(suelo_i2))
  )
)

print(tabla_estructura)


## 5. Valores faltantes por variable ---------------------------

tabla_faltantes <- tibble(
  Variable = names(suelo_i2),
  NA_n = map_int(
    suelo_i2,
    ~ sum(is.na(.x))
  ),
  NA_porcentaje = map_dbl(
    suelo_i2,
    ~ 100 * mean(is.na(.x))
  )
)

print(tabla_faltantes)


## 6. Duplicados y distancias entre puntos ---------------------

duplicados_nodo <- suelo_i2 |>
  filter(
    duplicated(Nodo) |
      duplicated(Nodo, fromLast = TRUE)
  )

duplicados_xy <- suelo_i2 |>
  group_by(X, Y) |>
  filter(n() > 1) |>
  ungroup()

coords <- suelo_i2 |>
  select(X, Y) |>
  as.matrix()

matriz_distancias <- as.matrix(
  dist(coords)
)

diag(matriz_distancias) <- NA_real_

distancia_vecino_mas_cercano <- apply(
  matriz_distancias,
  1,
  min,
  na.rm = TRUE
)

tabla_distancias <- tibble(
  Nodo = suelo_i2$Nodo,
  X = suelo_i2$X,
  Y = suelo_i2$Y,
  Distancia_vecino_mas_cercano_m =
    distancia_vecino_mas_cercano
)

resumen_distancias <- tibble(
  Indicador = c(
    "Distancia mínima entre nodos (m)",
    "Mediana distancia al vecino más cercano (m)",
    "Media distancia al vecino más cercano (m)",
    "Distancia máxima al vecino más cercano (m)"
  ),
  Valor = c(
    min(
      distancia_vecino_mas_cercano,
      na.rm = TRUE
    ),
    median(
      distancia_vecino_mas_cercano,
      na.rm = TRUE
    ),
    mean(
      distancia_vecino_mas_cercano,
      na.rm = TRUE
    ),
    max(
      distancia_vecino_mas_cercano,
      na.rm = TRUE
    )
  )
)

print(resumen_distancias)


## 7. Estadística descriptiva ----------------------------------

variables_suelo <- c(
  "N",
  "P",
  "K",
  "pH",
  "T",
  "Hv",
  "CE",
  "CE_corregida_OLS"
)

describir_variable <- function(x, nombre) {

  x_valido <- x[
    is.finite(x)
  ]

  media <- mean(x_valido)
  desviacion <- sd(x_valido)

  tibble(
    Variable = nombre,
    n = length(x_valido),
    Faltantes = sum(!is.finite(x)),
    Media = media,
    Mediana = median(x_valido),
    DE = desviacion,
    CV_porcentaje =
      ifelse(
        media == 0,
        NA_real_,
        100 * desviacion / abs(media)
      ),
    Minimo = min(x_valido),
    Q1 = quantile(
      x_valido,
      0.25,
      names = FALSE
    ),
    Q3 = quantile(
      x_valido,
      0.75,
      names = FALSE
    ),
    Maximo = max(x_valido),
    Rango = diff(range(x_valido)),
    Asimetria = e1071::skewness(
      x_valido,
      type = 2,
      na.rm = TRUE
    ),
    Curtosis_exceso = e1071::kurtosis(
      x_valido,
      type = 2,
      na.rm = TRUE
    )
  )
}

tabla_descriptiva <- map_dfr(
  variables_suelo,
  function(v) {
    describir_variable(
      suelo_i2[[v]],
      v
    )
  }
)

print(tabla_descriptiva)


## 8. Detección exploratoria de valores extremos ---------------

detectar_extremos_iqr <- function(datos, variable) {

  x <- datos[[variable]]

  q1 <- quantile(
    x,
    0.25,
    na.rm = TRUE
  )

  q3 <- quantile(
    x,
    0.75,
    na.rm = TRUE
  )

  iqr <- q3 - q1

  limite_inferior <- q1 - 1.5 * iqr
  limite_superior <- q3 + 1.5 * iqr

  datos |>
    transmute(
      Nodo,
      X,
      Y,
      Variable = variable,
      Valor = .data[[variable]],
      Limite_inferior = limite_inferior,
      Limite_superior = limite_superior,
      Extremo_IQR =
        Valor < limite_inferior |
        Valor > limite_superior
    ) |>
    filter(Extremo_IQR)
}

tabla_extremos <- map_dfr(
  variables_suelo,
  function(v) {
    detectar_extremos_iqr(
      suelo_i2,
      v
    )
  }
)

print(tabla_extremos)

# Nota:
# Una observación marcada mediante IQR NO se elimina automáticamente.
# Debe revisarse en contexto espacial y experimental.


## 9. Correlaciones --------------------------------------------

datos_cor <- suelo_i2 |>
  select(
    all_of(variables_suelo)
  )

cor_pearson <- cor(
  datos_cor,
  use = "pairwise.complete.obs",
  method = "pearson"
)

cor_spearman <- cor(
  datos_cor,
  use = "pairwise.complete.obs",
  method = "spearman"
)

cor_kendall <- cor(
  datos_cor,
  use = "pairwise.complete.obs",
  method = "kendall"
)

convertir_cor_tabla <- function(matriz, metodo) {

  as.data.frame(matriz) |>
    mutate(
      Variable_1 = rownames(matriz),
      .before = 1
    ) |>
    pivot_longer(
      cols = -Variable_1,
      names_to = "Variable_2",
      values_to = "Correlacion"
    ) |>
    mutate(
      Metodo = metodo,
      .before = 1
    )
}

tabla_correlaciones <- bind_rows(
  convertir_cor_tabla(
    cor_pearson,
    "Pearson"
  ),
  convertir_cor_tabla(
    cor_spearman,
    "Spearman"
  ),
  convertir_cor_tabla(
    cor_kendall,
    "Kendall"
  )
)

print(
  tabla_correlaciones |>
    filter(
      Variable_1 != Variable_2
    ) |>
    arrange(
      Metodo,
      desc(abs(Correlacion))
    )
)


## 10. Figuras de distribución --------------------------------

datos_largos <- suelo_i2 |>
  select(
    Nodo,
    all_of(variables_suelo)
  ) |>
  pivot_longer(
    cols = all_of(variables_suelo),
    names_to = "Variable",
    values_to = "Valor"
  )

g_histogramas <- ggplot(
  datos_largos,
  aes(
    x = Valor
  )
) +
  geom_histogram(
    bins = 8,
    color = "white"
  ) +
  facet_wrap(
    ~ Variable,
    scales = "free",
    ncol = 3
  ) +
  labs(
    title =
      "Distribución de las variables de suelo en I2",
    x = "Valor",
    y = "Frecuencia"
  ) +
  theme_minimal()

g_boxplots <- ggplot(
  datos_largos,
  aes(
    x = Variable,
    y = Valor
  )
) +
  geom_boxplot() +
  facet_wrap(
    ~ Variable,
    scales = "free",
    ncol = 3
  ) +
  labs(
    title =
      "Diagramas de caja de las variables de suelo en I2",
    x = NULL,
    y = "Valor"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

g_qq <- ggplot(
  datos_largos,
  aes(
    sample = Valor
  )
) +
  stat_qq() +
  stat_qq_line() +
  facet_wrap(
    ~ Variable,
    scales = "free",
    ncol = 3
  ) +
  labs(
    title =
      "Gráficos cuantiles normales de las variables de suelo en I2",
    x = "Cuantil teórico",
    y = "Cuantil observado"
  ) +
  theme_minimal()


## 11. Mapas exploratorios de puntos ---------------------------

g_puntos <- ggplot(
  datos_largos |>
    left_join(
      suelo_i2 |>
        select(
          Nodo,
          X,
          Y
        ),
      by = "Nodo"
    ),
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
    scales = "free",
    ncol = 3
  ) +
  coord_equal() +
  labs(
    title =
      "Distribución espacial de los valores observados en I2",
    x = "X (m)",
    y = "Y (m)",
    fill = "Valor"
  ) +
  theme_minimal()


## 12. Mapas de calor de correlación ---------------------------

preparar_heatmap <- function(matriz, metodo) {

  as.data.frame(matriz) |>
    mutate(
      Variable_1 = rownames(matriz),
      .before = 1
    ) |>
    pivot_longer(
      cols = -Variable_1,
      names_to = "Variable_2",
      values_to = "Correlacion"
    ) |>
    mutate(
      Metodo = metodo
    )
}

heatmap_datos <- bind_rows(
  preparar_heatmap(
    cor_pearson,
    "Pearson"
  ),
  preparar_heatmap(
    cor_spearman,
    "Spearman"
  )
)

g_correlaciones <- ggplot(
  heatmap_datos,
  aes(
    x = Variable_1,
    y = Variable_2,
    fill = Correlacion
  )
) +
  geom_tile() +
  geom_text(
    aes(
      label = sprintf(
        "%.2f",
        Correlacion
      )
    ),
    size = 3
  ) +
  facet_wrap(
    ~ Metodo
  ) +
  scale_fill_gradient2(
    limits = c(-1, 1),
    midpoint = 0
  ) +
  labs(
    title =
      "Correlaciones entre variables de suelo en I2",
    x = NULL,
    y = NULL,
    fill = "r / ρ"
  ) +
  theme_minimal() +
  theme(
    axis.text.x =
      element_text(
        angle = 45,
        hjust = 1
      )
  )


## 13. Comparación CE original y corregida ---------------------

g_ce_comparacion <- ggplot(
  suelo_i2,
  aes(
    x = CE,
    y = CE_corregida_OLS
  )
) +
  geom_point(
    size = 2.5
  ) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed"
  ) +
  labs(
    title =
      "CE original y CE corregida mediante OLS",
    subtitle =
      "CE corregida = 266.6932 + 0.9303 × CE original",
    x = "CE original del sensor (µS/cm)",
    y = "CE corregida OLS (µS/cm)"
  ) +
  theme_minimal()


## 14. Exportación de figuras ---------------------------------

ggsave(
  file.path(
    ruta_figuras,
    "01_histogramas_variables_I2.png"
  ),
  g_histogramas,
  width = 12,
  height = 9,
  dpi = 400
)

ggsave(
  file.path(
    ruta_figuras,
    "02_boxplots_variables_I2.png"
  ),
  g_boxplots,
  width = 12,
  height = 9,
  dpi = 400
)

ggsave(
  file.path(
    ruta_figuras,
    "03_QQplots_variables_I2.png"
  ),
  g_qq,
  width = 12,
  height = 9,
  dpi = 400
)

ggsave(
  file.path(
    ruta_figuras,
    "04_distribucion_espacial_puntos_I2.png"
  ),
  g_puntos,
  width = 12,
  height = 9,
  dpi = 400
)

ggsave(
  file.path(
    ruta_figuras,
    "05_correlaciones_I2.png"
  ),
  g_correlaciones,
  width = 13,
  height = 7,
  dpi = 400
)

ggsave(
  file.path(
    ruta_figuras,
    "06_CE_original_corregida_I2.png"
  ),
  g_ce_comparacion,
  width = 7,
  height = 6,
  dpi = 400
)


## 15. Exportación de tablas y base depurada -------------------

write_xlsx(
  list(
    Estructura =
      tabla_estructura,

    Valores_faltantes =
      tabla_faltantes,

    Resumen_distancias =
      resumen_distancias,

    Distancias_por_nodo =
      tabla_distancias,

    Descriptiva =
      tabla_descriptiva,

    Valores_extremos_IQR =
      tabla_extremos,

    Correlaciones =
      tabla_correlaciones,

    Duplicados_nodo =
      duplicados_nodo,

    Duplicados_XY =
      duplicados_xy
  ),
  path = file.path(
    ruta_tablas,
    "Auditoria_datos_suelo_I2.xlsx"
  )
)

write_csv(
  suelo_i2,
  file.path(
    ruta_datos_procesados,
    "Suelo_I2_auditado_CE_corregida.csv"
  )
)


## 16. Registro de sesión --------------------------------------

capture.output(
  sessionInfo(),
  file = file.path(
    ruta_salida,
    "sessionInfo_Script01.txt"
  )
)


## 17. Resumen final en consola --------------------------------

cat(
  "\n============================================\n"
)

cat(
  "SCRIPT 01 FINALIZADO\n"
)

cat(
  "Datos de entrada:\n",
  ruta_datos,
  "\n\n"
)

cat(
  "Número de nodos:",
  nrow(suelo_i2),
  "\n"
)

cat(
  "Base auditada guardada en:\n",
  file.path(
    ruta_datos_procesados,
    "Suelo_I2_auditado_CE_corregida.csv"
  ),
  "\n\n"
)

cat(
  "Resultados guardados en:\n",
  ruta_salida,
  "\n"
)

cat(
  "============================================\n"
)
