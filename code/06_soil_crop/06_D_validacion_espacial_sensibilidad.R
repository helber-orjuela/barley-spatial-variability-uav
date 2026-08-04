
## ============================================================
## 06_SUELO_CULTIVO – SCRIPT 06-D
## VALIDACIÓN ESPACIAL, INFLUENCIA Y SENSIBILIDAD DE MODELOS
## ============================================================
##
## Entrada:
##   06_SUELO_CULTIVO/03_MODELOS_INTEGRADOS/objetos_R/
##   objetos_modelos_integrados_06C.rds
##
## Salida:
##   06_SUELO_CULTIVO/04_VALIDACION_SENSIBILIDAD/
##
## PROPÓSITO
##   1. Evaluar autocorrelación espacial de residuos LOOCV.
##   2. Comparar M1 frente al mejor M2 de cada respuesta.
##   3. Validar dejando una repetición completa por fuera.
##   4. Identificar nodos influyentes y estabilidad de coeficientes.
##   5. Comparar CE original y CE corregida OLS.
##   6. Determinar si las mejoras predictivas son robustas.
##
## MÉTODOS
##   - Moran global de residuos con k vecinos más cercanos.
##   - Semivariograma experimental de residuos.
##   - Leave-one-replication-out cross-validation.
##   - Distancia de Cook y leverage del ajuste completo.
##   - Sensibilidad CE original frente a CE corregida.
##
## PRINCIPIOS DE INTERPRETACIÓN
##   - Se priorizan magnitud del error, estabilidad y estructura
##     residual; no se usa el p-valor como criterio central.
##   - CE original y corregida pueden producir resultados idénticos
##     por tratarse de una transformación lineal.
##   - Las variables de suelo proceden de Kriging.
##   - N, P y K son índices instrumentales nominales.
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
  "patchwork",
  "sf",
  "spdep",
  "gstat",
  "tibble"
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
library(patchwork)
library(sf)
library(spdep)
library(gstat)
library(tibble)


## 1. Parámetros ------------------------------------------------

set.seed(2026)

k_vecinos <- 4
nsim_moran <- 999
epsilon <- 1e-12
umbral_cook <- NULL
umbral_leverage_factor <- 2


## 2. Rutas ----------------------------------------------------

ruta_base <- paste0(
  "D:/UDCA_2017_CLASE/46_ARTICULOS_2025/",
  "ARTICULO_CEBADA_DOCTORADO/AJUSTES_ARTÍCULO"
)

ruta_objeto <- file.path(
  ruta_base,
  "06_SUELO_CULTIVO",
  "03_MODELOS_INTEGRADOS",
  "objetos_R",
  "objetos_modelos_integrados_06C.rds"
)

ruta_salida <- file.path(
  ruta_base,
  "06_SUELO_CULTIVO",
  "04_VALIDACION_SENSIBILIDAD"
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

if (!file.exists(ruta_objeto)) {
  stop(
    "No existe el objeto del Script 06-C:\n",
    ruta_objeto
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


## 3. Cargar objeto --------------------------------------------

obj <- readRDS(
  ruta_objeto
)

base <- as_tibble(
  obj$base
)

biblioteca_modelos <- obj$biblioteca_modelos
objetos_modelos <- obj$objetos_modelos
tabla_comparacion <- obj$tabla_comparacion
tabla_predicciones <- obj$tabla_predicciones
mejores_modelos <- obj$mejores_modelos
sensibilidad_ce_06C <- obj$sensibilidad_ce

columnas_requeridas <- c(
  "Nodo",
  "X",
  "Y",
  "Replication"
)

faltantes <- setdiff(
  columnas_requeridas,
  names(base)
)

if (length(faltantes) > 0) {
  stop(
    "Faltan columnas requeridas: ",
    paste(
      faltantes,
      collapse = ", "
    )
  )
}

base <- base |>
  mutate(
    Replication = factor(
      Replication
    )
  )


## 4. Seleccionar modelos de validación ------------------------

modelos_validacion <- bind_rows(
  tabla_comparacion |>
    filter(
      Modelo_ID == "M1_Tratamientos"
    ) |>
    select(
      Respuesta,
      Modelo_ID,
      Formula
    ),
  mejores_modelos |>
    select(
      Respuesta,
      Modelo_ID,
      Formula
    )
) |>
  distinct(
    Respuesta,
    Modelo_ID,
    .keep_all = TRUE
  ) |>
  mutate(
    Tipo_modelo = case_when(
      Modelo_ID == "M1_Tratamientos" ~
        "M1",
      TRUE ~
        "Mejor_M2"
    )
  )


## 5. Funciones auxiliares -------------------------------------

correlacion_segura <- function(
  x,
  y,
  metodo = "pearson"
) {

  validos <- is.finite(x) &
    is.finite(y)

  if (sum(validos) < 4) {
    return(
      NA_real_
    )
  }

  if (
    sd(x[validos]) <= epsilon |
      sd(y[validos]) <= epsilon
  ) {
    return(
      NA_real_
    )
  }

  suppressWarnings(
    cor(
      x[validos],
      y[validos],
      method = metodo
    )
  )
}


calcular_metricas <- function(
  observado,
  predicho
) {

  validos <- is.finite(observado) &
    is.finite(predicho)

  observado <- observado[validos]
  predicho <- predicho[validos]

  if (length(observado) < 4) {
    return(
      tibble(
        n = length(observado),
        RMSE = NA_real_,
        MAE = NA_real_,
        Sesgo = NA_real_,
        Cor_Pearson = NA_real_,
        Cor_Spearman = NA_real_
      )
    )
  }

  error <- predicho - observado

  tibble(
    n = length(observado),
    RMSE = sqrt(
      mean(
        error^2
      )
    ),
    MAE = mean(
      abs(
        error
      )
    ),
    Sesgo = mean(
      error
    ),
    Cor_Pearson = correlacion_segura(
      observado,
      predicho,
      "pearson"
    ),
    Cor_Spearman = correlacion_segura(
      observado,
      predicho,
      "spearman"
    )
  )
}


ajustar_y_predecir <- function(
  formula,
  entrenamiento,
  prueba
) {

  modelo <- try(
    lm(
      formula,
      data = entrenamiento,
      na.action = na.exclude
    ),
    silent = TRUE
  )

  if (inherits(modelo, "try-error")) {
    return(
      rep(
        NA_real_,
        nrow(prueba)
      )
    )
  }

  prediccion <- try(
    predict(
      modelo,
      newdata = prueba
    ),
    silent = TRUE
  )

  if (inherits(prediccion, "try-error")) {
    return(
      rep(
        NA_real_,
        nrow(prueba)
      )
    )
  }

  as.numeric(
    prediccion
  )
}


extraer_residuos_loocv <- function(
  respuesta,
  modelo_id
) {

  tabla_predicciones |>
    filter(
      Respuesta == respuesta,
      Modelo_ID == modelo_id
    ) |>
    select(
      Nodo,
      Respuesta,
      Modelo_ID,
      Observado,
      Predicho_LOOCV,
      Residuo_LOOCV
    ) |>
    left_join(
      base |>
        select(
          Nodo,
          X,
          Y,
          Replication
        ),
      by = "Nodo"
    )
}


## 6. Pesos espaciales -----------------------------------------

coords <- as.matrix(
  base |>
    select(
      X,
      Y
    )
)

if (anyDuplicated(coords) > 0) {
  stop(
    "Existen coordenadas duplicadas."
  )
}

k_real <- min(
  k_vecinos,
  nrow(coords) - 1
)

vecinos <- knearneigh(
  coords,
  k = k_real
) |>
  knn2nb()

pesos <- nb2listw(
  vecinos,
  style = "W",
  zero.policy = TRUE
)


## 7. Moran global de residuos ---------------------------------

calcular_moran_residual <- function(
  respuesta,
  modelo_id,
  tipo_modelo
) {

  df <- extraer_residuos_loocv(
    respuesta,
    modelo_id
  )

  residuo <- df$Residuo_LOOCV

  validos <- is.finite(
    residuo
  )

  if (sum(validos) != nrow(df)) {
    return(
      tibble(
        Respuesta = respuesta,
        Modelo_ID = modelo_id,
        Tipo_modelo = tipo_modelo,
        Moran_I = NA_real_,
        Esperanza = NA_real_,
        Varianza = NA_real_,
        Sim_Q025 = NA_real_,
        Sim_Mediana = NA_real_,
        Sim_Q975 = NA_real_,
        Clasificacion = "No evaluable"
      )
    )
  }

  prueba <- moran.test(
    residuo,
    listw = pesos,
    zero.policy = TRUE,
    randomisation = TRUE
  )

  simulacion <- moran.mc(
    residuo,
    listw = pesos,
    nsim = nsim_moran,
    zero.policy = TRUE
  )

  sim <- as.numeric(
    simulacion$res
  )

  q <- quantile(
    sim,
    probs = c(
      0.025,
      0.50,
      0.975
    ),
    na.rm = TRUE,
    names = FALSE
  )

  observado <- as.numeric(
    prueba$estimate[
      "Moran I statistic"
    ]
  )

  clasificacion <- case_when(
    observado > q[3] ~
      "Autocorrelación residual positiva",
    observado < q[1] ~
      "Autocorrelación residual negativa",
    TRUE ~
      "Compatible con referencia permutada"
  )

  tibble(
    Respuesta = respuesta,
    Modelo_ID = modelo_id,
    Tipo_modelo = tipo_modelo,
    Moran_I = observado,
    Esperanza = as.numeric(
      prueba$estimate[
        "Expectation"
      ]
    ),
    Varianza = as.numeric(
      prueba$estimate[
        "Variance"
      ]
    ),
    Sim_Q025 = q[1],
    Sim_Mediana = q[2],
    Sim_Q975 = q[3],
    Clasificacion = clasificacion
  )
}


tabla_moran_residuos <- pmap_dfr(
  modelos_validacion,
  function(
    Respuesta,
    Modelo_ID,
    Formula,
    Tipo_modelo
  ) {
    calcular_moran_residual(
      respuesta = Respuesta,
      modelo_id = Modelo_ID,
      tipo_modelo = Tipo_modelo
    )
  }
)


## 8. Semivariogramas de residuos ------------------------------

calcular_variograma_residual <- function(
  respuesta,
  modelo_id,
  tipo_modelo
) {

  df <- extraer_residuos_loocv(
    respuesta,
    modelo_id
  ) |>
    filter(
      is.finite(
        Residuo_LOOCV
      )
    )

  if (nrow(df) < 8) {
    return(
      tibble()
    )
  }

  sf_residuo <- st_as_sf(
    df,
    coords = c(
      "X",
      "Y"
    ),
    crs = 32618,
    remove = FALSE
  )

  formula_variograma <- Residuo_LOOCV ~ 1

  variograma <- try(
    variogram(
      formula_variograma,
      sf_residuo,
      cutoff = max(
        dist(
          as.matrix(
            df |>
              select(
                X,
                Y
              )
          )
        )
      ) *
        0.60,
      width = max(
        dist(
          as.matrix(
            df |>
              select(
                X,
                Y
              )
          )
        )
      ) *
        0.60 /
        6
    ),
    silent = TRUE
  )

  if (inherits(variograma, "try-error")) {
    return(
      tibble()
    )
  }

  as_tibble(
    variograma
  ) |>
    transmute(
      Respuesta = respuesta,
      Modelo_ID = modelo_id,
      Tipo_modelo = tipo_modelo,
      Distancia_m = dist,
      Semivarianza = gamma,
      Pares = np
    )
}


tabla_variogramas <- pmap_dfr(
  modelos_validacion,
  function(
    Respuesta,
    Modelo_ID,
    Formula,
    Tipo_modelo
  ) {
    calcular_variograma_residual(
      respuesta = Respuesta,
      modelo_id = Modelo_ID,
      tipo_modelo = Tipo_modelo
    )
  }
)


## 9. Leave-one-replication-out CV -----------------------------

validacion_repeticion <- function(
  respuesta,
  modelo_id,
  formula_texto,
  tipo_modelo
) {

  formula_modelo <- as.formula(
    formula_texto
  )

  niveles_repeticion <- levels(
    base$Replication
  )

  predicciones <- map_dfr(
    niveles_repeticion,
    function(rep_excluida) {

      entrenamiento <- base |>
        filter(
          Replication != rep_excluida
        )

      prueba <- base |>
        filter(
          Replication == rep_excluida
        )

      predicho <- ajustar_y_predecir(
        formula = formula_modelo,
        entrenamiento = entrenamiento,
        prueba = prueba
      )

      tibble(
        Nodo = prueba$Nodo,
        Respuesta = respuesta,
        Modelo_ID = modelo_id,
        Tipo_modelo = tipo_modelo,
        Replication_excluida = rep_excluida,
        Observado = prueba[[respuesta]],
        Predicho = predicho,
        Residuo = prueba[[respuesta]] -
          predicho
      )
    }
  )

  metricas <- calcular_metricas(
    predicciones$Observado,
    predicciones$Predicho
  ) |>
    mutate(
      Respuesta = respuesta,
      Modelo_ID = modelo_id,
      Tipo_modelo = tipo_modelo,
      .before = 1
    )

  metricas_por_repeticion <- predicciones |>
    group_by(
      Respuesta,
      Modelo_ID,
      Tipo_modelo,
      Replication_excluida
    ) |>
    group_modify(
      function(.x, .y) {
        calcular_metricas(
          .x$Observado,
          .x$Predicho
        )
      }
    ) |>
    ungroup()

  list(
    predicciones = predicciones,
    metricas = metricas,
    metricas_por_repeticion = metricas_por_repeticion
  )
}


resultados_loro <- pmap(
  modelos_validacion,
  function(
    Respuesta,
    Modelo_ID,
    Formula,
    Tipo_modelo
  ) {
    validacion_repeticion(
      respuesta = Respuesta,
      modelo_id = Modelo_ID,
      formula_texto = Formula,
      tipo_modelo = Tipo_modelo
    )
  }
)

tabla_loro_predicciones <- bind_rows(
  map(
    resultados_loro,
    "predicciones"
  )
)

tabla_loro_metricas <- bind_rows(
  map(
    resultados_loro,
    "metricas"
  )
)

tabla_loro_por_repeticion <- bind_rows(
  map(
    resultados_loro,
    "metricas_por_repeticion"
  )
)


## 10. Comparación LORO M1 vs M2 -------------------------------

comparacion_loro <- tabla_loro_metricas |>
  select(
    Respuesta,
    Tipo_modelo,
    Modelo_ID,
    RMSE,
    MAE,
    Sesgo,
    Cor_Pearson,
    Cor_Spearman
  ) |>
  pivot_wider(
    names_from = Tipo_modelo,
    values_from = c(
      Modelo_ID,
      RMSE,
      MAE,
      Sesgo,
      Cor_Pearson,
      Cor_Spearman
    )
  ) |>
  mutate(
    Mejora_RMSE_porcentaje = ifelse(
      is.finite(RMSE_M1) &
        RMSE_M1 > epsilon,
      100 *
        (
          RMSE_M1 -
            RMSE_Mejor_M2
        ) /
        RMSE_M1,
      NA_real_
    ),
    Mejora_MAE_porcentaje = ifelse(
      is.finite(MAE_M1) &
        MAE_M1 > epsilon,
      100 *
        (
          MAE_M1 -
            MAE_Mejor_M2
        ) /
        MAE_M1,
      NA_real_
    ),
    Estado_LORO = case_when(
      Mejora_RMSE_porcentaje > 0 &
        Mejora_MAE_porcentaje > 0 ~
        "Mejora robusta por repetición",
      Mejora_RMSE_porcentaje > 0 ~
        "Mejora solo en RMSE",
      Mejora_MAE_porcentaje > 0 ~
        "Mejora solo en MAE",
      TRUE ~
        "Sin mejora por repetición"
    )
  )


## 11. Influencia y leverage ------------------------------------

extraer_influencia <- function(
  respuesta,
  modelo_id,
  tipo_modelo
) {

  clave <- paste(
    respuesta,
    modelo_id,
    sep = "__"
  )

  modelo <- objetos_modelos[[clave]]

  if (is.null(modelo)) {
    return(
      tibble()
    )
  }

  cooks <- cooks.distance(
    modelo
  )

  leverage <- hatvalues(
    modelo
  )

  resid_est <- rstandard(
    modelo
  )

  n <- nobs(
    modelo
  )

  p <- modelo$rank

  umbral_cook_local <- ifelse(
    is.null(umbral_cook),
    4 / n,
    umbral_cook
  )

  umbral_leverage <- umbral_leverage_factor *
    p /
    n

  filas_modelo <- as.integer(
    rownames(
      model.frame(
        modelo
      )
    )
  )

  tibble(
    Fila_base = filas_modelo,
    Nodo = base$Nodo[filas_modelo],
    Respuesta = respuesta,
    Modelo_ID = modelo_id,
    Tipo_modelo = tipo_modelo,
    Cook = as.numeric(
      cooks
    ),
    Leverage = as.numeric(
      leverage
    ),
    Residuo_estandarizado = as.numeric(
      resid_est
    ),
    Umbral_Cook = umbral_cook_local,
    Umbral_Leverage = umbral_leverage,
    Influyente_Cook = Cook >
      Umbral_Cook,
    Influyente_Leverage = Leverage >
      Umbral_Leverage,
    Residuo_extremo = abs(
      Residuo_estandarizado
    ) > 2,
    Estado_influencia = case_when(
      Influyente_Cook |
        Influyente_Leverage |
        Residuo_extremo ~
        "Revisar",
      TRUE ~
        "Sin señal fuerte"
    )
  )
}


tabla_influencia <- pmap_dfr(
  modelos_validacion,
  function(
    Respuesta,
    Modelo_ID,
    Formula,
    Tipo_modelo
  ) {
    extraer_influencia(
      respuesta = Respuesta,
      modelo_id = Modelo_ID,
      tipo_modelo = Tipo_modelo
    )
  }
)

resumen_influencia <- tabla_influencia |>
  group_by(
    Respuesta,
    Modelo_ID,
    Tipo_modelo
  ) |>
  summarise(
    n_Cook = sum(
      Influyente_Cook,
      na.rm = TRUE
    ),
    n_Leverage = sum(
      Influyente_Leverage,
      na.rm = TRUE
    ),
    n_Residuos_extremos = sum(
      Residuo_extremo,
      na.rm = TRUE
    ),
    Cook_maximo = max(
      Cook,
      na.rm = TRUE
    ),
    Leverage_maximo = max(
      Leverage,
      na.rm = TRUE
    ),
    .groups = "drop"
  )


## 12. Sensibilidad CE original vs corregida -------------------

sensibilidad_ce <- sensibilidad_ce_06C |>
  mutate(
    Equivalencia_RMSE = abs(
      Diferencia_RMSE_corregida_menos_original
    ) < 1e-8,
    Equivalencia_MAE = abs(
      Diferencia_MAE_corregida_menos_original
    ) < 1e-8,
    Interpretacion = case_when(
      Equivalencia_RMSE &
        Equivalencia_MAE ~
        "Resultados numéricamente equivalentes",
      TRUE ~
        "Resultados diferentes"
    )
  )


## 13. Síntesis integrada --------------------------------------

moran_ancho <- tabla_moran_residuos |>
  select(
    Respuesta,
    Tipo_modelo,
    Moran_I,
    Clasificacion
  ) |>
  pivot_wider(
    names_from = Tipo_modelo,
    values_from = c(
      Moran_I,
      Clasificacion
    )
  )

influencia_ancha <- resumen_influencia |>
  select(
    Respuesta,
    Tipo_modelo,
    n_Cook,
    n_Leverage,
    n_Residuos_extremos
  ) |>
  pivot_wider(
    names_from = Tipo_modelo,
    values_from = c(
      n_Cook,
      n_Leverage,
      n_Residuos_extremos
    )
  )

tabla_sintesis <- comparacion_loro |>
  left_join(
    moran_ancho,
    by = "Respuesta"
  ) |>
  left_join(
    influencia_ancha,
    by = "Respuesta"
  ) |>
  mutate(
    Reduccion_abs_Moran = abs(
      Moran_I_M1
    ) -
      abs(
        Moran_I_Mejor_M2
      ),
    Robustez_integrada = case_when(
      Estado_LORO ==
        "Mejora robusta por repetición" &
        Reduccion_abs_Moran > 0 ~
        "Alta",
      Estado_LORO ==
        "Mejora robusta por repetición" |
        Reduccion_abs_Moran > 0 ~
        "Moderada",
      TRUE ~
        "Baja"
    ),
    Advertencia = case_when(
      str_detect(
        Modelo_ID_Mejor_M2,
        "CE_corregida"
      ) ~
        "CE corregida equivalente linealmente a CE original",
      str_detect(
        Modelo_ID_Mejor_M2,
        "Indice_"
      ) ~
        "Índice instrumental nominal",
      str_detect(
        Modelo_ID_Mejor_M2,
        "Exploratorio"
      ) ~
        "pH y Hv con alcance exploratorio",
      TRUE ~
        "Sin advertencia adicional"
    )
  )


## 14. Figuras -------------------------------------------------

g_moran <- tabla_moran_residuos |>
  ggplot(
    aes(
      x = Tipo_modelo,
      y = Moran_I,
      shape = Clasificacion
    )
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  geom_point(
    size = 3
  ) +
  facet_wrap(
    ~ Respuesta,
    scales = "free_y",
    ncol = 3
  ) +
  labs(
    title = "Autocorrelación espacial de residuos LOOCV",
    x = NULL,
    y = "Moran I",
    shape = "Clasificación"
  ) +
  theme_minimal()

ggsave(
  filename = file.path(
    ruta_figuras,
    "01_Moran_residuos_M1_M2.png"
  ),
  plot = g_moran,
  width = 13,
  height = 10,
  dpi = 400
)


g_variogramas <- tabla_variogramas |>
  ggplot(
    aes(
      x = Distancia_m,
      y = Semivarianza,
      shape = Tipo_modelo
    )
  ) +
  geom_line(
    aes(
      group = Tipo_modelo
    )
  ) +
  geom_point(
    aes(
      size = Pares
    )
  ) +
  facet_wrap(
    ~ Respuesta,
    scales = "free_y",
    ncol = 3
  ) +
  labs(
    title = "Semivariogramas experimentales de residuos LOOCV",
    x = "Distancia (m)",
    y = "Semivarianza",
    shape = "Modelo",
    size = "Pares"
  ) +
  theme_minimal()

ggsave(
  filename = file.path(
    ruta_figuras,
    "02_semivariogramas_residuos_M1_M2.png"
  ),
  plot = g_variogramas,
  width = 14,
  height = 11,
  dpi = 400
)


g_loro <- comparacion_loro |>
  select(
    Respuesta,
    Mejora_RMSE_porcentaje,
    Mejora_MAE_porcentaje
  ) |>
  pivot_longer(
    cols = c(
      Mejora_RMSE_porcentaje,
      Mejora_MAE_porcentaje
    ),
    names_to = "Metrica",
    values_to = "Mejora"
  ) |>
  mutate(
    Metrica = case_when(
      Metrica ==
        "Mejora_RMSE_porcentaje" ~
        "RMSE",
      Metrica ==
        "Mejora_MAE_porcentaje" ~
        "MAE",
      TRUE ~
        Metrica
    )
  ) |>
  ggplot(
    aes(
      x = Respuesta,
      y = Mejora,
      shape = Metrica
    )
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  geom_point(
    size = 3,
    position = position_dodge(
      width = 0.35
    )
  ) +
  labs(
    title = "Validación dejando una repetición por fuera",
    subtitle = "Valores positivos favorecen el mejor modelo M2",
    x = NULL,
    y = "Mejora frente a M1 (%)",
    shape = "Métrica"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(
      angle = 35,
      hjust = 1
    )
  )

ggsave(
  filename = file.path(
    ruta_figuras,
    "03_validacion_leave_one_replication_out.png"
  ),
  plot = g_loro,
  width = 11,
  height = 7,
  dpi = 400
)


g_influencia <- tabla_influencia |>
  ggplot(
    aes(
      x = Leverage,
      y = Cook,
      shape = Estado_influencia
    )
  ) +
  geom_point() +
  facet_grid(
    Tipo_modelo ~ Respuesta,
    scales = "free"
  ) +
  labs(
    title = "Diagnóstico de influencia de los modelos",
    x = "Leverage",
    y = "Distancia de Cook",
    shape = "Estado"
  ) +
  theme_minimal()

ggsave(
  filename = file.path(
    ruta_figuras,
    "04_influencia_Cook_leverage.png"
  ),
  plot = g_influencia,
  width = 16,
  height = 8,
  dpi = 400
)


residuos_mapa <- tabla_predicciones |>
  semi_join(
    modelos_validacion |>
      select(
        Respuesta,
        Modelo_ID
      ),
    by = c(
      "Respuesta",
      "Modelo_ID"
    )
  ) |>
  left_join(
    base |>
      select(
        Nodo,
        X,
        Y
      ),
    by = "Nodo"
  ) |>
  left_join(
    modelos_validacion |>
      select(
        Respuesta,
        Modelo_ID,
        Tipo_modelo
      ),
    by = c(
      "Respuesta",
      "Modelo_ID"
    )
  )

g_mapas_residuos <- ggplot(
  residuos_mapa,
  aes(
    x = X,
    y = Y,
    size = abs(
      Residuo_LOOCV
    ),
    shape = Tipo_modelo
  )
) +
  geom_point() +
  facet_wrap(
    ~ Respuesta,
    scales = "fixed",
    ncol = 3
  ) +
  coord_equal() +
  labs(
    title = "Distribución espacial de residuos LOOCV",
    x = "X (m)",
    y = "Y (m)",
    size = "|Residuo|",
    shape = "Modelo"
  ) +
  theme_minimal()

ggsave(
  filename = file.path(
    ruta_figuras,
    "05_mapas_residuos_LOOCV.png"
  ),
  plot = g_mapas_residuos,
  width = 14,
  height = 11,
  dpi = 400
)


panel_validacion <- g_moran /
  g_loro +
  plot_annotation(
    title = "Validación espacial y por repetición de modelos integrados"
  )

ggsave(
  filename = file.path(
    ruta_figuras,
    "06_panel_validacion_integrada.png"
  ),
  plot = panel_validacion,
  width = 13,
  height = 18,
  dpi = 400
)


## 15. Exportar resultados -------------------------------------

write_csv(
  tabla_moran_residuos,
  file.path(
    ruta_datos,
    "Moran_residuos_modelos.csv"
  )
)

write_csv(
  tabla_variogramas,
  file.path(
    ruta_datos,
    "Semivariogramas_residuos.csv"
  )
)

write_csv(
  tabla_loro_predicciones,
  file.path(
    ruta_datos,
    "Predicciones_leave_one_replication_out.csv"
  )
)

write_csv(
  tabla_loro_metricas,
  file.path(
    ruta_datos,
    "Metricas_leave_one_replication_out.csv"
  )
)

write_csv(
  tabla_influencia,
  file.path(
    ruta_datos,
    "Diagnostico_influencia_nodos.csv"
  )
)

write_csv(
  tabla_sintesis,
  file.path(
    ruta_datos,
    "Sintesis_validacion_sensibilidad.csv"
  )
)


write_xlsx(
  list(
    Sintesis =
      tabla_sintesis,

    Moran_residuos =
      tabla_moran_residuos,

    Variogramas_residuos =
      tabla_variogramas,

    LORO_metricas =
      tabla_loro_metricas,

    LORO_por_repeticion =
      tabla_loro_por_repeticion,

    LORO_predicciones =
      tabla_loro_predicciones,

    Influencia_resumen =
      resumen_influencia,

    Influencia_nodos =
      tabla_influencia,

    Sensibilidad_CE =
      sensibilidad_ce
  ),
  path = file.path(
    ruta_tablas,
    "Resultados_validacion_sensibilidad_modelos.xlsx"
  )
)


## 16. Guardar objeto para Script 06-E -------------------------

saveRDS(
  list(
    base = base,
    modelos_validacion = modelos_validacion,
    tabla_moran_residuos = tabla_moran_residuos,
    tabla_variogramas = tabla_variogramas,
    tabla_loro_metricas = tabla_loro_metricas,
    tabla_loro_por_repeticion = tabla_loro_por_repeticion,
    tabla_loro_predicciones = tabla_loro_predicciones,
    tabla_influencia = tabla_influencia,
    resumen_influencia = resumen_influencia,
    sensibilidad_ce = sensibilidad_ce,
    tabla_sintesis = tabla_sintesis
  ),
  file = file.path(
    ruta_objetos,
    "objetos_validacion_sensibilidad_06D.rds"
  )
)


## 17. Registro de sesión --------------------------------------

capture.output(
  sessionInfo(),
  file = file.path(
    ruta_salida,
    "sessionInfo_06_D.txt"
  )
)


## 18. Resumen en consola --------------------------------------

cat(
  "\n============================================\n"
)

cat(
  "SCRIPT 06-D FINALIZADO\n\n"
)

cat(
  "Síntesis de validación:\n"
)

print(
  tabla_sintesis |>
    select(
      Respuesta,
      Modelo_ID_Mejor_M2,
      Mejora_RMSE_porcentaje,
      Mejora_MAE_porcentaje,
      Moran_I_M1,
      Moran_I_Mejor_M2,
      Reduccion_abs_Moran,
      Estado_LORO,
      Robustez_integrada
    )
)

cat(
  "\nSensibilidad CE original vs corregida:\n"
)

print(
  sensibilidad_ce |>
    select(
      Respuesta,
      Diferencia_RMSE_corregida_menos_original,
      Diferencia_MAE_corregida_menos_original,
      Interpretacion
    )
)

cat(
  "\nObjeto para Script 06-E:\n",
  file.path(
    ruta_objetos,
    "objetos_validacion_sensibilidad_06D.rds"
  ),
  "\n"
)

cat(
  "============================================\n"
)
