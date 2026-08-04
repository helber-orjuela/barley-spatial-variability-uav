
## ============================================================
## 06_SUELO_CULTIVO – SCRIPT 06-C
## MODELOS INTEGRADOS SUELO–CULTIVO Y VALIDACIÓN LOOCV
## ============================================================
##
## Entrada:
##   06_SUELO_CULTIVO/02_SELECCION_COVARIABLES/objetos_R/
##   objetos_seleccion_covariables_06B.rds
##
## Salida:
##   06_SUELO_CULTIVO/03_MODELOS_INTEGRADOS/
##
## PROPÓSITO
##   1. Comparar modelos de diseño experimental con modelos que
##      incorporan covariables de suelo.
##   2. Evaluar si las superficies edáficas aportan capacidad
##      predictiva adicional sobre las respuestas del cultivo.
##   3. Comparar CE original y CE corregida OLS.
##   4. Evaluar por separado los índices instrumentales N, P y K.
##   5. Priorizar RMSE, MAE, correlación observado–predicho,
##      tamaño del efecto y estabilidad, no p-valores aislados.
##
## RESPUESTAS
##   Macollas, Altura, Diametro, NDVI, Nitratos,
##   Potasio_savia y Calcio_savia.
##
## MODELOS
##   M0: respuesta ~ Replication
##   M1: respuesta ~ Replication + Dose_percent +
##                   Omission_description
##   M2: M1 + conjunto candidato de covariables de suelo
##
## VALIDACIÓN
##   - LOOCV sobre los 56 nodos.
##   - La matriz de diseño se construye una sola vez por modelo.
##   - Se conserva la codificación completa de los factores.
##   - Los modelos con deficiencia de rango se identifican.
##
## ADVERTENCIAS
##   - Las variables de suelo son predicciones de Kriging.
##   - N, P y K son índices instrumentales nominales.
##   - CE corregida OLS procede de calibración exploratoria.
##   - pH y Hv se consideran covariables exploratorias.
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
library(tibble)


## 1. Parámetros ------------------------------------------------

set.seed(2026)

epsilon <- 1e-12
umbral_mejora_relevante <- 5
umbral_cor_predicha <- 0.30


## 2. Rutas ----------------------------------------------------

ruta_base <- paste0(
  "D:/UDCA_2017_CLASE/46_ARTICULOS_2025/",
  "ARTICULO_CEBADA_DOCTORADO/AJUSTES_ARTÍCULO"
)

ruta_objeto <- file.path(
  ruta_base,
  "06_SUELO_CULTIVO",
  "02_SELECCION_COVARIABLES",
  "objetos_R",
  "objetos_seleccion_covariables_06B.rds"
)

ruta_salida <- file.path(
  ruta_base,
  "06_SUELO_CULTIVO",
  "03_MODELOS_INTEGRADOS"
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
    "No existe el objeto del Script 06-B:\n",
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
  obj$base_modelo
)

variables_cultivo <- intersect(
  obj$variables_cultivo,
  names(base)
)

recomendacion_modelos <- obj$recomendacion_modelos

if (length(variables_cultivo) == 0) {
  stop(
    "No se encontraron variables de cultivo."
  )
}

if (
  is.null(recomendacion_modelos) ||
    nrow(recomendacion_modelos) == 0
) {
  stop(
    "No existe la tabla de escenarios recomendados de 06-B."
  )
}


## 4. Preparar factores y variables ----------------------------

columnas_diseno_requeridas <- c(
  "Replication",
  "Dose_percent",
  "Omission_description"
)

faltantes_diseno <- setdiff(
  columnas_diseno_requeridas,
  names(base)
)

if (length(faltantes_diseno) > 0) {
  stop(
    "Faltan columnas del diseño experimental: ",
    paste(
      faltantes_diseno,
      collapse = ", "
    )
  )
}

base <- base |>
  mutate(
    Replication = factor(
      Replication
    ),
    Dose_percent = factor(
      as.character(
        Dose_percent
      )
    ),
    Omission_description = factor(
      Omission_description
    )
  )

for (variable in variables_cultivo) {
  base[[variable]] <- as.numeric(
    base[[variable]]
  )
}


## 5. Funciones auxiliares -------------------------------------

separar_covariables <- function(texto) {

  if (
    is.na(texto) ||
      texto == ""
  ) {
    return(
      character(0)
    )
  }

  str_split(
    texto,
    "\\s*\\+\\s*"
  )[[1]]
}


crear_formula <- function(
  respuesta,
  predictores
) {

  if (length(predictores) == 0) {
    return(
      as.formula(
        paste0(
          respuesta,
          " ~ 1"
        )
      )
    )
  }

  as.formula(
    paste0(
      respuesta,
      " ~ ",
      paste(
        predictores,
        collapse = " + "
      )
    )
  )
}


correlacion_segura <- function(
  x,
  y,
  metodo = "pearson"
) {

  validos <- is.finite(
    x
  ) &
    is.finite(
      y
    )

  if (sum(validos) < 4) {
    return(
      NA_real_
    )
  }

  if (
    sd(
      x[validos]
    ) <= epsilon ||
      sd(
        y[validos]
      ) <= epsilon
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


ajustar_matriz <- function(
  formula,
  datos
) {

  marco <- model.frame(
    formula,
    data = datos,
    na.action = na.pass,
    drop.unused.levels = FALSE
  )

  y <- model.response(
    marco
  )

  X <- model.matrix(
    formula,
    data = marco
  )

  validos <- is.finite(
    y
  ) &
    apply(
      X,
      1,
      function(fila) {
        all(
          is.finite(
            fila
          )
        )
      }
    )

  list(
    y = as.numeric(
      y
    ),
    X = X,
    validos = validos,
    marco = marco
  )
}


loocv_matriz <- function(
  formula,
  datos
) {

  componentes <- ajustar_matriz(
    formula,
    datos
  )

  y <- componentes$y
  X <- componentes$X
  validos <- componentes$validos

  indices_validos <- which(
    validos
  )

  predicciones <- rep(
    NA_real_,
    length(
      y
    )
  )

  rangos <- rep(
    NA_integer_,
    length(
      y
    )
  )

  if (length(indices_validos) < 4) {
    return(
      list(
        predicciones = predicciones,
        rangos = rangos,
        n = length(
          indices_validos
        )
      )
    )
  }

  for (indice_prueba in indices_validos) {

    indices_entrenamiento <- setdiff(
      indices_validos,
      indice_prueba
    )

    X_train <- X[
      indices_entrenamiento,
      ,
      drop = FALSE
    ]

    y_train <- y[
      indices_entrenamiento
    ]

    ajuste <- lm.fit(
      x = X_train,
      y = y_train
    )

    coeficientes <- ajuste$coefficients

    coeficientes[
      !is.finite(
        coeficientes
      )
    ] <- 0

    predicciones[
      indice_prueba
    ] <- as.numeric(
      X[
        indice_prueba,
        ,
        drop = FALSE
      ] %*%
        coeficientes
    )

    rangos[
      indice_prueba
    ] <- ajuste$rank
  }

  list(
    predicciones = predicciones,
    rangos = rangos,
    n = length(
      indices_validos
    )
  )
}


calcular_metricas <- function(
  observado,
  predicho
) {

  validos <- is.finite(
    observado
  ) &
    is.finite(
      predicho
    )

  observado <- observado[
    validos
  ]

  predicho <- predicho[
    validos
  ]

  if (length(observado) < 4) {
    return(
      tibble(
        n = length(
          observado
        ),
        RMSE = NA_real_,
        MAE = NA_real_,
        Sesgo = NA_real_,
        RRMSE_porcentaje = NA_real_,
        Cor_Pearson = NA_real_,
        Cor_Spearman = NA_real_
      )
    )
  }

  errores <- predicho -
    observado

  media_observada <- mean(
    observado
  )

  tibble(
    n = length(
      observado
    ),
    RMSE = sqrt(
      mean(
        errores^2
      )
    ),
    MAE = mean(
      abs(
        errores
      )
    ),
    Sesgo = mean(
      errores
    ),
    RRMSE_porcentaje = ifelse(
      abs(
        media_observada
      ) > epsilon,
      100 *
        sqrt(
          mean(
            errores^2
          )
        ) /
        abs(
          media_observada
        ),
      NA_real_
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


ajuste_final_seguro <- function(
  formula,
  datos
) {

  modelo <- try(
    lm(
      formula,
      data = datos,
      na.action = na.exclude
    ),
    silent = TRUE
  )

  if (inherits(modelo, "try-error")) {
    return(
      NULL
    )
  }

  modelo
}


extraer_coeficientes <- function(
  modelo,
  respuesta,
  modelo_id,
  escenario,
  familia
) {

  if (is.null(modelo)) {
    return(
      tibble()
    )
  }

  coeficientes <- coef(
    modelo
  )

  errores_estandar <- try(
    sqrt(
      diag(
        vcov(
          modelo
        )
      )
    ),
    silent = TRUE
  )

  if (inherits(errores_estandar, "try-error")) {
    errores_estandar <- rep(
      NA_real_,
      length(
        coeficientes
      )
    )
  }

  tibble(
    Respuesta = respuesta,
    Modelo_ID = modelo_id,
    Familia = familia,
    Escenario = escenario,
    Termino = names(
      coeficientes
    ),
    Estimacion = as.numeric(
      coeficientes
    ),
    Error_estandar = as.numeric(
      errores_estandar
    ),
    IC95_inferior = Estimacion -
      1.96 *
        Error_estandar,
    IC95_superior = Estimacion +
      1.96 *
        Error_estandar,
    Coeficiente_alias = !is.finite(
      Estimacion
    )
  )
}


calcular_resumen_ajuste <- function(
  modelo
) {

  if (is.null(modelo)) {
    return(
      tibble(
        R2 = NA_real_,
        R2_ajustado = NA_real_,
        AIC = NA_real_,
        BIC = NA_real_,
        Sigma_residual = NA_real_,
        Rango_modelo = NA_integer_,
        n_coeficientes = NA_integer_,
        n_alias = NA_integer_
      )
    )
  }

  resumen <- summary(
    modelo
  )

  tibble(
    R2 = resumen$r.squared,
    R2_ajustado = resumen$adj.r.squared,
    AIC = AIC(
      modelo
    ),
    BIC = BIC(
      modelo
    ),
    Sigma_residual = resumen$sigma,
    Rango_modelo = modelo$rank,
    n_coeficientes = length(
      coef(
        modelo
      )
    ),
    n_alias = sum(
      !is.finite(
        coef(
          modelo
        )
      )
    )
  )
}


## 6. Definir biblioteca de modelos ----------------------------

predictores_M0 <- c(
  "Replication"
)

predictores_M1 <- c(
  "Replication",
  "Dose_percent",
  "Omission_description"
)

biblioteca_modelos <- bind_rows(
  tibble(
    Modelo_ID = "M0_Repeticion",
    Familia = "M0",
    Escenario = "Diseño_reducido",
    Covariables = "",
    Predictores = list(
      predictores_M0
    ),
    Advertencia = "Modelo de referencia con repetición"
  ),
  tibble(
    Modelo_ID = "M1_Tratamientos",
    Familia = "M1",
    Escenario = "Diseño_completo",
    Covariables = "",
    Predictores = list(
      predictores_M1
    ),
    Advertencia = "Modelo base para evaluar aporte adicional del suelo"
  ),
  recomendacion_modelos |>
    transmute(
      Modelo_ID = Escenario,
      Familia = "M2",
      Escenario = Escenario,
      Covariables = Covariables,
      Predictores = map(
        Covariables,
        function(texto) {
          unique(
            c(
              predictores_M1,
              separar_covariables(
                texto
              )
            )
          )
        }
      ),
      Advertencia = Advertencia
    )
)


## 7. Ejecutar modelos -----------------------------------------

resultados_metricas <- list()
resultados_predicciones <- list()
resultados_coeficientes <- list()
resultados_ajuste <- list()
objetos_modelos <- list()

contador <- 1

for (respuesta in variables_cultivo) {

  for (
    indice_modelo in seq_len(
      nrow(
        biblioteca_modelos
      )
    )
  ) {

    modelo_id <- biblioteca_modelos$Modelo_ID[
      indice_modelo
    ]

    familia <- biblioteca_modelos$Familia[
      indice_modelo
    ]

    escenario <- biblioteca_modelos$Escenario[
      indice_modelo
    ]

    predictores <- biblioteca_modelos$Predictores[[indice_modelo]]

    predictores_disponibles <- intersect(
      predictores,
      names(
        base
      )
    )

    formula_modelo <- crear_formula(
      respuesta = respuesta,
      predictores = predictores_disponibles
    )

    validacion <- loocv_matriz(
      formula = formula_modelo,
      datos = base
    )

    observado <- base[[respuesta]]
    predicho <- validacion$predicciones

    metricas <- calcular_metricas(
      observado = observado,
      predicho = predicho
    ) |>
      mutate(
        Respuesta = respuesta,
        Modelo_ID = modelo_id,
        Familia = familia,
        Escenario = escenario,
        Formula = paste(
          deparse(
            formula_modelo
          ),
          collapse = ""
        ),
        n_predictores_nominales = length(
          predictores_disponibles
        ),
        Rango_LOOCV_min = ifelse(
          all(
            is.na(
              validacion$rangos
            )
          ),
          NA_integer_,
          min(
            validacion$rangos,
            na.rm = TRUE
          )
        ),
        Rango_LOOCV_max = ifelse(
          all(
            is.na(
              validacion$rangos
            )
          ),
          NA_integer_,
          max(
            validacion$rangos,
            na.rm = TRUE
          )
        ),
        .before = 1
      )

    modelo_final <- ajuste_final_seguro(
      formula = formula_modelo,
      datos = base
    )

    resumen_ajuste <- calcular_resumen_ajuste(
      modelo_final
    ) |>
      mutate(
        Respuesta = respuesta,
        Modelo_ID = modelo_id,
        Familia = familia,
        Escenario = escenario,
        .before = 1
      )

    coeficientes <- extraer_coeficientes(
      modelo = modelo_final,
      respuesta = respuesta,
      modelo_id = modelo_id,
      escenario = escenario,
      familia = familia
    )

    predicciones <- tibble(
      Nodo = base$Nodo,
      Respuesta = respuesta,
      Modelo_ID = modelo_id,
      Familia = familia,
      Escenario = escenario,
      Observado = observado,
      Predicho_LOOCV = predicho,
      Residuo_LOOCV = observado -
        predicho
    )

    clave <- paste(
      respuesta,
      modelo_id,
      sep = "__"
    )

    objetos_modelos[[clave]] <- modelo_final

    resultados_metricas[[contador]] <- metricas

    resultados_predicciones[[contador]] <- predicciones

    resultados_coeficientes[[contador]] <- coeficientes

    resultados_ajuste[[contador]] <- resumen_ajuste

    contador <- contador + 1
  }
}

tabla_metricas <- bind_rows(
  resultados_metricas
)

tabla_predicciones <- bind_rows(
  resultados_predicciones
)

tabla_coeficientes <- bind_rows(
  resultados_coeficientes
)

tabla_ajuste <- bind_rows(
  resultados_ajuste
)


## 8. Comparar contra M1 ---------------------------------------

referencia_M1 <- tabla_metricas |>
  filter(
    Modelo_ID == "M1_Tratamientos"
  ) |>
  select(
    Respuesta,
    RMSE_M1 = RMSE,
    MAE_M1 = MAE,
    Cor_Pearson_M1 = Cor_Pearson,
    Cor_Spearman_M1 = Cor_Spearman
  )

tabla_comparacion <- tabla_metricas |>
  left_join(
    referencia_M1,
    by = "Respuesta"
  ) |>
  mutate(
    Delta_RMSE = RMSE -
      RMSE_M1,
    Mejora_RMSE_porcentaje = ifelse(
      is.finite(
        RMSE_M1
      ) &
        RMSE_M1 > epsilon,
      100 *
        (
          RMSE_M1 -
            RMSE
        ) /
        RMSE_M1,
      NA_real_
    ),
    Delta_MAE = MAE -
      MAE_M1,
    Mejora_MAE_porcentaje = ifelse(
      is.finite(
        MAE_M1
      ) &
        MAE_M1 > epsilon,
      100 *
        (
          MAE_M1 -
            MAE
        ) /
        MAE_M1,
      NA_real_
    ),
    Cambio_cor_Pearson = Cor_Pearson -
      Cor_Pearson_M1,
    Cambio_cor_Spearman = Cor_Spearman -
      Cor_Spearman_M1
  ) |>
  left_join(
    tabla_ajuste,
    by = c(
      "Respuesta",
      "Modelo_ID",
      "Familia",
      "Escenario"
    )
  )


## 9. Tamaño del efecto incremental ----------------------------

sse_por_modelo <- tabla_predicciones |>
  filter(
    is.finite(
      Residuo_LOOCV
    )
  ) |>
  group_by(
    Respuesta,
    Modelo_ID
  ) |>
  summarise(
    SSE_LOOCV = sum(
      Residuo_LOOCV^2
    ),
    .groups = "drop"
  )

sse_M1 <- sse_por_modelo |>
  filter(
    Modelo_ID == "M1_Tratamientos"
  ) |>
  select(
    Respuesta,
    SSE_M1 = SSE_LOOCV
  )

tabla_efecto <- sse_por_modelo |>
  left_join(
    sse_M1,
    by = "Respuesta"
  ) |>
  mutate(
    R2_incremental_LOOCV = ifelse(
      is.finite(
        SSE_M1
      ) &
        SSE_M1 > epsilon,
      (
        SSE_M1 -
          SSE_LOOCV
      ) /
        SSE_M1,
      NA_real_
    )
  )


tabla_comparacion <- tabla_comparacion |>
  left_join(
    tabla_efecto |>
      select(
        Respuesta,
        Modelo_ID,
        SSE_LOOCV,
        SSE_M1,
        R2_incremental_LOOCV
      ),
    by = c(
      "Respuesta",
      "Modelo_ID"
    )
  )


## 10. Clasificar desempeño ------------------------------------

tabla_comparacion <- tabla_comparacion |>
  mutate(
    Estado_modelo = case_when(
      n_alias > 0 ~
        "Revisar: coeficientes alias",

      !is.finite(
        RMSE
      ) ~
        "No evaluable",

      Modelo_ID == "M1_Tratamientos" ~
        "Referencia",

      Familia == "M0" ~
        "Modelo reducido",

      Mejora_RMSE_porcentaje >=
        umbral_mejora_relevante &
        Mejora_MAE_porcentaje >= 0 &
        Cor_Pearson >=
          umbral_cor_predicha ~
        "Mejora predictiva relevante",

      Mejora_RMSE_porcentaje > 0 &
        Mejora_MAE_porcentaje > 0 ~
        "Mejora predictiva limitada",

      Mejora_RMSE_porcentaje > 0 ~
        "Mejora solo en RMSE",

      TRUE ~
        "Sin mejora frente a M1"
    ),
    Evidencia_integrada = case_when(
      Estado_modelo ==
        "Mejora predictiva relevante" &
        R2_incremental_LOOCV >= 0.10 ~
        "Alta",

      Estado_modelo %in% c(
        "Mejora predictiva relevante",
        "Mejora predictiva limitada"
      ) ~
        "Moderada",

      TRUE ~
        "Baja"
    )
  )


## 11. Selección del mejor modelo por respuesta ----------------

mejores_modelos <- tabla_comparacion |>
  filter(
    Familia == "M2",
    is.finite(
      RMSE
    )
  ) |>
  group_by(
    Respuesta
  ) |>
  arrange(
    RMSE,
    MAE,
    desc(
      Cor_Pearson
    )
  ) |>
  slice(
    1
  ) |>
  ungroup() |>
  mutate(
    Supera_M1 = Mejora_RMSE_porcentaje > 0 &
      Mejora_MAE_porcentaje > 0
  )


## 12. Sensibilidad CE original vs corregida -------------------

sensibilidad_ce <- tabla_comparacion |>
  filter(
    Modelo_ID %in% c(
      "M2_CE_original",
      "M2_CE_corregida"
    )
  ) |>
  select(
    Respuesta,
    Modelo_ID,
    RMSE,
    MAE,
    Cor_Pearson,
    Cor_Spearman,
    Mejora_RMSE_porcentaje,
    Mejora_MAE_porcentaje,
    R2_incremental_LOOCV
  ) |>
  pivot_wider(
    names_from = Modelo_ID,
    values_from = c(
      RMSE,
      MAE,
      Cor_Pearson,
      Cor_Spearman,
      Mejora_RMSE_porcentaje,
      Mejora_MAE_porcentaje,
      R2_incremental_LOOCV
    )
  ) |>
  mutate(
    Diferencia_RMSE_corregida_menos_original =
      RMSE_M2_CE_corregida -
      RMSE_M2_CE_original,
    Diferencia_MAE_corregida_menos_original =
      MAE_M2_CE_corregida -
      MAE_M2_CE_original
  )


## 13. Figuras -------------------------------------------------

g_rmse <- tabla_comparacion |>
  filter(
    Familia %in% c(
      "M1",
      "M2"
    )
  ) |>
  ggplot(
    aes(
      x = Modelo_ID,
      y = Mejora_RMSE_porcentaje,
      shape = Estado_modelo
    )
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  geom_point(
    size = 2.5
  ) +
  facet_wrap(
    ~ Respuesta,
    scales = "free_y",
    ncol = 3
  ) +
  labs(
    title = "Cambio del RMSE frente al modelo de tratamientos",
    subtitle = "Valores positivos indican mejora mediante covariables de suelo",
    x = NULL,
    y = "Mejora del RMSE (%)",
    shape = "Clasificación"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 7
    )
  )

ggsave(
  filename = file.path(
    ruta_figuras,
    "01_mejora_RMSE_modelos_integrados.png"
  ),
  plot = g_rmse,
  width = 15,
  height = 11,
  dpi = 400
)


g_mae <- tabla_comparacion |>
  filter(
    Familia %in% c(
      "M1",
      "M2"
    )
  ) |>
  ggplot(
    aes(
      x = Modelo_ID,
      y = Mejora_MAE_porcentaje,
      shape = Estado_modelo
    )
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  geom_point(
    size = 2.5
  ) +
  facet_wrap(
    ~ Respuesta,
    scales = "free_y",
    ncol = 3
  ) +
  labs(
    title = "Cambio del MAE frente al modelo de tratamientos",
    subtitle = "Valores positivos indican mejora mediante covariables de suelo",
    x = NULL,
    y = "Mejora del MAE (%)",
    shape = "Clasificación"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 7
    )
  )

ggsave(
  filename = file.path(
    ruta_figuras,
    "02_mejora_MAE_modelos_integrados.png"
  ),
  plot = g_mae,
  width = 15,
  height = 11,
  dpi = 400
)


pred_mejores <- tabla_predicciones |>
  semi_join(
    mejores_modelos |>
      select(
        Respuesta,
        Modelo_ID
      ),
    by = c(
      "Respuesta",
      "Modelo_ID"
    )
  )

g_obs_pred <- ggplot(
  pred_mejores,
  aes(
    x = Observado,
    y = Predicho_LOOCV
  )
) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed"
  ) +
  geom_point() +
  facet_wrap(
    ~ Respuesta,
    scales = "free",
    ncol = 3
  ) +
  labs(
    title = "Observado frente a predicho para el mejor modelo M2",
    subtitle = "Predicciones obtenidas mediante LOOCV",
    x = "Observado",
    y = "Predicho LOOCV"
  ) +
  theme_minimal()

ggsave(
  filename = file.path(
    ruta_figuras,
    "03_observado_predicho_mejores_modelos.png"
  ),
  plot = g_obs_pred,
  width = 13,
  height = 11,
  dpi = 400
)


g_ce <- sensibilidad_ce |>
  select(
    Respuesta,
    starts_with(
      "Mejora_RMSE_porcentaje_"
    )
  ) |>
  pivot_longer(
    cols = -Respuesta,
    names_to = "Escenario",
    values_to = "Mejora_RMSE"
  ) |>
  mutate(
    Escenario = case_when(
      str_detect(
        Escenario,
        "corregida"
      ) ~
        "CE corregida OLS",

      TRUE ~
        "CE original"
    )
  ) |>
  ggplot(
    aes(
      x = Respuesta,
      y = Mejora_RMSE,
      shape = Escenario
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
    title = "Sensibilidad a CE original y CE corregida",
    x = NULL,
    y = "Mejora del RMSE frente a M1 (%)",
    shape = "Escenario"
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
    "04_sensibilidad_CE_original_corregida.png"
  ),
  plot = g_ce,
  width = 11,
  height = 7,
  dpi = 400
)


panel_metricas <- g_rmse /
  g_mae +
  plot_annotation(
    title = "Desempeño predictivo de los modelos integrados suelo–cultivo"
  )

ggsave(
  filename = file.path(
    ruta_figuras,
    "05_panel_metricas_modelos_integrados.png"
  ),
  plot = panel_metricas,
  width = 15,
  height = 22,
  dpi = 400
)


## 14. Exportar resultados -------------------------------------

write_csv(
  tabla_comparacion,
  file.path(
    ruta_datos,
    "Comparacion_modelos_integrados.csv"
  )
)

write_csv(
  tabla_predicciones,
  file.path(
    ruta_datos,
    "Predicciones_LOOCV_modelos_integrados.csv"
  )
)

write_csv(
  tabla_coeficientes,
  file.path(
    ruta_datos,
    "Coeficientes_modelos_integrados.csv"
  )
)

write_csv(
  mejores_modelos,
  file.path(
    ruta_datos,
    "Mejores_modelos_por_respuesta.csv"
  )
)

write_csv(
  sensibilidad_ce,
  file.path(
    ruta_datos,
    "Sensibilidad_CE_original_corregida.csv"
  )
)


write_xlsx(
  list(
    Mejores_modelos =
      mejores_modelos,

    Comparacion_modelos =
      tabla_comparacion,

    Sensibilidad_CE =
      sensibilidad_ce,

    Metricas_LOOCV =
      tabla_metricas,

    Ajuste_completo =
      tabla_ajuste,

    Coeficientes =
      tabla_coeficientes,

    Predicciones_LOOCV =
      tabla_predicciones,

    Biblioteca_modelos =
      biblioteca_modelos |>
        mutate(
          Predictores = map_chr(
            Predictores,
            function(x) {
              paste(
                x,
                collapse = " + "
              )
            }
          )
        )
  ),
  path = file.path(
    ruta_tablas,
    "Resultados_modelos_integrados_suelo_cultivo.xlsx"
  )
)


## 15. Guardar objeto para Script 06-D -------------------------

saveRDS(
  list(
    base = base,
    biblioteca_modelos = biblioteca_modelos,
    objetos_modelos = objetos_modelos,
    tabla_metricas = tabla_metricas,
    tabla_ajuste = tabla_ajuste,
    tabla_comparacion = tabla_comparacion,
    tabla_predicciones = tabla_predicciones,
    tabla_coeficientes = tabla_coeficientes,
    mejores_modelos = mejores_modelos,
    sensibilidad_ce = sensibilidad_ce
  ),
  file = file.path(
    ruta_objetos,
    "objetos_modelos_integrados_06C.rds"
  )
)


## 16. Registro de sesión --------------------------------------

capture.output(
  sessionInfo(),
  file = file.path(
    ruta_salida,
    "sessionInfo_06_C.txt"
  )
)


## 17. Resumen en consola --------------------------------------

cat(
  "\n============================================\n"
)

cat(
  "SCRIPT 06-C FINALIZADO\n\n"
)

cat(
  "Respuestas evaluadas:\n",
  paste(
    variables_cultivo,
    collapse = ", "
  ),
  "\n\n"
)

cat(
  "Mejores modelos M2 por respuesta:\n"
)

print(
  mejores_modelos |>
    select(
      Respuesta,
      Modelo_ID,
      RMSE,
      MAE,
      Cor_Pearson,
      Mejora_RMSE_porcentaje,
      Mejora_MAE_porcentaje,
      R2_incremental_LOOCV,
      Estado_modelo
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
      Diferencia_MAE_corregida_menos_original
    )
)

cat(
  "\nObjeto para Script 06-D:\n",
  file.path(
    ruta_objetos,
    "objetos_modelos_integrados_06C.rds"
  ),
  "\n"
)

cat(
  "============================================\n"
)
