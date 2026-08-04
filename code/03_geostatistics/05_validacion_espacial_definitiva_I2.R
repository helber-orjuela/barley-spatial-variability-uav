
## ============================================================
## SCRIPT 05. VALIDACIÓN ESPACIAL DEFINITIVA – SUELO I2
## Artículo cebada: revisión mayor
##
## Entrada:
##   AJUSTES_ARTÍCULO/04_KRIGING_I2/objetos_R/
##   objetos_kriging_I2.rds
##
## También reutiliza:
##   AJUSTES_ARTÍCULO/03_ANISOTROPIA_I2/objetos_R/
##   objetos_anisotropia_I2.rds
##
## Salidas:
##   AJUSTES_ARTÍCULO/05_VALIDACION_ESPACIAL_I2/
##
## Objetivos:
##   1. Ejecutar validación cruzada LOOCV para cada variable.
##   2. Calcular ME, MAE, RMSE, RMSE normalizado, RPD, NSE,
##      Pearson, Spearman y CCC de Lin, sin usar p-valores.
##   3. Generar diagnósticos observado-predicho, residuos,
##      QQ-plot, histograma y mapa de residuos.
##   4. Comparar CE original frente a CE corregida OLS.
##   5. Consolidar ranking de desempeño y texto interpretativo.
##
## NOTAS:
##   - N, P y K son índices instrumentales nominales.
##   - pH se conserva solo como diagnóstico exploratorio.
##   - CE corregida OLS es específica de este conjunto de datos.
##   - La varianza de Kriging no incluye la incertidumbre adicional
##     de la calibración OLS.
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
  "DescTools",
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
library(DescTools)
library(readr)


## 1. Rutas ----------------------------------------------------

ruta_base <- paste0(
  "D:/UDCA_2017_CLASE/46_ARTICULOS_2025/",
  "ARTICULO_CEBADA_DOCTORADO/AJUSTES_ARTÍCULO"
)

ruta_objeto_kriging <- file.path(
  ruta_base,
  "04_KRIGING_I2",
  "objetos_R",
  "objetos_kriging_I2.rds"
)

ruta_objeto_anisotropia <- file.path(
  ruta_base,
  "03_ANISOTROPIA_I2",
  "objetos_R",
  "objetos_anisotropia_I2.rds"
)

ruta_salida <- file.path(
  ruta_base,
  "05_VALIDACION_ESPACIAL_I2"
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
  "datos_validacion"
)

if (!file.exists(ruta_objeto_kriging)) {
  stop(
    "No existe el objeto del Script 04:\n",
    ruta_objeto_kriging
  )
}

if (!file.exists(ruta_objeto_anisotropia)) {
  stop(
    "No existe el objeto del Script 03:\n",
    ruta_objeto_anisotropia
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
  ruta_datos,
  recursive = TRUE,
  showWarnings = FALSE
)


## 2. Cargar objetos -------------------------------------------

obj_k <- readRDS(
  ruta_objeto_kriging
)

obj_a <- readRDS(
  ruta_objeto_anisotropia
)

suelo_sf <- obj_k$suelo_sf
lote <- obj_k$lote
resultados_kriging <- obj_k$resultados_kriging
tabla_resumen_superficies <- obj_k$tabla_resumen_superficies
tabla_modelos_finales <- obj_k$tabla_modelos_finales
clasificacion_variables <- obj_k$clasificacion_variables

modelos_finales <- obj_a$modelos_finales
tabla_comparacion_iso_anis <- obj_a$tabla_comparacion_iso_anis

variables_modelar <- intersect(
  names(modelos_finales),
  names(resultados_kriging)
)

if (length(variables_modelar) == 0) {
  stop(
    "No se encontraron variables comunes entre los modelos y ",
    "las superficies de Kriging."
  )
}

lote <- st_transform(
  lote,
  st_crs(suelo_sf)
)


## 3. Funciones métricas ---------------------------------------

calcular_ccc <- function(observado, predicho) {

  resultado <- try(
    DescTools::CCC(
      observado,
      predicho,
      ci = "z-transform"
    ),
    silent = TRUE
  )

  if (inherits(resultado, "try-error")) {
    return(
      tibble(
        CCC = NA_real_,
        CCC_LI95 = NA_real_,
        CCC_LS95 = NA_real_
      )
    )
  }

  tibble(
    CCC = as.numeric(
      resultado$rho.c[1]
    ),
    CCC_LI95 = as.numeric(
      resultado$rho.c[2]
    ),
    CCC_LS95 = as.numeric(
      resultado$rho.c[3]
    )
  )
}


calcular_metricas <- function(observado, predicho, var_prediccion) {

  residuo <- observado - predicho

  validos <- is.finite(observado) &
    is.finite(predicho)

  observado_v <- observado[validos]
  predicho_v <- predicho[validos]
  residuo_v <- residuo[validos]

  if (length(observado_v) < 3) {
    return(
      tibble(
        n = length(observado_v),
        ME = NA_real_,
        MAE = NA_real_,
        RMSE = NA_real_,
        NRMSE_media_porcentaje = NA_real_,
        NRMSE_rango_porcentaje = NA_real_,
        RPD = NA_real_,
        NSE = NA_real_,
        Pearson = NA_real_,
        Spearman = NA_real_,
        RMSSE = NA_real_,
        MSNE = NA_real_
      )
    )
  }

  media_obs <- mean(
    observado_v
  )

  rango_obs <- diff(
    range(
      observado_v
    )
  )

  de_obs <- sd(
    observado_v
  )

  rmse <- sqrt(
    mean(
      residuo_v^2
    )
  )

  denominador_nse <- sum(
    (
      observado_v -
        media_obs
    )^2
  )

  nse <- ifelse(
    denominador_nse > 0,
    1 -
      sum(
        residuo_v^2
      ) /
      denominador_nse,
    NA_real_
  )

  var_valida <- validos &
    is.finite(var_prediccion) &
    var_prediccion > 0

  residuo_est <- rep(
    NA_real_,
    length(observado)
  )

  residuo_est[var_valida] <-
    residuo[var_valida] /
    sqrt(
      var_prediccion[var_valida]
    )

  tibble(
    n = length(observado_v),
    ME = mean(
      residuo_v
    ),
    MAE = mean(
      abs(
        residuo_v
      )
    ),
    RMSE = rmse,
    NRMSE_media_porcentaje = ifelse(
      media_obs != 0,
      100 * rmse / abs(media_obs),
      NA_real_
    ),
    NRMSE_rango_porcentaje = ifelse(
      rango_obs > 0,
      100 * rmse / rango_obs,
      NA_real_
    ),
    RPD = ifelse(
      rmse > 0,
      de_obs / rmse,
      NA_real_
    ),
    NSE = nse,
    Pearson = suppressWarnings(
      cor(
        observado_v,
        predicho_v,
        method = "pearson"
      )
    ),
    Spearman = suppressWarnings(
      cor(
        observado_v,
        predicho_v,
        method = "spearman"
      )
    ),
    RMSSE = ifelse(
      any(var_valida),
      sqrt(
        mean(
          residuo_est^2,
          na.rm = TRUE
        )
      ),
      NA_real_
    ),
    MSNE = ifelse(
      any(var_valida),
      mean(
        residuo_est^2,
        na.rm = TRUE
      ),
      NA_real_
    )
  )
}


clasificar_desempeno <- function(
  nrmse_rango,
  rpd,
  nse
) {

  if (
    !is.finite(nrmse_rango) ||
      !is.finite(rpd)
  ) {
    return(
      "No clasificable"
    )
  }

  if (
    nrmse_rango <= 20 &&
      rpd >= 1.5 &&
      is.finite(nse) &&
      nse > 0
  ) {
    return(
      "Aceptable exploratorio"
    )
  }

  if (
    nrmse_rango <= 30 &&
      rpd >= 1.0
  ) {
    return(
      "Limitado"
    )
  }

  "Deficiente"
}


## 4. Validación cruzada LOOCV -------------------------------

validaciones <- list()
metricas_variables <- list()
errores_validacion <- list()

for (variable in variables_modelar) {

  message(
    "Validación LOOCV de ",
    variable
  )

  modelo <- modelos_finales[[variable]]

  formula_variable <- as.formula(
    paste(
      variable,
      "~ 1"
    )
  )

  cv <- try(
    krige.cv(
      formula = formula_variable,
      locations = suelo_sf,
      model = modelo,
      nfold = nrow(suelo_sf),
      verbose = FALSE
    ),
    silent = TRUE
  )

  if (inherits(cv, "try-error")) {

    errores_validacion[[variable]] <- tibble(
      Variable = variable,
      Error = as.character(cv)
    )

    next
  }

  coords <- st_coordinates(
    cv
  )

  cv_df <- cv |>
    st_drop_geometry() |>
    mutate(
      Nodo = suelo_sf$Nodo,
      X = coords[, 1],
      Y = coords[, 2],
      Variable = variable,
      Observado = observed,
      Predicho = var1.pred,
      Residuo = observed - var1.pred,
      Residuo_absoluto = abs(
        observed - var1.pred
      ),
      Residuo_estandarizado = ifelse(
        is.finite(var1.var) &
          var1.var > 0,
        (
          observed -
            var1.pred
        ) /
          sqrt(
            var1.var
          ),
        NA_real_
      ),
      .before = 1
    )

  metricas <- calcular_metricas(
    observado = cv_df$Observado,
    predicho = cv_df$Predicho,
    var_prediccion = cv_df$var1.var
  )

  ccc <- calcular_ccc(
    cv_df$Observado,
    cv_df$Predicho
  )

  fila_modelo <- tabla_modelos_finales |>
    filter(
      Variable == variable
    )

  fila_superficie <- tabla_resumen_superficies |>
    filter(
      Variable == variable
    )

  resumen <- bind_cols(
    tibble(
      Variable = variable,
      Modelo_final = ifelse(
        nrow(fila_modelo) > 0,
        fila_modelo$Modelo_codigo[1],
        NA_character_
      ),
      Tipo_modelo_final = ifelse(
        nrow(fila_modelo) > 0,
        fila_modelo$Tipo_modelo_final[1],
        NA_character_
      ),
      Angulo = ifelse(
        nrow(fila_modelo) > 0,
        fila_modelo$Angulo_final[1],
        NA_real_
      ),
      Relacion_anisotropia = ifelse(
        nrow(fila_modelo) > 0,
        fila_modelo$Relacion_anisotropia_final[1],
        NA_real_
      )
    ),
    metricas,
    ccc
  ) |>
    mutate(
      Clasificacion_desempeno =
        clasificar_desempeno(
          NRMSE_rango_porcentaje,
          RPD,
          NSE
        ),
      Error_estandar_kriging_superficie = ifelse(
        nrow(fila_superficie) > 0,
        fila_superficie$Error_estandar_medio[1],
        NA_real_
      )
    )

  validaciones[[variable]] <- cv_df
  metricas_variables[[variable]] <- resumen

  write_csv(
    cv_df,
    file.path(
      ruta_datos,
      paste0(
        "Validacion_LOOCV_",
        variable,
        "_I2.csv"
      )
    )
  )
}


tabla_metricas <- bind_rows(
  metricas_variables
)

tabla_errores_validacion <- bind_rows(
  errores_validacion
)

if (nrow(tabla_metricas) == 0) {
  stop(
    "No se obtuvo ninguna validación válida."
  )
}

print(
  tabla_metricas
)


## 5. Ranking de desempeño ------------------------------------

tabla_ranking <- tabla_metricas |>
  mutate(
    Rank_RMSE_relativo = min_rank(
      NRMSE_rango_porcentaje
    ),
    Rank_RPD = min_rank(
      desc(
        RPD
      )
    ),
    Rank_NSE = min_rank(
      desc(
        NSE
      )
    ),
    Rank_CCC = min_rank(
      desc(
        CCC
      )
    ),
    Rank_MSNE = min_rank(
      abs(
        MSNE - 1
      )
    ),
    Score_global =
      0.30 * Rank_RMSE_relativo +
      0.20 * Rank_RPD +
      0.15 * Rank_NSE +
      0.15 * Rank_CCC +
      0.20 * Rank_MSNE
  ) |>
  arrange(
    Score_global
  ) |>
  mutate(
    Posicion = row_number(),
    .before = 1
  )

print(
  tabla_ranking
)


## 6. Comparación CE original vs corregida ---------------------

if (
  all(
    c(
      "CE",
      "CE_corregida_OLS"
    ) %in% tabla_metricas$Variable
  )
) {

  met_ce <- tabla_metricas |>
    filter(
      Variable == "CE"
    )

  met_ce_corr <- tabla_metricas |>
    filter(
      Variable == "CE_corregida_OLS"
    )

  tabla_comparacion_ce <- tibble(
    Indicador = c(
      "ME",
      "MAE",
      "RMSE",
      "NRMSE_media_porcentaje",
      "NRMSE_rango_porcentaje",
      "RPD",
      "NSE",
      "Pearson",
      "Spearman",
      "CCC",
      "RMSSE",
      "MSNE"
    ),
    CE_original = c(
      met_ce$ME,
      met_ce$MAE,
      met_ce$RMSE,
      met_ce$NRMSE_media_porcentaje,
      met_ce$NRMSE_rango_porcentaje,
      met_ce$RPD,
      met_ce$NSE,
      met_ce$Pearson,
      met_ce$Spearman,
      met_ce$CCC,
      met_ce$RMSSE,
      met_ce$MSNE
    ),
    CE_corregida = c(
      met_ce_corr$ME,
      met_ce_corr$MAE,
      met_ce_corr$RMSE,
      met_ce_corr$NRMSE_media_porcentaje,
      met_ce_corr$NRMSE_rango_porcentaje,
      met_ce_corr$RPD,
      met_ce_corr$NSE,
      met_ce_corr$Pearson,
      met_ce_corr$Spearman,
      met_ce_corr$CCC,
      met_ce_corr$RMSSE,
      met_ce_corr$MSNE
    )
  ) |>
    mutate(
      Diferencia =
        CE_corregida -
        CE_original,
      Cambio_porcentual = ifelse(
        CE_original != 0,
        100 *
          (
            CE_corregida -
              CE_original
          ) /
          abs(
            CE_original
          ),
        NA_real_
      )
    )

  mejora_rmse_ce <- 100 * (
    met_ce$RMSE -
      met_ce_corr$RMSE
  ) / met_ce$RMSE

  mejora_mae_ce <- 100 * (
    met_ce$MAE -
      met_ce_corr$MAE
  ) / met_ce$MAE

  conservacion_rango <- tabla_comparacion_iso_anis |>
    filter(
      Variable %in% c(
        "CE",
        "CE_corregida_OLS"
      )
    )

} else {

  tabla_comparacion_ce <- tibble(
    Indicador = "Comparación no disponible",
    CE_original = NA_real_,
    CE_corregida = NA_real_,
    Diferencia = NA_real_,
    Cambio_porcentual = NA_real_
  )

  mejora_rmse_ce <- NA_real_
  mejora_mae_ce <- NA_real_
  conservacion_rango <- tibble()
}


## 7. Diagnósticos gráficos ------------------------------------

lote_borde <- st_boundary(
  st_union(
    lote
  )
)

crear_graficos_validacion <- function(
  variable,
  datos
) {

  rango_total <- range(
    c(
      datos$Observado,
      datos$Predicho
    ),
    na.rm = TRUE
  )

  g_obs_pred <- ggplot(
    datos,
    aes(
      x = Observado,
      y = Predicho
    )
  ) +
    geom_point(
      size = 2.4
    ) +
    geom_abline(
      intercept = 0,
      slope = 1,
      linetype = "dashed"
    ) +
    coord_equal(
      xlim = rango_total,
      ylim = rango_total
    ) +
    labs(
      title =
        "Observado vs predicho",
      x =
        "Observado",
      y =
        "Predicho LOOCV"
    ) +
    theme_minimal()

  g_res_pred <- ggplot(
    datos,
    aes(
      x = Predicho,
      y = Residuo
    )
  ) +
    geom_point(
      size = 2.4
    ) +
    geom_hline(
      yintercept = 0,
      linetype = "dashed"
    ) +
    labs(
      title =
        "Residuos vs predicción",
      x =
        "Predicho LOOCV",
      y =
        "Residuo"
    ) +
    theme_minimal()

  g_hist <- ggplot(
    datos,
    aes(
      x = Residuo
    )
  ) +
    geom_histogram(
      bins = 8,
      color = "white"
    ) +
    geom_vline(
      xintercept = 0,
      linetype = "dashed"
    ) +
    labs(
      title =
        "Distribución de residuos",
      x =
        "Residuo",
      y =
        "Frecuencia"
    ) +
    theme_minimal()

  g_qq <- ggplot(
    datos,
    aes(
      sample = Residuo
    )
  ) +
    stat_qq() +
    stat_qq_line() +
    labs(
      title =
        "QQ-plot de residuos",
      x =
        "Cuantil teórico",
      y =
        "Cuantil observado"
    ) +
    theme_minimal()

  g_mapa <- ggplot() +
    geom_sf(
      data = lote_borde,
      fill = NA,
      linewidth = 0.45
    ) +
    geom_point(
      data = datos,
      aes(
        x = X,
        y = Y,
        size = Residuo_absoluto,
        fill = Residuo
      ),
      shape = 21
    ) +
    coord_sf(
      expand = FALSE
    ) +
    labs(
      title =
        "Distribución espacial de residuos",
      x =
        "X (m)",
      y =
        "Y (m)",
      size =
        "|Residuo|",
      fill =
        "Residuo"
    ) +
    theme_minimal()

  panel <- (
    g_obs_pred |
      g_res_pred
  ) / (
    g_hist |
      g_qq
  ) / g_mapa +
    plot_annotation(
      title = paste0(
        "Validación espacial LOOCV de ",
        variable
      )
    )

  list(
    panel = panel,
    observado_predicho = g_obs_pred,
    residuos_prediccion = g_res_pred,
    histograma = g_hist,
    qq = g_qq,
    mapa = g_mapa
  )
}


graficos_validacion <- list()

for (variable in names(validaciones)) {

  graficos <- crear_graficos_validacion(
    variable,
    validaciones[[variable]]
  )

  graficos_validacion[[variable]] <- graficos

  ggsave(
    filename = file.path(
      ruta_figuras,
      paste0(
        "Validacion_LOOCV_",
        variable,
        "_I2.png"
      )
    ),
    plot = graficos$panel,
    width = 12,
    height = 15,
    dpi = 400
  )
}


## 8. Panel comparativo observado-predicho ---------------------

graficos_obs_pred <- map(
  graficos_validacion,
  "observado_predicho"
)

panel_obs_pred <- wrap_plots(
  graficos_obs_pred,
  ncol = 2
) +
  plot_annotation(
    title =
      "Comparación observado-predicho mediante LOOCV en I2"
  )

ggsave(
  filename = file.path(
    ruta_figuras,
    "Panel_observado_predicho_LOOCV_I2.png"
  ),
  plot = panel_obs_pred,
  width = 14,
  height = 18,
  dpi = 400
)


## 9. Panel de mapas de residuos -------------------------------

graficos_mapas <- map(
  graficos_validacion,
  "mapa"
)

panel_mapas_residuos <- wrap_plots(
  graficos_mapas,
  ncol = 2
) +
  plot_annotation(
    title =
      "Distribución espacial de los residuos LOOCV en I2"
  )

ggsave(
  filename = file.path(
    ruta_figuras,
    "Panel_mapas_residuos_LOOCV_I2.png"
  ),
  plot = panel_mapas_residuos,
  width = 14,
  height = 18,
  dpi = 400
)


## 10. Interpretación automática -------------------------------

lineas_interpretacion <- c(
  "VALIDACIÓN ESPACIAL EXPLORATORIA DE SUELO – I2",
  "================================================",
  "",
  paste0(
    "Se validaron ",
    nrow(tabla_metricas),
    " variables mediante leave-one-out cross-validation (LOOCV)."
  ),
  paste0(
    "No se utilizaron p-valores como criterio de desempeño. ",
    "La interpretación se fundamentó en la magnitud del error, ",
    "la concordancia, la correlación, el RPD, el NSE y la ",
    "calibración de la incertidumbre mediante MSNE."
  ),
  ""
)

for (i in seq_len(nrow(tabla_ranking))) {

  fila <- tabla_ranking[i, ]

  linea <- paste0(
    fila$Variable,
    ": RMSE=",
    round(
      fila$RMSE,
      3
    ),
    "; MAE=",
    round(
      fila$MAE,
      3
    ),
    "; NRMSE(rango)=",
    round(
      fila$NRMSE_rango_porcentaje,
      1
    ),
    "%; RPD=",
    round(
      fila$RPD,
      2
    ),
    "; NSE=",
    round(
      fila$NSE,
      2
    ),
    "; Spearman=",
    round(
      fila$Spearman,
      2
    ),
    "; CCC=",
    round(
      fila$CCC,
      2
    ),
    "; MSNE=",
    round(
      fila$MSNE,
      2
    ),
    ". Clasificación: ",
    fila$Clasificacion_desempeno,
    "."
  )

  lineas_interpretacion <- c(
    lineas_interpretacion,
    linea
  )
}

lineas_interpretacion <- c(
  lineas_interpretacion,
  "",
  "COMPARACIÓN DE CE",
  "-----------------"
)

if (
  is.finite(
    mejora_rmse_ce
  )
) {

  lineas_interpretacion <- c(
    lineas_interpretacion,
    paste0(
      "La CE corregida modificó el RMSE LOOCV en ",
      round(
        mejora_rmse_ce,
        1
      ),
      "% respecto a la CE original."
    ),
    paste0(
      "La modificación del MAE fue de ",
      round(
        mejora_mae_ce,
        1
      ),
      "%."
    ),
    paste0(
      "Como la corrección OLS es lineal, la geometría del patrón ",
      "espacial se conserva; cambian principalmente la escala, ",
      "la magnitud del error y la varianza expresada en las nuevas unidades."
    )
  )

} else {

  lineas_interpretacion <- c(
    lineas_interpretacion,
    "No fue posible realizar la comparación completa entre CE original y corregida."
  )
}

lineas_interpretacion <- c(
  lineas_interpretacion,
  "",
  "ADVERTENCIAS",
  "------------",
  "N, P y K deben interpretarse como índices instrumentales nominales y no como determinaciones químicas independientes.",
  "El pH se conserva como diagnóstico exploratorio debido a su falta de concordancia frente al laboratorio.",
  "La incertidumbre de la calibración OLS no está incorporada dentro de la varianza de Kriging.",
  "Los resultados son específicos del lote, el sensor, la campaña I2 y el tamaño muestral disponible."
)

writeLines(
  lineas_interpretacion,
  con = file.path(
    ruta_salida,
    "Interpretacion_validacion_espacial_I2.txt"
  )
)


## 11. Exportación de tablas -----------------------------------

write_xlsx(
  list(
    Metricas_LOOCV =
      tabla_metricas,

    Ranking_desempeno =
      tabla_ranking,

    Comparacion_CE =
      tabla_comparacion_ce,

    Modelos_finales =
      tabla_modelos_finales,

    Comparacion_iso_anis =
      tabla_comparacion_iso_anis,

    Resumen_superficies =
      tabla_resumen_superficies,

    Errores_validacion =
      tabla_errores_validacion
  ),
  path = file.path(
    ruta_tablas,
    "Resultados_validacion_espacial_I2.xlsx"
  )
)


## 12. Base consolidada de validación --------------------------

tabla_validacion_completa <- bind_rows(
  validaciones
)

write_csv(
  tabla_validacion_completa,
  file.path(
    ruta_datos,
    "Validacion_LOOCV_todas_variables_I2.csv"
  )
)


## 13. Registro de sesión --------------------------------------

capture.output(
  sessionInfo(),
  file = file.path(
    ruta_salida,
    "sessionInfo_Script05.txt"
  )
)


## 14. Resumen en consola --------------------------------------

cat(
  "\n============================================\n"
)

cat(
  "SCRIPT 05 FINALIZADO\n\n"
)

cat(
  "Ranking de desempeño:\n"
)

print(
  tabla_ranking |>
    select(
      Posicion,
      Variable,
      RMSE,
      NRMSE_rango_porcentaje,
      RPD,
      NSE,
      Spearman,
      CCC,
      MSNE,
      Clasificacion_desempeno
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
