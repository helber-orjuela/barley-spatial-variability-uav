
## ============================================================
## 06_SUELO_CULTIVO – SCRIPT 06-E
## SÍNTESIS FINAL PARA PUBLICACIÓN
## ============================================================
##
## Entradas:
##   06_SUELO_CULTIVO/03_MODELOS_INTEGRADOS/objetos_R/
##   objetos_modelos_integrados_06C.rds
##
##   06_SUELO_CULTIVO/04_VALIDACION_SENSIBILIDAD/objetos_R/
##   objetos_validacion_sensibilidad_06D.rds
##
## Salida:
##   06_SUELO_CULTIVO/05_SINTESIS_PUBLICACION/
##
## PROPÓSITO
##   1. Integrar desempeño LOOCV, estructura espacial residual,
##      influencia y sensibilidad CE.
##   2. Separar evidencia predictiva, evidencia espacial y
##      advertencias metodológicas.
##   3. Generar tablas compactas para manuscrito y suplemento.
##   4. Crear figuras finales de alta resolución.
##   5. Evitar interpretar como válida una validación por repetición
##      que no produjo predicciones evaluables.
##
## NOTA METODOLÓGICA IMPORTANTE
##   En los resultados de 06-D, la validación dejando una repetición
##   completa por fuera produjo n = 0 para todos los modelos. Esto
##   ocurre porque Replication forma parte de la fórmula y el nivel
##   excluido no está presente durante el ajuste. Por tanto:
##
##   - LORO no se utilizará para seleccionar modelos.
##   - Se reportará como "no evaluable con esta parametrización".
##   - La síntesis se apoyará en LOOCV, Moran residual, influencia
##     y sensibilidad entre escenarios.
##
## PRINCIPIOS DE INTERPRETACIÓN
##   - No se usa el p-valor como criterio principal.
##   - CE original y corregida OLS son linealmente equivalentes
##     para correlación y desempeño cuando se usan como predictor.
##   - N, P y K son índices instrumentales nominales.
##   - pH y Hv se conservan como variables exploratorias.
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


## 1. Rutas ----------------------------------------------------

ruta_base <- paste0(
  "D:/UDCA_2017_CLASE/46_ARTICULOS_2025/",
  "ARTICULO_CEBADA_DOCTORADO/AJUSTES_ARTÍCULO"
)

ruta_06C <- file.path(
  ruta_base,
  "06_SUELO_CULTIVO",
  "03_MODELOS_INTEGRADOS",
  "objetos_R",
  "objetos_modelos_integrados_06C.rds"
)

ruta_06D <- file.path(
  ruta_base,
  "06_SUELO_CULTIVO",
  "04_VALIDACION_SENSIBILIDAD",
  "objetos_R",
  "objetos_validacion_sensibilidad_06D.rds"
)

ruta_salida <- file.path(
  ruta_base,
  "06_SUELO_CULTIVO",
  "05_SINTESIS_PUBLICACION"
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

faltantes <- c(
  ruta_06C,
  ruta_06D
)[
  !file.exists(
    c(
      ruta_06C,
      ruta_06D
    )
  )
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

obj_C <- readRDS(
  ruta_06C
)

obj_D <- readRDS(
  ruta_06D
)

tabla_comparacion <- obj_C$tabla_comparacion
tabla_predicciones <- obj_C$tabla_predicciones
tabla_coeficientes <- obj_C$tabla_coeficientes
mejores_modelos <- obj_C$mejores_modelos
sensibilidad_ce <- obj_D$sensibilidad_ce
tabla_moran <- obj_D$tabla_moran_residuos
tabla_variogramas <- obj_D$tabla_variogramas
tabla_influencia <- obj_D$tabla_influencia
resumen_influencia <- obj_D$resumen_influencia
tabla_loro_metricas <- obj_D$tabla_loro_metricas
tabla_loro_por_repeticion <- obj_D$tabla_loro_por_repeticion


## 3. Validar estructura mínima --------------------------------

objetos_requeridos <- list(
  tabla_comparacion = tabla_comparacion,
  mejores_modelos = mejores_modelos,
  tabla_moran = tabla_moran,
  resumen_influencia = resumen_influencia,
  sensibilidad_ce = sensibilidad_ce
)

nulos <- names(
  objetos_requeridos
)[
  map_lgl(
    objetos_requeridos,
    is.null
  )
]

if (length(nulos) > 0) {
  stop(
    "Faltan objetos requeridos: ",
    paste(
      nulos,
      collapse = ", "
    )
  )
}


## 4. Diagnóstico de validación por repetición -----------------

resumen_loro_estado <- tabla_loro_metricas |>
  summarise(
    n_modelos = n(),
    n_modelos_evaluables = sum(
      n >= 4 &
        is.finite(
          RMSE
        )
    ),
    n_modelos_no_evaluables = sum(
      n < 4 |
        !is.finite(
          RMSE
        )
    )
  ) |>
  mutate(
    Estado_LORO = case_when(
      n_modelos_evaluables == 0 ~
        "No evaluable con la parametrización actual",
      n_modelos_no_evaluables > 0 ~
        "Parcialmente evaluable",
      TRUE ~
        "Evaluable"
    ),
    Explicacion = case_when(
      n_modelos_evaluables == 0 ~
        paste0(
          "Replication forma parte del modelo y el nivel excluido ",
          "no está presente en el entrenamiento; por ello no se ",
          "obtuvieron predicciones válidas."
        ),
      TRUE ~
        "Se obtuvieron predicciones para al menos algunos modelos."
    )
  )


## 5. Seleccionar referencia M1 y mejores M2 -------------------

tabla_M1 <- tabla_comparacion |>
  filter(
    Modelo_ID == "M1_Tratamientos"
  ) |>
  select(
    Respuesta,
    RMSE_M1 = RMSE,
    MAE_M1 = MAE,
    Sesgo_M1 = Sesgo,
    Cor_Pearson_M1 = Cor_Pearson,
    Cor_Spearman_M1 = Cor_Spearman
  )

tabla_M2 <- mejores_modelos |>
  select(
    Respuesta,
    Modelo_ID_M2 = Modelo_ID,
    RMSE_M2 = RMSE,
    MAE_M2 = MAE,
    Sesgo_M2 = Sesgo,
    Cor_Pearson_M2 = Cor_Pearson,
    Cor_Spearman_M2 = Cor_Spearman,
    Mejora_RMSE_porcentaje,
    Mejora_MAE_porcentaje,
    R2_incremental_LOOCV,
    Estado_modelo,
    Evidencia_integrada
  )


## 6. Integrar Moran residual ----------------------------------

moran_ancho <- tabla_moran |>
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


## 7. Integrar influencia --------------------------------------

influencia_ancha <- resumen_influencia |>
  select(
    Respuesta,
    Tipo_modelo,
    n_Cook,
    n_Leverage,
    n_Residuos_extremos,
    Cook_maximo,
    Leverage_maximo
  ) |>
  pivot_wider(
    names_from = Tipo_modelo,
    values_from = c(
      n_Cook,
      n_Leverage,
      n_Residuos_extremos,
      Cook_maximo,
      Leverage_maximo
    )
  )


## 8. Tabla integrada final ------------------------------------

tabla_final <- tabla_M1 |>
  left_join(
    tabla_M2,
    by = "Respuesta"
  ) |>
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
    Mejora_error_conjunta = (
      Mejora_RMSE_porcentaje > 0
    ) &
      (
        Mejora_MAE_porcentaje > 0
      ),
    Reduce_autocorrelacion = Reduccion_abs_Moran > 0,
    Aumenta_Cook = n_Cook_Mejor_M2 > n_Cook_M1,
    Aumenta_residuos_extremos =
      n_Residuos_extremos_Mejor_M2 >
      n_Residuos_extremos_M1,
    Robustez_publicacion = case_when(
      Mejora_error_conjunta &
        Reduce_autocorrelacion &
        !Aumenta_residuos_extremos ~
        "Moderada-alta",
      Mejora_error_conjunta &
        Reduce_autocorrelacion ~
        "Moderada",
      Mejora_error_conjunta ~
        "Limitada",
      TRUE ~
        "Sin evidencia de mejora"
    ),
    Interpretacion = case_when(
      Respuesta == "NDVI" ~
        paste0(
          "El modelo con índice K y temperatura redujo RMSE y MAE, ",
          "pero mantuvo autocorrelación residual positiva; la mejora ",
          "predictiva debe considerarse exploratoria."
        ),
      Respuesta == "Potasio_savia" ~
        paste0(
          "El modelo con índice P y temperatura redujo RMSE y MAE ",
          "y disminuyó la magnitud de Moran residual; representa la ",
          "señal integrada más consistente, aunque de magnitud limitada."
        ),
      Mejora_error_conjunta ~
        paste0(
          "El modelo M2 mejoró ambos errores, pero requiere revisión ",
          "de estructura residual e influencia."
        ),
      TRUE ~
        paste0(
          "Las covariables de suelo no mejoraron simultáneamente ",
          "RMSE y MAE frente al modelo de tratamientos."
        )
    ),
    Advertencia_modelo = case_when(
      str_detect(
        Modelo_ID_M2,
        "Indice_"
      ) ~
        "Índice instrumental nominal",
      str_detect(
        Modelo_ID_M2,
        "CE_corregida"
      ) ~
        "CE corregida linealmente equivalente a CE original",
      str_detect(
        Modelo_ID_M2,
        "Exploratorio"
      ) ~
        "pH y Hv con alcance exploratorio",
      TRUE ~
        "Sin advertencia adicional"
    )
  ) |>
  arrange(
    desc(
      Mejora_error_conjunta
    ),
    desc(
      Mejora_RMSE_porcentaje
    )
  )


## 9. Sensibilidad CE ------------------------------------------

tabla_ce_compacta <- sensibilidad_ce |>
  select(
    Respuesta,
    Diferencia_RMSE_corregida_menos_original,
    Diferencia_MAE_corregida_menos_original,
    Equivalencia_RMSE,
    Equivalencia_MAE,
    Interpretacion
  )


## 10. Tabla compacta para manuscrito --------------------------

tabla_manuscrito <- tabla_final |>
  select(
    Respuesta,
    Modelo_ID_M2,
    RMSE_M1,
    RMSE_M2,
    Mejora_RMSE_porcentaje,
    MAE_M1,
    MAE_M2,
    Mejora_MAE_porcentaje,
    Cor_Pearson_M2,
    R2_incremental_LOOCV,
    Moran_I_M1,
    Moran_I_Mejor_M2,
    Clasificacion_Mejor_M2,
    n_Cook_Mejor_M2,
    n_Residuos_extremos_Mejor_M2,
    Robustez_publicacion,
    Advertencia_modelo
  )


## 11. Textos base prudentes -----------------------------------

texto_general <- tibble(
  Seccion = c(
    "Aporte general del suelo",
    "NDVI",
    "Potasio en savia",
    "CE original y corregida",
    "Validación por repetición",
    "Alcance metodológico"
  ),
  Texto = c(
    paste0(
      "The addition of kriged soil covariates did not improve all ",
      "crop responses. Improvements were selective and generally ",
      "small, indicating that the treatment structure remained the ",
      "main source of predictive information."
    ),
    paste0(
      "For NDVI, the model including the nominal soil K index and ",
      "soil temperature reduced LOOCV RMSE and MAE. However, positive ",
      "residual spatial autocorrelation remained, so the improvement ",
      "should be interpreted as exploratory rather than definitive."
    ),
    paste0(
      "For sap potassium, the model including the nominal soil P ",
      "index and soil temperature reduced both LOOCV errors and the ",
      "absolute magnitude of residual Moran's I, providing the most ",
      "internally consistent integrated signal."
    ),
    paste0(
      "Models based on original and OLS-corrected EC were numerically ",
      "equivalent because the correction was a linear transformation. ",
      "The correction changes measurement scale but not spatial order."
    ),
    paste0(
      "Leave-one-replication-out validation was not evaluable under ",
      "the current parameterization because replication was included ",
      "as a fixed predictor and the excluded level was absent from the ",
      "training data."
    ),
    paste0(
      "Soil N, P, and K values are nominal instrumental indices and ",
      "must not be interpreted as validated nutrient concentrations. ",
      "All soil covariates at crop nodes were obtained from kriging."
    )
  )
)


## 12. Figuras finales -----------------------------------------

fig_metricas <- tabla_final |>
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
      Metrica == "Mejora_RMSE_porcentaje" ~
        "RMSE",
      Metrica == "Mejora_MAE_porcentaje" ~
        "MAE",
      TRUE ~
        Metrica
    )
  )

g_metricas <- ggplot(
  fig_metricas,
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
    title = "Aporte predictivo de las covariables de suelo",
    subtitle = "Valores positivos indican reducción del error frente a M1",
    x = NULL,
    y = "Mejora (%)",
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
    "01_mejora_predictiva_final.png"
  ),
  plot = g_metricas,
  width = 11,
  height = 7,
  dpi = 500
)


fig_moran <- tabla_final |>
  select(
    Respuesta,
    Moran_I_M1,
    Moran_I_Mejor_M2
  ) |>
  pivot_longer(
    cols = c(
      Moran_I_M1,
      Moran_I_Mejor_M2
    ),
    names_to = "Modelo",
    values_to = "Moran_I"
  ) |>
  mutate(
    Modelo = case_when(
      Modelo == "Moran_I_M1" ~
        "M1 tratamientos",
      Modelo == "Moran_I_Mejor_M2" ~
        "Mejor M2",
      TRUE ~
        Modelo
    )
  )

g_moran <- ggplot(
  fig_moran,
  aes(
    x = Respuesta,
    y = Moran_I,
    shape = Modelo
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
    title = "Autocorrelación espacial residual",
    x = NULL,
    y = "Moran I",
    shape = "Modelo"
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
    "02_Moran_residual_final.png"
  ),
  plot = g_moran,
  width = 11,
  height = 7,
  dpi = 500
)


fig_influencia <- tabla_final |>
  select(
    Respuesta,
    n_Cook_M1,
    n_Cook_Mejor_M2,
    n_Residuos_extremos_M1,
    n_Residuos_extremos_Mejor_M2
  ) |>
  pivot_longer(
    cols = -Respuesta,
    names_to = "Indicador",
    values_to = "Nodos"
  ) |>
  mutate(
    Modelo = case_when(
      str_detect(
        Indicador,
        "Mejor_M2"
      ) ~
        "Mejor M2",
      TRUE ~
        "M1"
    ),
    Diagnostico = case_when(
      str_detect(
        Indicador,
        "Cook"
      ) ~
        "Cook",
      TRUE ~
        "Residuo extremo"
    )
  )

g_influencia <- ggplot(
  fig_influencia,
  aes(
    x = Respuesta,
    y = Nodos,
    shape = Modelo
  )
) +
  geom_point(
    size = 3,
    position = position_dodge(
      width = 0.35
    )
  ) +
  facet_wrap(
    ~ Diagnostico,
    scales = "free_y"
  ) +
  labs(
    title = "Diagnóstico de influencia",
    x = NULL,
    y = "Número de nodos",
    shape = "Modelo"
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
    "03_influencia_final.png"
  ),
  plot = g_influencia,
  width = 12,
  height = 7,
  dpi = 500
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
    title = "Observado frente a predicho en los mejores modelos M2",
    subtitle = "Predicción obtenida mediante LOOCV",
    x = "Observado",
    y = "Predicho"
  ) +
  theme_minimal()

ggsave(
  filename = file.path(
    ruta_figuras,
    "04_observado_predicho_final.png"
  ),
  plot = g_obs_pred,
  width = 13,
  height = 10,
  dpi = 500
)


panel_final <- (
  g_metricas +
    g_moran
) /
  (
    g_influencia +
      g_obs_pred
  ) +
  plot_annotation(
    title = "Síntesis de modelos integrados suelo–cultivo"
  )

ggsave(
  filename = file.path(
    ruta_figuras,
    "05_panel_sintesis_suelo_cultivo.png"
  ),
  plot = panel_final,
  width = 18,
  height = 16,
  dpi = 500
)


## 13. Exportar resultados -------------------------------------

write_csv(
  tabla_final,
  file.path(
    ruta_datos,
    "Sintesis_final_modelos_suelo_cultivo.csv"
  )
)

write_csv(
  tabla_manuscrito,
  file.path(
    ruta_datos,
    "Tabla_compacta_modelos_para_manuscrito.csv"
  )
)

write_csv(
  texto_general,
  file.path(
    ruta_datos,
    "Textos_base_resultados_discusion.csv"
  )
)

write_csv(
  tabla_ce_compacta,
  file.path(
    ruta_datos,
    "Sensibilidad_CE_compacta.csv"
  )
)


write_xlsx(
  list(
    Resumen_ejecutivo =
      texto_general,

    Sintesis_final =
      tabla_final,

    Tabla_manuscrito =
      tabla_manuscrito,

    Sensibilidad_CE =
      tabla_ce_compacta,

    Estado_LORO =
      resumen_loro_estado,

    Moran_residual =
      tabla_moran,

    Influencia_resumen =
      resumen_influencia,

    Mejores_modelos_06C =
      mejores_modelos,

    Comparacion_completa =
      tabla_comparacion
  ),
  path = file.path(
    ruta_tablas,
    "Sintesis_final_modelos_suelo_cultivo_publicacion.xlsx"
  )
)


## 14. Guardar objeto final ------------------------------------

saveRDS(
  list(
    tabla_final = tabla_final,
    tabla_manuscrito = tabla_manuscrito,
    texto_general = texto_general,
    tabla_ce_compacta = tabla_ce_compacta,
    resumen_loro_estado = resumen_loro_estado
  ),
  file = file.path(
    ruta_objetos,
    "objetos_sintesis_final_suelo_cultivo_06E.rds"
  )
)


## 15. Registro de sesión --------------------------------------

capture.output(
  sessionInfo(),
  file = file.path(
    ruta_salida,
    "sessionInfo_06_E.txt"
  )
)


## 16. Resumen en consola --------------------------------------

cat(
  "\n============================================\n"
)

cat(
  "SCRIPT 06-E FINALIZADO\n\n"
)

cat(
  "Estado de validación por repetición:\n"
)

print(
  resumen_loro_estado
)

cat(
  "\nSíntesis final:\n"
)

print(
  tabla_final |>
    select(
      Respuesta,
      Modelo_ID_M2,
      Mejora_RMSE_porcentaje,
      Mejora_MAE_porcentaje,
      Moran_I_M1,
      Moran_I_Mejor_M2,
      Robustez_publicacion,
      Advertencia_modelo
    )
)

cat(
  "\nObjeto final:\n",
  file.path(
    ruta_objetos,
    "objetos_sintesis_final_suelo_cultivo_06E.rds"
  ),
  "\n"
)

cat(
  "============================================\n"
)
