
## ============================================================
## 04_TRATAMIENTOS – SCRIPT 04
## DIAGNÓSTICO ESPACIAL DE RESIDUOS Y FIGURAS FINALES
## ============================================================
##
## Entrada:
##   04_TRATAMIENTOS/03_MODELOS_TRATAMIENTOS/objetos_R/
##   objetos_modelos_tratamientos.rds
##
## Salida:
##   04_TRATAMIENTOS/04_DIAGNOSTICO_RESIDUAL/
##
## PROPÓSITO
##   1. Evaluar si los residuos de los modelos seleccionados
##      conservan estructura espacial.
##   2. Calcular Moran I con vecinos k-NN.
##   3. Construir semivariogramas experimentales de residuos.
##   4. Identificar observaciones influyentes.
##   5. Generar figuras y tablas consolidadas para el manuscrito.
##
## INTERPRETACIÓN
##   - Moran I cercano a 0: ausencia de patrón espacial residual claro.
##   - Moran I positivo: residuos semejantes tienden a agruparse.
##   - Moran I negativo: alternancia local de residuos.
##   - El p-valor se exporta únicamente como información secundaria.
##   - La interpretación principal se basa en magnitud, patrón gráfico,
##     intervalos de simulación y estabilidad de los modelos.
## ============================================================


## 0. Paquetes -------------------------------------------------

paquetes <- c(
  "dplyr",
  "tidyr",
  "purrr",
  "stringr",
  "ggplot2",
  "writexl",
  "readr",
  "sf",
  "spdep",
  "gstat",
  "patchwork",
  "broom"
)

instalar <- paquetes[!paquetes %in% rownames(installed.packages())]

if (length(instalar) > 0) {
  install.packages(instalar)
}

library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(ggplot2)
library(writexl)
library(readr)
library(sf)
library(spdep)
library(gstat)
library(patchwork)
library(broom)


## 1. Rutas ----------------------------------------------------

ruta_base <- paste0(
  "D:/UDCA_2017_CLASE/46_ARTICULOS_2025/",
  "ARTICULO_CEBADA_DOCTORADO/AJUSTES_ARTÍCULO"
)

ruta_objeto <- file.path(
  ruta_base,
  "04_TRATAMIENTOS",
  "03_MODELOS_TRATAMIENTOS",
  "objetos_R",
  "objetos_modelos_tratamientos.rds"
)

ruta_salida <- file.path(
  ruta_base,
  "04_TRATAMIENTOS",
  "04_DIAGNOSTICO_RESIDUAL"
)

ruta_tablas <- file.path(ruta_salida, "tablas")
ruta_figuras <- file.path(ruta_salida, "figuras")
ruta_datos <- file.path(ruta_salida, "datos_procesados")
ruta_objetos <- file.path(ruta_salida, "objetos_R")

if (!file.exists(ruta_objeto)) {
  stop(
    "No existe el objeto del Script 03:\n",
    ruta_objeto
  )
}

dir.create(ruta_tablas, recursive = TRUE, showWarnings = FALSE)
dir.create(ruta_figuras, recursive = TRUE, showWarnings = FALSE)
dir.create(ruta_datos, recursive = TRUE, showWarnings = FALSE)
dir.create(ruta_objetos, recursive = TRUE, showWarnings = FALSE)


## 2. Cargar resultados ----------------------------------------

obj <- readRDS(ruta_objeto)

datos <- obj$datos
resultados_modelos <- obj$resultados_modelos
modelos_seleccionados <- obj$modelos_seleccionados
tabla_metricas <- obj$tabla_metricas
tabla_tamanos_efecto <- obj$tabla_tamanos_efecto
tabla_diagnostico <- obj$tabla_diagnostico
catalogo_respuestas <- obj$catalogo_respuestas

columnas_espaciales <- c("Nodo", "X", "Y")

faltantes_espaciales <- setdiff(
  columnas_espaciales,
  names(datos)
)

if (length(faltantes_espaciales) > 0) {
  stop(
    "Faltan columnas espaciales en el objeto: ",
    paste(faltantes_espaciales, collapse = ", ")
  )
}

tabla_diagnostico <- tabla_diagnostico |>
  left_join(
    datos |>
      select(Nodo, X, Y),
    by = "Nodo"
  )

variables_validas <- intersect(
  modelos_seleccionados$Variable,
  unique(tabla_diagnostico$Variable)
)

if (length(variables_validas) == 0) {
  stop("No hay variables válidas para diagnóstico residual.")
}


## 3. Funciones auxiliares -------------------------------------

crear_pesos_knn <- function(coords, k = 4) {

  n <- nrow(coords)

  if (n < 3) {
    return(NULL)
  }

  k_real <- min(k, n - 1)

  vecinos <- knearneigh(
    coords,
    k = k_real
  ) |>
    knn2nb()

  nb2listw(
    vecinos,
    style = "W",
    zero.policy = TRUE
  )
}


diagnostico_moran <- function(df, variable, nsim = 999, k = 4) {

  df <- df |>
    filter(
      Variable == variable,
      is.finite(Residuo),
      is.finite(X),
      is.finite(Y)
    )

  if (nrow(df) < 5) {
    return(
      list(
        resumen = tibble(
          Variable = variable,
          n = nrow(df),
          k_vecinos = NA_integer_,
          Moran_I = NA_real_,
          Esperanza_I = NA_real_,
          Varianza_I = NA_real_,
          Z_Moran = NA_real_,
          p_valor_informativo = NA_real_,
          Sim_Q025 = NA_real_,
          Sim_Mediana = NA_real_,
          Sim_Q975 = NA_real_,
          Diagnostico = "Número insuficiente de observaciones"
        ),
        simulaciones = tibble()
      )
    )
  }

  coords <- as.matrix(
    df |>
      select(X, Y)
  )

  listw <- crear_pesos_knn(
    coords,
    k = k
  )

  prueba <- moran.test(
    df$Residuo,
    listw,
    zero.policy = TRUE,
    randomisation = TRUE
  )

  simulacion <- moran.mc(
    df$Residuo,
    listw,
    nsim = nsim,
    zero.policy = TRUE,
    return_boot = TRUE
  )

  valores_sim <- as.numeric(simulacion$t)

  moran_i <- as.numeric(
    prueba$estimate[["Moran I statistic"]]
  )

  esperanza <- as.numeric(
    prueba$estimate[["Expectation"]]
  )

  varianza <- as.numeric(
    prueba$estimate[["Variance"]]
  )

  z_moran <- ifelse(
    is.finite(varianza) && varianza > 0,
    (moran_i - esperanza) / sqrt(varianza),
    NA_real_
  )

  q_sim <- quantile(
    valores_sim,
    probs = c(0.025, 0.50, 0.975),
    na.rm = TRUE,
    names = FALSE
  )

  diagnostico <- case_when(
    !is.finite(moran_i) ~
      "No interpretable",

    moran_i > q_sim[3] ~
      "Autocorrelación positiva mayor que el intervalo simulado",

    moran_i < q_sim[1] ~
      "Autocorrelación negativa menor que el intervalo simulado",

    TRUE ~
      "Sin evidencia gráfica fuerte de estructura espacial residual"
  )

  list(
    resumen = tibble(
      Variable = variable,
      n = nrow(df),
      k_vecinos = min(k, nrow(df) - 1),
      Moran_I = moran_i,
      Esperanza_I = esperanza,
      Varianza_I = varianza,
      Z_Moran = z_moran,
      p_valor_informativo = prueba$p.value,
      Sim_Q025 = q_sim[1],
      Sim_Mediana = q_sim[2],
      Sim_Q975 = q_sim[3],
      Diagnostico = diagnostico
    ),
    simulaciones = tibble(
      Variable = variable,
      Moran_simulado = valores_sim
    )
  )
}


crear_variograma_residual <- function(df, variable) {

  df <- df |>
    filter(
      Variable == variable,
      is.finite(Residuo),
      is.finite(X),
      is.finite(Y)
    )

  if (nrow(df) < 5) {
    return(NULL)
  }

  sf_res <- st_as_sf(
    df,
    coords = c("X", "Y"),
    crs = 32618,
    remove = FALSE
  )

  distancia_maxima <- max(
    as.vector(
      dist(
        st_coordinates(sf_res)
      )
    ),
    na.rm = TRUE
  )

  cutoff <- 0.50 * distancia_maxima
  numero_lags <- 7
  width <- cutoff / numero_lags

  vexp <- variogram(
    Residuo ~ 1,
    sf_res,
    cutoff = cutoff,
    width = width
  )

  as_tibble(vexp) |>
    transmute(
      Variable = variable,
      Numero_pares = np,
      Distancia_media_m = dist,
      Semivarianza = gamma
    )
}


## 4. Moran I para todos los modelos ---------------------------

resultados_moran <- map(
  variables_validas,
  function(variable) {
    diagnostico_moran(
      df = tabla_diagnostico,
      variable = variable,
      nsim = 999,
      k = 4
    )
  }
)

names(resultados_moran) <- variables_validas

tabla_moran <- bind_rows(
  map(resultados_moran, "resumen")
) |>
  left_join(
    catalogo_respuestas,
    by = "Variable"
  )

simulaciones_moran <- bind_rows(
  map(resultados_moran, "simulaciones")
)


## 5. Semivariogramas de residuos ------------------------------

variogramas_residuales <- map(
  variables_validas,
  function(variable) {
    crear_variograma_residual(
      df = tabla_diagnostico,
      variable = variable
    )
  }
)

names(variogramas_residuales) <- variables_validas

variogramas_residuales <- variogramas_residuales[
  !map_lgl(variogramas_residuales, is.null)
]

tabla_variogramas <- bind_rows(
  variogramas_residuales
)


## 6. Observaciones influyentes --------------------------------

tabla_influencia <- tabla_diagnostico |>
  group_by(Variable) |>
  mutate(
    n_modelo = n(),
    umbral_cook = 4 / n_modelo,
    umbral_leverage = 2 * mean(Leverage, na.rm = TRUE),
    Flag_Cook = Cook > umbral_cook,
    Flag_Leverage = Leverage > umbral_leverage,
    Flag_Residuo = abs(Residuo_estudentizado) > 2,
    Influyente = Flag_Cook | Flag_Leverage | Flag_Residuo
  ) |>
  ungroup()

resumen_influencia <- tabla_influencia |>
  group_by(Variable) |>
  summarise(
    n = n(),
    n_Cook = sum(Flag_Cook, na.rm = TRUE),
    n_Leverage = sum(Flag_Leverage, na.rm = TRUE),
    n_Residuo_extremo = sum(Flag_Residuo, na.rm = TRUE),
    n_Influyentes = sum(Influyente, na.rm = TRUE),
    Cook_max = max(Cook, na.rm = TRUE),
    Residuo_estudentizado_max = max(
      abs(Residuo_estudentizado),
      na.rm = TRUE
    ),
    .groups = "drop"
  )


## 7. Síntesis de modelos ---------------------------------------

tabla_sintesis <- modelos_seleccionados |>
  left_join(
    tabla_moran |>
      select(
        Variable,
        Moran_I,
        Sim_Q025,
        Sim_Q975,
        Diagnostico
      ),
    by = "Variable"
  ) |>
  left_join(
    resumen_influencia |>
      select(
        Variable,
        n_Influyentes,
        Cook_max,
        Residuo_estudentizado_max
      ),
    by = "Variable"
  ) |>
  mutate(
    Estabilidad_predictiva = case_when(
      Spearman_LOOCV >= 0.70 ~ "Alta",
      Spearman_LOOCV >= 0.40 ~ "Moderada",
      Spearman_LOOCV >= 0.20 ~ "Limitada",
      TRUE ~ "Muy limitada"
    ),
    Estado_residual = case_when(
      str_detect(
        Diagnostico,
        "Sin evidencia"
      ) ~ "Adecuado espacialmente",
      is.na(Diagnostico) ~ "No evaluado",
      TRUE ~ "Requiere discusión"
    )
  )


## 8. Figuras: mapas de residuos -------------------------------

crear_mapa_residuos <- function(variable) {

  df <- tabla_diagnostico |>
    filter(
      Variable == variable
    )

  ggplot(
    df,
    aes(
      x = X,
      y = Y
    )
  ) +
    geom_point(
      aes(
        fill = Residuo,
        size = abs(Residuo)
      ),
      shape = 21
    ) +
    coord_equal() +
    labs(
      title = variable,
      x = "X (m)",
      y = "Y (m)",
      fill = "Residuo",
      size = "|Residuo|"
    ) +
    theme_minimal()
}

mapas_residuos <- map(
  variables_validas,
  crear_mapa_residuos
)

names(mapas_residuos) <- variables_validas

panel_mapas <- wrap_plots(
  mapas_residuos,
  ncol = 3
) +
  plot_annotation(
    title = "Distribución espacial de residuos de los modelos de tratamientos"
  )

ggsave(
  filename = file.path(
    ruta_figuras,
    "Panel_mapas_residuos_tratamientos.png"
  ),
  plot = panel_mapas,
  width = 16,
  height = 18,
  dpi = 400
)


## 9. Figuras: semivariogramas residuales ----------------------

if (nrow(tabla_variogramas) > 0) {

  g_variogramas <- ggplot(
    tabla_variogramas,
    aes(
      x = Distancia_media_m,
      y = Semivarianza,
      size = Numero_pares
    )
  ) +
    geom_point() +
    geom_line(
      aes(group = Variable),
      linewidth = 0.35
    ) +
    facet_wrap(
      ~ Variable,
      scales = "free",
      ncol = 3
    ) +
    labs(
      title = "Semivariogramas experimentales de residuos",
      x = "Distancia (m)",
      y = "Semivarianza",
      size = "Pares"
    ) +
    theme_minimal()

  ggsave(
    filename = file.path(
      ruta_figuras,
      "Semivariogramas_residuales_tratamientos.png"
    ),
    plot = g_variogramas,
    width = 15,
    height = 17,
    dpi = 400
  )
}


## 10. Figuras: Moran observado y simulaciones -----------------

if (nrow(simulaciones_moran) > 0) {

  g_moran <- ggplot(
    simulaciones_moran,
    aes(x = Moran_simulado)
  ) +
    geom_histogram(
      bins = 30,
      color = "white"
    ) +
    geom_vline(
      data = tabla_moran,
      aes(xintercept = Moran_I),
      linetype = "dashed",
      linewidth = 0.7
    ) +
    facet_wrap(
      ~ Variable,
      scales = "free_y",
      ncol = 3
    ) +
    labs(
      title = "Moran I observado frente a distribución simulada",
      subtitle = "Línea discontinua: Moran I observado",
      x = "Moran I",
      y = "Frecuencia simulada"
    ) +
    theme_minimal()

  ggsave(
    filename = file.path(
      ruta_figuras,
      "Moran_residuos_distribuciones_simuladas.png"
    ),
    plot = g_moran,
    width = 15,
    height = 17,
    dpi = 400
  )
}


## 11. Figura: observado, ajustado y residuos ------------------

g_obs_ajustado <- ggplot(
  tabla_diagnostico,
  aes(
    x = Observado,
    y = Ajustado
  )
) +
  geom_point() +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed"
  ) +
  facet_wrap(
    ~ Variable,
    scales = "free",
    ncol = 3
  ) +
  labs(
    title = "Valores observados frente a ajustados",
    x = "Observado",
    y = "Ajustado"
  ) +
  theme_minimal()

ggsave(
  filename = file.path(
    ruta_figuras,
    "Panel_observado_ajustado_tratamientos.png"
  ),
  plot = g_obs_ajustado,
  width = 15,
  height = 17,
  dpi = 400
)


## 12. Exportar datos y tablas ---------------------------------

write_csv(
  tabla_influencia,
  file.path(
    ruta_datos,
    "Diagnostico_influencia_espacial_residuos.csv"
  )
)

write_csv(
  tabla_variogramas,
  file.path(
    ruta_datos,
    "Semivariogramas_residuales.csv"
  )
)

write_csv(
  tabla_moran,
  file.path(
    ruta_datos,
    "Moran_residuos.csv"
  )
)

write_xlsx(
  list(
    Sintesis_modelos = tabla_sintesis,
    Moran_residuos = tabla_moran,
    Semivariogramas_residuos = tabla_variogramas,
    Resumen_influencia = resumen_influencia,
    Observaciones_influyentes = tabla_influencia |>
      filter(Influyente),
    Diagnosticos_completos = tabla_influencia,
    Modelos_seleccionados = modelos_seleccionados,
    Tamanos_efecto = tabla_tamanos_efecto
  ),
  path = file.path(
    ruta_tablas,
    "Resultados_diagnostico_residual_tratamientos.xlsx"
  )
)


## 13. Guardar objeto final ------------------------------------

saveRDS(
  list(
    tabla_sintesis = tabla_sintesis,
    tabla_moran = tabla_moran,
    tabla_variogramas = tabla_variogramas,
    tabla_influencia = tabla_influencia,
    resumen_influencia = resumen_influencia,
    simulaciones_moran = simulaciones_moran
  ),
  file = file.path(
    ruta_objetos,
    "objetos_diagnostico_residual_tratamientos.rds"
  )
)


## 14. Registro de sesión --------------------------------------

capture.output(
  sessionInfo(),
  file = file.path(
    ruta_salida,
    "sessionInfo_04_04.txt"
  )
)


## 15. Resumen en consola --------------------------------------

cat(
  "\n============================================\n"
)

cat(
  "SCRIPT 04-04 FINALIZADO\n\n"
)

cat(
  "Variables evaluadas:\n",
  paste(variables_validas, collapse = ", "),
  "\n\n"
)

cat(
  "Resumen de autocorrelación residual:\n"
)

print(
  tabla_moran |>
    select(
      Variable,
      Moran_I,
      Sim_Q025,
      Sim_Q975,
      Diagnostico
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
