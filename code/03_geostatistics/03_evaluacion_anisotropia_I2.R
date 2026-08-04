
## ============================================================
## SCRIPT 03. EVALUACIÓN DE ANISOTROPÍA – SUELO I2
## Artículo cebada: revisión mayor
##
## Entrada:
##   AJUSTES_ARTÍCULO/02_SEMIVARIOGRAMAS_I2/objetos_R/
##   objetos_semivariogramas_I2.rds
##
## Salidas:
##   AJUSTES_ARTÍCULO/03_ANISOTROPIA_I2/
##
## Objetivos:
##   1. Construir semivariogramas direccionales a 0°, 45°, 90° y 135°.
##   2. Documentar el número de pares por dirección e intervalo.
##   3. Ajustar el modelo teórico seleccionado en cada dirección.
##   4. Probar modelos anisotrópicos geométricos mediante una rejilla
##      limitada de direcciones y relaciones de anisotropía.
##   5. Comparar el modelo isotrópico seleccionado frente a candidatos
##      anisotrópicos mediante LOOCV.
##   6. Recomendar anisotropía solo cuando exista evidencia direccional
##      suficiente y una mejora predictiva relevante.
##
## IMPORTANTE:
##   - El Script 03 NO genera todavía los mapas finales.
##   - El Kriging y los mapas de varianza se generan en el Script 04.
##   - pH se conserva únicamente como diagnóstico exploratorio.
##   - N, P y K son índices instrumentales nominales y redundantes con CE.
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
library(patchwork)


## 1. Rutas ----------------------------------------------------

ruta_base <- paste0(
  "D:/UDCA_2017_CLASE/46_ARTICULOS_2025/",
  "ARTICULO_CEBADA_DOCTORADO/AJUSTES_ARTÍCULO"
)

ruta_objeto <- file.path(
  ruta_base,
  "02_SEMIVARIOGRAMAS_I2",
  "objetos_R",
  "objetos_semivariogramas_I2.rds"
)

ruta_salida <- file.path(
  ruta_base,
  "03_ANISOTROPIA_I2"
)

ruta_tablas <- file.path(
  ruta_salida,
  "tablas"
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
    "No existe el objeto del Script 02:\n",
    ruta_objeto,
    "\nEjecuta primero el Script 02."
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
  ruta_objetos,
  recursive = TRUE,
  showWarnings = FALSE
)


## 2. Cargar resultados del Script 02 --------------------------

obj <- readRDS(
  ruta_objeto
)

suelo_sf <- obj$suelo_sf
lote <- obj$lote
tabla_modelos_seleccionados <- obj$tabla_modelos_seleccionados
modelos_ajustados <- obj$modelos_ajustados
cutoff_global <- obj$cutoff
width_global <- obj$width
variables_modelar <- obj$variables

if (nrow(tabla_modelos_seleccionados) == 0) {
  stop(
    "El Script 02 no dejó modelos seleccionados."
  )
}


## 3. Configuración direccional --------------------------------

angulos <- c(
  0,
  45,
  90,
  135
)

tolerancia_horizontal <- 22.5

# Para evitar conclusiones apoyadas en muy pocos pares.
min_pares_intervalo <- 5
min_pares_direccion_total <- 20
min_intervalos_utiles <- 3

# Rejilla limitada para anisotropía geométrica.
# gstat usa anis = c(ángulo, relación eje menor/eje mayor).
relaciones_anisotropia <- c(
  0.35,
  0.50,
  0.65,
  0.80
)

# Criterios para recomendar anisotropía:
mejora_minima_rmse_porcentaje <- 5
mejora_minima_msne <- 0.05


## 4. Funciones auxiliares -------------------------------------

crear_variograma_direccional <- function(
  variable,
  datos_sf,
  angulo,
  cutoff,
  width,
  tol_hor
) {

  formula_variable <- as.formula(
    paste(
      variable,
      "~ 1"
    )
  )

  variogram(
    formula_variable,
    datos_sf,
    alpha = angulo,
    tol.hor = tol_hor,
    cutoff = cutoff,
    width = width,
    cressie = FALSE
  )
}


extraer_componentes <- function(modelo) {

  nugget <- modelo |>
    filter(model == "Nug") |>
    pull(psill)

  if (length(nugget) == 0) {
    nugget <- 0
  }

  estructura <- modelo |>
    filter(model != "Nug")

  if (nrow(estructura) == 0) {
    return(
      tibble(
        Pepita = nugget[1],
        Psill = NA_real_,
        Silla = NA_real_,
        Rango = NA_real_
      )
    )
  }

  psill <- estructura$psill[1]
  rango <- estructura$range[1]

  tibble(
    Pepita = nugget[1],
    Psill = psill,
    Silla = nugget[1] + psill,
    Rango = rango
  )
}


calcular_cv <- function(
  variable,
  datos_sf,
  modelo
) {

  formula_variable <- as.formula(
    paste(
      variable,
      "~ 1"
    )
  )

  cv <- try(
    krige.cv(
      formula_variable,
      datos_sf,
      model = modelo,
      nfold = nrow(datos_sf),
      verbose = FALSE
    ),
    silent = TRUE
  )

  if (inherits(cv, "try-error")) {
    return(
      list(
        metricas = tibble(
          ME = NA_real_,
          MAE = NA_real_,
          RMSE = NA_real_,
          RMSSE = NA_real_,
          MSNE = NA_real_,
          Correlacion_obs_pred = NA_real_
        ),
        datos = NULL
      )
    )
  }

  residuo <- cv$observed - cv$var1.pred

  var_valida <- is.finite(cv$var1.var) &
    cv$var1.var > 0

  residuo_est <- rep(
    NA_real_,
    length(residuo)
  )

  residuo_est[var_valida] <-
    residuo[var_valida] /
    sqrt(cv$var1.var[var_valida])

  metricas <- tibble(
    ME = mean(
      residuo,
      na.rm = TRUE
    ),
    MAE = mean(
      abs(residuo),
      na.rm = TRUE
    ),
    RMSE = sqrt(
      mean(
        residuo^2,
        na.rm = TRUE
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
    ),
    Correlacion_obs_pred = suppressWarnings(
      cor(
        cv$observed,
        cv$var1.pred,
        use = "complete.obs"
      )
    )
  )

  list(
    metricas = metricas,
    datos = cv
  )
}


ajustar_modelo_direccional <- function(
  vexp_dir,
  modelo_iso
) {

  if (
    is.null(vexp_dir) ||
      nrow(vexp_dir) < 3
  ) {
    return(NULL)
  }

  modelo_inicial <- modelo_iso

  ajuste <- try(
    fit.variogram(
      vexp_dir,
      modelo_inicial,
      fit.method = 6
    ),
    silent = TRUE
  )

  if (
    inherits(ajuste, "try-error") ||
      isTRUE(attr(ajuste, "singular"))
  ) {
    return(NULL)
  }

  ajuste
}


crear_modelo_anisotropico <- function(
  modelo_iso,
  angulo,
  relacion
) {

  fila_nugget <- modelo_iso |>
    filter(model == "Nug")

  fila_estructura <- modelo_iso |>
    filter(model != "Nug")

  if (nrow(fila_estructura) == 0) {
    return(NULL)
  }

  nugget <- ifelse(
    nrow(fila_nugget) > 0,
    fila_nugget$psill[1],
    0
  )

  vgm(
    psill = fila_estructura$psill[1],
    model = as.character(
      fila_estructura$model[1]
    ),
    range = fila_estructura$range[1],
    nugget = nugget,
    anis = c(
      angulo,
      relacion
    )
  )
}


clasificar_evidencia_direccional <- function(
  tabla_direcciones
) {

  direcciones_validas <- tabla_direcciones |>
    filter(
      Ajuste_direccional_valido,
      Pares_totales >= min_pares_direccion_total,
      Intervalos_utiles >= min_intervalos_utiles,
      is.finite(Rango_direccional)
    )

  if (nrow(direcciones_validas) < 2) {
    return(
      tibble(
        Evidencia_direccional = "Insuficiente",
        Direccion_mayor_continuidad = NA_real_,
        Rango_mayor = NA_real_,
        Rango_menor = NA_real_,
        Relacion_rangos = NA_real_
      )
    )
  }

  rango_mayor <- max(
    direcciones_validas$Rango_direccional,
    na.rm = TRUE
  )

  rango_menor <- min(
    direcciones_validas$Rango_direccional,
    na.rm = TRUE
  )

  direccion_mayor <- direcciones_validas |>
    filter(
      Rango_direccional == rango_mayor
    ) |>
    slice(1) |>
    pull(Angulo)

  relacion <- rango_menor / rango_mayor

  evidencia <- case_when(
    !is.finite(relacion) ~ "Insuficiente",
    relacion <= 0.50 ~ "Fuerte",
    relacion <= 0.75 ~ "Moderada",
    TRUE ~ "Débil"
  )

  tibble(
    Evidencia_direccional = evidencia,
    Direccion_mayor_continuidad = direccion_mayor,
    Rango_mayor = rango_mayor,
    Rango_menor = rango_menor,
    Relacion_rangos = relacion
  )
}


## 5. Variogramas direccionales y ajustes ----------------------

variogramas_direccionales <- list()
modelos_direccionales <- list()
tablas_direccionales <- list()

for (variable in variables_modelar) {

  message(
    "Variogramas direccionales: ",
    variable
  )

  fila_sel <- tabla_modelos_seleccionados |>
    filter(
      Variable == variable
    )

  if (nrow(fila_sel) == 0) {
    next
  }

  codigo_modelo <- fila_sel$Modelo_codigo[1]
  modelo_iso <- modelos_ajustados[[variable]][[codigo_modelo]]

  lista_vdir <- list()
  lista_mdir <- list()
  lista_resumen <- list()

  for (angulo in angulos) {

    vdir <- crear_variograma_direccional(
      variable = variable,
      datos_sf = suelo_sf,
      angulo = angulo,
      cutoff = cutoff_global,
      width = width_global,
      tol_hor = tolerancia_horizontal
    )

    ajuste_dir <- ajustar_modelo_direccional(
      vdir,
      modelo_iso
    )

    pares_totales <- sum(
      vdir$np,
      na.rm = TRUE
    )

    intervalos_utiles <- sum(
      vdir$np >= min_pares_intervalo,
      na.rm = TRUE
    )

    if (is.null(ajuste_dir)) {

      resumen_dir <- tibble(
        Variable = variable,
        Angulo = angulo,
        Pares_totales = pares_totales,
        Intervalos = nrow(vdir),
        Intervalos_utiles = intervalos_utiles,
        Ajuste_direccional_valido = FALSE,
        Pepita_direccional = NA_real_,
        Silla_direccional = NA_real_,
        Rango_direccional = NA_real_,
        SSErr_direccional = NA_real_
      )

    } else {

      componentes <- extraer_componentes(
        ajuste_dir
      )

      ss_err <- attr(
        ajuste_dir,
        "SSErr"
      )

      resumen_dir <- tibble(
        Variable = variable,
        Angulo = angulo,
        Pares_totales = pares_totales,
        Intervalos = nrow(vdir),
        Intervalos_utiles = intervalos_utiles,
        Ajuste_direccional_valido = TRUE,
        Pepita_direccional = componentes$Pepita,
        Silla_direccional = componentes$Silla,
        Rango_direccional = componentes$Rango,
        SSErr_direccional = ifelse(
          length(ss_err) == 0,
          NA_real_,
          as.numeric(ss_err)
        )
      )
    }

    lista_vdir[[as.character(angulo)]] <- vdir
    lista_mdir[[as.character(angulo)]] <- ajuste_dir
    lista_resumen[[as.character(angulo)]] <- resumen_dir
  }

  variogramas_direccionales[[variable]] <- lista_vdir
  modelos_direccionales[[variable]] <- lista_mdir
  tablas_direccionales[[variable]] <- bind_rows(
    lista_resumen
  )
}


tabla_direccional_completa <- bind_rows(
  tablas_direccionales
)


## 6. Evidencia direccional por variable -----------------------

tabla_evidencia_direccional <- map_dfr(
  variables_modelar,
  function(variable) {

    tabla_var <- tabla_direccional_completa |>
      filter(
        Variable == variable
      )

    if (nrow(tabla_var) == 0) {
      return(
        tibble(
          Variable = variable,
          Evidencia_direccional = "Sin resultados",
          Direccion_mayor_continuidad = NA_real_,
          Rango_mayor = NA_real_,
          Rango_menor = NA_real_,
          Relacion_rangos = NA_real_
        )
      )
    }

    bind_cols(
      tibble(
        Variable = variable
      ),
      clasificar_evidencia_direccional(
        tabla_var
      )
    )
  }
)

print(
  tabla_evidencia_direccional
)


## 7. Rejilla de candidatos anisotrópicos ----------------------

resultados_anisotropia <- list()
modelos_anisotropicos <- list()
cv_anisotropicos <- list()

for (variable in variables_modelar) {

  message(
    "Comparación isotrópica/anisotrópica: ",
    variable
  )

  fila_iso <- tabla_modelos_seleccionados |>
    filter(
      Variable == variable
    )

  if (nrow(fila_iso) == 0) {
    next
  }

  codigo_modelo <- fila_iso$Modelo_codigo[1]
  modelo_iso <- modelos_ajustados[[variable]][[codigo_modelo]]

  cv_iso <- calcular_cv(
    variable,
    suelo_sf,
    modelo_iso
  )

  tabla_iso <- bind_cols(
    tibble(
      Variable = variable,
      Tipo_modelo = "Isotrópico",
      Angulo = NA_real_,
      Relacion_anisotropia = 1,
      Modelo_codigo = codigo_modelo
    ),
    cv_iso$metricas
  )

  candidatos_var <- list()
  modelos_var <- list()
  cv_var <- list()

  indice <- 1

  for (angulo in angulos) {

    for (relacion in relaciones_anisotropia) {

      modelo_anis <- crear_modelo_anisotropico(
        modelo_iso,
        angulo,
        relacion
      )

      if (is.null(modelo_anis)) {
        next
      }

      cv_anis <- calcular_cv(
        variable,
        suelo_sf,
        modelo_anis
      )

      fila <- bind_cols(
        tibble(
          Variable = variable,
          Tipo_modelo = "Anisotrópico",
          Angulo = angulo,
          Relacion_anisotropia = relacion,
          Modelo_codigo = codigo_modelo
        ),
        cv_anis$metricas
      )

      candidatos_var[[indice]] <- fila
      modelos_var[[indice]] <- modelo_anis
      cv_var[[indice]] <- cv_anis$datos

      indice <- indice + 1
    }
  }

  tabla_anis_var <- bind_rows(
    candidatos_var
  )

  if (nrow(tabla_anis_var) > 0) {

    mejor_anis <- tabla_anis_var |>
      filter(
        is.finite(RMSE),
        is.finite(MSNE)
      ) |>
      mutate(
        Score =
          rank(RMSE, ties.method = "min") +
          rank(
            abs(MSNE - 1),
            ties.method = "min"
          ) +
          0.5 * rank(
            MAE,
            ties.method = "min"
          )
      ) |>
      arrange(
        Score,
        RMSE,
        abs(MSNE - 1)
      ) |>
      slice(1)

  } else {

    mejor_anis <- tibble()
  }

  evidencia_var <- tabla_evidencia_direccional |>
    filter(
      Variable == variable
    )

  if (nrow(mejor_anis) == 1) {

    mejora_rmse <- 100 * (
      tabla_iso$RMSE -
        mejor_anis$RMSE
    ) / tabla_iso$RMSE

    mejora_msne <- abs(
      tabla_iso$MSNE - 1
    ) - abs(
      mejor_anis$MSNE - 1
    )

    evidencia_suficiente <-
      evidencia_var$Evidencia_direccional %in%
      c(
        "Moderada",
        "Fuerte"
      )

    recomendar <- (
      evidencia_suficiente &&
        is.finite(mejora_rmse) &&
        mejora_rmse >= mejora_minima_rmse_porcentaje &&
        is.finite(mejora_msne) &&
        mejora_msne >= mejora_minima_msne
    )

    comparacion <- bind_cols(
      tabla_iso,
      tibble(
        Mejor_angulo_anisotropico =
          mejor_anis$Angulo,
        Mejor_relacion_anisotropia =
          mejor_anis$Relacion_anisotropia,
        RMSE_mejor_anisotropico =
          mejor_anis$RMSE,
        MAE_mejor_anisotropico =
          mejor_anis$MAE,
        MSNE_mejor_anisotropico =
          mejor_anis$MSNE,
        Mejora_RMSE_porcentaje =
          mejora_rmse,
        Mejora_MSNE =
          mejora_msne,
        Evidencia_direccional =
          evidencia_var$Evidencia_direccional,
        Recomendar_anisotropia =
          recomendar,
        Motivo_decision = case_when(
          !evidencia_suficiente ~
            "Evidencia direccional insuficiente o débil",
          mejora_rmse <
            mejora_minima_rmse_porcentaje ~
            "Mejora de RMSE inferior al umbral",
          mejora_msne <
            mejora_minima_msne ~
            "MSNE no mejora suficientemente",
          TRUE ~
            "Evidencia direccional y mejora predictiva"
        )
      )
    )

    indice_mejor <- which(
      tabla_anis_var$Angulo ==
        mejor_anis$Angulo &
        tabla_anis_var$Relacion_anisotropia ==
        mejor_anis$Relacion_anisotropia
    )[1]

    modelos_anisotropicos[[variable]] <-
      modelos_var[[indice_mejor]]

    cv_anisotropicos[[variable]] <-
      cv_var[[indice_mejor]]

  } else {

    comparacion <- bind_cols(
      tabla_iso,
      tibble(
        Mejor_angulo_anisotropico =
          NA_real_,
        Mejor_relacion_anisotropia =
          NA_real_,
        RMSE_mejor_anisotropico =
          NA_real_,
        MAE_mejor_anisotropico =
          NA_real_,
        MSNE_mejor_anisotropico =
          NA_real_,
        Mejora_RMSE_porcentaje =
          NA_real_,
        Mejora_MSNE =
          NA_real_,
        Evidencia_direccional =
          evidencia_var$Evidencia_direccional,
        Recomendar_anisotropia =
          FALSE,
        Motivo_decision =
          "No se obtuvo candidato anisotrópico válido"
      )
    )

    modelos_anisotropicos[[variable]] <- NULL
    cv_anisotropicos[[variable]] <- NULL
  }

  resultados_anisotropia[[variable]] <- list(
    isotropico = tabla_iso,
    candidatos = tabla_anis_var,
    comparacion = comparacion
  )
}


tabla_candidatos_anisotropicos <- bind_rows(
  map(
    resultados_anisotropia,
    "candidatos"
  )
)

tabla_comparacion_iso_anis <- bind_rows(
  map(
    resultados_anisotropia,
    "comparacion"
  )
)

print(
  tabla_comparacion_iso_anis |>
    select(
      Variable,
      RMSE,
      MSNE,
      Mejor_angulo_anisotropico,
      Mejor_relacion_anisotropia,
      RMSE_mejor_anisotropico,
      MSNE_mejor_anisotropico,
      Mejora_RMSE_porcentaje,
      Evidencia_direccional,
      Recomendar_anisotropia,
      Motivo_decision
    )
)


## 8. Modelo final recomendado para Script 04 -----------------

modelos_finales <- list()
tabla_modelos_finales <- list()

for (variable in variables_modelar) {

  fila_sel <- tabla_modelos_seleccionados |>
    filter(
      Variable == variable
    )

  fila_comp <- tabla_comparacion_iso_anis |>
    filter(
      Variable == variable
    )

  if (
    nrow(fila_sel) == 0 ||
      nrow(fila_comp) == 0
  ) {
    next
  }

  codigo_modelo <- fila_sel$Modelo_codigo[1]
  modelo_iso <- modelos_ajustados[[variable]][[codigo_modelo]]

  if (
    isTRUE(
      fila_comp$Recomendar_anisotropia[1]
    ) &&
      !is.null(
        modelos_anisotropicos[[variable]]
      )
  ) {

    modelo_final <- modelos_anisotropicos[[variable]]
    tipo_final <- "Anisotrópico"
    angulo_final <- fila_comp$Mejor_angulo_anisotropico[1]
    relacion_final <- fila_comp$Mejor_relacion_anisotropia[1]

  } else {

    modelo_final <- modelo_iso
    tipo_final <- "Isotrópico"
    angulo_final <- NA_real_
    relacion_final <- 1
  }

  modelos_finales[[variable]] <- modelo_final

  tabla_modelos_finales[[variable]] <- tibble(
    Variable = variable,
    Modelo_codigo = codigo_modelo,
    Tipo_modelo_final = tipo_final,
    Angulo_final = angulo_final,
    Relacion_anisotropia_final = relacion_final,
    Recomendar_anisotropia =
      fila_comp$Recomendar_anisotropia[1],
    Motivo_decision =
      fila_comp$Motivo_decision[1]
  )
}

tabla_modelos_finales <- bind_rows(
  tabla_modelos_finales
)

print(
  tabla_modelos_finales
)


## 9. Figuras direccionales ------------------------------------

crear_panel_direccional <- function(
  variable
) {

  lista_v <- variogramas_direccionales[[variable]]
  lista_m <- modelos_direccionales[[variable]]

  if (is.null(lista_v)) {
    return(NULL)
  }

  graficos <- map(
    as.character(angulos),
    function(angulo_txt) {

      vdir <- lista_v[[angulo_txt]]
      mdir <- lista_m[[angulo_txt]]

      g <- ggplot(
        as_tibble(vdir),
        aes(
          x = dist,
          y = gamma
        )
      ) +
        geom_point(
          aes(
            size = np
          )
        ) +
        labs(
          title = paste0(
            angulo_txt,
            "°"
          ),
          x = "Distancia (m)",
          y = "Semivarianza",
          size = "Pares"
        ) +
        theme_minimal()

      if (!is.null(mdir)) {

        curva <- variogramLine(
          mdir,
          maxdist = max(
            vdir$dist,
            na.rm = TRUE
          ),
          n = 200
        )

        g <- g +
          geom_line(
            data = curva,
            aes(
              x = dist,
              y = gamma
            ),
            inherit.aes = FALSE
          )
      }

      g
    }
  )

  wrap_plots(
    graficos,
    ncol = 2
  ) +
    plot_annotation(
      title = paste0(
        "Semivariogramas direccionales de ",
        variable
      )
    )
}


paneles_direccionales <- list()

for (variable in variables_modelar) {

  panel <- crear_panel_direccional(
    variable
  )

  if (is.null(panel)) {
    next
  }

  paneles_direccionales[[variable]] <- panel

  ggsave(
    filename = file.path(
      ruta_figuras,
      paste0(
        "Variogramas_direccionales_",
        variable,
        ".png"
      )
    ),
    plot = panel,
    width = 11,
    height = 8,
    dpi = 400
  )
}


## 10. Figura resumen de rangos direccionales ------------------

g_rangos <- tabla_direccional_completa |>
  filter(
    Ajuste_direccional_valido,
    is.finite(Rango_direccional)
  ) |>
  ggplot(
    aes(
      x = factor(Angulo),
      y = Rango_direccional,
      group = Variable
    )
  ) +
  geom_point(
    size = 2.5
  ) +
  geom_line() +
  facet_wrap(
    ~ Variable,
    scales = "free_y",
    ncol = 3
  ) +
  labs(
    title =
      "Rangos ajustados por dirección en I2",
    subtitle =
      "Interpretar únicamente direcciones con suficiente número de pares",
    x =
      "Dirección (°)",
    y =
      "Rango del modelo direccional (m)"
  ) +
  theme_minimal()

ggsave(
  filename = file.path(
    ruta_figuras,
    "Resumen_rangos_direccionales_I2.png"
  ),
  plot = g_rangos,
  width = 13,
  height = 10,
  dpi = 400
)


## 11. Exportación de tablas -----------------------------------

write_xlsx(
  list(
    Configuracion =
      tibble(
        Parametro = c(
          "Ángulos",
          "Tolerancia horizontal (°)",
          "Mínimo pares por intervalo",
          "Mínimo pares totales por dirección",
          "Mínimo intervalos útiles",
          "Mejora mínima RMSE (%)",
          "Mejora mínima MSNE"
        ),
        Valor = c(
          paste(
            angulos,
            collapse = ", "
          ),
          tolerancia_horizontal,
          min_pares_intervalo,
          min_pares_direccion_total,
          min_intervalos_utiles,
          mejora_minima_rmse_porcentaje,
          mejora_minima_msne
        )
      ),

    Direccionales =
      tabla_direccional_completa,

    Evidencia_direccional =
      tabla_evidencia_direccional,

    Candidatos_anisotropicos =
      tabla_candidatos_anisotropicos,

    Comparacion_iso_anis =
      tabla_comparacion_iso_anis,

    Modelos_finales =
      tabla_modelos_finales
  ),
  path = file.path(
    ruta_tablas,
    "Resultados_anisotropia_I2.xlsx"
  )
)


## 12. Guardar objetos para Kriging – Script 04 ----------------

saveRDS(
  list(
    suelo_sf =
      suelo_sf,

    lote =
      lote,

    variables =
      variables_modelar,

    tabla_modelos_seleccionados =
      tabla_modelos_seleccionados,

    tabla_direccional_completa =
      tabla_direccional_completa,

    tabla_evidencia_direccional =
      tabla_evidencia_direccional,

    tabla_comparacion_iso_anis =
      tabla_comparacion_iso_anis,

    tabla_modelos_finales =
      tabla_modelos_finales,

    modelos_finales =
      modelos_finales,

    modelos_isotropicos =
      modelos_ajustados,

    modelos_anisotropicos =
      modelos_anisotropicos,

    cutoff =
      cutoff_global,

    width =
      width_global
  ),
  file = file.path(
    ruta_objetos,
    "objetos_anisotropia_I2.rds"
  )
)


## 13. Registro de sesión --------------------------------------

capture.output(
  sessionInfo(),
  file = file.path(
    ruta_salida,
    "sessionInfo_Script03.txt"
  )
)


## 14. Resumen en consola --------------------------------------

cat(
  "\n============================================\n"
)

cat(
  "SCRIPT 03 FINALIZADO\n\n"
)

cat(
  "Decisión isotropía/anisotropía:\n"
)

print(
  tabla_modelos_finales
)

cat(
  "\nResultados guardados en:\n",
  ruta_salida,
  "\n"
)

cat(
  "\nObjeto para Script 04 (Kriging):\n",
  file.path(
    ruta_objetos,
    "objetos_anisotropia_I2.rds"
  ),
  "\n"
)

cat(
  "============================================\n"
)
