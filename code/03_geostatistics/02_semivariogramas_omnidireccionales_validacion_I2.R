
## ============================================================
## SCRIPT 02. SEMIVARIOGRAMAS OMNIDIRECCIONALES Y VALIDACIÓN – I2
## Artículo cebada: revisión mayor
##
## Entrada principal:
##   AJUSTES_ARTÍCULO/01_AUDITORIA_I2/datos_procesados/
##   Suelo_I2_auditado_CE_corregida.csv
##
## Polígono:
##   AJUSTES_ARTÍCULO/DATOS BASE/Lote_2.shp
##
## Salidas:
##   AJUSTES_ARTÍCULO/02_SEMIVARIOGRAMAS_I2/
##
## Alcance:
##   1. Construir semivariogramas empíricos omnidireccionales.
##   2. Documentar número de pares por intervalo.
##   3. Ajustar modelos Sph, Exp, Gau y Bes.
##   4. Detectar ajustes singulares o físicamente inconsistentes.
##   5. Validar por LOOCV con ME, MAE, RMSE, RMSSE y MSNE.
##   6. Seleccionar modelos candidatos con criterio multivariable.
##   7. Evaluar sensibilidad de Hv con y sin el nodo 10.
##
## IMPORTANTE:
##   - Este script NO evalúa anisotropía. Eso corresponde al Script 03.
##   - Este script NO genera mapas finales de Kriging. Eso corresponde
##     al Script 04.
##   - pH se procesa únicamente como diagnóstico exploratorio.
##   - No se usan p-valores como criterio de selección.
##   - La selección combina ajuste, magnitud del error, calibración de la
##     incertidumbre y consistencia física del modelo.
## ============================================================


## 0. Paquetes -------------------------------------------------

paquetes <- c(
  "sf",
  "gstat",
  "readr",
  "dplyr",
  "tidyr",
  "purrr",
  "ggplot2",
  "writexl",
  "patchwork",
  "scales"
)

instalar <- paquetes[
  !paquetes %in% rownames(installed.packages())
]

if (length(instalar) > 0) {
  install.packages(instalar)
}

library(sf)
library(gstat)
library(readr)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(writexl)
library(patchwork)
library(scales)


## 1. Rutas ----------------------------------------------------

ruta_base <- paste0(
  "D:/UDCA_2017_CLASE/46_ARTICULOS_2025/",
  "ARTICULO_CEBADA_DOCTORADO/AJUSTES_ARTÍCULO"
)

ruta_datos <- file.path(
  ruta_base,
  "01_AUDITORIA_I2",
  "datos_procesados",
  "Suelo_I2_auditado_CE_corregida.csv"
)

ruta_lote <- file.path(
  ruta_base,
  "DATOS BASE",
  "Lote_2.shp"
)

ruta_salida <- file.path(
  ruta_base,
  "02_SEMIVARIOGRAMAS_I2"
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

if (!file.exists(ruta_datos)) {
  stop(
    "No existe la base auditada:\n",
    ruta_datos,
    "\nEjecuta primero el Script 01."
  )
}

if (!file.exists(ruta_lote)) {
  stop(
    "No existe el shapefile del lote:\n",
    ruta_lote
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


## 2. Lectura de datos y polígono -------------------------------

suelo <- read_csv(
  ruta_datos,
  show_col_types = FALSE
)

lote <- st_read(
  ruta_lote,
  quiet = TRUE
)

# CRS objetivo confirmado para este flujo.
crs_objetivo <- 32618

if (is.na(st_crs(lote))) {
  warning(
    "El shapefile no tenía CRS. Se asignará EPSG:32618 ",
    "sin reproyectar. Verifica que esto sea correcto."
  )
  st_crs(lote) <- crs_objetivo
} else {
  lote <- st_transform(
    lote,
    crs_objetivo
  )
}

# Crear puntos desde X y Y.
suelo_sf <- st_as_sf(
  suelo,
  coords = c("X", "Y"),
  crs = crs_objetivo,
  remove = FALSE
)

# Verificar inclusión espacial.
dentro_lote <- lengths(
  st_intersects(
    suelo_sf,
    lote
  )
) > 0

tabla_inclusion <- suelo |>
  transmute(
    Nodo,
    X,
    Y,
    Dentro_lote = dentro_lote
  )

if (any(!dentro_lote)) {
  warning(
    sum(!dentro_lote),
    " puntos se encuentran fuera del lote. ",
    "No se eliminarán automáticamente."
  )
}


## 3. Variables y clasificación --------------------------------

variables_modelar <- c(
  "N",
  "P",
  "K",
  "pH",
  "T",
  "Hv",
  "CE",
  "CE_corregida_OLS"
)

clasificacion_variables <- tibble(
  Variable = variables_modelar,
  Grupo = c(
    "Índice instrumental nominal",
    "Índice instrumental nominal",
    "Índice instrumental nominal",
    "Diagnóstico exploratorio",
    "Variable física",
    "Variable física",
    "Variable proximal principal",
    "Estimación exploratoria calibrada"
  ),
  Uso_interpretativo = c(
    "Suplementario",
    "Suplementario",
    "Suplementario",
    "No usar para inferencia agronómica principal",
    "Principal con cautela",
    "Principal",
    "Principal",
    "Principal comparativo"
  )
)


## 4. Parámetros del semivariograma empírico --------------------

coords <- st_coordinates(
  suelo_sf
)

distancias <- as.vector(
  dist(coords)
)

distancia_maxima <- max(
  distancias,
  na.rm = TRUE
)

# Cutoff conservador: 50% de la distancia máxima.
# Reduce el dominio de distancias con pocos pares.
cutoff_global <- 0.50 * distancia_maxima

# Se eligen 7 intervalos para 28 nodos.
# Si la geometría produce bins vacíos, gstat los omitirá.
numero_lags <- 7
width_global <- cutoff_global / numero_lags

tabla_parametros_variograma <- tibble(
  Parametro = c(
    "Número de nodos",
    "Distancia máxima entre puntos (m)",
    "Cutoff usado (m)",
    "Número objetivo de intervalos",
    "Ancho de intervalo (m)"
  ),
  Valor = c(
    nrow(suelo_sf),
    distancia_maxima,
    cutoff_global,
    numero_lags,
    width_global
  )
)

print(tabla_parametros_variograma)


## 5. Modelos candidatos ---------------------------------------

modelos_candidatos <- c(
  "Sph",
  "Exp",
  "Gau",
  "Bes"
)

nombres_modelos <- c(
  Sph = "Esférico",
  Exp = "Exponencial",
  Gau = "Gaussiano",
  Bes = "Bessel"
)


## 6. Funciones auxiliares -------------------------------------

extraer_componentes_modelo <- function(vfit) {

  nugget <- vfit |>
    filter(model == "Nug") |>
    pull(psill)

  if (length(nugget) == 0) {
    nugget <- 0
  }

  estructuras <- vfit |>
    filter(model != "Nug")

  if (nrow(estructuras) == 0) {
    return(
      tibble(
        Pepita = nugget[1],
        Varianza_estructural = NA_real_,
        Silla_total = NA_real_,
        Rango_modelo = NA_real_,
        GDE = NA_real_
      )
    )
  }

  estructura <- sum(
    estructuras$psill,
    na.rm = TRUE
  )

  silla <- nugget[1] + estructura

  # Para este flujo se usa una sola estructura.
  rango <- estructuras$range[1]

  gde <- ifelse(
    silla > 0,
    estructura / silla,
    NA_real_
  )

  tibble(
    Pepita = nugget[1],
    Varianza_estructural = estructura,
    Silla_total = silla,
    Rango_modelo = rango,
    GDE = gde
  )
}


calcular_metricas_cv <- function(
  formula_modelo,
  datos_sf,
  modelo
) {

  cv <- try(
    krige.cv(
      formula_modelo,
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
        tabla = tibble(
          ME = NA_real_,
          MAE = NA_real_,
          RMSE = NA_real_,
          RMSSE = NA_real_,
          MSNE = NA_real_,
          Correlacion_observado_predicho = NA_real_
        ),
        datos = NULL,
        error = as.character(cv)
      )
    )
  }

  residuo <- cv$observed - cv$var1.pred

  varianza_valida <- is.finite(cv$var1.var) &
    cv$var1.var > 0

  residuo_estandarizado <- rep(
    NA_real_,
    length(residuo)
  )

  residuo_estandarizado[varianza_valida] <-
    residuo[varianza_valida] /
    sqrt(cv$var1.var[varianza_valida])

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
      any(varianza_valida),
      sqrt(
        mean(
          residuo_estandarizado^2,
          na.rm = TRUE
        )
      ),
      NA_real_
    ),
    MSNE = ifelse(
      any(varianza_valida),
      mean(
        residuo_estandarizado^2,
        na.rm = TRUE
      ),
      NA_real_
    ),
    Correlacion_observado_predicho = suppressWarnings(
      cor(
        cv$observed,
        cv$var1.pred,
        use = "complete.obs"
      )
    )
  )

  list(
    tabla = metricas,
    datos = cv,
    error = NA_character_
  )
}


evaluar_consistencia_modelo <- function(
  parametros,
  cutoff,
  varianza_datos,
  ajuste_singular
) {

  pepita <- parametros$Pepita
  silla <- parametros$Silla_total
  rango <- parametros$Rango_modelo

  razones <- c()

  if (isTRUE(ajuste_singular)) {
    razones <- c(
      razones,
      "Ajuste singular"
    )
  }

  if (
    !is.finite(pepita) ||
      pepita < 0
  ) {
    razones <- c(
      razones,
      "Pepita negativa/no finita"
    )
  }

  if (
    !is.finite(silla) ||
      silla <= 0
  ) {
    razones <- c(
      razones,
      "Silla no positiva/no finita"
    )
  }

  if (
    !is.finite(rango) ||
      rango <= 0
  ) {
    razones <- c(
      razones,
      "Rango no positivo/no finito"
    )
  }

  if (
    is.finite(rango) &&
      rango > 3 * cutoff
  ) {
    razones <- c(
      razones,
      "Rango muy superior al dominio observado"
    )
  }

  if (
    is.finite(silla) &&
      varianza_datos > 0 &&
      silla > 20 * varianza_datos
  ) {
    razones <- c(
      razones,
      "Silla excesiva frente a la varianza muestral"
    )
  }

  if (length(razones) == 0) {
    tibble(
      Modelo_valido = TRUE,
      Motivo_descarte = ""
    )
  } else {
    tibble(
      Modelo_valido = FALSE,
      Motivo_descarte = paste(
        razones,
        collapse = "; "
      )
    )
  }
}


ajustar_modelo_individual <- function(
  vexp,
  variable,
  datos_sf,
  modelo_codigo,
  formula_modelo,
  cutoff,
  width
) {

  varianza_variable <- var(
    datos_sf[[variable]],
    na.rm = TRUE
  )

  gamma_min <- min(
    vexp$gamma,
    na.rm = TRUE
  )

  gamma_max <- max(
    vexp$gamma,
    na.rm = TRUE
  )

  nugget_ini <- max(
    0,
    min(
      gamma_min,
      0.25 * varianza_variable
    )
  )

  psill_ini <- max(
    varianza_variable - nugget_ini,
    0.25 * varianza_variable,
    .Machine$double.eps
  )

  range_ini <- max(
    cutoff / 3,
    width
  )

  modelo_inicial <- vgm(
    psill = psill_ini,
    model = modelo_codigo,
    range = range_ini,
    nugget = nugget_ini
  )

  fit_obj <- try(
    fit.variogram(
      vexp,
      modelo_inicial,
      fit.method = 6
    ),
    silent = TRUE
  )

  if (inherits(fit_obj, "try-error")) {
    return(
      list(
        resumen = tibble(
          Variable = variable,
          Modelo_codigo = modelo_codigo,
          Modelo = nombres_modelos[[modelo_codigo]],
          SSErr_ajuste = NA_real_,
          Singular = NA,
          Pepita = NA_real_,
          Varianza_estructural = NA_real_,
          Silla_total = NA_real_,
          Rango_modelo = NA_real_,
          GDE = NA_real_,
          ME = NA_real_,
          MAE = NA_real_,
          RMSE = NA_real_,
          RMSSE = NA_real_,
          MSNE = NA_real_,
          Correlacion_observado_predicho = NA_real_,
          Modelo_valido = FALSE,
          Motivo_descarte = paste(
            "Error de ajuste:",
            as.character(fit_obj)
          )
        ),
        modelo = NULL,
        cv = NULL
      )
    )
  }

  singular <- isTRUE(
    attr(
      fit_obj,
      "singular"
    )
  )

  ss_err <- attr(
    fit_obj,
    "SSErr"
  )

  componentes <- extraer_componentes_modelo(
    fit_obj
  )

  consistencia <- evaluar_consistencia_modelo(
    componentes,
    cutoff,
    varianza_variable,
    singular
  )

  cv_resultado <- calcular_metricas_cv(
    formula_modelo,
    datos_sf,
    fit_obj
  )

  resumen <- bind_cols(
    tibble(
      Variable = variable,
      Modelo_codigo = modelo_codigo,
      Modelo = nombres_modelos[[modelo_codigo]],
      SSErr_ajuste = ifelse(
        length(ss_err) == 0,
        NA_real_,
        as.numeric(ss_err)
      ),
      Singular = singular
    ),
    componentes,
    cv_resultado$tabla,
    consistencia
  )

  list(
    resumen = resumen,
    modelo = fit_obj,
    cv = cv_resultado$datos
  )
}


crear_variograma_empirico <- function(
  variable,
  datos_sf,
  cutoff,
  width
) {

  formula_modelo <- as.formula(
    paste(
      variable,
      "~ 1"
    )
  )

  variogram(
    formula_modelo,
    datos_sf,
    cutoff = cutoff,
    width = width,
    cressie = FALSE
  )
}


seleccionar_modelo <- function(tabla_variable) {

  candidatos <- tabla_variable |>
    filter(
      Modelo_valido,
      is.finite(RMSE),
      is.finite(MSNE),
      is.finite(SSErr_ajuste)
    )

  if (nrow(candidatos) == 0) {
    return(
      tabla_variable |>
        slice(0) |>
        mutate(
          Score_seleccion = numeric(),
          Seleccionado = logical()
        )
    )
  }

  candidatos <- candidatos |>
    mutate(
      Rank_RMSE = min_rank(RMSE),
      Rank_MAE = min_rank(MAE),
      Rank_MSNE = min_rank(
        abs(MSNE - 1)
      ),
      Rank_ME = min_rank(
        abs(ME)
      ),
      Rank_SSErr = min_rank(
        SSErr_ajuste
      ),
      Score_seleccion =
        0.35 * Rank_RMSE +
        0.20 * Rank_MAE +
        0.20 * Rank_MSNE +
        0.10 * Rank_ME +
        0.15 * Rank_SSErr
    ) |>
    arrange(
      Score_seleccion,
      RMSE,
      abs(MSNE - 1)
    ) |>
    mutate(
      Seleccionado =
        row_number() == 1
    )

  candidatos
}


## 7. Ajuste completo por variable ------------------------------

resultados_variables <- list()
variogramas_empiricos <- list()
modelos_ajustados <- list()
validaciones_cv <- list()

for (variable in variables_modelar) {

  message(
    "Procesando variable: ",
    variable
  )

  formula_variable <- as.formula(
    paste(
      variable,
      "~ 1"
    )
  )

  vexp <- crear_variograma_empirico(
    variable,
    suelo_sf,
    cutoff_global,
    width_global
  )

  variogramas_empiricos[[variable]] <- vexp

  resultados_modelos <- map(
    modelos_candidatos,
    function(modelo_codigo) {

      ajustar_modelo_individual(
        vexp = vexp,
        variable = variable,
        datos_sf = suelo_sf,
        modelo_codigo = modelo_codigo,
        formula_modelo = formula_variable,
        cutoff = cutoff_global,
        width = width_global
      )
    }
  )

  tabla_variable <- bind_rows(
    map(
      resultados_modelos,
      "resumen"
    )
  )

  seleccion <- seleccionar_modelo(
    tabla_variable
  )

  tabla_variable <- tabla_variable |>
    left_join(
      seleccion |>
        select(
          Modelo_codigo,
          Score_seleccion,
          Seleccionado
        ),
      by = "Modelo_codigo"
    ) |>
    mutate(
      Score_seleccion = ifelse(
        is.na(Score_seleccion),
        Inf,
        Score_seleccion
      ),
      Seleccionado = ifelse(
        is.na(Seleccionado),
        FALSE,
        Seleccionado
      )
    )

  resultados_variables[[variable]] <- tabla_variable

  modelos_ajustados[[variable]] <- setNames(
    map(
      resultados_modelos,
      "modelo"
    ),
    modelos_candidatos
  )

  validaciones_cv[[variable]] <- setNames(
    map(
      resultados_modelos,
      "cv"
    ),
    modelos_candidatos
  )
}


tabla_todos_modelos <- bind_rows(
  resultados_variables
)

tabla_modelos_seleccionados <- tabla_todos_modelos |>
  filter(
    Seleccionado
  ) |>
  arrange(
    match(
      Variable,
      variables_modelar
    )
  )

print(tabla_modelos_seleccionados)


## 8. Tabla de pares por intervalo -----------------------------

tabla_pares_intervalo <- imap_dfr(
  variogramas_empiricos,
  function(vexp, variable) {

    as_tibble(vexp) |>
      transmute(
        Variable = variable,
        Intervalo = row_number(),
        Numero_pares = np,
        Distancia_media_m = dist,
        Semivarianza = gamma
      )
  }
)

print(tabla_pares_intervalo)


## 9. Sensibilidad de Hv con y sin nodo 10 ---------------------

suelo_hv_sin_10 <- suelo_sf |>
  filter(
    Nodo != 10
  )

vexp_hv_completo <- variogramas_empiricos[["Hv"]]

vexp_hv_sin_10 <- crear_variograma_empirico(
  "Hv",
  suelo_hv_sin_10,
  cutoff_global,
  width_global
)

ajustes_hv_completo <- map(
  modelos_candidatos,
  function(modelo_codigo) {

    ajustar_modelo_individual(
      vexp = vexp_hv_completo,
      variable = "Hv",
      datos_sf = suelo_sf,
      modelo_codigo = modelo_codigo,
      formula_modelo = Hv ~ 1,
      cutoff = cutoff_global,
      width = width_global
    )
  }
)

ajustes_hv_sin_10 <- map(
  modelos_candidatos,
  function(modelo_codigo) {

    ajustar_modelo_individual(
      vexp = vexp_hv_sin_10,
      variable = "Hv",
      datos_sf = suelo_hv_sin_10,
      modelo_codigo = modelo_codigo,
      formula_modelo = Hv ~ 1,
      cutoff = cutoff_global,
      width = width_global
    )
  }
)

tabla_hv_completo <- bind_rows(
  map(
    ajustes_hv_completo,
    "resumen"
  )
) |>
  mutate(
    Escenario = "Hv_completo_28",
    .before = 1
  )

tabla_hv_sin_10 <- bind_rows(
  map(
    ajustes_hv_sin_10,
    "resumen"
  )
) |>
  mutate(
    Escenario = "Hv_sin_nodo_10",
    .before = 1
  )

tabla_sensibilidad_hv <- bind_rows(
  tabla_hv_completo,
  tabla_hv_sin_10
)

seleccion_hv_completo <- seleccionar_modelo(
  tabla_hv_completo
)

seleccion_hv_sin_10 <- seleccionar_modelo(
  tabla_hv_sin_10
)

tabla_hv_seleccionados <- bind_rows(
  seleccion_hv_completo |>
    mutate(
      Escenario = "Hv_completo_28",
      .before = 1
    ),
  seleccion_hv_sin_10 |>
    mutate(
      Escenario = "Hv_sin_nodo_10",
      .before = 1
    )
) |>
  filter(
    Seleccionado
  )

print(tabla_hv_seleccionados)


## 10. Figuras: modelos por variable ---------------------------

crear_dataframe_curva <- function(
  variable,
  modelo_obj,
  vexp
) {

  if (is.null(modelo_obj)) {
    return(NULL)
  }

  distancia_max <- max(
    vexp$dist,
    na.rm = TRUE
  )

  curva <- variogramLine(
    modelo_obj,
    maxdist = distancia_max,
    n = 250
  )

  as_tibble(curva) |>
    mutate(
      Variable = variable
    )
}


graficos_variogramas <- list()

for (variable in variables_modelar) {

  seleccionado <- tabla_modelos_seleccionados |>
    filter(
      Variable == variable
    )

  vexp <- variogramas_empiricos[[variable]]

  if (nrow(seleccionado) == 0) {

    g <- ggplot(
      as_tibble(vexp),
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
        title = variable,
        subtitle = "Sin modelo candidato válido",
        x = "Distancia (m)",
        y = "Semivarianza",
        size = "Pares"
      ) +
      theme_minimal()

    graficos_variogramas[[variable]] <- g
    next
  }

  modelo_codigo <- seleccionado$Modelo_codigo[1]
  
  modelo_obj <-
    modelos_ajustados[[variable]][[modelo_codigo]]

  curva <- crear_dataframe_curva(
    variable,
    modelo_obj,
    vexp
  )

  if (is.null(curva)) {

    g <- ggplot(
      as_tibble(vexp),
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
        title = variable,
        subtitle = "El modelo seleccionado no conservó una curva válida",
        x = "Distancia (m)",
        y = "Semivarianza",
        size = "Pares"
      ) +
      theme_minimal()

    graficos_variogramas[[variable]] <- g
    next
  }

  subtitulo <- paste0(
    seleccionado$Modelo[1],
    " | Pepita=",
    round(
      seleccionado$Pepita[1],
      3
    ),
    " | Silla=",
    round(
      seleccionado$Silla_total[1],
      3
    ),
    " | Rango=",
    round(
      seleccionado$Rango_modelo[1],
      2
    ),
    " m | RMSE=",
    round(
      seleccionado$RMSE[1],
      3
    ),
    " | MSNE=",
    round(
      seleccionado$MSNE[1],
      3
    )
  )

  g <- ggplot(
    as_tibble(vexp),
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
    geom_line(
      data = curva,
      aes(
        x = dist,
        y = gamma
      ),
      linewidth = 0.8,
      inherit.aes = FALSE
    ) +
    scale_size_continuous(
      range = c(2, 6)
    ) +
    labs(
      title = variable,
      subtitle = subtitulo,
      x = "Distancia (m)",
      y = "Semivarianza",
      size = "Pares"
    ) +
    theme_minimal()

  graficos_variogramas[[variable]] <- g
}


panel_variogramas <- wrap_plots(
  graficos_variogramas,
  ncol = 2
)

print(panel_variogramas)


## 11. Figura comparativa CE original y corregida --------------

g_ce_original <- graficos_variogramas[["CE"]]
g_ce_corregida <- graficos_variogramas[["CE_corregida_OLS"]]

panel_ce <- g_ce_original |
  g_ce_corregida

print(panel_ce)


## 12. Figura sensibilidad Hv ----------------------------------

crear_grafico_hv <- function(
  vexp,
  ajustes,
  tabla_seleccion,
  titulo
) {

  if (
    nrow(tabla_seleccion) == 0
  ) {
    return(
      ggplot(
        as_tibble(vexp),
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
          title = titulo,
          subtitle = "Sin modelo válido",
          x = "Distancia (m)",
          y = "Semivarianza"
        ) +
        theme_minimal()
    )
  }

  codigo <- tabla_seleccion$Modelo_codigo[1]

  indice_modelo <- match(
    codigo,
    modelos_candidatos
  )

  modelo <- ajustes[[indice_modelo]]$modelo

  curva <- variogramLine(
    modelo,
    maxdist = max(vexp$dist),
    n = 250
  )

  ggplot(
    as_tibble(vexp),
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
    geom_line(
      data = curva,
      aes(
        x = dist,
        y = gamma
      ),
      inherit.aes = FALSE
    ) +
    labs(
      title = titulo,
      subtitle = paste0(
        tabla_seleccion$Modelo[1],
        " | Rango=",
        round(
          tabla_seleccion$Rango_modelo[1],
          2
        ),
        " m | RMSE=",
        round(
          tabla_seleccion$RMSE[1],
          3
        ),
        " | MSNE=",
        round(
          tabla_seleccion$MSNE[1],
          3
        )
      ),
      x = "Distancia (m)",
      y = "Semivarianza",
      size = "Pares"
    ) +
    theme_minimal()
}


g_hv_completo <- crear_grafico_hv(
  vexp_hv_completo,
  ajustes_hv_completo,
  seleccion_hv_completo |>
    filter(
      Seleccionado
    ),
  "Hv con 28 nodos"
)

g_hv_sin_10 <- crear_grafico_hv(
  vexp_hv_sin_10,
  ajustes_hv_sin_10,
  seleccion_hv_sin_10 |>
    filter(
      Seleccionado
    ),
  "Hv sin nodo 10"
)

panel_hv <- g_hv_completo |
  g_hv_sin_10

print(panel_hv)


## 13. Exportación de figuras ----------------------------------

ggsave(
  filename = file.path(
    ruta_figuras,
    "01_semivariogramas_omnidireccionales_I2.png"
  ),
  plot = panel_variogramas,
  width = 14,
  height = 18,
  dpi = 400
)

ggsave(
  filename = file.path(
    ruta_figuras,
    "01_semivariogramas_omnidireccionales_I2.tiff"
  ),
  plot = panel_variogramas,
  width = 14,
  height = 18,
  dpi = 600,
  compression = "lzw"
)

ggsave(
  filename = file.path(
    ruta_figuras,
    "02_comparacion_semivariogramas_CE.png"
  ),
  plot = panel_ce,
  width = 12,
  height = 6,
  dpi = 400
)

ggsave(
  filename = file.path(
    ruta_figuras,
    "03_sensibilidad_semivariograma_Hv.png"
  ),
  plot = panel_hv,
  width = 12,
  height = 6,
  dpi = 400
)


## 14. Exportación de tablas -----------------------------------

write_xlsx(
  list(
    Parametros_generales =
      tabla_parametros_variograma,

    Clasificacion_variables =
      clasificacion_variables,

    Inclusion_lote =
      tabla_inclusion,

    Pares_por_intervalo =
      tabla_pares_intervalo,

    Todos_modelos =
      tabla_todos_modelos,

    Modelos_seleccionados =
      tabla_modelos_seleccionados,

    Sensibilidad_Hv_todos =
      tabla_sensibilidad_hv,

    Sensibilidad_Hv_seleccion =
      tabla_hv_seleccionados
  ),
  path = file.path(
    ruta_tablas,
    "Resultados_semivariogramas_I2.xlsx"
  )
)


## 15. Guardar objetos para scripts siguientes -----------------

saveRDS(
  list(
    suelo_sf = suelo_sf,
    lote = lote,
    variogramas_empiricos =
      variogramas_empiricos,
    modelos_ajustados =
      modelos_ajustados,
    validaciones_cv =
      validaciones_cv,
    tabla_todos_modelos =
      tabla_todos_modelos,
    tabla_modelos_seleccionados =
      tabla_modelos_seleccionados,
    cutoff =
      cutoff_global,
    width =
      width_global,
    variables =
      variables_modelar
  ),
  file = file.path(
    ruta_objetos,
    "objetos_semivariogramas_I2.rds"
  )
)


## 16. Registro de sesión --------------------------------------

capture.output(
  sessionInfo(),
  file = file.path(
    ruta_salida,
    "sessionInfo_Script02.txt"
  )
)


## 17. Resumen en consola --------------------------------------

cat(
  "\n============================================\n"
)

cat(
  "SCRIPT 02 FINALIZADO\n\n"
)

cat(
  "Modelos seleccionados:\n"
)

print(
  tabla_modelos_seleccionados |>
    select(
      Variable,
      Modelo,
      Pepita,
      Silla_total,
      Rango_modelo,
      GDE,
      RMSE,
      MSNE,
      Correlacion_observado_predicho
    )
)

cat(
  "\nResultados guardados en:\n",
  ruta_salida,
  "\n"
)

cat(
  "\nObjeto para Script 03:\n",
  file.path(
    ruta_objetos,
    "objetos_semivariogramas_I2.rds"
  ),
  "\n"
)

cat(
  "============================================\n"
)
