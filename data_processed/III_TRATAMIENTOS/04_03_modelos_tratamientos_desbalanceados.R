
## ============================================================
## 04_TRATAMIENTOS – SCRIPT 03
## MODELOS DE TRATAMIENTOS CON DISEÑO DESBALANCEADO
## ============================================================
##
## Entrada:
##   04_TRATAMIENTOS/02_INTEGRACION_SUELO_CULTIVO/objetos_R/
##   objetos_integracion_suelo_cultivo_tratamientos.rds
##
## Salida:
##   04_TRATAMIENTOS/03_MODELOS_TRATAMIENTOS/
##
## PROPÓSITO
##   1. Evaluar los efectos de dosis, omisión y repetición sobre
##      las variables observadas del cultivo.
##   2. Analizar de forma exploratoria las variables de suelo
##      predichas mediante Kriging en los 56 nodos.
##   3. Comparar tres formulaciones:
##        M0: respuesta ~ repetición
##        M1: respuesta ~ repetición + dosis + omisión
##        M2: respuesta ~ repetición + dosis * omisión
##   4. Priorizar magnitud del efecto, intervalos de confianza,
##      concordancia y desempeño predictivo; no usar p-valores
##      como criterio central.
##
## DECISIONES METODOLÓGICAS
##   - Replication se usa como factor fijo de bloqueo porque solo
##     presenta tres niveles y el diseño está desbalanceado.
##   - Se emplean errores estándar robustos HC3.
##   - La interacción dosis × omisión se conserva cuando mejora
##     el desempeño predictivo o la plausibilidad agronómica.
##   - Las variables de suelo interpoladas no se consideran
##     observaciones independientes de laboratorio; sus modelos
##     son exploratorios.
##   - N, P y K son índices instrumentales nominales.
##   - pH y Hv requieren interpretación cautelosa.
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
  "broom",
  "car",
  "effectsize",
  "emmeans",
  "sandwich",
  "performance",
  "patchwork"
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
library(writexl)
library(readr)
library(broom)
library(car)
library(effectsize)
library(emmeans)
library(sandwich)
library(performance)
library(patchwork)


## 1. Rutas ----------------------------------------------------

ruta_base <- paste0(
  "D:/UDCA_2017_CLASE/46_ARTICULOS_2025/",
  "ARTICULO_CEBADA_DOCTORADO/AJUSTES_ARTÍCULO"
)

ruta_objeto <- file.path(
  ruta_base,
  "04_TRATAMIENTOS",
  "02_INTEGRACION_SUELO_CULTIVO",
  "objetos_R",
  "objetos_integracion_suelo_cultivo_tratamientos.rds"
)

ruta_salida <- file.path(
  ruta_base,
  "04_TRATAMIENTOS",
  "03_MODELOS_TRATAMIENTOS"
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
  "datos_modelos"
)

ruta_objetos <- file.path(
  ruta_salida,
  "objetos_R"
)

if (!file.exists(ruta_objeto)) {
  stop(
    "No existe el objeto del Script 02:\n",
    ruta_objeto
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

dir.create(
  ruta_objetos,
  recursive = TRUE,
  showWarnings = FALSE
)


## 2. Cargar y preparar base -----------------------------------

obj <- readRDS(
  ruta_objeto
)

datos <- obj$base_integrada

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
  "Suelo_N",
  "Suelo_P",
  "Suelo_K",
  "Suelo_pH",
  "Suelo_T",
  "Suelo_Hv",
  "Suelo_CE",
  "Suelo_CE_corregida_OLS"
)

variables_cultivo <- intersect(
  variables_cultivo,
  names(
    datos
  )
)

variables_suelo <- intersect(
  variables_suelo,
  names(
    datos
  )
)

variables_respuesta <- c(
  variables_cultivo,
  variables_suelo
)

if (length(variables_respuesta) == 0) {
  stop(
    "No se encontraron variables de respuesta."
  )
}

datos <- datos |>
  mutate(
    Replication = factor(
      Replication
    ),
    Dose_percent = factor(
      Dose_percent,
      levels = c(
        75,
        100
      ),
      labels = c(
        "D75",
        "D100"
      )
    ),
    Omission_description = factor(
      Omission_description,
      levels = c(
        "Control",
        "Completo",
        "Sin N",
        "Sin P",
        "Sin K",
        "Sin Mg y S",
        "Sin Fe",
        "Sin Zn",
        "Sin Mn",
        "Sin B",
        "Sin Cu",
        "Sin Mo"
      )
    )
  )

# Contrastes suma-cero para interpretar modelos desbalanceados.
contrasts(
  datos$Dose_percent
) <- contr.sum(
  nlevels(
    datos$Dose_percent
  )
)

contrasts(
  datos$Omission_description
) <- contr.sum(
  nlevels(
    datos$Omission_description
  )
)

contrasts(
  datos$Replication
) <- contr.sum(
  nlevels(
    datos$Replication
  )
)


## 3. Catálogo de respuestas -----------------------------------

catalogo_respuestas <- tibble(
  Variable = variables_respuesta,
  Grupo = c(
    rep(
      "Cultivo_observado",
      length(
        variables_cultivo
      )
    ),
    rep(
      "Suelo_predicho_exploratorio",
      length(
        variables_suelo
      )
    )
  )
) |>
  mutate(
    Interpretacion = case_when(
      Variable %in% c(
        "Suelo_N",
        "Suelo_P",
        "Suelo_K"
      ) ~
        "Índice instrumental nominal; no equivale a análisis químico independiente",

      Variable == "Suelo_pH" ~
        "Exploratorio; sensor de pH con validación limitada",

      Variable == "Suelo_Hv" ~
        "Exploratorio; estructura espacial débil o efecto pepita dominante",

      Variable == "Suelo_CE_corregida_OLS" ~
        "CE calibrada exploratoriamente; la varianza no incluye incertidumbre de calibración",

      Variable == "Suelo_CE" ~
        "Lectura proximal original de CE",

      TRUE ~
        "Respuesta observada del cultivo"
    )
  )


## 4. Funciones auxiliares -------------------------------------

metricas_prediccion <- function(
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

  residuo <- observado -
    predicho

  tibble(
    n = length(
      observado
    ),
    ME = mean(
      residuo
    ),
    MAE = mean(
      abs(
        residuo
      )
    ),
    RMSE = sqrt(
      mean(
        residuo^2
      )
    ),
    Pearson = suppressWarnings(
      cor(
        observado,
        predicho,
        method = "pearson"
      )
    ),
    Spearman = suppressWarnings(
      cor(
        observado,
        predicho,
        method = "spearman"
      )
    )
  )
}


calcular_aicc <- function(
  modelo
) {

  k <- attr(
    logLik(
      modelo
    ),
    "df"
  )

  n <- nobs(
    modelo
  )

  aic <- AIC(
    modelo
  )

  if (
    n -
      k -
      1 <= 0
  ) {
    return(
      NA_real_
    )
  }

  aic +
    (
      2 *
        k *
        (
          k +
            1
        )
    ) /
    (
      n -
        k -
        1
    )
}


loocv_lm <- function(
  formula_modelo,
  datos_modelo,
  variable
) {

  predicciones <- rep(
    NA_real_,
    nrow(
      datos_modelo
    )
  )

  observado <- datos_modelo[[variable]]

  for (
    i in seq_len(
      nrow(
        datos_modelo
      )
    )
  ) {

    entrenamiento <- datos_modelo[
      -i,
      ,
      drop = FALSE
    ]

    prueba <- datos_modelo[
      i,
      ,
      drop = FALSE
    ]

    ajuste <- try(
      lm(
        formula_modelo,
        data = entrenamiento
      ),
      silent = TRUE
    )

    if (
      inherits(
        ajuste,
        "try-error"
      )
    ) {
      next
    }

    pred <- try(
      predict(
        ajuste,
        newdata = prueba
      ),
      silent = TRUE
    )

    if (
      !inherits(
        pred,
        "try-error"
      )
    ) {
      predicciones[i] <- as.numeric(
        pred
      )
    }
  }

  metricas_prediccion(
    observado = observado,
    predicho = predicciones
  )
}


extraer_anova_robusta <- function(
  modelo,
  variable,
  codigo_modelo
) {

  resultado <- try(
    car::Anova(
      modelo,
      type = 3,
      white.adjust = "hc3"
    ),
    silent = TRUE
  )

  if (
    inherits(
      resultado,
      "try-error"
    )
  ) {
    return(
      tibble(
        Variable = variable,
        Modelo = codigo_modelo,
        Termino = NA_character_,
        Df = NA_real_,
        Estadistico_F = NA_real_,
        p_valor_informativo = NA_real_
      )
    )
  }

  tabla <- as.data.frame(
    resultado
  )

  tabla$Termino <- rownames(
    tabla
  )

  rownames(
    tabla
  ) <- NULL

  nombres <- names(
    tabla
  )

  columna_df <- nombres[
    str_detect(
      nombres,
      "^Df$|Df"
    )
  ][1]

  columna_f <- nombres[
    str_detect(
      nombres,
      "F value|F"
    )
  ][1]

  columna_p <- nombres[
    str_detect(
      nombres,
      "Pr\\(>F\\)"
    )
  ][1]

  tibble(
    Variable = variable,
    Modelo = codigo_modelo,
    Termino = tabla$Termino,
    Df = if (
      !is.na(
        columna_df
      )
    ) {
      tabla[[columna_df]]
    } else {
      NA_real_
    },
    Estadistico_F = if (
      !is.na(
        columna_f
      )
    ) {
      tabla[[columna_f]]
    } else {
      NA_real_
    },
    p_valor_informativo = if (
      !is.na(
        columna_p
      )
    ) {
      tabla[[columna_p]]
    } else {
      NA_real_
    }
  )
}


extraer_tamanos_efecto <- function(
  modelo,
  variable,
  codigo_modelo
) {

  eta <- try(
    effectsize::eta_squared(
      modelo,
      partial = TRUE,
      ci = 0.95
    ),
    silent = TRUE
  )

  omega <- try(
    effectsize::omega_squared(
      modelo,
      partial = TRUE,
      ci = 0.95
    ),
    silent = TRUE
  )

  eta_tabla <- if (
    inherits(
      eta,
      "try-error"
    )
  ) {
    tibble()
  } else {
    as.data.frame(
      eta
    ) |>
      as_tibble() |>
      mutate(
        Medida = "Eta2_parcial"
      )
  }

  omega_tabla <- if (
    inherits(
      omega,
      "try-error"
    )
  ) {
    tibble()
  } else {
    as.data.frame(
      omega
    ) |>
      as_tibble() |>
      mutate(
        Medida = "Omega2_parcial"
      )
  }

  bind_rows(
    eta_tabla,
    omega_tabla
  ) |>
    mutate(
      Variable = variable,
      Modelo = codigo_modelo,
      .before = 1
    )
}


obtener_emmeans <- function(
  modelo,
  variable
) {

  vcov_robusta <- sandwich::vcovHC(
    modelo,
    type = "HC3"
  )

  medias_dosis <- try(
    emmeans(
      modelo,
      ~ Dose_percent,
      vcov. = vcov_robusta
    ),
    silent = TRUE
  )

  medias_omision <- try(
    emmeans(
      modelo,
      ~ Omission_description,
      vcov. = vcov_robusta
    ),
    silent = TRUE
  )

  medias_interaccion <- try(
    emmeans(
      modelo,
      ~ Dose_percent *
        Omission_description,
      vcov. = vcov_robusta
    ),
    silent = TRUE
  )

  contraste_dosis_por_omision <- try(
    contrast(
      emmeans(
        modelo,
        ~ Dose_percent |
          Omission_description,
        vcov. = vcov_robusta
      ),
      method = "revpairwise",
      adjust = "none"
    ),
    silent = TRUE
  )

  convertir <- function(
    objeto,
    tipo
  ) {

    if (
      inherits(
        objeto,
        "try-error"
      )
    ) {
      return(
        tibble()
      )
    }

    as.data.frame(
      objeto
    ) |>
      as_tibble() |>
      mutate(
        Variable = variable,
        Tipo = tipo,
        .before = 1
      )
  }

  list(
    dosis = convertir(
      medias_dosis,
      "Medias_ajustadas_dosis"
    ),
    omision = convertir(
      medias_omision,
      "Medias_ajustadas_omision"
    ),
    interaccion = convertir(
      medias_interaccion,
      "Medias_ajustadas_dosis_omision"
    ),
    contrastes = convertir(
      contraste_dosis_por_omision,
      "Contraste_D100_vs_D75_por_omision"
    )
  )
}


## 5. Ajustar modelos por variable -----------------------------

resultados_modelos <- list()
metricas_modelos <- list()
anova_robusta <- list()
tamanos_efecto <- list()
emmeans_dosis <- list()
emmeans_omision <- list()
emmeans_interaccion <- list()
contrastes_dosis <- list()
datos_diagnostico <- list()
errores_modelos <- list()

for (
  variable in variables_respuesta
) {

  message(
    "Procesando ",
    variable
  )

  datos_variable <- datos |>
    filter(
      is.finite(
        .data[[variable]]
      )
    )

  if (
    nrow(
      datos_variable
    ) < 20
  ) {

    errores_modelos[[variable]] <- tibble(
      Variable = variable,
      Error = "Menos de 20 observaciones válidas"
    )

    next
  }

  formula_m0 <- as.formula(
    paste0(
      variable,
      " ~ Replication"
    )
  )

  formula_m1 <- as.formula(
    paste0(
      variable,
      " ~ Replication + Dose_percent + Omission_description"
    )
  )

  formula_m2 <- as.formula(
    paste0(
      variable,
      " ~ Replication + Dose_percent * Omission_description"
    )
  )

  modelos <- list(
    M0_Bloque = lm(
      formula_m0,
      data = datos_variable
    ),
    M1_Aditivo = lm(
      formula_m1,
      data = datos_variable
    ),
    M2_Interaccion = lm(
      formula_m2,
      data = datos_variable
    )
  )

  resultados_modelos[[variable]] <- modelos

  formulas <- list(
    M0_Bloque = formula_m0,
    M1_Aditivo = formula_m1,
    M2_Interaccion = formula_m2
  )

  tabla_metricas_variable <- imap_dfr(
    modelos,
    function(
      modelo,
      codigo
    ) {

      resumen <- summary(
        modelo
      )

      cv <- loocv_lm(
        formula_modelo = formulas[[codigo]],
        datos_modelo = datos_variable,
        variable = variable
      )

      tibble(
        Variable = variable,
        Modelo = codigo,
        n = nobs(
          modelo
        ),
        Parametros = attr(
          logLik(
            modelo
          ),
          "df"
        ),
        R2 = resumen$r.squared,
        R2_ajustado = resumen$adj.r.squared,
        Sigma_residual = sigma(
          modelo
        ),
        AIC = AIC(
          modelo
        ),
        AICc = calcular_aicc(
          modelo
        ),
        BIC = BIC(
          modelo
        )
      ) |>
        bind_cols(
          cv |>
            rename(
              n_LOOCV = n,
              ME_LOOCV = ME,
              MAE_LOOCV = MAE,
              RMSE_LOOCV = RMSE,
              Pearson_LOOCV = Pearson,
              Spearman_LOOCV = Spearman
            )
        )
    }
  )

  # Selección basada en RMSE LOOCV y complejidad.
  tabla_metricas_variable <- tabla_metricas_variable |>
    mutate(
      Rank_RMSE = min_rank(
        RMSE_LOOCV
      ),
      Rank_AICc = min_rank(
        AICc
      ),
      Penalizacion_complejidad = case_when(
        Modelo == "M0_Bloque" ~ 0,
        Modelo == "M1_Aditivo" ~ 0.5,
        Modelo == "M2_Interaccion" ~ 1,
        TRUE ~ 0
      ),
      Score_seleccion =
        0.65 *
          Rank_RMSE +
        0.25 *
          Rank_AICc +
        0.10 *
          Penalizacion_complejidad
    ) |>
    arrange(
      Score_seleccion,
      RMSE_LOOCV
    ) |>
    mutate(
      Seleccionado =
        row_number() == 1
    )

  metricas_modelos[[variable]] <- tabla_metricas_variable

  codigo_seleccionado <- tabla_metricas_variable |>
    filter(
      Seleccionado
    ) |>
    pull(
      Modelo
    )

  modelo_seleccionado <- modelos[[codigo_seleccionado]]

  anova_robusta[[variable]] <- extraer_anova_robusta(
    modelo = modelo_seleccionado,
    variable = variable,
    codigo_modelo = codigo_seleccionado
  )

  tamanos_efecto[[variable]] <- extraer_tamanos_efecto(
    modelo = modelo_seleccionado,
    variable = variable,
    codigo_modelo = codigo_seleccionado
  )

  emms <- obtener_emmeans(
    modelo = modelo_seleccionado,
    variable = variable
  )

  emmeans_dosis[[variable]] <- emms$dosis

  emmeans_omision[[variable]] <- emms$omision

  emmeans_interaccion[[variable]] <- emms$interaccion

  contrastes_dosis[[variable]] <- emms$contrastes

  diagnostico <- augment(
    modelo_seleccionado,
    data = datos_variable
  ) |>
    transmute(
      Nodo,
      Treatment_ID,
      Dose_percent,
      Omission_description,
      Replication,
      Observado = .data[[variable]],
      Ajustado = .fitted,
      Residuo = .resid,
      Residuo_estudentizado = .std.resid,
      Leverage = .hat,
      Cook = .cooksd,
      Variable = variable,
      Modelo = codigo_seleccionado
    )

  datos_diagnostico[[variable]] <- diagnostico
}


## 6. Consolidar resultados ------------------------------------

tabla_metricas <- bind_rows(
  metricas_modelos
) |>
  left_join(
    catalogo_respuestas,
    by = "Variable"
  )

tabla_anova_robusta <- bind_rows(
  anova_robusta
) |>
  left_join(
    catalogo_respuestas,
    by = "Variable"
  )

tabla_tamanos_efecto <- bind_rows(
  tamanos_efecto
) |>
  left_join(
    catalogo_respuestas,
    by = "Variable"
  )

tabla_emmeans_dosis <- bind_rows(
  emmeans_dosis
)

tabla_emmeans_omision <- bind_rows(
  emmeans_omision
)

tabla_emmeans_interaccion <- bind_rows(
  emmeans_interaccion
)

tabla_contrastes_dosis <- bind_rows(
  contrastes_dosis
)

tabla_diagnostico <- bind_rows(
  datos_diagnostico
)

tabla_errores <- bind_rows(
  errores_modelos
)

modelos_seleccionados <- tabla_metricas |>
  filter(
    Seleccionado
  ) |>
  select(
    Variable,
    Grupo,
    Modelo,
    n,
    R2,
    R2_ajustado,
    RMSE_LOOCV,
    MAE_LOOCV,
    Pearson_LOOCV,
    Spearman_LOOCV,
    AICc,
    Interpretacion
  ) |>
  arrange(
    Grupo,
    RMSE_LOOCV
  )

print(
  modelos_seleccionados
)


## 7. Figuras de medias ajustadas ------------------------------

crear_figura_emmeans <- function(
  variable
) {

  datos_emm <- tabla_emmeans_interaccion |>
    filter(
      Variable == variable
    )

  if (
    nrow(
      datos_emm
    ) == 0
  ) {
    return(
      NULL
    )
  }

  columna_media <- intersect(
    c(
      "emmean",
      "response"
    ),
    names(
      datos_emm
    )
  )[1]

  columna_inferior <- intersect(
    c(
      "lower.CL",
      "asymp.LCL"
    ),
    names(
      datos_emm
    )
  )[1]

  columna_superior <- intersect(
    c(
      "upper.CL",
      "asymp.UCL"
    ),
    names(
      datos_emm
    )
  )[1]

  if (
    is.na(
      columna_media
    ) ||
      is.na(
        columna_inferior
      ) ||
      is.na(
        columna_superior
      )
  ) {
    return(
      NULL
    )
  }

  ggplot(
    datos_emm,
    aes(
      x = Omission_description,
      y = .data[[columna_media]],
      group = Dose_percent,
      shape = Dose_percent
    )
  ) +
    geom_point(
      position = position_dodge(
        width = 0.35
      ),
      size = 2.2
    ) +
    geom_errorbar(
      aes(
        ymin = .data[[columna_inferior]],
        ymax = .data[[columna_superior]]
      ),
      position = position_dodge(
        width = 0.35
      ),
      width = 0.18
    ) +
    labs(
      title = variable,
      subtitle =
        "Medias ajustadas e intervalos robustos HC3",
      x = "Tratamiento de omisión",
      y = "Media ajustada",
      shape = "Dosis"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 1
      )
    )
}


figuras_emmeans <- map(
  variables_respuesta,
  crear_figura_emmeans
)

names(
  figuras_emmeans
) <- variables_respuesta

figuras_emmeans <- figuras_emmeans[
  !map_lgl(
    figuras_emmeans,
    is.null
  )
]

for (
  variable in names(
    figuras_emmeans
  )
) {

  ggsave(
    filename = file.path(
      ruta_figuras,
      paste0(
        "Medias_ajustadas_",
        variable,
        ".png"
      )
    ),
    plot = figuras_emmeans[[variable]],
    width = 10,
    height = 7,
    dpi = 400
  )
}


## 8. Panel cultivo observado ----------------------------------

figuras_cultivo <- figuras_emmeans[
  intersect(
    variables_cultivo,
    names(
      figuras_emmeans
    )
  )
]

if (
  length(
    figuras_cultivo
  ) > 0
) {

  panel_cultivo <- wrap_plots(
    figuras_cultivo,
    ncol = 2
  ) +
    plot_annotation(
      title =
        "Efectos ajustados de dosis y omisión sobre el cultivo"
    )

  ggsave(
    filename = file.path(
      ruta_figuras,
      "Panel_efectos_tratamientos_cultivo.png"
    ),
    plot = panel_cultivo,
    width = 15,
    height = 20,
    dpi = 400
  )
}


## 9. Panel suelo exploratorio ---------------------------------

figuras_suelo <- figuras_emmeans[
  intersect(
    variables_suelo,
    names(
      figuras_emmeans
    )
  )
]

if (
  length(
    figuras_suelo
  ) > 0
) {

  panel_suelo <- wrap_plots(
    figuras_suelo,
    ncol = 2
  ) +
    plot_annotation(
      title =
        "Efectos exploratorios sobre variables de suelo predichas",
      subtitle =
        "Las predicciones de Kriging no sustituyen mediciones independientes"
    )

  ggsave(
    filename = file.path(
      ruta_figuras,
      "Panel_efectos_tratamientos_suelo_exploratorio.png"
    ),
    plot = panel_suelo,
    width = 15,
    height = 21,
    dpi = 400
  )
}


## 10. Figuras de diagnóstico ----------------------------------

g_residuos <- ggplot(
  tabla_diagnostico,
  aes(
    x = Ajustado,
    y = Residuo
  )
) +
  geom_point() +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  facet_wrap(
    ~ Variable,
    scales = "free",
    ncol = 3
  ) +
  labs(
    title =
      "Residuos frente a valores ajustados",
    x = "Valor ajustado",
    y = "Residuo"
  ) +
  theme_minimal()

ggsave(
  filename = file.path(
    ruta_figuras,
    "Diagnostico_residuos_modelos_tratamientos.png"
  ),
  plot = g_residuos,
  width = 14,
  height = 16,
  dpi = 400
)

g_influencia <- ggplot(
  tabla_diagnostico,
  aes(
    x = Leverage,
    y = Cook
  )
) +
  geom_point() +
  facet_wrap(
    ~ Variable,
    scales = "free",
    ncol = 3
  ) +
  labs(
    title =
      "Influencia de observaciones en los modelos",
    x = "Leverage",
    y = "Distancia de Cook"
  ) +
  theme_minimal()

ggsave(
  filename = file.path(
    ruta_figuras,
    "Diagnostico_influencia_modelos_tratamientos.png"
  ),
  plot = g_influencia,
  width = 14,
  height = 16,
  dpi = 400
)


## 11. Exportar tablas -----------------------------------------

write_csv(
  tabla_diagnostico,
  file.path(
    ruta_datos,
    "Diagnosticos_modelos_tratamientos.csv"
  )
)

write_csv(
  tabla_contrastes_dosis,
  file.path(
    ruta_datos,
    "Contrastes_D100_vs_D75_por_omision.csv"
  )
)

write_xlsx(
  list(
    Catalogo_respuestas =
      catalogo_respuestas,

    Modelos_seleccionados =
      modelos_seleccionados,

    Comparacion_modelos =
      tabla_metricas,

    ANOVA_robusta_HC3 =
      tabla_anova_robusta,

    Tamanos_efecto =
      tabla_tamanos_efecto,

    Medias_dosis =
      tabla_emmeans_dosis,

    Medias_omision =
      tabla_emmeans_omision,

    Medias_dosis_omision =
      tabla_emmeans_interaccion,

    Contrastes_dosis_omision =
      tabla_contrastes_dosis,

    Diagnosticos =
      tabla_diagnostico,

    Errores =
      tabla_errores
  ),
  path = file.path(
    ruta_tablas,
    "Resultados_modelos_tratamientos.xlsx"
  )
)


## 12. Guardar objeto para Script 04 ---------------------------

saveRDS(
  list(
    datos =
      datos,

    catalogo_respuestas =
      catalogo_respuestas,

    resultados_modelos =
      resultados_modelos,

    modelos_seleccionados =
      modelos_seleccionados,

    tabla_metricas =
      tabla_metricas,

    tabla_anova_robusta =
      tabla_anova_robusta,

    tabla_tamanos_efecto =
      tabla_tamanos_efecto,

    tabla_emmeans_interaccion =
      tabla_emmeans_interaccion,

    tabla_contrastes_dosis =
      tabla_contrastes_dosis,

    tabla_diagnostico =
      tabla_diagnostico
  ),
  file = file.path(
    ruta_objetos,
    "objetos_modelos_tratamientos.rds"
  )
)


## 13. Registro de sesión --------------------------------------

capture.output(
  sessionInfo(),
  file = file.path(
    ruta_salida,
    "sessionInfo_04_03.txt"
  )
)


## 14. Resumen en consola --------------------------------------

cat(
  "\n============================================\n"
)

cat(
  "SCRIPT 04-03 FINALIZADO\n\n"
)

cat(
  "Modelos seleccionados:\n"
)

print(
  modelos_seleccionados |>
    select(
      Variable,
      Grupo,
      Modelo,
      R2_ajustado,
      RMSE_LOOCV,
      Spearman_LOOCV
    )
)

cat(
  "\nObjeto para Script 04:\n",
  file.path(
    ruta_objetos,
    "objetos_modelos_tratamientos.rds"
  ),
  "\n"
)

cat(
  "============================================\n"
)
