
## ============================================================
## 05_PATRON_PUNTOS – SCRIPT 05-C
## ASOCIACIÓN ESPACIAL ENTRE VARIABLES DE SUELO Y CULTIVO
## ============================================================
##
## Entrada:
##   05_PATRON_PUNTOS/02_PATRONES_POR_GRUPO/objetos_R/
##   objetos_patrones_por_grupo.rds
##
## Salida:
##   05_PATRON_PUNTOS/03_ASOCIACION_SUELO_CULTIVO/
##
## PREGUNTA
##   ¿Las variables del cultivo presentan correspondencia espacial
##   con las condiciones del suelo predichas en los mismos 56 nodos?
##
## ENFOQUE
##   1. Correlaciones Pearson y Spearman en valores originales.
##   2. Correlaciones sobre residuos ajustados por:
##        repetición + dosis + omisión.
##   3. Moran bivariado global:
##        suelo_i frente al rezago espacial del cultivo_j.
##   4. Referencia Monte Carlo para Moran bivariado.
##   5. Asociación por clases de distancia mediante productos
##      cruzados estandarizados.
##   6. Comparación específica entre CE original y CE corregida OLS.
##
## PRECAUCIONES
##   - Los valores de suelo fueron predichos mediante Kriging en las
##     coordenadas del cultivo; no son mediciones independientes.
##   - N, P y K son índices instrumentales nominales.
##   - pH y Hv tienen alcance exploratorio.
##   - CE corregida OLS incorpora una calibración exploratoria.
##   - La interpretación se centra en magnitud, estabilidad,
##     dirección y coherencia espacial, no en p-valores aislados.
##   - Moran bivariado no implica causalidad.
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
  "spdep",
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
library(readr)
library(writexl)
library(sf)
library(spdep)
library(patchwork)


## 1. Parámetros reproducibles ---------------------------------

set.seed(2026)

nsim_moran <- 999
k_vecinos <- 4
n_clases_distancia <- 6


## 2. Rutas ----------------------------------------------------

ruta_base <- paste0(
  "D:/UDCA_2017_CLASE/46_ARTICULOS_2025/",
  "ARTICULO_CEBADA_DOCTORADO/AJUSTES_ARTÍCULO"
)

ruta_objeto <- file.path(
  ruta_base,
  "05_PATRON_PUNTOS",
  "02_PATRONES_POR_GRUPO",
  "objetos_R",
  "objetos_patrones_por_grupo.rds"
)

ruta_salida <- file.path(
  ruta_base,
  "05_PATRON_PUNTOS",
  "03_ASOCIACION_SUELO_CULTIVO"
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

if (!file.exists(ruta_objeto)) {
  stop(
    "No existe el objeto del Script 05-B:\n",
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


## 3. Cargar y preparar datos ----------------------------------

obj <- readRDS(
  ruta_objeto
)

datos <- obj$datos
cultivo_sf <- obj$cultivo_sf
lote <- obj$lote

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

variables_cultivo <- intersect(
  variables_cultivo,
  names(datos)
)

variables_suelo <- intersect(
  variables_suelo,
  names(datos)
)

if (length(variables_cultivo) == 0) {
  stop(
    "No se encontraron variables de cultivo."
  )
}

if (length(variables_suelo) == 0) {
  stop(
    "No se encontraron variables de suelo."
  )
}

columnas_requeridas <- c(
  "Nodo",
  "X",
  "Y",
  "Replication",
  "Dose_percent",
  "Omission_description"
)

faltantes <- setdiff(
  columnas_requeridas,
  names(datos)
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

datos <- datos |>
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


## 4. Pesos espaciales k-NN ------------------------------------

coords <- as.matrix(
  datos |>
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

tabla_vecinos <- tibble(
  Nodo = datos$Nodo,
  Numero_vecinos = card(
    vecinos
  )
)


## 5. Funciones auxiliares -------------------------------------

estandarizar <- function(x) {

  x <- as.numeric(
    x
  )

  if (
    all(
      !is.finite(
        x
      )
    )
  ) {
    return(
      rep(
        NA_real_,
        length(
          x
        )
      )
    )
  }

  desviacion <- sd(
    x,
    na.rm = TRUE
  )

  if (
    !is.finite(
      desviacion
    ) ||
      desviacion == 0
  ) {
    return(
      rep(
        0,
        length(
          x
        )
      )
    )
  }

  as.numeric(
    scale(
      x
    )
  )
}


residualizar <- function(
  variable,
  df
) {

  formula_modelo <- as.formula(
    paste0(
      variable,
      " ~ Replication + Dose_percent + Omission_description"
    )
  )

  modelo <- try(
    lm(
      formula_modelo,
      data = df,
      na.action = na.exclude
    ),
    silent = TRUE
  )

  if (
    inherits(
      modelo,
      "try-error"
    )
  ) {
    return(
      rep(
        NA_real_,
        nrow(
          df
        )
      )
    )
  }

  residuales <- residuals(
    modelo
  )

  salida <- rep(
    NA_real_,
    nrow(
      df
    )
  )

  indices <- as.integer(
    names(
      residuales
    )
  )

  salida[
    indices
  ] <- as.numeric(
    residuales
  )

  salida
}


correlacion_segura <- function(
  x,
  y,
  metodo
) {

  validos <- is.finite(
    x
  ) &
    is.finite(
      y
    )

  if (
    sum(
      validos
    ) < 4
  ) {
    return(
      NA_real_
    )
  }

  suppressWarnings(
    cor(
      x[
        validos
      ],
      y[
        validos
      ],
      method = metodo
    )
  )
}


moran_bivariado <- function(
  x,
  y,
  listw
) {

  zx <- estandarizar(
    x
  )

  zy <- estandarizar(
    y
  )

  validos <- is.finite(
    zx
  ) &
    is.finite(
      zy
    )

  if (
    !all(
      validos
    )
  ) {
    return(
      NA_real_
    )
  }

  rezago_y <- lag.listw(
    listw,
    zy,
    zero.policy = TRUE
  )

  mean(
    zx *
      rezago_y
  )
}


simular_moran_bivariado <- function(
  x,
  y,
  listw,
  nsim
) {

  observado <- moran_bivariado(
    x,
    y,
    listw
  )

  simulados <- map_dbl(
    seq_len(
      nsim
    ),
    function(i) {

      y_perm <- sample(
        y,
        replace = FALSE
      )

      moran_bivariado(
        x,
        y_perm,
        listw
      )
    }
  )

  q <- quantile(
    simulados,
    probs = c(
      0.025,
      0.50,
      0.975
    ),
    na.rm = TRUE,
    names = FALSE
  )

  p_mc <- (
    1 +
      sum(
        abs(
          simulados -
            median(
              simulados,
              na.rm = TRUE
            )
        ) >=
          abs(
            observado -
              median(
                simulados,
                na.rm = TRUE
              )
          ),
        na.rm = TRUE
      )
  ) /
    (
      nsim +
        1
    )

  clasificacion <- case_when(
    !is.finite(
      observado
    ) ~
      "No interpretable",

    observado >
      q[3] ~
      "Asociación espacial positiva mayor que la referencia",

    observado <
      q[1] ~
      "Asociación espacial negativa menor que la referencia",

    TRUE ~
      "Compatible con la referencia permutada"
  )

  list(
    observado = observado,
    q025 = q[1],
    mediana = q[2],
    q975 = q[3],
    p_mc = p_mc,
    clasificacion = clasificacion,
    simulados = simulados
  )
}


## 6. Crear residuos ajustados ---------------------------------

datos_residuales <- datos |>
  select(
    Nodo,
    X,
    Y,
    Replication,
    Dose_percent,
    Omission_description,
    all_of(
      variables_suelo
    ),
    all_of(
      variables_cultivo
    )
  )

for (variable in c(
  variables_suelo,
  variables_cultivo
)) {

  datos_residuales[[paste0(
      "RES_",
      variable
    )]] <- residualizar(
    variable = variable,
    df = datos
  )
}


## 7. Correlaciones globales -----------------------------------

pares <- crossing(
  Variable_suelo = variables_suelo,
  Variable_cultivo = variables_cultivo
)

tabla_correlaciones <- pmap_dfr(
  pares,
  function(
    Variable_suelo,
    Variable_cultivo
  ) {

    x <- datos[[Variable_suelo]]

    y <- datos[[Variable_cultivo]]

    xr <- datos_residuales[[paste0(
        "RES_",
        Variable_suelo
      )]]

    yr <- datos_residuales[[paste0(
        "RES_",
        Variable_cultivo
      )]]

    tibble(
      Variable_suelo = Variable_suelo,
      Variable_cultivo = Variable_cultivo,
      n_original = sum(
        is.finite(
          x
        ) &
          is.finite(
            y
          )
      ),
      Pearson_original = correlacion_segura(
        x,
        y,
        "pearson"
      ),
      Spearman_original = correlacion_segura(
        x,
        y,
        "spearman"
      ),
      n_residual = sum(
        is.finite(
          xr
        ) &
          is.finite(
            yr
          )
      ),
      Pearson_residual = correlacion_segura(
        xr,
        yr,
        "pearson"
      ),
      Spearman_residual = correlacion_segura(
        xr,
        yr,
        "spearman"
      )
    )
  }
) |>
  mutate(
    Cambio_abs_Spearman =
      abs(
        Spearman_residual
      ) -
      abs(
        Spearman_original
      ),
    Estabilidad_direccion = case_when(
      sign(
        Spearman_original
      ) ==
        sign(
          Spearman_residual
        ) ~
        "Dirección estable",

      TRUE ~
        "Dirección cambia tras ajuste"
    ),
    Magnitud_residual = case_when(
      abs(
        Spearman_residual
      ) >= 0.70 ~
        "Alta",

      abs(
        Spearman_residual
      ) >= 0.40 ~
        "Moderada",

      abs(
        Spearman_residual
      ) >= 0.20 ~
        "Limitada",

      TRUE ~
        "Muy limitada"
    )
  )


## 8. Moran bivariado global -----------------------------------

resultados_moran <- list()
simulaciones_moran <- list()

for (
  i in seq_len(
    nrow(
      pares
    )
  )
) {

  variable_suelo <- pares$Variable_suelo[
    i
  ]

  variable_cultivo <- pares$Variable_cultivo[
    i
  ]

  clave <- paste0(
    variable_suelo,
    "__",
    variable_cultivo
  )

  x_original <- datos[[variable_suelo]]

  y_original <- datos[[variable_cultivo]]

  x_residual <- datos_residuales[[paste0(
      "RES_",
      variable_suelo
    )]]

  y_residual <- datos_residuales[[paste0(
      "RES_",
      variable_cultivo
    )]]

  moran_original <- simular_moran_bivariado(
    x = x_original,
    y = y_original,
    listw = pesos,
    nsim = nsim_moran
  )

  moran_residual <- simular_moran_bivariado(
    x = x_residual,
    y = y_residual,
    listw = pesos,
    nsim = nsim_moran
  )

  resultados_moran[[clave]] <- bind_rows(
    tibble(
      Variable_suelo = variable_suelo,
      Variable_cultivo = variable_cultivo,
      Escenario = "Original",
      Moran_bivariado = moran_original$observado,
      Sim_Q025 = moran_original$q025,
      Sim_Mediana = moran_original$mediana,
      Sim_Q975 = moran_original$q975,
      p_Monte_Carlo_informativo = moran_original$p_mc,
      Clasificacion = moran_original$clasificacion
    ),
    tibble(
      Variable_suelo = variable_suelo,
      Variable_cultivo = variable_cultivo,
      Escenario = "Residual_ajustado",
      Moran_bivariado = moran_residual$observado,
      Sim_Q025 = moran_residual$q025,
      Sim_Mediana = moran_residual$mediana,
      Sim_Q975 = moran_residual$q975,
      p_Monte_Carlo_informativo = moran_residual$p_mc,
      Clasificacion = moran_residual$clasificacion
    )
  )

  simulaciones_moran[[paste0(
      clave,
      "__Original"
    )]] <- tibble(
    Variable_suelo = variable_suelo,
    Variable_cultivo = variable_cultivo,
    Escenario = "Original",
    Moran_simulado = moran_original$simulados
  )

  simulaciones_moran[[paste0(
      clave,
      "__Residual"
    )]] <- tibble(
    Variable_suelo = variable_suelo,
    Variable_cultivo = variable_cultivo,
    Escenario = "Residual_ajustado",
    Moran_simulado = moran_residual$simulados
  )
}

tabla_moran <- bind_rows(
  resultados_moran
)

tabla_simulaciones_moran <- bind_rows(
  simulaciones_moran
)


## 9. Asociación por clases de distancia -----------------------
##
## Versión vectorizada:
##   - La tabla de pares espaciales se construye una sola vez.
##   - Se reutiliza para los 56 pares suelo–cultivo.
##   - Evita crear miles de tibbles dentro de bucles anidados.
##   - Conserva la definición simétrica del producto cruzado:
##       0.5 * (z_xi * z_yj + z_xj * z_yi)

matriz_distancias <- as.matrix(
  dist(
    coords
  )
)

indices_pares <- which(
  upper.tri(
    matriz_distancias
  ),
  arr.ind = TRUE
)

tabla_pares_espaciales <- tibble(
  i = indices_pares[, 1],
  j = indices_pares[, 2],
  Distancia_m = matriz_distancias[
    indices_pares
  ]
)

cortes <- unique(
  quantile(
    tabla_pares_espaciales$Distancia_m,
    probs = seq(
      0,
      1,
      length.out = n_clases_distancia + 1
    ),
    na.rm = TRUE,
    names = FALSE
  )
)

if (length(cortes) < 3) {
  stop(
    "No fue posible construir clases de distancia."
  )
}

tabla_pares_espaciales <- tabla_pares_espaciales |>
  mutate(
    Clase_distancia = cut(
      Distancia_m,
      breaks = cortes,
      include.lowest = TRUE,
      right = TRUE
    )
  )

calcular_asociacion_distancia <- function(
  x,
  y,
  variable_suelo,
  variable_cultivo,
  escenario,
  tabla_pares
) {

  zx <- estandarizar(
    x
  )

  zy <- estandarizar(
    y
  )

  if (
    all(
      !is.finite(
        zx
      )
    ) ||
      all(
        !is.finite(
          zy
        )
      )
  ) {
    return(
      tibble()
    )
  }

  productos <- 0.5 *
    (
      zx[
        tabla_pares$i
      ] *
        zy[
          tabla_pares$j
        ] +
      zx[
        tabla_pares$j
      ] *
        zy[
          tabla_pares$i
        ]
    )

  tabla_pares |>
    mutate(
      Producto_cruzado = productos
    ) |>
    filter(
      is.finite(
        Producto_cruzado
      ),
      !is.na(
        Clase_distancia
      )
    ) |>
    group_by(
      Clase_distancia
    ) |>
    summarise(
      n_pares = n(),
      Distancia_media_m = mean(
        Distancia_m,
        na.rm = TRUE
      ),
      Asociacion_cruzada_media = mean(
        Producto_cruzado,
        na.rm = TRUE
      ),
      DE_asociacion = sd(
        Producto_cruzado,
        na.rm = TRUE
      ),
      .groups = "drop"
    ) |>
    mutate(
      Variable_suelo = variable_suelo,
      Variable_cultivo = variable_cultivo,
      Escenario = escenario,
      Clase_distancia = as.character(
        Clase_distancia
      ),
      .before = 1
    )
}

asociacion_distancia <- pmap_dfr(
  pares,
  function(
    Variable_suelo,
    Variable_cultivo
  ) {

    original <- calcular_asociacion_distancia(
      x = datos[[Variable_suelo]],
      y = datos[[Variable_cultivo]],
      variable_suelo = Variable_suelo,
      variable_cultivo = Variable_cultivo,
      escenario = "Original",
      tabla_pares = tabla_pares_espaciales
    )

    residual <- calcular_asociacion_distancia(
      x = datos_residuales[[paste0(
          "RES_",
          Variable_suelo
        )]
      ],
      y = datos_residuales[[paste0(
          "RES_",
          Variable_cultivo
        )]
      ],
      variable_suelo = Variable_suelo,
      variable_cultivo = Variable_cultivo,
      escenario = "Residual_ajustado",
      tabla_pares = tabla_pares_espaciales
    )

    bind_rows(
      original,
      residual
    )
  }
)


## 10. Comparación CE original y corregida ---------------------

tabla_ce <- tabla_correlaciones |>
  filter(
    Variable_suelo %in% c(
      "Suelo_CE",
      "Suelo_CE_corregida_OLS"
    )
  ) |>
  select(
    Variable_suelo,
    Variable_cultivo,
    Pearson_original,
    Spearman_original,
    Pearson_residual,
    Spearman_residual,
    Magnitud_residual
  )

comparacion_ce <- tabla_ce |>
  pivot_longer(
    cols = c(
      Pearson_original,
      Spearman_original,
      Pearson_residual,
      Spearman_residual
    ),
    names_to = "Metrica",
    values_to = "Valor"
  ) |>
  pivot_wider(
    names_from = Variable_suelo,
    values_from = Valor
  ) |>
  mutate(
    Diferencia_corregida_menos_original =
      Suelo_CE_corregida_OLS -
      Suelo_CE
  )


## 11. Síntesis de asociaciones -------------------------------

tabla_sintesis <- tabla_correlaciones |>
  left_join(
    tabla_moran |>
      filter(
        Escenario == "Residual_ajustado"
      ) |>
      select(
        Variable_suelo,
        Variable_cultivo,
        Moran_bivariado_residual =
          Moran_bivariado,
        Moran_Q025 =
          Sim_Q025,
        Moran_Q975 =
          Sim_Q975,
        Clasificacion_Moran =
          Clasificacion
      ),
    by = c(
      "Variable_suelo",
      "Variable_cultivo"
    )
  ) |>
  mutate(
    Coherencia_correlacion_espacial = case_when(
      sign(
        Spearman_residual
      ) ==
        sign(
          Moran_bivariado_residual
        ) ~
        "Dirección coherente",

      TRUE ~
        "Dirección no coincidente"
    ),
    Prioridad_interpretacion = case_when(
      abs(
        Spearman_residual
      ) >= 0.40 &
        str_detect(
          Clasificacion_Moran,
          "mayor que"
        ) ~
        "Alta: asociación residual y espacial concordantes",

      abs(
        Spearman_residual
      ) >= 0.20 ~
        "Media: asociación residual limitada o moderada",

      TRUE ~
        "Baja: evidencia débil"
    )
  ) |>
  arrange(
    desc(
      abs(
        Spearman_residual
      )
    ),
    desc(
      abs(
        Moran_bivariado_residual
      )
    )
  )


## 12. Figuras: mapas de calor ---------------------------------

cor_heatmap <- tabla_correlaciones |>
  dplyr::select(
    Variable_suelo,
    Variable_cultivo,
    Spearman_original,
    Spearman_residual
  ) |>
  tidyr::pivot_longer(
    cols = c(
      Spearman_original,
      Spearman_residual
    ),
    names_to = "Escenario",
    values_to = "Spearman"
  ) |>
  mutate(
    Escenario = case_when(
      Escenario == "Spearman_original" ~ "Original",
      Escenario == "Spearman_residual" ~ "Residual ajustado",
      TRUE ~ as.character(Escenario)
    )
  )

g_correlaciones <- ggplot(
  cor_heatmap,
  aes(
    x = Variable_cultivo,
    y = Variable_suelo,
    fill = Spearman
  )
) +
  geom_tile() +
  geom_text(
    aes(
      label = sprintf(
        "%.2f",
        Spearman
      )
    ),
    size = 2.7
  ) +
  facet_wrap(
    ~ Escenario
  ) +
  labs(
    title =
      "Correlaciones de Spearman entre suelo y cultivo",
    subtitle =
      "Comparación entre valores originales y residuos ajustados",
    x = "Variable de cultivo",
    y = "Variable de suelo",
    fill = "ρ"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

ggsave(
  filename = file.path(
    ruta_figuras,
    "01_heatmap_correlaciones_suelo_cultivo.png"
  ),
  plot = g_correlaciones,
  width = 14,
  height = 9,
  dpi = 400
)


moran_heatmap <- tabla_moran |>
  filter(
    Escenario == "Residual_ajustado"
  )

g_moran <- ggplot(
  moran_heatmap,
  aes(
    x = Variable_cultivo,
    y = Variable_suelo,
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
    size = 2.7
  ) +
  labs(
    title =
      "Moran bivariado residual suelo–cultivo",
    subtitle = paste0(
      "Pesos k-NN con k = ",
      k_real
    ),
    x = "Variable de cultivo",
    y = "Variable de suelo",
    fill = "I bivariado"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

ggsave(
  filename = file.path(
    ruta_figuras,
    "02_heatmap_moran_bivariado_residual.png"
  ),
  plot = g_moran,
  width = 11,
  height = 8,
  dpi = 400
)


## 13. Figuras: asociación por distancia -----------------------

top_pares <- tabla_sintesis |>
  slice_head(
    n = min(
      12,
      nrow(tabla_sintesis)
    )
  ) |>
  transmute(
    Clave = paste(
      Variable_suelo,
      Variable_cultivo,
      sep = "__"
    )
  ) |>
  pull(
    Clave
  )

distancia_fig <- asociacion_distancia |>
  mutate(
    Clave = paste(
      Variable_suelo,
      Variable_cultivo,
      sep = "__"
    )
  ) |>
  filter(
    Clave %in% top_pares
  )

g_distancia <- ggplot(
  distancia_fig,
  aes(
    x = Distancia_media_m,
    y = Asociacion_cruzada_media,
    linetype = Escenario
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dotted"
  ) +
  geom_line() +
  geom_point() +
  facet_wrap(
    Variable_suelo ~ Variable_cultivo,
    scales = "free_y"
  ) +
  labs(
    title =
      "Asociación suelo–cultivo por clases de distancia",
    subtitle =
      "Se muestran los pares con mayor asociación residual",
    x = "Distancia media (m)",
    y = "Producto cruzado estandarizado medio",
    linetype = "Escenario"
  ) +
  theme_minimal()

ggsave(
  filename = file.path(
    ruta_figuras,
    "03_asociacion_cruzada_por_distancia.png"
  ),
  plot = g_distancia,
  width = 16,
  height = 18,
  dpi = 400
)


## 14. Figura CE original vs corregida -------------------------

g_ce <- ggplot(
  tabla_ce |>
    pivot_longer(
      cols = c(
        Spearman_original,
        Spearman_residual
      ),
      names_to = "Escenario",
      values_to = "Spearman"
    ),
  aes(
    x = Variable_cultivo,
    y = Spearman,
    shape = Variable_suelo
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dotted"
  ) +
  geom_point(
    size = 2.8,
    position = position_dodge(
      width = 0.35
    )
  ) +
  facet_wrap(
    ~ Escenario
  ) +
  labs(
    title =
      "Relación del cultivo con CE original y CE corregida",
    x = "Variable de cultivo",
    y = "Spearman",
    shape = "Variable de CE"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

ggsave(
  filename = file.path(
    ruta_figuras,
    "04_comparacion_CE_original_corregida.png"
  ),
  plot = g_ce,
  width = 12,
  height = 7,
  dpi = 400
)


## 15. Figura de síntesis --------------------------------------

panel_sintesis <- (
  g_correlaciones /
    g_moran
) +
  plot_annotation(
    title =
      "Síntesis de asociación espacial suelo–cultivo"
  )

ggsave(
  filename = file.path(
    ruta_figuras,
    "05_panel_sintesis_asociacion_suelo_cultivo.png"
  ),
  plot = panel_sintesis,
  width = 14,
  height = 17,
  dpi = 400
)


## 16. Exportar datos ------------------------------------------

write_csv(
  tabla_correlaciones,
  file.path(
    ruta_datos,
    "Correlaciones_suelo_cultivo.csv"
  )
)

write_csv(
  tabla_moran,
  file.path(
    ruta_datos,
    "Moran_bivariado_suelo_cultivo.csv"
  )
)

write_csv(
  asociacion_distancia,
  file.path(
    ruta_datos,
    "Asociacion_suelo_cultivo_por_distancia.csv"
  )
)

write_csv(
  datos_residuales,
  file.path(
    ruta_datos,
    "Datos_residuales_ajustados_suelo_cultivo.csv"
  )
)


## 17. Exportar libro ------------------------------------------

write_xlsx(
  list(
    Sintesis =
      tabla_sintesis,

    Correlaciones =
      tabla_correlaciones,

    Moran_bivariado =
      tabla_moran,

    Asociacion_distancia =
      asociacion_distancia,

    Comparacion_CE =
      comparacion_ce,

    Vecinos =
      tabla_vecinos
  ),
  path = file.path(
    ruta_tablas,
    "Resultados_asociacion_espacial_suelo_cultivo.xlsx"
  )
)


## 18. Guardar objeto para Script 05-D -------------------------

saveRDS(
  list(
    datos = datos,
    datos_residuales = datos_residuales,
    cultivo_sf = cultivo_sf,
    lote = lote,
    pesos = pesos,
    vecinos = vecinos,
    tabla_correlaciones = tabla_correlaciones,
    tabla_moran = tabla_moran,
    tabla_sintesis = tabla_sintesis,
    asociacion_distancia = asociacion_distancia,
    comparacion_ce = comparacion_ce
  ),
  file = file.path(
    ruta_objetos,
    "objetos_asociacion_espacial_suelo_cultivo.rds"
  )
)


## 19. Registro de sesión --------------------------------------

capture.output(
  sessionInfo(),
  file = file.path(
    ruta_salida,
    "sessionInfo_05_C.txt"
  )
)


## 20. Resumen en consola --------------------------------------

cat(
  "\n============================================\n"
)

cat(
  "SCRIPT 05-C FINALIZADO\n\n"
)

cat(
  "Pares suelo–cultivo evaluados: ",
  nrow(
    pares
  ),
  "\n"
)

cat(
  "Vecinos k-NN utilizados: ",
  k_real,
  "\n\n"
)

cat(
  "Principales asociaciones residuales:\n"
)

print(
  tabla_sintesis |>
    select(
      Variable_suelo,
      Variable_cultivo,
      Spearman_residual,
      Moran_bivariado_residual,
      Prioridad_interpretacion
    ) |>
    slice_head(
      n = 12
    )
)

cat(
  "\nObjeto para Script 05-D:\n",
  file.path(
    ruta_objetos,
    "objetos_asociacion_espacial_suelo_cultivo.rds"
  ),
  "\n"
)

cat(
  "============================================\n"
)
