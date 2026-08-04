
## ============================================================
## 05_PATRON_PUNTOS – SCRIPT 05-D
## HOTSPOTS, COLDSPOTS Y ASOCIACIÓN ESPACIAL LOCAL
## ============================================================
##
## Entrada:
##   05_PATRON_PUNTOS/03_ASOCIACION_SUELO_CULTIVO/objetos_R/
##   objetos_asociacion_espacial_suelo_cultivo.rds
##
## Salida:
##   05_PATRON_PUNTOS/04_HOTSPOTS_LISA/
##
## PREGUNTAS
##   1. ¿Dónde se concentran valores localmente altos y bajos?
##   2. ¿Qué nodos forman clústeres Alto–Alto o Bajo–Bajo?
##   3. ¿Qué nodos se comportan como atípicos espaciales?
##   4. ¿Los patrones locales permanecen después de controlar
##      repetición, dosis y omisión?
##
## MÉTODOS
##   - Moran local (LISA) con pesos k-NN.
##   - Getis–Ord Gi* con vecindad k-NN e inclusión del propio nodo.
##   - Referencia Monte Carlo por permutación.
##   - Escenarios:
##       a) valores originales;
##       b) residuos ajustados por repetición + dosis + omisión.
##
## PRECAUCIONES
##   - Los mapas locales son exploratorios.
##   - No se usa el p-valor como criterio principal.
##   - La clasificación se basa en la posición del estadístico
##     observado frente al intervalo central de permutación.
##   - Las variables de suelo proceden de Kriging y no constituyen
##     mediciones independientes en cada nodo de cultivo.
##   - N, P y K son índices instrumentales nominales.
##   - pH y Hv requieren interpretación especialmente prudente.
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

nsim_local <- 499
k_vecinos <- 4
prob_inferior <- 0.025
prob_superior <- 0.975


## 2. Rutas ----------------------------------------------------

ruta_base <- paste0(
  "D:/UDCA_2017_CLASE/46_ARTICULOS_2025/",
  "ARTICULO_CEBADA_DOCTORADO/AJUSTES_ARTÍCULO"
)

ruta_objeto <- file.path(
  ruta_base,
  "05_PATRON_PUNTOS",
  "03_ASOCIACION_SUELO_CULTIVO",
  "objetos_R",
  "objetos_asociacion_espacial_suelo_cultivo.rds"
)

ruta_salida <- file.path(
  ruta_base,
  "05_PATRON_PUNTOS",
  "04_HOTSPOTS_LISA"
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
    "No existe el objeto del Script 05-C:\n",
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
datos_residuales <- obj$datos_residuales
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

variables_analisis <- c(
  variables_cultivo,
  variables_suelo
)

if (length(variables_analisis) == 0) {
  stop(
    "No se encontraron variables para el análisis local."
  )
}

columnas_requeridas <- c(
  "Nodo",
  "X",
  "Y"
)

faltantes <- setdiff(
  columnas_requeridas,
  names(datos)
)

if (length(faltantes) > 0) {
  stop(
    "Faltan columnas: ",
    paste(
      faltantes,
      collapse = ", "
    )
  )
}


## 4. Construir pesos espaciales -------------------------------

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

pesos_W <- nb2listw(
  vecinos,
  style = "W",
  zero.policy = TRUE
)

vecinos_con_self <- include.self(
  vecinos
)

pesos_G <- nb2listw(
  vecinos_con_self,
  style = "B",
  zero.policy = TRUE
)

W_moran <- listw2mat(
  pesos_W
)

W_gi <- listw2mat(
  pesos_G
)


## 5. Funciones auxiliares -------------------------------------

estandarizar <- function(x) {

  x <- as.numeric(
    x
  )

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


obtener_vector <- function(
  variable,
  escenario
) {

  if (escenario == "Original") {
    return(
      as.numeric(
        datos[[variable]]
      )
    )
  }

  nombre_residual <- paste0(
    "RES_",
    variable
  )

  if (!nombre_residual %in% names(datos_residuales)) {
    return(
      rep(
        NA_real_,
        nrow(
          datos
        )
      )
    )
  }

  as.numeric(
    datos_residuales[[nombre_residual]]
  )
}


estadistico_moran_local <- function(
  x,
  W
) {

  z <- estandarizar(
    x
  )

  rezago <- as.numeric(
    W %*% z
  )

  list(
    z = z,
    rezago = rezago,
    I_local = z * rezago
  )
}


estadistico_gi_local <- function(
  x,
  W
) {

  x <- as.numeric(
    x
  )

  suma_local <- as.numeric(
    W %*% x
  )

  media_global <- mean(
    x,
    na.rm = TRUE
  )

  desviacion_global <- sd(
    x,
    na.rm = TRUE
  )

  suma_pesos <- rowSums(
    W
  )

  suma_pesos_cuadrados <- rowSums(
    W^2
  )

  n <- length(
    x
  )

  denominador <- desviacion_global *
    sqrt(
      (
        n *
          suma_pesos_cuadrados -
          suma_pesos^2
      ) /
        max(
          n - 1,
          1
        )
    )

  gi <- (
    suma_local -
      media_global *
        suma_pesos
  ) /
    denominador

  gi[
    !is.finite(
      gi
    )
  ] <- NA_real_

  gi
}


clasificar_cuadrante <- function(
  z,
  lag_z,
  I_obs,
  q025,
  q975
) {

  extremo <- case_when(
    I_obs > q975 ~ "Positivo_extremo",
    I_obs < q025 ~ "Negativo_extremo",
    TRUE ~ "Dentro_referencia"
  )

  case_when(
    extremo == "Positivo_extremo" &
      z > 0 &
      lag_z > 0 ~
      "Alto-Alto",

    extremo == "Positivo_extremo" &
      z < 0 &
      lag_z < 0 ~
      "Bajo-Bajo",

    extremo == "Negativo_extremo" &
      z > 0 &
      lag_z < 0 ~
      "Alto-Bajo",

    extremo == "Negativo_extremo" &
      z < 0 &
      lag_z > 0 ~
      "Bajo-Alto",

    TRUE ~
      "Sin_clasificar"
  )
}


analizar_variable_local <- function(
  variable,
  escenario
) {

  x <- obtener_vector(
    variable,
    escenario
  )

  validos <- is.finite(
    x
  )

  if (sum(validos) != length(x)) {
    warning(
      "La variable ",
      variable,
      " en escenario ",
      escenario,
      " contiene valores no finitos. Se omite."
    )

    return(
      list(
        lisa = tibble(),
        gi = tibble(),
        simulaciones_lisa = tibble(),
        simulaciones_gi = tibble()
      )
    )
  }

  moran_obs <- estadistico_moran_local(
    x,
    W_moran
  )

  gi_obs <- estadistico_gi_local(
    x,
    W_gi
  )

  I_sim <- matrix(
    NA_real_,
    nrow = length(x),
    ncol = nsim_local
  )

  Gi_sim <- matrix(
    NA_real_,
    nrow = length(x),
    ncol = nsim_local
  )

  for (s in seq_len(nsim_local)) {

    x_perm <- sample(
      x,
      replace = FALSE
    )

    I_sim[, s] <- estadistico_moran_local(
      x_perm,
      W_moran
    )$I_local

    Gi_sim[, s] <- estadistico_gi_local(
      x_perm,
      W_gi
    )
  }

  I_q025 <- apply(
    I_sim,
    1,
    quantile,
    probs = prob_inferior,
    na.rm = TRUE
  )

  I_q975 <- apply(
    I_sim,
    1,
    quantile,
    probs = prob_superior,
    na.rm = TRUE
  )

  Gi_q025 <- apply(
    Gi_sim,
    1,
    quantile,
    probs = prob_inferior,
    na.rm = TRUE
  )

  Gi_q975 <- apply(
    Gi_sim,
    1,
    quantile,
    probs = prob_superior,
    na.rm = TRUE
  )

  p_lisa <- (
    1 +
      rowSums(
        abs(
          sweep(
            I_sim,
            1,
            apply(
              I_sim,
              1,
              median,
              na.rm = TRUE
            ),
            "-"
          )
        ) >=
          abs(
            moran_obs$I_local -
              apply(
                I_sim,
                1,
                median,
                na.rm = TRUE
              )
          ),
        na.rm = TRUE
      )
  ) /
    (
      nsim_local +
        1
    )

  p_gi <- (
    1 +
      rowSums(
        abs(
          sweep(
            Gi_sim,
            1,
            apply(
              Gi_sim,
              1,
              median,
              na.rm = TRUE
            ),
            "-"
          )
        ) >=
          abs(
            gi_obs -
              apply(
                Gi_sim,
                1,
                median,
                na.rm = TRUE
              )
          ),
        na.rm = TRUE
      )
  ) /
    (
      nsim_local +
        1
    )

  cuadrante <- clasificar_cuadrante(
    z = moran_obs$z,
    lag_z = moran_obs$rezago,
    I_obs = moran_obs$I_local,
    q025 = I_q025,
    q975 = I_q975
  )

  clase_gi <- case_when(
    gi_obs > Gi_q975 ~ "Hotspot",
    gi_obs < Gi_q025 ~ "Coldspot",
    TRUE ~ "Sin_clasificar"
  )

  tabla_lisa <- tibble(
    Nodo = datos$Nodo,
    X = datos$X,
    Y = datos$Y,
    Variable = variable,
    Escenario = escenario,
    Valor = x,
    Z = moran_obs$z,
    Rezago_Z = moran_obs$rezago,
    Moran_local = moran_obs$I_local,
    Perm_Q025 = I_q025,
    Perm_Q975 = I_q975,
    p_Monte_Carlo_informativo = p_lisa,
    Clase_LISA = cuadrante
  )

  tabla_gi <- tibble(
    Nodo = datos$Nodo,
    X = datos$X,
    Y = datos$Y,
    Variable = variable,
    Escenario = escenario,
    Valor = x,
    Gi_estrella = gi_obs,
    Perm_Q025 = Gi_q025,
    Perm_Q975 = Gi_q975,
    p_Monte_Carlo_informativo = p_gi,
    Clase_Gi = clase_gi
  )

  simulaciones_lisa <- tibble(
    Variable = variable,
    Escenario = escenario,
    Nodo = rep(
      datos$Nodo,
      times = nsim_local
    ),
    Simulacion = rep(
      seq_len(
        nsim_local
      ),
      each = nrow(
        datos
      )
    ),
    Moran_local_simulado = as.vector(
      I_sim
    )
  )

  simulaciones_gi <- tibble(
    Variable = variable,
    Escenario = escenario,
    Nodo = rep(
      datos$Nodo,
      times = nsim_local
    ),
    Simulacion = rep(
      seq_len(
        nsim_local
      ),
      each = nrow(
        datos
      )
    ),
    Gi_simulado = as.vector(
      Gi_sim
    )
  )

  list(
    lisa = tabla_lisa,
    gi = tabla_gi,
    simulaciones_lisa = simulaciones_lisa,
    simulaciones_gi = simulaciones_gi
  )
}


## 6. Ejecutar análisis local ----------------------------------

escenarios <- c(
  "Original",
  "Residual_ajustado"
)

combinaciones <- crossing(
  Variable = variables_analisis,
  Escenario = escenarios
)

resultados_locales <- pmap(
  combinaciones,
  function(
    Variable,
    Escenario
  ) {

    message(
      "Procesando ",
      Variable,
      " – ",
      Escenario
    )

    analizar_variable_local(
      variable = Variable,
      escenario = Escenario
    )
  }
)

tabla_lisa <- bind_rows(
  map(
    resultados_locales,
    "lisa"
  )
)

tabla_gi <- bind_rows(
  map(
    resultados_locales,
    "gi"
  )
)

simulaciones_lisa <- bind_rows(
  map(
    resultados_locales,
    "simulaciones_lisa"
  )
)

simulaciones_gi <- bind_rows(
  map(
    resultados_locales,
    "simulaciones_gi"
  )
)


## 7. Resúmenes por variable -----------------------------------

resumen_lisa <- tabla_lisa |>
  count(
    Variable,
    Escenario,
    Clase_LISA,
    name = "n_nodos"
  ) |>
  group_by(
    Variable,
    Escenario
  ) |>
  mutate(
    Porcentaje = 100 *
      n_nodos /
      sum(
        n_nodos
      )
  ) |>
  ungroup()

resumen_gi <- tabla_gi |>
  count(
    Variable,
    Escenario,
    Clase_Gi,
    name = "n_nodos"
  ) |>
  group_by(
    Variable,
    Escenario
  ) |>
  mutate(
    Porcentaje = 100 *
      n_nodos /
      sum(
        n_nodos
      )
  ) |>
  ungroup()


## 8. Estabilidad entre escenarios -----------------------------

estabilidad_lisa <- tabla_lisa |>
  select(
    Nodo,
    Variable,
    Escenario,
    Clase_LISA
  ) |>
  pivot_wider(
    names_from = Escenario,
    values_from = Clase_LISA
  ) |>
  mutate(
    Clase_estable =
      Original ==
      Residual_ajustado
  )

estabilidad_gi <- tabla_gi |>
  select(
    Nodo,
    Variable,
    Escenario,
    Clase_Gi
  ) |>
  pivot_wider(
    names_from = Escenario,
    values_from = Clase_Gi
  ) |>
  mutate(
    Clase_estable =
      Original ==
      Residual_ajustado
  )

resumen_estabilidad <- estabilidad_lisa |>
  group_by(
    Variable
  ) |>
  summarise(
    LISA_n_estables = sum(
      Clase_estable,
      na.rm = TRUE
    ),
    LISA_porcentaje_estable = 100 *
      mean(
        Clase_estable,
        na.rm = TRUE
      ),
    .groups = "drop"
  ) |>
  left_join(
    estabilidad_gi |>
      group_by(
        Variable
      ) |>
      summarise(
        Gi_n_estables = sum(
          Clase_estable,
          na.rm = TRUE
        ),
        Gi_porcentaje_estable = 100 *
          mean(
            Clase_estable,
            na.rm = TRUE
          ),
        .groups = "drop"
      ),
    by = "Variable"
  )


## 9. Síntesis interpretativa ----------------------------------

tabla_sintesis <- resumen_estabilidad |>
  left_join(
    resumen_lisa |>
      filter(
        Escenario == "Residual_ajustado",
        Clase_LISA != "Sin_clasificar"
      ) |>
      group_by(
        Variable
      ) |>
      summarise(
        Nodos_LISA_clasificados = sum(
          n_nodos
        ),
        .groups = "drop"
      ),
    by = "Variable"
  ) |>
  left_join(
    resumen_gi |>
      filter(
        Escenario == "Residual_ajustado",
        Clase_Gi != "Sin_clasificar"
      ) |>
      group_by(
        Variable
      ) |>
      summarise(
        Nodos_Gi_clasificados = sum(
          n_nodos
        ),
        .groups = "drop"
      ),
    by = "Variable"
  ) |>
  mutate(
    Nodos_LISA_clasificados = replace_na(
      Nodos_LISA_clasificados,
      0L
    ),
    Nodos_Gi_clasificados = replace_na(
      Nodos_Gi_clasificados,
      0L
    ),
    Evidencia_local_residual = case_when(
      Nodos_LISA_clasificados >= 5 |
        Nodos_Gi_clasificados >= 5 ~
        "Estructura local residual relevante",

      Nodos_LISA_clasificados > 0 |
        Nodos_Gi_clasificados > 0 ~
        "Estructura local residual limitada",

      TRUE ~
        "Sin clasificación local residual"
    ),
    Advertencia = case_when(
      Variable %in% c(
        "Suelo_N",
        "Suelo_P",
        "Suelo_K"
      ) ~
        "Índice instrumental nominal",

      Variable == "Suelo_pH" ~
        "pH exploratorio por validación limitada",

      Variable == "Suelo_Hv" ~
        "Hv con estructura espacial débil",

      Variable == "Suelo_CE_corregida_OLS" ~
        "CE corregida mediante calibración exploratoria",

      str_detect(
        Variable,
        "^Suelo_"
      ) ~
        "Variable interpolada mediante Kriging",

      TRUE ~
        "Variable observada del cultivo"
    )
  ) |>
  arrange(
    desc(
      Nodos_LISA_clasificados +
        Nodos_Gi_clasificados
    )
  )


## 10. Preparar geometrías para mapas --------------------------

lote_sf <- st_as_sf(
  lote
)

lote_sf <- st_transform(
  lote_sf,
  st_crs(
    cultivo_sf
  )
)

lote_borde <- st_boundary(
  st_union(
    lote_sf
  )
)


## 11. Funciones de mapas --------------------------------------

crear_mapa_lisa <- function(
  grupo_variables,
  escenario,
  titulo
) {

  df <- tabla_lisa |>
    filter(
      Variable %in% grupo_variables,
      Escenario == escenario
    )

  ggplot() +
    geom_sf(
      data = lote_borde,
      inherit.aes = FALSE,
      fill = NA,
      linewidth = 0.45
    ) +
    geom_point(
      data = df,
      aes(
        x = X,
        y = Y,
        shape = Clase_LISA
      ),
      size = 2
    ) +
    facet_wrap(
      ~ Variable,
      ncol = 3
    ) +
    coord_sf(
      expand = FALSE
    ) +
    labs(
      title = titulo,
      subtitle = escenario,
      x = "X (m)",
      y = "Y (m)",
      shape = "Clase LISA"
    ) +
    theme_minimal()
}


crear_mapa_gi <- function(
  grupo_variables,
  escenario,
  titulo
) {

  df <- tabla_gi |>
    filter(
      Variable %in% grupo_variables,
      Escenario == escenario
    )

  ggplot() +
    geom_sf(
      data = lote_borde,
      inherit.aes = FALSE,
      fill = NA,
      linewidth = 0.45
    ) +
    geom_point(
      data = df,
      aes(
        x = X,
        y = Y,
        shape = Clase_Gi
      ),
      size = 2
    ) +
    facet_wrap(
      ~ Variable,
      ncol = 3
    ) +
    coord_sf(
      expand = FALSE
    ) +
    labs(
      title = titulo,
      subtitle = escenario,
      x = "X (m)",
      y = "Y (m)",
      shape = "Gi*"
    ) +
    theme_minimal()
}


## 12. Figuras LISA --------------------------------------------

g_lisa_cultivo_original <- crear_mapa_lisa(
  grupo_variables = variables_cultivo,
  escenario = "Original",
  titulo = "Clústeres y atípicos locales del cultivo"
)

g_lisa_cultivo_residual <- crear_mapa_lisa(
  grupo_variables = variables_cultivo,
  escenario = "Residual_ajustado",
  titulo = "LISA del cultivo después de controlar el diseño"
)

g_lisa_suelo_original <- crear_mapa_lisa(
  grupo_variables = variables_suelo,
  escenario = "Original",
  titulo = "Clústeres y atípicos locales del suelo predicho"
)

g_lisa_suelo_residual <- crear_mapa_lisa(
  grupo_variables = variables_suelo,
  escenario = "Residual_ajustado",
  titulo = "LISA residual de las variables de suelo"
)

ggsave(
  filename = file.path(
    ruta_figuras,
    "01_LISA_cultivo_original.png"
  ),
  plot = g_lisa_cultivo_original,
  width = 14,
  height = 11,
  dpi = 400
)

ggsave(
  filename = file.path(
    ruta_figuras,
    "02_LISA_cultivo_residual.png"
  ),
  plot = g_lisa_cultivo_residual,
  width = 14,
  height = 11,
  dpi = 400
)

ggsave(
  filename = file.path(
    ruta_figuras,
    "03_LISA_suelo_original.png"
  ),
  plot = g_lisa_suelo_original,
  width = 14,
  height = 12,
  dpi = 400
)

ggsave(
  filename = file.path(
    ruta_figuras,
    "04_LISA_suelo_residual.png"
  ),
  plot = g_lisa_suelo_residual,
  width = 14,
  height = 12,
  dpi = 400
)


## 13. Figuras Gi* ---------------------------------------------

g_gi_cultivo_original <- crear_mapa_gi(
  grupo_variables = variables_cultivo,
  escenario = "Original",
  titulo = "Hotspots y coldspots del cultivo"
)

g_gi_cultivo_residual <- crear_mapa_gi(
  grupo_variables = variables_cultivo,
  escenario = "Residual_ajustado",
  titulo = "Gi* del cultivo después de controlar el diseño"
)

g_gi_suelo_original <- crear_mapa_gi(
  grupo_variables = variables_suelo,
  escenario = "Original",
  titulo = "Hotspots y coldspots del suelo predicho"
)

g_gi_suelo_residual <- crear_mapa_gi(
  grupo_variables = variables_suelo,
  escenario = "Residual_ajustado",
  titulo = "Gi* residual de las variables de suelo"
)

ggsave(
  filename = file.path(
    ruta_figuras,
    "05_Gi_cultivo_original.png"
  ),
  plot = g_gi_cultivo_original,
  width = 14,
  height = 11,
  dpi = 400
)

ggsave(
  filename = file.path(
    ruta_figuras,
    "06_Gi_cultivo_residual.png"
  ),
  plot = g_gi_cultivo_residual,
  width = 14,
  height = 11,
  dpi = 400
)

ggsave(
  filename = file.path(
    ruta_figuras,
    "07_Gi_suelo_original.png"
  ),
  plot = g_gi_suelo_original,
  width = 14,
  height = 12,
  dpi = 400
)

ggsave(
  filename = file.path(
    ruta_figuras,
    "08_Gi_suelo_residual.png"
  ),
  plot = g_gi_suelo_residual,
  width = 14,
  height = 12,
  dpi = 400
)


## 14. Figura de conteos ---------------------------------------

conteos_figura <- bind_rows(
  resumen_lisa |>
    transmute(
      Variable,
      Escenario,
      Metodo = "LISA",
      Clase = Clase_LISA,
      n_nodos
    ),
  resumen_gi |>
    transmute(
      Variable,
      Escenario,
      Metodo = "Gi*",
      Clase = Clase_Gi,
      n_nodos
    )
) |>
  filter(
    Clase != "Sin_clasificar"
  )

g_conteos <- ggplot(
  conteos_figura,
  aes(
    x = Variable,
    y = n_nodos,
    fill = Clase
  )
) +
  geom_col() +
  facet_grid(
    Metodo ~ Escenario,
    scales = "free_y"
  ) +
  labs(
    title = "Número de nodos con clasificación espacial local",
    x = "Variable",
    y = "Número de nodos",
    fill = "Clase"
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
    "09_conteo_clasificaciones_locales.png"
  ),
  plot = g_conteos,
  width = 15,
  height = 10,
  dpi = 400
)


## 15. Panel de síntesis ---------------------------------------

panel_sintesis <- (
  g_lisa_cultivo_residual /
    g_gi_cultivo_residual
) +
  plot_annotation(
    title = "Síntesis de estructura espacial local residual del cultivo"
  )

ggsave(
  filename = file.path(
    ruta_figuras,
    "10_panel_sintesis_local_cultivo_residual.png"
  ),
  plot = panel_sintesis,
  width = 14,
  height = 21,
  dpi = 400
)


## 16. Exportar datos ------------------------------------------

write_csv(
  tabla_lisa,
  file.path(
    ruta_datos,
    "Resultados_LISA_local.csv"
  )
)

write_csv(
  tabla_gi,
  file.path(
    ruta_datos,
    "Resultados_Getis_Ord_Gi.csv"
  )
)

write_csv(
  estabilidad_lisa,
  file.path(
    ruta_datos,
    "Estabilidad_LISA_original_residual.csv"
  )
)

write_csv(
  estabilidad_gi,
  file.path(
    ruta_datos,
    "Estabilidad_Gi_original_residual.csv"
  )
)


## 17. Exportar libro ------------------------------------------

write_xlsx(
  list(
    Sintesis =
      tabla_sintesis,

    Resumen_LISA =
      resumen_lisa,

    Resumen_Gi =
      resumen_gi,

    Estabilidad =
      resumen_estabilidad,

    LISA_nodos =
      tabla_lisa,

    Gi_nodos =
      tabla_gi,

    Estabilidad_LISA =
      estabilidad_lisa,

    Estabilidad_Gi =
      estabilidad_gi
  ),
  path = file.path(
    ruta_tablas,
    "Resultados_hotspots_LISA_Gi.xlsx"
  )
)


## 18. Guardar objeto para Script 05-E -------------------------

saveRDS(
  list(
    datos = datos,
    datos_residuales = datos_residuales,
    cultivo_sf = cultivo_sf,
    lote = lote,
    vecinos = vecinos,
    pesos_W = pesos_W,
    pesos_G = pesos_G,
    tabla_lisa = tabla_lisa,
    tabla_gi = tabla_gi,
    resumen_lisa = resumen_lisa,
    resumen_gi = resumen_gi,
    resumen_estabilidad = resumen_estabilidad,
    tabla_sintesis = tabla_sintesis
  ),
  file = file.path(
    ruta_objetos,
    "objetos_hotspots_LISA_Gi.rds"
  )
)


## 19. Registro de sesión --------------------------------------

capture.output(
  sessionInfo(),
  file = file.path(
    ruta_salida,
    "sessionInfo_05_D.txt"
  )
)


## 20. Resumen en consola --------------------------------------

cat(
  "\n============================================\n"
)

cat(
  "SCRIPT 05-D FINALIZADO\n\n"
)

cat(
  "Variables analizadas: ",
  length(
    variables_analisis
  ),
  "\n"
)

cat(
  "Escenarios: Original y Residual_ajustado\n"
)

cat(
  "Permutaciones por variable y escenario: ",
  nsim_local,
  "\n\n"
)

cat(
  "Síntesis de estructura local residual:\n"
)

print(
  tabla_sintesis |>
    select(
      Variable,
      Nodos_LISA_clasificados,
      Nodos_Gi_clasificados,
      LISA_porcentaje_estable,
      Gi_porcentaje_estable,
      Evidencia_local_residual
    )
)

cat(
  "\nObjeto para Script 05-E:\n",
  file.path(
    ruta_objetos,
    "objetos_hotspots_LISA_Gi.rds"
  ),
  "\n"
)

cat(
  "============================================\n"
)
