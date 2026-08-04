
## ============================================================
## 05_PATRON_PUNTOS – SCRIPT 05-A
## DIAGNÓSTICO GLOBAL DEL PATRÓN DE PUNTOS
## ============================================================
##
## Entrada principal:
##   04_TRATAMIENTOS/02_INTEGRACION_SUELO_CULTIVO/objetos_R/
##   objetos_integracion_suelo_cultivo_tratamientos.rds
##
## Salida:
##   05_PATRON_PUNTOS/01_DIAGNOSTICO_GLOBAL/
##
## PREGUNTA
##   ¿La disposición de los 56 nodos de cultivo es compatible con
##   un patrón regular, aleatorio o agregado dentro del lote?
##
## ALCANCE
##   - Este script caracteriza la geometría de los NODOS MUESTREADOS.
##   - No interpreta cada nodo como una planta individual.
##   - La regularidad puede reflejar el diseño de muestreo o el arreglo
##     experimental, no un proceso biológico espontáneo.
##   - Las envolventes CSR se usan como referencia exploratoria.
##   - La interpretación no se centra en p-valores.
##
## ANÁLISIS
##   1. Auditoría geométrica del patrón.
##   2. Distancias al vecino más cercano.
##   3. Índice de Clark–Evans con referencia Monte Carlo.
##   4. Densidad Kernel.
##   5. Funciones G, F y J.
##   6. Funciones K de Ripley y L de Besag.
##   7. Envolventes CSR mediante simulaciones.
##   8. Tabla de síntesis e interpretación automática prudente.
## ============================================================


## 0. Paquetes -------------------------------------------------

paquetes <- c(
  "sf",
  "dplyr",
  "tidyr",
  "purrr",
  "ggplot2",
  "readr",
  "writexl",
  "patchwork",
  "spatstat.geom",
  "spatstat.explore"
)

instalar <- paquetes[
  !paquetes %in% rownames(installed.packages())
]

if (length(instalar) > 0) {
  install.packages(instalar)
}

library(sf)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(readr)
library(writexl)
library(patchwork)
library(spatstat.geom)
library(spatstat.explore)


## 1. Parámetros reproducibles ---------------------------------

set.seed(2026)

nsim_envolventes <- 199
nsim_clark_evans <- 999
correccion_borde <- "border"


## 2. Rutas ----------------------------------------------------

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
  "05_PATRON_PUNTOS",
  "01_DIAGNOSTICO_GLOBAL"
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
    "No existe el objeto del Script 04-02:\n",
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


## 3. Cargar datos ---------------------------------------------

obj <- readRDS(
  ruta_objeto
)

cultivo_sf <- obj$cultivo_sf
lote <- obj$lote
base_integrada <- obj$base_integrada

if (is.null(cultivo_sf) || is.null(lote)) {
  stop(
    "El objeto no contiene cultivo_sf o lote."
  )
}

if (is.na(st_crs(cultivo_sf))) {
  stop(
    "Los nodos de cultivo no tienen CRS."
  )
}

lote <- st_as_sf(
  lote
)

lote <- st_transform(
  lote,
  st_crs(cultivo_sf)
)

lote <- st_make_valid(
  lote
)

lote_union <- st_union(
  lote
)


## 4. Auditoría geométrica -------------------------------------

coordenadas <- st_coordinates(
  cultivo_sf
)

if (nrow(coordenadas) != nrow(cultivo_sf)) {
  stop(
    "La geometría de cultivo no corresponde a puntos simples."
  )
}

duplicados_xy <- duplicated(
  data.frame(
    X = coordenadas[, 1],
    Y = coordenadas[, 2]
  )
)

dentro_lote <- lengths(
  st_intersects(
    cultivo_sf,
    lote_union
  )
) > 0

tabla_auditoria <- tibble(
  Indicador = c(
    "Número de nodos",
    "Coordenadas duplicadas",
    "Nodos dentro o sobre el lote",
    "Nodos fuera del lote",
    "Área de la ventana (m²)",
    "CRS",
    "Simulaciones Clark-Evans",
    "Simulaciones envolventes CSR"
  ),
  Valor = c(
    nrow(cultivo_sf),
    sum(duplicados_xy),
    sum(dentro_lote),
    sum(!dentro_lote),
    round(
      as.numeric(
        st_area(lote_union)
      ),
      4
    ),
    st_crs(cultivo_sf)$input,
    nsim_clark_evans,
    nsim_envolventes
  )
)

if (any(duplicados_xy)) {
  stop(
    "Existen coordenadas duplicadas. Revise los nodos antes de continuar."
  )
}

if (any(!dentro_lote)) {
  stop(
    "Existen nodos fuera del lote. Corrija la coherencia espacial."
  )
}


## 5. Convertir sf a patrón spatstat ---------------------------

ventana <- as.owin(
  lote_union
)

patron <- ppp(
  x = coordenadas[, 1],
  y = coordenadas[, 2],
  window = ventana,
  checkdup = TRUE
)

if (patron$n != nrow(cultivo_sf)) {
  stop(
    "El número de puntos del objeto ppp no coincide con los nodos."
  )
}

area_ventana <- area.owin(
  ventana
)

intensidad_global <- patron$n /
  area_ventana


## 6. Distancia al vecino más cercano --------------------------

dist_nn <- nndist(
  patron
)

tabla_vecino <- base_integrada |>
  select(
    Nodo,
    Experimental_Unit,
    Treatment_ID,
    Dose_percent,
    Replication,
    X,
    Y
  ) |>
  mutate(
    Distancia_vecino_mas_cercano_m = dist_nn
  )

resumen_vecino <- tibble(
  Indicador = c(
    "Distancia mínima (m)",
    "Distancia media (m)",
    "Distancia mediana (m)",
    "Desviación estándar (m)",
    "Distancia máxima (m)",
    "Intensidad global (nodos/m²)"
  ),
  Valor = c(
    min(dist_nn),
    mean(dist_nn),
    median(dist_nn),
    sd(dist_nn),
    max(dist_nn),
    intensidad_global
  )
)


## 7. Clark–Evans con referencia Monte Carlo -------------------

media_nn_observada <- mean(
  dist_nn
)

media_nn_csr_teorica <- 1 /
  (
    2 *
      sqrt(
        intensidad_global
      )
  )

R_clark_evans <- media_nn_observada /
  media_nn_csr_teorica

simular_media_nn <- function(i) {

  sim <- runifpoint(
    n = patron$n,
    win = ventana
  )

  mean(
    nndist(sim)
  )
}

medias_nn_simuladas <- map_dbl(
  seq_len(nsim_clark_evans),
  simular_media_nn
)

q_nn <- quantile(
  medias_nn_simuladas,
  probs = c(
    0.025,
    0.50,
    0.975
  ),
  na.rm = TRUE,
  names = FALSE
)

p_monte_carlo_bilateral <- (
  1 +
    sum(
      abs(
        medias_nn_simuladas -
          median(
            medias_nn_simuladas
          )
      ) >=
        abs(
          media_nn_observada -
            median(
              medias_nn_simuladas
            )
        )
    )
) /
  (
    nsim_clark_evans +
      1
  )

clasificacion_clark <- case_when(
  media_nn_observada >
    q_nn[3] ~
    "Más regular que la referencia CSR",

  media_nn_observada <
    q_nn[1] ~
    "Más agregado que la referencia CSR",

  TRUE ~
    "Compatible con la referencia CSR"
)

tabla_clark_evans <- tibble(
  n = patron$n,
  Area_m2 = area_ventana,
  Intensidad_nodos_m2 = intensidad_global,
  Media_NN_observada_m = media_nn_observada,
  Media_NN_CSR_teorica_m = media_nn_csr_teorica,
  R_Clark_Evans = R_clark_evans,
  Sim_Q025_m = q_nn[1],
  Sim_Mediana_m = q_nn[2],
  Sim_Q975_m = q_nn[3],
  p_Monte_Carlo_informativo = p_monte_carlo_bilateral,
  Clasificacion = clasificacion_clark
)

tabla_simulaciones_clark <- tibble(
  Simulacion = seq_len(
    nsim_clark_evans
  ),
  Media_NN_CSR_m = medias_nn_simuladas
)


## 8. Ancho de banda y densidad Kernel -------------------------

sigma_diggle <- try(
  bw.diggle(
    patron
  ),
  silent = TRUE
)

if (
  inherits(
    sigma_diggle,
    "try-error"
  ) ||
    !is.finite(
      sigma_diggle
    ) ||
    sigma_diggle <= 0
) {

  sigma_ppl <- try(
    bw.ppl(
      patron
    ),
    silent = TRUE
  )

  if (
    inherits(
      sigma_ppl,
      "try-error"
    ) ||
      !is.finite(
        sigma_ppl
      ) ||
      sigma_ppl <= 0
  ) {

    sigma_kernel <- 1.5 *
      median(
        dist_nn
      )

    metodo_sigma <- "1.5 × mediana de vecino más cercano"

  } else {

    sigma_kernel <- as.numeric(
      sigma_ppl
    )

    metodo_sigma <- "bw.ppl"
  }

} else {

  sigma_kernel <- as.numeric(
    sigma_diggle
  )

  metodo_sigma <- "bw.diggle"
}

densidad_kernel <- density.ppp(
  patron,
  sigma = sigma_kernel,
  edge = TRUE,
  at = "pixels",
  eps = min(
    diff(
      ventana$xrange
    ),
    diff(
      ventana$yrange
    )
  ) /
    300
)

kernel_df <- as.data.frame(
  densidad_kernel
)

names(kernel_df) <- c(
  "X",
  "Y",
  "Densidad_nodos_m2"
)

tabla_kernel <- tibble(
  Metodo_sigma = metodo_sigma,
  Sigma_m = sigma_kernel,
  Densidad_min = min(
    kernel_df$Densidad_nodos_m2,
    na.rm = TRUE
  ),
  Densidad_media = mean(
    kernel_df$Densidad_nodos_m2,
    na.rm = TRUE
  ),
  Densidad_max = max(
    kernel_df$Densidad_nodos_m2,
    na.rm = TRUE
  )
)


## 9. Funciones espaciales y envolventes CSR -------------------

# Cada función de resumen espacial admite correcciones de borde y
# rangos de distancia diferentes. En particular, G, F y J no deben
# forzarse con la corrección "border" usada para K y L. Además, una
# ventana estrecha y rotada puede producir valores NA si se impone un
# rango común demasiado amplio. Por ello, spatstat selecciona aquí el
# rango válido de r para cada función y se emplea la corrección
# apropiada en cada caso.

simulador_csr <- expression(
  runifpoint(
    n = patron$n,
    win = ventana
  )
)

crear_envolvente <- function(
  fun,
  correccion,
  ...
) {

  resultado <- try(
    envelope(
      Y = patron,
      fun = fun,
      nsim = nsim_envolventes,
      simulate = simulador_csr,
      correction = correccion,
      savefuns = TRUE,
      savepatterns = FALSE,
      verbose = FALSE,
      ...
    ),
    silent = TRUE
  )

  if (inherits(resultado, "try-error")) {
    stop(
      "No fue posible calcular la envolvente para ",
      deparse(substitute(fun)),
      ". Detalle: ",
      as.character(resultado)
    )
  }

  resultado
}

# G, F y J: estimadores Kaplan-Meier, adecuados para funciones de
# distancia y compatibles con ventanas irregulares.
envolvente_G <- crear_envolvente(
  fun = Gest,
  correccion = "km"
)

envolvente_F <- crear_envolvente(
  fun = Fest,
  correccion = "km"
)

envolvente_J <- crear_envolvente(
  fun = Jest,
  correccion = "km"
)

# K y L: corrección traslacional, estable para ventanas poligonales y
# sin la erosión extrema que puede ocasionar la corrección border.
envolvente_K <- crear_envolvente(
  fun = Kest,
  correccion = "translation"
)

envolvente_L <- crear_envolvente(
  fun = Lest,
  correccion = "translation"
)


## 10. Convertir envolventes a tablas --------------------------

envolvente_a_tabla <- function(
  env,
  nombre_funcion
) {

  df <- as.data.frame(
    env
  )

  columnas_requeridas <- c(
    "r",
    "obs",
    "lo",
    "hi"
  )

  faltantes <- setdiff(
    columnas_requeridas,
    names(df)
  )

  if (length(faltantes) > 0) {
    stop(
      "La envolvente ",
      nombre_funcion,
      " no contiene: ",
      paste(
        faltantes,
        collapse = ", "
      )
    )
  }

  df |>
    transmute(
      Funcion = nombre_funcion,
      r = r,
      Observado = obs,
      Limite_inferior = lo,
      Limite_superior = hi,
      Posicion = case_when(
        Observado >
          Limite_superior ~
          "Sobre_envolvente",

        Observado <
          Limite_inferior ~
          "Bajo_envolvente",

        TRUE ~
          "Dentro_envolvente"
      )
    )
}

tabla_G <- envolvente_a_tabla(
  envolvente_G,
  "G"
)

tabla_F <- envolvente_a_tabla(
  envolvente_F,
  "F"
)

tabla_J <- envolvente_a_tabla(
  envolvente_J,
  "J"
)

tabla_K <- envolvente_a_tabla(
  envolvente_K,
  "K"
)

tabla_L <- envolvente_a_tabla(
  envolvente_L,
  "L"
)

tabla_envolventes <- bind_rows(
  tabla_G,
  tabla_F,
  tabla_J,
  tabla_K,
  tabla_L
)


## 11. Síntesis de desviaciones de las envolventes ------------

resumir_envolvente <- function(
  tabla,
  funcion
) {

  valida <- tabla |>
    filter(
      is.finite(
        Observado
      ),
      is.finite(
        Limite_inferior
      ),
      is.finite(
        Limite_superior
      )
    )

  n_total <- nrow(
    valida
  )

  if (n_total == 0) {
    return(
      tibble(
        Funcion = funcion,
        Porcentaje_dentro = NA_real_,
        Porcentaje_sobre = NA_real_,
        Porcentaje_bajo = NA_real_,
        r_primera_salida_m = NA_real_,
        Direccion_principal = "No interpretable"
      )
    )
  }

  salidas <- valida |>
    filter(
      Posicion != "Dentro_envolvente"
    )

  direccion <- case_when(
    mean(
      valida$Posicion ==
        "Sobre_envolvente"
    ) >
      mean(
        valida$Posicion ==
          "Bajo_envolvente"
      ) ~
      "Predominio sobre la envolvente",

    mean(
      valida$Posicion ==
        "Bajo_envolvente"
    ) >
      mean(
        valida$Posicion ==
          "Sobre_envolvente"
      ) ~
      "Predominio bajo la envolvente",

    TRUE ~
      "Sin dirección dominante"
  )

  tibble(
    Funcion = funcion,
    Porcentaje_dentro = 100 *
      mean(
        valida$Posicion ==
          "Dentro_envolvente"
      ),
    Porcentaje_sobre = 100 *
      mean(
        valida$Posicion ==
          "Sobre_envolvente"
      ),
    Porcentaje_bajo = 100 *
      mean(
        valida$Posicion ==
          "Bajo_envolvente"
      ),
    r_primera_salida_m = ifelse(
      nrow(salidas) > 0,
      min(
        salidas$r,
        na.rm = TRUE
      ),
      NA_real_
    ),
    Direccion_principal = direccion
  )
}

tabla_resumen_envolventes <- map_dfr(
  c(
    "G",
    "F",
    "J",
    "K",
    "L"
  ),
  function(fun) {

    resumir_envolvente(
      tabla_envolventes |>
        filter(
          Funcion == fun
        ),
      fun
    )
  }
)


## 12. Interpretación combinada prudente -----------------------

interpretacion_KL <- tabla_resumen_envolventes |>
  filter(
    Funcion %in% c(
      "K",
      "L"
    )
  )

predominio_KL <- case_when(
  mean(
    interpretacion_KL$Porcentaje_sobre,
    na.rm = TRUE
  ) >
    mean(
      interpretacion_KL$Porcentaje_bajo,
      na.rm = TRUE
    ) ~
    "Tendencia de agregación a algunas escalas",

  mean(
    interpretacion_KL$Porcentaje_bajo,
    na.rm = TRUE
  ) >
    mean(
      interpretacion_KL$Porcentaje_sobre,
      na.rm = TRUE
    ) ~
    "Tendencia de regularidad a algunas escalas",

  TRUE ~
    "Sin tendencia multiescala dominante"
)

interpretacion_global <- tibble(
  Componente = c(
    "Clark-Evans",
    "K y L",
    "Conclusión global",
    "Advertencia metodológica"
  ),
  Interpretacion = c(
    clasificacion_clark,
    predominio_KL,
    paste0(
      "El patrón debe interpretarse como la geometría de los ",
      "56 nodos de muestreo. La conclusión final debe combinar ",
      "Clark-Evans, G/F/J y K/L, evitando atribuir la regularidad ",
      "del diseño a un proceso biológico espontáneo."
    ),
    paste0(
      "Los puntos son unidades observadas del experimento y no ",
      "la localización exhaustiva de plantas individuales."
    )
  )
)


## 13. Figuras -------------------------------------------------

lote_borde <- st_boundary(
  lote_union
)

g_patron <- ggplot() +
  geom_sf(
    data = lote_borde,
    fill = NA,
    linewidth = 0.65
  ) +
  geom_sf(
    data = cultivo_sf,
    size = 1.8
  ) +
  coord_sf(
    expand = FALSE
  ) +
  labs(
    title = "Patrón espacial de los 56 nodos de cultivo",
    subtitle = paste0(
      "Intensidad global = ",
      round(
        intensidad_global,
        4
      ),
      " nodos/m²"
    ),
    x = "X (m)",
    y = "Y (m)"
  ) +
  theme_minimal()

ggsave(
  filename = file.path(
    ruta_figuras,
    "01_patron_global_nodos.png"
  ),
  plot = g_patron,
  width = 8,
  height = 7,
  dpi = 400
)


g_kernel <- ggplot() +
  geom_raster(
    data = kernel_df,
    aes(
      x = X,
      y = Y,
      fill = Densidad_nodos_m2
    )
  ) +
  geom_sf(
    data = lote_borde,
    fill = NA,
    linewidth = 0.55
  ) +
  geom_sf(
    data = cultivo_sf,
    size = 0.8
  ) +
  coord_sf(
    expand = FALSE
  ) +
  labs(
    title = "Densidad Kernel de los nodos de cultivo",
    subtitle = paste0(
      metodo_sigma,
      "; sigma = ",
      round(
        sigma_kernel,
        3
      ),
      " m"
    ),
    x = "X (m)",
    y = "Y (m)",
    fill = "Nodos/m²"
  ) +
  theme_minimal()

ggsave(
  filename = file.path(
    ruta_figuras,
    "02_densidad_kernel_nodos.png"
  ),
  plot = g_kernel,
  width = 8,
  height = 7,
  dpi = 400
)


g_vecino <- ggplot(
  tabla_vecino,
  aes(
    x = Distancia_vecino_mas_cercano_m
  )
) +
  geom_histogram(
    bins = 12,
    color = "white"
  ) +
  geom_vline(
    xintercept = media_nn_observada,
    linetype = "solid"
  ) +
  geom_vline(
    xintercept = q_nn[c(1, 3)],
    linetype = "dashed"
  ) +
  labs(
    title = "Distancia al vecino más cercano",
    subtitle = paste0(
      "Media observada = ",
      round(
        media_nn_observada,
        3
      ),
      " m; intervalo CSR simulado = ",
      round(
        q_nn[1],
        3
      ),
      "–",
      round(
        q_nn[3],
        3
      ),
      " m"
    ),
    x = "Distancia (m)",
    y = "Frecuencia"
  ) +
  theme_minimal()

ggsave(
  filename = file.path(
    ruta_figuras,
    "03_distancia_vecino_mas_cercano.png"
  ),
  plot = g_vecino,
  width = 8,
  height = 6,
  dpi = 400
)


g_clark <- ggplot(
  tabla_simulaciones_clark,
  aes(
    x = Media_NN_CSR_m
  )
) +
  geom_histogram(
    bins = 30,
    color = "white"
  ) +
  geom_vline(
    xintercept = media_nn_observada,
    linetype = "solid",
    linewidth = 0.8
  ) +
  geom_vline(
    xintercept = q_nn[c(1, 3)],
    linetype = "dashed"
  ) +
  labs(
    title = "Clark–Evans con referencia Monte Carlo",
    subtitle = clasificacion_clark,
    x = "Media de distancia al vecino más cercano (m)",
    y = "Frecuencia simulada"
  ) +
  theme_minimal()

ggsave(
  filename = file.path(
    ruta_figuras,
    "04_clark_evans_monte_carlo.png"
  ),
  plot = g_clark,
  width = 8,
  height = 6,
  dpi = 400
)


crear_grafico_envolvente <- function(
  tabla,
  funcion
) {

  ggplot(
    tabla,
    aes(
      x = r
    )
  ) +
    geom_ribbon(
      aes(
        ymin = Limite_inferior,
        ymax = Limite_superior
      ),
      alpha = 0.35
    ) +
    geom_line(
      aes(
        y = Observado
      ),
      linewidth = 0.65
    ) +
    labs(
      title = paste0(
        "Función ",
        funcion,
        " con envolvente CSR"
      ),
      x = "r (m)",
      y = funcion
    ) +
    theme_minimal()
}

g_G <- crear_grafico_envolvente(
  tabla_G,
  "G"
)

g_F <- crear_grafico_envolvente(
  tabla_F,
  "F"
)

g_J <- crear_grafico_envolvente(
  tabla_J,
  "J"
)

g_K <- crear_grafico_envolvente(
  tabla_K,
  "K"
)

g_L <- crear_grafico_envolvente(
  tabla_L,
  "L"
)

panel_GFJ <- (
  g_G +
    g_F +
    g_J
) +
  plot_annotation(
    title = "Funciones de distancia G, F y J"
  )

ggsave(
  filename = file.path(
    ruta_figuras,
    "05_panel_funciones_G_F_J.png"
  ),
  plot = panel_GFJ,
  width = 15,
  height = 5.5,
  dpi = 400
)

panel_KL <- (
  g_K +
    g_L
) +
  plot_annotation(
    title = "Funciones K de Ripley y L de Besag"
  )

ggsave(
  filename = file.path(
    ruta_figuras,
    "06_panel_funciones_K_L.png"
  ),
  plot = panel_KL,
  width = 11,
  height = 5.5,
  dpi = 400
)


panel_sintesis <- (
  g_patron +
    g_kernel
) /
  (
    g_clark +
      g_L
  ) +
  plot_annotation(
    title = "Diagnóstico global del patrón de puntos"
  )

ggsave(
  filename = file.path(
    ruta_figuras,
    "07_panel_sintesis_patron_global.png"
  ),
  plot = panel_sintesis,
  width = 13,
  height = 12,
  dpi = 400
)


## 14. Exportar datos ------------------------------------------

write_csv(
  tabla_vecino,
  file.path(
    ruta_datos,
    "Distancias_vecino_mas_cercano.csv"
  )
)

write_csv(
  kernel_df,
  file.path(
    ruta_datos,
    "Densidad_kernel_nodos.csv"
  )
)

write_csv(
  tabla_envolventes,
  file.path(
    ruta_datos,
    "Funciones_espaciales_envolventes_CSR.csv"
  )
)

write_csv(
  tabla_simulaciones_clark,
  file.path(
    ruta_datos,
    "Simulaciones_Clark_Evans.csv"
  )
)


## 15. Exportar libro de resultados ----------------------------

write_xlsx(
  list(
    Auditoria = tabla_auditoria,
    Resumen_vecino = resumen_vecino,
    Clark_Evans = tabla_clark_evans,
    Kernel = tabla_kernel,
    Resumen_envolventes = tabla_resumen_envolventes,
    Interpretacion_global = interpretacion_global,
    Distancias_nodos = tabla_vecino,
    Funciones_envolventes = tabla_envolventes
  ),
  path = file.path(
    ruta_tablas,
    "Resultados_diagnostico_global_patron_puntos.xlsx"
  )
)


## 16. Guardar objeto para Script 05-B -------------------------

saveRDS(
  list(
    patron = patron,
    ventana = ventana,
    cultivo_sf = cultivo_sf,
    lote = lote_union,
    base_integrada = base_integrada,
    sigma_kernel = sigma_kernel,
    densidad_kernel = densidad_kernel,
    tabla_clark_evans = tabla_clark_evans,
    tabla_envolventes = tabla_envolventes,
    tabla_resumen_envolventes = tabla_resumen_envolventes,
    interpretacion_global = interpretacion_global
  ),
  file = file.path(
    ruta_objetos,
    "objetos_diagnostico_global_patron_puntos.rds"
  )
)


## 17. Registro de sesión --------------------------------------

capture.output(
  sessionInfo(),
  file = file.path(
    ruta_salida,
    "sessionInfo_05_A.txt"
  )
)


## 18. Resumen en consola --------------------------------------

cat(
  "\n============================================\n"
)

cat(
  "SCRIPT 05-A FINALIZADO\n\n"
)

cat(
  "Número de nodos: ",
  patron$n,
  "\n"
)

cat(
  "Área de la ventana: ",
  round(
    area_ventana,
    3
  ),
  " m²\n"
)

cat(
  "Clark–Evans R: ",
  round(
    R_clark_evans,
    4
  ),
  "\n"
)

cat(
  "Clasificación Clark–Evans: ",
  clasificacion_clark,
  "\n"
)

cat(
  "Síntesis K/L: ",
  predominio_KL,
  "\n\n"
)

cat(
  "Objeto para Script 05-B:\n",
  file.path(
    ruta_objetos,
    "objetos_diagnostico_global_patron_puntos.rds"
  ),
  "\n"
)

cat(
  "============================================\n"
)
