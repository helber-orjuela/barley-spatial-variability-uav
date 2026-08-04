
## ============================================================
## 05_PATRON_PUNTOS – SCRIPT 05-B
## PATRONES ESPACIALES SEGÚN DOSIS, REPETICIÓN Y TRATAMIENTO
## ============================================================
##
## Entrada:
##   05_PATRON_PUNTOS/01_DIAGNOSTICO_GLOBAL/objetos_R/
##   objetos_diagnostico_global_patron_puntos.rds
##
## Salida:
##   05_PATRON_PUNTOS/02_PATRONES_POR_GRUPO/
##
## CONTEXTO DEL SCRIPT 05-A
##   El patrón global presentó:
##     - 56 nodos dentro del lote.
##     - Clark–Evans R > 1.
##     - Mayor regularidad que la referencia CSR.
##     - Tendencia de regularidad en K y L a algunas escalas.
##
## PREGUNTAS DEL SCRIPT 05-B
##   1. ¿La regularidad global se conserva dentro de cada dosis?
##   2. ¿Las repeticiones presentan configuraciones espaciales distintas?
##   3. ¿Qué tratamientos tienen cobertura espacial suficiente para
##      una descripción geométrica?
##
## DECISIONES METODOLÓGICAS
##   - Análisis multiescala formal: dosis y repeticiones con n >= 8.
##   - Tratamientos individuales: descripción geométrica únicamente,
##     porque normalmente contienen 2–3 nodos y no permiten estimar
##     de manera estable funciones K/L ni envolventes CSR.
##   - Se utiliza la MISMA ventana del lote para todos los grupos.
##   - La densidad Kernel usa el mismo sigma global del Script 05-A,
##     permitiendo comparación visual entre dosis.
##   - Las simulaciones mantienen el número de puntos de cada grupo.
##   - Los p-valores Monte Carlo se exportan solo como información
##     complementaria; la interpretación prioriza magnitud y gráficos.
## ============================================================


## 0. Paquetes -------------------------------------------------

paquetes <- c(
  "sf",
  "dplyr",
  "tidyr",
  "purrr",
  "stringr",
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
library(stringr)
library(ggplot2)
library(readr)
library(writexl)
library(patchwork)
library(spatstat.geom)
library(spatstat.explore)


## 1. Parámetros reproducibles ---------------------------------

set.seed(2026)

nsim_clark <- 499
nsim_envolventes <- 99
n_min_multiescala <- 8
n_min_geometria <- 2


## 2. Rutas ----------------------------------------------------

ruta_base <- paste0(
  "D:/UDCA_2017_CLASE/46_ARTICULOS_2025/",
  "ARTICULO_CEBADA_DOCTORADO/AJUSTES_ARTÍCULO"
)

ruta_objeto <- file.path(
  ruta_base,
  "05_PATRON_PUNTOS",
  "01_DIAGNOSTICO_GLOBAL",
  "objetos_R",
  "objetos_diagnostico_global_patron_puntos.rds"
)

ruta_salida <- file.path(
  ruta_base,
  "05_PATRON_PUNTOS",
  "02_PATRONES_POR_GRUPO"
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
    "No existe el objeto del Script 05-A:\n",
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


## 3. Cargar objeto --------------------------------------------

obj <- readRDS(
  ruta_objeto
)

patron_global <- obj$patron
ventana <- obj$ventana
cultivo_sf <- obj$cultivo_sf
lote <- obj$lote
base_integrada <- obj$base_integrada
sigma_global <- obj$sigma_kernel
tabla_global_clark <- obj$tabla_clark_evans

columnas_requeridas <- c(
  "Nodo",
  "X",
  "Y",
  "Dose_percent",
  "Replication",
  "Treatment_ID",
  "Treatment_base",
  "Omission_description"
)

faltantes <- setdiff(
  columnas_requeridas,
  names(base_integrada)
)

if (length(faltantes) > 0) {
  stop(
    "Faltan columnas en base_integrada: ",
    paste(
      faltantes,
      collapse = ", "
    )
  )
}

datos <- base_integrada |>
  mutate(
    Grupo_dosis = case_when(
      as.character(Dose_percent) %in% c(
        "100",
        "D100"
      ) ~ "D100",

      as.character(Dose_percent) %in% c(
        "75",
        "D75"
      ) ~ "D75",

      TRUE ~ as.character(Dose_percent)
    ),
    Grupo_repeticion = paste0(
      "R",
      as.character(Replication)
    ),
    Grupo_tratamiento = as.character(
      Treatment_ID
    )
  )


## 4. Funciones auxiliares -------------------------------------

crear_ppp_grupo <- function(df) {

  ppp(
    x = df$X,
    y = df$Y,
    window = ventana,
    checkdup = TRUE
  )
}


resumen_geometrico_grupo <- function(
  df,
  tipo_grupo,
  grupo
) {

  n <- nrow(
    df
  )

  if (n < n_min_geometria) {
    return(
      tibble(
        Tipo_grupo = tipo_grupo,
        Grupo = grupo,
        n = n,
        X_centroide = NA_real_,
        Y_centroide = NA_real_,
        NN_min_m = NA_real_,
        NN_media_m = NA_real_,
        NN_mediana_m = NA_real_,
        NN_max_m = NA_real_,
        Dispersion_X_m = NA_real_,
        Dispersion_Y_m = NA_real_,
        Area_convexa_m2 = NA_real_,
        Intensidad_nodos_m2 = n /
          area.owin(
            ventana
          ),
        Analisis_multiescala = FALSE
      )
    )
  }

  patron <- crear_ppp_grupo(
    df
  )

  dist_nn <- nndist(
    patron
  )

  area_convexa <- if (n >= 3) {

    poligono <- st_convex_hull(
      st_union(
        st_as_sf(
          df,
          coords = c(
            "X",
            "Y"
          ),
          crs = st_crs(
            cultivo_sf
          ),
          remove = FALSE
        )
      )
    )

    as.numeric(
      st_area(
        poligono
      )
    )

  } else {

    NA_real_
  }

  tibble(
    Tipo_grupo = tipo_grupo,
    Grupo = grupo,
    n = n,
    X_centroide = mean(
      df$X
    ),
    Y_centroide = mean(
      df$Y
    ),
    NN_min_m = min(
      dist_nn
    ),
    NN_media_m = mean(
      dist_nn
    ),
    NN_mediana_m = median(
      dist_nn
    ),
    NN_max_m = max(
      dist_nn
    ),
    Dispersion_X_m = sd(
      df$X
    ),
    Dispersion_Y_m = sd(
      df$Y
    ),
    Area_convexa_m2 = area_convexa,
    Intensidad_nodos_m2 = n /
      area.owin(
        ventana
      ),
    Analisis_multiescala = n >=
      n_min_multiescala
  )
}


clark_evans_grupo <- function(
  df,
  tipo_grupo,
  grupo
) {

  n <- nrow(
    df
  )

  if (n < n_min_multiescala) {

    return(
      tibble(
        Tipo_grupo = tipo_grupo,
        Grupo = grupo,
        n = n,
        Media_NN_observada_m = NA_real_,
        Media_NN_CSR_teorica_m = NA_real_,
        R_Clark_Evans = NA_real_,
        Sim_Q025_m = NA_real_,
        Sim_Mediana_m = NA_real_,
        Sim_Q975_m = NA_real_,
        p_Monte_Carlo_informativo = NA_real_,
        Clasificacion = paste0(
          "No estimado: n < ",
          n_min_multiescala
        )
      )
    )
  }

  patron <- crear_ppp_grupo(
    df
  )

  area <- area.owin(
    ventana
  )

  intensidad <- n /
    area

  media_observada <- mean(
    nndist(
      patron
    )
  )

  media_csr_teorica <- 1 /
    (
      2 *
        sqrt(
          intensidad
        )
    )

  R <- media_observada /
    media_csr_teorica

  medias_simuladas <- map_dbl(
    seq_len(
      nsim_clark
    ),
    function(i) {

      sim <- runifpoint(
        n = n,
        win = ventana
      )

      mean(
        nndist(
          sim
        )
      )
    }
  )

  q <- quantile(
    medias_simuladas,
    probs = c(
      0.025,
      0.50,
      0.975
    ),
    names = FALSE,
    na.rm = TRUE
  )

  p_mc <- (
    1 +
      sum(
        abs(
          medias_simuladas -
            median(
              medias_simuladas
            )
        ) >=
          abs(
            media_observada -
              median(
                medias_simuladas
              )
          )
      )
  ) /
    (
      nsim_clark +
        1
    )

  clasificacion <- case_when(
    media_observada >
      q[3] ~
      "Más regular que CSR",

    media_observada <
      q[1] ~
      "Más agregado que CSR",

    TRUE ~
      "Compatible con CSR"
  )

  tibble(
    Tipo_grupo = tipo_grupo,
    Grupo = grupo,
    n = n,
    Media_NN_observada_m = media_observada,
    Media_NN_CSR_teorica_m = media_csr_teorica,
    R_Clark_Evans = R,
    Sim_Q025_m = q[1],
    Sim_Mediana_m = q[2],
    Sim_Q975_m = q[3],
    p_Monte_Carlo_informativo = p_mc,
    Clasificacion = clasificacion
  )
}


envolvente_grupo <- function(
  df,
  tipo_grupo,
  grupo,
  funcion = c(
    "K",
    "L"
  )
) {

  funcion <- match.arg(
    funcion
  )

  n <- nrow(
    df
  )

  if (n < n_min_multiescala) {
    return(
      tibble()
    )
  }

  patron <- crear_ppp_grupo(
    df
  )

  simulador <- expression(
    runifpoint(
      patron$n,
      win = ventana
    )
  )

  fun_obj <- if (
    funcion == "K"
  ) {
    Kest
  } else {
    Lest
  }

  resultado <- try(
    envelope(
      Y = patron,
      fun = fun_obj,
      nsim = nsim_envolventes,
      simulate = simulador,
      correction = "translation",
      savefuns = TRUE,
      savepatterns = FALSE,
      verbose = FALSE
    ),
    silent = TRUE
  )

  if (
    inherits(
      resultado,
      "try-error"
    )
  ) {

    warning(
      "No se pudo calcular ",
      funcion,
      " para ",
      tipo_grupo,
      " = ",
      grupo
    )

    return(
      tibble()
    )
  }

  df_env <- as.data.frame(
    resultado
  )

  columnas <- c(
    "r",
    "obs",
    "lo",
    "hi"
  )

  if (
    !all(
      columnas %in%
        names(
          df_env
        )
    )
  ) {
    return(
      tibble()
    )
  }

  df_env |>
    transmute(
      Tipo_grupo = tipo_grupo,
      Grupo = grupo,
      n = n,
      Funcion = funcion,
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


analizar_conjunto_grupos <- function(
  datos,
  columna_grupo,
  tipo_grupo
) {

  grupos <- split(
    datos,
    datos[[columna_grupo]]
  )

  resumen_geometrico <- imap_dfr(
    grupos,
    function(
      df,
      grupo
    ) {

      resumen_geometrico_grupo(
        df = df,
        tipo_grupo = tipo_grupo,
        grupo = grupo
      )
    }
  )

  clark <- imap_dfr(
    grupos,
    function(
      df,
      grupo
    ) {

      clark_evans_grupo(
        df = df,
        tipo_grupo = tipo_grupo,
        grupo = grupo
      )
    }
  )

  envolventes <- imap_dfr(
    grupos,
    function(
      df,
      grupo
    ) {

      bind_rows(
        envolvente_grupo(
          df = df,
          tipo_grupo = tipo_grupo,
          grupo = grupo,
          funcion = "K"
        ),
        envolvente_grupo(
          df = df,
          tipo_grupo = tipo_grupo,
          grupo = grupo,
          funcion = "L"
        )
      )
    }
  )

  list(
    grupos = grupos,
    resumen_geometrico = resumen_geometrico,
    clark = clark,
    envolventes = envolventes
  )
}


## 5. Análisis por dosis ---------------------------------------

analisis_dosis <- analizar_conjunto_grupos(
  datos = datos,
  columna_grupo = "Grupo_dosis",
  tipo_grupo = "Dosis"
)


## 6. Análisis por repetición ----------------------------------

analisis_repeticion <- analizar_conjunto_grupos(
  datos = datos,
  columna_grupo = "Grupo_repeticion",
  tipo_grupo = "Repeticion"
)


## 7. Descripción por tratamiento ------------------------------

grupos_tratamiento <- split(
  datos,
  datos$Grupo_tratamiento
)

resumen_tratamientos <- imap_dfr(
  grupos_tratamiento,
  function(
    df,
    grupo
  ) {

    resumen_geometrico_grupo(
      df = df,
      tipo_grupo = "Tratamiento",
      grupo = grupo
    )
  }
) |>
  left_join(
    datos |>
      distinct(
        Grupo_tratamiento,
        Treatment_base,
        Omission_description,
        Grupo_dosis
      ),
    by = c(
      "Grupo" =
        "Grupo_tratamiento"
    )
  ) |>
  mutate(
    Justificacion = case_when(
      n < n_min_multiescala ~
        paste0(
          "Descripción geométrica únicamente; n < ",
          n_min_multiescala
        ),

      TRUE ~
        "Apto para análisis multiescala"
    )
  )


## 8. Consolidar resultados ------------------------------------

tabla_geometria <- bind_rows(
  analisis_dosis$resumen_geometrico,
  analisis_repeticion$resumen_geometrico,
  resumen_tratamientos
)

tabla_clark <- bind_rows(
  analisis_dosis$clark,
  analisis_repeticion$clark
)

tabla_envolventes <- bind_rows(
  analisis_dosis$envolventes,
  analisis_repeticion$envolventes
)


## 9. Resumir K y L por grupo ----------------------------------

resumen_envolventes <- tabla_envolventes |>
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
  ) |>
  group_by(
    Tipo_grupo,
    Grupo,
    n,
    Funcion
  ) |>
  summarise(
    Porcentaje_dentro = 100 *
      mean(
        Posicion ==
          "Dentro_envolvente"
      ),
    Porcentaje_sobre = 100 *
      mean(
        Posicion ==
          "Sobre_envolvente"
      ),
    Porcentaje_bajo = 100 *
      mean(
        Posicion ==
          "Bajo_envolvente"
      ),
    r_primera_salida_m = ifelse(
      any(
        Posicion !=
          "Dentro_envolvente"
      ),
      min(
        r[
          Posicion !=
            "Dentro_envolvente"
        ],
        na.rm = TRUE
      ),
      NA_real_
    ),
    .groups = "drop"
  ) |>
  mutate(
    Interpretacion = case_when(
      Porcentaje_bajo >
        Porcentaje_sobre ~
        "Tendencia de regularidad",

      Porcentaje_sobre >
        Porcentaje_bajo ~
        "Tendencia de agregación",

      TRUE ~
        "Sin dirección dominante"
    )
  )


## 10. Densidad Kernel por dosis -------------------------------

crear_kernel_grupo <- function(
  df,
  grupo
) {

  patron <- crear_ppp_grupo(
    df
  )

  imagen <- density.ppp(
    patron,
    sigma = sigma_global,
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

  as.data.frame(
    imagen
  ) |>
    setNames(
      c(
        "X",
        "Y",
        "Densidad_nodos_m2"
      )
    ) |>
    mutate(
      Grupo = grupo
    )
}

kernel_dosis <- imap_dfr(
  analisis_dosis$grupos,
  crear_kernel_grupo
)


## 11. Tabla de síntesis ---------------------------------------

sintesis_grupos <- tabla_clark |>
  left_join(
    resumen_envolventes |>
      filter(
        Funcion == "L"
      ) |>
      select(
        Tipo_grupo,
        Grupo,
        L_Porcentaje_dentro =
          Porcentaje_dentro,
        L_Porcentaje_sobre =
          Porcentaje_sobre,
        L_Porcentaje_bajo =
          Porcentaje_bajo,
        L_Interpretacion =
          Interpretacion
      ),
    by = c(
      "Tipo_grupo",
      "Grupo"
    )
  ) |>
  left_join(
    tabla_geometria |>
      select(
        Tipo_grupo,
        Grupo,
        NN_media_m,
        Area_convexa_m2,
        Intensidad_nodos_m2
      ),
    by = c(
      "Tipo_grupo",
      "Grupo"
    )
  ) |>
  mutate(
    Conclusion = case_when(
      str_detect(
        Clasificacion,
        "regular"
      ) &
        str_detect(
          L_Interpretacion,
          "regularidad"
        ) ~
        "Regularidad consistente entre Clark–Evans y L",

      str_detect(
        Clasificacion,
        "agregado"
      ) &
        str_detect(
          L_Interpretacion,
          "agregación"
        ) ~
        "Agregación consistente entre Clark–Evans y L",

      is.na(
        L_Interpretacion
      ) ~
        Clasificacion,

      TRUE ~
        "Resultados dependientes de la escala; interpretar conjuntamente"
    )
  )


## 12. Figuras: distribución por dosis -------------------------

lote_borde <- st_boundary(
  st_union(
    lote
  )
)

g_dosis <- ggplot() +
  geom_sf(
    data = lote_borde,
    fill = NA,
    linewidth = 0.55
  ) +
  geom_point(
    data = datos,
    aes(
      x = X,
      y = Y,
      shape = Grupo_dosis
    ),
    size = 2.2
  ) +
  coord_sf(
    expand = FALSE
  ) +
  labs(
    title = "Distribución espacial de los nodos según dosis",
    subtitle = paste0(
      "D100: ",
      sum(
        datos$Grupo_dosis ==
          "D100"
      ),
      " nodos; D75: ",
      sum(
        datos$Grupo_dosis ==
          "D75"
      ),
      " nodos"
    ),
    x = "X (m)",
    y = "Y (m)",
    shape = "Dosis"
  ) +
  theme_minimal()

ggsave(
  filename = file.path(
    ruta_figuras,
    "01_distribucion_nodos_por_dosis.png"
  ),
  plot = g_dosis,
  width = 8,
  height = 7,
  dpi = 400
)


g_kernel_dosis <- ggplot(
  kernel_dosis,
  aes(
    x = X,
    y = Y,
    fill = Densidad_nodos_m2
  )
) +
  geom_raster() +
  geom_sf(
    data = lote_borde,
    inherit.aes = FALSE,
    fill = NA,
    linewidth = 0.45
  ) +
  geom_point(
    data = datos,
    aes(
      x = X,
      y = Y
    ),
    inherit.aes = FALSE,
    size = 0.7
  ) +
  facet_wrap(
    ~ Grupo,
    ncol = 2
  ) +
  coord_sf(
    expand = FALSE
  ) +
  labs(
    title = "Densidad Kernel por dosis",
    subtitle = paste0(
      "Sigma común = ",
      round(
        sigma_global,
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
    "02_densidad_kernel_por_dosis.png"
  ),
  plot = g_kernel_dosis,
  width = 12,
  height = 6,
  dpi = 400
)


## 13. Figuras: Clark–Evans por grupo --------------------------

g_clark <- ggplot(
  tabla_clark,
  aes(
    x = Grupo,
    y = R_Clark_Evans,
    shape = Tipo_grupo
  )
) +
  geom_hline(
    yintercept = 1,
    linetype = "dashed"
  ) +
  geom_point(
    size = 3
  ) +
  facet_wrap(
    ~ Tipo_grupo,
    scales = "free_x",
    ncol = 1
  ) +
  labs(
    title = "Índice de Clark–Evans por dosis y repetición",
    subtitle = "R > 1 indica mayor separación media que la referencia CSR",
    x = "Grupo",
    y = "R de Clark–Evans",
    shape = "Tipo"
  ) +
  theme_minimal()

ggsave(
  filename = file.path(
    ruta_figuras,
    "03_clark_evans_por_grupo.png"
  ),
  plot = g_clark,
  width = 9,
  height = 8,
  dpi = 400
)


## 14. Figuras: K y L por dosis --------------------------------

crear_grafico_envolvente <- function(
  tabla,
  funcion,
  tipo_grupo
) {

  ggplot(
    tabla |>
      filter(
        Funcion == funcion,
        Tipo_grupo == tipo_grupo
      ),
    aes(
      x = r
    )
  ) +
    geom_ribbon(
      aes(
        ymin = Limite_inferior,
        ymax = Limite_superior
      ),
      alpha = 0.30
    ) +
    geom_line(
      aes(
        y = Observado
      ),
      linewidth = 0.65
    ) +
    facet_wrap(
      ~ Grupo,
      scales = "free_y"
    ) +
    labs(
      title = paste0(
        "Función ",
        funcion,
        " por ",
        tolower(
          tipo_grupo
        )
      ),
      x = "r (m)",
      y = funcion
    ) +
    theme_minimal()
}

if (
  nrow(
    tabla_envolventes |>
      filter(
        Tipo_grupo == "Dosis"
      )
  ) > 0
) {

  g_K_dosis <- crear_grafico_envolvente(
    tabla = tabla_envolventes,
    funcion = "K",
    tipo_grupo = "Dosis"
  )

  g_L_dosis <- crear_grafico_envolvente(
    tabla = tabla_envolventes,
    funcion = "L",
    tipo_grupo = "Dosis"
  )

  panel_KL_dosis <- g_K_dosis +
    g_L_dosis +
    plot_annotation(
      title = "Comparación multiescala entre dosis"
    )

  ggsave(
    filename = file.path(
      ruta_figuras,
      "04_panel_K_L_por_dosis.png"
    ),
    plot = panel_KL_dosis,
    width = 13,
    height = 7,
    dpi = 400
  )
}


## 15. Figura descriptiva por tratamiento ----------------------

g_tratamientos <- ggplot(
  datos,
  aes(
    x = X,
    y = Y
  )
) +
  geom_point(
    size = 1.8
  ) +
  facet_wrap(
    ~ Treatment_ID,
    ncol = 4
  ) +
  coord_equal() +
  labs(
    title = "Localización de los nodos por tratamiento",
    subtitle = paste0(
      "Panel descriptivo; no se estiman K/L con grupos de ",
      "2–3 nodos"
    ),
    x = "X (m)",
    y = "Y (m)"
  ) +
  theme_minimal()

ggsave(
  filename = file.path(
    ruta_figuras,
    "05_panel_localizacion_por_tratamiento.png"
  ),
  plot = g_tratamientos,
  width = 14,
  height = 18,
  dpi = 400
)


## 16. Panel de síntesis ---------------------------------------

panel_sintesis <- (
  g_dosis +
    g_kernel_dosis
) /
  g_clark +
  plot_annotation(
    title = "Síntesis de patrones espaciales por dosis y repetición"
  )

ggsave(
  filename = file.path(
    ruta_figuras,
    "06_panel_sintesis_patrones_por_grupo.png"
  ),
  plot = panel_sintesis,
  width = 14,
  height = 14,
  dpi = 400
)


## 17. Exportar datos ------------------------------------------

write_csv(
  tabla_geometria,
  file.path(
    ruta_datos,
    "Resumen_geometrico_grupos.csv"
  )
)

write_csv(
  tabla_clark,
  file.path(
    ruta_datos,
    "Clark_Evans_por_grupo.csv"
  )
)

write_csv(
  tabla_envolventes,
  file.path(
    ruta_datos,
    "Envolventes_K_L_por_grupo.csv"
  )
)

write_csv(
  kernel_dosis,
  file.path(
    ruta_datos,
    "Kernel_por_dosis.csv"
  )
)


## 18. Exportar libro ------------------------------------------

write_xlsx(
  list(
    Contexto_global =
      tabla_global_clark,

    Sintesis_grupos =
      sintesis_grupos,

    Clark_Evans =
      tabla_clark,

    Resumen_K_L =
      resumen_envolventes,

    Geometria_dosis_repeticion =
      tabla_geometria |>
        filter(
          Tipo_grupo !=
            "Tratamiento"
        ),

    Tratamientos_descriptivos =
      resumen_tratamientos,

    Envolventes_K_L =
      tabla_envolventes
  ),
  path = file.path(
    ruta_tablas,
    "Resultados_patrones_por_dosis_repeticion_tratamiento.xlsx"
  )
)


## 19. Guardar objeto para Script 05-C -------------------------

saveRDS(
  list(
    datos = datos,
    cultivo_sf = cultivo_sf,
    lote = lote,
    ventana = ventana,
    sigma_global = sigma_global,
    tabla_geometria = tabla_geometria,
    tabla_clark = tabla_clark,
    tabla_envolventes = tabla_envolventes,
    resumen_envolventes = resumen_envolventes,
    sintesis_grupos = sintesis_grupos,
    kernel_dosis = kernel_dosis
  ),
  file = file.path(
    ruta_objetos,
    "objetos_patrones_por_grupo.rds"
  )
)


## 20. Registro de sesión --------------------------------------

capture.output(
  sessionInfo(),
  file = file.path(
    ruta_salida,
    "sessionInfo_05_B.txt"
  )
)


## 21. Resumen en consola --------------------------------------

cat(
  "\n============================================\n"
)

cat(
  "SCRIPT 05-B FINALIZADO\n\n"
)

cat(
  "Clark–Evans por dosis y repetición:\n"
)

print(
  tabla_clark |>
    select(
      Tipo_grupo,
      Grupo,
      n,
      R_Clark_Evans,
      Clasificacion
    )
)

cat(
  "\nTratamientos descritos: ",
  nrow(
    resumen_tratamientos
  ),
  "\n"
)

cat(
  "Tratamientos con K/L formal: 0, salvo que alguno alcance n >= ",
  n_min_multiescala,
  "\n\n"
)

cat(
  "Objeto para Script 05-C:\n",
  file.path(
    ruta_objetos,
    "objetos_patrones_por_grupo.rds"
  ),
  "\n"
)

cat(
  "============================================\n"
)
