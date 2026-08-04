
## ============================================================
## 06_SUELO_CULTIVO – SCRIPT 06-A
## CONSOLIDACIÓN Y AUDITORÍA DE LA BASE INTEGRADA
## ============================================================
##
## Entrada principal:
##   04_TRATAMIENTOS/02_INTEGRACION_SUELO_CULTIVO/objetos_R/
##   objetos_integracion_suelo_cultivo_tratamientos.rds
##
## Entradas opcionales:
##   04_KRIGING_I2/objetos_R/objetos_kriging_I2.rds
##   05_PATRON_PUNTOS/03_ASOCIACION_SUELO_CULTIVO/objetos_R/
##   objetos_asociacion_espacial_suelo_cultivo.rds
##
## Salida:
##   06_SUELO_CULTIVO/01_BASE_INTEGRADA_AUDITADA/
##
## PROPÓSITO
##   1. Consolidar en una sola tabla los 56 nodos de cultivo.
##   2. Integrar diseño experimental, cultivo y suelo predicho.
##   3. Incorporar incertidumbre de Kriging cuando esté disponible.
##   4. Verificar correspondencia uno a uno, CRS, duplicados,
##      valores faltantes, rangos y coherencia espacial.
##   5. Generar un diccionario de datos y una matriz de trazabilidad.
##
## PRECAUCIONES
##   - Las variables de suelo en los nodos de cultivo son predicciones
##     obtenidas mediante Kriging, no mediciones independientes.
##   - N, P y K se conservan como índices instrumentales nominales.
##   - CE corregida OLS proviene de una calibración exploratoria.
##   - pH y Hv requieren interpretación prudente.
## ============================================================


## 0. Paquetes -------------------------------------------------

paquetes <- c(
  "sf",
  "dplyr",
  "tidyr",
  "purrr",
  "stringr",
  "readr",
  "writexl",
  "ggplot2",
  "patchwork"
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
library(readr)
library(writexl)
library(ggplot2)
library(patchwork)


## 1. Rutas ----------------------------------------------------

ruta_base <- paste0(
  "D:/UDCA_2017_CLASE/46_ARTICULOS_2025/",
  "ARTICULO_CEBADA_DOCTORADO/AJUSTES_ARTÍCULO"
)

ruta_integracion <- file.path(
  ruta_base,
  "04_TRATAMIENTOS",
  "02_INTEGRACION_SUELO_CULTIVO",
  "objetos_R",
  "objetos_integracion_suelo_cultivo_tratamientos.rds"
)

ruta_kriging <- file.path(
  ruta_base,
  "04_KRIGING_I2",
  "objetos_R",
  "objetos_kriging_I2.rds"
)

ruta_asociacion <- file.path(
  ruta_base,
  "05_PATRON_PUNTOS",
  "03_ASOCIACION_SUELO_CULTIVO",
  "objetos_R",
  "objetos_asociacion_espacial_suelo_cultivo.rds"
)

ruta_salida <- file.path(
  ruta_base,
  "06_SUELO_CULTIVO",
  "01_BASE_INTEGRADA_AUDITADA"
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

if (!file.exists(ruta_integracion)) {
  stop(
    "No existe el objeto principal de integración:\n",
    ruta_integracion
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

obj_integracion <- readRDS(
  ruta_integracion
)

obj_kriging <- if (file.exists(ruta_kriging)) {
  readRDS(ruta_kriging)
} else {
  NULL
}

obj_asociacion <- if (file.exists(ruta_asociacion)) {
  readRDS(ruta_asociacion)
} else {
  NULL
}


## 3. Extraer componentes principales --------------------------

base_integrada <- obj_integracion$base_integrada
cultivo_sf <- obj_integracion$cultivo_sf
lote <- obj_integracion$lote

if (is.null(base_integrada)) {
  stop(
    "El objeto principal no contiene 'base_integrada'."
  )
}

if (is.null(cultivo_sf)) {
  stop(
    "El objeto principal no contiene 'cultivo_sf'."
  )
}

if (is.null(lote)) {
  stop(
    "El objeto principal no contiene 'lote'."
  )
}

base_integrada <- as_tibble(
  base_integrada
)

cultivo_sf <- st_as_sf(
  cultivo_sf
)

lote <- st_as_sf(
  lote
)


## 4. Verificación de identificadores y coordenadas ------------

columnas_minimas <- c(
  "Nodo",
  "X",
  "Y"
)

faltantes_minimos <- setdiff(
  columnas_minimas,
  names(base_integrada)
)

if (length(faltantes_minimos) > 0) {
  stop(
    "Faltan columnas mínimas en base_integrada: ",
    paste(
      faltantes_minimos,
      collapse = ", "
    )
  )
}

base_integrada <- base_integrada |>
  mutate(
    Nodo = as.integer(
      Nodo
    ),
    X = as.numeric(
      X
    ),
    Y = as.numeric(
      Y
    )
  )

if (anyDuplicated(base_integrada$Nodo) > 0) {
  stop(
    "Existen identificadores Nodo duplicados."
  )
}

if (
  anyDuplicated(
    base_integrada |>
      select(
        X,
        Y
      )
  ) > 0
) {
  stop(
    "Existen coordenadas X-Y duplicadas."
  )
}


## 5. Normalizar CRS y coherencia espacial ---------------------

if (is.na(st_crs(cultivo_sf))) {
  stop(
    "cultivo_sf no tiene CRS definido."
  )
}

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

cultivo_sf <- cultivo_sf |>
  arrange(
    match(
      Nodo,
      base_integrada$Nodo
    )
  )

coords_sf <- st_coordinates(
  cultivo_sf
)

if (nrow(coords_sf) != nrow(base_integrada)) {
  stop(
    "El número de geometrías no coincide con la base integrada."
  )
}

diferencia_X <- abs(
  coords_sf[, 1] -
    base_integrada$X
)

diferencia_Y <- abs(
  coords_sf[, 2] -
    base_integrada$Y
)

tolerancia_coordenadas_m <- 0.10

coordenadas_coherentes <- (
  diferencia_X <= tolerancia_coordenadas_m
) &
  (
    diferencia_Y <= tolerancia_coordenadas_m
  )

dentro_lote <- lengths(
  st_intersects(
    cultivo_sf,
    lote_union
  )
) > 0


## 6. Catálogo de variables esperadas --------------------------

variables_diseno <- intersect(
  c(
    "Experimental_Unit",
    "Treatment_ID",
    "Treatment_base",
    "Dose_percent",
    "Omission_description",
    "Replication"
  ),
  names(base_integrada)
)

variables_cultivo <- intersect(
  c(
    "Macollas",
    "Altura",
    "Diametro",
    "NDVI",
    "Nitratos",
    "Potasio_savia",
    "Calcio_savia"
  ),
  names(base_integrada)
)

variables_suelo <- intersect(
  c(
    "Suelo_CE",
    "Suelo_CE_corregida_OLS",
    "Suelo_N",
    "Suelo_P",
    "Suelo_K",
    "Suelo_pH",
    "Suelo_T",
    "Suelo_Hv"
  ),
  names(base_integrada)
)

variables_control <- intersect(
  c(
    "Distancia_nodo_suelo_m",
    "Nodo_suelo_cercano",
    "Dentro_lote"
  ),
  names(base_integrada)
)


## 7. Incorporar incertidumbre de Kriging ----------------------

extraer_incertidumbre_kriging <- function(
  objeto_kriging,
  nodos_sf
) {

  if (is.null(objeto_kriging)) {
    return(
      tibble()
    )
  }

  resultados <- objeto_kriging$resultados_kriging

  if (is.null(resultados)) {
    return(
      tibble()
    )
  }

  variables_disponibles <- names(
    resultados
  )

  if (length(variables_disponibles) == 0) {
    return(
      tibble()
    )
  }

  extraer_variable <- function(variable) {

    elemento <- resultados[[variable]]

    if (is.null(elemento)) {
      return(
        NULL
      )
    }

    modelo <- elemento$modelo_final

    if (is.null(modelo)) {
      modelo <- elemento$modelo
    }

    if (is.null(modelo)) {
      return(
        NULL
      )
    }

    fuente_sf <- objeto_kriging$i2_sf

    if (is.null(fuente_sf)) {
      fuente_sf <- objeto_kriging$suelo_sf
    }

    if (is.null(fuente_sf)) {
      return(
        NULL
      )
    }

    fuente_sf <- st_as_sf(
      fuente_sf
    )

    fuente_sf <- st_transform(
      fuente_sf,
      st_crs(nodos_sf)
    )

    nodos_objetivo <- st_transform(
      nodos_sf,
      st_crs(fuente_sf)
    )

    formula_variable <- as.formula(
      paste0(
        variable,
        " ~ 1"
      )
    )

    prediccion <- try(
      gstat::krige(
        formula = formula_variable,
        locations = fuente_sf,
        newdata = nodos_objetivo,
        model = modelo,
        debug.level = 0
      ),
      silent = TRUE
    )

    if (inherits(prediccion, "try-error")) {
      return(
        NULL
      )
    }

    tibble(
      Nodo = nodos_sf$Nodo,
      Variable = variable,
      Prediccion_recalculada = as.numeric(
        prediccion$var1.pred
      ),
      Varianza_kriging = as.numeric(
        prediccion$var1.var
      ),
      Error_estandar_kriging = sqrt(
        pmax(
          as.numeric(
            prediccion$var1.var
          ),
          0
        )
      )
    )
  }

  bind_rows(
    map(
      variables_disponibles,
      extraer_variable
    )
  )
}

incertidumbre_larga <- extraer_incertidumbre_kriging(
  objeto_kriging = obj_kriging,
  nodos_sf = cultivo_sf
)

if (nrow(incertidumbre_larga) > 0) {

  incertidumbre_ancha <- incertidumbre_larga |>
    select(
      Nodo,
      Variable,
      Varianza_kriging,
      Error_estandar_kriging
    ) |>
    pivot_wider(
      names_from = Variable,
      values_from = c(
        Varianza_kriging,
        Error_estandar_kriging
      ),
      names_sep = "__"
    )

} else {

  incertidumbre_ancha <- tibble(
    Nodo = base_integrada$Nodo
  )
}


## 8. Consolidar base final ------------------------------------

base_auditada <- base_integrada |>
  left_join(
    incertidumbre_ancha,
    by = "Nodo"
  ) |>
  mutate(
    CRS = st_crs(cultivo_sf)$input,
    Dentro_lote_auditado = dentro_lote,
    Diferencia_X_geometria_m = diferencia_X,
    Diferencia_Y_geometria_m = diferencia_Y,
    Coordenadas_coherentes = coordenadas_coherentes,
    Fila_completa_diseno = if (
      length(variables_diseno) > 0
    ) {
      complete.cases(
        across(
          all_of(
            variables_diseno
          )
        )
      )
    } else {
      TRUE
    },
    Fila_completa_cultivo = if (
      length(variables_cultivo) > 0
    ) {
      complete.cases(
        across(
          all_of(
            variables_cultivo
          )
        )
      )
    } else {
      TRUE
    },
    Fila_completa_suelo = if (
      length(variables_suelo) > 0
    ) {
      complete.cases(
        across(
          all_of(
            variables_suelo
          )
        )
      )
    } else {
      TRUE
    },
    Estado_auditoria = case_when(
      !Dentro_lote_auditado ~
        "Revisar: fuera del lote",

      !Coordenadas_coherentes ~
        "Revisar: discrepancia geométrica",

      !Fila_completa_diseno ~
        "Revisar: diseño incompleto",

      !Fila_completa_cultivo ~
        "Revisar: cultivo incompleto",

      !Fila_completa_suelo ~
        "Revisar: suelo incompleto",

      TRUE ~
        "Apto"
    )
  )


## 9. Auditoría de valores faltantes ---------------------------

tabla_faltantes <- tibble(
  Variable = names(
    base_auditada
  ),
  n_total = nrow(
    base_auditada
  ),
  n_faltantes = map_int(
    base_auditada,
    function(x) {
      sum(
        is.na(
          x
        )
      )
    }
  ),
  Porcentaje_faltantes = 100 *
    n_faltantes /
    n_total
) |>
  arrange(
    desc(
      n_faltantes
    ),
    Variable
  )


## 10. Resumen numérico y rangos -------------------------------

columnas_numericas <- names(
  base_auditada
)[
  map_lgl(
    base_auditada,
    is.numeric
  )
]

resumen_numerico <- map_dfr(
  columnas_numericas,
  function(variable) {

    x <- base_auditada[[variable]]

    tibble(
      Variable = variable,
      n = sum(
        is.finite(
          x
        )
      ),
      Minimo = ifelse(
        any(
          is.finite(
            x
          )
        ),
        min(
          x,
          na.rm = TRUE
        ),
        NA_real_
      ),
      Q1 = ifelse(
        any(
          is.finite(
            x
          )
        ),
        quantile(
          x,
          0.25,
          na.rm = TRUE,
          names = FALSE
        ),
        NA_real_
      ),
      Mediana = ifelse(
        any(
          is.finite(
            x
          )
        ),
        median(
          x,
          na.rm = TRUE
        ),
        NA_real_
      ),
      Media = ifelse(
        any(
          is.finite(
            x
          )
        ),
        mean(
          x,
          na.rm = TRUE
        ),
        NA_real_
      ),
      Q3 = ifelse(
        any(
          is.finite(
            x
          )
        ),
        quantile(
          x,
          0.75,
          na.rm = TRUE,
          names = FALSE
        ),
        NA_real_
      ),
      Maximo = ifelse(
        any(
          is.finite(
            x
          )
        ),
        max(
          x,
          na.rm = TRUE
        ),
        NA_real_
      ),
      DE = ifelse(
        sum(
          is.finite(
            x
          )
        ) > 1,
        sd(
          x,
          na.rm = TRUE
        ),
        NA_real_
      ),
      CV_porcentaje = ifelse(
        is.finite(
          mean(
            x,
            na.rm = TRUE
          )
        ) &&
          mean(
            x,
            na.rm = TRUE
          ) != 0,
        100 *
          sd(
            x,
            na.rm = TRUE
          ) /
          abs(
            mean(
              x,
              na.rm = TRUE
            )
          ),
        NA_real_
      )
    )
  }
)


## 11. Detección exploratoria de valores extremos --------------

variables_extremos <- unique(
  c(
    variables_cultivo,
    variables_suelo
  )
)

tabla_extremos <- map_dfr(
  variables_extremos,
  function(variable) {

    x <- as.numeric(
      base_auditada[[variable]]
    )

    q1 <- quantile(
      x,
      0.25,
      na.rm = TRUE,
      names = FALSE
    )

    q3 <- quantile(
      x,
      0.75,
      na.rm = TRUE,
      names = FALSE
    )

    iqr <- q3 - q1

    limite_inferior <- q1 -
      1.5 *
        iqr

    limite_superior <- q3 +
      1.5 *
        iqr

    tibble(
      Nodo = base_auditada$Nodo,
      Variable = variable,
      Valor = x,
      Limite_inferior = limite_inferior,
      Limite_superior = limite_superior,
      Extremo_IQR = (
        x < limite_inferior
      ) |
        (
          x > limite_superior
        )
    )
  }
)


## 12. Diccionario y trazabilidad ------------------------------

clasificar_origen <- function(variable) {

  case_when(
    variable %in% c(
      "Nodo",
      "X",
      "Y",
      "CRS"
    ) ~
      "Identificación y geometría",

    variable %in% variables_diseno ~
      "Diseño experimental",

    variable %in% variables_cultivo ~
      "Medición de cultivo",

    variable %in% variables_suelo ~
      "Predicción de suelo por Kriging",

    str_detect(
      variable,
      "^Varianza_kriging__"
    ) ~
      "Incertidumbre de Kriging",

    str_detect(
      variable,
      "^Error_estandar_kriging__"
    ) ~
      "Incertidumbre de Kriging",

    variable %in% variables_control ~
      "Control de integración espacial",

    TRUE ~
      "Variable auxiliar o auditoría"
  )
}

clasificar_unidad <- function(variable) {

  case_when(
    variable %in% c(
      "X",
      "Y",
      "Diferencia_X_geometria_m",
      "Diferencia_Y_geometria_m",
      "Distancia_nodo_suelo_m"
    ) ~
      "m",

    variable %in% c(
      "Altura"
    ) ~
      "cm",

    variable %in% c(
      "Diametro"
    ) ~
      "mm",

    variable %in% c(
      "NDVI"
    ) ~
      "adimensional",

    variable %in% c(
      "Suelo_CE",
      "Suelo_CE_corregida_OLS"
    ) ~
      "µS/cm",

    variable %in% c(
      "Suelo_pH"
    ) ~
      "unidades de pH",

    variable %in% c(
      "Suelo_T"
    ) ~
      "°C",

    variable %in% c(
      "Suelo_Hv"
    ) ~
      "índice instrumental",

    variable %in% c(
      "Suelo_N",
      "Suelo_P",
      "Suelo_K"
    ) ~
      "índice instrumental nominal",

    variable %in% c(
      "Macollas"
    ) ~
      "conteo",

    variable %in% c(
      "Nitratos",
      "Potasio_savia",
      "Calcio_savia"
    ) ~
      "unidad del equipo de savia",

    str_detect(
      variable,
      "^Error_estandar_kriging__"
    ) ~
      "misma unidad de la variable",

    str_detect(
      variable,
      "^Varianza_kriging__"
    ) ~
      "unidad²",

    TRUE ~
      NA_character_
  )
}

diccionario <- tibble(
  Variable = names(
    base_auditada
  ),
  Clase_R = map_chr(
    base_auditada,
    function(x) {
      class(
        x
      )[1]
    }
  )
) |>
  mutate(
    Origen = map_chr(
      Variable,
      clasificar_origen
    ),
    Unidad = map_chr(
      Variable,
      clasificar_unidad
    ),
    Interpretacion = case_when(
      Variable == "Suelo_CE_corregida_OLS" ~
        "CE corregida mediante calibración OLS exploratoria",

      Variable %in% c(
        "Suelo_N",
        "Suelo_P",
        "Suelo_K"
      ) ~
        "Índice instrumental nominal; no interpretar como concentración química validada",

      Variable == "Suelo_pH" ~
        "Variable exploratoria por validación limitada del sensor",

      Variable == "Suelo_Hv" ~
        "Variable con estructura espacial débil o efecto pepita dominante",

      str_detect(
        Variable,
        "^Suelo_"
      ) ~
        "Valor predicho en el nodo de cultivo mediante modelo geoestadístico",

      TRUE ~
        NA_character_
    )
  )


## 13. Matriz de trazabilidad ----------------------------------

trazabilidad <- tibble(
  Componente = c(
    "Diseño experimental",
    "Variables de cultivo",
    "Variables de suelo I2",
    "CE corregida",
    "Incertidumbre de Kriging",
    "Geometría de nodos",
    "Polígono del lote",
    "Base residual previa"
  ),
  Fuente = c(
    basename(
      ruta_integracion
    ),
    basename(
      ruta_integracion
    ),
    basename(
      ruta_integracion
    ),
    basename(
      ruta_integracion
    ),
    ifelse(
      file.exists(
        ruta_kriging
      ),
      basename(
        ruta_kriging
      ),
      "No disponible"
    ),
    basename(
      ruta_integracion
    ),
    basename(
      ruta_integracion
    ),
    ifelse(
      file.exists(
        ruta_asociacion
      ),
      basename(
        ruta_asociacion
      ),
      "No utilizada"
    )
  ),
  Estado = c(
    ifelse(
      length(
        variables_diseno
      ) > 0,
      "Disponible",
      "No disponible"
    ),
    ifelse(
      length(
        variables_cultivo
      ) > 0,
      "Disponible",
      "No disponible"
    ),
    ifelse(
      length(
        variables_suelo
      ) > 0,
      "Disponible",
      "No disponible"
    ),
    ifelse(
      "Suelo_CE_corregida_OLS" %in% names(
        base_auditada
      ),
      "Disponible",
      "No disponible"
    ),
    ifelse(
      nrow(
        incertidumbre_larga
      ) > 0,
      "Incorporada",
      "No incorporada"
    ),
    "Disponible",
    "Disponible",
    ifelse(
      file.exists(
        ruta_asociacion
      ),
      "Disponible como soporte",
      "No disponible"
    )
  )
)


## 14. Resumen de auditoría ------------------------------------

resumen_auditoria <- tibble(
  Indicador = c(
    "Número de filas",
    "Nodos únicos",
    "Coordenadas únicas",
    "Nodos dentro del lote",
    "Nodos con coordenadas coherentes",
    "Filas aptas",
    "Filas con observación",
    "Variables de diseño",
    "Variables de cultivo",
    "Variables de suelo",
    "Variables con incertidumbre de Kriging"
  ),
  Valor = c(
    nrow(
      base_auditada
    ),
    n_distinct(
      base_auditada$Nodo
    ),
    n_distinct(
      paste(
        base_auditada$X,
        base_auditada$Y
      )
    ),
    sum(
      base_auditada$Dentro_lote_auditado
    ),
    sum(
      base_auditada$Coordenadas_coherentes
    ),
    sum(
      base_auditada$Estado_auditoria ==
        "Apto"
    ),
    sum(
      base_auditada$Estado_auditoria !=
        "Apto"
    ),
    length(
      variables_diseno
    ),
    length(
      variables_cultivo
    ),
    length(
      variables_suelo
    ),
    length(
      grep(
        "^Error_estandar_kriging__",
        names(
          base_auditada
        )
      )
    )
  )
)


## 15. Figuras de auditoría ------------------------------------

estado_figura <- base_auditada |>
  count(
    Estado_auditoria,
    name = "n"
  )

g_estado <- ggplot(
  estado_figura,
  aes(
    x = Estado_auditoria,
    y = n
  )
) +
  geom_col() +
  labs(
    title = "Estado de auditoría de los 56 nodos",
    x = NULL,
    y = "Número de nodos"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(
      angle = 30,
      hjust = 1
    )
  )

faltantes_figura <- tabla_faltantes |>
  filter(
    n_faltantes > 0
  )

if (nrow(faltantes_figura) > 0) {

  g_faltantes <- ggplot(
    faltantes_figura,
    aes(
      x = Porcentaje_faltantes,
      y = reorder(
        Variable,
        Porcentaje_faltantes
      )
    )
  ) +
    geom_col() +
    labs(
      title = "Porcentaje de valores faltantes",
      x = "Valores faltantes (%)",
      y = NULL
    ) +
    theme_minimal()

} else {

  g_faltantes <- ggplot() +
    annotate(
      "text",
      x = 0,
      y = 0,
      label = "No se detectaron valores faltantes"
    ) +
    xlim(
      -1,
      1
    ) +
    ylim(
      -1,
      1
    ) +
    labs(
      title = "Auditoría de valores faltantes"
    ) +
    theme_void()
}

g_coordenadas <- ggplot(
  base_auditada,
  aes(
    x = X,
    y = Y,
    shape = Estado_auditoria
  )
) +
  geom_point(
    size = 2
  ) +
  coord_equal() +
  labs(
    title = "Coherencia espacial de la base integrada",
    x = "X (m)",
    y = "Y (m)",
    shape = "Estado"
  ) +
  theme_minimal()

panel_auditoria <- (
  g_estado +
    g_faltantes
) /
  g_coordenadas +
  plot_annotation(
    title = "Auditoría de la base integrada suelo–cultivo"
  )

ggsave(
  filename = file.path(
    ruta_figuras,
    "Panel_auditoria_base_integrada.png"
  ),
  plot = panel_auditoria,
  width = 13,
  height = 11,
  dpi = 400
)


## 16. Exportar base y tablas ----------------------------------

write_csv(
  base_auditada,
  file.path(
    ruta_datos,
    "Base_maestra_suelo_cultivo_auditada.csv"
  )
)

write_csv(
  incertidumbre_larga,
  file.path(
    ruta_datos,
    "Incertidumbre_kriging_por_nodo.csv"
  )
)

write_csv(
  tabla_extremos,
  file.path(
    ruta_datos,
    "Valores_extremos_exploratorios.csv"
  )
)

write_xlsx(
  list(
    Resumen_auditoria =
      resumen_auditoria,

    Base_maestra =
      base_auditada,

    Valores_faltantes =
      tabla_faltantes,

    Resumen_numerico =
      resumen_numerico,

    Valores_extremos =
      tabla_extremos,

    Diccionario =
      diccionario,

    Trazabilidad =
      trazabilidad,

    Incertidumbre_kriging =
      incertidumbre_larga
  ),
  path = file.path(
    ruta_tablas,
    "Resultados_auditoria_base_integrada_suelo_cultivo.xlsx"
  )
)


## 17. Guardar objeto para Script 06-B -------------------------

saveRDS(
  list(
    base_auditada = base_auditada,
    cultivo_sf = cultivo_sf,
    lote = lote_union,
    variables_diseno = variables_diseno,
    variables_cultivo = variables_cultivo,
    variables_suelo = variables_suelo,
    incertidumbre_kriging = incertidumbre_larga,
    resumen_auditoria = resumen_auditoria,
    diccionario = diccionario,
    trazabilidad = trazabilidad
  ),
  file = file.path(
    ruta_objetos,
    "objetos_base_integrada_auditada_06A.rds"
  )
)


## 18. Registro de sesión --------------------------------------

capture.output(
  sessionInfo(),
  file = file.path(
    ruta_salida,
    "sessionInfo_06_A.txt"
  )
)


## 19. Resumen en consola --------------------------------------

cat(
  "\n============================================\n"
)

cat(
  "SCRIPT 06-A FINALIZADO\n\n"
)

cat(
  "Resumen de auditoría:\n"
)

print(
  resumen_auditoria
)

cat(
  "\nVariables de diseño:\n",
  paste(
    variables_diseno,
    collapse = ", "
  ),
  "\n"
)

cat(
  "\nVariables de cultivo:\n",
  paste(
    variables_cultivo,
    collapse = ", "
  ),
  "\n"
)

cat(
  "\nVariables de suelo:\n",
  paste(
    variables_suelo,
    collapse = ", "
  ),
  "\n"
)

cat(
  "\nObjeto para Script 06-B:\n",
  file.path(
    ruta_objetos,
    "objetos_base_integrada_auditada_06A.rds"
  ),
  "\n"
)

cat(
  "============================================\n"
)
