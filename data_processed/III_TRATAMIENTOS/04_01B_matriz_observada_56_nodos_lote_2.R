
## ============================================================
## 04_TRATAMIENTOS – SCRIPT 01B CORREGIDO
## MATRIZ OBSERVADA DE TRATAMIENTOS PARA EL LOTE ANALIZADO
## ============================================================
##
## Entrada:
##   AJUSTES_ARTÍCULO/04_TRATAMIENTOS/DATOS_BASE/
##   Nodos_cultivo_lote_2_por_tratamiento.xlsx
##
## Hoja:
##   Nodos_cultivo_lote_2
##
## Alcance:
##   - Trabaja exclusivamente con los 56 nodos de cultivo medidos
##     en el lote analizado.
##   - No compara con un diseño teórico de 72 unidades.
##   - No clasifica los 16 registros restantes como faltantes,
##     porque pertenecen a otro lote.
##   - Reconoce 24 combinaciones tratamiento × dosis.
##   - Documenta el desbalance real de repeticiones:
##       A = 100%: 32 nodos
##       B = 75%: 24 nodos
##
## Propósito:
##   1. Limpiar y estandarizar códigos.
##   2. Construir la matriz observada del lote.
##   3. Auditar tratamientos, dosis y repeticiones.
##   4. Resumir variables de cultivo por unidad observada.
##   5. Generar figuras y objeto para los siguientes scripts.
## ============================================================


## 0. Paquetes -------------------------------------------------

paquetes <- c(
  "readxl",
  "readr",
  "writexl",
  "dplyr",
  "tidyr",
  "purrr",
  "stringr",
  "ggplot2",
  "sf"
)

instalar <- paquetes[
  !paquetes %in% rownames(installed.packages())
]

if (length(instalar) > 0) {
  install.packages(instalar)
}

library(readxl)
library(readr)
library(writexl)
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(ggplot2)
library(sf)


## 1. Rutas ----------------------------------------------------

ruta_base <- paste0(
  "D:/UDCA_2017_CLASE/46_ARTICULOS_2025/",
  "ARTICULO_CEBADA_DOCTORADO/AJUSTES_ARTÍCULO"
)

ruta_tratamientos <- file.path(
  ruta_base,
  "04_TRATAMIENTOS"
)

ruta_datos_base <- file.path(
  ruta_tratamientos,
  "DATOS_BASE"
)

archivos_xlsx <- list.files(
  ruta_datos_base,
  pattern = "Nodos_cultivo_lote_2_por_tratamiento.*\\.xlsx$",
  full.names = TRUE,
  ignore.case = TRUE
)

if (length(archivos_xlsx) == 0) {
  stop(
    "No se encontró el archivo de nodos por tratamiento en:\n",
    ruta_datos_base
  )
}

if (length(archivos_xlsx) > 1) {
  warning(
    "Se encontraron varios archivos compatibles. Se usará:\n",
    archivos_xlsx[1]
  )
}

ruta_nodos <- archivos_xlsx[1]

ruta_salida <- file.path(
  ruta_tratamientos,
  "01_MATRIZ_MAESTRA",
  "nodos_tratamientos_lote_2"
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


## 2. Leer archivo por posición --------------------------------

datos_crudos <- read_excel(
  ruta_nodos,
  sheet = "Nodos_cultivo_lote_2",
  .name_repair = "unique"
)

if (ncol(datos_crudos) < 13) {
  stop(
    "El archivo tiene menos de 13 columnas; no es posible ",
    "localizar tratamiento y repetición."
  )
}

# Extraer columnas por posición de forma segura.
# .data[[...]] exige nombres de columna, no índices numéricos.
nombres_columnas <- names(
  datos_crudos
)

nodos <- tibble(
  Nodo = as.integer(
    datos_crudos[[nombres_columnas[1]]]
  ),
  X = as.numeric(
    datos_crudos[[nombres_columnas[2]]]
  ),
  Y = as.numeric(
    datos_crudos[[nombres_columnas[3]]]
  ),
  Macollas = as.numeric(
    datos_crudos[[nombres_columnas[4]]]
  ),
  Altura = as.numeric(
    datos_crudos[[nombres_columnas[5]]]
  ),
  Diametro = as.numeric(
    datos_crudos[[nombres_columnas[6]]]
  ),
  NDVI = as.numeric(
    datos_crudos[[nombres_columnas[7]]]
  ),
  Nitratos = as.numeric(
    datos_crudos[[nombres_columnas[8]]]
  ),
  Potasio_savia = as.numeric(
    datos_crudos[[nombres_columnas[9]]]
  ),
  Calcio_savia = as.numeric(
    datos_crudos[[nombres_columnas[10]]]
  ),
  Codigo_original = as.character(
    datos_crudos[[nombres_columnas[12]]]
  ),
  Replication = as.integer(
    datos_crudos[[nombres_columnas[13]]]
  )
) |>
  filter(
    !is.na(Nodo)
  )


## 3. Catálogo y normalización ---------------------------------

normalizar_omision <- function(x) {

  x <- str_to_upper(
    str_trim(
      x
    )
  )

  case_when(
    x %in% c("CONTROL", "CONTR") ~ "Control",
    x %in% c("COMPLETO", "COMP", "COMPL") ~ "Completo",
    x == "N" ~ "Sin N",
    x == "P" ~ "Sin P",
    x == "K" ~ "Sin K",
    x %in% c("MG", "MG Y S", "MG+S", "MGS") ~ "Sin Mg y S",
    x == "FE" ~ "Sin Fe",
    x == "ZN" ~ "Sin Zn",
    x == "MN" ~ "Sin Mn",
    x %in% c("B", "BORO") ~ "Sin B",
    x == "CU" ~ "Sin Cu",
    x == "MO" ~ "Sin Mo",
    TRUE ~ NA_character_
  )
}

catalogo <- tribble(
  ~Treatment_base, ~Omission_description, ~Omission_public,
  "T0", "Control", "Control",
  "T1", "Completo", "Complete",
  "T2", "Sin N", "Minus_N",
  "T3", "Sin P", "Minus_P",
  "T4", "Sin K", "Minus_K",
  "T5", "Sin Mg y S", "Minus_MgS",
  "T6", "Sin Fe", "Minus_Fe",
  "T7", "Sin Zn", "Minus_Zn",
  "T8", "Sin Mn", "Minus_Mn",
  "T9", "Sin B", "Minus_B",
  "T10", "Sin Cu", "Minus_Cu",
  "T11", "Sin Mo", "Minus_Mo"
)

nodos_limpios <- nodos |>
  mutate(
    Codigo_original = str_to_upper(
      str_trim(
        Codigo_original
      )
    ),
    Dose_code = str_extract(
      Codigo_original,
      "^[AB]"
    ),
    Etiqueta_omision = str_remove(
      Codigo_original,
      "^[AB]-"
    ),
    Omission_description = normalizar_omision(
      Etiqueta_omision
    ),
    Dose_percent = case_when(
      Dose_code == "A" ~ 100,
      Dose_code == "B" ~ 75,
      TRUE ~ NA_real_
    )
  ) |>
  left_join(
    catalogo,
    by = "Omission_description"
  ) |>
  mutate(
    Treatment_ID = ifelse(
      !is.na(Treatment_base) &
        !is.na(Dose_code),
      paste0(
        Treatment_base,
        "-",
        Dose_code
      ),
      NA_character_
    ),
    Public_treatment_code = ifelse(
      !is.na(Omission_public) &
        !is.na(Dose_percent),
      paste0(
        "D",
        Dose_percent,
        "_",
        Omission_public
      ),
      NA_character_
    ),
    Public_unit_code = ifelse(
      !is.na(Public_treatment_code) &
        !is.na(Replication),
      paste0(
        Public_treatment_code,
        "_R",
        Replication
      ),
      NA_character_
    ),
    Experimental_Unit = paste0(
      "L2_NODE_",
      str_pad(
        Nodo,
        width = 2,
        pad = "0"
      )
    ),
    Geometry_status = case_when(
      is.finite(X) &
        is.finite(Y) ~ "Point_available",
      TRUE ~ "Missing_coordinates"
    )
  ) |>
  select(
    Nodo,
    Experimental_Unit,
    X,
    Y,
    Codigo_original,
    Treatment_ID,
    Treatment_base,
    Dose_code,
    Dose_percent,
    Omission_description,
    Omission_public,
    Replication,
    Public_treatment_code,
    Public_unit_code,
    Geometry_status,
    Macollas,
    Altura,
    Diametro,
    NDVI,
    Nitratos,
    Potasio_savia,
    Calcio_savia
  )


## 4. Auditoría del lote observado -----------------------------

codigos_no_reconocidos <- nodos_limpios |>
  filter(
    is.na(Treatment_ID) |
      is.na(Dose_percent) |
      is.na(Replication)
  )

duplicados_nodo <- nodos_limpios |>
  count(
    Nodo,
    name = "n"
  ) |>
  filter(
    n > 1
  )

duplicados_unidad <- nodos_limpios |>
  count(
    Treatment_ID,
    Replication,
    name = "n_nodos"
  ) |>
  filter(
    n_nodos > 1
  )

conteo_tratamientos <- nodos_limpios |>
  count(
    Treatment_ID,
    Treatment_base,
    Dose_code,
    Dose_percent,
    Omission_description,
    name = "n_nodos"
  ) |>
  arrange(
    Dose_percent,
    Treatment_base
  )

conteo_repeticiones <- nodos_limpios |>
  count(
    Dose_code,
    Dose_percent,
    Replication,
    name = "n_nodos"
  ) |>
  arrange(
    Dose_percent,
    Replication
  )

resumen_lote <- tibble(
  Indicador = c(
    "Nodos observados en el lote",
    "Tratamientos base observados",
    "Combinaciones tratamiento-dosis observadas",
    "Nodos con dosis 100%",
    "Nodos con dosis 75%",
    "Nodos de repetición 1",
    "Nodos de repetición 2",
    "Nodos de repetición 3",
    "Códigos no reconocidos",
    "Nodos duplicados"
  ),
  Valor = c(
    nrow(
      nodos_limpios
    ),
    n_distinct(
      nodos_limpios$Treatment_base
    ),
    n_distinct(
      nodos_limpios$Treatment_ID
    ),
    sum(
      nodos_limpios$Dose_percent == 100,
      na.rm = TRUE
    ),
    sum(
      nodos_limpios$Dose_percent == 75,
      na.rm = TRUE
    ),
    sum(
      nodos_limpios$Replication == 1,
      na.rm = TRUE
    ),
    sum(
      nodos_limpios$Replication == 2,
      na.rm = TRUE
    ),
    sum(
      nodos_limpios$Replication == 3,
      na.rm = TRUE
    ),
    nrow(
      codigos_no_reconocidos
    ),
    nrow(
      duplicados_nodo
    )
  )
)

print(
  resumen_lote
)


## 5. Matriz observada de 56 unidades --------------------------

matriz_observada <- nodos_limpios |>
  mutate(
    Area_m2 = 4,
    Lote = "Lote_2",
    Observation_status = "Observed"
  ) |>
  select(
    Nodo,
    Experimental_Unit,
    Lote,
    Treatment_ID,
    Treatment_base,
    Dose_code,
    Dose_percent,
    Omission_description,
    Omission_public,
    Replication,
    Area_m2,
    Public_treatment_code,
    Public_unit_code,
    X,
    Y,
    Geometry_status,
    Observation_status,
    Macollas,
    Altura,
    Diametro,
    NDVI,
    Nitratos,
    Potasio_savia,
    Calcio_savia
  )


## 6. Resumen de cultivo por tratamiento -----------------------

resumen_tratamiento <- matriz_observada |>
  group_by(
    Treatment_ID,
    Treatment_base,
    Dose_code,
    Dose_percent,
    Omission_description
  ) |>
  summarise(
    n = n(),
    Repeticiones = paste(
      sort(
        unique(
          Replication
        )
      ),
      collapse = ", "
    ),
    Macollas_media = mean(
      Macollas,
      na.rm = TRUE
    ),
    Altura_media = mean(
      Altura,
      na.rm = TRUE
    ),
    Diametro_media = mean(
      Diametro,
      na.rm = TRUE
    ),
    NDVI_media = mean(
      NDVI,
      na.rm = TRUE
    ),
    Nitratos_media = mean(
      Nitratos,
      na.rm = TRUE
    ),
    Potasio_savia_media = mean(
      Potasio_savia,
      na.rm = TRUE
    ),
    Calcio_savia_media = mean(
      Calcio_savia,
      na.rm = TRUE
    ),
    .groups = "drop"
  )


## 7. Figura espacial ------------------------------------------

g_distribucion <- ggplot(
  matriz_observada,
  aes(
    x = X,
    y = Y
  )
) +
  geom_point(
    aes(
      shape = factor(
        Dose_percent
      )
    ),
    size = 2.5
  ) +
  geom_text(
    aes(
      label = paste0(
        Treatment_base,
        "-R",
        Replication
      )
    ),
    size = 2.2,
    nudge_y = 0.35,
    check_overlap = TRUE
  ) +
  coord_equal() +
  labs(
    title =
      "Distribución espacial de las 56 unidades observadas",
    subtitle =
      "Lote 2; A=100% y B=75%",
    x = "X (m)",
    y = "Y (m)",
    shape = "Dosis (%)"
  ) +
  theme_minimal()

ggsave(
  filename = file.path(
    ruta_figuras,
    "Distribucion_56_nodos_por_tratamiento.png"
  ),
  plot = g_distribucion,
  width = 11,
  height = 8,
  dpi = 400
)


## 8. Figura de cobertura real por combinación ----------------

g_conteo <- ggplot(
  conteo_tratamientos,
  aes(
    x = Treatment_ID,
    y = n_nodos,
    fill = factor(
      Dose_percent
    )
  )
) +
  geom_col() +
  geom_hline(
    yintercept = 3,
    linetype = "dashed"
  ) +
  labs(
    title =
      "Número de observaciones por combinación de tratamiento",
    subtitle =
      "La línea discontinua representa tres repeticiones",
    x = "Tratamiento",
    y = "Número de nodos",
    fill = "Dosis (%)"
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
    "Conteo_observaciones_por_tratamiento.png"
  ),
  plot = g_conteo,
  width = 12,
  height = 7,
  dpi = 400
)


## 9. Exportar -------------------------------------------------

write_csv(
  matriz_observada,
  file.path(
    ruta_datos,
    "Matriz_observada_56_nodos_lote_2.csv"
  )
)

write_csv(
  resumen_tratamiento,
  file.path(
    ruta_datos,
    "Resumen_cultivo_por_tratamiento_lote_2.csv"
  )
)

write_xlsx(
  list(
    Resumen_lote =
      resumen_lote,

    Matriz_56_nodos =
      matriz_observada,

    Conteo_tratamientos =
      conteo_tratamientos,

    Conteo_repeticiones =
      conteo_repeticiones,

    Resumen_tratamiento =
      resumen_tratamiento,

    Codigos_no_reconocidos =
      codigos_no_reconocidos,

    Duplicados_nodo =
      duplicados_nodo,

    Unidades_con_submuestras =
      duplicados_unidad,

    Catalogo =
      catalogo
  ),
  path = file.path(
    ruta_tablas,
    "Resultados_matriz_observada_56_nodos.xlsx"
  )
)


## 10. Guardar objeto ------------------------------------------

nodos_sf <- st_as_sf(
  matriz_observada,
  coords = c(
    "X",
    "Y"
  ),
  crs = 32618,
  remove = FALSE
)

saveRDS(
  list(
    matriz_observada =
      matriz_observada,

    nodos_sf =
      nodos_sf,

    resumen_tratamiento =
      resumen_tratamiento,

    conteo_tratamientos =
      conteo_tratamientos,

    catalogo =
      catalogo
  ),
  file = file.path(
    ruta_objetos,
    "objetos_tratamientos_56_nodos.rds"
  )
)


## 11. Registro de sesión --------------------------------------

capture.output(
  sessionInfo(),
  file = file.path(
    ruta_salida,
    "sessionInfo_04_01B_56_nodos.txt"
  )
)


## 12. Resumen en consola --------------------------------------

cat(
  "\n============================================\n"
)

cat(
  "SCRIPT 04-01B CORREGIDO FINALIZADO\n\n"
)

cat(
  "El análisis se restringió al lote con 56 nodos.\n\n"
)

print(
  resumen_lote
)

cat(
  "\nObjeto para los scripts siguientes:\n",
  file.path(
    ruta_objetos,
    "objetos_tratamientos_56_nodos.rds"
  ),
  "\n"
)

cat(
  "============================================\n"
)
