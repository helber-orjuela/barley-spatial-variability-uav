
## ============================================================
## 06_SUELO_CULTIVO – SCRIPT 06-B
## REDUNDANCIA, COLINEALIDAD Y SELECCIÓN DE COVARIABLES
## ============================================================
##
## Entrada:
##   06_SUELO_CULTIVO/01_BASE_INTEGRADA_AUDITADA/objetos_R/
##   objetos_base_integrada_auditada_06A.rds
##
## Salida:
##   06_SUELO_CULTIVO/02_SELECCION_COVARIABLES/
##
## PROPÓSITO
##   1. Cuantificar redundancia entre las variables de suelo.
##   2. Evitar incorporar simultáneamente covariables equivalentes.
##   3. Evaluar correlaciones, VIF y número de condición.
##   4. Construir conjuntos candidatos parsimoniosos.
##   5. Identificar covariables aptas para los modelos integrados 06-C.
##
## DECISIONES METODOLÓGICAS
##   - CE y CE corregida OLS no se incluyen juntas.
##   - N, P y K no se incluyen simultáneamente si presentan
##     redundancia casi perfecta.
##   - pH y Hv se conservan como variables exploratorias.
##   - Se evalúan conjuntos alternativos, no un único modelo global.
##   - El PCA se usa como descripción de la estructura de redundancia;
##     no sustituye automáticamente las variables originales.
##
## CRITERIOS ORIENTATIVOS
##   - |r| o |rho| >= 0.90: redundancia muy alta.
##   - VIF >= 5: colinealidad relevante.
##   - VIF >= 10: colinealidad severa.
##   - Número de condición >= 30: inestabilidad importante.
##
## NOTA
##   Las variables N, P y K son índices instrumentales nominales.
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
  "patchwork",
  "tibble"
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
library(patchwork)
library(tibble)


## 1. Parámetros ------------------------------------------------

umbral_cor_alta <- 0.90
umbral_cor_moderada <- 0.70
umbral_vif_relevante <- 5
umbral_vif_severo <- 10
umbral_condicion <- 30


## 2. Rutas ----------------------------------------------------

ruta_base <- paste0(
  "D:/UDCA_2017_CLASE/46_ARTICULOS_2025/",
  "ARTICULO_CEBADA_DOCTORADO/AJUSTES_ARTÍCULO"
)

ruta_objeto <- file.path(
  ruta_base,
  "06_SUELO_CULTIVO",
  "01_BASE_INTEGRADA_AUDITADA",
  "objetos_R",
  "objetos_base_integrada_auditada_06A.rds"
)

ruta_salida <- file.path(
  ruta_base,
  "06_SUELO_CULTIVO",
  "02_SELECCION_COVARIABLES"
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

if (!file.exists(ruta_objeto)) {
  stop(
    "No existe el objeto del Script 06-A:\n",
    ruta_objeto
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


## 3. Cargar objeto --------------------------------------------

obj <- readRDS(
  ruta_objeto
)

base <- obj$base_auditada
variables_suelo_obj <- obj$variables_suelo
variables_cultivo <- obj$variables_cultivo
variables_diseno <- obj$variables_diseno
diccionario <- obj$diccionario

if (is.null(base)) {
  stop(
    "El objeto no contiene 'base_auditada'."
  )
}

base <- as_tibble(
  base
)

if (
  "Estado_auditoria" %in% names(base)
) {
  base_modelo <- base |>
    filter(
      Estado_auditoria == "Apto"
    )
} else {
  base_modelo <- base
}

if (nrow(base_modelo) == 0) {
  stop(
    "No existen filas aptas para el análisis."
  )
}


## 4. Variables de suelo disponibles ---------------------------

variables_suelo_esperadas <- c(
  "Suelo_CE",
  "Suelo_CE_corregida_OLS",
  "Suelo_N",
  "Suelo_P",
  "Suelo_K",
  "Suelo_pH",
  "Suelo_T",
  "Suelo_Hv"
)

variables_suelo <- intersect(
  variables_suelo_esperadas,
  names(base_modelo)
)

if (length(variables_suelo) < 2) {
  stop(
    "Se requieren al menos dos variables de suelo."
  )
}

variables_numericas_validas <- variables_suelo[
  map_lgl(
    base_modelo[
      variables_suelo
    ],
    is.numeric
  )
]

variables_suelo <- variables_numericas_validas

variacion <- map_dbl(
  base_modelo[
    variables_suelo
  ],
  function(x) {
    sd(
      x,
      na.rm = TRUE
    )
  }
)

variables_suelo <- variables_suelo[
  is.finite(
    variacion
  ) &
    variacion > 0
]

if (length(variables_suelo) < 2) {
  stop(
    "Menos de dos variables de suelo tienen variación."
  )
}


## 5. Funciones auxiliares -------------------------------------

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

  if (sum(validos) < 4) {
    return(
      NA_real_
    )
  }

  suppressWarnings(
    cor(
      x[validos],
      y[validos],
      method = metodo
    )
  )
}


matriz_a_largo <- function(
  matriz,
  metodo
) {

  as.data.frame(
    as.table(
      matriz
    ),
    stringsAsFactors = FALSE
  ) |>
    as_tibble() |>
    setNames(
      c(
        "Variable_1",
        "Variable_2",
        "Correlacion"
      )
    ) |>
    mutate(
      Metodo = metodo
    )
}


calcular_vif <- function(
  df,
  variables
) {

  if (length(variables) == 1) {
    return(
      tibble(
        Variable = variables,
        VIF = 1,
        Tolerancia = 1
      )
    )
  }

  map_dfr(
    variables,
    function(variable_respuesta) {

      variables_predictoras <- setdiff(
        variables,
        variable_respuesta
      )

      formula_auxiliar <- as.formula(
        paste0(
          variable_respuesta,
          " ~ ",
          paste(
            variables_predictoras,
            collapse = " + "
          )
        )
      )

      modelo <- try(
        lm(
          formula_auxiliar,
          data = df
        ),
        silent = TRUE
      )

      if (inherits(modelo, "try-error")) {
        return(
          tibble(
            Variable = variable_respuesta,
            VIF = Inf,
            Tolerancia = 0
          )
        )
      }

      r2 <- summary(
        modelo
      )$r.squared

      vif <- ifelse(
        is.finite(r2) &&
          r2 < 1,
        1 / (1 - r2),
        Inf
      )

      tibble(
        Variable = variable_respuesta,
        VIF = vif,
        Tolerancia = ifelse(
          is.finite(vif) &&
            vif != 0,
          1 / vif,
          0
        )
      )
    }
  )
}


calcular_numero_condicion <- function(
  df,
  variables
) {

  X <- df |>
    select(
      all_of(
        variables
      )
    ) |>
    mutate(
      across(
        everything(),
        as.numeric
      )
    ) |>
    as.data.frame()

  X <- X[
    complete.cases(X),
    ,
    drop = FALSE
  ]

  if (
    nrow(X) < 3 ||
      ncol(X) < 2
  ) {
    return(
      NA_real_
    )
  }

  X_escalada <- scale(
    X
  )

  singular_values <- svd(
    X_escalada,
    nu = 0,
    nv = 0
  )$d

  singular_values <- singular_values[
    singular_values > sqrt(
      .Machine$double.eps
    )
  ]

  if (length(singular_values) < 2) {
    return(
      Inf
    )
  }

  max(
    singular_values
  ) /
    min(
      singular_values
    )
}


evaluar_conjunto <- function(
  nombre,
  variables
) {

  variables_presentes <- intersect(
    variables,
    variables_suelo
  )

  if (length(variables_presentes) == 0) {
    return(
      list(
        resumen = tibble(
          Conjunto = nombre,
          Variables = NA_character_,
          n_variables = 0,
          VIF_maximo = NA_real_,
          VIF_medio = NA_real_,
          Numero_condicion = NA_real_,
          Estado = "No disponible"
        ),
        vif = tibble()
      )
    )
  }

  vif_tabla <- calcular_vif(
    df = base_modelo,
    variables = variables_presentes
  ) |>
    mutate(
      Conjunto = nombre,
      .before = 1
    )

  numero_condicion <- calcular_numero_condicion(
    df = base_modelo,
    variables = variables_presentes
  )

  vif_max <- max(
    vif_tabla$VIF,
    na.rm = TRUE
  )

  vif_medio <- mean(
    vif_tabla$VIF[
      is.finite(
        vif_tabla$VIF
      )
    ],
    na.rm = TRUE
  )

  estado <- case_when(
    any(
      !is.finite(
        vif_tabla$VIF
      )
    ) ~
      "No apto: singularidad",

    vif_max >= umbral_vif_severo ~
      "No apto: VIF severo",

    numero_condicion >= umbral_condicion ~
      "Revisar: condición alta",

    vif_max >= umbral_vif_relevante ~
      "Revisar: VIF relevante",

    TRUE ~
      "Apto"
  )

  resumen <- tibble(
    Conjunto = nombre,
    Variables = paste(
      variables_presentes,
      collapse = " + "
    ),
    n_variables = length(
      variables_presentes
    ),
    VIF_maximo = vif_max,
    VIF_medio = vif_medio,
    Numero_condicion = numero_condicion,
    Estado = estado
  )

  list(
    resumen = resumen,
    vif = vif_tabla
  )
}


## 6. Correlaciones entre variables de suelo -------------------

matriz_pearson <- cor(
  base_modelo[
    variables_suelo
  ],
  use = "pairwise.complete.obs",
  method = "pearson"
)

matriz_spearman <- cor(
  base_modelo[
    variables_suelo
  ],
  use = "pairwise.complete.obs",
  method = "spearman"
)

correlaciones_largas <- bind_rows(
  matriz_a_largo(
    matriz_pearson,
    "Pearson"
  ),
  matriz_a_largo(
    matriz_spearman,
    "Spearman"
  )
) |>
  filter(
    Variable_1 != Variable_2
  ) |>
  mutate(
    Par_ordenado = map2_chr(
      Variable_1,
      Variable_2,
      function(a, b) {
        paste(
          sort(
            c(
              a,
              b
            )
          ),
          collapse = "__"
        )
      }
    )
  ) |>
  group_by(
    Metodo,
    Par_ordenado
  ) |>
  slice(
    1
  ) |>
  ungroup() |>
  select(
    -Par_ordenado
  ) |>
  mutate(
    Magnitud = abs(
      Correlacion
    ),
    Categoria = case_when(
      Magnitud >= umbral_cor_alta ~
        "Redundancia muy alta",

      Magnitud >= umbral_cor_moderada ~
        "Asociación alta",

      Magnitud >= 0.40 ~
        "Asociación moderada",

      Magnitud >= 0.20 ~
        "Asociación limitada",

      TRUE ~
        "Asociación débil"
    )
  ) |>
  arrange(
    Metodo,
    desc(
      Magnitud
    )
  )


## 7. Pares redundantes ----------------------------------------

pares_redundantes <- correlaciones_largas |>
  filter(
    Metodo == "Spearman",
    Magnitud >= umbral_cor_alta
  ) |>
  mutate(
    Recomendacion = case_when(
      (
        Variable_1 == "Suelo_CE" &
          Variable_2 == "Suelo_CE_corregida_OLS"
      ) |
        (
          Variable_2 == "Suelo_CE" &
            Variable_1 == "Suelo_CE_corregida_OLS"
        ) ~
        "No incluir CE original y corregida en el mismo modelo",

      Variable_1 %in% c(
        "Suelo_N",
        "Suelo_P",
        "Suelo_K"
      ) &
        Variable_2 %in% c(
          "Suelo_N",
          "Suelo_P",
          "Suelo_K"
        ) ~
        "Evaluar N, P y K en modelos alternativos separados",

      TRUE ~
        "Evitar inclusión simultánea sin justificación"
    )
  )


## 8. Conjuntos candidatos -------------------------------------

conjuntos_candidatos <- list(
  C01_CE_original = c(
    "Suelo_CE",
    "Suelo_pH",
    "Suelo_T",
    "Suelo_Hv"
  ),

  C02_CE_corregida = c(
    "Suelo_CE_corregida_OLS",
    "Suelo_pH",
    "Suelo_T",
    "Suelo_Hv"
  ),

  C03_Indice_N = c(
    "Suelo_N",
    "Suelo_pH",
    "Suelo_T",
    "Suelo_Hv"
  ),

  C04_Indice_P = c(
    "Suelo_P",
    "Suelo_pH",
    "Suelo_T",
    "Suelo_Hv"
  ),

  C05_Indice_K = c(
    "Suelo_K",
    "Suelo_pH",
    "Suelo_T",
    "Suelo_Hv"
  ),

  C06_CE_T = c(
    "Suelo_CE",
    "Suelo_T"
  ),

  C07_CEcor_T = c(
    "Suelo_CE_corregida_OLS",
    "Suelo_T"
  ),

  C08_N_T = c(
    "Suelo_N",
    "Suelo_T"
  ),

  C09_P_T = c(
    "Suelo_P",
    "Suelo_T"
  ),

  C10_K_T = c(
    "Suelo_K",
    "Suelo_T"
  ),

  C11_pH_T_Hv = c(
    "Suelo_pH",
    "Suelo_T",
    "Suelo_Hv"
  )
)

evaluaciones <- imap(
  conjuntos_candidatos,
  function(
    variables,
    nombre
  ) {
    evaluar_conjunto(
      nombre = nombre,
      variables = variables
    )
  }
)

tabla_conjuntos <- bind_rows(
  map(
    evaluaciones,
    "resumen"
  )
)

tabla_vif <- bind_rows(
  map(
    evaluaciones,
    "vif"
  )
)

conjuntos_aptos <- tabla_conjuntos |>
  filter(
    Estado == "Apto"
  ) |>
  arrange(
    VIF_maximo,
    Numero_condicion
  )


## 9. Evaluación individual de covariables ---------------------

tabla_covariables <- tibble(
  Variable = variables_suelo
) |>
  mutate(
    Tipo = case_when(
      Variable %in% c(
        "Suelo_CE",
        "Suelo_CE_corregida_OLS"
      ) ~
        "Conductividad eléctrica",

      Variable %in% c(
        "Suelo_N",
        "Suelo_P",
        "Suelo_K"
      ) ~
        "Índice instrumental nominal",

      Variable == "Suelo_pH" ~
        "pH exploratorio",

      Variable == "Suelo_T" ~
        "Temperatura",

      Variable == "Suelo_Hv" ~
        "Humedad instrumental",

      TRUE ~
        "Otra"
    ),
    Prioridad = case_when(
      Variable == "Suelo_T" ~
        "Covariable complementaria",

      Variable == "Suelo_CE" ~
        "Covariable principal: escenario original",

      Variable == "Suelo_CE_corregida_OLS" ~
        "Covariable principal: escenario corregido",

      Variable %in% c(
        "Suelo_N",
        "Suelo_P",
        "Suelo_K"
      ) ~
        "Evaluar por separado",

      Variable == "Suelo_pH" ~
        "Exploratoria",

      Variable == "Suelo_Hv" ~
        "Exploratoria con estructura espacial débil",

      TRUE ~
        "Secundaria"
    ),
    Restriccion = case_when(
      Variable == "Suelo_CE" ~
        "No combinar con Suelo_CE_corregida_OLS",

      Variable == "Suelo_CE_corregida_OLS" ~
        "No combinar con Suelo_CE",

      Variable %in% c(
        "Suelo_N",
        "Suelo_P",
        "Suelo_K"
      ) ~
        "No combinar simultáneamente con los otros índices NPK",

      Variable == "Suelo_pH" ~
        "Validación limitada del sensor",

      Variable == "Suelo_Hv" ~
        "Efecto pepita o dependencia espacial débil",

      TRUE ~
        "Sin restricción adicional"
    )
  )


## 10. PCA descriptivo -----------------------------------------

base_pca <- base_modelo |>
  select(
    all_of(
      variables_suelo
    )
  ) |>
  drop_na()

pca_resultado <- NULL
pca_varianza <- tibble()
pca_cargas <- tibble()
pca_puntajes <- tibble()

if (
  nrow(base_pca) >= 5 &&
    ncol(base_pca) >= 2
) {

  pca_resultado <- prcomp(
    base_pca,
    center = TRUE,
    scale. = TRUE
  )

  varianza <- pca_resultado$sdev^2

  pca_varianza <- tibble(
    Componente = paste0(
      "PC",
      seq_along(
        varianza
      )
    ),
    Varianza = varianza,
    Porcentaje = 100 *
      varianza /
      sum(
        varianza
      ),
    Porcentaje_acumulado = cumsum(
      100 *
        varianza /
        sum(
          varianza
        )
    )
  )

  pca_cargas <- as.data.frame(
    pca_resultado$rotation
  ) |>
    tibble::rownames_to_column(
      var = "Variable"
    ) |>
    as_tibble()

  pca_puntajes <- as.data.frame(
    pca_resultado$x
  ) |>
    mutate(
      Nodo = base_modelo$Nodo[
        complete.cases(
          base_modelo[
            variables_suelo
          ]
        )
      ],
      .before = 1
    ) |>
    as_tibble()
}


## 11. Figuras -------------------------------------------------

cor_plot <- correlaciones_largas |>
  filter(
    Metodo == "Spearman"
  ) |>
  bind_rows(
    tibble(
      Variable_1 = variables_suelo,
      Variable_2 = variables_suelo,
      Correlacion = 1,
      Metodo = "Spearman",
      Magnitud = 1,
      Categoria = "Diagonal"
    )
  )

g_cor <- ggplot(
  cor_plot,
  aes(
    x = Variable_1,
    y = Variable_2,
    fill = Correlacion
  )
) +
  geom_tile() +
  geom_text(
    aes(
      label = sprintf(
        "%.2f",
        Correlacion
      )
    ),
    size = 3
  ) +
  labs(
    title = "Redundancia entre covariables de suelo",
    subtitle = "Correlación de Spearman",
    x = NULL,
    y = NULL,
    fill = "rho"
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
    "01_correlaciones_covariables_suelo.png"
  ),
  plot = g_cor,
  width = 11,
  height = 9,
  dpi = 400
)


g_vif <- ggplot(
  tabla_vif,
  aes(
    x = VIF,
    y = reorder(
      Variable,
      VIF
    )
  )
) +
  geom_vline(
    xintercept = umbral_vif_relevante,
    linetype = "dashed"
  ) +
  geom_vline(
    xintercept = umbral_vif_severo,
    linetype = "dotted"
  ) +
  geom_point() +
  facet_wrap(
    ~ Conjunto,
    scales = "free_y"
  ) +
  labs(
    title = "VIF por conjunto candidato",
    x = "VIF",
    y = NULL
  ) +
  theme_minimal()

ggsave(
  filename = file.path(
    ruta_figuras,
    "02_VIF_conjuntos_candidatos.png"
  ),
  plot = g_vif,
  width = 14,
  height = 12,
  dpi = 400
)


g_conjuntos <- ggplot(
  tabla_conjuntos,
  aes(
    x = Numero_condicion,
    y = VIF_maximo,
    shape = Estado
  )
) +
  geom_vline(
    xintercept = umbral_condicion,
    linetype = "dashed"
  ) +
  geom_hline(
    yintercept = umbral_vif_relevante,
    linetype = "dashed"
  ) +
  geom_point(
    size = 3
  ) +
  geom_text(
    aes(
      label = Conjunto
    ),
    nudge_y = 0.15,
    size = 3,
    check_overlap = TRUE
  ) +
  labs(
    title = "Estabilidad de los conjuntos candidatos",
    subtitle = "Menor VIF y menor número de condición indican mayor estabilidad",
    x = "Número de condición",
    y = "VIF máximo",
    shape = "Estado"
  ) +
  theme_minimal()

ggsave(
  filename = file.path(
    ruta_figuras,
    "03_estabilidad_conjuntos_candidatos.png"
  ),
  plot = g_conjuntos,
  width = 10,
  height = 7,
  dpi = 400
)


if (nrow(pca_varianza) > 0) {

  g_pca_var <- ggplot(
    pca_varianza,
    aes(
      x = Componente,
      y = Porcentaje
    )
  ) +
    geom_col() +
    geom_line(
      aes(
        y = Porcentaje_acumulado,
        group = 1
      )
    ) +
    geom_point(
      aes(
        y = Porcentaje_acumulado
      )
    ) +
    labs(
      title = "PCA descriptivo de covariables de suelo",
      x = "Componente",
      y = "Varianza explicada o acumulada (%)"
    ) +
    theme_minimal()

  ggsave(
    filename = file.path(
      ruta_figuras,
      "04_PCA_varianza_explicada.png"
    ),
    plot = g_pca_var,
    width = 9,
    height = 6,
    dpi = 400
  )
}


## 12. Recomendación final para 06-C ---------------------------

recomendacion_modelos <- tibble(
  Escenario = c(
    "M2_CE_original",
    "M2_CE_corregida",
    "M2_Indice_N",
    "M2_Indice_P",
    "M2_Indice_K",
    "M2_Exploratorio_pH_T_Hv"
  ),
  Covariables = c(
    paste(
      intersect(
        c(
          "Suelo_CE",
          "Suelo_T"
        ),
        variables_suelo
      ),
      collapse = " + "
    ),
    paste(
      intersect(
        c(
          "Suelo_CE_corregida_OLS",
          "Suelo_T"
        ),
        variables_suelo
      ),
      collapse = " + "
    ),
    paste(
      intersect(
        c(
          "Suelo_N",
          "Suelo_T"
        ),
        variables_suelo
      ),
      collapse = " + "
    ),
    paste(
      intersect(
        c(
          "Suelo_P",
          "Suelo_T"
        ),
        variables_suelo
      ),
      collapse = " + "
    ),
    paste(
      intersect(
        c(
          "Suelo_K",
          "Suelo_T"
        ),
        variables_suelo
      ),
      collapse = " + "
    ),
    paste(
      intersect(
        c(
          "Suelo_pH",
          "Suelo_T",
          "Suelo_Hv"
        ),
        variables_suelo
      ),
      collapse = " + "
    )
  ),
  Uso = c(
    "Escenario principal con CE sin calibrar",
    "Sensibilidad con CE corregida OLS",
    "Índice N evaluado por separado",
    "Índice P evaluado por separado",
    "Índice K evaluado por separado",
    "Escenario exploratorio"
  ),
  Advertencia = c(
    "No combinar con CE corregida",
    "No combinar con CE original",
    "Índice instrumental nominal",
    "Índice instrumental nominal",
    "Índice instrumental nominal",
    "pH y Hv con validación limitada"
  )
) |>
  filter(
    Covariables != ""
  )


## 13. Exportar resultados -------------------------------------

write_csv(
  correlaciones_largas,
  file.path(
    ruta_datos,
    "Correlaciones_covariables_suelo.csv"
  )
)

write_csv(
  pares_redundantes,
  file.path(
    ruta_datos,
    "Pares_redundantes_covariables.csv"
  )
)

write_csv(
  tabla_conjuntos,
  file.path(
    ruta_datos,
    "Evaluacion_conjuntos_candidatos.csv"
  )
)

write_csv(
  tabla_vif,
  file.path(
    ruta_datos,
    "VIF_conjuntos_candidatos.csv"
  )
)

write_csv(
  recomendacion_modelos,
  file.path(
    ruta_datos,
    "Conjuntos_recomendados_para_06C.csv"
  )
)

if (nrow(pca_varianza) > 0) {
  write_csv(
    pca_varianza,
    file.path(
      ruta_datos,
      "PCA_varianza.csv"
    )
  )

  write_csv(
    pca_cargas,
    file.path(
      ruta_datos,
      "PCA_cargas.csv"
    )
  )
}


write_xlsx(
  list(
    Resumen_conjuntos =
      tabla_conjuntos,

    Conjuntos_aptos =
      conjuntos_aptos,

    Recomendacion_06C =
      recomendacion_modelos,

    Pares_redundantes =
      pares_redundantes,

    Correlaciones =
      correlaciones_largas,

    VIF =
      tabla_vif,

    Covariables =
      tabla_covariables,

    PCA_varianza =
      pca_varianza,

    PCA_cargas =
      pca_cargas
  ),
  path = file.path(
    ruta_tablas,
    "Resultados_redundancia_seleccion_covariables.xlsx"
  )
)


## 14. Guardar objeto para Script 06-C -------------------------

saveRDS(
  list(
    base_modelo = base_modelo,
    variables_suelo = variables_suelo,
    variables_cultivo = variables_cultivo,
    variables_diseno = variables_diseno,
    correlaciones = correlaciones_largas,
    pares_redundantes = pares_redundantes,
    conjuntos_candidatos = conjuntos_candidatos,
    tabla_conjuntos = tabla_conjuntos,
    tabla_vif = tabla_vif,
    conjuntos_aptos = conjuntos_aptos,
    recomendacion_modelos = recomendacion_modelos,
    pca = pca_resultado,
    pca_varianza = pca_varianza,
    pca_cargas = pca_cargas
  ),
  file = file.path(
    ruta_objetos,
    "objetos_seleccion_covariables_06B.rds"
  )
)


## 15. Registro de sesión --------------------------------------

capture.output(
  sessionInfo(),
  file = file.path(
    ruta_salida,
    "sessionInfo_06_B.txt"
  )
)


## 16. Resumen en consola --------------------------------------

cat(
  "\n============================================\n"
)

cat(
  "SCRIPT 06-B FINALIZADO\n\n"
)

cat(
  "Variables de suelo evaluadas:\n",
  paste(
    variables_suelo,
    collapse = ", "
  ),
  "\n\n"
)

cat(
  "Pares con redundancia muy alta:\n"
)

print(
  pares_redundantes |>
    select(
      Variable_1,
      Variable_2,
      Correlacion,
      Recomendacion
    )
)

cat(
  "\nConjuntos aptos:\n"
)

print(
  conjuntos_aptos
)

cat(
  "\nConjuntos recomendados para 06-C:\n"
)

print(
  recomendacion_modelos
)

cat(
  "\nObjeto para Script 06-C:\n",
  file.path(
    ruta_objetos,
    "objetos_seleccion_covariables_06B.rds"
  ),
  "\n"
)

cat(
  "============================================\n"
)
