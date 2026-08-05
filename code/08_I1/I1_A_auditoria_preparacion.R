
## ============================================================
## I1 – SCRIPT A. AUDITORÍA Y PREPARACIÓN COMPACTA
## ============================================================
## Entrada:
##   AJUSTES_ARTÍCULO/DATOS BASE/Nodos_Lote_2.csv
##
## Salida:
##   AJUSTES_ARTÍCULO/06_I1_COMPACTO/01_AUDITORIA
##
## Propósito:
##   Preparar I1 para comparación temporal exploratoria con I2.
##   No usa p-valores como criterio principal.
## ============================================================

paquetes <- c(
  "readr", "dplyr", "tidyr", "purrr",
  "ggplot2", "e1071", "writexl"
)

instalar <- paquetes[!paquetes %in% rownames(installed.packages())]
if (length(instalar) > 0) install.packages(instalar)

library(readr)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(e1071)
library(writexl)

ruta_base <- paste0(
  "D:/UDCA_2017_CLASE/46_ARTICULOS_2025/",
  "ARTICULO_CEBADA_DOCTORADO/AJUSTES_ARTÍCULO"
)

ruta_entrada <- file.path(
  ruta_base,
  "DATOS BASE",
  "Nodos_Lote_2.csv"
)

ruta_salida <- file.path(
  ruta_base,
  "06_I1_COMPACTO",
  "01_AUDITORIA"
)

ruta_figuras <- file.path(ruta_salida, "figuras")
ruta_tablas <- file.path(ruta_salida, "tablas")
ruta_datos <- file.path(ruta_salida, "datos_procesados")

if (!file.exists(ruta_entrada)) {
  stop("No existe el archivo I1:\n", ruta_entrada)
}

dir.create(ruta_figuras, recursive = TRUE, showWarnings = FALSE)
dir.create(ruta_tablas, recursive = TRUE, showWarnings = FALSE)
dir.create(ruta_datos, recursive = TRUE, showWarnings = FALSE)

i1 <- read_csv(
  ruta_entrada,
  show_col_types = FALSE,
  trim_ws = TRUE
)

names(i1) <- trimws(names(i1))

columnas_requeridas <- c(
  "Nodo", "X", "Y", "N", "P", "K", "pH", "T", "Hv", "CE"
)

faltantes <- setdiff(columnas_requeridas, names(i1))

if (length(faltantes) > 0) {
  stop(
    "Faltan columnas en I1: ",
    paste(faltantes, collapse = ", "),
    "\nColumnas encontradas: ",
    paste(names(i1), collapse = ", ")
  )
}

i1 <- i1 |>
  select(all_of(columnas_requeridas)) |>
  mutate(
    Nodo = as.integer(Nodo),
    across(
      c(X, Y, N, P, K, pH, T, Hv, CE),
      as.numeric
    )
  )

variables <- c("N", "P", "K", "pH", "T", "Hv", "CE")

tabla_estructura <- tibble(
  Indicador = c(
    "Filas",
    "Columnas",
    "Nodos únicos",
    "Coordenadas únicas",
    "Duplicados de nodo",
    "Duplicados XY",
    "Valores faltantes"
  ),
  Valor = c(
    nrow(i1),
    ncol(i1),
    n_distinct(i1$Nodo),
    n_distinct(paste(i1$X, i1$Y, sep = "_")),
    sum(duplicated(i1$Nodo)),
    sum(duplicated(i1 |> select(X, Y))),
    sum(is.na(i1))
  )
)

describir <- function(x, variable) {
  x <- x[is.finite(x)]

  tibble(
    Variable = variable,
    n = length(x),
    Media = mean(x),
    Mediana = median(x),
    DE = sd(x),
    CV_porcentaje = ifelse(
      mean(x) == 0,
      NA_real_,
      100 * sd(x) / abs(mean(x))
    ),
    Minimo = min(x),
    Q1 = quantile(x, 0.25, names = FALSE),
    Q3 = quantile(x, 0.75, names = FALSE),
    Maximo = max(x),
    Asimetria = e1071::skewness(x, type = 2),
    Curtosis_exceso = e1071::kurtosis(x, type = 2)
  )
}

tabla_descriptiva <- map_dfr(
  variables,
  function(v) describir(i1[[v]], v)
)

tabla_faltantes <- tibble(
  Variable = names(i1),
  NA_n = map_int(i1, ~ sum(is.na(.x))),
  NA_porcentaje = map_dbl(i1, ~ 100 * mean(is.na(.x)))
)

coords <- i1 |> select(X, Y) |> as.matrix()
matriz_dist <- as.matrix(dist(coords))
diag(matriz_dist) <- NA_real_

dist_nn <- apply(
  matriz_dist,
  1,
  min,
  na.rm = TRUE
)

tabla_distancias <- tibble(
  Nodo = i1$Nodo,
  Distancia_vecino_mas_cercano_m = dist_nn
)

matriz_cor_spearman <- cor(
  i1 |> select(all_of(variables)),
  use = "pairwise.complete.obs",
  method = "spearman"
)

tabla_cor_spearman <- matriz_cor_spearman |>
  as.data.frame() |>
  tibble::rownames_to_column(
    var = "Variable_1"
  )

datos_largos <- i1 |>
  select(Nodo, X, Y, all_of(variables)) |>
  pivot_longer(
    cols = all_of(variables),
    names_to = "Variable",
    values_to = "Valor"
  )

g_hist <- ggplot(
  datos_largos,
  aes(x = Valor)
) +
  geom_histogram(bins = 8, color = "white") +
  facet_wrap(~ Variable, scales = "free", ncol = 3) +
  labs(
    title = "Distribución de variables de suelo en I1",
    x = "Valor",
    y = "Frecuencia"
  ) +
  theme_minimal()

g_puntos <- ggplot(
  datos_largos,
  aes(x = X, y = Y, fill = Valor)
) +
  geom_point(shape = 21, size = 3) +
  facet_wrap(~ Variable, ncol = 3) +
  coord_equal() +
  labs(
    title = "Valores observados de suelo en I1",
    x = "X (m)",
    y = "Y (m)",
    fill = "Valor"
  ) +
  theme_minimal()

ggsave(
  file.path(ruta_figuras, "Distribuciones_I1.png"),
  g_hist,
  width = 12,
  height = 9,
  dpi = 400
)

ggsave(
  file.path(ruta_figuras, "Puntos_observados_I1.png"),
  g_puntos,
  width = 12,
  height = 9,
  dpi = 400
)

write_xlsx(
  list(
    Estructura = tabla_estructura,
    Descriptiva = tabla_descriptiva,
    Faltantes = tabla_faltantes,
    Distancias = tabla_distancias,
    Correlacion_Spearman = tabla_cor_spearman
  ),
  path = file.path(
    ruta_tablas,
    "Auditoria_compacta_I1.xlsx"
  )
)

write_csv(
  i1,
  file.path(
    ruta_datos,
    "Suelo_I1_auditado.csv"
  )
)

capture.output(
  sessionInfo(),
  file = file.path(
    ruta_salida,
    "sessionInfo_I1_A.R.txt"
  )
)

cat(
  "\nSCRIPT I1-A FINALIZADO\n",
  "Base auditada:\n",
  file.path(ruta_datos, "Suelo_I1_auditado.csv"),
  "\n"
)
