
## ============================================================
## I1 – SCRIPT B. SEMIVARIOGRAMAS Y VALIDACIÓN COMPACTA
## ============================================================
## Entrada:
##   AJUSTES_ARTÍCULO/06_I1_COMPACTO/01_AUDITORIA/
##   datos_procesados/Suelo_I1_auditado.csv
##
## Salida:
##   AJUSTES_ARTÍCULO/06_I1_COMPACTO/02_SEMIVARIOGRAMAS
##
## Propósito:
##   Obtener evidencia geoestadística temporal complementaria de I1.
##   Se ajustan modelos Sph, Exp, Gau y Bes y se selecciona usando
##   RMSE, MAE, MSNE y consistencia física.
## ============================================================

paquetes <- c(
  "sf", "gstat", "readr", "dplyr", "tidyr",
  "purrr", "ggplot2", "writexl", "patchwork"
)

instalar <- paquetes[!paquetes %in% rownames(installed.packages())]
if (length(instalar) > 0) install.packages(instalar)

library(sf)
library(gstat)
library(readr)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(writexl)
library(patchwork)

ruta_base <- paste0(
  "D:/UDCA_2017_CLASE/46_ARTICULOS_2025/",
  "ARTICULO_CEBADA_DOCTORADO/AJUSTES_ARTÍCULO"
)

ruta_entrada <- file.path(
  ruta_base,
  "06_I1_COMPACTO",
  "01_AUDITORIA",
  "datos_procesados",
  "Suelo_I1_auditado.csv"
)

ruta_lote <- file.path(
  ruta_base,
  "DATOS BASE",
  "Lote_2.shp"
)

ruta_salida <- file.path(
  ruta_base,
  "06_I1_COMPACTO",
  "02_SEMIVARIOGRAMAS"
)

ruta_tablas <- file.path(ruta_salida, "tablas")
ruta_figuras <- file.path(ruta_salida, "figuras")
ruta_objetos <- file.path(ruta_salida, "objetos_R")

if (!file.exists(ruta_entrada)) {
  stop("No existe la base auditada I1:\n", ruta_entrada)
}

if (!file.exists(ruta_lote)) {
  stop("No existe Lote_2.shp:\n", ruta_lote)
}

dir.create(ruta_tablas, recursive = TRUE, showWarnings = FALSE)
dir.create(ruta_figuras, recursive = TRUE, showWarnings = FALSE)
dir.create(ruta_objetos, recursive = TRUE, showWarnings = FALSE)

i1 <- read_csv(
  ruta_entrada,
  show_col_types = FALSE
)

crs_objetivo <- 32618

i1_sf <- st_as_sf(
  i1,
  coords = c("X", "Y"),
  crs = crs_objetivo,
  remove = FALSE
)

lote <- st_read(ruta_lote, quiet = TRUE)

if (is.na(st_crs(lote))) {
  st_crs(lote) <- crs_objetivo
} else {
  lote <- st_transform(lote, crs_objetivo)
}

variables <- c("N", "P", "K", "pH", "T", "Hv", "CE")
modelos_candidatos <- c("Sph", "Exp", "Gau", "Bes")

coords <- st_coordinates(i1_sf)
dist_max <- max(as.vector(dist(coords)), na.rm = TRUE)
cutoff <- 0.50 * dist_max
numero_lags <- 7
width <- cutoff / numero_lags

extraer_parametros <- function(modelo) {
  nugget <- modelo |>
    filter(model == "Nug") |>
    pull(psill)

  if (length(nugget) == 0) nugget <- 0

  estructura <- modelo |>
    filter(model != "Nug")

  if (nrow(estructura) == 0) {
    return(
      tibble(
        Pepita = nugget[1],
        Psill = NA_real_,
        Silla = NA_real_,
        Rango = NA_real_,
        GDE = NA_real_
      )
    )
  }

  psill <- estructura$psill[1]
  silla <- nugget[1] + psill

  tibble(
    Pepita = nugget[1],
    Psill = psill,
    Silla = silla,
    Rango = estructura$range[1],
    GDE = ifelse(silla > 0, psill / silla, NA_real_)
  )
}

calcular_cv <- function(variable, datos_sf, modelo) {
  fml <- as.formula(paste(variable, "~ 1"))

  cv <- try(
    krige.cv(
      fml,
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
        metricas = tibble(
          ME = NA_real_,
          MAE = NA_real_,
          RMSE = NA_real_,
          RMSSE = NA_real_,
          MSNE = NA_real_,
          Spearman = NA_real_
        ),
        datos = NULL
      )
    )
  }

  residuo <- cv$observed - cv$var1.pred
  valida_var <- is.finite(cv$var1.var) & cv$var1.var > 0
  residuo_est <- rep(NA_real_, length(residuo))

  residuo_est[valida_var] <-
    residuo[valida_var] / sqrt(cv$var1.var[valida_var])

  list(
    metricas = tibble(
      ME = mean(residuo, na.rm = TRUE),
      MAE = mean(abs(residuo), na.rm = TRUE),
      RMSE = sqrt(mean(residuo^2, na.rm = TRUE)),
      RMSSE = ifelse(
        any(valida_var),
        sqrt(mean(residuo_est^2, na.rm = TRUE)),
        NA_real_
      ),
      MSNE = ifelse(
        any(valida_var),
        mean(residuo_est^2, na.rm = TRUE),
        NA_real_
      ),
      Spearman = suppressWarnings(
        cor(
          cv$observed,
          cv$var1.pred,
          method = "spearman",
          use = "complete.obs"
        )
      )
    ),
    datos = cv
  )
}

ajustar_variable <- function(variable) {
  fml <- as.formula(paste(variable, "~ 1"))

  vexp <- variogram(
    fml,
    i1_sf,
    cutoff = cutoff,
    width = width,
    cressie = FALSE
  )

  var_datos <- var(i1_sf[[variable]], na.rm = TRUE)
  gamma_min <- min(vexp$gamma, na.rm = TRUE)
  nugget_ini <- max(0, min(gamma_min, 0.25 * var_datos))
  psill_ini <- max(var_datos - nugget_ini, 0.25 * var_datos)
  rango_ini <- max(cutoff / 3, width)

  resultados <- map(
    modelos_candidatos,
    function(codigo) {
      inicial <- vgm(
        psill = psill_ini,
        model = codigo,
        range = rango_ini,
        nugget = nugget_ini
      )

      ajuste <- try(
        fit.variogram(
          vexp,
          inicial,
          fit.method = 6
        ),
        silent = TRUE
      )

      if (
        inherits(ajuste, "try-error") ||
        isTRUE(attr(ajuste, "singular"))
      ) {
        return(
          list(
            resumen = tibble(
              Variable = variable,
              Modelo = codigo,
              Valido = FALSE,
              Motivo = "Error o ajuste singular",
              SSErr = NA_real_,
              Pepita = NA_real_,
              Psill = NA_real_,
              Silla = NA_real_,
              Rango = NA_real_,
              GDE = NA_real_,
              ME = NA_real_,
              MAE = NA_real_,
              RMSE = NA_real_,
              RMSSE = NA_real_,
              MSNE = NA_real_,
              Spearman = NA_real_
            ),
            modelo = NULL,
            cv = NULL
          )
        )
      }

      pars <- extraer_parametros(ajuste)
      cv <- calcular_cv(variable, i1_sf, ajuste)

      valido <- is.finite(pars$Rango) &&
        pars$Rango > 0 &&
        is.finite(pars$Silla) &&
        pars$Silla > 0 &&
        pars$Rango <= 3 * cutoff

      resumen <- bind_cols(
        tibble(
          Variable = variable,
          Modelo = codigo,
          Valido = valido,
          Motivo = ifelse(
            valido,
            "",
            "Parámetros fuera del dominio aceptado"
          ),
          SSErr = as.numeric(attr(ajuste, "SSErr"))
        ),
        pars,
        cv$metricas
      )

      list(
        resumen = resumen,
        modelo = ajuste,
        cv = cv$datos
      )
    }
  )

  tabla <- bind_rows(map(resultados, "resumen"))

  candidatos <- tabla |>
    filter(
      Valido,
      is.finite(RMSE),
      is.finite(MSNE),
      is.finite(SSErr)
    ) |>
    mutate(
      Rank_RMSE = min_rank(RMSE),
      Rank_MAE = min_rank(MAE),
      Rank_MSNE = min_rank(abs(MSNE - 1)),
      Rank_ME = min_rank(abs(ME)),
      Rank_SSErr = min_rank(SSErr),
      Score =
        0.35 * Rank_RMSE +
        0.20 * Rank_MAE +
        0.20 * Rank_MSNE +
        0.10 * Rank_ME +
        0.15 * Rank_SSErr
    ) |>
    arrange(Score, RMSE, abs(MSNE - 1))

  if (nrow(candidatos) > 0) {
    seleccionado <- candidatos |> slice(1)
    indice <- match(seleccionado$Modelo, modelos_candidatos)
    modelo_sel <- resultados[[indice]]$modelo
    cv_sel <- resultados[[indice]]$cv
  } else {
    seleccionado <- tibble()
    modelo_sel <- NULL
    cv_sel <- NULL
  }

  list(
    variograma = vexp,
    tabla = tabla,
    seleccionado = seleccionado,
    modelo = modelo_sel,
    cv = cv_sel
  )
}

resultados <- map(
  variables,
  ajustar_variable
)

names(resultados) <- variables

tabla_modelos <- bind_rows(
  map(resultados, "tabla")
)

tabla_seleccionados <- bind_rows(
  map(resultados, "seleccionado")
)

tabla_pares <- imap_dfr(
  resultados,
  function(res, variable) {
    as_tibble(res$variograma) |>
      transmute(
        Variable = variable,
        Intervalo = row_number(),
        Numero_pares = np,
        Distancia_media_m = dist,
        Semivarianza = gamma
      )
  }
)

graficos <- list()

for (variable in variables) {
  res <- resultados[[variable]]
  vexp <- res$variograma

  g <- ggplot(
    as_tibble(vexp),
    aes(x = dist, y = gamma)
  ) +
    geom_point(aes(size = np)) +
    labs(
      title = variable,
      x = "Distancia (m)",
      y = "Semivarianza",
      size = "Pares"
    ) +
    theme_minimal()

  if (!is.null(res$modelo)) {
    curva <- variogramLine(
      res$modelo,
      maxdist = max(vexp$dist, na.rm = TRUE),
      n = 250
    )

    g <- g +
      geom_line(
        data = curva,
        aes(x = dist, y = gamma),
        inherit.aes = FALSE
      ) +
      labs(
        subtitle = paste0(
          res$seleccionado$Modelo,
          " | Rango=",
          round(res$seleccionado$Rango, 2),
          " m | RMSE=",
          round(res$seleccionado$RMSE, 3),
          " | MSNE=",
          round(res$seleccionado$MSNE, 3)
        )
      )
  } else {
    g <- g +
      labs(subtitle = "Sin modelo válido")
  }

  graficos[[variable]] <- g
}

panel <- wrap_plots(
  graficos,
  ncol = 2
) +
  plot_annotation(
    title = "Semivariogramas omnidireccionales compactos de I1"
  )

ggsave(
  file.path(
    ruta_figuras,
    "Semivariogramas_compactos_I1.png"
  ),
  panel,
  width = 14,
  height = 16,
  dpi = 400
)

write_xlsx(
  list(
    Parametros = tibble(
      Parametro = c(
        "Número de nodos",
        "Distancia máxima (m)",
        "Cutoff (m)",
        "Número de intervalos",
        "Ancho de intervalo (m)"
      ),
      Valor = c(
        nrow(i1_sf),
        dist_max,
        cutoff,
        numero_lags,
        width
      )
    ),
    Pares_por_intervalo = tabla_pares,
    Todos_modelos = tabla_modelos,
    Modelos_seleccionados = tabla_seleccionados
  ),
  path = file.path(
    ruta_tablas,
    "Resultados_semivariogramas_I1.xlsx"
  )
)

modelos_finales <- map(
  resultados,
  "modelo"
)

validaciones <- map(
  resultados,
  "cv"
)

saveRDS(
  list(
    i1_sf = i1_sf,
    lote = lote,
    variables = variables,
    modelos_finales = modelos_finales,
    validaciones = validaciones,
    tabla_modelos = tabla_modelos,
    tabla_seleccionados = tabla_seleccionados,
    cutoff = cutoff,
    width = width
  ),
  file = file.path(
    ruta_objetos,
    "objetos_semivariogramas_I1.rds"
  )
)

capture.output(
  sessionInfo(),
  file = file.path(
    ruta_salida,
    "sessionInfo_I1_B.R.txt"
  )
)

cat(
  "\nSCRIPT I1-B FINALIZADO\n",
  "Modelos seleccionados:\n"
)

print(
  tabla_seleccionados |>
    select(
      Variable, Modelo, Pepita, Silla,
      Rango, GDE, RMSE, MSNE, Spearman
    )
)
