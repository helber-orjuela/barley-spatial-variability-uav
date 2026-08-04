
## ============================================================
## ANÁLISIS DE SENSIBILIDAD DE LA CALIBRACIÓN DE CE
## Comparación con todas las muestras y exclusión temporal
## de las muestras 2 y 4
## ============================================================

## 0. Paquetes -------------------------------------------------

paquetes <- c(
  "readxl",
  "writexl",
  "dplyr",
  "ggplot2",
  "DescTools",
  "boot",
  "patchwork",
  "tidyr"
)

instalar <- paquetes[
  !paquetes %in% rownames(installed.packages())
]

if (length(instalar) > 0) {
  install.packages(instalar)
}

library(readxl)
library(writexl)
library(dplyr)
library(ggplot2)
library(DescTools)
library(boot)
library(patchwork)


## 1. Rutas ----------------------------------------------------

ruta_base <- paste0(
  "D:/UDCA_2017_CLASE/46_ARTICULOS_2025/",
  "ARTICULO_CEBADA_DOCTORADO/AJUSTES_ARTÍCULO"
)

ruta_datos <- file.path(
  ruta_base,
  "Datos_sensor_lab_CE.xlsx"
)

ruta_salida <- file.path(
  ruta_base,
  "02_SENSOR_CE",
  "sensibilidad"
)

if (!dir.exists(ruta_base)) {
  stop("No existe la carpeta base:\n", ruta_base)
}

if (!file.exists(ruta_datos)) {
  stop("No existe el archivo:\n", ruta_datos)
}

dir.create(
  ruta_salida,
  recursive = TRUE,
  showWarnings = FALSE
)


## 2. Lectura --------------------------------------------------

ce <- read_excel(
  ruta_datos,
  sheet = "Hoja1"
) |>
  select(
    Muestra,
    CE_sensor_uS_cm,
    CE_lab_uS_cm
  ) |>
  mutate(
    Muestra = as.integer(Muestra),
    CE_sensor_uS_cm = as.numeric(CE_sensor_uS_cm),
    CE_lab_uS_cm = as.numeric(CE_lab_uS_cm)
  ) |>
  filter(
    !is.na(Muestra),
    !is.na(CE_sensor_uS_cm),
    !is.na(CE_lab_uS_cm)
  )

print(ce)


## 3. Funciones auxiliares -------------------------------------

calcular_ccc <- function(x, y) {

  resultado <- DescTools::CCC(
    x,
    y,
    ci = "z-transform"
  )

  c(
    CCC = unname(resultado$rho.c[1]),
    CCC_LI95 = unname(resultado$rho.c[2]),
    CCC_LS95 = unname(resultado$rho.c[3])
  )
}


loocv_ols <- function(datos) {

  n <- nrow(datos)
  pred <- numeric(n)

  for (i in seq_len(n)) {

    modelo_i <- lm(
      CE_lab_uS_cm ~ CE_sensor_uS_cm,
      data = datos[-i, ]
    )

    pred[i] <- predict(
      modelo_i,
      newdata = datos[i, , drop = FALSE]
    )
  }

  residuos <- datos$CE_lab_uS_cm - pred

  c(
    RMSE_LOOCV = sqrt(mean(residuos^2)),
    MAE_LOOCV = mean(abs(residuos))
  )
}


bootstrap_coeficientes <- function(datos, R = 2000) {

  boot_fun <- function(data, indices) {

    d <- data[indices, ]

    ajuste <- lm(
      CE_lab_uS_cm ~ CE_sensor_uS_cm,
      data = d
    )

    coef(ajuste)
  }

  set.seed(2026)

  boot_fit <- boot(
    data = datos,
    statistic = boot_fun,
    R = R
  )

  ci_intercepto <- boot.ci(
    boot_fit,
    type = "perc",
    index = 1
  )

  ci_pendiente <- boot.ci(
    boot_fit,
    type = "perc",
    index = 2
  )

  tibble(
    Parametro = c(
      "Intercepto",
      "Pendiente"
    ),
    Estimacion = c(
      boot_fit$t0[1],
      boot_fit$t0[2]
    ),
    LI_95 = c(
      ci_intercepto$percent[4],
      ci_pendiente$percent[4]
    ),
    LS_95 = c(
      ci_intercepto$percent[5],
      ci_pendiente$percent[5]
    )
  )
}


evaluar_escenario <- function(datos, nombre_escenario) {

  modelo <- lm(
    CE_lab_uS_cm ~ CE_sensor_uS_cm,
    data = datos
  )

  datos_eval <- datos |>
    mutate(
      CE_pred_OLS = predict(modelo),
      residuo_OLS = CE_lab_uS_cm - CE_pred_OLS,
      diferencia_sensor_lab =
        CE_sensor_uS_cm - CE_lab_uS_cm,
      error_absoluto =
        abs(diferencia_sensor_lab),
      error_cuadrado =
        diferencia_sensor_lab^2
    )

  ccc <- calcular_ccc(
    datos_eval$CE_sensor_uS_cm,
    datos_eval$CE_lab_uS_cm
  )

  loocv <- loocv_ols(datos_eval)

  resumen <- tibble(
    Escenario = nombre_escenario,
    n = nrow(datos_eval),

    Media_sensor =
      mean(datos_eval$CE_sensor_uS_cm),

    Media_laboratorio =
      mean(datos_eval$CE_lab_uS_cm),

    Sesgo_sensor_lab =
      mean(datos_eval$diferencia_sensor_lab),

    MAE_sin_calibrar =
      mean(datos_eval$error_absoluto),

    RMSE_sin_calibrar =
      sqrt(mean(datos_eval$error_cuadrado)),

    Pearson =
      cor(
        datos_eval$CE_sensor_uS_cm,
        datos_eval$CE_lab_uS_cm,
        method = "pearson"
      ),

    Spearman =
      cor(
        datos_eval$CE_sensor_uS_cm,
        datos_eval$CE_lab_uS_cm,
        method = "spearman"
      ),

    CCC =
      ccc["CCC"],

    CCC_LI95 =
      ccc["CCC_LI95"],

    CCC_LS95 =
      ccc["CCC_LS95"],

    Intercepto_OLS =
      unname(coef(modelo)[1]),

    Pendiente_OLS =
      unname(coef(modelo)[2]),

    RMSE_ajuste =
      sqrt(mean(datos_eval$residuo_OLS^2)),

    MAE_ajuste =
      mean(abs(datos_eval$residuo_OLS)),

    RMSE_LOOCV =
      loocv["RMSE_LOOCV"],

    MAE_LOOCV =
      loocv["MAE_LOOCV"]
  )

  bootstrap <- bootstrap_coeficientes(
    datos_eval,
    R = 2000
  ) |>
    mutate(
      Escenario = nombre_escenario,
      .before = 1
    )

  list(
    modelo = modelo,
    datos = datos_eval,
    resumen = resumen,
    bootstrap = bootstrap
  )
}


## 4. Escenarios -----------------------------------------------

escenarios <- list(

  Completo_28 =
    ce,

  Sin_muestra_2 =
    ce |>
      filter(Muestra != 2),

  Sin_muestra_4 =
    ce |>
      filter(Muestra != 4),

  Sin_muestras_2_y_4 =
    ce |>
      filter(!Muestra %in% c(2, 4))
)


resultados <- lapply(
  names(escenarios),
  function(nombre) {
    evaluar_escenario(
      escenarios[[nombre]],
      nombre
    )
  }
)

names(resultados) <- names(escenarios)


## 5. Tablas consolidadas --------------------------------------

tabla_sensibilidad <- bind_rows(
  lapply(
    resultados,
    function(x) x$resumen
  )
)

tabla_bootstrap <- bind_rows(
  lapply(
    resultados,
    function(x) x$bootstrap
  )
)

print(tabla_sensibilidad)
print(tabla_bootstrap)


## 6. Cambios relativos frente al escenario completo ----------

referencia <- tabla_sensibilidad |>
  filter(Escenario == "Completo_28")

# Extraer valores escalares numéricos
ref_intercepto <- referencia |>
  pull(Intercepto_OLS) |>
  as.numeric()

ref_pendiente <- referencia |>
  pull(Pendiente_OLS) |>
  as.numeric()

ref_rmse_loocv <- referencia |>
  pull(RMSE_LOOCV) |>
  as.numeric()

ref_mae_loocv <- referencia |>
  pull(MAE_LOOCV) |>
  as.numeric()

ref_pearson <- referencia |>
  pull(Pearson) |>
  as.numeric()

ref_spearman <- referencia |>
  pull(Spearman) |>
  as.numeric()

ref_ccc <- referencia |>
  pull(CCC) |>
  as.numeric()


tabla_cambios <- tabla_sensibilidad |>
  mutate(
    Intercepto_OLS = as.numeric(Intercepto_OLS),
    Pendiente_OLS = as.numeric(Pendiente_OLS),
    RMSE_LOOCV = as.numeric(RMSE_LOOCV),
    MAE_LOOCV = as.numeric(MAE_LOOCV),
    Pearson = as.numeric(Pearson),
    Spearman = as.numeric(Spearman),
    CCC = as.numeric(CCC),
    
    Cambio_intercepto =
      Intercepto_OLS - ref_intercepto,
    
    Cambio_pendiente =
      Pendiente_OLS - ref_pendiente,
    
    Cambio_RMSE_LOOCV =
      RMSE_LOOCV - ref_rmse_loocv,
    
    Cambio_MAE_LOOCV =
      MAE_LOOCV - ref_mae_loocv,
    
    Cambio_Pearson =
      Pearson - ref_pearson,
    
    Cambio_Spearman =
      Spearman - ref_spearman,
    
    Cambio_CCC =
      CCC - ref_ccc
  )

print(tabla_cambios)


## 7. Diagnóstico de estabilidad -------------------------------

umbral_cambio_pendiente <- 0.15
umbral_cambio_rmse <- 0.10

tabla_estabilidad <- tabla_cambios |>
  mutate(
    Cambio_rel_pendiente =
      abs(
        Cambio_pendiente /
          ref_pendiente
      ),
    
    Cambio_rel_RMSE_LOOCV =
      abs(
        Cambio_RMSE_LOOCV /
          ref_rmse_loocv
      ),
    
    Estabilidad_pendiente =
      ifelse(
        Cambio_rel_pendiente <=
          umbral_cambio_pendiente,
        "Estable",
        "Sensible"
      ),
    
    Estabilidad_RMSE =
      ifelse(
        Cambio_rel_RMSE_LOOCV <=
          umbral_cambio_rmse,
        "Estable",
        "Sensible"
      )
  )

print(tabla_estabilidad)

## 8. Figuras comparativas -------------------------------------

datos_grafico <- bind_rows(
  lapply(
    names(resultados),
    function(nombre) {

      resultados[[nombre]]$datos |>
        mutate(
          Escenario = nombre
        )
    }
  )
)


coeficientes <- tabla_sensibilidad |>
  select(
    Escenario,
    Intercepto_OLS,
    Pendiente_OLS
  )


g_escenarios <- ggplot(
  datos_grafico,
  aes(
    x = CE_sensor_uS_cm,
    y = CE_lab_uS_cm
  )
) +
  geom_point(
    size = 2
  ) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed"
  ) +
  geom_smooth(
    method = "lm",
    se = FALSE
  ) +
  facet_wrap(
    ~ Escenario,
    ncol = 2
  ) +
  labs(
    title =
      "Sensibilidad de la calibración OLS",

    subtitle =
      "Comparación con y sin las muestras 2 y 4",

    x =
      "CE del sensor (µS/cm)",

    y =
      "CE de laboratorio (µS/cm)"
  ) +
  theme_minimal()


tabla_coef_larga <- tabla_sensibilidad |>
  select(
    Escenario,
    Intercepto_OLS,
    Pendiente_OLS
  ) |>
  tidyr::pivot_longer(
    cols = c(
      Intercepto_OLS,
      Pendiente_OLS
    ),
    names_to = "Parametro",
    values_to = "Valor"
  )


g_coeficientes <- ggplot(
  tabla_coef_larga,
  aes(
    x = Escenario,
    y = Valor,
    group = Parametro
  )
) +
  geom_point(
    size = 2.5
  ) +
  geom_line() +
  facet_wrap(
    ~ Parametro,
    scales = "free_y",
    ncol = 1
  ) +
  labs(
    title =
      "Estabilidad de los coeficientes OLS",

    x = NULL,
    y = "Valor"
  ) +
  theme_minimal() +
  theme(
    axis.text.x =
      element_text(
        angle = 30,
        hjust = 1
      )
  )


tabla_error_larga <- tabla_sensibilidad |>
  select(
    Escenario,
    RMSE_LOOCV,
    MAE_LOOCV
  ) |>
  tidyr::pivot_longer(
    cols = c(
      RMSE_LOOCV,
      MAE_LOOCV
    ),
    names_to = "Metrica",
    values_to = "Valor"
  )


g_error <- ggplot(
  tabla_error_larga,
  aes(
    x = Escenario,
    y = Valor,
    group = Metrica
  )
) +
  geom_point(
    size = 2.5
  ) +
  geom_line() +
  facet_wrap(
    ~ Metrica,
    scales = "free_y",
    ncol = 1
  ) +
  labs(
    title =
      "Sensibilidad del error de validación cruzada",

    x = NULL,
    y = "µS/cm"
  ) +
  theme_minimal() +
  theme(
    axis.text.x =
      element_text(
        angle = 30,
        hjust = 1
      )
  )


panel_sensibilidad <- (
  g_escenarios
) / (
  g_coeficientes |
    g_error
)

print(panel_sensibilidad)


## 9. Exportación ----------------------------------------------

ggsave(
  filename = file.path(
    ruta_salida,
    "Figura_sensibilidad_calibracion_CE.png"
  ),
  plot = panel_sensibilidad,
  width = 13,
  height = 12,
  dpi = 400
)

ggsave(
  filename = file.path(
    ruta_salida,
    "Figura_sensibilidad_calibracion_CE.tiff"
  ),
  plot = panel_sensibilidad,
  width = 13,
  height = 12,
  dpi = 600,
  compression = "lzw"
)


write_xlsx(
  list(
    Resumen_sensibilidad =
      tabla_sensibilidad,

    Cambios_relativos =
      tabla_cambios,

    Diagnostico_estabilidad =
      tabla_estabilidad,

    Bootstrap =
      tabla_bootstrap,

    Datos_completos =
      resultados$Completo_28$datos,

    Datos_sin_muestra_2 =
      resultados$Sin_muestra_2$datos,

    Datos_sin_muestra_4 =
      resultados$Sin_muestra_4$datos,

    Datos_sin_muestras_2_y_4 =
      resultados$Sin_muestras_2_y_4$datos
  ),

  path = file.path(
    ruta_salida,
    "Resultados_sensibilidad_calibracion_CE.xlsx"
  )
)


## 10. Ecuaciones por escenario --------------------------------

for (nombre in names(resultados)) {

  modelo <- resultados[[nombre]]$modelo

  cat(
    "\n------------------------------------\n"
  )

  cat(
    "Escenario: ",
    nombre,
    "\n",
    sep = ""
  )

  cat(
    "CE_lab = ",
    round(
      coef(modelo)[1],
      4
    ),
    " + ",
    round(
      coef(modelo)[2],
      4
    ),
    " × CE_sensor\n",
    sep = ""
  )
}


cat(
  "\nResultados guardados en:\n",
  ruta_salida,
  "\n"
)

cat(
  "\nEste análisis no elimina automáticamente las muestras 2 y 4. ",
  "Evalúa cuánto cambia la calibración cuando se excluyen temporalmente.\n"
)
