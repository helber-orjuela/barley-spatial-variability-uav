## ============================================================
## CALIBRACIÓN Y CONCORDANCIA DE CONDUCTIVIDAD ELÉCTRICA
## Sensor OEM frente a método de laboratorio en pasta saturada
## ============================================================

## 0. Paquetes -------------------------------------------------

paquetes <- c(
  "readxl",
  "writexl",
  "dplyr",
  "ggplot2",
  "DescTools",
  "mcr",
  "boot",
  "patchwork"
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
library(mcr)
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
  "02_SENSOR_CE"
)

# Verificar carpeta base
if (!dir.exists(ruta_base)) {
  stop(
    "No existe la carpeta base:\n",
    ruta_base
  )
}

# Verificar archivo Excel
if (!file.exists(ruta_datos)) {
  stop(
    "No existe el archivo:\n",
    ruta_datos,
    "\n\nArchivos encontrados en la carpeta:\n",
    paste(
      list.files(ruta_base),
      collapse = "\n"
    )
  )
}

# Crear carpeta de salida
dir.create(
  ruta_salida,
  recursive = TRUE,
  showWarnings = FALSE
)

cat(
  "\nCarpeta base:\n",
  ruta_base,
  "\n\nArchivo de entrada:\n",
  ruta_datos,
  "\n\nCarpeta de salida:\n",
  ruta_salida,
  "\n\n"
)


## 2. Lectura de datos -----------------------------------------

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
summary(ce)

cat(
  "\nNúmero de pares válidos:",
  nrow(ce),
  "\n"
)


## 3. Verificaciones básicas ----------------------------------

if (nrow(ce) < 10) {
  stop(
    "El número de pares es demasiado bajo para realizar ",
    "una calibración exploratoria."
  )
}

if (any(duplicated(ce$Muestra))) {
  warning(
    "Existen identificadores de muestra duplicados."
  )
}

if (any(ce$CE_sensor_uS_cm <= 0)) {
  warning(
    "Existen valores de CE del sensor menores o iguales a cero."
  )
}

if (any(ce$CE_lab_uS_cm <= 0)) {
  warning(
    "Existen valores de CE de laboratorio menores o iguales a cero."
  )
}


## 4. Estadística descriptiva ---------------------------------

tabla_descriptiva <- data.frame(
  Metodo = c(
    "Sensor",
    "Laboratorio"
  ),
  n = c(
    nrow(ce),
    nrow(ce)
  ),
  Media = c(
    mean(ce$CE_sensor_uS_cm),
    mean(ce$CE_lab_uS_cm)
  ),
  Mediana = c(
    median(ce$CE_sensor_uS_cm),
    median(ce$CE_lab_uS_cm)
  ),
  Desviacion_estandar = c(
    sd(ce$CE_sensor_uS_cm),
    sd(ce$CE_lab_uS_cm)
  ),
  Minimo = c(
    min(ce$CE_sensor_uS_cm),
    min(ce$CE_lab_uS_cm)
  ),
  Maximo = c(
    max(ce$CE_sensor_uS_cm),
    max(ce$CE_lab_uS_cm)
  ),
  CV_porcentaje = c(
    100 *
      sd(ce$CE_sensor_uS_cm) /
      mean(ce$CE_sensor_uS_cm),
    100 *
      sd(ce$CE_lab_uS_cm) /
      mean(ce$CE_lab_uS_cm)
  )
)

print(tabla_descriptiva)


## 5. Cálculo de errores ---------------------------------------

ce_eval <- ce |>
  mutate(
    diferencia_sensor_lab =
      CE_sensor_uS_cm -
      CE_lab_uS_cm,
    
    diferencia_lab_sensor =
      CE_lab_uS_cm -
      CE_sensor_uS_cm,
    
    error_absoluto =
      abs(diferencia_sensor_lab),
    
    error_cuadrado =
      diferencia_sensor_lab^2,
    
    media_metodos =
      (
        CE_sensor_uS_cm +
          CE_lab_uS_cm
      ) / 2,
    
    error_relativo_porcentaje =
      100 *
      error_absoluto /
      CE_lab_uS_cm
  )

n <- nrow(ce_eval)

sesgo <- mean(
  ce_eval$diferencia_sensor_lab
)

mae <- mean(
  ce_eval$error_absoluto
)

rmse <- sqrt(
  mean(
    ce_eval$error_cuadrado
  )
)

mape <- mean(
  ce_eval$error_relativo_porcentaje
)

pearson <- cor(
  ce_eval$CE_sensor_uS_cm,
  ce_eval$CE_lab_uS_cm,
  method = "pearson"
)

spearman <- cor(
  ce_eval$CE_sensor_uS_cm,
  ce_eval$CE_lab_uS_cm,
  method = "spearman"
)


## 6. Concordancia de Lin --------------------------------------

ccc <- DescTools::CCC(
  ce_eval$CE_sensor_uS_cm,
  ce_eval$CE_lab_uS_cm,
  ci = "z-transform"
)

ccc_valor <- unname(
  ccc$rho.c[1]
)

ccc_li <- unname(
  ccc$rho.c[2]
)

ccc_ls <- unname(
  ccc$rho.c[3]
)


## 7. Tabla de métricas originales -----------------------------

tabla_metricas <- data.frame(
  Indicador = c(
    "Número de pares",
    "Media sensor (µS/cm)",
    "Media laboratorio (µS/cm)",
    "Sesgo sensor - laboratorio (µS/cm)",
    "MAE (µS/cm)",
    "RMSE (µS/cm)",
    "MAPE (%)",
    "Correlación de Pearson",
    "Correlación de Spearman",
    "Coeficiente de concordancia de Lin",
    "Lin CCC - límite inferior 95%",
    "Lin CCC - límite superior 95%"
  ),
  Valor = c(
    n,
    mean(ce_eval$CE_sensor_uS_cm),
    mean(ce_eval$CE_lab_uS_cm),
    sesgo,
    mae,
    rmse,
    mape,
    pearson,
    spearman,
    ccc_valor,
    ccc_li,
    ccc_ls
  )
)

print(tabla_metricas)


## 8. Regresión lineal ordinaria OLS ----------------------------

modelo_ols <- lm(
  CE_lab_uS_cm ~ CE_sensor_uS_cm,
  data = ce_eval
)

print(
  summary(modelo_ols)
)

intercepto_ols <- unname(
  coef(modelo_ols)[1]
)

pendiente_ols <- unname(
  coef(modelo_ols)[2]
)

ce_eval <- ce_eval |>
  mutate(
    CE_pred_OLS =
      predict(modelo_ols),
    
    residuo_OLS =
      CE_lab_uS_cm -
      CE_pred_OLS
  )

rmse_ols <- sqrt(
  mean(
    ce_eval$residuo_OLS^2
  )
)

mae_ols <- mean(
  abs(
    ce_eval$residuo_OLS
  )
)


## 9. Regresión de Deming --------------------------------------

modelo_deming <- mcreg(
  x = ce_eval$CE_sensor_uS_cm,
  y = ce_eval$CE_lab_uS_cm,
  method.reg = "Deming",
  error.ratio = 1,
  method.ci = "analytical"
)

printSummary(
  modelo_deming
)

coef_deming <- getCoefficients(
  modelo_deming
)

intercepto_deming <- coef_deming[
  1,
  "EST"
]

pendiente_deming <- coef_deming[
  2,
  "EST"
]

ce_eval <- ce_eval |>
  mutate(
    CE_pred_Deming =
      intercepto_deming +
      pendiente_deming *
      CE_sensor_uS_cm,
    
    residuo_Deming =
      CE_lab_uS_cm -
      CE_pred_Deming
  )

rmse_deming <- sqrt(
  mean(
    ce_eval$residuo_Deming^2
  )
)

mae_deming <- mean(
  abs(
    ce_eval$residuo_Deming
  )
)


## 10. Regresión Passing-Bablok --------------------------------

modelo_pb <- mcreg(
  x = ce_eval$CE_sensor_uS_cm,
  y = ce_eval$CE_lab_uS_cm,
  method.reg = "PaBa",
  method.ci = "analytical"
)

printSummary(
  modelo_pb
)

coef_pb <- getCoefficients(
  modelo_pb
)

intercepto_pb <- coef_pb[
  1,
  "EST"
]

pendiente_pb <- coef_pb[
  2,
  "EST"
]

ce_eval <- ce_eval |>
  mutate(
    CE_pred_Passing_Bablok =
      intercepto_pb +
      pendiente_pb *
      CE_sensor_uS_cm,
    
    residuo_Passing_Bablok =
      CE_lab_uS_cm -
      CE_pred_Passing_Bablok
  )

rmse_pb <- sqrt(
  mean(
    ce_eval$residuo_Passing_Bablok^2
  )
)

mae_pb <- mean(
  abs(
    ce_eval$residuo_Passing_Bablok
  )
)


## 11. Comparación de modelos ---------------------------------

tabla_modelos <- data.frame(
  Modelo = c(
    "Sin calibración",
    "OLS",
    "Deming",
    "Passing-Bablok"
  ),
  
  Intercepto = c(
    0,
    intercepto_ols,
    intercepto_deming,
    intercepto_pb
  ),
  
  Pendiente = c(
    1,
    pendiente_ols,
    pendiente_deming,
    pendiente_pb
  ),
  
  RMSE = c(
    rmse,
    rmse_ols,
    rmse_deming,
    rmse_pb
  ),
  
  MAE = c(
    mae,
    mae_ols,
    mae_deming,
    mae_pb
  ),
  
  Mejora_RMSE_porcentaje = c(
    0,
    100 *
      (rmse - rmse_ols) /
      rmse,
    100 *
      (rmse - rmse_deming) /
      rmse,
    100 *
      (rmse - rmse_pb) /
      rmse
  )
)

print(tabla_modelos)


## 12. Bland–Altman --------------------------------------------

sesgo_bland <- mean(
  ce_eval$diferencia_lab_sensor
)

sd_bland <- sd(
  ce_eval$diferencia_lab_sensor
)

limite_inferior_bland <-
  sesgo_bland -
  1.96 *
  sd_bland

limite_superior_bland <-
  sesgo_bland +
  1.96 *
  sd_bland

tabla_bland_altman <- data.frame(
  Sesgo_laboratorio_sensor =
    sesgo_bland,
  
  Limite_inferior_95 =
    limite_inferior_bland,
  
  Limite_superior_95 =
    limite_superior_bland
)

print(
  tabla_bland_altman
)


## 13. Validación cruzada LOOCV para OLS -----------------------

pred_loo_ols <- numeric(n)

for (i in seq_len(n)) {
  
  ajuste_i <- lm(
    CE_lab_uS_cm ~ CE_sensor_uS_cm,
    data = ce_eval[-i, ]
  )
  
  pred_loo_ols[i] <- predict(
    ajuste_i,
    newdata = ce_eval[
      i,
      ,
      drop = FALSE
    ]
  )
}

ce_eval <- ce_eval |>
  mutate(
    CE_pred_LOOCV_OLS =
      pred_loo_ols,
    
    residuo_LOOCV_OLS =
      CE_lab_uS_cm -
      CE_pred_LOOCV_OLS
  )

rmse_loocv_ols <- sqrt(
  mean(
    ce_eval$residuo_LOOCV_OLS^2
  )
)

mae_loocv_ols <- mean(
  abs(
    ce_eval$residuo_LOOCV_OLS
  )
)

tabla_loocv <- data.frame(
  Modelo =
    "OLS con validación cruzada LOOCV",
  
  RMSE_LOOCV =
    rmse_loocv_ols,
  
  MAE_LOOCV =
    mae_loocv_ols
)

print(
  tabla_loocv
)


## 14. Bootstrap de coeficientes OLS ---------------------------

boot_fun <- function(
    data,
    indices
) {
  
  d <- data[
    indices,
  ]
  
  fit <- lm(
    CE_lab_uS_cm ~ CE_sensor_uS_cm,
    data = d
  )
  
  coef(fit)
}

set.seed(2026)

boot_ols <- boot(
  data = ce_eval,
  statistic = boot_fun,
  R = 2000
)

ci_intercepto <- boot.ci(
  boot_ols,
  type = "perc",
  index = 1
)

ci_pendiente <- boot.ci(
  boot_ols,
  type = "perc",
  index = 2
)

tabla_bootstrap <- data.frame(
  Parametro = c(
    "Intercepto OLS",
    "Pendiente OLS"
  ),
  
  Estimacion = c(
    intercepto_ols,
    pendiente_ols
  ),
  
  Limite_inferior_95 = c(
    ci_intercepto$percent[4],
    ci_pendiente$percent[4]
  ),
  
  Limite_superior_95 = c(
    ci_intercepto$percent[5],
    ci_pendiente$percent[5]
  )
)

print(
  tabla_bootstrap
)


## 15. Detección exploratoria de observaciones influyentes -----

diagnostico_ols <- data.frame(
  Muestra =
    ce_eval$Muestra,
  
  Residuo_estudentizado =
    rstudent(
      modelo_ols
    ),
  
  Distancia_Cook =
    cooks.distance(
      modelo_ols
    ),
  
  Leverage =
    hatvalues(
      modelo_ols
    )
)

diagnostico_ols <- diagnostico_ols |>
  mutate(
    Influyente_Cook =
      Distancia_Cook >
      4 / n,
    
    Residuo_extremo =
      abs(
        Residuo_estudentizado
      ) > 2
  )

print(
  diagnostico_ols
)


## 16. Figuras -------------------------------------------------

g_concordancia <- ggplot(
  ce_eval,
  aes(
    x = CE_sensor_uS_cm,
    y = CE_lab_uS_cm
  )
) +
  geom_point(
    size = 2.5
  ) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed"
  ) +
  geom_smooth(
    method = "lm",
    se = TRUE
  ) +
  labs(
    title =
      "Comparación de conductividad eléctrica",
    
    subtitle =
      "Sensor proximal frente a laboratorio en pasta saturada",
    
    x =
      "CE del sensor (µS/cm)",
    
    y =
      "CE de laboratorio (µS/cm)"
  ) +
  theme_minimal()


g_bland <- ggplot(
  ce_eval,
  aes(
    x = media_metodos,
    y = diferencia_lab_sensor
  )
) +
  geom_point(
    size = 2.5
  ) +
  geom_hline(
    yintercept =
      sesgo_bland
  ) +
  geom_hline(
    yintercept =
      c(
        limite_inferior_bland,
        limite_superior_bland
      ),
    linetype =
      "dashed"
  ) +
  labs(
    title =
      "Análisis de Bland–Altman",
    
    x =
      "Media de los dos métodos (µS/cm)",
    
    y =
      "Laboratorio – sensor (µS/cm)"
  ) +
  theme_minimal()


g_residuos <- ggplot(
  ce_eval,
  aes(
    x = CE_pred_Deming,
    y = residuo_Deming
  )
) +
  geom_point(
    size = 2.5
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  labs(
    title =
      "Residuos del modelo de Deming",
    
    x =
      "CE estimada mediante Deming (µS/cm)",
    
    y =
      "Residuo (µS/cm)"
  ) +
  theme_minimal()


g_calibracion <- ggplot(
  ce_eval,
  aes(
    x = CE_sensor_uS_cm,
    y = CE_lab_uS_cm
  )
) +
  geom_point(
    size = 2.5
  ) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed"
  ) +
  geom_abline(
    intercept =
      intercepto_ols,
    
    slope =
      pendiente_ols,
    
    linetype =
      "solid"
  ) +
  geom_abline(
    intercept =
      intercepto_deming,
    
    slope =
      pendiente_deming,
    
    linetype =
      "dotdash"
  ) +
  geom_abline(
    intercept =
      intercepto_pb,
    
    slope =
      pendiente_pb,
    
    linetype =
      "longdash"
  ) +
  labs(
    title =
      "Modelos candidatos de calibración de CE",
    
    subtitle =
      "1:1, OLS, Deming y Passing–Bablok",
    
    x =
      "CE del sensor (µS/cm)",
    
    y =
      "CE de laboratorio (µS/cm)"
  ) +
  theme_minimal()


panel_calibracion <- (
  g_concordancia |
    g_bland
) / (
  g_residuos |
    g_calibracion
)

print(
  panel_calibracion
)


## 17. Exportación de figuras ---------------------------------

ggsave(
  filename = file.path(
    ruta_salida,
    "Figura_panel_calibracion_CE.png"
  ),
  plot = panel_calibracion,
  width = 13,
  height = 10,
  dpi = 400
)

ggsave(
  filename = file.path(
    ruta_salida,
    "Figura_panel_calibracion_CE.tiff"
  ),
  plot = panel_calibracion,
  width = 13,
  height = 10,
  dpi = 600,
  compression = "lzw"
)

ggsave(
  filename = file.path(
    ruta_salida,
    "Figura_concordancia_CE.png"
  ),
  plot = g_concordancia,
  width = 7,
  height = 6,
  dpi = 400
)

ggsave(
  filename = file.path(
    ruta_salida,
    "Figura_Bland_Altman_CE.png"
  ),
  plot = g_bland,
  width = 7,
  height = 6,
  dpi = 400
)

ggsave(
  filename = file.path(
    ruta_salida,
    "Figura_residuos_Deming_CE.png"
  ),
  plot = g_residuos,
  width = 7,
  height = 6,
  dpi = 400
)

ggsave(
  filename = file.path(
    ruta_salida,
    "Figura_modelos_calibracion_CE.png"
  ),
  plot = g_calibracion,
  width = 7,
  height = 6,
  dpi = 400
)


## 18. Exportación de resultados -------------------------------

write_xlsx(
  list(
    Datos_originales =
      ce,
    
    Datos_calculados =
      ce_eval,
    
    Estadistica_descriptiva =
      tabla_descriptiva,
    
    Metricas_concordancia =
      tabla_metricas,
    
    Modelos_calibracion =
      tabla_modelos,
    
    Bland_Altman =
      tabla_bland_altman,
    
    Validacion_LOOCV =
      tabla_loocv,
    
    Bootstrap_OLS =
      tabla_bootstrap,
    
    Diagnostico_OLS =
      diagnostico_ols
  ),
  
  path = file.path(
    ruta_salida,
    "Resultados_calibracion_CE.xlsx"
  )
)


## 19. Exportación de resúmenes de modelos ---------------------

capture.output(
  summary(
    modelo_ols
  ),
  
  file = file.path(
    ruta_salida,
    "Resumen_modelo_OLS.txt"
  )
)

capture.output(
  printSummary(
    modelo_deming
  ),
  
  file = file.path(
    ruta_salida,
    "Resumen_modelo_Deming.txt"
  )
)

capture.output(
  printSummary(
    modelo_pb
  ),
  
  file = file.path(
    ruta_salida,
    "Resumen_modelo_Passing_Bablok.txt"
  )
)


## 20. Ecuaciones candidatas ----------------------------------

cat(
  "\n========================================\n"
)

cat(
  "Ecuación OLS candidata:\n"
)

cat(
  "CE_lab = ",
  round(
    intercepto_ols,
    4
  ),
  " + ",
  round(
    pendiente_ols,
    4
  ),
  " × CE_sensor\n",
  sep = ""
)

cat(
  "\nEcuación Deming candidata:\n"
)

cat(
  "CE_lab = ",
  round(
    intercepto_deming,
    4
  ),
  " + ",
  round(
    pendiente_deming,
    4
  ),
  " × CE_sensor\n",
  sep = ""
)

cat(
  "\nEcuación Passing–Bablok candidata:\n"
)

cat(
  "CE_lab = ",
  round(
    intercepto_pb,
    4
  ),
  " + ",
  round(
    pendiente_pb,
    4
  ),
  " × CE_sensor\n",
  sep = ""
)

cat(
  "\n========================================\n"
)

cat(
  "\nResultados guardados en:\n",
  ruta_salida,
  "\n"
)

cat(
  "\nNo aplicar todavía una ecuación al conjunto espacial ",
  "hasta revisar RMSE, LOOCV, Bland–Altman, residuos ",
  "e intervalos de confianza.\n"
)