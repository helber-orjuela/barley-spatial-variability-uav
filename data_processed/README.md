# Data directory

## Overview

This directory contains the datasets used throughout the reproducible workflow. The information is organized according to the analytical stages described in the manuscript rather than by file type, allowing each methodological module to be reproduced independently.

```
data/
│
├── I_FASE_I2
├── IA_I1_COMPACTO
├── II_SENSOR_CE
├── III_TRATAMIENTOS
├── IV_SUELO_CULTIVO
├── V_UAV
└── README.md
```

---

# Folder description

## I_FASE_I2

Contains the original field datasets used during the initial data auditing stage.

Typical contents include:

- Sampling nodes
- Crop measurements
- Original soil observations
- Experimental design information

These files constitute the primary inputs for the analytical workflow.

---

## IA_I1_COMPACTO

Contains the compact spatial dataset generated from Phase I.

This dataset is used as the reference for exploratory geostatistical analyses and spatial validation.

---

## II_SENSOR_CE

Contains the datasets used during the electrical conductivity calibration workflow.

Typical files include:

- laboratory EC measurements;
- sensor observations;
- exploratory calibration datasets;
- OLS-corrected conductivity tables.

---

## III_TRATAMIENTOS

Contains the datasets required for treatment analyses.

Typical information includes:

- fertilizer doses;
- omission treatments;
- replication structure;
- response variables used in statistical modelling.

---

## IV_SUELO_CULTIVO

Contains the integrated soil–crop datasets generated after geostatistical processing.

These datasets combine:

- kriging predictions;
- soil properties;
- crop measurements;
- integrated modelling tables.

These files constitute the main analytical datasets used in the manuscript.

---

## V_UAV

Contains the UAV-derived datasets.

Typical contents include:

- extracted GLI values;
- extracted GRVI values;
- temporal comparison tables;
- UAV indices extracted at crop sampling nodes.

Original orthomosaics are not stored directly in this repository because of their size. They are archived separately through Zenodo.

---

# Data policy

This repository contains only the datasets required to reproduce the analyses presented in the manuscript.

Intermediate files generated automatically during script execution are intentionally omitted to reduce repository size and improve navigation.

Large raster products, temporary objects and auxiliary processing files are distributed separately when required.

---

# Coordinate reference system

Spatial datasets use the same projected coordinate reference system adopted throughout the manuscript.

---

# Relationship with the repository

Each data folder corresponds directly to one or more analytical modules contained in the `code/` directory.

This parallel organization facilitates workflow reproducibility and allows users to execute each stage independently.