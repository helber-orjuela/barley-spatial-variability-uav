[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21798804.svg)](https://doi.org/10.5281/zenodo.21798804)
# Barley spatial variability and UAV analysis

Reproducible workflow for spatial variability analysis integrating soil properties, crop traits, geostatistics, treatment effects, point-pattern analysis, and UAV-derived vegetation indices in barley experimental plots.

---

## Overview

This repository contains the complete analytical workflow developed to evaluate spatial variability in a barley experimental field.

The workflow integrates:

- soil data auditing;
- electrical conductivity calibration;
- exploratory and directional semivariogram analysis;
- anisotropy assessment;
- ordinary kriging;
- spatial cross-validation;
- treatment-effect models;
- residual spatial diagnostics;
- point-pattern analysis;
- soil–crop integration;
- UAV orthomosaic processing;
- GLI and GRVI calculation;
- temporal comparison between UAV acquisition dates.

The repository was prepared to support the reproducibility of the revised manuscript and to provide reviewers and future users with access to the scripts, processed outputs, and methodological documentation.

---

## Objectives

The analytical workflow was designed to:

1. audit and standardize soil, crop, and spatial datasets;
2. evaluate the concordance between sensor-based and laboratory electrical conductivity measurements;
3. characterize the spatial dependence of soil variables;
4. evaluate anisotropy and fit geostatistical models;
5. generate kriging prediction surfaces and assess their predictive performance;
6. evaluate treatment and fertilizer-dose effects;
7. diagnose residual spatial autocorrelation;
8. characterize the spatial configuration of experimental sampling nodes;
9. integrate predicted soil properties with crop-response variables;
10. derive UAV vegetation indices from RGB orthomosaics;
11. compare GLI and GRVI between 9 and 14 November;
12. provide a reproducible workflow for the manuscript and supplementary material.

---

## Repository structure

```text
barley-spatial-variability-uav/
│
├── README.md
├── .gitignore
│
└── code/
    ├── 01_audit_i2/
    ├── 02_sensor_ec/
    ├── 03_geostatistics/
    ├── 04_treatments/
    ├── 05_point_pattern/
    ├── 06_soil_crop/
    └── 07_uav/
