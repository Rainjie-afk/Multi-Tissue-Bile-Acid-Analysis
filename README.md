# Multi-Tissue Bile Acid Analysis

This repository provides a workflow for **multi-compartment bile acid analysis in NAFLD**, integrating bile acid profiles from the **liver, plasma, ileum, cecum, and feces**.

The workflow incorporates multivariate analyses implemented using the **mixOmics** framework (https://mixomics.org/) together with complementary statistical analyses. It includes:

* **DIABLO** for multiblock integration of bile acid profiles across tissues and biological compartments
* **Circos plots** for visualization of cross-compartment associations
* **Correlation analysis** for identifying relationships among bile acids across compartments
* **Piecewise structural equation modeling (piecewise SEM)** for evaluating potential directional relationships among liver, plasma, ileum, cecum, and fecal bile acid profiles
* **Optional sPLS-DA analysis** for supervised discrimination and feature selection

The workflow is designed to characterize coordinated alterations in bile acid metabolism across multiple anatomical compartments and to investigate potential inter-compartment relationships in NAFLD.
