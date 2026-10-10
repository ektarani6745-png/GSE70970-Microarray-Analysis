# Microarray Differential Expression Analysis

An R-based pipeline utilizing `limma` to analyze GEO dataset GSE70970.

## Workflow
1. Extracts and formats the raw series matrix.
2. Parses clinical metadata to classify patients into `Metastasis` vs `No_Metastasis` groups.
3. Performs array normalization and fits a linear model (`lmFit`).
4. Generates data tables for all, upregulated, and downregulated genes.
5. Exports high-resolution Volcano Plots and Heatmaps for the top 50 differentially expressed genes.
