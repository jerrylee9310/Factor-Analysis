# Project Analysis

This repository contains the analysis pipeline and replication materials for the paper: `On the Fisher's Additive Model: Revisiting Heritability through Genotype Encoding`. The `notebooks/` directory includes Jupyter Notebooks that guide users through the process of replicating the GRM creation and analysis described in the paper.

---

## Key Directories

### `notebooks/`
The `notebooks/` directory is organized into two subdirectories: `real/` and `simulation/`. Each contains Jupyter Notebooks for specific aspects of the analysis.

#### `simulation/`
This folder contains notebooks for simulation-based analyses:

1. **`H2Comparison.ipynb`**  
This notebook explores how genotype encoding schemes (AD vs. Factor models) influence heritability estimates. It uses simulations under purely additive genetic architectures to compare twin-based and SNP-based heritability estimates. Results highlight discrepancies caused by encoding differences, replicating the twin-SNP heritability gap observed in empirical studies.

2. **`Twin.r_DZ.ipynb`**  
This notebook simulates genetic relatedness in dizygotic (DZ) twins under two encoding models: Additive-Dominance (AD) and Factor. It explores how minor allele frequency (MAF) influences genetic covariance components (e.g., additive, dominance, heterozygous, homozygous). Results highlight distinct patterns of genetic relatedness captured by each model, with the Factor model providing a more detailed breakdown of heterozygous and homozygous effects.

#### `real/`
This folder contains notebooks for analyzing real-world datasets:

1. **`REML.ipynb`**  
This notebook demonstrates the use of Restricted Maximum Likelihood (REML) analysis to estimate variance components for heritability. It generates multi-GRM (genetic relatedness matrix) files for both Additive-Dominance (AD) and Factor models, then runs REML using GCTA software on phenotypic data. The results provide insights into how different genotype encoding schemes influence variance decomposition.

2. **`snph2.plot.ipynb`**  
This notebook provides a comprehensive guide to analyzing and visualizing SNP-based heritability results. It includes parsing REML output, performing meta-analysis, and calculating heritability components for both Additive-Dominance (AD) and Factor models. The notebook also generates detailed plots, such as heritability estimates, significance patterns, and covariance analyses, to interpret the genetic architecture of phenotypes.

3. **`snpToH2.ipynb`**  
This notebook maps SNP-based variance components to twin-based heritability estimates under the Factor model. It simulates MAF-dependent genetic relatedness for dizygotic twins and combines these results with UK Biobank MAF distributions to compute expected twin correlations. Using Falconer’s formula, it bridges SNP-based and twin-based heritability estimates for direct comparison with empirical studies.

4. **`H2.snp-twin.ipynb`**  
This notebook estimates SNP-based heritability using twin data and compares it with twin-based heritability estimates. It preprocesses data, computes heritability components under Additive-Dominance (AD) and Factor models, and applies Falconer’s formula to bridge the two methods. Statistical tests and visualizations, including heritability ratios and confidence intervals, provide insights into the alignment between SNP-based and twin-based heritability estimates.

---

## License

This project is licensed under the [MIT License](LICENSE)