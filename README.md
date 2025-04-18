# sMRD Code

This code is used for analysis and visualization of data from the sMRD manuscript.

# Installation

Install [anaconda or miniconda.](https://www.anaconda.com/docs/getting-started/miniconda/install#quickstart-install-instructions)

Then, create and activate a conda environment:
```
conda env create -f environment.yml
conda activate smrd
```

Installation should complete within a couple of minutes on a standard desktop computer.

# Usage

Synthetic example data is provided under:
```
data/sample_info/ex_sample_info.csv -- Sample metadata
data/features/ex_cnv.csv -- CNV calls
data/features/ex_snv.csv -- SNV calls
data/features/ex_enrichment.csv -- Raw Variant Calls
```

To generate intermediate data and figures, run the following script:
```
./run_analysis.sh
```

To perform variant enrichment analysis, run the R Markdown notebook under `analysis/enrichment_analysis.Rmd`.

To reproduce the manuscript figures or to run a custom analysis, replace the provided synthetic data with a new dataset.

# Output

Tabular output, phylogenetic trees, and statistics are written to the `data/output` directory. Plots are written to the `data/plots` directory. Runtime for example data should not exceed a couple of minutes on a standard desktop computer.
