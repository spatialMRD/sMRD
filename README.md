# sMRD Code

This code is used for analysis and visualization of data from the sMRD manuscript.

# Usage

Synthetic example data is provided under:
```
data/sample_info/ex_sample_info.csv -- Sample metadata
data/features/ex_cnv.csv -- CNV calls
data/features/ex_snv.csv -- SNV calls
data/features/ex_enrichment.csv -- Raw Variant Calls
```

To generate intermediate data and figures, create a conda environment, activate, and run the following script:
```
conda env create -f environment.yml
conda activate smrd
./run_analysis.sh
```

To perform variant enrichment analysis, run the R Markdown notebook under `analysis/enrichment_analysis.Rmd`.
