#! /bin/bash
# Run the analysis
mkdir -p data/output
mkdir -p data/plots

# Run the preprocessing
echo "Running the preprocessing..."
python ./analysis/generate_phylo_analysis.py
python ./analysis/process_spatial_data.py
python ./analysis/jaccard_index.py
python ./analysis/spatial_correlation.py

# plots
echo "Generating plots..."
Rscript ./plots/jaccard_plot.R
Rscript ./plots/spatial_correlation_plot.R
Rscript ./plots/plot_2d_phylo_trees.R
Rscript ./plots/plot_2d_phylo_trees.R
Rscript ./plots/plot_3d_phylo_trees.R

# remove tmp files
rm Rplots.pdf

echo "Done"
