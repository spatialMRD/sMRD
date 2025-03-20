#!/usr/bin/env Rscript

# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyverse)
library(ape)  # For neighbor joining and tree manipulation
# library(ggtree)  # For tree visualization - removed due to installation issues
library(wesanderson)  # For color palettes

DATA_DIR <- "data/"
base_dir <- paste0(DATA_DIR, "output/phylo/")

# Get all patient directories
patient_dirs <- list.dirs(base_dir, full.names = TRUE, recursive = FALSE)

# Initialize a list to store all phylo distance files
all_phylo_files <- c()

# Loop through each patient directory to find phylo distance files
for (dir in patient_dirs) {
  # Extract patient ID from directory path
  patient_id <- basename(dir)
  
  # Look for the phylo distances file
  phylo_file <- file.path(dir, paste0(patient_id, "_phylo_distances.csv"))
  
  # Check if file exists and add to list
  if (file.exists(phylo_file)) {
    all_phylo_files <- c(all_phylo_files, phylo_file)
  } else {
    cat("File not found:", phylo_file, "\n")
  }
}

# Print the number of files found
cat("Found", length(all_phylo_files), "phylo distance files\n")

# Function to read and process a phylo distance file
read_phylo_file <- function(file_path) {
  # load the file
  df <- read.csv(file_path)
  # convert to a flattened dataframe with cols: sample_id1, sample_id2, distance
  df <- as.data.frame(df)
  # pivot longer setting first col as sample_id1 and other cols as sample_id2
  df <- pivot_longer(df, cols = -1, names_to = "sample_id2", values_to = "distance")
  # reset colname X to sample_id1
  colnames(df)[1] <- "sample_id1"
  # replace "." in sample_id2 with "-"
  df$sample_id2 <- gsub("\\.", "-", df$sample_id2)
  # add a patient_id column
  df$patient_id <- basename(dirname(file_path))
  return(df)
}

# Read all files into a list of dataframes
phylo_dfs <- lapply(all_phylo_files, read_phylo_file)

# Combine all dataframes into one
if (length(phylo_dfs) > 0) {
  all_phylo_data <- bind_rows(phylo_dfs)
  cat("Combined data has", nrow(all_phylo_data), "rows\n")
} else {
  all_phylo_data <- data.frame()
  cat("No phylo distance files were found\n")
}

# Load the spatial samples data
spatial_samples <- read.csv(paste0(DATA_DIR, "output/all_spatial_dfs.csv"))

# Load the sample info df for special case processing
sample_info <- read.csv(paste0(DATA_DIR, "sample_info/ex_sample_info.csv"))

# Get a list of all sample_id in spatial_samples
spatial_sample_ids <- unique(spatial_samples$sample_id)

# Filter phylo data to only include rows where sample_id1 and sample_id2 are in the spatial_sample_ids list
all_phylo_data <- all_phylo_data %>%
  filter(sample_id1 %in% spatial_sample_ids & sample_id2 %in% spatial_sample_ids)

# Add patient_id column to spatial_samples, split sample_id by "-" and take the first element
spatial_samples$patient_id <- spatial_samples$sample_id %>% 
  str_split("-") %>% 
  map_chr(1)

# Get the root sample for each patient_id
root_sample <- spatial_samples %>%
  group_by(patient_id) %>%
  arrange(genomic_dist) %>%
  slice(1) %>%
  ungroup()

# Get this as a map, key = patient_id, value = sample_id
root_sample_map <- setNames(root_sample$sample_id, root_sample$patient_id)

# print root_sample_map to terminal
# Print root sample map to terminal
cat("Root sample map:\n")
for (patient in names(root_sample_map)) {
  cat(patient, "->", root_sample_map[patient], "\n")
}

### CREATE PHYLO TREES IN 2D ###

# Function to create a phylogenetic tree for a patient using base R plotting
create_phylo_tree <- function(patient_id, phylo_data, root_sample_id) {
  # Create output directory for this patient
  patient_dir <- file.path(paste0(DATA_DIR, "output/phylo"), patient_id)
  dir.create(patient_dir, showWarnings = FALSE)
  
  # Filter data for this patient if patient_id doesn't contain "_left" or "_right"
  if (!grepl("_left$|_right$", patient_id)) {
    patient_data <- phylo_data %>%
      filter(patient_id == !!patient_id)
  } else {
    # For special cases (left/right), use the provided data directly
    patient_data <- phylo_data
    # Extract the base patient_id (without _left or _right)
    base_patient_id <- sub("_left$|_right$", "", patient_id)
  }
  
  # Get unique sample IDs for this patient
  sample_ids <- unique(c(patient_data$sample_id1, patient_data$sample_id2))
  
  # Check if we have at least 3 samples (required for NJ tree)
  if (length(sample_ids) < 3) {
    cat("Patient", patient_id, "has less than 3 samples. Skipping.\n")
    return(NULL)
  }
  
  # Create a distance matrix
  n <- length(sample_ids)
  dist_matrix <- matrix(0, nrow = n, ncol = n)
  rownames(dist_matrix) <- sample_ids
  colnames(dist_matrix) <- sample_ids
  
  # Fill the distance matrix
  for (i in 1:nrow(patient_data)) {
    row <- patient_data[i, ]
    dist_matrix[row$sample_id1, row$sample_id2] <- row$distance
    dist_matrix[row$sample_id2, row$sample_id1] <- row$distance  # Make it symmetric
  }
  
  # Convert to a dist object
  dist_obj <- as.dist(dist_matrix)
  
  # Apply neighbor joining algorithm
  nj_tree <- nj(dist_obj)
  
  # Root the tree at the specified root sample
  if (root_sample_id %in% sample_ids) {
    # First, make sure the root sample is at node 1 (leftmost position)
    rooted_tree <- root(nj_tree, which(nj_tree$tip.label == root_sample_id), resolve.root = TRUE)
    
    # Ensure the tree is properly ladderized to have a consistent appearance
    rooted_tree <- ladderize(rooted_tree, right = TRUE)
  } else {
    # If root sample not found, use the first sample as root
    cat("Root sample", root_sample_id, "not found for patient", patient_id, ". Using first sample as root.\n")
    rooted_tree <- ladderize(nj_tree, right = TRUE)
  }
  
  # Get genomic distances for coloring
  # For special cases, use the base patient_id to filter spatial_samples
  if (exists("base_patient_id")) {
    filter_patient_id <- base_patient_id
  } else {
    filter_patient_id <- patient_id
  }
  
  # Get the genomic distances from spatial_samples instead of calculating from the tree
  sample_genomic_dist <- spatial_samples %>%
    filter(patient_id == !!filter_patient_id) %>%
    select(sample_id, genomic_dist_scaled)
  
  # Create a node data frame with genomic distances from spatial_samples
  node_data <- data.frame(
    label = rooted_tree$tip.label,
    genomic_dist_scaled = sample_genomic_dist$genomic_dist_scaled[match(rooted_tree$tip.label, sample_genomic_dist$sample_id)]
  )
  
  # Save the tree data
  write.csv(node_data, file.path(patient_dir, paste0(patient_id, "_tree_data.csv")), row.names = FALSE)
  write.tree(rooted_tree, file.path(patient_dir, paste0(patient_id, "_tree.nwk")))
  
  # Define the color palette from wesanderson
  pal <- wes_palette("Zissou1", 100, type = "continuous")
  
  # Create a color mapping function
  color_map <- colorRampPalette(pal)
  
  # Map genomic distances to colors
  tip_colors <- color_map(100)[cut(node_data$genomic_dist_scaled, breaks = 100)]
  
  # Create a small lookup table for the legend
  legend_breaks <- seq(0, 1, length.out = 5)
  legend_colors <- color_map(5)
  
  # Save the plot as PNG
  png(file.path(patient_dir, paste0(patient_id, "_phylo_tree.png")), width = 1200, height = 800)
  
  # Set up the plotting area with space for a legend at the bottom
  layout(matrix(c(1, 2), 2, 1), heights = c(4, 1))
  
  # Plot the tree with the root on the left
  plot(rooted_tree, 
       direction = "rightwards",  # Make the tree grow from left to right
       show.tip.label = FALSE,    # Don't show the sample IDs
       main = paste("Phylogenetic Tree for Patient", patient_id),
       cex.main = 1.5,
       root.edge = TRUE)          # Show the root edge
  
  # Add colored dots at the tips (larger size)
  tiplabels(pch = 19, col = tip_colors, cex = 3.5)
  
  # Find the root sample index
  root_index <- which(rooted_tree$tip.label == root_sample_id)
  if (length(root_index) > 0) {
    # Color the root sample dot black
    tiplabels(pch = 19, col = "black", cex = 3.5, tip = root_index)
    # Add a "Root" label
    #tiplabels(text = "Root", tip = root_index, adj = c(-0.5, 0.5), frame = "none", cex = 0.8)
  }
  
  # Add a horizontal color bar at the bottom
  par(mar = c(4, 4, 2, 4))
  plot(0, 0, type = "n", axes = FALSE, xlab = "", ylab = "", xlim = c(0, 1), ylim = c(0, 1))
  title("Scaled Genomic Distance", cex.main = 1.2)
  
  # Draw the horizontal color scale
  x_positions <- seq(0.1, 0.9, length.out = 100)
  rect_width <- x_positions[2] - x_positions[1]
  
  for (i in 1:100) {
    rect(x_positions[i], 0.4, x_positions[i] + rect_width, 0.6, 
         col = color_map(100)[i], border = NA)
  }
  
  # Add labels to the color scale
  text_positions <- seq(0.1, 0.9, length.out = 5)
  text_labels <- seq(0, 1, length.out = 5)
  
  for (i in 1:5) {
    text(text_positions[i], 0.25, labels = text_labels[i], cex = 0.8)
  }
  
  dev.off()
  
  # Save the plot as PDF
  pdf(file.path(patient_dir, paste0(patient_id, "_phylo_tree.pdf")), width = 12, height = 8)
  
  # Set up the plotting area with space for a legend at the bottom
  layout(matrix(c(1, 2), 2, 1), heights = c(4, 1))
  
  # Plot the tree with the root on the left
  plot(rooted_tree, 
       direction = "rightwards",  # Make the tree grow from left to right
       show.tip.label = FALSE,    # Don't show the sample IDs
       main = paste("Phylogenetic Tree for Patient", patient_id),
       cex.main = 1.5,
       root.edge = TRUE)          # Show the root edge
  
  # Add colored dots at the tips (larger size)
  tiplabels(pch = 19, col = tip_colors, cex = 3)
  
  # Find the root sample index
  root_index <- which(rooted_tree$tip.label == root_sample_id)
  if (length(root_index) > 0) {
    # Color the root sample dot black
    tiplabels(pch = 19, col = "black", cex = 3, tip = root_index)
    # Add a "Root" label
    tiplabels(text = "Root", tip = root_index, adj = c(-0.5, 0.5), frame = "none", cex = 0.8)
  }
  
  # Add a horizontal color bar at the bottom
  par(mar = c(4, 4, 2, 4))
  plot(0, 0, type = "n", axes = FALSE, xlab = "", ylab = "", xlim = c(0, 1), ylim = c(0, 1))
  title("Scaled Genomic Distance", cex.main = 1.2)
  
  # Draw the horizontal color scale
  x_positions <- seq(0.1, 0.9, length.out = 100)
  rect_width <- x_positions[2] - x_positions[1]
  
  for (i in 1:100) {
    rect(x_positions[i], 0.4, x_positions[i] + rect_width, 0.6, 
         col = color_map(100)[i], border = NA)
  }
  
  # Add labels to the color scale
  text_positions <- seq(0.1, 0.9, length.out = 5)
  text_labels <- seq(0, 1, length.out = 5)
  
  for (i in 1:5) {
    text(text_positions[i], 0.25, labels = text_labels[i], cex = 0.8)
  }
  
  dev.off()
  
  # Return the tree object
  return(rooted_tree)
}

# Process each patient
patient_ids <- unique(all_phylo_data$patient_id)
for (pid in patient_ids) {
  cat("Processing patient:", pid, "\n")
  root_id <- root_sample_map[pid]
  if (is.na(root_id)) {
    cat("No root sample found for patient", pid, ". Skipping.\n")
    next
  }
  
  tryCatch({
    create_phylo_tree(pid, all_phylo_data, root_id)
  }, error = function(e) {
    cat("Error processing patient", pid, ":", e$message, "\n")
  })
}

# Debug: print structure of all_phylo_data and sample_info
cat("\nStructure of all_phylo_data:\n")
str(all_phylo_data)

cat("\nSample of all_phylo_data:\n")
print(head(all_phylo_data, 3))

cat("\nStructure of sample_info:\n")
str(sample_info)

cat("\nSample of sample_info:\n")
print(head(sample_info[, c("sample_id", "tumor_vs_lymphnode")], 3))

# print columns of sample_info to terminal
cat("Columns of sample_info:\n")
print(colnames(sample_info))