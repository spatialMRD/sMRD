#!/usr/bin/env Rscript

# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyverse)
library(ape) # For neighbor joining and tree manipulation
library(plotly)
library(wesanderson)
library(reticulate)
# Get the current Python executable path from the system
python_path <- Sys.which("python")
use_python(python_path)

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

# Print root sample map to terminal
cat("Root sample map:\n")
for (patient in names(root_sample_map)) {
  cat(patient, "->", root_sample_map[patient], "\n")
}

# Function to create a phylogenetic tree and extract connections
create_phylo_tree_and_connections <- function(patient_id, phylo_data, root_sample_id) {
  # Filter data for this patient
  patient_data <- phylo_data %>%
    filter(patient_id == !!patient_id)

  # Get unique sample IDs for this patient
  sample_ids <- unique(c(patient_data$sample_id1, patient_data$sample_id2))

  # Check if we have at least 3 samples (required for NJ tree)
  if (length(sample_ids) < 3) {
    cat("Patient", patient_id, "has less than 3 samples. Using direct MST approach.\n")
    # Fall back to the original approach for 2 samples
    if (length(sample_ids) == 2) {
      connections <- data.frame(
        sample1 = sample_ids[1],
        sample2 = sample_ids[2]
      )
      return(list(connections = connections))
    } else {
      cat("Patient", patient_id, "has less than 2 samples. Skipping.\n")
      return(NULL)
    }
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
    dist_matrix[row$sample_id2, row$sample_id1] <- row$distance # Make it symmetric
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
    root_sample_id <- rooted_tree$tip.label[1]
  }

  # Extract distances from the NJ tree
  # Get the cophenetic distances from the tree
  tree_dist_matrix <- cophenetic.phylo(rooted_tree)

  # Initialize the list of added samples with the root sample
  added_samples <- c(root_sample_id)

  # Initialize the list of connections
  connections <- list()

  # Iteratively add samples until all are added (MST algorithm)
  while (length(added_samples) < length(sample_ids)) {
    min_dist <- Inf
    min_sample <- NULL
    min_connection <- NULL

    # For each sample not yet added
    for (sample in sample_ids[!(sample_ids %in% added_samples)]) {
      # Find the closest already-added sample
      for (added_sample in added_samples) {
        dist <- tree_dist_matrix[sample, added_sample]
        if (dist < min_dist) {
          min_dist <- dist
          min_sample <- sample
          min_connection <- added_sample
        }
      }
    }

    # Add the sample with minimum distance
    if (!is.null(min_sample)) {
      added_samples <- c(added_samples, min_sample)
      connections[[length(connections) + 1]] <- list(sample1 = min_connection, sample2 = min_sample)
    } else {
      # This should not happen, but just in case
      break
    }
  }

  # Convert the list to a data frame
  connections_df <- do.call(rbind, lapply(connections, as.data.frame))

  # Return the connections
  return(list(connections = connections_df))
}

# Function to create a 3D network plot
create_3d_network_plot <- function(points_df, connections_df, plot_name) {
  # Define the custom color palette
  # Convert wes_palette colors to a Plotly-compatible colorscale
  pal <- wes_palette("Zissou1", 100, type = "continuous")
  colorscale <- lapply(seq_along(pal), function(i) {
    list((i - 1) / (length(pal) - 1), pal[i])
  })

  # Identify the point with the lowest genomic_dist_scaled
  min_point <- points_df[which.min(points_df$genomic_dist_scaled), ]
  # update all points so min_point is 0,0,0
  points_df$x <- points_df$x - min_point$x
  points_df$y <- points_df$y - min_point$y
  points_df$z <- points_df$z - min_point$z
  min_point <- points_df[which.min(points_df$genomic_dist_scaled), ]

  # Create the basic scatter plot
  p <- plot_ly() %>%
    # Add points colored by genomic_dist_scaled with size scaled by genomic_dist_scaled
    add_trace(
      data = points_df,
      x = ~x, y = ~y, z = ~z,
      type = "scatter3d",
      mode = "markers",
      marker = list(
        size = ~ 6 + 14 * genomic_dist_scaled, # Base size of 6, scaling up to 20 based on genomic distance
        color = ~genomic_dist_scaled,
        colorscale = colorscale,
        colorbar = list(
          title = "Scaled Genomic Distance",
          titleside = "top",
          orientation = "h",
          xanchor = "center",
          x = 0.5,
          y = -.2, # Moved up from 0 to prevent overlap
          len = 0.6,
          thickness = 15,
          yanchor = "bottom"
        )
      ),
      showlegend = FALSE,
      hoverinfo = "text",
      text = ~ paste("Sample ID:", sample_id, "<br>Genomic Distance:", genomic_dist_scaled)
    ) %>%
    # Add the special point with lowest genomic_dist_scaled in black
    add_trace(
      data = min_point,
      x = ~x, y = ~y, z = ~z,
      type = "scatter3d",
      mode = "markers",
      marker = list(
        size = 6, # Smaller size for the root node
        color = "black"
      ),
      showlegend = FALSE,
      text = ~ paste("Root Sample ID:", sample_id, "<br>Genomic Distance:", genomic_dist_scaled),
      hoverinfo = "text"
    )

  # Add connections as lines
  if (!is.null(connections_df) && nrow(connections_df) > 0) {
    for (i in 1:nrow(connections_df)) {
      # Get coordinates for sample1
      point1 <- points_df[points_df$sample_id == connections_df$sample1[i], ]
      # Get coordinates for sample2
      point2 <- points_df[points_df$sample_id == connections_df$sample2[i], ]

      # Add line if both points exist in the points_df
      if (nrow(point1) > 0 && nrow(point2) > 0) {
        p <- p %>% add_trace(
          x = c(point1$x, point2$x),
          y = c(point1$y, point2$y),
          z = c(point1$z, point2$z),
          type = "scatter3d",
          mode = "lines",
          line = list(
            color = "black",
            width = 1.5,
            alpha = 0.7
          ),
          showlegend = FALSE,
          hoverinfo = "none"
        )
      }
    }
  }

  # Update layout with increased margins to prevent axis label truncation
  p <- p %>% layout(
    scene = list(
      xaxis = list(title = "X"),
      yaxis = list(title = "Y"),
      zaxis = list(title = "Z"),
      camera = list(
        eye = list(x = 1.6, y = 1.6, z = 1.6) # Zoomed out camera position
      )
    ),
    title = plot_name,
    margin = list(l = 50, r = 50, b = 50, t = 50) # Increased margins to prevent axis label truncation
  )

  # Save a png of the plot as a static image to the new directory
  save_image(p, paste0(DATA_DIR, "output/phylo/", plot_name, ".png"), scale = 2)
  # Save as a pdf to the new directory
  save_image(p, paste0(DATA_DIR, "output/phylo/", plot_name, ".pdf"), scale = 2)

  return(p)
}

# Process each patient
df <- spatial_samples
df$patient_id <- sapply(strsplit(df$sample_id, "-"), function(x) x[1])
all_pats <- unique(df$patient_id)

# Process regular patients (excluding special cases)
special_patients <- c("BC081", "BC089")
regular_patients <- all_pats

for (pat in regular_patients) {
  print(paste("Processing regular patient:", pat))

  # Get the root sample ID
  root_id <- root_sample_map[pat]
  if (is.na(root_id)) {
    cat("No root sample found for patient", pat, ". Skipping.\n")
    next
  }

  # Create phylogenetic tree and extract connections
  tryCatch(
    {
      result <- create_phylo_tree_and_connections(pat, all_phylo_data, root_id)
      if (!is.null(result)) {
        # Filter data for this patient
        df_patient <- df %>% filter(patient_id == pat)

        # Create 3D network plot
        create_3d_network_plot(df_patient, result$connections, pat)
      }
    },
    error = function(e) {
      cat("Error processing patient", pat, ":", e$message, "\n")
    }
  )
}
