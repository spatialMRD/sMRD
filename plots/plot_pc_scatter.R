library(ggplot2)
library(dplyr)
library(reticulate)
library("wesanderson")


DATA_DIR <- "data/"
df <- read.csv(paste0(DATA_DIR, "output/all_spatial_dfs.csv"))
conn_df <- read.csv(paste0(DATA_DIR, "output/all_spatial_connections.txt"), sep = "\t")



plot_pc_distances <- function(df) {
  # Load required library
  require(ggplot2)

  # Create a new dataframe with distance from origin for sorting
  df$dist_from_origin <- sqrt(df$PC1^2 + df$PC2^2)

  # Sort dataframe by distance from origin (descending)
  df_sorted <- df[order(-df$dist_from_origin), ]

  # Create the plot
  pal <- wes_palette("Zissou1", 100, type = "continuous")
  # Create the plot
  p <- ggplot() +
    # Add lines from origin to points
    geom_segment(
      data = df_sorted,
      aes(
        x = 0, y = 0,
        xend = PC1, yend = PC2,
        color = genomic_dist_scaled
      ),
      alpha = .6
    ) +
    # Add points colored by genomic_dist
    geom_point(
      data = df_sorted,
      aes(
        x = PC1, y = PC2,
        color = genomic_dist_scaled, size = 3 * genomic_dist_scaled
      ),
      alpha = .6
    ) +
    # Add origin point in blue
    geom_point(
      data = data.frame(x = 0, y = 0),
      aes(x = x, y = y),
      color = "black",
      size = 1, alpha = 1
    ) +
    scale_color_gradientn(colors = pal, name = "Genomic\nDistance") +
    scale_size_continuous(range = c(0.5, 3), guide = "none") + # Control size range and hide size legend
    labs(
      x = "Spatial PC1",
      y = "Spatial PC2"
    ) +
    theme_classic()
  # coord_cartesian(xlim = c(-4.5, 4.5), ylim = c(-4.5, 4.5))

  return(p)
}

df_progressors <- df %>% filter(progression == 1)
df_nonprogressors <- df %>% filter(progression == 0)

p <- plot_pc_distances(df_progressors)

# save this plot to file
ggsave(paste0(DATA_DIR, "plots/pc_distances_progressors.png"), p, width = 4, height = 4)
# save as a pdf
ggsave(paste0(DATA_DIR, "plots/pc_distances_progressors.pdf"), p, width = 4, height = 4)

p <- plot_pc_distances(df_nonprogressors)

ggsave(paste0(DATA_DIR, "plots/pc_distances_nonprogressors.png"), p, width = 4, height = 4)
ggsave(paste0(DATA_DIR, "plots/pc_distances_nonprogressors.pdf"), p, width = 4, height = 4)

library(dplyr)

compute_patient_stats <- function(data) {
  data %>%
    # Add a column for Euclidean distance from origin using PC1 and PC2
    mutate(pc_distance = sqrt(PC1^2 + PC2^2)) %>%
    group_by(patient_id) %>%
    summarise(
      median_pc_distance = mean(pc_distance), # Median PC distance for the patient
      mean_genomic_dist_below_median = mean(genomic_dist_scaled[pc_distance <= median_pc_distance]), # Mean genomic_dist_scaled below the median PC distance
      mean_genomic_dist_above_median = mean(genomic_dist_scaled[pc_distance > median_pc_distance]) # Mean genomic_dist_scaled above the median PC distance
    )
}
