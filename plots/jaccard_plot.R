library(dplyr)
library(ggplot2)

DATA_DIR <- "data/output"
PLOT_DIR <- "data/plots"
# get jaccard stats from DATA_DIR
df <- read.csv(file.path(DATA_DIR, "jaccard_stats.csv"))

create_boxplot <- function(data, metric, y_label) {
  # Perform a t-test between progression groups
  t_test <- t.test(data[[metric]] ~ data$progression, data = data, var.equal = TRUE)
  p_value <- t_test$p.value

  # Print the p-value to the CLI
  #cat(metric, "p-value:", p_value, "\n")
  # print the metric, t test stat and p value
  cat(metric, "t-test statistic:", t_test$statistic, "\n")
  cat(metric, "p-value:", p_value, "\n")

  # Calculate positions for significance bar
  max_y <- max(data[[metric]], na.rm = TRUE)
  y_bar <- max_y * 1.05

  # Create the boxplot
  plot <- ggplot(data, aes(x = factor(progression, labels = c("Non-Progressors", "Progressors")), y = .data[[metric]])) +
    geom_boxplot(outlier.shape = NA, width = 0.6) +
    geom_jitter(width = 0.2, alpha = 0.5) +
    labs(
      title = NULL,
      x = NULL,
      y = y_label
    ) +
    theme_classic(base_size = 14) +
    theme(
      panel.grid = element_blank(),
      axis.line.x.bottom = element_line(color = "black"),
      axis.line.y.left = element_line(color = "black"),
      axis.title.x = element_blank()
    ) +
    annotate("segment", x = 1, xend = 2, y = y_bar, yend = y_bar, color = "black") +
    annotate("text", x = 1.5, y = y_bar * 1.02, label = "*", size = 5)

  return(plot)
}

metrics <- list(
    mean = "Mean Jaccard Index",
    std = "Standard Deviation Jaccard Index"
)

plots <- lapply(names(metrics), function(metric) create_boxplot(df, metric, metrics[[metric]]))

# Arrange the plots side by side
library(gridExtra)
output_plot <- grid.arrange(grobs = plots, ncol = 2, bottom = "Progression Status")

# Save the plot as a smaller figure
ggsave(file.path(PLOT_DIR, "jaccard_boxplots.png"), plot = output_plot, width = 8, height = 4)
# also save as a vector graphic
ggsave(file.path(PLOT_DIR, "jaccard_boxplots.svg"), plot = output_plot, width = 8, height = 4)
