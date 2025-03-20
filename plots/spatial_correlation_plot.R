library(ggplot2)
library(dplyr)

DATA_DIR <- "data/"
df <- read.csv(paste0(DATA_DIR, "output/spatial_correlation.csv"))

create_boxplot <- function(data, metric, y_label) {
  # Perform a t-test between progression groups
  t_test <- t.test(data[[metric]] ~ data$progression, data = data, var.equal = TRUE)
  p_value <- t_test$p.value

  # Print the p-value to the CLI
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

# generate the plot
p <- create_boxplot(df, "corr", "Pearson Correlation")

ggsave(paste0(DATA_DIR, "plots/spatial_correlation.png"), plot = p, width = 4, height = 4)
ggsave(paste0(DATA_DIR, "plots/spatial_correlation.svg"), plot = p, width = 4, height = 4)
