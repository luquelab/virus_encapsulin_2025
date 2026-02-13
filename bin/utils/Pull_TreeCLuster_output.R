# Load necessary libraries
library(tidyverse)
library(gridExtra)
library(ggplot2)

# Define thresholds for phage categorization
small_threshold <- 11600
mini_threshold <- 5000

# Define file paths
clusters_file <- "/Users/abelardoaguilar/projects/github_repos/mini-devel/mini-devel/bin/Modules_refactoring/Clustering_module/results/Proteins_identified_as_MCPs_VIBRANT_efamHMMSEARCH.faa_avg_clade_clusters.txt"
metadata_file <- "~/projects/github_repos/mini-devel/mini-devel/bin/Modules_refactoring/1_VIBRANT/data/DTRs_20kb.csv"
output_dir <- "~/projects/github_repos/mini-devel/mini-devel/bin/Modules_refactoring/Clustering_module/results/Absolute_clustering_refactoring_5000/"
colname_identifiers_metadata <- "genome_id"
colname_identifiers_clusters <- "SequenceName"

# Make directory if it does not exist
dir.create(output_dir, showWarnings = FALSE)

# Read data
clusters_df <- read_delim(clusters_file, delim = "\t", col_names = TRUE)
metadata_df <- read_csv(metadata_file)

# Find equivalent identifiers
cluster_ids <- clusters_df[[colname_identifiers_clusters]]
metadata_ids <- metadata_df[[colname_identifiers_metadata]]

find_metadata_for_a_given_cluster_id <- function(cluster_id, metadata_df){
  # Extract metadata identifiers
  metadata_ids <- metadata_df[[colname_identifiers_metadata]]
  # Find which metadata_id is a substring of the cluster_id
  for (metadata_id in metadata_ids) {
    if (grepl(metadata_id, cluster_id)) {
      return(metadata_id)
    }
  }
}


# Process cluster data
clusters_df <- clusters_df %>%
  mutate(
    Genome_ID = sapply(cluster_ids, find_metadata_for_a_given_cluster_id, metadata_df),
    ABSOLUTE_CLUSTER = ifelse(grepl("-1", ClusterNumber),"Unclustered",ClusterNumber),
  )
# Transform ABSOLUTE_CLUSTER to character
clusters_df$ABSOLUTE_CLUSTER <- as.character(clusters_df$ABSOLUTE_CLUSTER)

# Summarize clusters
per_cluster_df <- clusters_df %>%
  group_by(ABSOLUTE_CLUSTER) %>%
  summarise(
    n = n(),
    max_length = 0, min_length = 0, max_length_ID = "", min_length_ID = "",
    length_difference = 0, median_length = 0, mean_length = 0,
    biomes1_n = 0, biomes2_n = 0, biomes3_n = 0,
    biomes1 = "", biomes2 = "", biomes3 = "",
    genomes_ids = "", is_new_MCP = FALSE,
    mini_contig_a = FALSE, mini_contig_b = FALSE, mini_contig_c = FALSE
  )

# Function to calculate cluster attributes
calculate_attributes <- function(cluster_row, clusters_df, metadata_df) {
  cur_genomes <- clusters_df %>% filter(ABSOLUTE_CLUSTER == cluster_row$ABSOLUTE_CLUSTER) %>% pull(Genome_ID)
  cur_metadata <- metadata_df %>% filter(genome_id %in% cur_genomes)
  cur_lengths <- cur_metadata$contig_length
  
  # Handle empty cur_lengths by setting default values
  cluster_row$max_length <- ifelse(length(cur_lengths) > 0, max(cur_lengths, na.rm = TRUE), NA)
  cluster_row$min_length <- ifelse(length(cur_lengths) > 0, min(cur_lengths, na.rm = TRUE), NA)
  cluster_row$length_difference <- ifelse(length(cur_lengths) > 0, cluster_row$max_length - cluster_row$min_length, NA)
  cluster_row$median_length <- ifelse(length(cur_lengths) > 0, median(cur_lengths, na.rm = TRUE), NA)
  cluster_row$mean_length <- ifelse(length(cur_lengths) > 0, mean(cur_lengths, na.rm = TRUE), NA)
  cluster_row$max_length_ID <- ifelse(length(cur_lengths) > 0, cur_metadata %>% filter(contig_length == cluster_row$max_length) %>% pull(genome_id) %>% first(), NA)
  cluster_row$min_length_ID <- ifelse(length(cur_lengths) > 0, cur_metadata %>% filter(contig_length == cluster_row$min_length) %>% pull(genome_id) %>% first(), NA)
  
  # Calculate biome percentages
  for (i in 1:3) {
    biome_col <- paste0("biome", i)
    biome_table <- table(cur_metadata[[biome_col]])
    biome_percentages <- ifelse(length(biome_table) > 0, paste(names(biome_table), "_percentage_", round(100 * biome_table / sum(biome_table)), collapse = ","), "")
    cluster_row[[paste0("biomes", i)]] <- biome_percentages
    cluster_row[[paste0("biomes", i, "_n")]] <- length(biome_table)
  }
  
  cluster_row$genomes_ids <- paste(cur_genomes, collapse = ",")
  cluster_row$is_new_MCP <- grepl("efam", cluster_row$ABSOLUTE_CLUSTER)
  
  # Check for mini contigs
  check_mini_contig <- function(ids) any(grepl(ids, cluster_row$genomes_ids))
  cluster_row$mini_contig_a <- check_mini_contig("DTR_29916|DTR_758375")
  cluster_row$mini_contig_b <- check_mini_contig("DTR_605894|DTR_610893")
  cluster_row$mini_contig_c <- check_mini_contig("DTR_204949|DTR_248616")
  
  return(as.data.frame(cluster_row))
}

# Apply function to each row of per_cluster_df
per_cluster_df <- per_cluster_df %>%
  rowwise() %>%
  do(calculate_attributes(., clusters_df, metadata_df)) %>%
  ungroup()

# Print interesting clusters
interesting_clusters <- per_cluster_df %>%
  filter(mini_contig_a | mini_contig_b | mini_contig_c) %>%
  pull(ABSOLUTE_CLUSTER) %>%
  unique()
print(paste("Interesting clusters:", paste(interesting_clusters, collapse = ", ")))

# Merge clusters with metadata
clans_jgi_merge_df <- merge(clusters_df, metadata_df, by.x = "Genome_ID", by.y = colname_identifiers_metadata, all = TRUE)
clans_jgi_merge_df$ABSOLUTE_CLUSTER[is.na(clans_jgi_merge_df$ABSOLUTE_CLUSTER)] <- "NOT_INCLUDED"

# Calculate summary statistics
summary_stats <- clans_jgi_merge_df %>%
  group_by(ABSOLUTE_CLUSTER) %>%
  summarise(
    N = n(),
    min_length = min(contig_length, na.rm = TRUE),
    max_length = max(contig_length, na.rm = TRUE),
    median_length = median(contig_length, na.rm = TRUE),
    mean_length = mean(contig_length, na.rm = TRUE),
    std_deviation = sd(contig_length, na.rm = TRUE),
    cv = std_deviation / mean_length
  )

# Find the largest cluster and exclude it from color determination
largest_cluster <- summary_stats %>% filter(N == max(N)) %>% pull(ABSOLUTE_CLUSTER)
summary_stats_filtered <- summary_stats %>% filter(ABSOLUTE_CLUSTER != largest_cluster)

# Determine color categories for small to medium-sized clusters
n_ranges <- cut(
  summary_stats_filtered$N, 
  breaks = c(0, 10, 20, 50, 100, max(summary_stats_filtered$N) + 1), 
  labels = c("Very Small", "Small", "Medium", "Large", "Very Large"), 
  right = FALSE
)

# Plotting functions
plot_boxplot <- function(df, summary_df, summary_df_filtered, n_ranges) {
  median_order <- summary_df %>% arrange(median_length) %>% pull(ABSOLUTE_CLUSTER)
  median_order <- median_order[!(median_order %in% c("NO_CLUSTER", "VIBRANT_efam_0","NOT_INCLUDED"))]
  median_order <- unique(c(median_order, "Unclustered", "NOT_INCLUDED", "VIBRANT_efam_0"))
  
  # Define a color palette
  color_palette <- c("Very Small" = "#E5E5E5", "Small" = "#BEBEBE", "Medium" = "#7F7F7F", "Large" = "#4D4D4D", "Very Large" = "#1A1A1A")
  
  # Create a color mapping based on ranges
  cluster_colors <- setNames(as.character(color_palette[n_ranges]), summary_stats_filtered$ABSOLUTE_CLUSTER)
  cluster_colors[largest_cluster] <- "black"  # Special color for the largest cluster
  
  df %>%
    mutate(ABSOLUTE_CLUSTER = factor(ABSOLUTE_CLUSTER, levels = median_order)) %>%
    ggplot(aes(x = ABSOLUTE_CLUSTER, y = contig_length, fill = ABSOLUTE_CLUSTER)) +
    geom_boxplot() +
    geom_text(data = summary_df, aes(angle = 90, x = ABSOLUTE_CLUSTER, label = sprintf("N: %d Avg Len: %.1f CV: %.2f", N, mean_length, cv),
                                     y = 28000)) +
    geom_hline(yintercept = small_threshold, linetype = "dashed", color = "blue", size = 1) +
    geom_hline(yintercept = mini_threshold, linetype = "dashed", color = "red", size = 1) +
    scale_fill_manual(values = cluster_colors) +
    theme_minimal() +
    theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Boxplot of Contig Length by Cluster (Sorted by Median)",
         x = "Cluster",
         y = "Contig Length") +
    ylim(NA, 35000) +
    theme(plot.margin = margin(20, 20, 20, 20))  # Add margin around the plot
}

plot_histogram <- function(df, title) {
  bin_width <- 500
  breaks_hist_Toni <- seq.int(from = 0, to = 20000, by = bin_width)
  
  labels_hist_Toni <- c(
    rep("1-3000", 6),
    rep("3001-6000", 6),
    rep("6001-8500", 5),
    rep("8501-10500", 4),
    rep("10501-15000", 9),
    rep("15001-20000", 10)
  )
  
  df$breaks_hist_Toni <- cut(df$contig_length,
                             breaks = breaks_hist_Toni,
                             labels = labels_hist_Toni,
                             include.lowest = TRUE, right = FALSE)
  
  unique_labels <- c("1-3000", "3001-6000", "6001-8500", "8501-10500", "10501-15000", "15001-20000")
  colors <- c("#8e8b8cff", "#a69bb7ff", "#cdaabbff", "#f2bfabff", "#fedcb0ff", "#f7f5b3ff")
  names(colors) <- unique_labels
  
  cur_n_plot <- nrow(df)
  
  gg <- ggplot(df, aes(x = contig_length, fill = breaks_hist_Toni)) +
    geom_histogram(aes(y = ..count..), binwidth = bin_width, colour = "black", boundary = 0) +
    scale_fill_manual(values = colors) +
    geom_density(aes(y = ..density.. * cur_n_plot * bin_width), alpha = 0.1, fill = "#FF6666", adjust = 0.9) +
    labs(title = paste(title, "\nN =", cur_n_plot), x = "Contig length (pb)", y = "Count") +
    xlim(0, 20000) +
    theme_minimal() +
    geom_vline(xintercept = small_threshold, linetype = "dashed", color = "blue", size = 1) +
    geom_vline(xintercept = mini_threshold, linetype = "dashed", color = "red", size = 1) +
    theme(plot.margin = margin(20, 20, 20, 20))  # Add margin around the plot
  
  ggsave(filename = paste0(title, ".svg"), plot = gg, path = output_dir)
  return(gg)
}

# Plot data
boxplot <- plot_boxplot(clans_jgi_merge_df, summary_stats, summary_stats_filtered, n_ranges)
print(boxplot)
ggsave(filename = "Boxplot.svg", plot = boxplot, path = output_dir)

original_plot <- plot_histogram(clans_jgi_merge_df, "All phages\nDistribution of mini phages")
print(original_plot)

# Prepare a list to store plots for the PDF
plots_list <- list(boxplot, original_plot)

# Plot by biome levels
plot_biome_levels <- function(df, plot_histogram_fn) {
  unique_biomes <- unique(na.omit(df$biome1))
  plots <- list()
  for (b1 in unique_biomes) {
    cur_df_b1 <- filter(df, biome1 == b1)
    b1_plot <- plot_histogram_fn(cur_df_b1, paste("Biome Level 1:", b1))
    plots <- append(plots, list(b1_plot))
    
    for (b2 in unique(na.omit(cur_df_b1$biome2))) {
      cur_df_b2 <- filter(cur_df_b1, biome2 == b2)
      b2_plot <- plot_histogram_fn(cur_df_b2, paste("Biome Levels:", b1, "&", b2))
      plots <- append(plots, list(b2_plot))
    }
  }
  return(plots)
}

biome_plots <- plot_biome_levels(clans_jgi_merge_df, plot_histogram)
plots_list <- c(plots_list, biome_plots)

# Filter and save plots after removing unclustered contigs
clans_jgi_merge_df <- clans_jgi_merge_df %>% filter(!(ABSOLUTE_CLUSTER %in% c("VIBRANT_efam_0", "NO_CLUSTER")))
filtered_plot <- plot_histogram(clans_jgi_merge_df, "All phages\nDistribution of mini phages")
print(filtered_plot)
plots_list <- c(plots_list, list(filtered_plot))

filtered_biome_plots <- plot_biome_levels(clans_jgi_merge_df, plot_histogram)
plots_list <- c(plots_list, filtered_biome_plots)

# Create an empty plot as a placeholder
empty_plot <- ggplot() + theme_void()

# Save all plots to a single PDF with vertical orientation and increased whitespace
pdf(file.path(output_dir, "All_Plots_Vertical.pdf"), width = 8.5, height = 11)
for (i in seq(1, length(plots_list), by = 2)) {
  plot1 <- ggplotGrob(plots_list[[i]])
  plot2 <- if (i + 1 <= length(plots_list)) ggplotGrob(plots_list[[i + 1]]) else ggplotGrob(empty_plot)
  
  # Use arrangeGrob with heights to increase space between plots
  grid.arrange(
    arrangeGrob(plot1, plot2, nrow = 2, heights = c(0.45, 0.45)),  # Adjust heights for space between plots
    padding = unit(2, "cm")  # Add padding around the plots
  )
}
dev.off()

# Export results as TSV files
write_tsv(clusters_df, file.path(output_dir, "clusters.tsv"))
write_tsv(per_cluster_df, file.path(output_dir, "per_cluster.tsv"))
write_tsv(clans_jgi_merge_df, file.path(output_dir, "CLANS_JGI_merge.tsv"))

# Additional analysis for biome level statistics
biome_stats_df <- metadata_df %>%
  group_by(biome1, biome2) %>%
  summarise(
    howmany = n(),
    howmanymini = sum(contig_length <= mini_threshold),
    howmanysmall = sum(contig_length <= small_threshold & contig_length > mini_threshold),
    howmanybig = sum(contig_length > small_threshold),
    relativefreqmini = (howmanymini / howmany) * 100,
    relativefreqsmall = (howmanysmall / howmany) * 100
  ) %>%
  filter(howmanymini > 0 | howmanysmall > 0) %>%
  mutate(result = paste(biome1, biome2, howmany, howmanymini, howmanysmall, relativefreqmini, relativefreqsmall))

# Print biome stats dataframe
print(biome_stats_df)
write_tsv(biome_stats_df, file.path(output_dir, "biome_stats_df.tsv"))


# change NA to zeros in clan_jgi_merge_df
clans_jgi_merge_df[is.na(clans_jgi_merge_df)] <- 0
