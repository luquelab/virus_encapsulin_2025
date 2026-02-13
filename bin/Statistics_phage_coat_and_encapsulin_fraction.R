#!/usr/bin/env Rscript

# ------------------------------------------------------------
# Group comparisons for two metrics:
#   1) percentage_phage_coat (positive-only)
#   2) enc_perc_fraction     (positive-only)
#
#   - Phage coat: Welch ANOVA + Games–Howell
#   - Encapsulin: Kruskal–Wallis + Dunn
#
# Outputs (core):
#   - results_phage_coat_pairwise.csv
#   - results_encapsulin_pairwise.csv
#   - summary_by_group_positive_only.csv
#
# Outputs (lookup):
#   - lookup_phage_coat_pairwise.csv
#   - lookup_encapsulin_pairwise.csv
#   - matrix_phage_coat_p.csv
#   - matrix_phage_coat_stars.csv
#   - matrix_encapsulin_p.csv
#   - matrix_encapsulin_stars.csv
# ------------------------------------------------------------

library(tidyverse)

set.seed(1) # This allows reproducibility of bootstrap and permutation results, but can be changed or removed

# -------------------------
# User settings
# -------------------------
csv_path <- "Statistics_phage_coat_and_encapsulin_fraction_data.csv"
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) {
  csv_path <- args[1]
}
# compare only positive values (as requested)
positive_only <- TRUE

# bootstrap/permutation settings (trade-off: speed vs stability)
B_boot <- 2000     # bootstrap replicates for CIs
B_perm <- 4999     # permutation replicates for median test
alpha  <- 0.05
min_n_infer <- 5  # minimum n per group for omnibus + post-hoc inference

# Group order (will keep only those present)
default_group_order <- c("Micro", "Mini", "Control", "Fringe", "Refseq", "PICIs",
                         "ENC_1", "ENC_2", "ENC_3", "ENC_4")

# -------------------------
# Helper: safe package availability
# -------------------------
has_pkg <- function(pkg) {
  requireNamespace(pkg, quietly = TRUE)
}

# -------------------------
# Helper: pretty section prints
# -------------------------
banner <- function(txt) {
  cat("\n", strrep("=", 78), "\n", sep = "")
  cat(txt, "\n")
  cat(strrep("=", 78), "\n", sep = "")
}

subbanner <- function(txt) {
  cat("\n", strrep("-", 78), "\n", sep = "")
  cat(txt, "\n")
  cat(strrep("-", 78), "\n", sep = "")
}

# -------------------------
# Helpers: p-value labels + delta magnitude
# -------------------------
p_to_stars <- function(p) {
  ifelse(is.na(p), NA_character_,
         ifelse(p < 1e-4, "****",
                ifelse(p < 1e-3, "***",
                       ifelse(p < 1e-2, "**",
                              ifelse(p < 0.05, "*", "ns")))))
}

# Common rule-of-thumb bins for Cliff's delta magnitude
# (Romano et al. style thresholds)
delta_magnitude <- function(delta) {
  if (is.na(delta)) return(NA_character_)
  ad <- abs(delta)
  if (ad < 0.147) return("negligible")
  if (ad < 0.330) return("small")
  if (ad < 0.474) return("medium")
  "large"
}

# Make an order-invariant pair_id so joins never fail due to swapped order
make_pair_id <- function(a, b) {
  paste(pmin(a, b), pmax(a, b), sep = "__")
}

# -------------------------
# Helper: Cliff's delta (effect size)
# delta = P(x > y) - P(x < y)
# using Mann-Whitney U:
# delta = (2U)/(n_x*n_y) - 1
# -------------------------
cliffs_delta <- function(x, y) {
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  nx <- length(x); ny <- length(y)
  if (nx < 1 || ny < 1) return(NA_real_)
  r <- rank(c(x, y), ties.method = "average")
  rx <- sum(r[1:nx])
  U  <- rx - nx * (nx + 1) / 2
  (2 * U) / (nx * ny) - 1
}

# Bootstrap CI for a statistic on two samples
boot_diff_ci <- function(x, y, stat = c("mean", "median", "trimmed_mean"),
                         trim = 0.2, B = 2000, conf = 0.95) {
  stat <- match.arg(stat)
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  nx <- length(x); ny <- length(y)
  if (nx < 2 || ny < 2) {
    return(tibble(est = NA_real_, lo = NA_real_, hi = NA_real_))
  }
  f <- switch(
    stat,
    mean = mean,
    median = median,
    trimmed_mean = function(z) mean(z, trim = trim)
  )
  est <- f(x) - f(y)
  boots <- replicate(B, {
    xb <- sample(x, nx, replace = TRUE)
    yb <- sample(y, ny, replace = TRUE)
    f(xb) - f(yb)
  })
  a <- (1 - conf) / 2
  lo <- unname(quantile(boots, probs = a, na.rm = TRUE))
  hi <- unname(quantile(boots, probs = 1 - a, na.rm = TRUE))
  tibble(est = est, lo = lo, hi = hi)
}

# Bootstrap CI for Cliff's delta
boot_delta_ci <- function(x, y, B = 2000, conf = 0.95) {
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  nx <- length(x); ny <- length(y)
  if (nx < 2 || ny < 2) {
    return(tibble(delta = NA_real_, lo = NA_real_, hi = NA_real_))
  }
  est <- cliffs_delta(x, y)
  boots <- replicate(B, {
    xb <- sample(x, nx, replace = TRUE)
    yb <- sample(y, ny, replace = TRUE)
    cliffs_delta(xb, yb)
  })
  a <- (1 - conf) / 2
  lo <- unname(quantile(boots, probs = a, na.rm = TRUE))
  hi <- unname(quantile(boots, probs = 1 - a, na.rm = TRUE))
  tibble(delta = est, lo = lo, hi = hi)
}

# Permutation test for difference in medians (two-sided)
perm_median_p <- function(x, y, B = 4999) {
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  nx <- length(x); ny <- length(y)
  if (nx < 2 || ny < 2) return(NA_real_)
  obs <- median(x) - median(y)
  z <- c(x, y)
  cnt <- 0L
  for (b in seq_len(B)) {
    idx <- sample.int(nx + ny, nx, replace = FALSE)
    xb <- z[idx]
    yb <- z[-idx]
    stat <- median(xb) - median(yb)
    if (abs(stat) >= abs(obs)) cnt <- cnt + 1L
  }
  (cnt + 1) / (B + 1)  # add-one correction
}

# -------------------------
# Read and validate data
# -------------------------
banner("1) Reading data and validating columns")

if (!file.exists(csv_path)) {
  stop("CSV not found: ", csv_path,
       "\nPut it next to this script, or change csv_path at the top.")
}

data <- read.csv(csv_path, stringsAsFactors = FALSE)
cat("Loaded rows:", nrow(data), " | columns:", ncol(data), "\n")

required_cols <- c("group", "percentage_phage_coat", "enc_perc_fraction")
missing <- setdiff(required_cols, names(data))
if (length(missing) > 0) {
  stop("Missing required columns: ", paste(missing, collapse = ", "))
}

# Clean group and enforce order
present_groups <- sort(unique(data$group))
cat("Groups present (raw):", paste(present_groups, collapse = ", "), "\n")

group_order <- default_group_order[default_group_order %in% present_groups]
unknown_groups <- setdiff(present_groups, default_group_order)
if (length(unknown_groups) > 0) {
  cat("NOTE: Found groups not in default_group_order:\n  ",
      paste(unknown_groups, collapse = ", "), "\n")
  cat("They will be kept but placed after the default ordered groups.\n")
  group_order <- c(group_order, unknown_groups)
}

data <- data %>%
  mutate(group = factor(group, levels = group_order))

cat("Final group order used:\n  ", paste(levels(data$group), collapse = "  |  "), "\n")

subbanner("Quick counts by group (all rows)")
print(data %>% count(group, name = "n_total") %>% arrange(desc(n_total)))

# -------------------------
# Summary function
# -------------------------
summarize_by_group <- function(df, metric, positive_only = TRUE) {
  metric_sym <- rlang::sym(metric)
  out <- df %>%
    filter(!is.na(!!metric_sym)) %>%
    { if (positive_only) filter(., !!metric_sym > 0) else . } %>%
    group_by(group) %>%
    summarise(
      n = n(),
      mean = mean(!!metric_sym),
      median = median(!!metric_sym),
      sd = sd(!!metric_sym),
      iqr = IQR(!!metric_sym),
      min = min(!!metric_sym),
      max = max(!!metric_sym),
      pct_at_100 = mean((!!metric_sym) == 100) * 100,
      .groups = "drop"
    ) %>%
    complete(group = factor(levels(df$group), levels = levels(df$group)),
             fill = list(n = 0, mean = NA_real_, median = NA_real_, sd = NA_real_,
                         iqr = NA_real_, min = NA_real_, max = NA_real_, pct_at_100 = NA_real_)) %>%
    arrange(group)
  out
}

# -------------------------
# Pairwise engine (effect sizes + CIs + permutation medians)
# -------------------------
pairwise_table <- function(df, metric, positive_only = TRUE,
                           do_perm_median = TRUE,
                           B_boot = 2000, B_perm = 4999) {
  
  metric_sym <- rlang::sym(metric)
  
  df2 <- df %>%
    filter(!is.na(!!metric_sym)) %>%
    { if (positive_only) filter(., !!metric_sym > 0) else . } %>%
    mutate(group = droplevels(group))
  
  # Diagnostics
  cat("Metric:", metric, " | positive_only:", positive_only, "\n")
  cat("Rows used after filtering:", nrow(df2), "\n")
  
  counts <- df2 %>% count(group, name = "n") %>% arrange(desc(n))
  print(counts)
  
  tiny <- counts %>% filter(n > 0 & n < 5)
  if (nrow(tiny) > 0) {
    cat("WARNING: Some groups have very small n (<5) after filtering:\n")
    print(tiny)
  }
  
  var_tbl <- df2 %>%
    group_by(group) %>%
    summarise(var = var(!!metric_sym), sd = sd(!!metric_sym), .groups = "drop")
  zv <- var_tbl %>% filter(is.na(var) | var == 0)
  if (nrow(zv) > 0) {
    cat("WARNING: Zero (or NA) variance groups detected for this metric:\n")
    print(zv)
    cat("This can break some parametric post-hoc procedures; rank/permutation methods remain OK.\n")
  }
  
  if (metric == "enc_perc_fraction") {
    frac100 <- mean(df2[[metric]] == 100) * 100
    cat(sprintf("Info: %.2f%% of positive observations are exactly 100 for enc_perc_fraction.\n", frac100))
  }
  
  grps <- levels(df2$group)
  pairs <- combn(grps, 2, simplify = FALSE)
  
  tab <- map_dfr(pairs, function(p) {
    g1 <- p[1]; g2 <- p[2]
    x <- df2 %>% filter(group == g1) %>% pull(!!metric_sym)
    y <- df2 %>% filter(group == g2) %>% pull(!!metric_sym)
    
    n1 <- length(x); n2 <- length(y)
    
    if (n1 == 0 || n2 == 0) {
      return(tibble(
        group1 = g1, group2 = g2,
        pair_id = make_pair_id(g1, g2),
        n1 = n1, n2 = n2,
        mean1 = NA_real_, mean2 = NA_real_,
        median1 = NA_real_, median2 = NA_real_,
        diff_mean = NA_real_, diff_mean_lo = NA_real_, diff_mean_hi = NA_real_,
        diff_median = NA_real_, diff_median_lo = NA_real_, diff_median_hi = NA_real_,
        cliffs_delta = NA_real_, delta_lo = NA_real_, delta_hi = NA_real_,
        prob_superiority = NA_real_,
        p_perm_median = NA_real_
      ))
    }
    
    ci_mean   <- boot_diff_ci(x, y, stat = "mean", B = B_boot)
    ci_median <- boot_diff_ci(x, y, stat = "median", B = B_boot)
    ci_delta  <- boot_delta_ci(x, y, B = B_boot)
    pperm     <- if (do_perm_median) perm_median_p(x, y, B = B_perm) else NA_real_
    
    delta <- ci_delta$delta
    prob_sup <- ifelse(is.na(delta), NA_real_, (delta + 1) / 2)
    
    tibble(
      group1 = g1, group2 = g2,
      pair_id = make_pair_id(g1, g2),
      n1 = n1, n2 = n2,
      mean1 = mean(x), mean2 = mean(y),
      median1 = median(x), median2 = median(y),
      diff_mean = ci_mean$est, diff_mean_lo = ci_mean$lo, diff_mean_hi = ci_mean$hi,
      diff_median = ci_median$est, diff_median_lo = ci_median$lo, diff_median_hi = ci_median$hi,
      cliffs_delta = ci_delta$delta, delta_lo = ci_delta$lo, delta_hi = ci_delta$hi,
      prob_superiority = prob_sup,
      p_perm_median = pperm
    )
  })
  
  tab %>%
    mutate(
      p_perm_median_BH   = p.adjust(p_perm_median, method = "BH"),
      p_perm_median_Holm = p.adjust(p_perm_median, method = "holm")
    )
}

# -------------------------
# Test wrappers: Welch ANOVA, Kruskal-Wallis, Post-hoc
# -------------------------
filter_for_inference <- function(df2, metric, min_n_infer) {
  counts <- df2 %>% count(group, name = "n") %>% arrange(desc(n))
  keep <- counts %>% filter(n >= min_n_infer) %>% pull(group)
  
  dropped <- counts %>% filter(n < min_n_infer)
  if (nrow(dropped) > 0) {
    cat("INFO: Excluding groups from inferential tests due to small n <", min_n_infer, ":\n", sep = "")
    print(dropped)
  }
  
  df2 %>% filter(group %in% keep) %>% droplevels()
}

run_welch_anova <- function(df, metric, positive_only = TRUE, min_n_infer = 5) {
  metric_sym <- rlang::sym(metric)
  df2 <- df %>%
    filter(!is.na(!!metric_sym)) %>%
    { if (positive_only) filter(., !!metric_sym > 0) else . } %>%
    mutate(group = droplevels(group))
  
  df2i <- filter_for_inference(df2, metric, min_n_infer)
  
  if (nlevels(df2i$group) < 2) {
    cat("Welch ANOVA skipped: <2 groups with n>=", min_n_infer, "\n", sep = "")
    return(NULL)
  }
  
  oneway.test(df2i[[metric]] ~ df2i$group, var.equal = FALSE)
}

run_kruskal <- function(df, metric, positive_only = TRUE, min_n_infer = 5) {
  metric_sym <- rlang::sym(metric)
  df2 <- df %>%
    filter(!is.na(!!metric_sym)) %>%
    { if (positive_only) filter(., !!metric_sym > 0) else . } %>%
    mutate(group = droplevels(group))
  
  df2i <- filter_for_inference(df2, metric, min_n_infer)
  
  if (nlevels(df2i$group) < 2) {
    cat("Kruskal-Wallis skipped: <2 groups with n>=", min_n_infer, "\n", sep = "")
    return(NULL)
  }
  
  kruskal.test(df2i[[metric]] ~ df2i$group)
}

posthoc_games_howell <- function(df, metric, positive_only = TRUE, min_n_infer = 5) {
  metric_sym <- rlang::sym(metric)
  df2 <- df %>%
    filter(!is.na(!!metric_sym)) %>%
    { if (positive_only) filter(., !!metric_sym > 0) else . } %>%
    mutate(group = droplevels(group))
  
  df2i <- filter_for_inference(df2, metric, min_n_infer)
  if (nlevels(df2i$group) < 2) return(NULL)
  
  if (has_pkg("rstatix")) {
    suppressPackageStartupMessages(library(rstatix))
    out <- df2i %>%
      rstatix::games_howell_test(as.formula(paste(metric, "~ group"))) %>%
      transmute(
        group1, group2,
        pair_id = make_pair_id(group1, group2),
        estimate = estimate,
        conf.low = conf.low,
        conf.high = conf.high,
        p_primary = p.adj,
        primary_method = "Games-Howell (rstatix) [p.adj returned]"
      )
    return(out)
  }
  
  cat("NOTE: rstatix not available -> using pairwise Welch t-tests (Holm) as fallback.\n")
  pw <- pairwise.t.test(df2i[[metric]], df2i$group, pool.sd = FALSE, p.adjust.method = "holm")
  m <- pw$p.value
  
  as.data.frame(as.table(m)) %>%
    as_tibble() %>%
    filter(!is.na(Freq)) %>%
    transmute(
      group1 = as.character(Var1),
      group2 = as.character(Var2),
      pair_id = make_pair_id(group1, group2),
      estimate = NA_real_, conf.low = NA_real_, conf.high = NA_real_,
      p_primary = Freq,
      primary_method = "Pairwise Welch t-test (Holm-adjusted) [fallback]"
    )
}

posthoc_dunn <- function(df, metric, positive_only = TRUE, p_adjust = "BH", min_n_infer = 5) {
  metric_sym <- rlang::sym(metric)
  df2 <- df %>%
    filter(!is.na(!!metric_sym)) %>%
    { if (positive_only) filter(., !!metric_sym > 0) else . } %>%
    mutate(group = droplevels(group))
  
  df2i <- filter_for_inference(df2, metric, min_n_infer)
  if (nlevels(df2i$group) < 2) return(NULL)
  
  if (has_pkg("rstatix")) {
    suppressPackageStartupMessages(library(rstatix))
    out <- df2i %>%
      rstatix::dunn_test(as.formula(paste(metric, "~ group")), p.adjust.method = p_adjust) %>%
      transmute(
        group1, group2,
        pair_id = make_pair_id(group1, group2),
        Z = statistic,
        p_primary = p.adj,
        primary_method = paste0("Dunn (rstatix) + ", p_adjust, " [p.adj returned]")
      )
    return(out)
  }
  
  if (has_pkg("FSA")) {
    suppressPackageStartupMessages(library(FSA))
    dt <- FSA::dunnTest(df2i[[metric]] ~ df2i$group, method = p_adjust)
    out <- dt$res %>%
      as_tibble() %>%
      separate(Comparison, into = c("group1", "group2"), sep = " - ") %>%
      transmute(
        group1, group2,
        pair_id = make_pair_id(group1, group2),
        Z = Z,
        p_primary = P.adj,
        primary_method = paste0("Dunn (FSA) + ", p_adjust, " [p.adj returned]")
      )
    return(out)
  }
  
  cat("NOTE: Dunn not available -> rely on permutation median tests.\n")
  NULL
}

# -------------------------
#Build easy "lookup" table + matrix exports
# -------------------------
make_lookup_table <- function(pw_out, metric, alpha = 0.05) {
  pw_out %>%
    mutate(
      metric = metric,
      p = p_primary,
      p_stars = p_to_stars(p),
      significant = !is.na(p) & p < alpha,
      delta_mag = map_chr(cliffs_delta, delta_magnitude),
      higher_mean_group   = ifelse(is.na(diff_mean), NA_character_, ifelse(diff_mean > 0, group1, group2)),
      higher_median_group = ifelse(is.na(diff_median), NA_character_, ifelse(diff_median > 0, group1, group2))
    ) %>%
    transmute(
      metric,
      group1, group2, n1, n2,
      mean1, mean2, median1, median2,
      higher_mean_group, higher_median_group,
      diff_mean, diff_mean_lo, diff_mean_hi,
      diff_median, diff_median_lo, diff_median_hi,
      cliffs_delta, delta_lo, delta_hi,
      prob_superiority, delta_mag,
      p_primary = p,
      p_stars,
      significant,
      primary_method,
      # keep permutation track as a sensitivity analysis column:
      p_perm_median,
      p_perm_median_BH,
      p_perm_median_Holm
    ) %>%
    arrange(p_primary)
}

write_matrix_outputs <- function(pw_out, groups, value_col, out_csv) {
  # symmetric long
  v <- rlang::sym(value_col)
  
  sym <- bind_rows(
    pw_out %>% transmute(row = group1, col = group2, value = !!v),
    pw_out %>% transmute(row = group2, col = group1, value = !!v)
  ) %>%
    mutate(
      row = factor(row, levels = groups),
      col = factor(col, levels = groups)
    )
  
  # include diagonal explicitly
  diag <- tibble(
    row = factor(groups, levels = groups),
    col = factor(groups, levels = groups),
    value = NA
  )
  
  mat <- bind_rows(sym, diag) %>%
    distinct(row, col, .keep_all = TRUE) %>%
    arrange(row, col) %>%
    pivot_wider(names_from = col, values_from = value)
  
  write.csv(mat, out_csv, row.names = FALSE)
}

# -------------------------
# Main analysis runner per metric
# -------------------------
analyze_metric <- function(df, metric, metric_label, primary_strategy,
                           positive_only = TRUE,
                           do_perm_median = TRUE,
                           B_boot = 2000, B_perm = 4999) {
  
  banner(paste0("2) Analysis for ", metric_label, "  [metric=", metric, "]"))
  
  subbanner("Summary by group (positive-only filter applied here)")
  summ <- summarize_by_group(df, metric, positive_only = positive_only)
  print(summ)
  
  subbanner("Omnibus tests (run both as diagnostics)")
  welch <- run_welch_anova(df, metric, positive_only = positive_only, min_n_infer = min_n_infer)
  kw    <- run_kruskal(df, metric, positive_only = positive_only, min_n_infer = min_n_infer)
  
  if (!is.null(welch)) {
    cat("\nWelch ANOVA:\n")
    print(welch)
  } else cat("Welch ANOVA: not run (need >=2 groups with data)\n")
  
  if (!is.null(kw)) {
    cat("\nKruskal-Wallis:\n")
    print(kw)
  } else cat("Kruskal-Wallis: not run (need >=2 groups with data)\n")
  
  subbanner("Pairwise effect sizes + bootstrap CIs + permutation median tests")
  pw_stats <- pairwise_table(df, metric,
                             positive_only = positive_only,
                             do_perm_median = do_perm_median,
                             B_boot = B_boot, B_perm = B_perm)
  
  cat("\nPreview (top 10 rows) of pairwise stats:\n")
  print(pw_stats %>% slice(1:10))
  
  subbanner(paste0("Primary post-hoc strategy: ", primary_strategy))
  
  # IMPORTANT: primary p-values returned by GH/Dunn are already adjusted (p.adj).
  # We do NOT re-adjust them (avoids double-correction).
  if (primary_strategy == "welch + games-howell") {
    gh <- posthoc_games_howell(df, metric, positive_only = positive_only, min_n_infer = min_n_infer)
    cat("\nPost-hoc results head:\n")
    print(head(gh, 12))
    
    pw_out <- pw_stats %>%
      left_join(gh %>% dplyr::select(pair_id, p_primary, primary_method), by = "pair_id")
  }
  
  if (primary_strategy == "kruskal + dunn") {
    dunn <- posthoc_dunn(df, metric, positive_only = positive_only, p_adjust = "BH", min_n_infer = min_n_infer)
    
    if (!is.null(dunn)) {
      cat("\nDunn post-hoc results head:\n")
      print(head(dunn, 12))
      
      pw_out <- pw_stats %>%
        left_join(dunn %>% dplyr::select(pair_id, p_primary, primary_method), by = "pair_id")
    } else {
      cat("\nNo Dunn available -> using permutation median p-values as primary.\n")
      pw_out <- pw_stats %>%
        mutate(
          p_primary = p_perm_median,
          primary_method = "Permutation test on median difference (primary)"
        )
    }
  }
  
  # Flag significant comparisons (using primary p-values)
  pw_out <- pw_out %>%
    mutate(
      sig_primary = !is.na(p_primary) & p_primary < alpha,
      p_primary_stars = p_to_stars(p_primary),
      delta_mag = map_chr(cliffs_delta, delta_magnitude)
    )
  
  subbanner("Top significant comparisons (if any) by primary p-value")
  top_sig <- pw_out %>%
    filter(sig_primary) %>%
    arrange(p_primary) %>%
    select(group1, group2, n1, n2,
           diff_mean, diff_mean_lo, diff_mean_hi,
           diff_median, diff_median_lo, diff_median_hi,
           cliffs_delta, delta_lo, delta_hi, prob_superiority, delta_mag,
           p_primary, p_primary_stars, primary_method,
           p_perm_median_BH) %>%
    slice(1:15)
  
  if (nrow(top_sig) == 0) {
    cat("No significant pairwise differences detected under primary method at alpha=", alpha, "\n", sep="")
  } else {
    print(top_sig)
  }
  
  list(summary = summ, pairwise = pw_out, welch = welch, kw = kw)
}

# -------------------------
# Run analyses (positive-only)
# -------------------------
banner("2) Running the full analysis pipeline")

res_phage <- analyze_metric(
  df = data,
  metric = "percentage_phage_coat",
  metric_label = "Phage coat (%)",
  primary_strategy = "welch + games-howell",
  positive_only = positive_only,
  do_perm_median = TRUE,
  B_boot = B_boot,
  B_perm = B_perm
)

res_enc <- analyze_metric(
  df = data,
  metric = "enc_perc_fraction",
  metric_label = "Encapsulin (%) among positives",
  primary_strategy = "kruskal + dunn",
  positive_only = positive_only,
  do_perm_median = TRUE,
  B_boot = B_boot,
  B_perm = B_perm
)

# -------------------------
# Save outputs
# -------------------------
banner("3) Writing output CSV files")

summary_out <- full_join(
  res_phage$summary %>% rename_with(~ paste0("phage_", .x), -group),
  res_enc$summary   %>% rename_with(~ paste0("enc_", .x), -group),
  by = "group"
)

write.csv(summary_out, "summary_by_group_positive_only.csv", row.names = FALSE)
write.csv(res_phage$pairwise, "results_phage_coat_pairwise.csv", row.names = FALSE)
write.csv(res_enc$pairwise, "results_encapsulin_pairwise.csv", row.names = FALSE)

# NEW: easy lookup tables
lookup_phage <- make_lookup_table(res_phage$pairwise, metric = "percentage_phage_coat", alpha = alpha)
lookup_enc   <- make_lookup_table(res_enc$pairwise, metric = "enc_perc_fraction", alpha = alpha)

write.csv(lookup_phage, "lookup_phage_coat_pairwise.csv", row.names = FALSE)
write.csv(lookup_enc,   "lookup_encapsulin_pairwise.csv", row.names = FALSE)

# NEW: quick p-value + stars matrices
groups_used <- levels(data$group)

# phage coat matrices
write_matrix_outputs(res_phage$pairwise, groups_used, "p_primary", "matrix_phage_coat_p.csv")
write_matrix_outputs(res_phage$pairwise, groups_used, "p_primary_stars", "matrix_phage_coat_stars.csv")

# encapsulin matrices
write_matrix_outputs(res_enc$pairwise, groups_used, "p_primary", "matrix_encapsulin_p.csv")
write_matrix_outputs(res_enc$pairwise, groups_used, "p_primary_stars", "matrix_encapsulin_stars.csv")

cat("Wrote:\n",
    " - summary_by_group_positive_only.csv\n",
    " - results_phage_coat_pairwise.csv\n",
    " - results_encapsulin_pairwise.csv\n",
    " - lookup_phage_coat_pairwise.csv\n",
    " - lookup_encapsulin_pairwise.csv\n",
    " - matrix_phage_coat_p.csv\n",
    " - matrix_phage_coat_stars.csv\n",
    " - matrix_encapsulin_p.csv\n",
    " - matrix_encapsulin_stars.csv\n", sep = "")

# -------------------------
# Final guidance prints (for your revision narrative)
# -------------------------
banner("4) Suggested Methods wording (draft)")

cat(
  "Draft Methods language (edit freely):\n\n",
  "- All analyses were restricted to proteins with positive metric values (>0).\n",
  "- For phage coat coverage, we tested for differences among groups using Welch’s one-way ANOVA\n",
  "  (heteroscedasticity-robust). Pairwise comparisons were performed using the Games–Howell\n",
  "  procedure when available; otherwise, pairwise Welch t-tests with Holm-adjusted p-values were used.\n",
  "- For encapsulin fraction among positives, we used Kruskal–Wallis rank-sum tests followed by\n",
  "  Dunn’s post-hoc comparisons with Benjamini–Hochberg-adjusted p-values when available; as a\n",
  "  distribution-free sensitivity analysis we also computed permutation tests for pairwise\n",
  "  median differences.\n",
  "- In addition to p-values, we report effect sizes (Cliff’s delta / probability of superiority)\n",
  "  and bootstrap confidence intervals for mean and median differences.\n\n",
  sep = ""
)

banner("DONE")
#END

