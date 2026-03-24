#!/usr/bin/env Rscript

# BHPMF phylogeny-aware gap filling script for data.csv
# - Uses AccSpeciesName1 + species_phylo_tree.nwk to build hierarchy info.
# - Completes LDMC and Leafph via BHPMF.
# - Adds uncertainty columns (posterior SD + 95% CI columns).
# - Outputs grouped (FROM) and ungrouped figures.

suppressPackageStartupMessages({
  library(ape)
  library(ggplot2)
})

source_local_bhpmf <- function(repo_root) {
  r_files <- c(
    "R/countNumTraits.R",
    "R/generatePermutation.R",
    "R/splitData.R",
    "R/findFold.R",
    "R/findPhyloInfo.R",
    "R/buildUpperLevelMat.R",
    "R/preprocess_cv.R",
    "R/plot_rmse_std.R",
    "R/tune_BHPMF.R",
    "R/gap_filling.R"
  )
  for (f in r_files) {
    source(file.path(repo_root, f))
  }
}

normalize_species <- function(x) {
  x <- tolower(x)
  x <- gsub("_", " ", x)
  x <- gsub("[^a-z0-9 ]", " ", x)
  x <- gsub("\\s+", " ", x)
  trimws(x)
}

choose_cut_levels <- function(n_species) {
  # 2 tree-derived clade levels + species level.
  k_high <- max(8, floor(sqrt(n_species)))
  k_high <- min(k_high, n_species)
  k_low <- max(3, floor(k_high / 3))
  c(k_low = k_low, k_high = k_high)
}

build_hierarchy_from_tree <- function(df, tree) {
  sp_data_raw <- as.character(df$AccSpeciesName1)
  sp_data_norm <- normalize_species(sp_data_raw)

  tip_raw <- tree$tip.label
  tip_norm <- normalize_species(tip_raw)

  # Map normalized species names to original tree tip label.
  tip_map <- setNames(tip_raw, tip_norm)

  matched <- sp_data_norm %in% names(tip_map)

  tip_labels_for_data <- ifelse(matched, tip_map[sp_data_norm], NA)

  # Build clade assignments from tree for matched species.
  hc <- as.hclust.phylo(tree)
  k <- choose_cut_levels(length(tree$tip.label))
  clade_low <- cutree(hc, k = k[["k_low"]])
  clade_high <- cutree(hc, k = k[["k_high"]])

  # Names are tree tip labels.
  clade_low_map <- setNames(clade_low, names(clade_low))
  clade_high_map <- setNames(clade_high, names(clade_high))

  # Prepare fallback genus-based clustering for unmatched species.
  genus <- sapply(strsplit(sp_data_norm, " "), function(v) if (length(v) > 0) v[1] else "unknown")
  genus[is.na(genus) | genus == ""] <- "unknown"

  # Map clades.
  clade_l1 <- character(nrow(df))
  clade_l2 <- character(nrow(df))
  species_level <- character(nrow(df))

  for (i in seq_len(nrow(df))) {
    if (!is.na(tip_labels_for_data[i])) {
      tip_i <- tip_labels_for_data[i]
      clade_l1[i] <- sprintf("tree_l1_%03d", clade_low_map[[tip_i]])
      clade_l2[i] <- sprintf("tree_l2_%03d", clade_high_map[[tip_i]])
      species_level[i] <- sprintf("sp_%s", tip_i)
    } else {
      # If species not found in tree, keep deterministic fallback levels.
      clade_l1[i] <- sprintf("fallback_genus_%s", genus[i])
      clade_l2[i] <- sprintf("fallback_species_%s", gsub(" ", "_", sp_data_norm[i]))
      species_level[i] <- sprintf("sp_%s", gsub(" ", "_", sp_data_norm[i]))
    }
  }

  hierarchy <- data.frame(
    ObsID = seq_len(nrow(df)),
    CladeLevel1 = clade_l1,
    CladeLevel2 = clade_l2,
    Species = species_level,
    stringsAsFactors = FALSE
  )

  list(
    hierarchy = as.matrix(hierarchy),
    matched_ratio = mean(matched),
    unmatched_species = unique(sp_data_raw[!matched])
  )
}

plot_outputs <- function(df_out, out_dir) {
  dir.create(file.path(out_dir, "plots"), showWarnings = FALSE, recursive = TRUE)

  # Ungrouped: observed vs imputed density (all rows)
  long_all <- rbind(
    data.frame(Trait = "LDMC", Type = "Observed", Value = df_out$LDMC),
    data.frame(Trait = "LDMC", Type = "ImputedMean", Value = df_out$LDMC_bhpmf_mean),
    data.frame(Trait = "Leafph", Type = "Observed", Value = df_out$Leafph),
    data.frame(Trait = "Leafph", Type = "ImputedMean", Value = df_out$Leafph_bhpmf_mean)
  )
  long_all <- long_all[!is.na(long_all$Value), ]

  p_density <- ggplot(long_all, aes(x = Value, color = Type, fill = Type)) +
    geom_density(alpha = 0.2) +
    facet_wrap(~Trait, scales = "free") +
    theme_bw() +
    labs(title = "Ungrouped distribution: observed vs BHPMF mean")
  ggsave(file.path(out_dir, "plots", "ungrouped_density_observed_vs_imputed.png"), p_density, width = 10, height = 5, dpi = 150)

  # Ungrouped uncertainty distribution
  unc_long <- rbind(
    data.frame(Trait = "LDMC", SD = df_out$LDMC_bhpmf_sd),
    data.frame(Trait = "Leafph", SD = df_out$Leafph_bhpmf_sd)
  )
  unc_long <- unc_long[!is.na(unc_long$SD), ]
  p_unc <- ggplot(unc_long, aes(x = SD)) +
    geom_histogram(bins = 40, fill = "#2c7fb8", alpha = 0.8) +
    facet_wrap(~Trait, scales = "free") +
    theme_bw() +
    labs(title = "Ungrouped uncertainty (posterior SD)", x = "Posterior SD")
  ggsave(file.path(out_dir, "plots", "ungrouped_uncertainty_hist.png"), p_unc, width = 10, height = 5, dpi = 150)

  # Grouped by FROM: imputed mean boxplot
  group_long <- rbind(
    data.frame(FROM = df_out$FROM, Trait = "LDMC", Value = df_out$LDMC_bhpmf_mean),
    data.frame(FROM = df_out$FROM, Trait = "Leafph", Value = df_out$Leafph_bhpmf_mean)
  )
  group_long <- group_long[!is.na(group_long$Value), ]

  p_group_mean <- ggplot(group_long, aes(x = FROM, y = Value, fill = FROM)) +
    geom_boxplot(outlier.alpha = 0.2) +
    facet_wrap(~Trait, scales = "free") +
    theme_bw() +
    labs(title = "Grouped by FROM: BHPMF imputed mean")
  ggsave(file.path(out_dir, "plots", "grouped_from_imputed_mean_boxplot.png"), p_group_mean, width = 9, height = 5, dpi = 150)

  # Grouped by FROM: uncertainty boxplot
  group_unc <- rbind(
    data.frame(FROM = df_out$FROM, Trait = "LDMC", SD = df_out$LDMC_bhpmf_sd),
    data.frame(FROM = df_out$FROM, Trait = "Leafph", SD = df_out$Leafph_bhpmf_sd)
  )
  group_unc <- group_unc[!is.na(group_unc$SD), ]

  p_group_unc <- ggplot(group_unc, aes(x = FROM, y = SD, fill = FROM)) +
    geom_boxplot(outlier.alpha = 0.2) +
    facet_wrap(~Trait, scales = "free") +
    theme_bw() +
    labs(title = "Grouped by FROM: uncertainty (posterior SD)")
  ggsave(file.path(out_dir, "plots", "grouped_from_uncertainty_boxplot.png"), p_group_unc, width = 9, height = 5, dpi = 150)

  # Grouped + ungrouped scatter on observed entries only.
  obs_idx_ldmc <- !is.na(df_out$LDMC)
  obs_idx_leafph <- !is.na(df_out$Leafph)

  scatter_df <- rbind(
    data.frame(Trait = "LDMC", FROM = df_out$FROM[obs_idx_ldmc],
               Observed = df_out$LDMC[obs_idx_ldmc],
               Imputed = df_out$LDMC_bhpmf_mean[obs_idx_ldmc]),
    data.frame(Trait = "Leafph", FROM = df_out$FROM[obs_idx_leafph],
               Observed = df_out$Leafph[obs_idx_leafph],
               Imputed = df_out$Leafph_bhpmf_mean[obs_idx_leafph])
  )

  p_scatter_group <- ggplot(scatter_df, aes(x = Observed, y = Imputed, color = FROM)) +
    geom_point(alpha = 0.4, size = 1) +
    geom_abline(slope = 1, intercept = 0, linetype = 2) +
    facet_wrap(~Trait, scales = "free") +
    theme_bw() +
    labs(title = "Observed vs imputed (grouped by FROM)")
  ggsave(file.path(out_dir, "plots", "grouped_from_observed_vs_imputed_scatter.png"), p_scatter_group, width = 10, height = 5, dpi = 150)

  p_scatter_all <- ggplot(scatter_df, aes(x = Observed, y = Imputed)) +
    geom_point(alpha = 0.35, color = "#1f78b4", size = 1) +
    geom_abline(slope = 1, intercept = 0, linetype = 2) +
    facet_wrap(~Trait, scales = "free") +
    theme_bw() +
    labs(title = "Observed vs imputed (ungrouped)")
  ggsave(file.path(out_dir, "plots", "ungrouped_observed_vs_imputed_scatter.png"), p_scatter_all, width = 10, height = 5, dpi = 150)
}

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  repo_root <- ifelse(length(args) >= 1, args[[1]], ".")
  out_dir <- ifelse(length(args) >= 2, args[[2]], file.path(repo_root, "outputs_bhpmf"))

  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  data_path <- file.path(repo_root, "data.csv")
  tree_path <- file.path(repo_root, "species_phylo_tree.nwk")

  if (!file.exists(data_path)) stop("data.csv not found: ", data_path)
  if (!file.exists(tree_path)) stop("species_phylo_tree.nwk not found: ", tree_path)

  source_local_bhpmf(repo_root)

  df <- read.csv(data_path, stringsAsFactors = FALSE, check.names = FALSE)

  # Traits to fill.
  trait_cols <- c("LDMC", "Leafph")
  for (tc in trait_cols) {
    df[[tc]] <- as.numeric(df[[tc]])
  }

  tree <- read.tree(tree_path)

  h <- build_hierarchy_from_tree(df, tree)
  hierarchy_info <- h$hierarchy

  # Save a diagnostic for species matching.
  writeLines(
    c(
      sprintf("Tree match ratio: %.4f", h$matched_ratio),
      sprintf("Unmatched species count: %d", length(h$unmatched_species)),
      "Unmatched species:",
      if (length(h$unmatched_species) > 0) paste0("- ", h$unmatched_species) else "- none"
    ),
    con = file.path(out_dir, "species_tree_match_report.txt")
  )

  X <- as.matrix(df[, trait_cols])

  mean_path <- file.path(out_dir, "bhpmf_mean_imputed.tsv")
  sd_path <- file.path(out_dir, "bhpmf_sd_imputed.tsv")
  tmp_dir <- file.path(out_dir, "tmp")
  dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)

  # Core BHPMF run
  GapFilling(
    X = X,
    hierarchy.info = hierarchy_info,
    prediction.level = ncol(hierarchy_info),
    used.num.hierarchy.levels = ncol(hierarchy_info) - 1,
    num.samples = 1200,
    burn = 200,
    gaps = 2,
    num.latent.feats = 10,
    tuning = FALSE,
    tmp.dir = tmp_dir,
    mean.gap.filled.output.path = mean_path,
    std.gap.filled.output.path = sd_path,
    rmse.plot.test.data = TRUE,
    verbose = TRUE
  )

  mean_imp <- read.table(mean_path, header = TRUE, sep = "\t", check.names = FALSE)
  sd_imp <- read.table(sd_path, header = TRUE, sep = "\t", check.names = FALSE)

  # New columns for uncertainty-aware completion result.
  df$LDMC_bhpmf_mean <- mean_imp$LDMC
  df$LDMC_bhpmf_sd <- sd_imp$LDMC
  df$LDMC_bhpmf_ci95_low <- df$LDMC_bhpmf_mean - 1.96 * df$LDMC_bhpmf_sd
  df$LDMC_bhpmf_ci95_high <- df$LDMC_bhpmf_mean + 1.96 * df$LDMC_bhpmf_sd

  df$Leafph_bhpmf_mean <- mean_imp$Leafph
  df$Leafph_bhpmf_sd <- sd_imp$Leafph
  df$Leafph_bhpmf_ci95_low <- df$Leafph_bhpmf_mean - 1.96 * df$Leafph_bhpmf_sd
  df$Leafph_bhpmf_ci95_high <- df$Leafph_bhpmf_mean + 1.96 * df$Leafph_bhpmf_sd

  # Output completed table.
  write.csv(df, file.path(out_dir, "data_bhpmf_completed_with_uncertainty.csv"), row.names = FALSE)

  # Grouped + ungrouped plots.
  plot_outputs(df, out_dir)

  message("Done. Outputs written to: ", normalizePath(out_dir))
}

main()
