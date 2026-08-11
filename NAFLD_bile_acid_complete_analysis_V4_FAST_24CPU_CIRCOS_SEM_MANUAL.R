# ==============================================================================
# COMPLETE BILE-ACID ANALYSIS WORKFLOW
# Liver, plasma, ileum, cecum and feces: NC, MCD and 3_5 groups
# ==============================================================================
# Input expected:
#   NAFLD_BA_Allcontaint(4).csv
#
# Main outputs:
#   00_QC/                  data checks and sample summaries
#   01_Pathway_indices/     conjugation and transformation indices
#   02_Univariate/          group summaries, tests and effect sizes
#   03_Figures/             boxplots and heat maps
#   04_Multivariate/        PCA, PERMANOVA and dispersion tests
#   05_Correlation/         within-compartment residual correlations
#   06_Optional_matched/    mixed models / DIABLO / SEM, when Animal_ID exists
#
# Important:
# 1. The default analysis compares groups separately within each compartment.
# 2. Do not infer Animal_ID from sample names. Fill sample_metadata_map(3).csv first.
# 3. Ratios are pathway-balance proxies, not direct metabolic-flux estimates.
# 4. Features beginning with "log2_" are already log2 ratios and are never
#    log-transformed a second time.
# 5. Revised metadata contains 24 complete matched animals (8/group) across 5 parts.
# ==============================================================================

# ------------------------------- 0. SETTINGS -----------------------------------

set.seed(20260807)
options(stringsAsFactors = FALSE)

input_candidates <- c(
  "NAFLD_BA_Allcontaint(1).csv",
  "NAFLD_BA_Allcontaint.csv"
)
existing_input <- input_candidates[file.exists(input_candidates)]
input_file <- if (length(existing_input) > 0) existing_input[1] else input_candidates[1]
out_dir <- "NAFLD_BA_R_results_V4_FAST_24CPU_CIRCOSFIX"

# Number of permutations for PERMANOVA. Use 9999 for final analysis.
n_permutations <- 9999

# Whether to create one PDF page for every individual BA/pathway feature.
make_individual_boxplots <- TRUE

# --------------------- MANUAL SEM FEATURE SELECTION ----------------------------
# Choose ONE feature for each compartment in the proposed
# Liver -> Plasma -> Ileum -> Cecum -> Feces chain.
#
# IMPORTANT naming:
#   Raw "T.b.MCA" is renamed earlier in this workflow to "TbetaMCA".
#
# Any numeric column present in matched_data can be used, including:
#   - individual bile acids in ba_columns, e.g. "TbetaMCA", "TCA", "DCA"
#   - pathway summaries, e.g. "Conjugated_fraction", "Secondary_fraction"
#   - totals/families/indices, e.g. "Total_BA", "CA_to_DCA_index"
#
# Edit ONLY the five values below to test another biological chain.
SEM_FEATURE_BY_PART <- c(
  Liver  = "TbetaMCA",
  Plasma = "TbetaMCA",
  Ileum  = "TbetaMCA",
  Cecum  = "TbetaMCA",
  Feces  = "TbetaMCA"
)

# Example of a mixed-feature chain:
# SEM_FEATURE_BY_PART <- c(
#   Liver  = "TbetaMCA",
#   Plasma = "Conjugated_fraction",
#   Ileum  = "Conjugated_fraction",
#   Cecum  = "Secondary_fraction",
#   Feces  = "Total_BA"
# )
# individual bile acids
# "CA"
# "alphaMCA"
# "betaMCA"
# "omegaMCA"
# "DCA"
# "UDCA"
# "HDCA"
# "CDCA"
# "TUDCA"
# "TCDCA"
# "TDCA"
# "TCA"
# "TbetaMCA"
# "TalphaMCA"
# "GCDCA"
# "GDCA"
# "LCA"
# "GCA"

# pathway indices
# "Total_BA"
# "Conjugated_BA"
# "Unconjugated_BA"
# "Taurine_conjugated_BA"
# "Glycine_conjugated_BA"
# "Primary_BA"
# "Secondary_BA"
# 
# "CA_family"
# "CDCA_family"
# "DCA_family"
# "UDCA_family"
# "alphaMCA_family"
# "betaMCA_family"
# "omegaMCA_family"
# 
# "Conjugated_fraction"
# "Taurine_fraction_total"
# "Glycine_fraction_total"
# "Secondary_fraction"
# "Fraction_12aOH"
# 
# "CA_conjugation"
# "CDCA_conjugation"
# "DCA_conjugation"
# "UDCA_conjugation"
# "alphaMCA_conjugation"
# "betaMCA_conjugation"
# 
# "CA_to_DCA_index"
# "CDCA_to_UDCA_index"
# "betaMCA_to_omegaMCA_index"
# 
# "log2_taurine_to_glycine"
# "log2_12aOH_to_non12aOH"
# "log2_primary_to_secondary"

# ---------------------- mixOmics speed / parallel settings --------------------
# Maximum number of workers used by mixOmics. DIABLO can use the full 24
# workers because its joint keepX search is the expensive step. Tissue-specific
# sPLS-DA is tiny (about 24 samples x 18 features), so it is run serially. This
# avoids unnecessary SOCK connections and their cleanup/serialization warnings.
MIXOMICS_MAX_WORKERS <- 24L
SPLSDA_MAX_WORKERS <- 1L   # intentionally serial: avoids 4-worker socket warnings
DIABLO_MAX_WORKERS <- 24L # retain full parallelism for expensive multiblock tuning

# Coarse keepX grid. With ~18 bile acids per block, three biologically useful
# sparsity levels are sufficient and avoid the combinatorial explosion caused
# by seven values per block in DIABLO (7^5 = 16,807 combinations).
MIXOMICS_KEEPX_GRID <- c(3L, 8L, 15L)

# Hyperparameter tuning should be relatively light. Final-model performance is
# evaluated more thoroughly after the optimal keepX values have been selected.
SPLSDA_TUNE_NREPEAT <- 3L
SPLSDA_PERF_NREPEAT <- 20L
DIABLO_TUNE_NREPEAT <- 3L
DIABLO_PERF_NREPEAT <- 20L

# Optional metadata file. The script creates a template if it does not exist.
metadata_candidates <- c(
  "sample_metadata_map(1).csv",
  "sample_metadata_map.csv"
)
existing_metadata <- metadata_candidates[file.exists(metadata_candidates)]
metadata_file <- if (length(existing_metadata) > 0) {
  existing_metadata[1]
} else {
  metadata_candidates[length(metadata_candidates)]
}

# Core packages used in the main workflow.
core_packages <- c(
  "tidyverse", "broom", "vegan", "pheatmap", "patchwork", "ggrepel"
)

missing_core <- core_packages[
  !vapply(core_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_core) > 0) {
  stop(
    "Install the following packages before running this script:\n",
    paste0("install.packages(c(",
           paste(sprintf('"%s"', missing_core), collapse = ", "),
           "))")
  )
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(broom)
  library(vegan)
  library(pheatmap)
  library(patchwork)
  library(ggrepel)
})

# Optional packages for advanced matched-animal analyses:
# install.packages(c("lme4", "lmerTest", "emmeans", "glmmTMB",
#                    "piecewiseSEM", "zCompositions", "compositions"))
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install(c("mixOmics", "MOFA2"))

subdirs <- c(
  "00_QC", "01_Pathway_indices", "02_Univariate", "03_Figures",
  "04_Multivariate", "05_Correlation", "06_Optional_matched"
)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
walk(file.path(out_dir, subdirs), dir.create, showWarnings = FALSE, recursive = TRUE)

# ---------------------------- 1. HELPER FUNCTIONS ------------------------------

safe_divide <- function(numerator, denominator) {
  ifelse(is.finite(denominator) & denominator > 0,
         numerator / denominator,
         NA_real_)
}

sum_columns <- function(data, columns) {
  missing_columns <- setdiff(columns, names(data))
  if (length(missing_columns) > 0) {
    stop("Missing columns: ", paste(missing_columns, collapse = ", "))
  }
  rowSums(data[, columns, drop = FALSE], na.rm = TRUE)
}

minimum_positive <- function(x, fallback = 1e-12) {
  positive <- x[is.finite(x) & x > 0]
  if (length(positive) == 0) fallback else min(positive) / 2
}

format_p <- function(p) {
  case_when(
    is.na(p) ~ NA_character_,
    p < 0.001 ~ "<0.001",
    TRUE ~ sprintf("%.3f", p)
  )
}

significance_symbol <- function(q) {
  case_when(
    is.na(q) ~ "",
    q < 0.001 ~ "***",
    q < 0.01 ~ "**",
    q < 0.05 ~ "*",
    TRUE ~ ""
  )
}

# Build a small, adaptive keepX grid. The preferred grid is 3/8/15, but the
# values are clipped to the number of available variables. If clipping produces
# fewer than three unique values, construct three approximately spaced values
# whenever at least three predictors are available.
make_fast_keepx_grid <- function(p, preferred = MIXOMICS_KEEPX_GRID) {
  p <- as.integer(p)
  if (!is.finite(p) || p < 1L) return(integer(0))

  grid <- sort(unique(as.integer(pmin(preferred, p))))
  grid <- grid[grid >= 1L & grid <= p]

  if (length(grid) < 3L && p >= 3L) {
    grid <- sort(unique(as.integer(round(seq(1, p, length.out = 3L)))))
  }

  grid
}

# Reusable BiocParallel backends -------------------------------------------------
# Creating a fresh SnowParam for every tuning/performance call opens many SOCK
# connections. When those temporary objects are garbage-collected, R prints
# repeated "closing unused connection" warnings. Cache one backend per worker
# count, reuse it for all mixOmics calls, and explicitly stop it once at the end.
.mixomics_bp_cache <- new.env(parent = emptyenv())

get_reusable_bpparam <- function(workers) {
  workers <- max(1L, as.integer(workers))
  if (workers <= 1L || !requireNamespace("BiocParallel", quietly = TRUE)) {
    return(NULL)
  }

  key <- paste0("workers_", workers)
  if (exists(key, envir = .mixomics_bp_cache, inherits = FALSE)) {
    return(get(key, envir = .mixomics_bp_cache, inherits = FALSE))
  }

  bp <- BiocParallel::SnowParam(
    workers = workers,
    type = "SOCK",
    progressbar = FALSE,
    RNGseed = 20260807,
    stop.on.error = TRUE,
    exportglobals = FALSE
  )
  bp <- tryCatch(
    {
      BiocParallel::bpstart(bp)
      bp
    },
    error = function(e) {
      warning(
        paste0(
          "Could not start reusable SnowParam with ", workers,
          " workers; mixOmics will fall back to its non-BPPARAM mode. Reason: ",
          conditionMessage(e)
        ),
        call. = FALSE
      )
      NULL
    }
  )

  if (!is.null(bp)) {
    assign(key, bp, envir = .mixomics_bp_cache)
  }
  bp
}

stop_mixomics_backends <- function() {
  keys <- ls(envir = .mixomics_bp_cache, all.names = TRUE)
  if (length(keys) == 0L) return(invisible(NULL))

  for (key in keys) {
    bp <- get(key, envir = .mixomics_bp_cache, inherits = FALSE)
    try(BiocParallel::bpstop(bp), silent = TRUE)
    rm(list = key, envir = .mixomics_bp_cache)
  }
  invisible(NULL)
}

# Finalizer is a backup if execution is interrupted before the explicit cleanup
# at the bottom of this script.
.parallel_cleanup_guard <- new.env(parent = emptyenv())
reg.finalizer(
  .parallel_cleanup_guard,
  function(e) {
    try(stop_mixomics_backends(), silent = TRUE)
  },
  onexit = TRUE
)

# mixOmics versions differ: recent versions may expose BPPARAM, whereas older
# releases use cpus. Supply only an argument supported by the installed method.
get_mixomics_parallel_args <- function(fun, workers) {
  arg_names <- names(formals(fun))
  workers <- max(1L, as.integer(workers))

  if ("BPPARAM" %in% arg_names && requireNamespace("BiocParallel", quietly = TRUE)) {
    bp <- get_reusable_bpparam(workers)
    if (!is.null(bp)) return(list(BPPARAM = bp))
  }

  if ("cpus" %in% arg_names) {
    return(list(cpus = workers))
  }

  list()
}

# perf() is an S3 generic in mixOmics, so inspect its class-specific method first.
get_mixomics_perf_parallel_args <- function(object, workers) {
  for (cls in class(object)) {
    method <- tryCatch(
      utils::getS3method("perf", cls, optional = TRUE),
      error = function(e) NULL
    )
    if (!is.null(method)) {
      args <- get_mixomics_parallel_args(method, workers)
      if (length(args) > 0L) return(args)
    }
  }

  get_mixomics_parallel_args(mixOmics::perf, workers)
}

# Detect logical cores once. Reserve the full available pool (up to 24) for
# DIABLO, while tissue sPLS-DA is intentionally serial to avoid SOCK overhead.
detected_mixomics_cores <- suppressWarnings(parallel::detectCores(logical = TRUE))
if (!is.finite(detected_mixomics_cores) || detected_mixomics_cores < 1L) {
  detected_mixomics_cores <- 1L
}

MIXOMICS_WORKERS <- min(MIXOMICS_MAX_WORKERS, as.integer(detected_mixomics_cores))
SPLSDA_WORKERS <- min(SPLSDA_MAX_WORKERS, MIXOMICS_WORKERS)
DIABLO_WORKERS <- min(DIABLO_MAX_WORKERS, MIXOMICS_WORKERS)

message(
  "mixOmics speed settings: sPLS-DA workers=", SPLSDA_WORKERS,
  ", DIABLO workers=", DIABLO_WORKERS,
  ", keepX grid=", paste(MIXOMICS_KEEPX_GRID, collapse = "/"),
  ", sPLS-DA tune repeats=", SPLSDA_TUNE_NREPEAT,
  ", sPLS-DA perf repeats=", SPLSDA_PERF_NREPEAT,
  ", DIABLO tune repeats=", DIABLO_TUNE_NREPEAT,
  ", DIABLO perf repeats=", DIABLO_PERF_NREPEAT
)

write_csv(
  tibble(
    Setting = c(
      "Detected_logical_cores", "Maximum_workers",
      "sPLSDA_workers", "DIABLO_workers",
      "keepX_grid", "sPLSDA_tuning_repeats", "sPLSDA_performance_repeats",
      "DIABLO_tuning_repeats", "DIABLO_performance_repeats",
      "parallel_backend_policy"
    ),
    Value = c(
      as.character(detected_mixomics_cores),
      as.character(MIXOMICS_MAX_WORKERS),
      as.character(SPLSDA_WORKERS),
      as.character(DIABLO_WORKERS),
      paste(MIXOMICS_KEEPX_GRID, collapse = ";"),
      as.character(SPLSDA_TUNE_NREPEAT),
      as.character(SPLSDA_PERF_NREPEAT),
      as.character(DIABLO_TUNE_NREPEAT),
      as.character(DIABLO_PERF_NREPEAT),
      "Serial sPLS-DA; reusable SnowParam up to 24 workers for DIABLO"
    )
  ),
  file.path(out_dir, "00_QC", "mixOmics_speed_settings.csv")
)

# Cliff's delta is calculated as P(group2 > group1) - P(group2 < group1).
# Positive values mean that group2 tends to be higher than group1.
cliffs_delta <- function(group2, group1) {
  group2 <- group2[is.finite(group2)]
  group1 <- group1[is.finite(group1)]
  if (length(group2) == 0 || length(group1) == 0) return(NA_real_)
  differences <- outer(group2, group1, FUN = "-")
  mean(sign(differences))
}

cliffs_magnitude <- function(delta) {
  abs_delta <- abs(delta)
  case_when(
    is.na(abs_delta) ~ NA_character_,
    abs_delta < 0.147 ~ "negligible",
    abs_delta < 0.330 ~ "small",
    abs_delta < 0.474 ~ "medium",
    TRUE ~ "large"
  )
}

safe_kruskal <- function(data) {
  data <- data %>% filter(is.finite(Value), !is.na(Group))
  if (n_distinct(data$Group) < 2 || n_distinct(data$Value) < 2) {
    return(tibble(statistic = NA_real_, df = NA_real_, p_value = NA_real_))
  }
  fit <- tryCatch(
    kruskal.test(Value ~ Group, data = data),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    tibble(statistic = NA_real_, df = NA_real_, p_value = NA_real_)
  } else {
    tibble(
      statistic = unname(fit$statistic),
      df = unname(fit$parameter),
      p_value = fit$p.value
    )
  }
}

pairwise_test_one <- function(data, group1, group2, feature_name) {
  x1 <- data %>% filter(Group == group1) %>% pull(Value)
  x2 <- data %>% filter(Group == group2) %>% pull(Value)
  x1 <- x1[is.finite(x1)]
  x2 <- x2[is.finite(x2)]

  is_prelogged <- stringr::str_detect(as.character(feature_name), "^log2_")
  effect_definition <- if (is_prelogged) {
    "difference_in_existing_log2_ratio"
  } else {
    "log2_ratio_of_group_medians"
  }

  if (length(x1) < 2 || length(x2) < 2) {
    return(tibble(
      n_group1 = length(x1), n_group2 = length(x2),
      median_group1 = NA_real_, median_group2 = NA_real_,
      pseudocount = if (is_prelogged) 0 else NA_real_,
      log2FC = NA_real_,
      Effect_definition = effect_definition,
      W = NA_real_, p_value = NA_real_,
      cliffs_delta = NA_real_, cliffs_magnitude = NA_character_
    ))
  }

  median1 <- median(x1, na.rm = TRUE)
  median2 <- median(x2, na.rm = TRUE)

  if (is_prelogged) {
    # This variable is already on a log2-ratio scale and may be negative.
    # Its group effect is an arithmetic difference, not another logarithm.
    pseudo <- 0
    effect_log2 <- median2 - median1
  } else {
    pseudo <- minimum_positive(c(x1, x2))
    effect_log2 <- log2((median2 + pseudo) / (median1 + pseudo))
  }

  test <- tryCatch(
    wilcox.test(x2, x1, exact = FALSE, conf.int = FALSE),
    error = function(e) NULL
  )

  delta <- cliffs_delta(x2, x1)

  tibble(
    n_group1 = length(x1),
    n_group2 = length(x2),
    median_group1 = median1,
    median_group2 = median2,
    pseudocount = pseudo,
    log2FC = effect_log2,
    Effect_definition = effect_definition,
    W = if (is.null(test)) NA_real_ else unname(test$statistic),
    p_value = if (is.null(test)) NA_real_ else test$p.value,
    cliffs_delta = delta,
    cliffs_magnitude = cliffs_magnitude(delta)
  )
}

# Simple centered log-ratio transform after zero replacement.
# zCompositions::cmultRepl is used when available; otherwise a documented
# small-value replacement is applied.
make_clr <- function(x) {
  x <- as.matrix(x)
  storage.mode(x) <- "double"

  if (any(rowSums(x, na.rm = TRUE) <= 0)) {
    stop("At least one sample has a total BA concentration of zero.")
  }

  relative <- x / rowSums(x)

  if (requireNamespace("zCompositions", quietly = TRUE)) {
    replaced <- zCompositions::cmultRepl(
      relative,
      label = 0,
      method = "CZM",
      output = "prop"
    )
  } else {
    positive <- relative[relative > 0 & is.finite(relative)]
    replacement <- if (length(positive) == 0) 1e-12 else min(positive) / 2
    replaced <- relative
    replaced[!is.finite(replaced) | replaced <= 0] <- replacement
    replaced <- replaced / rowSums(replaced)
  }

  log_replaced <- log(replaced)
  clr <- log_replaced - rowMeans(log_replaced)
  colnames(clr) <- colnames(x)
  rownames(clr) <- rownames(x)
  clr
}

safe_scale_matrix <- function(x) {
  x <- as.matrix(x)
  variable <- apply(x, 2, sd, na.rm = TRUE) > 0
  x <- x[, variable, drop = FALSE]
  scale(x)
}

# ------------------------------ 2. READ DATA -----------------------------------

if (!file.exists(input_file)) {
  stop("Input file not found: ", normalizePath(input_file, mustWork = FALSE))
}

raw_data <- readr::read_csv(input_file, show_col_types = FALSE)

required_original <- c(
  "Sample.Name", "Group", "Part", "CA", "a.MCA", "b.MCA", "w.MCA",
  "DCA", "UDCA", "HDCA", "CDCA", "TUDCA", "TCDCA.1", "TDCA",
  "TCA", "T.b.MCA", "T.a.MCA", "GCDCA", "GDCA", "LCA", "GCA.2"
)

missing_original <- setdiff(required_original, names(raw_data))
if (length(missing_original) > 0) {
  stop("The input table is missing: ", paste(missing_original, collapse = ", "))
}

data <- raw_data %>%
  rename(
    Sample_ID = Sample.Name,
    alphaMCA = a.MCA,
    betaMCA = b.MCA,
    omegaMCA = w.MCA,
    TCDCA = TCDCA.1,
    TbetaMCA = T.b.MCA,
    TalphaMCA = T.a.MCA,
    GCA = GCA.2
  ) %>%
  mutate(
    Group = as.character(Group),
    Part = case_when(
      str_to_lower(str_trim(Part)) == "liver" ~ "Liver",
      str_to_lower(str_trim(Part)) == "plasma" ~ "Plasma",
      str_to_lower(str_trim(Part)) == "ileum" ~ "Ileum",
      str_to_lower(str_trim(Part)) == "cecum" ~ "Cecum",
      str_to_lower(str_trim(Part)) == "feces" ~ "Feces",
      TRUE ~ str_trim(as.character(Part))
    )
  )

ba_columns <- c(
  "CA", "alphaMCA", "betaMCA", "omegaMCA", "DCA", "UDCA", "HDCA",
  "CDCA", "TUDCA", "TCDCA", "TDCA", "TCA", "TbetaMCA",
  "TalphaMCA", "GCDCA", "GDCA", "LCA", "GCA"
)

data <- data %>% mutate(across(all_of(ba_columns), as.numeric))

if (any(data[ba_columns] < 0, na.rm = TRUE)) {
  stop("Negative bile-acid values were detected. Check the input data.")
}

if (anyDuplicated(data$Sample_ID) > 0) {
  data %>%
    filter(duplicated(Sample_ID) | duplicated(Sample_ID, fromLast = TRUE)) %>%
    write_csv(file.path(out_dir, "00_QC", "duplicated_sample_ID.csv"))
  stop("Duplicated Sample_ID values exist in the bile-acid table. See 00_QC/duplicated_sample_ID.csv")
}

# Create a metadata template only when no metadata file is available.
if (!file.exists(metadata_file)) {
  metadata_template <- data %>%
    select(Sample_ID, Group, Part) %>%
    mutate(
      Treatment_detail = Group,
      Animal_ID = NA_character_,
      Replicate_Set = NA_character_,
      Batch = NA_character_,
      Notes = NA_character_
    )
  write_csv(metadata_template, metadata_file, na = "")
  message("Created metadata template: ", metadata_file)
}

# Read metadata as character so mixed numeric/alphanumeric Animal_ID values
# (e.g. 1, 182, A80) are preserved exactly.
metadata <- readr::read_csv(
  metadata_file,
  show_col_types = FALSE,
  col_types = readr::cols(.default = readr::col_character())
)

if (!"Sample_ID" %in% names(metadata)) {
  stop("metadata_file must contain Sample_ID")
}
if (anyDuplicated(metadata$Sample_ID) > 0) {
  metadata %>%
    filter(duplicated(Sample_ID) | duplicated(Sample_ID, fromLast = TRUE)) %>%
    write_csv(file.path(out_dir, "00_QC", "duplicated_metadata_Sample_ID.csv"))
  stop("Duplicated Sample_ID values exist in metadata. See 00_QC/duplicated_metadata_Sample_ID.csv")
}

# Guarantee expected metadata columns exist.
for (nm in c("Group", "Part", "Treatment_detail", "Animal_ID",
             "Replicate_Set", "Batch", "Notes")) {
  if (!nm %in% names(metadata)) metadata[[nm]] <- NA_character_
}

metadata <- metadata %>%
  mutate(
    Group = if_else(is.na(Group) | Group == "", NA_character_, str_trim(Group)),
    Part = case_when(
      is.na(Part) | Part == "" ~ NA_character_,
      str_to_lower(str_trim(Part)) == "liver" ~ "Liver",
      str_to_lower(str_trim(Part)) == "plasma" ~ "Plasma",
      str_to_lower(str_trim(Part)) == "ileum" ~ "Ileum",
      str_to_lower(str_trim(Part)) == "cecum" ~ "Cecum",
      str_to_lower(str_trim(Part)) == "feces" ~ "Feces",
      TRUE ~ str_trim(Part)
    ),
    Treatment_detail = if_else(
      is.na(Treatment_detail) | Treatment_detail == "",
      Group,
      str_trim(Treatment_detail)
    ),
    Animal_ID = na_if(str_trim(Animal_ID), ""),
    Replicate_Set = na_if(str_trim(Replicate_Set), ""),
    Batch = na_if(str_trim(Batch), ""),
    Notes = na_if(str_trim(Notes), "")
  )

# Validate the two uploaded files before metadata is allowed to control grouping.
metadata_validation <- full_join(
  data %>%
    transmute(
      Sample_ID,
      BA_Group = as.character(Group),
      BA_Part = as.character(Part),
      In_BA = TRUE
    ),
  metadata %>%
    transmute(
      Sample_ID,
      Metadata_Group = Group,
      Metadata_Part = Part,
      In_Metadata = TRUE
    ),
  by = "Sample_ID"
) %>%
  mutate(
    In_BA = replace_na(In_BA, FALSE),
    In_Metadata = replace_na(In_Metadata, FALSE),
    Group_match = In_BA & In_Metadata & !is.na(Metadata_Group) &
      BA_Group == Metadata_Group,
    Part_match = In_BA & In_Metadata & !is.na(Metadata_Part) &
      BA_Part == Metadata_Part,
    Validation_status = case_when(
      !In_BA ~ "metadata_only_sample",
      !In_Metadata ~ "BA_only_sample",
      is.na(Metadata_Group) ~ "missing_metadata_group",
      is.na(Metadata_Part) ~ "missing_metadata_part",
      !Group_match ~ "group_mismatch",
      !Part_match ~ "part_mismatch",
      TRUE ~ "matched"
    )
  )

write_csv(
  metadata_validation,
  file.path(out_dir, "00_QC", "BA_metadata_sample_validation.csv")
)

bad_metadata_rows <- metadata_validation %>%
  filter(Validation_status != "matched")

if (nrow(bad_metadata_rows) > 0) {
  stop(
    "Bile-acid data and metadata are inconsistent. See ",
    file.path(out_dir, "00_QC", "BA_metadata_sample_validation.csv")
  )
}

# Current revised dataset design QC.
current_design_qc <- metadata %>%
  count(Group, Part, name = "n_samples") %>%
  arrange(Part, Group)

write_csv(
  current_design_qc,
  file.path(out_dir, "00_QC", "current_design_group_by_part_counts.csv")
)

current_animal_qc <- metadata %>%
  filter(!is.na(Animal_ID), Animal_ID != "") %>%
  group_by(Animal_ID) %>%
  summarise(
    n_groups = n_distinct(Group),
    Group = paste(sort(unique(as.character(Group))), collapse = ";"),
    n_parts = n_distinct(Part),
    Parts = paste(sort(unique(as.character(Part))), collapse = ";"),
    n_rows = n(),
    .groups = "drop"
  )

write_csv(
  current_animal_qc,
  file.path(out_dir, "00_QC", "current_metadata_AnimalID_design_QC.csv")
)

if (nrow(metadata) == 120L &&
    n_distinct(metadata$Animal_ID[!is.na(metadata$Animal_ID)]) == 24L &&
    all(current_design_qc$n_samples == 8L) &&
    all(current_animal_qc$n_groups == 1L) &&
    all(current_animal_qc$n_parts == 5L) &&
    all(current_animal_qc$n_rows == 5L)) {
  message(
    "Revised metadata design confirmed: 120 samples, 24 complete animals, ",
    "8 animals per group, five compartments per animal."
  )
} else {
  warning(
    "The loaded files do not exactly match the expected revised design. ",
    "Inspect 00_QC/current_design_group_by_part_counts.csv and ",
    "00_QC/current_metadata_AnimalID_design_QC.csv."
  )
}

# Metadata is now authoritative for Group, Part and Treatment_detail after the
# validation above. This avoids inference from sample-name strings.
data <- data %>%
  select(-any_of(c("Treatment_detail", "Animal_ID", "Replicate_Set", "Batch", "Notes"))) %>%
  left_join(
    metadata %>%
      select(Sample_ID, Group, Part, Treatment_detail,
             Animal_ID, Replicate_Set, Batch, Notes) %>%
      rename(
        Group_metadata = Group,
        Part_metadata = Part
      ),
    by = "Sample_ID"
  ) %>%
  mutate(
    Group = factor(Group_metadata, levels = c("NC", "MCD", "3_5")),
    Part = factor(Part_metadata, levels = c("Liver", "Plasma", "Ileum", "Cecum", "Feces")),
    Treatment_detail = coalesce(Treatment_detail, as.character(Group))
  ) %>%
  select(-Group_metadata, -Part_metadata)

if (any(is.na(data$Group))) {
  stop("Metadata contains Group values outside NC, MCD, 3_5.")
}
if (any(is.na(data$Part))) {
  stop("Metadata contains Part values outside Liver, Plasma, Ileum, Cecum, Feces.")
}

# Keep raw rows by default. Technical replicates are NOT automatically merged.
write_csv(data, file.path(out_dir, "00_QC", "cleaned_input_data.csv"))

# ------------------------------- 3. QC -----------------------------------------

sample_counts <- data %>% count(Part, Group, name = "n_samples")
write_csv(sample_counts, file.path(out_dir, "00_QC", "sample_counts.csv"))

missing_summary <- data %>%
  pivot_longer(all_of(ba_columns), names_to = "BA", values_to = "Value") %>%
  group_by(Part, Group, BA) %>%
  summarise(
    n = n(),
    n_missing = sum(is.na(Value)),
    percent_missing = 100 * mean(is.na(Value)),
    n_zero = sum(Value == 0, na.rm = TRUE),
    percent_zero = 100 * mean(Value == 0, na.rm = TRUE),
    .groups = "drop"
  )
write_csv(missing_summary, file.path(out_dir, "00_QC", "missing_and_zero_summary.csv"))

sample_qc <- data %>%
  mutate(
    Total_measured_BA = sum_columns(., ba_columns),
    Detected_BA_number = rowSums(across(all_of(ba_columns), ~ .x > 0), na.rm = TRUE)
  ) %>%
  select(Sample_ID, Group, Part, Treatment_detail, Animal_ID, Batch,
         Total_measured_BA, Detected_BA_number)
write_csv(sample_qc, file.path(out_dir, "00_QC", "sample_level_QC.csv"))

p_total <- ggplot(sample_qc, aes(Group, Total_measured_BA)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.12, height = 0, size = 2) +
  facet_wrap(~Part, scales = "free_y") +
  scale_y_continuous(trans = scales::pseudo_log_trans(base = 10)) +
  labs(x = NULL, y = "Total measured BA (pseudo-log scale)") +
  theme_bw(base_size = 12)

ggsave(
  file.path(out_dir, "00_QC", "total_BA_by_part_group.pdf"),
  p_total, width = 10, height = 6
)

p_detected <- ggplot(sample_qc, aes(Group, Detected_BA_number)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.12, height = 0, size = 2) +
  facet_wrap(~Part) +
  labs(x = NULL, y = "Number of detected BA species") +
  theme_bw(base_size = 12)

ggsave(
  file.path(out_dir, "00_QC", "detected_species_by_part_group.pdf"),
  p_detected, width = 10, height = 6
)

# ------------------------- 4. PATHWAY DEFINITIONS ------------------------------

unconjugated_ba <- c(
  "CA", "alphaMCA", "betaMCA", "omegaMCA", "DCA", "UDCA", "HDCA",
  "CDCA", "LCA"
)

taurine_conjugated_ba <- c(
  "TCA", "TCDCA", "TDCA", "TUDCA", "TalphaMCA", "TbetaMCA"
)

glycine_conjugated_ba <- c("GCA", "GCDCA", "GDCA")
conjugated_ba <- c(taurine_conjugated_ba, glycine_conjugated_ba)

# Parent-family definitions.
ca_family <- c("CA", "TCA", "GCA")
cdca_family <- c("CDCA", "TCDCA", "GCDCA")
dca_family <- c("DCA", "TDCA", "GDCA")
udca_family <- c("UDCA", "TUDCA")
alpha_mca_family <- c("alphaMCA", "TalphaMCA")
beta_mca_family <- c("betaMCA", "TbetaMCA")
omega_mca_family <- c("omegaMCA")

# Primary versus secondary classification used here.
# UDCA, DCA, LCA, HDCA and omegaMCA are treated as secondary/microbial products.
primary_ba <- c(
  ca_family, cdca_family, alpha_mca_family, beta_mca_family
)
secondary_ba <- c(
  dca_family, udca_family, "LCA", "HDCA", omega_mca_family
)

# Common 12alpha-hydroxylated pool.
ba_12aOH <- c(ca_family, dca_family)
non_12aOH <- setdiff(ba_columns, ba_12aOH)

# ----------------------- 5. CALCULATE PATHWAY INDICES --------------------------

analysis_data <- data

analysis_data$Total_BA <- sum_columns(analysis_data, ba_columns)
analysis_data$Conjugated_BA <- sum_columns(analysis_data, conjugated_ba)
analysis_data$Unconjugated_BA <- sum_columns(analysis_data, unconjugated_ba)
analysis_data$Taurine_conjugated_BA <- sum_columns(analysis_data, taurine_conjugated_ba)
analysis_data$Glycine_conjugated_BA <- sum_columns(analysis_data, glycine_conjugated_ba)
analysis_data$Primary_BA <- sum_columns(analysis_data, primary_ba)
analysis_data$Secondary_BA <- sum_columns(analysis_data, secondary_ba)
analysis_data$BA_12aOH <- sum_columns(analysis_data, ba_12aOH)
analysis_data$BA_non12aOH <- sum_columns(analysis_data, non_12aOH)

analysis_data$CA_family <- sum_columns(analysis_data, ca_family)
analysis_data$CDCA_family <- sum_columns(analysis_data, cdca_family)
analysis_data$DCA_family <- sum_columns(analysis_data, dca_family)
analysis_data$UDCA_family <- sum_columns(analysis_data, udca_family)
analysis_data$alphaMCA_family <- sum_columns(analysis_data, alpha_mca_family)
analysis_data$betaMCA_family <- sum_columns(analysis_data, beta_mca_family)
analysis_data$omegaMCA_family <- sum_columns(analysis_data, omega_mca_family)

analysis_data$Conjugated_fraction <- safe_divide(
  analysis_data$Conjugated_BA, analysis_data$Total_BA
)
analysis_data$Taurine_fraction_total <- safe_divide(
  analysis_data$Taurine_conjugated_BA, analysis_data$Total_BA
)
analysis_data$Glycine_fraction_total <- safe_divide(
  analysis_data$Glycine_conjugated_BA, analysis_data$Total_BA
)
analysis_data$Secondary_fraction <- safe_divide(
  analysis_data$Secondary_BA,
  analysis_data$Primary_BA + analysis_data$Secondary_BA
)
analysis_data$Fraction_12aOH <- safe_divide(
  analysis_data$BA_12aOH, analysis_data$Total_BA
)

analysis_data$CA_conjugation <- safe_divide(
  analysis_data$TCA + analysis_data$GCA,
  analysis_data$CA + analysis_data$TCA + analysis_data$GCA
)
analysis_data$CDCA_conjugation <- safe_divide(
  analysis_data$TCDCA + analysis_data$GCDCA,
  analysis_data$CDCA + analysis_data$TCDCA + analysis_data$GCDCA
)
analysis_data$DCA_conjugation <- safe_divide(
  analysis_data$TDCA + analysis_data$GDCA,
  analysis_data$DCA + analysis_data$TDCA + analysis_data$GDCA
)
analysis_data$UDCA_conjugation <- safe_divide(
  analysis_data$TUDCA,
  analysis_data$UDCA + analysis_data$TUDCA
)
analysis_data$alphaMCA_conjugation <- safe_divide(
  analysis_data$TalphaMCA,
  analysis_data$alphaMCA + analysis_data$TalphaMCA
)
analysis_data$betaMCA_conjugation <- safe_divide(
  analysis_data$TbetaMCA,
  analysis_data$betaMCA + analysis_data$TbetaMCA
)

# Product/(substrate + product) balance indices.
analysis_data$CA_to_DCA_index <- safe_divide(
  analysis_data$DCA_family,
  analysis_data$CA_family + analysis_data$DCA_family
)
analysis_data$CDCA_to_UDCA_index <- safe_divide(
  analysis_data$UDCA_family,
  analysis_data$CDCA_family + analysis_data$UDCA_family
)
analysis_data$betaMCA_to_omegaMCA_index <- safe_divide(
  analysis_data$omegaMCA,
  analysis_data$betaMCA + analysis_data$omegaMCA
)

# Log-ratio indices use a data-derived pseudocount only to avoid log(0).
global_pseudocount <- minimum_positive(unlist(analysis_data[ba_columns]))
analysis_data$log2_taurine_to_glycine <- log2(
  (analysis_data$Taurine_conjugated_BA + global_pseudocount) /
    (analysis_data$Glycine_conjugated_BA + global_pseudocount)
)
analysis_data$log2_12aOH_to_non12aOH <- log2(
  (analysis_data$BA_12aOH + global_pseudocount) /
    (analysis_data$BA_non12aOH + global_pseudocount)
)
analysis_data$log2_primary_to_secondary <- log2(
  (analysis_data$Primary_BA + global_pseudocount) /
    (analysis_data$Secondary_BA + global_pseudocount)
)

pathway_features <- c(
  "Total_BA", "Conjugated_BA", "Unconjugated_BA",
  "Taurine_conjugated_BA", "Glycine_conjugated_BA",
  "Primary_BA", "Secondary_BA",
  "CA_family", "CDCA_family", "DCA_family", "UDCA_family",
  "alphaMCA_family", "betaMCA_family", "omegaMCA_family",
  "Conjugated_fraction", "Taurine_fraction_total", "Glycine_fraction_total",
  "Secondary_fraction", "Fraction_12aOH",
  "CA_conjugation", "CDCA_conjugation", "DCA_conjugation",
  "UDCA_conjugation", "alphaMCA_conjugation", "betaMCA_conjugation",
  "CA_to_DCA_index", "CDCA_to_UDCA_index",
  "betaMCA_to_omegaMCA_index",
  "log2_taurine_to_glycine", "log2_12aOH_to_non12aOH",
  "log2_primary_to_secondary"
)

write_csv(
  analysis_data,
  file.path(out_dir, "01_Pathway_indices", "sample_level_BA_and_indices.csv")
)

# Long tables used throughout the script.
ba_long <- analysis_data %>%
  select(Sample_ID, Group, Part, Treatment_detail, Animal_ID, Batch,
         all_of(ba_columns)) %>%
  pivot_longer(all_of(ba_columns), names_to = "Feature", values_to = "Value") %>%
  mutate(Feature_type = "Individual_BA")

pathway_long <- analysis_data %>%
  select(Sample_ID, Group, Part, Treatment_detail, Animal_ID, Batch,
         all_of(pathway_features)) %>%
  pivot_longer(all_of(pathway_features), names_to = "Feature", values_to = "Value") %>%
  mutate(Feature_type = "Pathway_index")

all_long <- bind_rows(ba_long, pathway_long)

# ----------------------- 6. DESCRIPTIVE STATISTICS -----------------------------

group_summary <- all_long %>%
  group_by(Part, Group, Feature_type, Feature) %>%
  summarise(
    n_total = n(),
    n_nonmissing = sum(is.finite(Value)),
    mean = mean(Value, na.rm = TRUE),
    sd = sd(Value, na.rm = TRUE),
    median = median(Value, na.rm = TRUE),
    Q1 = quantile(Value, 0.25, na.rm = TRUE),
    Q3 = quantile(Value, 0.75, na.rm = TRUE),
    minimum = min(Value, na.rm = TRUE),
    maximum = max(Value, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(
  group_summary,
  file.path(out_dir, "02_Univariate", "group_descriptive_statistics.csv")
)

# --------------------- 7. KRUSKAL-WALLIS TESTS --------------------------------

kruskal_results <- all_long %>%
  group_by(Part, Feature_type, Feature) %>%
  nest() %>%
  mutate(test = map(data, safe_kruskal)) %>%
  select(-data) %>%
  unnest(test) %>%
  ungroup() %>%
  group_by(Part) %>%
  mutate(q_BH_within_part = p.adjust(p_value, method = "BH")) %>%
  ungroup() %>%
  mutate(
    q_BH_global = p.adjust(p_value, method = "BH"),
    p_display = format_p(p_value),
    q_symbol = significance_symbol(q_BH_within_part)
  ) %>%
  arrange(Part, q_BH_within_part, p_value)

write_csv(
  kruskal_results,
  file.path(out_dir, "02_Univariate", "kruskal_wallis_all_features.csv")
)

# ----------------- 8. PAIRWISE TESTS + EFFECT SIZES ----------------------------

contrasts <- tribble(
  ~Group1, ~Group2, ~Contrast,
  "NC",  "MCD", "MCD_vs_NC",
  "MCD", "3_5", "3_5_vs_MCD",
  "NC",  "3_5", "3_5_vs_NC"
)

pairwise_results <- all_long %>%
  group_by(Part, Feature_type, Feature) %>%
  nest() %>%
  crossing(contrasts) %>%
  mutate(
    result = pmap(
      list(data, Group1, Group2, Feature),
      ~pairwise_test_one(..1, ..2, ..3, ..4)
    )
  ) %>%
  select(-data) %>%
  unnest(result) %>%
  ungroup() %>%
  group_by(Part, Contrast) %>%
  mutate(q_BH_within_part_contrast = p.adjust(p_value, method = "BH")) %>%
  ungroup() %>%
  mutate(
    q_BH_global = p.adjust(p_value, method = "BH"),
    p_display = format_p(p_value),
    q_symbol = significance_symbol(q_BH_within_part_contrast)
  ) %>%
  arrange(Part, Contrast, q_BH_within_part_contrast, p_value)

write_csv(
  pairwise_results,
  file.path(out_dir, "02_Univariate", "pairwise_wilcoxon_effect_sizes.csv")
)

# ------------------ 9. PARAMETRIC LOG2 LINEAR MODELS ---------------------------
# These models estimate mean log2 differences. Use together with, not instead of,
# the rank-based tests above.

prepare_model_feature <- function(data, feature_name) {
  data <- data %>%
    filter(is.finite(Value), !is.na(Group)) %>%
    mutate(
      Group = factor(as.character(Group), levels = c("NC", "MCD", "3_5"))
    )

  is_prelogged <- stringr::str_detect(as.character(feature_name), "^log2_")

  if (is_prelogged) {
    # These pathway variables are already log2 ratios and may legitimately be
    # negative. Applying log2() again would generate NaNs.
    data <- data %>%
      mutate(
        Pseudocount = 0,
        Log2Value = Value,
        Transformation = "already_log2_ratio"
      )
  } else {
    pseudo <- minimum_positive(data$Value)
    data <- data %>%
      mutate(
        Pseudocount = pseudo,
        Log2Value = log2(Value + pseudo),
        Transformation = "log2_value_plus_pseudocount"
      )
  }

  data %>%
    filter(is.finite(Log2Value)) %>%
    droplevels()
}

empty_lm_contrasts <- function(data, status) {
  pseudo <- if (nrow(data) > 0) first(data$Pseudocount) else NA_real_
  transformation <- if (nrow(data) > 0) {
    first(data$Transformation)
  } else {
    NA_character_
  }

  tibble(
    Contrast = c("MCD_vs_NC", "3_5_vs_MCD", "3_5_vs_NC"),
    estimate_log2FC = NA_real_,
    SE = NA_real_,
    df = NA_real_,
    statistic = NA_real_,
    p_value = NA_real_,
    Pseudocount = pseudo,
    Transformation = transformation,
    Model_status = status
  )
}

fit_log2_lm <- function(data, feature_name) {
  data <- prepare_model_feature(data, feature_name)

  required_groups <- c("NC", "MCD", "3_5")
  observed_groups <- unique(as.character(data$Group))
  group_counts <- table(factor(data$Group, levels = required_groups))

  if (!all(required_groups %in% observed_groups) ||
      any(group_counts < 2) ||
      n_distinct(data$Log2Value) < 2) {
    return(empty_lm_contrasts(
      data,
      "skipped_insufficient_groups_or_variation"
    ))
  }

  # Reapply the intended reference and coefficient order after filtering.
  data <- data %>%
    mutate(Group = factor(as.character(Group), levels = required_groups))

  fit <- tryCatch(
    lm(Log2Value ~ Group, data = data, na.action = na.fail),
    error = function(e) NULL
  )

  if (is.null(fit)) {
    return(empty_lm_contrasts(data, "skipped_lm_error"))
  }

  manual_contrasts <- function(fit_object) {
    model_coef <- coef(fit_object)
    vc <- vcov(fit_object)
    needed <- c("GroupMCD", "Group3_5")

    if (!all(needed %in% names(model_coef)) ||
        !all(needed %in% rownames(vc)) ||
        !all(needed %in% colnames(vc))) {
      return(NULL)
    }

    b_mcd <- unname(model_coef["GroupMCD"])
    b_35 <- unname(model_coef["Group3_5"])

    get_manual <- function(name, estimate, variance) {
      if (!is.finite(variance) || variance < 0) {
        return(tibble(
          Contrast = name,
          estimate_log2FC = estimate,
          SE = NA_real_,
          df = df.residual(fit_object),
          statistic = NA_real_,
          p_value = NA_real_
        ))
      }

      se <- sqrt(variance)

      if (!is.finite(se) || se == 0) {
        t_value <- NA_real_
        p <- NA_real_
      } else {
        t_value <- estimate / se
        p <- 2 * pt(
          abs(t_value),
          df = df.residual(fit_object),
          lower.tail = FALSE
        )
      }

      tibble(
        Contrast = name,
        estimate_log2FC = estimate,
        SE = se,
        df = df.residual(fit_object),
        statistic = t_value,
        p_value = p
      )
    }

    bind_rows(
      get_manual(
        "MCD_vs_NC",
        b_mcd,
        vc["GroupMCD", "GroupMCD"]
      ),
      get_manual(
        "3_5_vs_MCD",
        b_35 - b_mcd,
        vc["Group3_5", "Group3_5"] +
          vc["GroupMCD", "GroupMCD"] -
          2 * vc["Group3_5", "GroupMCD"]
      ),
      get_manual(
        "3_5_vs_NC",
        b_35,
        vc["Group3_5", "Group3_5"]
      )
    )
  }

  out <- NULL

  if (requireNamespace("emmeans", quietly = TRUE)) {
    out <- tryCatch({
      emm <- emmeans::emmeans(fit, ~Group)
      con <- emmeans::contrast(
        emm,
        method = list(
          MCD_vs_NC = c(-1, 1, 0),
          `3_5_vs_MCD` = c(0, -1, 1),
          `3_5_vs_NC` = c(-1, 0, 1)
        )
      )

      as.data.frame(summary(con, infer = c(TRUE, TRUE))) %>%
        transmute(
          Contrast = as.character(contrast),
          estimate_log2FC = estimate,
          SE = SE,
          df = df,
          statistic = t.ratio,
          p_value = p.value
        )
    }, error = function(e) NULL)
  }

  if (is.null(out)) {
    out <- manual_contrasts(fit)
  }

  if (is.null(out)) {
    return(empty_lm_contrasts(data, "skipped_contrast_error"))
  }

  out %>%
    mutate(
      Pseudocount = first(data$Pseudocount),
      Transformation = first(data$Transformation),
      Model_status = "fitted"
    )
}

lm_results <- all_long %>%
  group_by(Part, Feature_type, Feature) %>%
  nest() %>%
  mutate(result = map2(data, Feature, fit_log2_lm)) %>%
  select(-data) %>%
  unnest(result) %>%
  ungroup() %>%
  group_by(Part, Contrast) %>%
  mutate(q_BH_within_part_contrast = p.adjust(p_value, method = "BH")) %>%
  ungroup() %>%
  mutate(q_BH_global = p.adjust(p_value, method = "BH"))

write_csv(
  lm_results,
  file.path(out_dir, "02_Univariate", "log2_linear_model_contrasts.csv")
)

# ------------------ 10. DISEASE/REGULATION PATTERN ANALYSIS --------------------
# Restoration_fraction on log2 medians:
#   1 = complete return from MCD to NC
#   0 = no movement away from MCD
#  >1 = overshoot beyond NC
#  <0 = further change in the disease direction

summarise_feature_profile <- function(data, feature_name) {
  transformed <- prepare_model_feature(data, feature_name)

  transformed %>%
    group_by(Group) %>%
    summarise(
      Median_raw = median(Value, na.rm = TRUE),
      Median_log2 = median(Log2Value, na.rm = TRUE),
      Pseudocount = first(Pseudocount),
      Transformation = first(Transformation),
      .groups = "drop"
    )
}

median_profiles <- all_long %>%
  group_by(Part, Feature_type, Feature) %>%
  nest() %>%
  mutate(
    profile = map2(data, Feature, summarise_feature_profile)
  ) %>%
  select(-data) %>%
  unnest(profile) %>%
  pivot_wider(
    names_from = Group,
    values_from = c(Median_raw, Median_log2),
    names_sep = "__"
  ) %>%
  mutate(
    Disease_effect = Median_log2__MCD - Median_log2__NC,
    Regulation_effect = Median_log2__3_5 - Median_log2__MCD,
    Residual_difference_from_NC = Median_log2__3_5 - Median_log2__NC,
    Restoration_fraction = ifelse(
      abs(Disease_effect) > 1e-12,
      (Median_log2__MCD - Median_log2__3_5) /
        (Median_log2__MCD - Median_log2__NC),
      NA_real_
    ),
    MCD_distance_from_NC = abs(Median_log2__MCD - Median_log2__NC),
    Group35_distance_from_NC = abs(Median_log2__3_5 - Median_log2__NC),
    Pattern = case_when(
      is.na(Disease_effect) ~ "Insufficient data",
      abs(Disease_effect) < 0.25 ~ "Small disease change",
      Group35_distance_from_NC < MCD_distance_from_NC &
        sign(Regulation_effect) == -sign(Disease_effect) &
        Restoration_fraction >= 0 & Restoration_fraction <= 1.25 ~ "Restored toward NC",
      sign(Regulation_effect) == -sign(Disease_effect) &
        Restoration_fraction > 1.25 ~ "Overshoot beyond NC",
      sign(Regulation_effect) == sign(Disease_effect) ~ "Further shifted from NC",
      abs(Regulation_effect) < 0.25 ~ "Little regulation effect",
      TRUE ~ "Complex/intermediate"
    )
  )

write_csv(
  median_profiles,
  file.path(out_dir, "02_Univariate", "disease_regulation_restoration_patterns.csv")
)

# Robust sample-level distance from the NC profile within each compartment.
# Each BA is log2 transformed, then standardized using the NC median and NC MAD.
calculate_nc_distance <- function(part_data) {
  x <- as.matrix(part_data[, ba_columns, drop = FALSE])
  pseudo <- apply(x, 2, minimum_positive)
  xlog <- sweep(x, 2, pseudo, "+")
  xlog <- log2(xlog)

  nc_rows <- part_data$Group == "NC"
  nc_median <- apply(xlog[nc_rows, , drop = FALSE], 2, median, na.rm = TRUE)
  nc_mad <- apply(xlog[nc_rows, , drop = FALSE], 2, mad, constant = 1, na.rm = TRUE)

  # Fall back to NC SD, then 1, when MAD is zero.
  nc_sd <- apply(xlog[nc_rows, , drop = FALSE], 2, sd, na.rm = TRUE)
  nc_mad[!is.finite(nc_mad) | nc_mad == 0] <- nc_sd[!is.finite(nc_mad) | nc_mad == 0]
  nc_mad[!is.finite(nc_mad) | nc_mad == 0] <- 1

  z <- sweep(xlog, 2, nc_median, "-")
  z <- sweep(z, 2, nc_mad, "/")
  distance <- sqrt(rowMeans(z^2, na.rm = TRUE))

  part_data %>%
    transmute(
      Sample_ID, Group, Part, Treatment_detail,
      Distance_from_NC_profile = distance
    )
}

nc_distance <- analysis_data %>%
  group_split(Part, .keep = TRUE) %>%
  map_dfr(calculate_nc_distance)

write_csv(
  nc_distance,
  file.path(out_dir, "02_Univariate", "sample_distance_from_NC_profile.csv")
)

nc_distance_test <- nc_distance %>%
  filter(Group %in% c("MCD", "3_5")) %>%
  group_by(Part) %>%
  summarise(
    median_MCD = median(Distance_from_NC_profile[Group == "MCD"], na.rm = TRUE),
    median_3_5 = median(Distance_from_NC_profile[Group == "3_5"], na.rm = TRUE),
    p_value = tryCatch(
      wilcox.test(
        Distance_from_NC_profile[Group == "3_5"],
        Distance_from_NC_profile[Group == "MCD"],
        exact = FALSE
      )$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) %>%
  mutate(q_BH = p.adjust(p_value, "BH"))

write_csv(
  nc_distance_test,
  file.path(out_dir, "02_Univariate", "MCD_vs_3_5_distance_to_NC_test.csv")
)

p_nc_distance <- ggplot(nc_distance, aes(Group, Distance_from_NC_profile)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.12, size = 2) +
  facet_wrap(~Part, scales = "free_y") +
  labs(x = NULL, y = "Robust multivariate distance from NC profile") +
  theme_bw(base_size = 12)

ggsave(
  file.path(out_dir, "03_Figures", "distance_from_NC_profile.pdf"),
  p_nc_distance, width = 10, height = 6
)

# ----------------------------- 11. BOXPLOTS ------------------------------------

# All individual BAs in a compact faceted PDF.
ba_plot_data <- ba_long %>%
  group_by(Part, Feature) %>%
  mutate(Pseudocount = minimum_positive(Value), Log10_value = log10(Value + Pseudocount)) %>%
  ungroup()

p_all_ba <- ggplot(ba_plot_data, aes(Group, Log10_value)) +
  geom_boxplot(outlier.shape = NA, width = 0.65) +
  geom_jitter(width = 0.13, size = 0.8, alpha = 0.8) +
  facet_grid(Part ~ Feature, scales = "free_y") +
  labs(x = NULL, y = "log10(concentration + feature-specific pseudocount)") +
  theme_bw(base_size = 8) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text.x = element_text(size = 7),
    strip.text.y = element_text(size = 8)
  )

ggsave(
  file.path(out_dir, "03_Figures", "all_individual_BA_boxplots.pdf"),
  p_all_ba, width = 24, height = 13, limitsize = FALSE
)

fraction_features <- pathway_features[str_detect(pathway_features, "fraction|conjugation|index")]
p_fraction <- pathway_long %>%
  filter(Feature %in% fraction_features) %>%
  ggplot(aes(Group, Value)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.13, size = 1) +
  facet_grid(Part ~ Feature, scales = "free_y") +
  labs(x = NULL, y = "Pathway balance") +
  theme_bw(base_size = 8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  file.path(out_dir, "03_Figures", "all_pathway_fraction_boxplots.pdf"),
  p_fraction, width = 22, height = 13, limitsize = FALSE
)

if (make_individual_boxplots) {
  individual_pdf <- file.path(out_dir, "03_Figures", "one_page_per_feature.pdf")
  pdf(individual_pdf, width = 10, height = 6)
  for (feature_name in unique(all_long$Feature)) {
    plot_data <- all_long %>% filter(Feature == feature_name)
    p <- ggplot(plot_data, aes(Group, Value)) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(width = 0.13, size = 2) +
      facet_wrap(~Part, scales = "free_y") +
      scale_y_continuous(trans = scales::pseudo_log_trans(base = 10)) +
      labs(title = feature_name, x = NULL, y = "Value (pseudo-log scale)") +
      theme_bw(base_size = 12)
    print(p)
  }
  dev.off()
}

# ---------------------------- 12. HEAT MAPS ------------------------------------

make_effect_heatmap <- function(feature_type, output_name) {
  # pairwise_results may retain dplyr grouping attributes. If Part is kept as
  # an implicit ID column during pivot_wider(), the same BA name occurs once
  # per tissue and cannot be used directly as row names. Explicitly ungroup,
  # construct a unique Feature x (Part__Contrast) key, and aggregate any
  # accidental duplicate key before widening.
  heat <- pairwise_results %>%
    ungroup() %>%
    filter(Feature_type == feature_type) %>%
    transmute(
      Feature = as.character(Feature),
      Column = paste(as.character(Part), as.character(Contrast), sep = "__"),
      log2FC = as.numeric(log2FC),
      q_value = as.numeric(q_BH_within_part_contrast)
    ) %>%
    group_by(Feature, Column) %>%
    summarise(
      log2FC = if (all(!is.finite(log2FC))) {
        NA_real_
      } else {
        median(log2FC[is.finite(log2FC)], na.rm = TRUE)
      },
      q_value = if (all(is.na(q_value))) {
        NA_real_
      } else {
        min(q_value, na.rm = TRUE)
      },
      .groups = "drop"
    ) %>%
    mutate(q_symbol = significance_symbol(q_value))

  if (nrow(heat) == 0L) {
    warning("No data available for heat map: ", feature_type)
    return(invisible(NULL))
  }

  effect_df <- heat %>%
    select(Feature, Column, log2FC) %>%
    pivot_wider(
      id_cols = Feature,
      names_from = Column,
      values_from = log2FC
    ) %>%
    arrange(Feature)

  symbol_df <- heat %>%
    select(Feature, Column, q_symbol) %>%
    pivot_wider(
      id_cols = Feature,
      names_from = Column,
      values_from = q_symbol,
      values_fill = ""
    ) %>%
    arrange(Feature)

  # Convert only after Feature is guaranteed to be unique.
  stopifnot(!anyDuplicated(effect_df$Feature))
  stopifnot(!anyDuplicated(symbol_df$Feature))

  effect_matrix <- effect_df %>%
    column_to_rownames("Feature") %>%
    as.matrix()

  symbol_matrix <- symbol_df %>%
    column_to_rownames("Feature") %>%
    as.matrix()

  symbol_matrix <- symbol_matrix[
    rownames(effect_matrix),
    colnames(effect_matrix),
    drop = FALSE
  ]

  pdf(
    file.path(out_dir, "03_Figures", output_name),
    width = max(11, ncol(effect_matrix) * 0.7),
    height = max(7, nrow(effect_matrix) * 0.35)
  )
  pheatmap::pheatmap(
    effect_matrix,
    cluster_rows = nrow(effect_matrix) > 1,
    cluster_cols = FALSE,
    scale = "none",
    display_numbers = symbol_matrix,
    number_color = "black",
    fontsize_number = 10,
    main = paste0(
      feature_type,
      ": median log2 fold changes\n",
      "* BH q<0.05; ** q<0.01; *** q<0.001"
    )
  )
  dev.off()

  invisible(list(effect_matrix = effect_matrix, symbol_matrix = symbol_matrix))
}

make_effect_heatmap("Individual_BA", "individual_BA_log2FC_heatmap.pdf")
make_effect_heatmap("Pathway_index", "pathway_indices_log2FC_heatmap.pdf")

# Restoration heat map.
restoration_df <- median_profiles %>%
  ungroup() %>%
  filter(Feature_type == "Individual_BA") %>%
  transmute(
    Feature = as.character(Feature),
    Column = as.character(Part),
    Restoration_fraction = as.numeric(Restoration_fraction)
  ) %>%
  group_by(Feature, Column) %>%
  summarise(
    Restoration_fraction = if (all(!is.finite(Restoration_fraction))) {
      NA_real_
    } else {
      median(Restoration_fraction[is.finite(Restoration_fraction)], na.rm = TRUE)
    },
    .groups = "drop"
  ) %>%
  pivot_wider(
    id_cols = Feature,
    names_from = Column,
    values_from = Restoration_fraction
  ) %>%
  arrange(Feature)

stopifnot(!anyDuplicated(restoration_df$Feature))

restoration_matrix <- restoration_df %>%
  column_to_rownames("Feature") %>%
  as.matrix()

pdf(
  file.path(out_dir, "03_Figures", "BA_restoration_fraction_heatmap.pdf"),
  width = 8, height = 9
)
pheatmap::pheatmap(
  restoration_matrix,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  main = "3_5 restoration fraction\n1 = return to NC; 0 = no reversal; <0 = further shifted"
)
dev.off()

# -------------------- 13. PCA + PERMANOVA: ABSOLUTE DATA -----------------------

run_absolute_multivariate <- function(part_name) {
  part_data <- analysis_data %>% filter(Part == part_name)
  x <- as.matrix(part_data[, ba_columns, drop = FALSE])
  rownames(x) <- part_data$Sample_ID

  pseudocounts <- apply(x, 2, minimum_positive)
  log_x <- sweep(x, 2, pseudocounts, "+")
  log_x <- log10(log_x)
  scaled_x <- safe_scale_matrix(log_x)

  pca <- prcomp(scaled_x, center = FALSE, scale. = FALSE)
  explained <- 100 * pca$sdev^2 / sum(pca$sdev^2)

  scores <- as_tibble(pca$x[, 1:min(3, ncol(pca$x)), drop = FALSE],
                      rownames = "Sample_ID") %>%
    left_join(part_data %>% select(Sample_ID, Group, Treatment_detail), by = "Sample_ID")

  loadings <- as_tibble(pca$rotation, rownames = "BA")
  write_csv(scores,
            file.path(out_dir, "04_Multivariate",
                      paste0(part_name, "_absolute_PCA_scores.csv")))
  write_csv(loadings,
            file.path(out_dir, "04_Multivariate",
                      paste0(part_name, "_absolute_PCA_loadings.csv")))

  p <- ggplot(scores, aes(PC1, PC2, shape = Group, label = Sample_ID)) +
    geom_point(size = 3) +
    ggrepel::geom_text_repel(size = 2.2, max.overlaps = 20) +
    labs(
      title = paste(part_name, "PCA: log-scaled absolute BA"),
      x = sprintf("PC1 (%.1f%%)", explained[1]),
      y = sprintf("PC2 (%.1f%%)", explained[2])
    ) +
    theme_bw(base_size = 12)

  ggsave(
    file.path(out_dir, "04_Multivariate",
              paste0(part_name, "_absolute_PCA.pdf")),
    p, width = 8, height = 6
  )

  distance <- dist(scaled_x, method = "euclidean")
  meta <- part_data %>% select(Group, Treatment_detail)

  permanova <- vegan::adonis2(
    distance ~ Group,
    data = meta,
    permutations = n_permutations,
    by = "margin"
  )
  permanova_out <- as.data.frame(permanova) %>%
    rownames_to_column("Term") %>%
    mutate(Part = part_name, Analysis = "Absolute_log_scaled")

  dispersion <- vegan::betadisper(distance, part_data$Group)
  dispersion_test <- vegan::permutest(dispersion, permutations = n_permutations)
  dispersion_out <- as.data.frame(dispersion_test$tab) %>%
    rownames_to_column("Term") %>%
    mutate(Part = part_name, Analysis = "Absolute_log_scaled")

  list(permanova = permanova_out, dispersion = dispersion_out)
}

absolute_results <- map(
  levels(droplevels(analysis_data$Part)),
  run_absolute_multivariate
)

absolute_permanova <- map_dfr(absolute_results, "permanova")
absolute_dispersion <- map_dfr(absolute_results, "dispersion")

write_csv(
  absolute_permanova,
  file.path(out_dir, "04_Multivariate", "absolute_PERMANOVA.csv")
)
write_csv(
  absolute_dispersion,
  file.path(out_dir, "04_Multivariate", "absolute_beta_dispersion_tests.csv")
)

# ----------------------- 14. COMPOSITIONAL CLR ANALYSIS -------------------------

run_clr_multivariate <- function(part_name) {
  part_data <- analysis_data %>% filter(Part == part_name)
  x <- as.matrix(part_data[, ba_columns, drop = FALSE])
  rownames(x) <- part_data$Sample_ID

  clr_x <- make_clr(x)
  variable <- apply(clr_x, 2, sd, na.rm = TRUE) > 0
  clr_x <- clr_x[, variable, drop = FALSE]

  pca <- prcomp(clr_x, center = TRUE, scale. = FALSE)
  explained <- 100 * pca$sdev^2 / sum(pca$sdev^2)

  scores <- as_tibble(pca$x[, 1:min(3, ncol(pca$x)), drop = FALSE],
                      rownames = "Sample_ID") %>%
    left_join(part_data %>% select(Sample_ID, Group, Treatment_detail), by = "Sample_ID")
  loadings <- as_tibble(pca$rotation, rownames = "BA")

  write_csv(scores,
            file.path(out_dir, "04_Multivariate",
                      paste0(part_name, "_CLR_PCA_scores.csv")))
  write_csv(loadings,
            file.path(out_dir, "04_Multivariate",
                      paste0(part_name, "_CLR_PCA_loadings.csv")))

  p <- ggplot(scores, aes(PC1, PC2, shape = Group, label = Sample_ID)) +
    geom_point(size = 3) +
    ggrepel::geom_text_repel(size = 2.2, max.overlaps = 20) +
    labs(
      title = paste(part_name, "PCA: CLR bile-acid composition"),
      x = sprintf("PC1 (%.1f%%)", explained[1]),
      y = sprintf("PC2 (%.1f%%)", explained[2])
    ) +
    theme_bw(base_size = 12)

  ggsave(
    file.path(out_dir, "04_Multivariate",
              paste0(part_name, "_CLR_PCA.pdf")),
    p, width = 8, height = 6
  )

  distance <- dist(clr_x, method = "euclidean")  # Aitchison distance
  meta <- part_data %>% select(Group, Treatment_detail)

  permanova <- vegan::adonis2(
    distance ~ Group,
    data = meta,
    permutations = n_permutations,
    by = "margin"
  )
  permanova_out <- as.data.frame(permanova) %>%
    rownames_to_column("Term") %>%
    mutate(Part = part_name, Analysis = "CLR_composition")

  dispersion <- vegan::betadisper(distance, part_data$Group)
  dispersion_test <- vegan::permutest(dispersion, permutations = n_permutations)
  dispersion_out <- as.data.frame(dispersion_test$tab) %>%
    rownames_to_column("Term") %>%
    mutate(Part = part_name, Analysis = "CLR_composition")

  list(permanova = permanova_out, dispersion = dispersion_out)
}

clr_results <- map(
  levels(droplevels(analysis_data$Part)),
  run_clr_multivariate
)

clr_permanova <- map_dfr(clr_results, "permanova")
clr_dispersion <- map_dfr(clr_results, "dispersion")

write_csv(
  clr_permanova,
  file.path(out_dir, "04_Multivariate", "CLR_PERMANOVA.csv")
)
write_csv(
  clr_dispersion,
  file.path(out_dir, "04_Multivariate", "CLR_beta_dispersion_tests.csv")
)

# ------------------ 15. GROUP-CENTROID TRAJECTORY IN PCA -----------------------
# This illustrates the NC -> MCD -> 3_5 direction for each compartment.

make_centroid_trajectory <- function(part_name) {
  part_data <- analysis_data %>% filter(Part == part_name)
  x <- as.matrix(part_data[, ba_columns, drop = FALSE])
  pseudo <- apply(x, 2, minimum_positive)
  xlog <- log10(sweep(x, 2, pseudo, "+"))
  xscaled <- safe_scale_matrix(xlog)
  pca <- prcomp(xscaled, center = FALSE, scale. = FALSE)
  explained <- 100 * pca$sdev^2 / sum(pca$sdev^2)

  scores <- as_tibble(pca$x[, 1:2, drop = FALSE]) %>%
    bind_cols(part_data %>% select(Group))

  centroids <- scores %>%
    group_by(Group) %>%
    summarise(PC1 = mean(PC1), PC2 = mean(PC2), .groups = "drop") %>%
    mutate(Group_order = match(Group, c("NC", "MCD", "3_5"))) %>%
    arrange(Group_order)

  p <- ggplot(scores, aes(PC1, PC2, shape = Group)) +
    geom_point(alpha = 0.55, size = 2.5) +
    geom_path(
      data = centroids,
      aes(PC1, PC2, group = 1),
      inherit.aes = FALSE,
      linewidth = 0.8,
      arrow = arrow(length = grid::unit(0.18, "cm"), type = "closed")
    ) +
    geom_point(
      data = centroids,
      aes(PC1, PC2, shape = Group),
      inherit.aes = FALSE,
      size = 4
    ) +
    geom_text_repel(
      data = centroids,
      aes(PC1, PC2, label = Group),
      inherit.aes = FALSE
    ) +
    labs(
      title = paste(part_name, "group-centroid trajectory"),
      subtitle = "NC -> MCD -> 3_5",
      x = sprintf("PC1 (%.1f%%)", explained[1]),
      y = sprintf("PC2 (%.1f%%)", explained[2])
    ) +
    theme_bw(base_size = 12)

  ggsave(
    file.path(out_dir, "04_Multivariate",
              paste0(part_name, "_PCA_centroid_trajectory.pdf")),
    p, width = 7, height = 6
  )
}

walk(levels(droplevels(analysis_data$Part)), make_centroid_trajectory)

# ---------------- 16. RESIDUAL SPEARMAN CORRELATION NETWORKS ------------------
# Group effects are removed feature-by-feature before correlation. This reduces
# correlations that are caused only by NC/MCD/3_5 mean differences.

pairwise_cor_test <- function(matrix_data) {
  features <- colnames(matrix_data)
  pairs <- combn(features, 2, simplify = FALSE)

  map_dfr(pairs, function(pair) {
    x <- matrix_data[, pair[1]]
    y <- matrix_data[, pair[2]]
    keep <- is.finite(x) & is.finite(y)
    if (sum(keep) < 5 || sd(x[keep]) == 0 || sd(y[keep]) == 0) {
      return(tibble(
        Feature1 = pair[1], Feature2 = pair[2],
        n = sum(keep), rho = NA_real_, p_value = NA_real_
      ))
    }
    test <- suppressWarnings(cor.test(x[keep], y[keep], method = "spearman", exact = FALSE))
    tibble(
      Feature1 = pair[1], Feature2 = pair[2],
      n = sum(keep), rho = unname(test$estimate), p_value = test$p.value
    )
  }) %>%
    mutate(q_BH = p.adjust(p_value, "BH"))
}

run_residual_correlation <- function(part_name) {
  part_data <- analysis_data %>% filter(Part == part_name)
  x <- part_data %>% select(all_of(ba_columns))

  log_x <- map_dfc(x, function(column) {
    pseudo <- minimum_positive(column)
    log2(column + pseudo)
  })
  names(log_x) <- ba_columns

  residual_matrix <- map_dfc(ba_columns, function(feature) {
    fit <- lm(log_x[[feature]] ~ part_data$Group)
    tibble(value = residuals(fit))
  })
  names(residual_matrix) <- ba_columns

  correlation_table <- pairwise_cor_test(as.matrix(residual_matrix))
  write_csv(
    correlation_table,
    file.path(out_dir, "05_Correlation",
              paste0(part_name, "_group_adjusted_spearman_correlations.csv"))
  )

  cor_matrix <- cor(residual_matrix, method = "spearman", use = "pairwise.complete.obs")

  # Keep all correlations in the heat map. The CSV supplies FDR values.
  pdf(
    file.path(out_dir, "05_Correlation",
              paste0(part_name, "_group_adjusted_correlation_heatmap.pdf")),
    width = 10, height = 9
  )
  pheatmap::pheatmap(
    cor_matrix,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    main = paste(part_name, "group-adjusted Spearman correlations")
  )
  dev.off()
}

walk(levels(droplevels(analysis_data$Part)), run_residual_correlation)

# ---------------- 17. OPTIONAL PLASMA TREATMENT-DETAIL ANALYSIS -----------------
# This analysis is metadata-driven. It runs only when the plasma 3_5 group has
# at least two distinct Treatment_detail labels. In the current metadata all
# plasma 3_5 samples are simply labelled "3_5", so this section is skipped.

plasma_35_details <- analysis_data %>%
  filter(Part == "Plasma", Group == "3_5") %>%
  distinct(Treatment_detail) %>%
  filter(!is.na(Treatment_detail), Treatment_detail != "") %>%
  pull(Treatment_detail)

run_plasma_detail_analysis <- length(plasma_35_details) >= 2

write_csv(
  tibble(
    Analysis = "Plasma_3_5_Treatment_detail",
    n_distinct_details = length(plasma_35_details),
    Details = paste(plasma_35_details, collapse = ";"),
    Status = if_else(run_plasma_detail_analysis, "run", "skipped_no_subgroups_in_metadata")
  ),
  file.path(out_dir, "02_Univariate", "plasma_treatment_detail_status.csv")
)

if (run_plasma_detail_analysis) {
  plasma_sensitivity <- analysis_data %>%
    filter(Part == "Plasma") %>%
    mutate(
      Plasma_group = case_when(
        Group == "NC" ~ "NC",
        Group == "MCD" ~ "MCD",
        Group == "3_5" ~ as.character(Treatment_detail),
        TRUE ~ as.character(Group)
      ),
      Plasma_group = factor(Plasma_group)
    )

  plasma_sensitivity_long <- plasma_sensitivity %>%
    select(Sample_ID, Plasma_group, all_of(ba_columns), all_of(pathway_features)) %>%
    pivot_longer(
      all_of(c(ba_columns, pathway_features)),
      names_to = "Feature", values_to = "Value"
    )

  plasma_sensitivity_summary <- plasma_sensitivity_long %>%
    group_by(Plasma_group, Feature) %>%
    summarise(
      n = sum(is.finite(Value)),
      median = median(Value, na.rm = TRUE),
      Q1 = quantile(Value, 0.25, na.rm = TRUE),
      Q3 = quantile(Value, 0.75, na.rm = TRUE),
      .groups = "drop"
    )

  write_csv(
    plasma_sensitivity_summary,
    file.path(out_dir, "02_Univariate", "plasma_treatment_detail_summary.csv")
  )

  plasma_kruskal <- plasma_sensitivity_long %>%
    rename(Group = Plasma_group) %>%
    group_by(Feature) %>%
    nest() %>%
    mutate(test = map(data, safe_kruskal)) %>%
    select(-data) %>%
    unnest(test) %>%
    mutate(q_BH = p.adjust(p_value, "BH"))

  write_csv(
    plasma_kruskal,
    file.path(out_dir, "02_Univariate", "plasma_treatment_detail_kruskal.csv")
  )
}

# ---------------------- 18. OPTIONAL sPLS-DA -----------------------------------
# With only 8 samples per group, interpret supervised classification cautiously.
# Repeated cross-validation is mandatory. The code runs only if mixOmics exists.
#
# IMPORTANT:
# tune.splsda() expects test.keepX to be ONE numeric vector containing the
# candidate numbers of variables. It must not be a list by component.
# A named list is required only for tune.block.splsda()/DIABLO, where each list
# element corresponds to one data block.

if (requireNamespace("mixOmics", quietly = TRUE)) {

  run_splsda <- function(part_name) {
    tryCatch({
      part_data <- analysis_data %>%
        filter(Part == part_name) %>%
        droplevels()

      y <- droplevels(part_data$Group)

      if (nlevels(y) < 2) {
        stop("fewer than two outcome groups")
      }

      group_sizes <- table(y)
      if (min(group_sizes) < 2) {
        stop("at least one group has fewer than two samples")
      }

      x_raw <- as.matrix(part_data[, ba_columns, drop = FALSE])
      storage.mode(x_raw) <- "double"

      # Log-transform each bile acid using a feature-specific positive offset.
      pseudo <- apply(x_raw, 2, minimum_positive)
      x <- log10(sweep(x_raw, 2, pseudo, "+"))

      # Remove features that cannot contribute to the model. This also protects
      # mixOmics from all-NA, non-finite, or zero-variance columns.
      valid_feature <- apply(x, 2, function(z) {
        z <- z[is.finite(z)]
        length(z) == nrow(x) && length(unique(z)) > 1
      })

      x <- x[, valid_feature, drop = FALSE]

      if (ncol(x) < 2) {
        stop("fewer than two nonconstant bile-acid features remain")
      }

      # For three groups, two latent components are ordinarily sufficient.
      ncomp_use <- min(2L, nlevels(y) - 1L, ncol(x))
      if (ncomp_use < 1L) {
        stop("no valid latent component can be fitted")
      }

      # Fast coarse tuning grid: normally c(3, 8, 15).
      # A single-tissue sPLS-DA has only three candidate values rather than seven.
      candidates <- make_fast_keepx_grid(ncol(x))

      # tune.splsda() requires test.keepX to contain more than two entries.
      if (length(candidates) < 3) {
        stop(paste0(
          "test.keepX needs >=3 distinct candidates after filtering; ncol(X) = ",
          ncol(x)
        ))
      }

      folds_use <- min(4L, as.integer(min(group_sizes)))
      if (folds_use < 2L) {
        stop("fewer than two cross-validation folds are possible")
      }

      set.seed(20260807)
      tune_parallel <- get_mixomics_parallel_args(
        mixOmics::tune.splsda, SPLSDA_WORKERS
      )
      tuned <- do.call(
        mixOmics::tune.splsda,
        c(
          list(
            X = x,
            Y = y,
            ncomp = ncomp_use,
            validation = "Mfold",
            folds = folds_use,
            nrepeat = SPLSDA_TUNE_NREPEAT,
            dist = "centroids.dist",
            measure = "BER",
            test.keepX = candidates,
            near.zero.var = TRUE,
            progressBar = FALSE
          ),
          tune_parallel
        )
      )

      optimal_keepX <- as.integer(tuned$choice.keepX)
      optimal_keepX <- rep(optimal_keepX, length.out = ncomp_use)
      optimal_keepX <- pmax(1L, pmin(optimal_keepX, ncol(x)))

      fit <- mixOmics::splsda(
        X = x,
        Y = y,
        ncomp = ncomp_use,
        keepX = optimal_keepX,
        near.zero.var = TRUE
      )

      set.seed(20260807)
      perf_parallel <- get_mixomics_perf_parallel_args(
        fit, SPLSDA_WORKERS
      )
      perf <- do.call(
        mixOmics::perf,
        c(
          list(
            object = fit,
            validation = "Mfold",
            folds = folds_use,
            nrepeat = SPLSDA_PERF_NREPEAT,
            dist = "centroids.dist",
            progressBar = FALSE
          ),
          perf_parallel
        )
      )

      saveRDS(
        list(
          tuning = tuned,
          model = fit,
          performance = perf,
          retained_features = colnames(x),
          removed_features = colnames(x_raw)[!valid_feature],
          test.keepX = candidates,
          optimal.keepX = optimal_keepX,
          folds = folds_use,
          tuning_nrepeat = SPLSDA_TUNE_NREPEAT,
          performance_nrepeat = SPLSDA_PERF_NREPEAT,
          workers = SPLSDA_WORKERS
        ),
        file.path(
          out_dir,
          "04_Multivariate",
          paste0(part_name, "_sPLSDA.rds")
        )
      )

      if (ncomp_use >= 2L) {
        pdf(
          file.path(
            out_dir,
            "04_Multivariate",
            paste0(part_name, "_sPLSDA_samples.pdf")
          ),
          width = 8,
          height = 7
        )
        mixOmics::plotIndiv(
          fit,
          comp = c(1, 2),
          group = y,
          ind.names = FALSE,
          ellipse = TRUE,
          legend = TRUE,
          title = paste(part_name, "sPLS-DA")
        )
        dev.off()
      }

      pdf(
        file.path(
          out_dir,
          "04_Multivariate",
          paste0(part_name, "_sPLSDA_loadings.pdf")
        ),
        width = 9,
        height = 7
      )
      mixOmics::plotLoadings(
        fit,
        comp = 1,
        method = "mean",
        contrib = "max"
      )
      dev.off()

      tibble(
        Part = as.character(part_name),
        Status = "completed",
        Samples = nrow(x),
        Retained_features = ncol(x),
        Removed_features = sum(!valid_feature),
        Components = ncomp_use,
        Folds = folds_use,
        Tuning_repeats = SPLSDA_TUNE_NREPEAT,
        Performance_repeats = SPLSDA_PERF_NREPEAT,
        Workers = SPLSDA_WORKERS,
        test_keepX = paste(candidates, collapse = ";"),
        optimal_keepX = paste(optimal_keepX, collapse = ";"),
        Message = NA_character_
      )

    }, error = function(e) {
      warning(
        paste0("sPLS-DA skipped for ", part_name, ": ", conditionMessage(e)),
        call. = FALSE
      )

      tibble(
        Part = as.character(part_name),
        Status = "skipped_error",
        Samples = NA_integer_,
        Retained_features = NA_integer_,
        Removed_features = NA_integer_,
        Components = NA_integer_,
        Folds = NA_integer_,
        Tuning_repeats = SPLSDA_TUNE_NREPEAT,
        Performance_repeats = SPLSDA_PERF_NREPEAT,
        Workers = SPLSDA_WORKERS,
        test_keepX = NA_character_,
        optimal_keepX = NA_character_,
        Message = conditionMessage(e)
      )
    })
  }

  splsda_status <- map_dfr(
    levels(droplevels(analysis_data$Part)),
    run_splsda
  )

  write_csv(
    splsda_status,
    file.path(out_dir, "04_Multivariate", "sPLSDA_run_status.csv")
  )

} else {
  message("mixOmics is not installed; sPLS-DA and DIABLO sections were skipped.")
}

# --------------------- 19. OPTIONAL MATCHED-ANIMAL MODELS ----------------------
# These sections run only when Animal_ID has been filled in metadata_file.

has_animal_id <- "Animal_ID" %in% names(analysis_data) &&
  sum(!is.na(analysis_data$Animal_ID) & analysis_data$Animal_ID != "") > 0

if (has_animal_id) {
  matched_data <- analysis_data %>%
    filter(!is.na(Animal_ID), Animal_ID != "") %>%
    mutate(
      Animal_ID = factor(Animal_ID),
      Group = droplevels(Group),
      Part = droplevels(Part)
    )

  # CRITICAL QC: one Animal_ID must belong to exactly one experimental group.
  # For sample_metadata_map(3).csv, all 24 Animal_IDs are now group-consistent
  # and each is represented once in all five compartments. Therefore the current
  # dataset should retain all 24 animals (8 NC, 8 MCD, 8 3_5) for matched analyses.
  # The generic exclusion code below is retained as a safety check for future files.
  animal_id_group_qc <- matched_data %>%
    mutate(
      Animal_ID_chr = as.character(Animal_ID),
      Group_chr = as.character(Group)
    ) %>%
    group_by(Animal_ID_chr) %>%
    summarise(
      n_groups = n_distinct(Group_chr),
      Groups = paste(sort(unique(Group_chr)), collapse = ";"),
      n_parts = n_distinct(Part),
      Parts = paste(sort(unique(as.character(Part))), collapse = ";"),
      n_rows = n(),
      .groups = "drop"
    ) %>%
    rename(Animal_ID = Animal_ID_chr)

  write_csv(
    animal_id_group_qc,
    file.path(out_dir, "06_Optional_matched", "AnimalID_group_consistency_QC.csv")
  )

  inconsistent_animal_ids <- animal_id_group_qc %>%
    filter(n_groups != 1) %>%
    pull(Animal_ID)

  if (length(inconsistent_animal_ids) > 0) {
    excluded_matched_ids <- animal_id_group_qc %>%
      filter(Animal_ID %in% inconsistent_animal_ids) %>%
      mutate(Reason = "Animal_ID_has_multiple_Group_assignments_across_tissues")

    write_csv(
      excluded_matched_ids,
      file.path(out_dir, "06_Optional_matched", "excluded_from_matched_analyses.csv")
    )

    message(
      paste0(
        "Matched analyses use only group-consistent Animal_IDs. Excluding ",
        length(inconsistent_animal_ids), " IDs: ",
        paste(inconsistent_animal_ids, collapse = ", "), "."
      )
    )
  }

  matched_data <- matched_data %>%
    filter(!as.character(Animal_ID) %in% inconsistent_animal_ids) %>%
    mutate(
      Animal_ID = droplevels(factor(Animal_ID)),
      Group = droplevels(Group),
      Part = droplevels(Part)
    )

  # Retain only IDs represented in all five compartments for cross-tissue models.
  required_matched_parts <- c("Liver", "Plasma", "Ileum", "Cecum", "Feces")
  complete_matched_ids <- matched_data %>%
    mutate(Part_chr = as.character(Part), Animal_ID_chr = as.character(Animal_ID)) %>%
    group_by(Animal_ID_chr) %>%
    summarise(
      n_parts = n_distinct(Part_chr),
      all_required_parts = all(required_matched_parts %in% unique(Part_chr)),
      .groups = "drop"
    ) %>%
    filter(n_parts == length(required_matched_parts), all_required_parts) %>%
    pull(Animal_ID_chr)

  matched_data <- matched_data %>%
    filter(as.character(Animal_ID) %in% complete_matched_ids) %>%
    mutate(
      Animal_ID = droplevels(factor(Animal_ID)),
      Group = droplevels(Group),
      Part = droplevels(Part)
    )

  matched_analysis_counts <- matched_data %>%
    distinct(Animal_ID, Group) %>%
    count(Group, name = "n_animals")

  write_csv(
    matched_analysis_counts,
    file.path(out_dir, "06_Optional_matched", "matched_analysis_animal_counts.csv")
  )

  # Expected design for the revised data/metadata supplied with this script:
  # 24 complete animals = 8 NC + 8 MCD + 8 3_5, each with five compartments.
  expected_group_counts <- tibble(
    Group = factor(c("NC", "MCD", "3_5"), levels = c("NC", "MCD", "3_5")),
    expected_n_animals = c(8L, 8L, 8L)
  )

  matched_design_check <- matched_analysis_counts %>%
    mutate(Group = factor(as.character(Group), levels = c("NC", "MCD", "3_5"))) %>%
    full_join(expected_group_counts, by = "Group") %>%
    mutate(
      n_animals = replace_na(n_animals, 0L),
      expected_n_animals = replace_na(expected_n_animals, 0L),
      matches_expected = n_animals == expected_n_animals
    ) %>%
    arrange(Group)

  write_csv(
    matched_design_check,
    file.path(out_dir, "06_Optional_matched", "matched_design_expected_vs_observed.csv")
  )

  if (n_distinct(matched_data$Animal_ID) != 24L ||
      any(!matched_design_check$matches_expected)) {
    warning(
      "Matched-data design differs from the expected revised dataset ",
      "(24 animals; 8 NC, 8 MCD, 8 3_5). See matched-design QC outputs."
    )
  }

  matched_summary <- tibble(
    n_input_samples = nrow(analysis_data),
    n_unique_Animal_ID_input = n_distinct(analysis_data$Animal_ID[!is.na(analysis_data$Animal_ID)]),
    n_group_inconsistent_ID = length(inconsistent_animal_ids),
    n_complete_group_consistent_ID = n_distinct(matched_data$Animal_ID),
    n_matched_rows = nrow(matched_data),
    Expected_parts_per_animal = length(required_matched_parts),
    Analysis_note = "Current data: all 24 Animal_IDs should be complete and group-consistent; QC still filters any future inconsistency"
  )
  write_csv(
    matched_summary,
    file.path(out_dir, "06_Optional_matched", "matched_analysis_summary.csv")
  )

  # Warn about duplicated tissue samples per animal. These may be technical
  # replicates and should be resolved before longitudinal/multiblock analysis.
  duplicated_animal_part <- matched_data %>%
    count(Animal_ID, Part) %>%
    filter(n > 1)

  write_csv(
    duplicated_animal_part,
    file.path(out_dir, "06_Optional_matched", "duplicated_AnimalID_by_Part.csv")
  )

  # ---------------- Mixed-effects models for each individual BA ----------------
  if (requireNamespace("lme4", quietly = TRUE) &&
      requireNamespace("lmerTest", quietly = TRUE) &&
      requireNamespace("emmeans", quietly = TRUE)) {

    fit_mixed_one_feature <- function(feature_name) {
      model_data <- matched_data %>%
        select(Animal_ID, Group, Part, Batch, all_of(feature_name)) %>%
        rename(Value = all_of(feature_name)) %>%
        filter(is.finite(Value))

      pseudo <- minimum_positive(model_data$Value)
      model_data <- model_data %>% mutate(Log2Value = log2(Value + pseudo))

      # Batch is included only when it contains at least two nonmissing levels.
      use_batch <- "Batch" %in% names(model_data) &&
        n_distinct(na.omit(model_data$Batch)) >= 2

      formula_used <- if (use_batch) {
        Log2Value ~ Group * Part + Batch + (1 | Animal_ID)
      } else {
        Log2Value ~ Group * Part + (1 | Animal_ID)
      }

      fit <- tryCatch(
        suppressMessages(
          suppressWarnings(
            lmerTest::lmer(
              formula_used,
              data = model_data,
              REML = FALSE,
              control = lme4::lmerControl(
                check.conv.singular = list(
                  action = "ignore",
                  tol = 1e-4
                )
              )
            )
          )
        ),
        error = function(e) NULL
      )
      if (is.null(fit)) return(tibble())

      singular_fit <- lme4::isSingular(fit, tol = 1e-4)
      model_status <- if (singular_fit) {
        "lmer_singular_random_effect_variance_near_zero"
      } else {
        "lmer_ok"
      }

      # The singular-fit convergence message is deliberately suppressed at
      # model-fitting time via lmerControl().  We still diagnose it explicitly
      # with isSingular() and retain that information in Model_status.

      emm <- emmeans::emmeans(fit, ~Group | Part)
      contrasts_out <- emmeans::contrast(
        emm,
        method = list(
          MCD_vs_NC = c(-1, 1, 0),
          `3_5_vs_MCD` = c(0, -1, 1),
          `3_5_vs_NC` = c(-1, 0, 1)
        )
      ) %>%
        as.data.frame() %>%
        as_tibble() %>%
        transmute(
          Feature = feature_name,
          Part,
          Contrast = contrast,
          estimate_log2FC = estimate,
          SE, df,
          statistic = t.ratio,
          p_value = p.value,
          Pseudocount = pseudo,
          Model_status = model_status
        )
      contrasts_out
    }

    mixed_results <- map_dfr(ba_columns, fit_mixed_one_feature) %>%
      group_by(Part, Contrast) %>%
      mutate(q_BH = p.adjust(p_value, "BH")) %>%
      ungroup()

    write_csv(
      mixed_results,
      file.path(out_dir, "06_Optional_matched", "mixed_effects_BA_contrasts.csv")
    )
  }

  # ---------------- Beta mixed models for fractions -----------------------------
  if (requireNamespace("glmmTMB", quietly = TRUE) &&
      requireNamespace("emmeans", quietly = TRUE)) {

    proportion_features <- c(
      "Conjugated_fraction", "Secondary_fraction", "Fraction_12aOH",
      "CA_conjugation", "CDCA_conjugation", "DCA_conjugation",
      "UDCA_conjugation", "alphaMCA_conjugation", "betaMCA_conjugation",
      "CA_to_DCA_index", "CDCA_to_UDCA_index",
      "betaMCA_to_omegaMCA_index"
    )

    fit_beta_mixed <- function(feature_name) {
      d <- matched_data %>%
        select(Animal_ID, Group, Part, all_of(feature_name)) %>%
        rename(Proportion = all_of(feature_name)) %>%
        filter(is.finite(Proportion), Proportion >= 0, Proportion <= 1)

      # Smithson-Verkuilen adjustment moves exact 0/1 into the open interval.
      n <- nrow(d)
      d <- d %>% mutate(Proportion_beta = (Proportion * (n - 1) + 0.5) / n)

      fit <- tryCatch(
        glmmTMB::glmmTMB(
          Proportion_beta ~ Group * Part + (1 | Animal_ID),
          family = glmmTMB::beta_family(link = "logit"),
          data = d
        ),
        error = function(e) NULL
      )
      if (is.null(fit)) return(tibble())

      emm <- emmeans::emmeans(fit, ~Group | Part, type = "response")
      pairs <- emmeans::contrast(
        emm,
        method = list(
          MCD_vs_NC = c(-1, 1, 0),
          `3_5_vs_MCD` = c(0, -1, 1),
          `3_5_vs_NC` = c(-1, 0, 1)
        ),
        type = "response"
      ) %>%
        as.data.frame() %>%
        as_tibble() %>%
        mutate(Feature = feature_name)
      pairs
    }

    beta_mixed_results <- map_dfr(proportion_features, fit_beta_mixed)
    write_csv(
      beta_mixed_results,
      file.path(out_dir, "06_Optional_matched", "beta_mixed_pathway_indices.csv")
    )
  }

  # ---------------- Cross-tissue correlations ----------------------------------
  # Only animals represented once in both compared compartments are used.
  cross_tissue_features <- c(
    "Conjugated_fraction", "Secondary_fraction", "CA_to_DCA_index",
    "CDCA_to_UDCA_index", "betaMCA_to_omegaMCA_index", "Total_BA"
  )

  matched_index_long <- matched_data %>%
    select(Animal_ID, Group, Part, all_of(c(ba_columns, cross_tissue_features))) %>%
    pivot_longer(
      all_of(c(ba_columns, cross_tissue_features)),
      names_to = "Feature", values_to = "Value"
    )

  part_pairs <- combn(levels(droplevels(matched_data$Part)), 2, simplify = FALSE)

  cross_tissue_correlations <- map_dfr(part_pairs, function(parts_pair) {
    p1 <- parts_pair[1]
    p2 <- parts_pair[2]

    wide <- matched_index_long %>%
      filter(Part %in% c(p1, p2)) %>%
      group_by(Animal_ID, Group, Part, Feature) %>%
      summarise(Value = mean(Value, na.rm = TRUE), .groups = "drop") %>%
      mutate(Part_feature = paste(Part, Feature, sep = "__")) %>%
      select(Animal_ID, Group, Part_feature, Value) %>%
      pivot_wider(names_from = Part_feature, values_from = Value)

    map_dfr(c(ba_columns, cross_tissue_features), function(feature_name) {
      col1 <- paste(p1, feature_name, sep = "__")
      col2 <- paste(p2, feature_name, sep = "__")
      if (!all(c(col1, col2) %in% names(wide))) return(tibble())

      x <- wide[[col1]]
      y <- wide[[col2]]
      keep <- is.finite(x) & is.finite(y)
      if (sum(keep) < 5) return(tibble())

      test <- suppressWarnings(cor.test(x[keep], y[keep], method = "spearman", exact = FALSE))
      tibble(
        Part1 = p1, Part2 = p2, Feature = feature_name,
        n_pairs = sum(keep), rho = unname(test$estimate), p_value = test$p.value
      )
    })
  }) %>%
    group_by(Part1, Part2) %>%
    mutate(q_BH = p.adjust(p_value, "BH")) %>%
    ungroup()

  write_csv(
    cross_tissue_correlations,
    file.path(out_dir, "06_Optional_matched", "cross_tissue_spearman_correlations.csv")
  )

  # ---------------- DIABLO multiblock sPLS-DA ----------------------------------
  if (requireNamespace("mixOmics", quietly = TRUE)) {
    # DIABLO requires exactly the same animals, in exactly the same order,
    # in every block and in Y. With sample_metadata_map(3).csv the expected
    # design is 24 animals (8/group) present in all five blocks.

    summarized <- matched_data %>%
      group_by(Animal_ID, Group, Part) %>%
      summarise(
        across(all_of(ba_columns), ~ mean(.x, na.rm = TRUE)),
        .groups = "drop"
      )

    required_parts <- required_matched_parts
    n_required_parts <- length(required_parts)

    complete_id_qc <- summarized %>%
      mutate(Animal_ID_chr = as.character(Animal_ID)) %>%
      group_by(Animal_ID_chr) %>%
      summarise(
        n_parts = n_distinct(Part),
        n_groups = n_distinct(Group),
        n_rows = n(),
        .groups = "drop"
      ) %>%
      rename(Animal_ID = Animal_ID_chr)

    write_csv(
      complete_id_qc,
      file.path(out_dir, "06_Optional_matched", "DIABLO_complete_ID_QC.csv")
    )

    complete_ids <- complete_id_qc %>%
      filter(
        n_parts == n_required_parts,
        n_groups == 1,
        n_rows == n_required_parts
      ) %>%
      pull(Animal_ID)

    complete_data <- summarized %>%
      filter(as.character(Animal_ID) %in% complete_ids) %>%
      mutate(
        Animal_ID = as.character(Animal_ID),
        Group = droplevels(Group),
        Part = droplevels(Part)
      )

    group_counts_diablo <- complete_data %>%
      distinct(Animal_ID, Group) %>%
      count(Group, name = "n")

    write_csv(
      group_counts_diablo,
      file.path(out_dir, "06_Optional_matched", "DIABLO_group_counts.csv")
    )

    can_run_diablo <- length(complete_ids) >= 12 &&
      n_distinct(complete_data$Group) == 3 &&
      nrow(group_counts_diablo) == 3 &&
      min(group_counts_diablo$n) >= 3

    if (can_run_diablo) {
      # Build Y first; this becomes the master order for every block.
      y_table <- complete_data %>%
        distinct(Animal_ID, Group) %>%
        arrange(Animal_ID)

      if (anyDuplicated(y_table$Animal_ID)) {
        stop("DIABLO aborted: Animal_ID is not unique in Y after QC filtering.")
      }

      id_order <- y_table$Animal_ID
      y <- droplevels(factor(y_table$Group, levels = c("NC", "MCD", "3_5")))
      block_names <- required_parts
      removed_diablo_features <- list()

      blocks <- set_names(map(block_names, function(part_name) {
        d <- complete_data %>%
          filter(Part == part_name) %>%
          arrange(match(Animal_ID, id_order))

        if (nrow(d) != length(id_order)) {
          stop(paste0(
            "DIABLO alignment error for ", part_name, ": block has ",
            nrow(d), " rows but Y has ", length(id_order), "."
          ))
        }

        if (!identical(as.character(d$Animal_ID), as.character(id_order))) {
          stop(paste0(
            "DIABLO alignment error for ", part_name,
            ": Animal_ID order does not match Y."
          ))
        }

        x_raw <- as.matrix(d[, ba_columns, drop = FALSE])
        storage.mode(x_raw) <- "numeric"
        pseudo <- apply(x_raw, 2, minimum_positive)
        x <- log10(sweep(x_raw, 2, pseudo, "+"))

        valid_feature <- vapply(seq_len(ncol(x)), function(j) {
          z <- x[, j]
          all(is.finite(z)) && length(unique(z)) > 1 &&
            is.finite(sd(z)) && sd(z) > 0
        }, logical(1))

        removed_diablo_features[[part_name]] <<- colnames(x)[!valid_feature]
        x <- x[, valid_feature, drop = FALSE]

        if (ncol(x) < 2) {
          stop(paste0(
            "DIABLO skipped: fewer than two usable features in ", part_name, "."
          ))
        }

        rownames(x) <- d$Animal_ID
        x
      }), block_names)

      block_n <- vapply(blocks, nrow, integer(1))
      block_order_ok <- vapply(
        blocks,
        function(z) identical(rownames(z), as.character(id_order)),
        logical(1)
      )

      if (any(block_n != length(y)) || !all(block_order_ok)) {
        stop(paste0(
          "DIABLO aborted after final alignment check. Block sizes: ",
          paste(names(block_n), block_n, sep = "=", collapse = ", "),
          "; Y=", length(y), "."
        ))
      }

      write_csv(
        y_table,
        file.path(out_dir, "06_Optional_matched", "DIABLO_samples_used.csv")
      )

      removed_diablo_tbl <- imap_dfr(
        removed_diablo_features,
        ~ tibble(Part = .y, Feature = .x)
      )
      write_csv(
        removed_diablo_tbl,
        file.path(out_dir, "06_Optional_matched", "DIABLO_removed_features.csv")
      )

      design <- matrix(
        0.1,
        nrow = length(blocks),
        ncol = length(blocks),
        dimnames = list(names(blocks), names(blocks))
      )
      diag(design) <- 0

      # Fast DIABLO grid: three keepX values per block, normally c(3, 8, 15).
      # With five blocks this means 3^5 = 243 combinations per component/repeat,
      # instead of 7^5 = 16,807 combinations with the previous grid.
      test_keepX <- imap(blocks, function(x, block_name) {
        candidates <- make_fast_keepx_grid(ncol(x))
        if (length(candidates) < 3) {
          stop(
            "DIABLO test.keepX needs >=3 candidate values for block ",
            block_name, "; only ", length(candidates), " available."
          )
        }
        candidates
      })

      folds_diablo <- min(4L, as.integer(min(group_counts_diablo$n)))
      if (folds_diablo < 2L) {
        stop("DIABLO skipped: fewer than two CV folds are possible.")
      }

      # Report the size of the tuning problem before it starts.
      combinations_per_component <- prod(vapply(test_keepX, length, integer(1)))
      message(
        "DIABLO fast tuning: ", combinations_per_component,
        " keepX combinations/component x ", DIABLO_TUNE_NREPEAT,
        " repeat(s), using up to ", DIABLO_WORKERS, " worker(s)."
      )

      set.seed(20260807)
      diablo_tune_parallel <- get_mixomics_parallel_args(
        mixOmics::tune.block.splsda, DIABLO_WORKERS
      )
      tuned_diablo <- do.call(
        mixOmics::tune.block.splsda,
        c(
          list(
            X = blocks,
            Y = y,
            ncomp = 2,
            test.keepX = test_keepX,
            design = design,
            validation = "Mfold",
            folds = folds_diablo,
            nrepeat = DIABLO_TUNE_NREPEAT,
            dist = "centroids.dist",
            progressBar = FALSE
          ),
          diablo_tune_parallel
        )
      )

      list_keepX <- tuned_diablo$choice.keepX

      diablo_fit <- mixOmics::block.splsda(
        X = blocks,
        Y = y,
        ncomp = 2,
        keepX = list_keepX,
        design = design
      )

      # More thorough repeated CV is reserved for the final selected model.
      set.seed(20260807)
      diablo_perf_parallel <- get_mixomics_perf_parallel_args(
        diablo_fit, DIABLO_WORKERS
      )
      diablo_perf <- tryCatch(
        do.call(
          mixOmics::perf,
          c(
            list(
              object = diablo_fit,
              validation = "Mfold",
              folds = folds_diablo,
              nrepeat = DIABLO_PERF_NREPEAT,
              dist = "centroids.dist",
              progressBar = FALSE
            ),
            diablo_perf_parallel
          )
        ),
        error = function(e) {
          warning(
            paste0("DIABLO final performance assessment skipped: ", conditionMessage(e)),
            call. = FALSE
          )
          NULL
        }
      )

      saveRDS(
        list(
          tuning = tuned_diablo,
          model = diablo_fit,
          performance = diablo_perf,
          Animal_ID_order = id_order,
          Y = y,
          test.keepX = test_keepX,
          optimal.keepX = list_keepX,
          removed_features = removed_diablo_features,
          folds = folds_diablo,
          tuning_nrepeat = DIABLO_TUNE_NREPEAT,
          performance_nrepeat = DIABLO_PERF_NREPEAT,
          workers = DIABLO_WORKERS,
          combinations_per_component = combinations_per_component
        ),
        file.path(out_dir, "06_Optional_matched", "DIABLO_model.rds")
      )

      pdf(
        file.path(out_dir, "06_Optional_matched", "DIABLO_samples.pdf"),
        width = 9, height = 7
      )
      mixOmics::plotIndiv(
        diablo_fit,
        ind.names = FALSE,
        legend = TRUE,
        title = "DIABLO multiblock BA profile"
      )
      dev.off()

      # -----------------------------------------------------------------------
      # Robust DIABLO circos export
      #
      # IMPORTANT: circosPlot() expects color.blocks to be a character vector
      # (one colour per block), not TRUE/FALSE.  The older code used
      # color.blocks = TRUE inside try(..., silent = TRUE); therefore a plotting
      # error could be hidden while an empty/unreadable PDF was still created.
      # -----------------------------------------------------------------------
      circos_pdf <- file.path(
        out_dir, "06_Optional_matched", "DIABLO_circos.pdf"
      )
      circos_png <- file.path(
        out_dir, "06_Optional_matched", "DIABLO_circos.png"
      )
      circos_diag <- file.path(
        out_dir, "06_Optional_matched", "DIABLO_circos_diagnostic.txt"
      )

      block_cols <- grDevices::hcl.colors(
        length(blocks), palette = "Dark 3"
      )
      group_cols <- grDevices::hcl.colors(
        length(unique(as.character(y))), palette = "Set 2"
      )

      draw_diablo_circos <- function() {
        mixOmics::circosPlot(
          diablo_fit,
          comp = seq_len(min(2L, diablo_fit$ncomp)),
          cutoff = 0.7,
          line = TRUE,
          color.blocks = block_cols,
          color.Y = group_cols,
          color.cor = c("firebrick3", "navy"),
          size.labels = 1.1,
          size.variables = 0.5,
          legend = TRUE
        )
      }

      circos_error <- NULL
      circos_result <- NULL

      # Draw first to a temporary PDF.  Only copy it to the final filename if
      # circosPlot finishes successfully.
      tmp_circos_pdf <- tempfile(fileext = ".pdf")
      grDevices::pdf(tmp_circos_pdf, width = 10, height = 10, onefile = TRUE)
      tryCatch(
        {
          circos_result <- draw_diablo_circos()
        },
        error = function(e) {
          circos_error <<- conditionMessage(e)
        },
        finally = {
          grDevices::dev.off()
        }
      )

      if (is.null(circos_error) &&
          file.exists(tmp_circos_pdf) &&
          is.finite(file.info(tmp_circos_pdf)$size) &&
          file.info(tmp_circos_pdf)$size > 1000) {

        file.copy(tmp_circos_pdf, circos_pdf, overwrite = TRUE)

        # Also create a PNG fallback, which is often easier to open on Windows.
        png_error <- NULL
        grDevices::png(
          circos_png, width = 3000, height = 3000, res = 300
        )
        tryCatch(
          {
            draw_diablo_circos()
          },
          error = function(e) {
            png_error <<- conditionMessage(e)
          },
          finally = {
            grDevices::dev.off()
          }
        )

        diag_lines <- c(
          "DIABLO circos export completed.",
          "cutoff = 0.7",
          paste0("blocks = ", paste(names(blocks), collapse = ", ")),
          paste0("PDF = ", circos_pdf),
          paste0("PNG = ", circos_png),
          if (!is.null(png_error)) paste0("PNG warning: ", png_error) else "PNG export completed."
        )
        writeLines(diag_lines, circos_diag)

        # circosPlot returns useful plotting information as attributes/results
        # in many mixOmics versions; preserve it for later inspection.
        saveRDS(
          circos_result,
          file.path(out_dir, "06_Optional_matched", "DIABLO_circos_result.rds")
        )

      } else {
        # Never leave behind a corrupted/empty final PDF.  Instead create a
        # valid diagnostic PDF and write the actual error to a text file.
        if (file.exists(circos_pdf)) unlink(circos_pdf)
        if (file.exists(circos_png)) unlink(circos_png)

        grDevices::pdf(circos_pdf, width = 10, height = 10)
        graphics::plot.new()
        graphics::text(
          0.5, 0.58,
          "DIABLO circos plot could not be generated",
          cex = 1.4, font = 2
        )
        graphics::text(
          0.5, 0.48,
          paste0("mixOmics error: ",
                 ifelse(is.null(circos_error),
                        "unknown plotting/device error",
                        circos_error)),
          cex = 0.9
        )
        graphics::text(
          0.5, 0.39,
          "See DIABLO_circos_diagnostic.txt for details.",
          cex = 0.9
        )
        grDevices::dev.off()

        writeLines(
          c(
            "DIABLO circos export FAILED.",
            "cutoff = 0.7",
            paste0(
              "Error: ",
              ifelse(is.null(circos_error),
                     "unknown plotting/device error",
                     circos_error)
            ),
            paste0("mixOmics version: ",
                   as.character(utils::packageVersion("mixOmics"))),
            paste0("R version: ", R.version.string)
          ),
          circos_diag
        )

        warning(
          "DIABLO circos plotting failed. See: ", circos_diag,
          call. = FALSE
        )
      }

      if (file.exists(tmp_circos_pdf)) unlink(tmp_circos_pdf)
    } else {
      message(
        "DIABLO skipped: need >=12 group-consistent complete animals, all 3 groups, and >=3 animals per group."
      )
    }
  }

  # ---------------- Piecewise structural-equation model ------------------------
  # This section ALWAYS creates a summary/diagnostic file.  The older version
  # only created piecewise_SEM_summary.txt when piecewiseSEM was installed and
  # the model ran successfully, so a skipped/error state could look like a
  # "missing result".  We now:
  #   1) always write the SEM input data and diagnostics;
  #   2) always fit/export the four component LM path models as a fallback;
  #   3) run piecewiseSEM::psem() when the package is available;
  #   4) save standardized path coefficients and the fitted psem object;
  #   5) always create piecewise_SEM_summary.txt, even on failure/skip.

  sem_dir <- file.path(out_dir, "06_Optional_matched")
  sem_summary_file <- file.path(sem_dir, "piecewise_SEM_summary.txt")
  sem_diag_file <- file.path(sem_dir, "piecewise_SEM_diagnostic.txt")
  sem_input_file <- file.path(sem_dir, "piecewise_SEM_input_data.csv")
  sem_lm_coef_file <- file.path(sem_dir, "piecewise_SEM_component_LM_coefficients.csv")
  sem_lm_std_file <- file.path(sem_dir, "piecewise_SEM_component_LM_standardized_coefficients.csv")
  sem_psem_coef_file <- file.path(sem_dir, "piecewise_SEM_path_coefficients.csv")
  sem_model_file <- file.path(sem_dir, "piecewise_SEM_model.rds")

  # Create a placeholder immediately so the expected file can never be absent.
  writeLines(
    c(
      "Piecewise SEM analysis initialized.",
      paste0("Date/time: ", Sys.time()),
      "The file will be overwritten with the final result below."
    ),
    sem_summary_file
  )

  # Manual feature selection is defined near the top of the script in
  # SEM_FEATURE_BY_PART.  The five compartments must be present exactly once.
  sem_parts <- c("Liver", "Plasma", "Ileum", "Cecum", "Feces")

  if (!identical(sort(names(SEM_FEATURE_BY_PART)), sort(sem_parts))) {
    stop(
      "SEM_FEATURE_BY_PART must have exactly these names: ",
      paste(sem_parts, collapse = ", ")
    )
  }

  sem_selected_features <- unique(unname(SEM_FEATURE_BY_PART[sem_parts]))

  missing_source_features <- setdiff(sem_selected_features, names(matched_data))
  if (length(missing_source_features) > 0L) {
    stop(
      "Selected SEM feature(s) not found in matched_data: ",
      paste(missing_source_features, collapse = ", "),
      "\nCheck spelling/capitalization. Raw T.b.MCA is called TbetaMCA in this script."
    )
  }

  nonnumeric_sem_features <- sem_selected_features[
    !vapply(matched_data[sem_selected_features], is.numeric, logical(1))
  ]
  if (length(nonnumeric_sem_features) > 0L) {
    stop(
      "Selected SEM feature(s) must be numeric: ",
      paste(nonnumeric_sem_features, collapse = ", ")
    )
  }

  sem_wide <- matched_data %>%
    group_by(Animal_ID, Group, Part) %>%
    summarise(
      across(
        all_of(sem_selected_features),
        ~ mean(.x, na.rm = TRUE)
      ),
      .groups = "drop"
    ) %>%
    pivot_wider(
      names_from = Part,
      values_from = all_of(sem_selected_features),
      names_glue = "{.value}__{Part}"
    ) %>%
    mutate(
      Animal_ID = as.character(Animal_ID),
      Group = droplevels(factor(Group, levels = c("NC", "MCD", "3_5")))
    )

  # Convert each manual Part -> Feature choice to its corresponding wide-column
  # name, e.g. Liver="TbetaMCA" becomes "TbetaMCA__Liver".
  sem_node_vars <- setNames(
    paste0(unname(SEM_FEATURE_BY_PART[sem_parts]), "__", sem_parts),
    sem_parts
  )
  required_sem_vars <- unname(sem_node_vars)

  # Export the exact manual selection used for this SEM run.
  sem_selection_table <- tibble(
    Part = sem_parts,
    Feature = unname(SEM_FEATURE_BY_PART[sem_parts]),
    SEM_variable = required_sem_vars
  )
  readr::write_csv(
    sem_selection_table,
    file.path(sem_dir, "piecewise_SEM_manual_feature_selection.csv")
  )

  missing_sem_vars <- setdiff(required_sem_vars, names(sem_wide))

  if (length(missing_sem_vars) == 0L) {
    sem_wide <- sem_wide %>%
      filter(if_all(all_of(required_sem_vars), is.finite)) %>%
      droplevels()
  }

  readr::write_csv(sem_wide, sem_input_file)

  sem_group_counts <- if (nrow(sem_wide) > 0L && "Group" %in% names(sem_wide)) {
    sem_wide %>% count(Group, name = "n")
  } else {
    tibble(Group = character(), n = integer())
  }
  readr::write_csv(
    sem_group_counts,
    file.path(sem_dir, "piecewise_SEM_group_counts.csv")
  )

  sem_diagnostics <- c(
    paste0("piecewiseSEM installed: ", requireNamespace("piecewiseSEM", quietly = TRUE)),
    paste0(
      "piecewiseSEM version: ",
      if (requireNamespace("piecewiseSEM", quietly = TRUE)) {
        as.character(utils::packageVersion("piecewiseSEM"))
      } else {
        "not installed"
      }
    ),
    paste0("Complete animals available for SEM: ", nrow(sem_wide)),
    paste0(
      "Manual SEM chain: ",
      paste(
        paste0(sem_parts, "[", unname(SEM_FEATURE_BY_PART[sem_parts]), "]"),
        collapse = " -> "
      )
    ),
    paste0(
      "Group counts: ",
      if (nrow(sem_group_counts) > 0L) {
        paste(sem_group_counts$Group, sem_group_counts$n, sep = "=", collapse = ", ")
      } else {
        "none"
      }
    ),
    paste0(
      "Missing required SEM variables: ",
      if (length(missing_sem_vars) == 0L) "none" else paste(missing_sem_vars, collapse = ", ")
    ),
    paste0("R version: ", R.version.string)
  )

  # Safe component models. Formulas are generated automatically from the
  # manual feature selection above. Each path still adjusts for experimental Group.
  make_sem_formula <- function(response_part, predictor_part) {
    stats::reformulate(
      termlabels = c(sem_node_vars[[predictor_part]], "Group"),
      response = sem_node_vars[[response_part]]
    )
  }

  sem_formulas <- list(
    Liver_to_Plasma = make_sem_formula("Plasma", "Liver"),
    Plasma_to_Ileum = make_sem_formula("Ileum", "Plasma"),
    Ileum_to_Cecum = make_sem_formula("Cecum", "Ileum"),
    Cecum_to_Feces = make_sem_formula("Feces", "Cecum")
  )

  sem_lm_models <- list()
  sem_lm_errors <- character()

  if (length(missing_sem_vars) == 0L && nrow(sem_wide) >= 6L && nlevels(sem_wide$Group) >= 2L) {
    for (nm in names(sem_formulas)) {
      fit_tmp <- tryCatch(
        stats::lm(sem_formulas[[nm]], data = sem_wide),
        error = function(e) e
      )
      if (inherits(fit_tmp, "error")) {
        sem_lm_errors[nm] <- conditionMessage(fit_tmp)
      } else {
        sem_lm_models[[nm]] <- fit_tmp
      }
    }
  } else {
    sem_lm_errors["component_models"] <- paste0(
      "Insufficient valid SEM data: n=", nrow(sem_wide),
      ", group levels=", if ("Group" %in% names(sem_wide)) nlevels(sem_wide$Group) else 0L,
      if (length(missing_sem_vars) > 0L) {
        paste0("; missing variables: ", paste(missing_sem_vars, collapse = ", "))
      } else ""
    )
  }

  sem_lm_coef <- purrr::imap_dfr(sem_lm_models, function(model, path_name) {
    broom::tidy(model, conf.int = TRUE) %>%
      mutate(
        Path = path_name,
        n = stats::nobs(model),
        r_squared = summary(model)$r.squared,
        adjusted_r_squared = summary(model)$adj.r.squared,
        .before = 1
      )
  })

  if (nrow(sem_lm_coef) == 0L) {
    sem_lm_coef <- tibble(
      Path = character(), term = character(), estimate = numeric(),
      std.error = numeric(), statistic = numeric(), p.value = numeric(),
      conf.low = numeric(), conf.high = numeric(), n = integer(),
      r_squared = numeric(), adjusted_r_squared = numeric()
    )
  }
  readr::write_csv(sem_lm_coef, sem_lm_coef_file)

  # Standardized component path models: numeric SEM variables are z-scored,
  # while Group remains a factor. For the biological chain predictor in each
  # equation, the estimate is therefore a standardized path coefficient.
  sem_std <- sem_wide
  if (length(missing_sem_vars) == 0L && nrow(sem_std) > 0L) {
    sem_std <- sem_std %>%
      mutate(
        across(
          all_of(required_sem_vars),
          ~ {
            sx <- stats::sd(.x, na.rm = TRUE)
            if (is.finite(sx) && sx > 0) as.numeric(scale(.x)) else .x
          }
        )
      )
  }

  sem_lm_std_models <- list()
  if (length(missing_sem_vars) == 0L && nrow(sem_std) >= 6L && nlevels(sem_std$Group) >= 2L) {
    for (nm in names(sem_formulas)) {
      fit_tmp <- tryCatch(
        stats::lm(sem_formulas[[nm]], data = sem_std),
        error = function(e) NULL
      )
      if (!is.null(fit_tmp)) sem_lm_std_models[[nm]] <- fit_tmp
    }
  }

  sem_lm_std_coef <- purrr::imap_dfr(sem_lm_std_models, function(model, path_name) {
    broom::tidy(model, conf.int = TRUE) %>%
      mutate(
        Path = path_name,
        n = stats::nobs(model),
        .before = 1
      )
  })
  if (nrow(sem_lm_std_coef) == 0L) {
    sem_lm_std_coef <- tibble(
      Path = character(), term = character(), estimate = numeric(),
      std.error = numeric(), statistic = numeric(), p.value = numeric(),
      conf.low = numeric(), conf.high = numeric(), n = integer()
    )
  }
  readr::write_csv(sem_lm_std_coef, sem_lm_std_file)

  psem_status <- "not_run"
  psem_error <- NULL
  sem_model <- NULL
  sem_summary_lines <- character()

  can_run_psem <- length(missing_sem_vars) == 0L &&
    nrow(sem_wide) >= 15L &&
    nlevels(sem_wide$Group) == 3L &&
    length(sem_lm_models) == length(sem_formulas)

  if (!requireNamespace("piecewiseSEM", quietly = TRUE)) {
    psem_status <- "piecewiseSEM_not_installed_component_LM_fallback_completed"
    psem_error <- paste0(
      "The package 'piecewiseSEM' is not installed. The four component path ",
      "regressions were still fitted and exported. To obtain Fisher's C and ",
      "piecewiseSEM-specific global fit statistics, install piecewiseSEM and ",
      "rerun this SEM section."
    )
  } else if (!can_run_psem) {
    psem_status <- "piecewiseSEM_skipped_insufficient_or_invalid_input"
    psem_error <- paste0(
      "Requirements not met. Need >=15 complete animals, all 3 groups, all ",
      "required variables, and four successfully fitted component models. ",
      "Observed n=", nrow(sem_wide),
      ", group levels=", if ("Group" %in% names(sem_wide)) nlevels(sem_wide$Group) else 0L,
      ", component models fitted=", length(sem_lm_models), "."
    )
  } else {
    sem_model <- tryCatch(
      do.call(piecewiseSEM::psem, unname(sem_lm_models)),
      error = function(e) {
        psem_error <<- conditionMessage(e)
        NULL
      }
    )

    if (!is.null(sem_model)) {
      sem_summary_obj <- tryCatch(
        summary(sem_model),
        error = function(e) {
          psem_error <<- conditionMessage(e)
          NULL
        }
      )

      if (!is.null(sem_summary_obj)) {
        psem_status <- "piecewiseSEM_completed"
        sem_summary_lines <- capture.output(print(sem_summary_obj))
        saveRDS(sem_model, sem_model_file)

        psem_coef <- tryCatch(
          piecewiseSEM::coefs(sem_model, standardize = "scale"),
          error = function(e) {
            sem_diagnostics <<- c(
              sem_diagnostics,
              paste0("piecewiseSEM::coefs warning: ", conditionMessage(e))
            )
            NULL
          }
        )

        if (!is.null(psem_coef)) {
          readr::write_csv(as.data.frame(psem_coef), sem_psem_coef_file)
        } else {
          readr::write_csv(
            tibble(Status = "piecewiseSEM_coefficient_extraction_failed"),
            sem_psem_coef_file
          )
        }
      } else {
        psem_status <- "piecewiseSEM_summary_failed_component_LM_fallback_available"
      }
    } else {
      psem_status <- "piecewiseSEM_fit_failed_component_LM_fallback_available"
    }
  }

  # Ensure the path-coefficient file exists even when psem was not run.
  if (!file.exists(sem_psem_coef_file)) {
    readr::write_csv(
      tibble(Status = psem_status),
      sem_psem_coef_file
    )
  }

  sem_diagnostics <- c(
    sem_diagnostics,
    paste0("SEM status: ", psem_status),
    if (!is.null(psem_error)) paste0("SEM message/error: ", psem_error) else "SEM message/error: none",
    paste0("Component LM models fitted: ", length(sem_lm_models), "/", length(sem_formulas)),
    if (length(sem_lm_errors) > 0L) {
      paste0(
        "Component LM errors: ",
        paste(names(sem_lm_errors), sem_lm_errors, sep = " -> ", collapse = " | ")
      )
    } else {
      "Component LM errors: none"
    },
    paste0("Input data: ", sem_input_file),
    paste0("Component coefficients: ", sem_lm_coef_file),
    paste0("Standardized component coefficients: ", sem_lm_std_file),
    paste0("piecewiseSEM coefficients: ", sem_psem_coef_file)
  )
  writeLines(sem_diagnostics, sem_diag_file)

  # ALWAYS overwrite the summary placeholder with a useful final report.
  summary_header <- c(
    "============================================================",
    "PIECEWISE STRUCTURAL-EQUATION MODEL: BILE-ACID AXIS",
    "============================================================",
    paste0("Status: ", psem_status),
    paste0("Complete matched animals: ", nrow(sem_wide)),
    paste0(
      "Groups: ",
      if (nrow(sem_group_counts) > 0L) {
        paste(sem_group_counts$Group, sem_group_counts$n, sep = "=", collapse = ", ")
      } else "none"
    ),
    "",
    "Proposed biological chain (manual feature selection):",
    paste0(
      "  ",
      paste(
        paste0(
          sem_parts, "[",
          unname(SEM_FEATURE_BY_PART[sem_parts]),
          "]"
        ),
        collapse = " -> "
      )
    ),
    "",
    "Each component model also adjusts for experimental Group (NC/MCD/3_5).",
    ""
  )

  fallback_lines <- c(
    "------------------------------------------------------------",
    "COMPONENT PATH REGRESSIONS (always generated)",
    "------------------------------------------------------------",
    if (length(sem_lm_models) > 0L) {
      unlist(purrr::imap(sem_lm_models, function(model, nm) {
        c(
          paste0("\n### ", nm),
          capture.output(summary(model))
        )
      }), use.names = FALSE)
    } else {
      "No component path regression could be fitted."
    },
    "",
    paste0("Raw component coefficients CSV: ", sem_lm_coef_file),
    paste0("Standardized component coefficients CSV: ", sem_lm_std_file),
    ""
  )

  psem_lines <- if (length(sem_summary_lines) > 0L) {
    c(
      "------------------------------------------------------------",
      "piecewiseSEM::psem GLOBAL MODEL SUMMARY",
      "------------------------------------------------------------",
      sem_summary_lines,
      "",
      paste0("piecewiseSEM path coefficients CSV: ", sem_psem_coef_file),
      paste0("Saved model object: ", sem_model_file)
    )
  } else {
    c(
      "------------------------------------------------------------",
      "piecewiseSEM GLOBAL MODEL",
      "------------------------------------------------------------",
      "A piecewiseSEM global model was not available.",
      paste0("Reason/status: ", psem_status),
      if (!is.null(psem_error)) paste0("Message: ", psem_error) else "",
      "The component path regressions above remain valid regression outputs,",
      "but Fisher's C / directed-separation global-fit statistics require",
      "a successful piecewiseSEM::psem model.",
      "",
      paste0("Diagnostic file: ", sem_diag_file)
    )
  }

  writeLines(
    c(summary_header, fallback_lines, psem_lines),
    sem_summary_file
  )
}

# --------------------- 20. OPTIONAL TECHNICAL-REPLICATE COLLAPSE ---------------
# This is deliberately not applied automatically. After filling Replicate_Set in
# sample_metadata_map(3).csv, uncomment and run this block before the analyses above.
#
# data_collapsed <- data %>%
#   mutate(
#     Analysis_unit = case_when(
#       !is.na(Replicate_Set) & Replicate_Set != "" ~ Replicate_Set,
#       !is.na(Animal_ID) & Animal_ID != "" ~ paste(Animal_ID, Part, sep = "__"),
#       TRUE ~ Sample_ID
#     )
#   ) %>%
#   group_by(Analysis_unit, Group, Part, Animal_ID, Batch, Treatment_detail) %>%
#   summarise(
#     Sample_ID = paste(Sample_ID, collapse = ";"),
#     across(all_of(ba_columns), ~ mean(.x, na.rm = TRUE)),
#     .groups = "drop"
#   )
#
# Then replace `data` with `data_collapsed` and rerun from Section 3.

# ----------------------------- 21. SESSION INFO --------------------------------

# Explicitly stop cached SOCK workers before collecting session information.
stop_mixomics_backends()

capture.output(
  sessionInfo(),
  file = file.path(out_dir, "sessionInfo.txt")
)

message("Analysis complete. Results are in: ", normalizePath(out_dir))
