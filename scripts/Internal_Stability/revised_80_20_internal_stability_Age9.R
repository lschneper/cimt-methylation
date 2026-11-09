#!/usr/bin/env Rscript

# Repeated 80/20 internal stability analysis for primary EWAS findings
# Uses the same model as the original EWAS:
#   log_Mean_Mean_CIMT ~ CpG_M_value + covariates
#
#
# parallel is part of base R.

Sys.setenv(
  OMP_NUM_THREADS = 1,
  MKL_NUM_THREADS = 1,
  BLAS_NUM_THREADS = 1,
  LAPACK_NUM_THREADS = 1
)

library(parallel)
library(qs2)
OUTDIR="path to output directory"
# Create sig probes file
sighits<-read.delim("path to age9_specific_meta_analyzed_ewas_results")

write.table(sighits[sighits$bacon.fdr_fixed<0.05,],file=paste0(OUTDIR,"/y9_meta_sensitivity_model1_metalistdf2_output_chipseeker_sigprobes.txt"),sep="\t",quote=F,row.names=F)
setwd(OUTDIR)
# =============================================================================
# 1. USER SETTINGS: EDIT THIS SECTION FOR EACH AGE
# =============================================================================

AGE_LABEL <- "y9"
N_REPEATS <- 30L
TRAIN_FRAC <- 0.80
SEED <- 123L
NCORES <- 8L
TEST_P_CUTOFF <- 0.05

OUTPUT_FOLDER <- "path to output folder/80_20_internal_stability_revised/y9"
dir.create(OUTPUT_FOLDER, recursive = TRUE, showWarnings = FALSE)

SURVEY_FILE <- "path to file containing dependent variable and covariate data with array ID"
BETA_FILE <- "path to BetaMatrix.saliva.n2627.qcfilter.020425.qs"
MANIFEST_FILE <- "path to EPIC.hg19.manifest.dataframe.rds"


# File containing the CpGs that constitute the PRIMARY manuscript findings.
# For ages 9 and 15, use the significant EPIC+450K meta-analysis results.
# For age 22, use the significant primary EPIC results.
# The file may be an RDS data.frame or a tab-delimited text file.
PRIMARY_RESULTS_FILE <- paste0(OUTDIR,"y9_meta_sensitivity_model1_metalistdf2_output_chipseeker_sigprobes.txt")

# Column names in PRIMARY_RESULTS_FILE.
PRIMARY_PROBE_COLUMN <- "probe"
PRIMARY_BETA_COLUMN <- "TE_fixed"  # meta-analysis beta at ages 9/15; EPIC beta at age 22

OUTCOME_COLUMN <- "log_Mean_Mean_CIMT"
SAMPLE_ID_COLUMN <- "Array_ID"

# Use the exact covariates from the primary EWAS model being assessed.
COVARIATES <- c(
  "sexFemale", "Leukocytes", "SVPC1", "SVPC2", "SVPC3", "SVPC4", "SVPC5", "SVPC6", "SVPC7", "age"
)

# Exclude a primary CpG only if more than this fraction of its M-values are
# non-finite/missing among analysis participants.
MAX_CPG_MISSING_FRACTION <- 0.10

# =============================================================================
# 2. HELPER FUNCTIONS
# =============================================================================

read_results_file <- function(path) {
  if (!file.exists(path)) stop("PRIMARY_RESULTS_FILE does not exist: ", path)
  if (grepl("\\.rds$", path, ignore.case = TRUE)) {
    readRDS(path)
  } else {
    read.delim(path, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  }
}

# Fit the original regression separately for each CpG in one sample subset.
# model.matrix() handles numeric and factor covariates exactly as lm() does.
fit_cpg_models <- function(m_matrix, outcome, covariate_design, sample_indices,
                           probe_names, ncores = 1L) {
  y_sub <- outcome[sample_indices]
  m_sub <- m_matrix[sample_indices, , drop = FALSE]

  # The full-cohort covariate design is created once and then row-subset.
  # This preserves the same factor/dummy columns in every random split.
  cov_design <- covariate_design[sample_indices, , drop = FALSE]

  fit_one <- function(j) {
    cpg <- m_sub[, j]
    keep <- is.finite(cpg) & is.finite(y_sub) & complete.cases(cov_design)
    n_obs <- sum(keep)

    if (n_obs == 0L) {
      return(c(beta = NA_real_, se = NA_real_, t = NA_real_, p = NA_real_, n = 0,
               df = NA_real_))
    }

    X <- cbind(cov_design[keep, , drop = FALSE], CpG_M_value = cpg[keep])
    y <- y_sub[keep]

    fit <- lm.fit(x = X, y = y)
    rank <- fit$rank
    df_resid <- n_obs - rank

    # The CpG is deliberately the final design-matrix column.
    cpg_col <- ncol(X)
    beta <- unname(fit$coefficients[cpg_col])

    if (!is.finite(beta) || df_resid <= 0L) {
      return(c(beta = beta, se = NA_real_, t = NA_real_, p = NA_real_,
               n = n_obs, df = df_resid))
    }

    rss <- sum(fit$residuals^2)
    sigma2 <- rss / df_resid

    # vcov for estimable coefficients, respecting pivoting in lm.fit().
    R <- qr.R(fit$qr)[seq_len(rank), seq_len(rank), drop = FALSE]
    cov_unscaled_rank <- chol2inv(R)
    pivot_rank <- fit$qr$pivot[seq_len(rank)]

    cpg_rank_position <- match(cpg_col, pivot_rank)
    if (is.na(cpg_rank_position)) {
      return(c(beta = beta, se = NA_real_, t = NA_real_, p = NA_real_,
               n = n_obs, df = df_resid))
    }

    se <- unname(
      sqrt(sigma2 * cov_unscaled_rank[cpg_rank_position, cpg_rank_position])
    )
    t_value <- unname(beta / se)
    p_value <- unname(2 * pt(abs(t_value), df = df_resid, lower.tail = FALSE))

    return(c(beta = beta, se = se, t = t_value, p = p_value,
      n = n_obs, df = df_resid))
  }

  result_list <- if (.Platform$OS.type == "unix" && ncores > 1L) {
    mclapply(seq_along(probe_names), fit_one, mc.cores = ncores,
             mc.preschedule = TRUE)
  } else {
    lapply(seq_along(probe_names), fit_one)
  }

  result <- as.data.frame(do.call(rbind, result_list), stringsAsFactors = FALSE)
  result$probe <- probe_names
  result[, c("probe", "beta", "se", "t", "p", "n", "df")]
}

safe_rate <- function(x, denominator) {
  if (denominator == 0L) return(NA_real_)
  x / denominator
}

# =============================================================================
# 3. LOAD AND PREPARE THE SAME ANALYSIS DATA USED BY THE PRIMARY EWAS
# =============================================================================

set.seed(SEED)

survey <- readRDS(SURVEY_FILE)
if (!SAMPLE_ID_COLUMN %in% names(survey)) {
  stop("SAMPLE_ID_COLUMN was not found in the survey data: ", SAMPLE_ID_COLUMN)
}
rownames(survey) <- survey[[SAMPLE_ID_COLUMN]]

required_survey_columns <- c(OUTCOME_COLUMN, COVARIATES)
missing_survey_columns <- setdiff(required_survey_columns, names(survey))
if (length(missing_survey_columns) > 0L) {
  stop("Missing survey columns: ", paste(missing_survey_columns, collapse = ", "))
}

beta <- qs_read(BETA_FILE)

primary_results <- read_results_file(PRIMARY_RESULTS_FILE)
required_primary_columns <- c(PRIMARY_PROBE_COLUMN, PRIMARY_BETA_COLUMN)
missing_primary_columns <- setdiff(required_primary_columns, names(primary_results))
if (length(missing_primary_columns) > 0L) {
  stop("Missing primary-results columns: ",
       paste(missing_primary_columns, collapse = ", "))
}

primary_results <- primary_results[, required_primary_columns, drop = FALSE]
names(primary_results) <- c("probe", "primary_beta")
primary_results$probe <- as.character(primary_results$probe)
primary_results$primary_beta <- as.numeric(primary_results$primary_beta)
primary_results <- primary_results[!duplicated(primary_results$probe), , drop = FALSE]

# Keep participants with complete outcome and covariates and with EPIC data.
analysis_samples <- intersect(rownames(survey), colnames(beta))
survey_analysis <- survey[analysis_samples, c(OUTCOME_COLUMN, COVARIATES), drop = FALSE]
complete_subjects <- complete.cases(survey_analysis)
survey_analysis <- survey_analysis[complete_subjects, , drop = FALSE]
analysis_samples <- rownames(survey_analysis)

if (length(analysis_samples) < 20L) {
  stop("Too few complete participants after matching and covariate filtering.")
}

# Remove non-cg probes and sex/mitochondrial probes as in the EWAS.
manifest <- readRDS(MANIFEST_FILE)
manifest$ID <- rownames(manifest)
excluded_probes <- manifest[
  manifest$seqnames %in% c("chrX", "chrY", "chrM") |
    substr(manifest$ID, 1, 2) != "cg",
  "ID"
]

eligible_primary <- primary_results$probe[
  primary_results$probe %in% rownames(beta) &
    !primary_results$probe %in% excluded_probes
]

not_on_epic_or_excluded <- setdiff(primary_results$probe, eligible_primary)

if (length(eligible_primary) == 0L) {
  stop("None of the primary CpGs were found in the eligible EPIC probes.")
}

# CpGs x samples, then convert beta values to M-values.
beta_primary <- beta[eligible_primary, analysis_samples, drop = FALSE]
rm(beta)
gc()

m_primary <- log2(beta_primary) - log2(1 - beta_primary)
rm(beta_primary)
gc()

# Analysis matrix must be samples x CpGs.
m_matrix <- t(m_primary)
rm(m_primary)
gc()

# Preserve requested primary-results order.
m_matrix <- m_matrix[, eligible_primary, drop = FALSE]

missing_fraction <- colMeans(!is.finite(m_matrix))
keep_missing <- missing_fraction <= MAX_CPG_MISSING_FRACTION
removed_for_missingness <- names(missing_fraction)[!keep_missing]
m_matrix <- m_matrix[, keep_missing, drop = FALSE]
probe_names <- colnames(m_matrix)

primary_lookup <- primary_results[match(probe_names, primary_results$probe), , drop = FALSE]
stopifnot(identical(primary_lookup$probe, probe_names))

outcome <- as.numeric(survey_analysis[[OUTCOME_COLUMN]])
covariate_df <- survey_analysis[, COVARIATES, drop = FALSE]
covariate_design <- model.matrix(~ ., data = covariate_df)

# Save an audit table giving the exact denominator and exclusions.
probe_audit <- primary_results
probe_audit$found_on_EPIC <- probe_audit$probe %in% eligible_primary
probe_audit$excluded_manifest_or_absent <- probe_audit$probe %in% not_on_epic_or_excluded
probe_audit$missing_fraction <- missing_fraction[match(probe_audit$probe, names(missing_fraction))]
probe_audit$included_in_stability_analysis <- probe_audit$probe %in% probe_names

write.table(
  probe_audit,
  file = file.path(OUTPUT_FOLDER, paste0(AGE_LABEL, "_primary_probe_audit.txt")),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# =============================================================================
# 4. CREATE REPRODUCIBLE SPLITS AND RUN BOTH TRAINING AND TEST MODELS
# =============================================================================

n_samples <- nrow(m_matrix)
n_train <- floor(TRAIN_FRAC * n_samples)

split_assignments <- vector("list", N_REPEATS)
repeat_results <- vector("list", N_REPEATS)

for (r in seq_len(N_REPEATS)) {
  cat("Starting repeat", r, "of", N_REPEATS, "at", format(Sys.time()), "\n")

  train_idx <- sort(sample.int(n_samples, size = n_train, replace = FALSE))
  test_idx <- setdiff(seq_len(n_samples), train_idx)

  split_assignments[[r]] <- data.frame(
    repeat_id = r,
    sample_id = rownames(m_matrix),
    subset = ifelse(seq_len(n_samples) %in% train_idx, "train", "test"),
    stringsAsFactors = FALSE
  )

  train_fit <- fit_cpg_models(
    m_matrix = m_matrix,
    outcome = outcome,
    covariate_design = covariate_design,
    sample_indices = train_idx,
    probe_names = probe_names,
    ncores = NCORES
  )

  test_fit <- fit_cpg_models(
    m_matrix = m_matrix,
    outcome = outcome,
    covariate_design = covariate_design,
    sample_indices = test_idx,
    probe_names = probe_names,
    ncores = NCORES
  )

  stopifnot(identical(train_fit$probe, test_fit$probe))

  rr <- data.frame(
    repeat_id = r,
    probe = probe_names,
    primary_beta = primary_lookup$primary_beta,
    train_beta = train_fit$beta,
    train_se = train_fit$se,
    train_t = train_fit$t,
    train_p = train_fit$p,
    train_n = train_fit$n,
    test_beta = test_fit$beta,
    test_se = test_fit$se,
    test_t = test_fit$t,
    test_p = test_fit$p,
    test_n = test_fit$n,
    stringsAsFactors = FALSE
  )

  rr$train_same_direction_primary <-
    is.finite(rr$train_beta) & is.finite(rr$primary_beta) &
    sign(rr$train_beta) == sign(rr$primary_beta)

  rr$test_same_direction_primary <-
    is.finite(rr$test_beta) & is.finite(rr$primary_beta) &
    sign(rr$test_beta) == sign(rr$primary_beta)

  rr$test_same_direction_train <-
    is.finite(rr$test_beta) & is.finite(rr$train_beta) &
    sign(rr$test_beta) == sign(rr$train_beta)

  rr$train_nominal <- is.finite(rr$train_p) & rr$train_p < TEST_P_CUTOFF
  rr$test_nominal <- is.finite(rr$test_p) & rr$test_p < TEST_P_CUTOFF

  rr$test_nominal_same_primary <-
    rr$test_nominal & rr$test_same_direction_primary

  rr$test_nominal_same_train <-
    rr$test_nominal & rr$test_same_direction_train

  repeat_results[[r]] <- rr

  # Checkpoint after every repeat so completed work is retained if a job stops.
  saveRDS(
    repeat_results[seq_len(r)],
    file = file.path(
      OUTPUT_FOLDER,
      paste0(AGE_LABEL, "_repeat_results_checkpoint_through_", r, ".RDS")
    )
  )
}

repeat_results_df <- do.call(rbind, repeat_results)
split_assignments_df <- do.call(rbind, split_assignments)

# =============================================================================
# 5. SUMMARIZE USING A CONSISTENT DENOMINATOR OF 30 REPEATS
# =============================================================================

split_by_probe <- split(repeat_results_df, repeat_results_df$probe)

summary_list <- lapply(split_by_probe, function(d) {
  total_repeats <- N_REPEATS

  data.frame(
    probe = d$probe[1],
    primary_beta = d$primary_beta[1],
    repeats_expected = total_repeats,
    repeats_with_test_model = sum(is.finite(d$test_p)),
    train_nominal_n = sum(d$train_nominal, na.rm = TRUE),
    train_same_primary_n = sum(d$train_same_direction_primary, na.rm = TRUE),
    test_same_primary_n = sum(d$test_same_direction_primary, na.rm = TRUE),
    test_same_train_n = sum(d$test_same_direction_train, na.rm = TRUE),
    test_nominal_n = sum(d$test_nominal, na.rm = TRUE),
    test_nominal_same_primary_n = sum(d$test_nominal_same_primary, na.rm = TRUE),
    test_nominal_same_train_n = sum(d$test_nominal_same_train, na.rm = TRUE),
    train_nominal_fraction = safe_rate(sum(d$train_nominal, na.rm = TRUE), total_repeats),
    test_same_primary_fraction = safe_rate(sum(d$test_same_direction_primary, na.rm = TRUE), total_repeats),
    test_same_train_fraction = safe_rate(sum(d$test_same_direction_train, na.rm = TRUE), total_repeats),
    test_nominal_fraction = safe_rate(sum(d$test_nominal, na.rm = TRUE), total_repeats),
    test_nominal_same_primary_fraction = safe_rate(sum(d$test_nominal_same_primary, na.rm = TRUE), total_repeats),
    test_nominal_same_train_fraction = safe_rate(sum(d$test_nominal_same_train, na.rm = TRUE), total_repeats),
    median_train_beta = median(d$train_beta, na.rm = TRUE),
    median_test_beta = median(d$test_beta, na.rm = TRUE),
    median_train_p = median(d$train_p, na.rm = TRUE),
    median_test_p = median(d$test_p, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
})

summary_df <- do.call(rbind, summary_list)
rownames(summary_df) <- NULL
summary_df <- summary_df[match(probe_names, summary_df$probe), , drop = FALSE]

# Optional descriptive thresholds. These are not labeled as validated/failed.
summary_df$test_same_primary_at_least_24_of_30 <- summary_df$test_same_primary_n >= 24L
summary_df$test_nominal_same_primary_at_least_1_of_30 <- summary_df$test_nominal_same_primary_n >= 1L
summary_df$test_nominal_same_primary_at_least_15_of_30 <- summary_df$test_nominal_same_primary_n >= 15L
summary_df$test_nominal_same_primary_at_least_21_of_30 <- summary_df$test_nominal_same_primary_n >= 21L
summary_df$test_nominal_same_primary_at_least_24_of_30 <- summary_df$test_nominal_same_primary_n >= 24L

# Age-level counts for easy manuscript/rebuttal reporting.
age_summary <- data.frame(
  age = AGE_LABEL,
  primary_CpGs_input = nrow(primary_results),
  primary_CpGs_evaluable = nrow(summary_df),
  excluded_absent_or_manifest = length(not_on_epic_or_excluded),
  excluded_missing_over_threshold = length(removed_for_missingness),
  same_primary_direction_at_least_24_of_30_n =
    sum(summary_df$test_same_primary_at_least_24_of_30, na.rm = TRUE),
  same_primary_direction_at_least_24_of_30_percent =
    100 * mean(summary_df$test_same_primary_at_least_24_of_30, na.rm = TRUE),
  nominal_same_primary_at_least_1_of_30_n =
    sum(summary_df$test_nominal_same_primary_at_least_1_of_30, na.rm = TRUE),
  nominal_same_primary_at_least_1_of_30_percent =
    100 * mean(summary_df$test_nominal_same_primary_at_least_1_of_30, na.rm = TRUE),
  nominal_same_primary_at_least_15_of_30_n =
    sum(summary_df$test_nominal_same_primary_at_least_15_of_30, na.rm = TRUE),
  nominal_same_primary_at_least_15_of_30_percent =
    100 * mean(summary_df$test_nominal_same_primary_at_least_15_of_30, na.rm = TRUE),
  nominal_same_primary_at_least_21_of_30_n =
    sum(summary_df$test_nominal_same_primary_at_least_21_of_30, na.rm = TRUE),
  nominal_same_primary_at_least_21_of_30_percent =
    100 * mean(summary_df$test_nominal_same_primary_at_least_21_of_30, na.rm = TRUE),
  nominal_same_primary_at_least_24_of_30_n =
    sum(summary_df$test_nominal_same_primary_at_least_24_of_30, na.rm = TRUE),
  nominal_same_primary_at_least_24_of_30_percent =
    100 * mean(summary_df$test_nominal_same_primary_at_least_24_of_30, na.rm = TRUE),
  stringsAsFactors = FALSE
)

# =============================================================================
# 6. SAVE ALL MATERIAL NEEDED FOR THE MANUSCRIPT, SUPPLEMENT, AND AUDIT TRAIL
# =============================================================================

saveRDS(repeat_results_df,
        file.path(OUTPUT_FOLDER, paste0(AGE_LABEL, "_80_20_repeat_level_results.RDS")))
saveRDS(summary_df,
        file.path(OUTPUT_FOLDER, paste0(AGE_LABEL, "_80_20_CpG_summary.RDS")))
saveRDS(split_assignments_df,
        file.path(OUTPUT_FOLDER, paste0(AGE_LABEL, "_80_20_split_assignments.RDS")))

write.table(
  repeat_results_df,
  file = file.path(OUTPUT_FOLDER, paste0(AGE_LABEL, "_80_20_repeat_level_results.txt")),
  sep = "\t", quote = FALSE, row.names = FALSE
)

write.table(
  summary_df,
  file = file.path(OUTPUT_FOLDER, paste0(AGE_LABEL, "_80_20_CpG_summary.txt")),
  sep = "\t", quote = FALSE, row.names = FALSE
)

write.table(
  age_summary,
  file = file.path(OUTPUT_FOLDER, paste0(AGE_LABEL, "_80_20_age_summary.txt")),
  sep = "\t", quote = FALSE, row.names = FALSE
)

write.table(
  split_assignments_df,
  file = file.path(OUTPUT_FOLDER, paste0(AGE_LABEL, "_80_20_split_assignments.txt")),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(as.data.frame(missing_fraction[removed_for_missingness]),
            file = file.path(OUTPUT_FOLDER, paste0(AGE_LABEL, "_80_20_removedCpGs_missing_fraction.txt")),
            sep = "\t", quote = FALSE, row.names = TRUE)
            
analysis_settings <- list(
  age_label = AGE_LABEL,
  n_repeats = N_REPEATS,
  train_fraction = TRAIN_FRAC,
  seed = SEED,
  ncores = NCORES,
  nominal_p_cutoff = TEST_P_CUTOFF,
  outcome = OUTCOME_COLUMN,
  covariates = COVARIATES,
  max_CpG_missing_fraction = MAX_CPG_MISSING_FRACTION,
  survey_file = SURVEY_FILE,
  beta_file = BETA_FILE,
  manifest_file = MANIFEST_FILE,
  primary_results_file = PRIMARY_RESULTS_FILE,
  primary_probe_column = PRIMARY_PROBE_COLUMN,
  primary_beta_column = PRIMARY_BETA_COLUMN,
  analysis_sample_n = n_samples,
  primary_CpG_n_input = nrow(primary_results),
  primary_CpG_n_evaluable = length(probe_names),
  session_info = sessionInfo()
)

saveRDS(analysis_settings,
        file.path(OUTPUT_FOLDER, paste0(AGE_LABEL, "_80_20_analysis_settings_and_sessionInfo.RDS")))

capture.output(
  print(analysis_settings),
  file = file.path(OUTPUT_FOLDER, paste0(AGE_LABEL, "_80_20_analysis_settings_and_sessionInfo.txt"))
)

cat("\nCompleted analysis.\n")
print(age_summary)

sig_counts<-tapply(repeat_results_df$test_nominal_same_primary,repeat_results_df$probe,sum)
summary_stats<-c(
  Age=AGE_LABEL,
  Median=summary(sig_counts)[3],
  Mean=summary(sig_counts)[4],
  Q1=summary(sig_counts)[2],
  Q3=summary(sig_counts)[5],
  Min=summary(sig_counts)[1],
  Max=summary(sig_counts)[6]
)

write.table(
  as.data.frame(summary_stats),
  file = paste0(OUTPUT_FOLDER, "/",AGE_LABEL, "_80_20_summary_stats.txt"),
  sep = "\t", quote = FALSE, row.names = TRUE
)
