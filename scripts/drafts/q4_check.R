# Verify saved Q4 results without rerunning cross-validation or team scripts.
# Run AFTER q4_analysis.R, or Source this file to check existing saved results.
stopifnot(requireNamespace("glmnet", quietly = TRUE),
          requireNamespace("rsample", quietly = TRUE))
project_root <- local({
  candidates <- c(".", "..", "../..")
  matches <- candidates[file.exists(file.path(candidates, "data", "biomarker-clean.RData"))]
  if (!length(matches)) stop("Open the team project before running this script.")
  normalizePath(matches[[1]], mustWork = TRUE)
})
result_dir <- file.path(project_root, "results", "sou-q4-final")
audit <- new.env(parent = globalenv())
load(file.path(result_dir, "sou-q4-results.RData"), envir = audit)
with(audit, {
  stopifnot(nrow(configurations) == 24L, nrow(folds) == 15L,
            nrow(cv_results) == 360L, !anyNA(cv_results),
            nrow(final_results) == 3L,
            !anyDuplicated(analysis_data$subject_id),
            !length(intersect(train_data$subject_id, test_data$subject_id)),
            length(simple_panel) <= 3,
            !any(c("ados", "subject_id", "group") %in% protein_names))
  for (i in seq_len(nrow(folds))) {
    a <- rsample::analysis(folds$splits[[i]])$subject_id
    b <- rsample::assessment(folds$splits[[i]])$subject_id
    stopifnot(!length(intersect(a, b)),
              all(c(a, b) %in% train_data$subject_id),
              !length(intersect(c(a, b), test_data$subject_id)))
  }
  for (result in list(baseline_fit, alternative_fit, simple_fit)) {
    stopifnot(length(result$probabilities) == nrow(test_data))
    if (result$config$method == "glm") {
      regenerated <- as.numeric(predict(result$fit,
        newdata = as.data.frame(final_data$x_evaluation), type = "response"))
    } else {
      panel <- if (result$config$pool == "all") colnames(final_data$x_train) else
        head(final_rankings$fuzzy, 25)
      regenerated <- as.numeric(predict(result$fit,
        newx = final_data$x_evaluation[, panel, drop = FALSE],
        s = result$lambda, type = "response"))
    }
    stopifnot(isTRUE(all.equal(regenerated, result$probabilities, tolerance = 1e-12)))
  }
  probabilities <- read.csv(file.path(result_dir, "exploratory-test-predictions.csv"))
  for (role in c("inclass", "primary", "simpler")) {
    column <- switch(role, inclass = "inclass_probability",
      primary = "alternative_probability",
      simpler = "simpler_probability")
    p <- probabilities[[column]]
    truth <- probabilities$group == "ASD"
    selected <- p > 0.5
    count <- confusion_counts[confusion_counts$role == role, ]
    stopifnot(all(is.finite(p)), all(p >= 0 & p <= 1),
      sum(selected == truth) == count$correct,
      sum(selected & truth) == count$true_ASD,
      sum(!selected & truth) == count$missed_ASD,
      sum(!selected & !truth) == count$true_TD,
      sum(selected & !truth) == count$false_ASD)
  }
  stopifnot(length(alternative_panel) == alternative_fit$metrics$n_proteins,
            length(simple_panel) == simple_fit$metrics$n_proteins)
})
cat("PASS: all 360 validation fits, disjoint partitions, saved model predictions,\n",
    "exported classification counts, three-row baseline table, and panel sizes checked.\n", sep = "")
