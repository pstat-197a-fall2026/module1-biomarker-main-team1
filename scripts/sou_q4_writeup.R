# Sou's report generator. Called by sou_q.4.R after the analysis finishes.
# To inspect results later, use load("results/sou-q4-final/sou-q4-results.RData").

pct <- function(x) sprintf("%.1f%%", 100 * x)
num <- function(x) sprintf("%.3f", x)
md_table <- function(x) {
  x <- as.data.frame(x)
  x[] <- lapply(x, function(v) ifelse(is.na(v), "Not recorded", as.character(v)))
  c(paste0("| ", paste(names(x), collapse = " | "), " |"),
    paste0("| ", paste(rep("---", ncol(x)), collapse = " | "), " |"),
    apply(x, 1, function(row) paste0("| ", paste(row, collapse = " | "), " |")))
}
html_escape <- function(x) as.character(htmltools::htmlEscape(as.character(x)))
html_table <- function(x) {
  x <- as.data.frame(x)
  x[] <- lapply(x, function(v) ifelse(is.na(v), "Not recorded", as.character(v)))
  paste0("<div class='table-wrap'><table><thead><tr>",
         paste0("<th>", html_escape(names(x)), "</th>", collapse = ""),
         "</tr></thead><tbody>",
         paste(apply(x, 1, function(row) paste0("<tr>",
           paste0("<td>", html_escape(row), "</td>", collapse = ""), "</tr>")),
           collapse = ""), "</tbody></table></div>")
}

cv_display <- cv_summary %>%
  transmute(Approach = model, `Mean CV accuracy` = pct(mean_accuracy),
            `Mean CV AUROC` = num(mean_auc),
            `Mean panel size` = sprintf("%.1f", mean_panel_size))
test_display <- benchmark_comparison %>%
  transmute(Approach = model, Proteins = n_proteins,
            Sensitivity = pct(sensitivity), Specificity = pct(specificity),
            Accuracy = pct(accuracy), AUROC = num(roc_auc))
count_display <- confusion_counts %>%
  transmute(Approach = role, `ASD correctly identified` = true_ASD,
            `ASD missed` = missed_ASD, `TD correctly identified` = true_TD,
            `TD called ASD` = false_ASD, `Total correct` = paste0(correct, "/", total))
primary_cv <- cv_summary %>% filter(model == chosen_model)
simple_cv <- cv_summary %>% filter(model == simple_model)
reference_cv <- cv_summary %>% filter(model == "hard10")
baseline_accuracy <- baseline_fit$metrics$accuracy
primary_accuracy <- sou_fit$metrics$accuracy
simple_accuracy <- simple_fit$metrics$accuracy
simple_matches <- length(simple_panel) < length(inclass_panel) &&
  simple_accuracy >= baseline_accuracy - accuracy_margin
primary_improves <- primary_accuracy > baseline_accuracy
simple_conclusion <- if (simple_matches) {
  paste0("The simpler panel meets our descriptive two-percentage-point accuracy ",
         "criterion against the rerun class baseline on this partition, using fewer proteins. ",
         "This is an observed comparison, not proof of statistical equivalence.")
} else {
  paste0("The simpler panel does not meet the descriptive two-percentage-point accuracy ",
         "criterion against the rerun class baseline on this partition. ",
         "It should not be presented as a demonstrated equivalent replacement.")
}
primary_conclusion <- if (primary_improves) {
  "The primary alternative has higher observed test accuracy than the rerun class baseline; the small reused test set cannot establish a reliable generalization gain."
} else if (abs(primary_accuracy - baseline_accuracy) < 1e-12) {
  "The primary alternative matches the rerun class baseline's observed accuracy; it does not improve it."
} else {
  "The primary alternative has lower observed accuracy than the rerun class baseline."
}
verdict <- paste(simple_conclusion, primary_conclusion)

sections <- list()
add_text <- function(title, paragraphs) {
  sections[[length(sections) + 1L]] <<- list(type = "text", title = title, body = paragraphs)
}
add_table <- function(title, data) {
  sections[[length(sections) + 1L]] <<- list(type = "table", title = title, body = data)
}

add_text("Question 4: objective and conclusion", c(
  "The assignment asks for either a simpler protein panel with comparable classification accuracy or an alternative panel with improved accuracy, benchmarked against the in-class analysis. A panel is the set of protein measurements used by the final classifier.",
  verdict))
add_text("How Question 3 informed this analysis", c(
  "Q3a showed why the test participants must not influence feature selection. Here, the original 80/20 participant partition and seed are retained. All new-model scaling and protein selection are fitted on training participants, including within every validation fold. No participant is deleted because of outlying measurements.",
  "Q3b's recorded results did not show an accuracy increase from taking the top 25 instead of top 10 proteins per method. This motivates comparing both candidate sets rather than assuming more proteins will help. Lasso is added to penalize large coefficients and set some to zero.",
  "Q3c motivates softer combination of t-test and random-forest rankings. The fuzzy score is the minimum of two rank-based membership scores. Unlike Q3c's 0.95 cutoff, this analysis uses a fixed number of top fuzzy-ranked proteins to control panel size. It also tests each ranking alone for small panels."))
add_text("Methods and choices fixed before this extended run", c(
  sprintf("The raw CSV contains %d participants and %d protein measurement columns. ADOS and subject ID are excluded as predictors. The seed is %d; %d participants train the models and %d form the reused test partition.",
          nrow(raw_data), length(protein_names), seed, nrow(train_data), nrow(test_data)),
  "Protein abundance is transformed with log10. Each protein's mean and standard deviation are estimated from the current training partition. Those values standardize both training and evaluation measurements; standardized values are clipped at -3 and +3. This uses a fixed clipping rule, not learned outlier removal.",
  "Within the 123 training participants, stratified five-fold cross-validation is repeated three times. The 24 fixed configurations are: top-10 and top-25 hard intersections; fuzzy top 1, 2, 3 and 5; t-test top 1, 2 and 3; random-forest top 1, 2 and 3; and six lasso penalty fractions (0.5, 0.2, 0.1, 0.05, 0.02, 0.01) for each of two pools (fuzzy top 25 and all proteins). This yields 360 validation fits. There is no threshold search: probability greater than 0.5 always predicts ASD.",
  "Welch t-tests and the group's existing adjusted-score formula rank proteins. That formula is retained for benchmark consistency; its numerical scores are not used to claim significance. Random forest uses 1,000 trees and MeanDecreaseGini importance. Ordinary logistic regression combines selected measurements; lasso logistic regression instead shrinks coefficients and selects predictors jointly. Lasso penalties are fractions of the largest penalty on each training fit's path; coefficients are not refitted with unpenalized logistic regression.",
  sprintf("The primary approach maximizes mean validation accuracy, with mean panel size and AUROC as tie-breakers. The simpler approach is restricted to fixed panels of 1-3 proteins; among models within two percentage points of the best small-panel validation accuracy, the smallest is chosen. These choices are saved before computing final test predictions. Coefficients of absolute magnitude at most %g are numerical zero for panel counting.", coefficient_tolerance),
  "These are exploratory extensions after earlier test scores were seen. The search is not expanded after this run to chase a test improvement. The CV maximum is a tuning score, not an unbiased estimate of the selected procedure's accuracy."))
add_table("Training cross-validation: all approaches", cv_display)
add_text("Locked choices", c(
  sprintf("The primary approach is %s: mean validation accuracy %s, mean AUROC %s and mean panel size %.1f. Its final training fit uses %d proteins.",
          chosen_model, pct(primary_cv$mean_accuracy), num(primary_cv$mean_auc),
          primary_cv$mean_panel_size, length(sou_panel)),
  sprintf("The simpler approach is %s: mean validation accuracy %s, mean AUROC %s and final panel size %d. The training-only hard-intersection reference has mean validation accuracy %s. Protein identities may vary across folds; protein-selection-frequency.csv records how often each was selected.",
          simple_model, pct(simple_cv$mean_accuracy), num(simple_cv$mean_auc),
          length(simple_panel), pct(reference_cv$mean_accuracy))))
add_table("Exploratory test results and class benchmarks", test_display)
add_table("What the classifications mean in children", count_display)
add_text("Protein panels", c(
  paste0("Rerun in-class panel (", length(inclass_panel), "): ", paste(inclass_panel, collapse = ", "), "."),
  paste0("Training-only reference panel (", length(reference_panel), "): ", paste(reference_panel, collapse = ", "), "."),
  paste0("Simpler panel (", length(simple_panel), "): ", paste(simple_panel, collapse = ", "), "."),
  paste0("Primary panel (", length(sou_panel), "): ", paste(sou_panel, collapse = ", "), ".")))
add_text("Benchmark interpretation", c(
  "The in-class procedure is reproduced in Sou's script using global preprocessing and selection, retaining those limitations solely for the requested benchmark. Teammates' scripts are neither sourced nor edited. The rerun uses the same outer participant partition as the new methods. The training-only reference is a distinct benchmark and is labeled separately.",
  sprintf("The original Q3 baseline comments record accuracy 83.9%% and AUROC 0.908. The reproduced class procedure in this environment gives accuracy %s and AUROC %s. The cause of the discrepancy is unresolved. A result must not be described as exceeding the historical 83.9%% unless it actually does.",
          pct(baseline_accuracy), num(baseline_fit$metrics$roc_auc)),
  sprintf("Against the rerun class baseline, the primary alternative changes accuracy by %+.1f percentage points and the simpler alternative by %+.1f points. Sensitivity and specificity must also be considered; a higher overall accuracy or AUROC does not imply more ASD participants are correctly identified at the 0.5 threshold.",
          100 * (primary_accuracy - baseline_accuracy), 100 * (simple_accuracy - baseline_accuracy))))
add_text("Limitations and final interpretation", c(
  verdict,
  "There are only 31 test participants, so one correct classification changes accuracy by about 3.2 percentage points. The same test subjects were viewed in earlier work; these results are exploratory. Repeated training CV supports model development but model selection can make its best score optimistic. A new independent sample or a fully nested evaluation would be needed for stronger claims. Repeated folds share participants and are not independent replicates.",
  "Protein selection frequency describes stability within these training splits. It does not demonstrate biological causation or establish a clinical ASD test. This report evaluates a statistical assignment, not a validated diagnostic tool."))
add_text("How to run and inspect the work", c(
  'Open the team RStudio project, restart R, and Source scripts/sou_q.4.R. Run no Q1-Q3 scripts first. The only helper sourced is Sou\'s own scripts/sou_q4_writeup.R. Wait for all 15 folds. This creates results/sou-q4-final/ and the named objects in the RStudio Environment.',
  'Inspect View(cv_summary), View(benchmark_comparison), simple_panel, sou_panel, and View(confusion_counts). To inspect without recomputing, load("results/sou-q4-final/sou-q4-results.RData"). A separate terminal run does not populate the RStudio Environment.',
  "The earlier scripts/drafts/sou_q4_starter.R is historical scratch work. Keep it in drafts as requested by the assignment; it is not an input or prerequisite. Older results/sou-q4/ and results/run-output/ belong to earlier attempts and are not the current findings."))
add_text("Files to share with the team", c(
  "Share scripts/sou_q.4.R, scripts/sou_q4_writeup.R, scripts/sou_q4_check.R, scripts/sou_q4_plan.md, scripts/drafts/sou_q4_starter.R, and results/sou-q4-final/. The current report is sou-q4-report.html; sou-q4-report.Rmd and sou-q4-report.md contain the same generated text, and sou-q4-team-section.qmd is a short contribution for the team's report. Source scripts/sou_q4_check.R to check saved results without rerunning training. No existing group report or analysis script is overwritten.",
  "The CSVs include all validation scores, selected proteins, training selection frequencies, locked choices, benchmark metrics, confusion counts, final coefficients, subject partitions and test probabilities. sou-q4-results.RData stores fitted model objects and results; session-info.txt and input-fingerprint.csv record the software and raw input used. The team still needs to integrate its final analysis into analysis-main.R and its shared report before submission."))

md <- c("# Question 4 — protein-panel comparison", "", "Author: Sou", "")
html_body <- c("<h1>Question 4 — protein-panel comparison</h1>", "<p class='byline'>Sou · exploratory analysis</p>")
for (section in sections) {
  md <- c(md, paste0("## ", section$title), "")
  html_body <- c(html_body, paste0("<h2>", html_escape(section$title), "</h2>"))
  if (section$type == "text") {
    for (paragraph in section$body) {
      md <- c(md, paragraph, "")
      html_body <- c(html_body, paste0("<p>", html_escape(paragraph), "</p>"))
    }
  } else {
    md <- c(md, md_table(section$body), "")
    html_body <- c(html_body, html_table(section$body))
  }
  if (section$title == "Training cross-validation: all approaches") {
    md <- c(md, "![Cross-validation accuracy and panel size](cv-comparison.png)", "")
    html_body <- c(html_body, "<img src='cv-comparison.png' alt='Cross-validation accuracy and panel size'>")
  }
}
writeLines(md, file.path(output_dir, "sou-q4-report.md"))
# A short standalone contribution for the shared Quarto report. No shared file
# is changed; the team can integrate this section during its final assembly.
team_section <- c(
  "### Question 4: alternative and simpler protein panels", "",
  "Motivated by Q3, we split participants before preprocessing and feature selection, compared both hard and fuzzy rankings, and investigated lasso logistic regression. We compared 24 fixed configurations using three repeats of stratified five-fold cross-validation within the 123 training participants. Scaling and selection were repeated within each fold. A primary approach maximized mean validation accuracy; a second approach favored a one- to three-protein panel within two percentage points of the best small-panel validation score. All classifiers used a fixed 0.5 threshold.", "",
  sprintf("Cross-validation selected %s as the primary approach (mean accuracy %s) and %s as the simpler approach (mean accuracy %s). The final primary panel contained %d proteins; the simpler panel contained %s. The original class procedure was independently reproduced within Sou's new script; no teammate scripts were edited or required as inputs.",
          chosen_model, pct(primary_cv$mean_accuracy), simple_model,
          pct(simple_cv$mean_accuracy), length(sou_panel),
          paste(simple_panel, collapse = ", ")), "",
  md_table(test_display), "",
  verdict, "",
  sprintf("The primary alternative's sensitivity was %s, compared with %s for the rerun class baseline; therefore overall accuracy and ASD identification must be interpreted separately. The historical baseline comments report 83.9%% accuracy, but the same procedure in the current environment produced %s; the cause of this discrepancy remains unresolved.",
          pct(sou_fit$metrics$sensitivity), pct(baseline_fit$metrics$sensitivity),
          pct(baseline_accuracy)), "",
  "The 31 test participants were examined in earlier work, so the final comparison is exploratory. One child changes test accuracy by about 3.2 percentage points. CV scores were used for selection and may be optimistic. Comparable observed accuracy is not proof of equivalence, and these results do not establish a clinical biomarker panel. Full protein lists, coefficients, confusion counts and reproducible code accompany this section.", "")
writeLines(team_section, file.path(output_dir, "sou-q4-team-section.qmd"))
# Static R Markdown is deliberate: the analysis script creates authoritative
# results once; rendering this report does not secretly rerun model selection.
writeLines(c("---", 'title: "Question 4: protein-panel comparison"',
             'author: "Sou"', "output: html_document", "---", "", md[-(1:4)]),
           file.path(output_dir, "sou-q4-report.Rmd"))
css <- "body{font:16px/1.65 system-ui,sans-serif;color:#233343;max-width:1100px;margin:48px auto;padding:0 24px}h1{font-size:32px;line-height:1.2;color:#174768}h2{margin-top:38px;font-size:22px;color:#174768}.byline{color:#63778a}table{border-collapse:collapse;width:100%;font-size:14px}th,td{padding:9px 12px;text-align:left;border-bottom:1px solid #dce4ea}th{background:#edf3f7}.table-wrap{overflow-x:auto}img{max-width:100%;height:auto}@media print{body{margin:0;font-size:11pt}h2{break-after:avoid}tr{break-inside:avoid}}"
writeLines(c("<!doctype html><html lang='en'><head><meta charset='utf-8'>",
             "<meta name='viewport' content='width=device-width,initial-scale=1'>",
             "<title>Sou — Question 4 findings</title>", paste0("<style>", css, "</style>"),
             "</head><body>", html_body, "</body></html>"),
           file.path(output_dir, "sou-q4-report.html"))
writeLines(c("# Current Q4 results", "", verdict, "",
             "Open sou-q4-report.html for the full findings and limitations.",
             "Run scripts/sou_q.4.R from the project root to regenerate all results.",
             "Keep the plot next to the HTML report when sharing it.",
             "The .Rmd and .md report sources can be copied into the team's report.",
             "Do not substitute the older results/sou-q4/ tables for these results."),
           file.path(output_dir, "README.md"))
