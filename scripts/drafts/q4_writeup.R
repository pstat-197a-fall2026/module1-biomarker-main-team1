# Q4 report generator, called by analysis-main.R or q4_analysis.R.
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

# One current-run comparison table. The historical numbers are a separate note.
cv_display <- cv_summary %>%
  transmute(Approach = model, `Mean CV accuracy` = pct(mean_accuracy),
            `Mean CV AUROC` = num(mean_auc), `Mean panel size` = sprintf("%.1f", mean_panel_size))
test_display <- final_results %>%
  transmute(Approach = model, Proteins = n_proteins, Sensitivity = pct(sensitivity),
            Specificity = pct(specificity), Accuracy = pct(accuracy), AUROC = num(roc_auc))
counts_display <- confusion_counts %>%
  transmute(Approach = role, `Correct ASD` = true_ASD, `Missed ASD` = missed_ASD,
            `Correct TD` = true_TD, `TD called ASD` = false_ASD,
            `Total correct` = paste0(correct, "/", total))
primary_cv <- cv_summary %>% filter(model == chosen_model)
simple_cv <- cv_summary %>% filter(model == simple_model)
baseline_accuracy <- baseline_fit$metrics$accuracy
primary_accuracy <- alternative_fit$metrics$accuracy
simple_accuracy <- simple_fit$metrics$accuracy
simple_matches <- length(simple_panel) < length(inclass_panel) &&
  simple_accuracy >= baseline_accuracy - accuracy_margin
primary_conclusion <- if (primary_accuracy > baseline_accuracy + 1e-12) {
  "The alternative has higher observed accuracy than the rerun class baseline on this partition, but the reused test set does not establish a reliable improvement."
} else if (abs(primary_accuracy - baseline_accuracy) < 1e-12) {
  "The alternative matches the rerun class baseline's observed accuracy; it does not improve accuracy."
} else {
  "The alternative has lower observed accuracy than the rerun class baseline."
}
simple_conclusion <- if (simple_matches) {
  "The smaller panel meets the descriptive two-percentage-point comparability criterion against the rerun class baseline, using fewer proteins. This is not proof of statistical equivalence."
} else {
  "The smaller panel does not meet the descriptive two-percentage-point comparability criterion against the rerun class baseline."
}
verdict <- paste(primary_conclusion, simple_conclusion)
methods <- c(
  "Question 4 investigates an alternative panel to improve classification accuracy, with a smaller panel as a secondary option. A panel is the set of protein measurements used by the classifier.",
  "To align with the team, all comparisons use the supplied biomarker-clean.RData directly: the same log-transformed, globally standardized and clipped values used in Q3 and the original class model. The random forest receives factor labels ASD and TD in that order. The seed is 101422, the train/test split is 80/20, ASD is the positive logistic outcome, and predicted probability greater than 0.5 means ASD. ADOS and subject IDs are excluded as predictors.",
  "The class baseline uses the top-10 t-test and random-forest lists selected on the full cleaned dataset, their hard intersection, and logistic regression fitted on the training participants. That original full-data selection is retained only for the requested class benchmark. Q4 selects proteins using training participants and uses cross-validation there to choose its approach. No teammate script is sourced or changed.",
  "Q3a motivates selecting proteins on training data. Q3b motivates comparing broader lists without assuming they improve accuracy. Q3c motivates rank-based fuzzy selection, with fixed top-k panels rather than its 0.95 cutoff. Lasso logistic regression adds a penalty that can shrink protein coefficients to zero.",
  sprintf("The existing 24-configuration search is unchanged: hard top-10/top-25 intersections; fuzzy top 1/2/3/5; t-test top 1/2/3; forest top 1/2/3; and six penalty fractions for lasso using fuzzy top 25 or all proteins. Three repeats of stratified five-fold CV within %d training participants yield 360 fits. The primary choice maximizes mean CV accuracy. The simpler choice is restricted to 1-3 proteins and prefers the smallest panel within two percentage points of the best small-panel CV accuracy. Choices are saved before test evaluation. This is a tuning rule, not an equivalence test.", nrow(train_data)),
  "Unlike the previous Q4 revision, this analysis does not recalculate preprocessing on training subjects or recode forest labels to FALSE/TRUE. It has one class benchmark in the final table. The correction changes training-CV results and can change selected panels; prior Q4 tables should not be mixed with this revision. The test set was seen before this correction, so all final test results remain exploratory.")
results_text <- c(
  sprintf("Training cross-validation selects %s as the alternative (mean accuracy %s) and %s as the smaller approach (mean accuracy %s). Fitting those rules to all training subjects produces %d and %d proteins respectively.", chosen_model, pct(primary_cv$mean_accuracy), simple_model, pct(simple_cv$mean_accuracy), length(alternative_panel), length(simple_panel)),
  verdict,
  sprintf("The alternative's sensitivity is %s versus %s for the class baseline. Overall accuracy, ASD detection and AUROC therefore need separate interpretation; a higher AUROC is not automatically higher accuracy at the fixed 0.5 threshold.", pct(alternative_fit$metrics$sensitivity), pct(baseline_fit$metrics$sensitivity)))
baseline_note <- sprintf("The team's q3_base.R comments record accuracy 83.9%% and AUROC 0.908. The class baseline reproduced in this run gives accuracy %s and AUROC %s. The historical numbers are not inserted into the current-run table as if they came from this run. The discrepancy cannot be attributed to the metric syntax alone: direct counting and yardstick agree on the current predictions. To resolve differences across computers, compare selected proteins, test participant IDs, predicted probabilities, input-file hash and sessionInfo().", pct(baseline_accuracy), num(baseline_fit$metrics$roc_auc))
limitations <- c(
  "The supplied clean data were standardized using all subjects. Reusing them aligns the comparison with the team, but preprocessing is not fully isolated from validation/test participants. The original class benchmark also selects features before splitting. These limitations are retained and disclosed, rather than describing this as a fully leakage-free pipeline.",
  sprintf("There are only %d test participants, so one classification changes accuracy by about %.1f percentage points. Test subjects were viewed in earlier work. CV scores were used to select models and their maximum may be optimistic; repeated folds share subjects and are not independent replications. No seed, panel or threshold was chosen to maximize the final test score.",nrow(test_data),100/nrow(test_data)),
  "The two-percentage-point margin is descriptive. Neither comparable observed accuracy nor a small improvement establishes statistical equivalence or a dependable generalization gain. These exploratory results do not establish biological causation or a clinical ASD test.")
sections <- list(
  list(title="Objective and alignment",type="text",body=methods),
  list(title="Current-run comparison with the class baseline",type="table",body=test_display),
  list(title="Findings",type="text",body=results_text),
  list(title="Classifications in individual children",type="table",body=counts_display),
  list(title="Historical baseline discrepancy",type="text",body=baseline_note),
  list(title="Selected proteins",type="text",body=c(
    paste0("Class baseline (",length(inclass_panel),"): ",paste(inclass_panel,collapse=", "),"."),
    paste0("Smaller panel (",length(simple_panel),"): ",paste(simple_panel,collapse=", "),"."),
    paste0("Alternative panel (",length(alternative_panel),"): ",paste(alternative_panel,collapse=", "),"."))),
  list(title="All training cross-validation results",type="table",body=cv_display),
  list(title="Limitations",type="text",body=limitations),
  list(title="Running and sharing this work",type="text",body=c(
    'Open the team project and Source scripts/drafts/q4_analysis.R. It reads data/biomarker-clean.RData and calls the Q4 report helper scripts/drafts/q4_writeup.R. Do not run teammate scripts first. Wait for all 15 folds, then inspect cv_summary, final_results, inclass_panel, alternative_panel and simple_panel.',
    'Use sou-q4-team-section.qmd as the contribution to the shared report. Keep q4_analysis.R and q4_writeup.R together in scripts/drafts/. The starter in scripts/drafts/ is historical scratch work. Current outputs are regenerated under results/sou-q4-final/; older results/sou-q4/ and results/run-output/ are not current.',
    'baseline-test-predictions.csv and final-panels.csv let teammates compare the actual baseline rather than rounded summary numbers. input-fingerprint.csv and session-info.txt record inputs and software. Other CSVs record validation results, choices, coefficients and confusion counts. The saved RData supports inspection without rerunning.'))
)
md <- c("# Question 4 — aligned class-baseline comparison", "", "Author: Sou", "")
html <- c("<h1>Question 4 — aligned class-baseline comparison</h1>","<p>Sou · exploratory analysis</p>")
for(section in sections) {
  md <- c(md,paste0("## ",section$title),"")
  html <- c(html,paste0("<h2>",html_escape(section$title),"</h2>"))
  if(section$type=="table") {
    md <- c(md,md_table(section$body),"");html <- c(html,html_table(section$body))
  } else for(paragraph in section$body) {
    md <- c(md,paragraph,"");html <- c(html,paste0("<p>",html_escape(paragraph),"</p>"))
  }
}
md <- c(md,"![Training cross-validation comparison](cv-comparison.png)")
html <- c(html,"<img src='cv-comparison.png' alt='Training cross-validation comparison'>")
writeLines(md,file.path(output_dir,"sou-q4-report.md"))
writeLines(c("---",'title: "Question 4: aligned class-baseline comparison"',
             'author: "Sou"',"output: html_document","---","",md[-(1:4)]),
           file.path(output_dir,"sou-q4-report.Rmd"))
css <- "body{font:16px/1.6 system-ui,sans-serif;max-width:1100px;margin:40px auto;padding:0 24px;color:#243747}h1,h2{color:#174768}h2{margin-top:32px}.table-wrap{overflow-x:auto}table{border-collapse:collapse;width:100%;font-size:14px}td,th{padding:9px;border-bottom:1px solid #dce4ea;text-align:left}th{background:#edf3f7}img{max-width:100%;height:auto}"
writeLines(c("<!doctype html><html lang='en'><head><meta charset='utf-8'>",
             "<meta name='viewport' content='width=device-width,initial-scale=1'>",
             "<title>Sou: aligned Q4 comparison</title>",paste0("<style>",css,"</style>"),
             "</head><body>",html,"</body></html>"),file.path(output_dir,"sou-q4-report.html"))
team_section <- c("### Question 4: alternative and simpler panels","",
  "For a consistent comparison with Q3 and the class analysis, we used the supplied cleaned dataset, ASD/TD random-forest labels, seed 101422, the same 80/20 split and a fixed 0.5 classification threshold. The class benchmark retains the original full-data t-test/random-forest top-10 intersection. Q4 selected proteins on training subjects and compared the same 24 configurations using three repeats of stratified five-fold training cross-validation. We investigated lasso as the primary alternative and a one- to three-protein panel as a secondary option.","",
  results_text[1],"",md_table(test_display),"",verdict,"",
  sprintf("The simpler panel contained %s. The class panel contained %s. The alternative's full protein list and coefficients can be regenerated by the accompanying Q4 script.",paste(simple_panel,collapse=", "),paste(inclass_panel,collapse=", ")),"",
  results_text[3],"",baseline_note,"",
  "The clean data retain global preprocessing, and the class benchmark retains full-data feature selection; this comparison is aligned with the team but is not fully leakage-free. The test set was seen earlier, making final scores exploratory. With only 31 test children, one classification changes accuracy by approximately 3.2 percentage points. Observed comparability is not proof of equivalence.","")
writeLines(team_section,file.path(output_dir,"sou-q4-team-section.qmd"))
writeLines(c("# Current Q4: aligned comparison","",verdict,"",
             "Use sou-q4-team-section.qmd for the shared report.",
             "Run scripts/drafts/q4_analysis.R to regenerate all current outputs.",
             "There is ONE class baseline in the current-run table.",
             "The team's historical 83.9% is a separate note until reproduced.",
             "Earlier Q4 outputs used different preprocessing and label coding; do not mix them with this revision."),
           file.path(output_dir,"README.md"))
