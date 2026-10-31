# Week 9 publication captions

## Newly rendered figures

### F1_distributions_main

F1. Empirical cumulative distributions of all eligible saved E1 estimates at nominal CV 0, 1 and 3. All three parameters are estimated jointly; truth is (1,0,0.5), q=0.5, n=1000 and R=2000 per cell. Dotted vertical lines denote truth. Counts show point-valid/planned fits; no confidence bands are claimed. Regular-grid exact/PFML estimates coincide, so their ECDFs overlap; line styles distinguish coincident methods. Selection is post-result presentation. All seven CV settings are supplied in supporting figures and f1_ecdf.csv.gz.

### S_F1_distributions_cv0_to_1

Supporting F1, first page. Complete saved E1 ECDFs at nominal CV 0, 0.25, 0.5 and 1, with the same definitions as F1. Every eligible estimate is retained, including extremes.

### S_F1_distributions_cv1p5_to_3

Supporting F1, second page. Complete saved E1 ECDFs at nominal CV 1.5, 2 and 3. These input Gamma gap laws have the registered floor; realized spans are not held fixed.

### F7_classification_n250

F7, n=250. Complete registered E2 pointwise classifications over the 4 by 5 nominal q/CV grid for each method and theta/sigma. W means within the broad practical criteria; B:b, B:c and B:bc denote supported bias, coverage or both breaches; U remains unresolved. Classification uses unrounded bias Monte Carlo intervals and operational-coverage Wilson bounds. No interpolation, continuous boundary, monotonicity or simultaneous guarantee is asserted. The full 360 eligible E2 rows are in T1_e2.csv; mu's zero-truth relative-bias criterion is inapplicable.

### S_calibration_e2_n250

Supporting E2 calibration, n=250. Complete-grid empirical SD divided by RMS reported SE on the same interval-valid subset. Cell labels show the ratio; the common color scale is symmetric in log2 and saturates outside [0.5,2], without clipping the numerical labels. This is descriptive calibration, not a robust variance formula or a coverage explanation by itself. All counts and unrounded values are in calibration.csv.

### F7_classification_n1000

F7, n=1000. Complete registered E2 pointwise classifications over the 4 by 5 nominal q/CV grid for each method and theta/sigma. W means within the broad practical criteria; B:b, B:c and B:bc denote supported bias, coverage or both breaches; U remains unresolved. Classification uses unrounded bias Monte Carlo intervals and operational-coverage Wilson bounds. No interpolation, continuous boundary, monotonicity or simultaneous guarantee is asserted. The full 360 eligible E2 rows are in T1_e2.csv; mu's zero-truth relative-bias criterion is inapplicable.

### S_calibration_e2_n1000

Supporting E2 calibration, n=1000. Complete-grid empirical SD divided by RMS reported SE on the same interval-valid subset. Cell labels show the ratio; the common color scale is symmetric in log2 and saturates outside [0.5,2], without clipping the numerical labels. This is descriptive calibration, not a robust variance formula or a coverage explanation by itself. All counts and unrounded values are in calibration.csv.

### F7_classification_n4000

F7, n=4000. Complete registered E2 pointwise classifications over the 4 by 5 nominal q/CV grid for each method and theta/sigma. W means within the broad practical criteria; B:b, B:c and B:bc denote supported bias, coverage or both breaches; U remains unresolved. Classification uses unrounded bias Monte Carlo intervals and operational-coverage Wilson bounds. No interpolation, continuous boundary, monotonicity or simultaneous guarantee is asserted. The full 360 eligible E2 rows are in T1_e2.csv; mu's zero-truth relative-bias criterion is inapplicable.

### S_calibration_e2_n4000

Supporting E2 calibration, n=4000. Complete-grid empirical SD divided by RMS reported SE on the same interval-valid subset. Cell labels show the ratio; the common color scale is symmetric in log2 and saturates outside [0.5,2], without clipping the numerical labels. This is descriptive calibration, not a robust variance formula or a coverage explanation by itself. All counts and unrounded values are in calibration.csv.

### F8_matched_mean

F8. Exploratory E4 matched-marginal MMPP/IID comparisons for mu across r=1,5,20 and all methods. Ratios are MMPP/IID RMSE, differences are MMPP minus IID operational coverage in percentage points. Bars are pointwise 95% Monte Carlo intervals using independent-cell formulas; RMSE-ratio limits use the log-delta method. n=1000, q=0.5 and R=2000 per independent source cell. The contrast includes realized-span consequences; PFML/Euler Wald covariance is model-based. Displaying mu prominently is post-result presentation; all parameters remain supplied.

### S_calibration_e4_mean

Supporting E4 calibration. Empirical SD and RMS reported SE for mu on the identical interval-valid subsets, with the same y-axis limits. These saved descriptive quantities do not have a derived Monte Carlo interval here. All methods, r values and both matched arms are retained. At r=20 both arms use the declared 5e-7 floor; MMPP labels denote floored MMPP-derived gaps. They do not establish a long-run covariance theorem.

## Byte-identical source figures

These figures retain their original titles, definitions, uncertainty, source manifests and historical review identities. Their bytes are copied unchanged. They are supporting evidence and the existing F2–F6 figure families; each original review is copied beside the source figures. The source study's sealed report defines its original inspection scope.

- `figures/source/e1/F2_availability.png` ← `artifacts/week5_e1/figures/F2_availability.png` (e1); original review: `artifacts/week5_e1/visual_review.json`.
- `figures/source/e1/F2_coverage.png` ← `artifacts/week5_e1/figures/F2_coverage.png` (e1); original review: `artifacts/week5_e1/visual_review.json`.
- `figures/source/e1/F3_mu_absolute.png` ← `artifacts/week5_e1/figures/F3_mu_absolute.png` (e1); original review: `artifacts/week5_e1/visual_review.json`.
- `figures/source/e1/F3_relative_bias.png` ← `artifacts/week5_e1/figures/F3_relative_bias.png` (e1); original review: `artifacts/week5_e1/visual_review.json`.
- `figures/source/e1/F4_rmse_ratios.png` ← `artifacts/week5_e1/figures/F4_rmse_ratios.png` (e1); original review: `artifacts/week5_e1/visual_review.json`.
- `figures/source/e2/F5_operational_n1000.png` ← `artifacts/week6_e2/figures/F5_operational_n1000.png` (e2); original review: `artifacts/week6_e2/visual_review.json`.
- `figures/source/e2/F5_operational_n250.png` ← `artifacts/week6_e2/figures/F5_operational_n250.png` (e2); original review: `artifacts/week6_e2/visual_review.json`.
- `figures/source/e2/F5_operational_n4000.png` ← `artifacts/week6_e2/figures/F5_operational_n4000.png` (e2); original review: `artifacts/week6_e2/visual_review.json`.
- `figures/source/e2/F6_availability_n1000.png` ← `artifacts/week6_e2/figures/F6_availability_n1000.png` (e2); original review: `artifacts/week6_e2/visual_review.json`.
- `figures/source/e2/F6_availability_n250.png` ← `artifacts/week6_e2/figures/F6_availability_n250.png` (e2); original review: `artifacts/week6_e2/visual_review.json`.
- `figures/source/e2/F6_availability_n4000.png` ← `artifacts/week6_e2/figures/F6_availability_n4000.png` (e2); original review: `artifacts/week6_e2/visual_review.json`.
- `figures/source/e2/F6_conditional_n1000.png` ← `artifacts/week6_e2/figures/F6_conditional_n1000.png` (e2); original review: `artifacts/week6_e2/visual_review.json`.
- `figures/source/e2/F6_conditional_n250.png` ← `artifacts/week6_e2/figures/F6_conditional_n250.png` (e2); original review: `artifacts/week6_e2/visual_review.json`.
- `figures/source/e2/F6_conditional_n4000.png` ← `artifacts/week6_e2/figures/F6_conditional_n4000.png` (e2); original review: `artifacts/week6_e2/visual_review.json`.
- `figures/source/e3_floor/W7_floor_paired_bias.png` ← `artifacts/week7_floor/figures/W7_floor_paired_bias.png` (e3_floor); original review: `artifacts/week7_floor/visual_review.json`.
- `figures/source/e3_floor/W7_floor_paired_coverage.png` ← `artifacts/week7_floor/figures/W7_floor_paired_coverage.png` (e3_floor); original review: `artifacts/week7_floor/visual_review.json`.
- `figures/source/e3_floor/W7_floor_rmse_ratio.png` ← `artifacts/week7_floor/figures/W7_floor_rmse_ratio.png` (e3_floor); original review: `artifacts/week7_floor/visual_review.json`.
- `figures/source/e3_floor/W7_floor_spacing.png` ← `artifacts/week7_floor/figures/W7_floor_spacing.png` (e3_floor); original review: `artifacts/week7_floor/visual_review.json`.
- `figures/source/e3_followup/W7_followup_availability.png` ← `artifacts/week7_followup/figures/W7_followup_availability.png` (e3_followup); original review: `artifacts/week7_followup/visual_review.json`.
- `figures/source/e3_followup/W7_followup_coverage.png` ← `artifacts/week7_followup/figures/W7_followup_coverage.png` (e3_followup); original review: `artifacts/week7_followup/visual_review.json`.
- `figures/source/e4/W8_availability.png` ← `artifacts/week8_e4/figures/W8_availability.png` (e4); original review: `artifacts/week8_e4/visual_review.json`.
- `figures/source/e4/W8_coverage_mu.png` ← `artifacts/week8_e4/figures/W8_coverage_mu.png` (e4); original review: `artifacts/week8_e4/visual_review.json`.
- `figures/source/e4/W8_coverage_sigma.png` ← `artifacts/week8_e4/figures/W8_coverage_sigma.png` (e4); original review: `artifacts/week8_e4/visual_review.json`.
- `figures/source/e4/W8_coverage_theta.png` ← `artifacts/week8_e4/figures/W8_coverage_theta.png` (e4); original review: `artifacts/week8_e4/visual_review.json`.
- `figures/source/e4/W8_matched_bias_mu.png` ← `artifacts/week8_e4/figures/W8_matched_bias_mu.png` (e4); original review: `artifacts/week8_e4/visual_review.json`.
- `figures/source/e4/W8_matched_bias_sigma.png` ← `artifacts/week8_e4/figures/W8_matched_bias_sigma.png` (e4); original review: `artifacts/week8_e4/visual_review.json`.
- `figures/source/e4/W8_matched_bias_theta.png` ← `artifacts/week8_e4/figures/W8_matched_bias_theta.png` (e4); original review: `artifacts/week8_e4/visual_review.json`.
- `figures/source/e4/W8_matched_coverage_mu.png` ← `artifacts/week8_e4/figures/W8_matched_coverage_mu.png` (e4); original review: `artifacts/week8_e4/visual_review.json`.
- `figures/source/e4/W8_matched_coverage_sigma.png` ← `artifacts/week8_e4/figures/W8_matched_coverage_sigma.png` (e4); original review: `artifacts/week8_e4/visual_review.json`.
- `figures/source/e4/W8_matched_coverage_theta.png` ← `artifacts/week8_e4/figures/W8_matched_coverage_theta.png` (e4); original review: `artifacts/week8_e4/visual_review.json`.
- `figures/source/e4/W8_matched_rmse_ratio_mu.png` ← `artifacts/week8_e4/figures/W8_matched_rmse_ratio_mu.png` (e4); original review: `artifacts/week8_e4/visual_review.json`.
- `figures/source/e4/W8_matched_rmse_ratio_sigma.png` ← `artifacts/week8_e4/figures/W8_matched_rmse_ratio_sigma.png` (e4); original review: `artifacts/week8_e4/visual_review.json`.
- `figures/source/e4/W8_matched_rmse_ratio_theta.png` ← `artifacts/week8_e4/figures/W8_matched_rmse_ratio_theta.png` (e4); original review: `artifacts/week8_e4/visual_review.json`.
- `figures/source/e4/W8_source_adjacent_product.png` ← `artifacts/week8_e4/figures/W8_source_adjacent_product.png` (e4); original review: `artifacts/week8_e4/visual_review.json`.
- `figures/source/e4/W8_source_clipping.png` ← `artifacts/week8_e4/figures/W8_source_clipping.png` (e4); original review: `artifacts/week8_e4/visual_review.json`.
- `figures/source/e4/W8_source_cv.png` ← `artifacts/week8_e4/figures/W8_source_cv.png` (e4); original review: `artifacts/week8_e4/visual_review.json`.
- `figures/source/e4/W8_source_mean.png` ← `artifacts/week8_e4/figures/W8_source_mean.png` (e4); original review: `artifacts/week8_e4/visual_review.json`.
- `figures/source/e4/W8_source_span_sd.png` ← `artifacts/week8_e4/figures/W8_source_span_sd.png` (e4); original review: `artifacts/week8_e4/visual_review.json`.
- `figures/source/e3_selected/W7_selected_E2_cv0_rep001424.png` ← `artifacts/week7_selected_figures/W7_selected_E2_cv0_rep001424.png` (e3_selected); original review: `artifacts/week7_evidence/selected_visual_review.json`.
- `figures/source/e3_selected/W7_selected_E2_cv0p5_rep000465.png` ← `artifacts/week7_selected_figures/W7_selected_E2_cv0p5_rep000465.png` (e3_selected); original review: `artifacts/week7_evidence/selected_visual_review.json`.
