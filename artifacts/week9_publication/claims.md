# Week 9 claims ledger

Generated from the pinned saved evidence. Numerical wording follows the source row selections and unrounded values in claims.json; selected comparisons are pointwise.

## e1_classifications [supported]

E1 has 19 within-criteria, 22 supported-breach, and 1 unresolved theta/sigma classifications; 21 mu classifications are inapplicable at zero truth.

Qualifications: Pointwise grid statements under the broad 80%-coverage/10%-relative-bias criteria; within criteria does not establish nominal 95% calibration.

Uncertainty: saved pointwise Monte Carlo uncertainty

Comparison: within_cell_summary

Independent unit: source_path; paired floor contrasts use independent coupling groups

## e1_availability [supported]

E1 retains 14,000 source paths and 42,000 estimator records, including 0 invalid points and 0 additional interval rejections (0 unavailable intervals in total).

Qualifications: Point-invalid records are already included in interval-unavailable records. Counts are estimator records, not independent source failures.

Uncertainty: complete retained-record accounting

Comparison: within_cell_summary

Independent unit: source_path; paired floor contrasts use independent coupling groups

## e2_classifications [supported]

E2 has 160 within-criteria, 195 supported-breach, and 5 unresolved theta/sigma classifications; 180 mu classifications are inapplicable at zero truth.

Qualifications: Pointwise grid statements under the broad 80%-coverage/10%-relative-bias criteria; within criteria does not establish nominal 95% calibration.

Uncertainty: saved pointwise Monte Carlo uncertainty

Comparison: within_cell_summary

Independent unit: source_path; paired floor contrasts use independent coupling groups

## e2_availability [supported]

E2 retains 120,000 source paths and 360,000 estimator records, including 85 invalid points and 1 additional interval rejections (86 unavailable intervals in total).

Qualifications: Point-invalid records are already included in interval-unavailable records. Counts are estimator records, not independent source failures.

Uncertainty: complete retained-record accounting

Comparison: within_cell_summary

Independent unit: source_path; paired floor contrasts use independent coupling groups

## e4_classifications [supported]

E4 has 16 within-criteria, 26 supported-breach, and 0 unresolved theta/sigma classifications; 21 mu classifications are inapplicable at zero truth.

Qualifications: Pointwise grid statements under the broad 80%-coverage/10%-relative-bias criteria; within criteria does not establish nominal 95% calibration.

Uncertainty: saved pointwise Monte Carlo uncertainty

Comparison: within_cell_summary

Independent unit: source_path; paired floor contrasts use independent coupling groups

## e4_availability [supported]

E4 retains 14,000 source paths and 42,000 estimator records, including 0 invalid points and 0 additional interval rejections (0 unavailable intervals in total).

Qualifications: Point-invalid records are already included in interval-unavailable records. Counts are estimator records, not independent source failures.

Uncertainty: complete retained-record accounting

Comparison: within_cell_summary

Independent unit: source_path; paired floor contrasts use independent coupling groups

## e2_short_span_exact_theta [supported]

At n=250 and q=0.1, all five exact-theta breaches are bias-driven: relative bias ranges from 17.1% to 22.6%, while operational coverage ranges from 93.6% to 95.3%.

Qualifications: The regular reference is included. The anticipated coverage collapse did not occur; finite-span effects cannot all be attributed to irregularity.

Uncertainty: saved pointwise Monte Carlo uncertainty

Comparison: within_cell_summary

Independent unit: source_path; paired floor contrasts use independent coupling groups

## e2_nonmonotone_pfml_theta [supported]

At n=250 and q=0.1, PFML theta classifications across ascending CV are supported_breach, supported_breach, supported_breach, within_criteria, unresolved.

Qualifications: A single maximum acceptable CV or a monotone interpolated operating boundary would misrepresent these observed classifications.

Uncertainty: saved pointwise classifications

Comparison: within_cell_summary

Independent unit: source_path; paired floor contrasts use independent coupling groups

## e2_nonmonotone_euler_sigma [supported]

At n=250 and q=0.1, EULER sigma classifications across ascending CV are supported_breach, supported_breach, supported_breach, unresolved, within_criteria.

Qualifications: A single maximum acceptable CV or a monotone interpolated operating boundary would misrepresent these observed classifications.

Uncertainty: saved pointwise classifications

Comparison: within_cell_summary

Independent unit: source_path; paired floor contrasts use independent coupling groups

## headline_unresolved [supported]

All 6 unresolved E1/E2 parameter-specific classifications are retained separately from supported breaches and within-criteria results.

Qualifications: Overlapping E1/E2 nominal cells use independent seeds and are not pooled or selected between to obtain a preferred classification.

Uncertainty: saved pointwise Monte Carlo uncertainty

Comparison: within_cell_summary

Independent unit: source_path; paired floor contrasts use independent coupling groups

## e2_operational_conditional_example [supported]

For regular exact theta at n=250,q=2, 1832 covering intervals give 91.60% operational coverage over 2000 planned paths versus 93.42% conditional coverage over 1961 available intervals.

Qualifications: Unavailable intervals remain noncoverage operationally; bias/RMSE use point-valid fits.

Uncertainty: saved pointwise Monte Carlo uncertainty

Comparison: within_cell_summary

Independent unit: source_path; paired floor contrasts use independent coupling groups

## e2_irregularity_mixed_rmse [supported]

At n=4000,q=2,CV=3, exact-estimator RMSE ratios to the regular reference are 0.815 for theta, 0.445 for sigma and 1.682 for mu.

Qualifications: The nominal design is matched; realized spans differ. This is evidence against a blanket assertion that irregularity worsens every parameter.

Uncertainty: saved pointwise Monte Carlo uncertainty

Comparison: independent_cells_nominal_design_reference

Independent unit: source_path; paired floor contrasts use independent coupling groups

## e2_exact_not_uniform_lowest_rmse [descriptive]

At n=250,q=0.1,CV=3, observed theta RMSE is 0.474796 for exact, 0.353144 for PFML and 0.315932 for Euler.

Qualifications: This descriptive point-error ranking does not establish a general ordering. Within-source inferential method contrasts require the saved paired-comparison uncertainty; practical classification assesses bias and coverage separately.

Uncertainty: saved per-method RMSE MCSE; this sentence asserts no paired significance test

Comparison: three_estimators_on_shared_source_paths

Independent unit: source_path; paired floor contrasts use independent coupling groups

## e3_followup_accounting [supported]

The independent exploratory follow-up retains 20,000 paths and 60,000 fits, with 421 invalid points and 23 additional interval rejections (444 unavailable intervals).

Qualifications: These exploratory results do not overwrite E1/E2 classifications. One exact point remains boundary_suspected.

Uncertainty: complete retained-record accounting

Comparison: within_cell_summary

Independent unit: source_path; paired floor contrasts use independent coupling groups

## e3_selected_replay_scope [supported]

The selected-path diagnostic archive evaluated 53 sources and verified 159 original fit replays.

Qualifications: This selected audit does not independently recalculate all diagnostic mathematics, prove global optimality, or replace headline outcomes.

Uncertainty: complete fixed-selection accounting; no prevalence inference

Comparison: original_fit_replay_at_fixed_selected_sources

Independent unit: selected source; selection is not a random frequency study

## e3_floor_sensitivity [supported]

Across the specified paired floor contrasts, the largest estimated RMSE-ratio departure from one is 1.00023261, with pointwise MC interval [0.99974475, 1.00072070]; the largest absolute observed coverage difference is 0.05 percentage points.

Qualifications: Small effects within this grid do not establish equivalence or universal floor invariance. The physical realizations are not all independent.

Uncertainty: saved pointwise Monte Carlo uncertainty

Comparison: paired_floor_shared_draws_and_standardized_innovations

Independent unit: source_path; paired floor contrasts use independent coupling groups

## e4_pfml_mu_r5 [supported]

At r=5, PFML mu has MMPP/IID RMSE ratio 1.083874 [95% MC 1.037545, 1.132272] and an operational-coverage difference of -2.70 percentage points [-4.158, -1.242].

Qualifications: Exploratory fixed design including realized-span consequences. This contrast does not establish a mean-bias shift or a universal clustering effect. PFML intervals use the registered working-model covariance.

Uncertainty: pointwise normal MC interval for coverage difference; log-delta MC interval for RMSE ratio

Comparison: independent_MMPP_minus_full_marginal_IID_cells

Independent unit: source_path; paired floor contrasts use independent coupling groups

## e4_pfml_mu_r20 [supported]

At r=20, PFML mu has MMPP/IID RMSE ratio 1.157863 [95% MC 1.107265, 1.210774] and an operational-coverage difference of -3.85 percentage points [-5.466, -2.234].

Qualifications: Exploratory fixed design including realized-span consequences. This contrast does not establish a mean-bias shift or a universal clustering effect. PFML intervals use the registered working-model covariance.

Uncertainty: pointwise normal MC interval for coverage difference; log-delta MC interval for RMSE ratio

Comparison: independent_MMPP_minus_full_marginal_IID_cells

Independent unit: source_path; paired floor contrasts use independent coupling groups

## e4_null_contrasts [supported]

All 27 of 27 matched bias intervals contain zero; all 18 of 18 theta/sigma RMSE-ratio intervals contain one.

Qualifications: Intervals containing the null do not demonstrate equivalence. These comparisons are pointwise, not a familywise guarantee.

Uncertainty: saved pointwise Monte Carlo uncertainty

Comparison: independent_MMPP_minus_full_marginal_IID_cells

Independent unit: source_path; paired floor contrasts use independent coupling groups

## e4_isolated_exact_sigma_coverage [supported]

At r=20, the exact-sigma operational-coverage contrast is +1.40 percentage points, with pointwise MC interval [+0.02457, +2.77543].

Qualifications: This isolated pointwise result among multiple comparisons does not establish a general beneficial clustering effect.

Uncertainty: saved pointwise Monte Carlo uncertainty

Comparison: independent_MMPP_minus_full_marginal_IID_cells

Independent unit: source_path; paired floor contrasts use independent coupling groups

## e4_pfml_mu_calibration [descriptive]

At r=20, PFML mu empirical SD / RMS reported SE is 1.145 under MMPP and 0.994 under IID sampling.

Qualifications: Both quantities use the same interval-valid fits. This descriptive comparison is not a derived dependence-robust covariance benchmark.

Uncertainty: descriptive saved variability estimates; no new MC interval claimed

Comparison: within_cell_summary

Independent unit: source_path; paired floor contrasts use independent coupling groups

## shared_expected_objective_targets [theoretical_under_stated_assumptions]

Under the stated exogeneity and stationary initialization assumptions, matched full one-gap marginals share the derived expected-objective minimizers; this does not imply equal finite-sample variance or coverage.

Qualifications: Population targets are interpretation benchmarks; all simulation error summaries continue to use the true OU parameters.

Uncertainty: no covariance or coverage theorem claimed

Comparison: stationary_expected_objective_calculation

Independent unit: not applicable

## beyond_first_order_theory [unsupported_claim_excluded]

A demonstrated departure from a fully matched first-order coverage theory is not established by the current experiments or population-target calculation.

Qualifications: Population targets are interpretation benchmarks; all simulation error summaries continue to use the true OU parameters.

Uncertainty: no covariance or coverage theorem claimed

Comparison: stationary_expected_objective_calculation

Independent unit: not applicable

## general_dependent_gap_covariance [unsupported_claim_excluded]

A general dependent-gap CLT and a dependence-robust PFML/Euler covariance formula are not established in this study.

Qualifications: Population targets are interpretation benchmarks; all simulation error summaries continue to use the true OU parameters.

Uncertainty: no covariance or coverage theorem claimed

Comparison: stationary_expected_objective_calculation

Independent unit: not applicable
