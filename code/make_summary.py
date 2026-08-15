"""Writes RESULTS_SUMMARY.md for the public repository.

Every figure is read from the analysis outputs rather than typed, so the summary cannot drift
from the results it describes.
"""
import json

J = json.load(open("RESULTS.json"))
K = json.load(open("SENS3.json"))
CL = json.load(open("CLOSURE.json"))
S = json.load(open("SENSITIVITY.json"))
P, X, Q, G, R = J["panel"], J["sens"], J["off"], J["gen"], J["round"]
sp = lambda n: f"{int(n):,}".replace(",", " ")

T = f"""# Results summary

Malignant central nervous system neoplasms and premature mortality in Europe: lifespan
variation and the reference dependence of years of life lost.

Aaron Lawson McLean, Christian Senft, Colin W. McLean.

This file records the design and every headline result produced by the code in this
repository. It stands in place of the manuscript, which is not distributed here. All figures
below are read directly from the analysis outputs in `results/` and `data/derived/`.

---

## Question

Life disparity, denoted e†, is the mean remaining life expectancy at death within a life table.
It measures variation in length of life and is distinct from life expectancy, which measures
average longevity. Removing a cause of death compresses the distribution of ages at death only
when those deaths fall below a threshold age a†, and widens it otherwise. The effect of a cause
on lifespan variation therefore has to be computed rather than inferred from its age profile.

Two questions are addressed:

1. What are the changes in life expectancy at birth and in life disparity implied by removing
   mortality from ICD-10 C70 to C72, and what age structure produces them?
2. How far do cross-country comparisons of years of life lost from this cause depend on the
   reference life table chosen?

## Design

Descriptive cross-national analysis of aggregate routinely collected mortality data. Eurostat
`hlth_cd_aro` supplies deaths coded to C70 to C72 by country of residence, sex and five-year
age group. The Human Mortality Database supplies exposures, all-cause deaths and period life
tables. A country-year-sex stratum is retained only when its five-year age cells sum exactly to
the separately reported all-age total, which establishes that any absent cell contributes zero.

Primary frame: the common period 2012 to 2018 for all {J['frame']['countries']} countries,
comprising {J['frame']['country_years']} of {P['panel_potential']} potentially eligible
both-sex country-years, {sp(J['frame']['deaths'])} C70 to C72 deaths and
{sp(J['frame']['allcause_deaths'])} all-cause deaths. Six country-years fail the reconciliation
rule, all in Iceland and Luxembourg.

Cause-deleted life tables are built from mₓ multiplied by Rₓ, one minus the C70 to C72 share of
all deaths in the age group, using the same function as the observed tables. Threshold ages are
obtained numerically by perturbing mₓ in each age group in turn and are identified as band
boundaries, since the data are grouped in five-year intervals. The effect on e† is decomposed
by age by cumulative removal, reported as the symmetric average of the forward and reverse
orderings. Uncertainty is assessed by resampling age-specific death counts from Poisson
distributions over 10 000 replicates; the resulting ranges are percentile intervals describing
the stability of the estimates and are not confidence intervals.

## Principal results

| Quantity | Value |
|---|---|
| Gain in life expectancy at birth on cause deletion | {J['gain_e0']['mean']} years, about {R['gain_days']} days (range {J['gain_e0']['lo']} to {J['gain_e0']['hi']}) |
| Reduction in life disparity e† | {J['d_edag']['mean']} years, about {R['ded_days']} days (range {J['d_edag']['lo']} to {J['d_edag']['hi']}) |
| Observed e† for comparison | about {J['d_edag']['edag_mean']} years |
| Countries with a positive reduction in e† | {J['d_edag']['positive']} of {J['frame']['countries']} |
| Countries whose 95% percentile interval excludes zero | {J['d_edag']['ci_excl_zero']} of {J['frame']['countries']} |
| Median threshold boundary a† | {J['threshold']['median']} years |
| Share of C70 to C72 deaths below a† | {J['threshold']['below_lo']}% to {J['threshold']['below_hi']}% |
| Mean age at C70 to C72 death | {J['threshold']['mad_lo']} to {J['threshold']['mad_hi']} years |
| Age-standardised mortality rate | {J['asmr']['lo']} to {J['asmr']['hi']} per 100 000 person-years |

Age-specific contributions to the reduction in e† are positive from childhood to about age 75,
peak between 50 and 65 years, and are slightly negative above age 80. The positive
contributions dominate in every country, which is why the net effect is compressive.

## Reference dependence of years of life lost

Pooled across the frame, C70 to C72 accounts for {R['deaths_pct']}% of all deaths and
{R['yll_pct']}% of all-cause years of life lost under national life tables, a ratio of
{J['shares']['ratio']} that exceeds 1.4 in every country (range {J['shares']['ratio_min']} to
{J['shares']['ratio_max']}).

| Quantity | Value |
|---|---|
| Age-standardised YLL rate, national reference | {J['asyr']['nat_lo']} to {J['asyr']['nat_hi']} per 100 000 |
| Age-standardised YLL rate, GBD 2023 reference | {Q['asyr_off_lo']} to {Q['asyr_off_hi']} per 100 000 |
| Increase in YLL per death, GBD 2023 reference | {Q['off_pct_mean']}% mean ({Q['off_pct_lo']}% to {Q['off_pct_hi']}%) |
| Increase in YLL per death, empirical minimum-mortality reference | {X['bp_pct_mean']}% mean |
| Countries moving two or more rank places, GBD reference | {Q['rank_shift_ge2_off']} of {J['frame']['countries']}, largest shift {Q['rank_shift_max_off']} |
| Countries moving two or more rank places, minimum-mortality reference | {X['rank_shift_ge2_bp']} of {J['frame']['countries']} |
| Spearman correlation between the two common references | {Q['rho_off_bp']} |

Life expectancy at birth is {Q['off_e0']} years in the official GBD 2023 schedule and
{X['bp_e0']} years in the empirical minimum-mortality schedule. Ranks were computed from
unrounded rates throughout.

Across nine causes spanning mean ages at death from {G['mad_lo']} to {G['mad_hi']} years and
{G['total_deaths'] / 1e6:.1f} million deaths, the reference-induced increase rises with the mean
age at death, from {G['min_infl']}% for {G['min_infl_cause'].lower()} to {G['max_infl']}% for
{G['max_infl_cause'].lower()} (Spearman ρ = {G['spearman_mad_infl']}, with one inversion).

## Robustness

| Analysis | Result |
|---|---|
| Balanced panel, {K['bal_n']} countries with all seven years | gain {K['bal_gain_e0']}, reduction in e† {K['bal_d_edag']}, against {K['full_gain_e0']} and {K['full_d_edag']} in the full frame; ordering unchanged |
| Later period, {K['lat_n']} countries {K['lat_years']} ({K['lat_country_years']} country-years, {sp(K['lat_deaths'])} deaths) | reduction in e† {K['lat_d_edag']} against {K['latprim_d_edag']} for the same countries in the primary window; positive in all {K['lat_d_edag_positive']}; Spearman {K['lat_rho']} between windows |
| Reference dependence in the later-period frame | {K['lat_ref_shift_ge2']} of {K['lat_n']} countries move two or more places; mean increase {K['lat_ref_pct_mean']}% |
| Observation window, all countries | Spearman {J['frames']['rho_prim_max']} against all available country-years and {J['frames']['rho_prim_excl']} excluding 2020 and 2021; largest shift {J['frames']['max_shift']} places |
| Cause fraction formed with Eurostat all-cause deaths | reduction in e† {X['src_d_edag_es']} against {X['src_d_edag_hmd']}; largest country difference {X['src_max_abs_diff_edag']} |
| Restricted to {X['src_lowdisc_n']} countries with sources agreeing within 1% | reduction in e† {X['src_lowdisc_d_edag']}, positive in all {X['src_lowdisc_positive']} |
| Threshold perturbation of 1%, 0.1% and 0.01% | identical band in every country-year |
| Decomposition ordering reversed | totals identical; signs agree in {S['rev_sign_agreement_pct']}% of age cells; largest single difference {S['rev_max_abs_diff']:.1e} years |
| Open-interval age set to 96 or 100 years | country mean ages at death change by at most {S['openage_max_change_96']} and {S['openage_max_change_100']} years |
| Open interval closed at the published remaining expectancy for exact age 95 | level of e† {CL['edag_alt_mean']} against {CL['edag_reported_mean']} years; mean reduction {CL['ded_alt_scaled']} to {CL['ded_alt_fixed']} against {CL['ded_reported']}, largest country difference {CL['ded_max_country_diff']} years, positive in all {CL['ded_positive_alt_fixed']}; share of all-cause YLL {CL['share_yll_alt_pct']}% against {CL['share_yll_reported_pct']}% |

Sex-specific results: gain in life expectancy {J['male']['gain']} years in males and
{J['female']['gain']} in females; reduction in e† {J['male']['ded']} and {J['female']['ded']};
mean age at death {J['male']['mad']} and {J['female']['mad']} years; median threshold boundary
{J['male']['thresh']} and {J['female']['thresh']} years. Sex-specific strata are more heavily
affected by Eurostat disclosure control, in the extreme leaving a single retained year for one
country in each sex, so country-level sex-specific estimates should be read with that in mind.

Agreement between the two data sources on all-cause totals across retained country-years:
median absolute difference {J['reconciliation']['median_abs']}%, 95th percentile
{J['reconciliation']['p95']}%, maximum {J['reconciliation']['max']}%.

## Sensitivity to the closure of the open age interval

The reconstruction collapses ages above 95 into one open interval and closes it with L = l / m.
The published Human Mortality Database series extends to 110 years and closes there, so
remaining expectancy at exact age 95 differs between the two, and because T(x) is cumulated
downward that difference reaches every younger age. Life expectancy at birth in the
reconstruction exceeds the published value by a mean of {CL['e0_recon_minus_published_mean']}
years and at most {CL['e0_recon_minus_published_max']} years, and life disparity by
{CL['edag_reported_mean'] - CL['edag_alt_mean']:.3f} years.

That difference exceeds the cause-deletion effect being measured. It does not affect the
reported results, because every reported quantity is either a difference between two tables
built by the same procedure or a ratio in which the closure enters numerator and denominator
alike. Recomputing with the open interval closed at the published value moves the mean
reduction in e† from {CL['ded_reported']} to between {CL['ded_alt_scaled']} and
{CL['ded_alt_fixed']} years, with a largest country difference of {CL['ded_max_country_diff']}
years, and the reduction stays positive in all {CL['ded_positive_alt_fixed']} countries. The
gain in life expectancy moves by at most {CL['gain_max_country_diff']} years. The share of
all-cause years of life lost moves from {CL['share_yll_reported_pct']}% to
{CL['share_yll_alt_pct']}%, and the age-standardised years of life lost rate by at most
{CL['asyr_max_country_pct_diff']}% in any country, changing the rank of
{CL['asyr_rank_changes']} of 28 countries by one place.

{CL['cns_deaths_open_pct']}% of C70 to C72 deaths fall in the open interval, against
{CL['allcause_deaths_open_pct']}% of all-cause deaths, which is why the cause-specific results
are less exposed than the life-table levels.

## What the analysis does not establish

These are life-table quantities and not intervention effects. The cause-deleted counterfactual
assumes independence between C70 to C72 mortality and competing causes within age intervals,
makes no allowance for replacement mortality, and does not represent achievable elimination of
these tumours. The assumption cannot be tested in these data.

The analysis covers malignant tumours only; the finest neoplasm breakdown available in the
Eurostat dataset places non-malignant tumours in an undifferentiated aggregate. Cause-of-death
coding for CNS tumours varies internationally, particularly in the assignment of deaths from
brain metastases to the primary site. Threshold ages are band boundaries rather than point
estimates. The primary frame ends in 2018; the later-period frame shows the direction and
magnitude unchanged through 2021 where that comparison is possible.

The resampling intervals describe instability arising from the size and age allocation of the
cause-specific death counts. They do not represent uncertainty in coding practice, in
reconciliation between the two data sources, or in the Human Mortality Database series.
Resampling was performed for the primary frame only; results for other frames are point
estimates.
"""

open("RESULTS_SUMMARY.md", "w").write(T)
print("summary written, %d words" % len(T.split()))
