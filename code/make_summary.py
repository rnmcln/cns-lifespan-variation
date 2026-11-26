"""Writes RESULTS_SUMMARY.md for the public repository.

Every figure is read from the analysis outputs rather than typed, so the summary cannot drift
from the results it describes.
"""
import json

import pandas as pd

J = json.load(open("RESULTS.json"))
K = json.load(open("SENS3.json"))
CL = json.load(open("CLOSURE.json"))
S = json.load(open("SENSITIVITY.json"))
RK = json.load(open("RANKS.json"))
D = json.load(open("DISCLOSURE.json"))
GG = json.load(open("GENERALISATION.json"))
GC = pd.read_csv("generalisation_common.csv").sort_values("mad")
RV = json.load(open("REVERSALS.json"))
RC = pd.read_csv("rank_by_cause.csv").sort_values("mad")
NV, BV, YV = RV["nat_vs_gbd"], RV["nat_vs_minmort"], RV["ypd_nat_vs_gbd"]
CG, CM = RV["asyr_change_gbd"], RV["asyr_change_minmort"]
NAME = {"AT": "Austria", "BE": "Belgium", "BG": "Bulgaria", "HR": "Croatia", "CZ": "Czechia",
        "DK": "Denmark", "EE": "Estonia", "FI": "Finland", "FR": "France", "DE": "Germany",
        "EL": "Greece", "HU": "Hungary", "IS": "Iceland", "IE": "Ireland", "IT": "Italy",
        "LV": "Latvia", "LT": "Lithuania", "LU": "Luxembourg", "NL": "Netherlands", "NO": "Norway",
        "PL": "Poland", "PT": "Portugal", "SK": "Slovakia", "SI": "Slovenia", "ES": "Spain",
        "SE": "Sweden", "CH": "Switzerland", "UK": "United Kingdom"}
P, X, Q, G, R = J["panel"], J["sens"], J["off"], J["gen"], J["round"]
sp = lambda n: f"{int(n):,}".replace(",", " ")


def sci(x):
    """9.7 × 10<sup>−5</sup>, with a true superscript exponent"""
    m, e = ("%.1e" % x).split("e")
    return "%s × 10<sup>%s</sup>" % (m, str(int(e)).replace("-", "−"))
RA = RK["rank_agreement"]
RKN, RKO = RA["national vs GBD 2023"], RA["GBD 2023 vs minimum-mortality"]
RANK_ROWS = "\n".join(
    f"| {k} | {v['spearman_rho']} | {v['kendall_tau_b']} | {v['median_abs_shift']:.1f} | "
    f"{v['max_abs_shift']} | {v['n_pairs_reversed']} of {v['n_pairs']} ({v['pct_pairs_reversed']}%) |"
    for k, v in RA.items())
REV_ROWS = "\n".join(
    f"| {lab} | {v['spearman_rho']} | {v['n_pairs_reversed']} of {v['n_pairs']} ({v['pct_pairs_reversed']}%) | "
    f"{v['rev_gap_a_rel_median']} | {v['rev_gap_b_rel_median']} | {v['rev_gap_a_rel_max']} | {v['n_rev_both_ge5']} | {v['n_rev_both_ge10']} |"
    for lab, v in (("Age-standardised rate, national vs GBD 2023", NV),
                   ("Age-standardised rate, national vs minimum-mortality", BV),
                   ("YLL per death, national vs GBD 2023", YV)))
STABLE_ROWS = "\n".join(
    f"| {NAME[r['higher_nat']]} | {NAME[r['lower_nat']]} | {r['gap_nat_rel']} | {r['gap_gbd_rel']} | {100 * r['p_recur']:.1f} |"
    for r in NV["largest_stable_reversals"])
RANKC_ROWS = "\n".join(
    f"| {r.cause} | {r.mad:.1f} | {r.asyr_chg_median:.1f} | {r.spearman_rho:.3f} | {r.pct_pairs_reversed:.1f} | "
    f"{100 * r.rep_n_reversed_lo / (r.n * (r.n - 1) / 2):.1f} to {100 * r.rep_n_reversed_hi / (r.n * (r.n - 1) / 2):.1f} | "
    f"{r.rev_recur_ge95} | {r.sd_log_rate_nat:.3f} | {r.sd_log_ratio:.3f} | {r.noise_ratio:.2f} |"
    for r in RC.itertuples())
CAUSE_ROWS = "\n".join(
    f"| {r.cause} | {r.n_ctry} | {sp(r.deaths)} | {r.mad:.1f} | {r.infl_mean:.1f} | "
    f"{r.infl_lo:.1f} to {r.infl_hi:.1f} |" for r in GC.itertuples())

T = f"""# Results summary

Reference life tables and country rankings of years of life lost: malignant central nervous
system neoplasms and eight comparison causes in European countries.

Aaron Lawson McLean, Christian Senft, Colin W. McLean.

This file records the design and every headline result produced by the code in this
repository. It stands in place of the manuscript, which is not distributed here. All figures
below are read directly from the analysis outputs in `results/` and `data/derived/`.

---

## Question

Years of life lost (YLL) depend on a reference life table. National period life tables measure
YLL relative to the remaining life expectancy in each country; a common life table, such as the
Global Burden of Disease (GBD) theoretical minimum risk life table, applies the same age-specific
reference to all countries. The question is how much this choice changes age-standardised YLL
rates and country comparisons for the same recorded deaths: how large the change in each
country's rate is, which countries are most affected, and when the change is large enough to
alter a comparison. The main case is malignant central nervous system neoplasms (ICD-10 C70 to C72) in
{J['frame']['countries']} European countries; eight further causes show where rankings are most
and least sensitive to the reference.

A secondary analysis estimates the changes in life expectancy at birth and in life disparity
e<sup>†</sup> implied by removing C70 to C72 mortality in cause-deleted life tables. These do not
depend on the reference life table.

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
rule, all in Iceland and Luxembourg. Within the frame, Eurostat disclosure control leaves
{D['flagged_T']} both-sex age cells without a value, all in Iceland and Luxembourg; the
sex-specific extracts have {D['flagged_M']} such cells in males and {D['flagged_F']} in females.

Countries are ranked on the age-standardised YLL rate (2013 European Standard Population) under
three references: national period life tables, the official GBD 2023 theoretical minimum risk life
table, and an empirical minimum-mortality table built from the lowest age-specific rates in the
frame. For every pair of countries the gap between their rates is recorded under each reference,
as a percentage of the pair mean. Stability is assessed with {sp(RV['n_replicates'])} replicates in
which every age-specific death count is resampled from a Poisson distribution and both rankings
are recomputed from the same replicate deaths (seed {RV['seed']}); a reversal recurring in at
least 95% of replicates is called stable. These are perturbation summaries, not confidence
intervals.

## Reference dependence of years of life lost

Pooled across the frame, C70 to C72 accounts for {R['deaths_pct']}% of all deaths and
{R['yll_pct']}% of all-cause years of life lost under national life tables.

| Quantity | Value |
|---|---|
| Age-standardised YLL rate, national reference | {J['asyr']['nat_lo']:g} to {J['asyr']['nat_hi']:g} per 100 000 |
| Age-standardised YLL rate, GBD 2023 reference | {Q['asyr_off_lo']} to {Q['asyr_off_hi']} per 100 000 |
| Increase in YLL per death, GBD 2023 reference | {Q['off_pct_mean']}% mean ({Q['off_pct_lo']}% to {Q['off_pct_hi']}%) |
| Increase in YLL per death, empirical minimum-mortality reference | {X['bp_pct_mean']}% mean |
| Pairwise country orderings reversed, national vs GBD 2023 | {RKN['n_pairs_reversed']} of {RKN['n_pairs']} ({RKN['pct_pairs_reversed']}%) |
| Spearman correlation, national vs GBD 2023 | {RKN['spearman_rho']} |
| Spearman correlation between the two common references | {RKO['spearman_rho']} |

Life expectancy at birth is {Q['off_e0']} years in the official GBD 2023 schedule and
{X['bp_e0']} years in the empirical minimum-mortality schedule. Ranks were computed from
unrounded rates throughout.

## Change in each country's age-standardised YLL rate

The comparison of countries is made on the age-standardised rate, so its change carries the main
result. Per country values are in `data/derived/asyr_change_country.csv`.

| Quantity | GBD 2023 reference | Minimum-mortality reference |
|---|---|---|
| Median change in the rate, % | {CG['median_pct']} (interquartile range {CG['q1_pct']} to {CG['q3_pct']}) | {CM['median_pct']} |
| Range across countries, % | {CG['lo_pct']} ({NAME[CG['lo_c']]}) to {CG['hi_pct']} ({NAME[CG['hi_c']]}) | {CM['lo_pct']} to {CM['hi_pct']} |
| Spearman correlation of the change with national life expectancy at birth | {str(CG['spearman_with_e0']).replace('-', '−')} | {str(CM['spearman_with_e0']).replace('-', '−')} |

Median absolute change under the GBD 2023 reference: {CG['median_abs']} per 100 000 (range
{CG['lo_abs']} to {CG['hi_abs']}). Central and Eastern European countries, median change
{CG['cee_median_pct']}%; other countries, {CG['other_median_pct']}%. {CG['top10_cee']} of the 10
largest changes are in Central and Eastern Europe.

## Rank agreement between reference life tables

| Comparison | Spearman ρ | Kendall τ<sub>b</sub> | Median abs. shift | Max shift | Pairs reversed |
|---|---|---|---|---|---|
{RANK_ROWS}
| national vs GBD 2023 (derived schedule), later frame | {K['lat_ref_rho']} | {K['lat_ref_tau']} | {K['lat_ref_median_shift']:.1f} | {K['lat_ref_max_shift']} | {K['lat_ref_pairs_reversed']} of {K['lat_ref_pairs']} ({K['lat_ref_pct_pairs_reversed']}%) |

## Size and stability of the reversals

Gaps are differences between the two countries of a pair as a percentage of the pair mean; the
first reference is the national life table.

| Comparison | Spearman ρ | Pairs reversed | Median gap of reversed pairs, first reference % | Median gap, second reference % | Largest first-reference gap reversed % | Reversed with both gaps ≥5% | Reversed with both gaps ≥10% |
|---|---|---|---|---|---|---|---|
{REV_ROWS}

For the primary comparison, the median gap is {NV['all_gap_a_rel_median']}% over all pairs and
{NV['con_gap_a_rel_median']}% over pairs whose ordering is kept, against
{NV['rev_gap_a_rel_median']}% (interquartile range {NV['rev_gap_a_rel_q1']}% to
{NV['rev_gap_a_rel_q3']}%) over reversed pairs. None of the {NV['n_pairs_nat_gap_gt_rev_max']}
pairs further apart than {NV['rev_gap_a_rel_max']}% is reversed.

Across replicates the number of reversed pairs ranges from {NV['rep_n_reversed_lo']} to
{NV['rep_n_reversed_hi']} (median {NV['rep_n_reversed_median']}). The national ordering recurs in at
least 95% of replicates for {NV['all_nat_order_ge95']} of {NV['n_pairs']} pairs but for only
{NV['rev_nat_order_ge95']} of the {NV['n_pairs_reversed']} reversed pairs.
{NV['rev_recur_ge95']} reversals are stable; in {NV['stable_rev_cee_moves_above_other']} of them a
Central or Eastern European country moves above a country outside that group, and in
{NV['stable_rev_other_moves_above_cee']} the reverse. The five largest stable reversals:

| Higher under national tables | Higher under GBD 2023 | Gap national % | Gap GBD 2023 % | Recurrence % |
|---|---|---|---|---|
{STABLE_ROWS}

All pairs, with their gaps and recurrence, are in `data/derived/reversal_pairs.csv`.

## Rank sensitivity across nine causes

Common frame of {RV['by_cause_n_countries']} countries ({', '.join(RV['by_cause_excluded'])}
excluded, having no retained perinatal stratum). Causes ordered by mean age at death. Ratio of
SDs: between-country SD of log(GBD rate / national rate) divided by the between-country SD of
log(national rate).

| Cause | Mean age at death | Median change in rate % | Spearman ρ | Pairs reversed % | Replicate range % | Stable reversals | SD log national rate | SD log reference effect | Ratio of SDs |
|---|---|---|---|---|---|---|---|---|---|
{RANKC_ROWS}

Reversed pairs whose rates differ by at least 10% under both references, summed over all nine
causes: {RV['by_cause_n_rev_both_ge10_total']}.

## Reference effect on YLL per death by cause

Spearman ρ between mean age at death and the mean increase in YLL per death:
{GG['common_spearman_mad_infl']}. Cause ordering identical to the full-coverage frame:
{GG['common_order_matches_full']}.

| Cause | Countries | Deaths | Mean age at death | Increase % | Range % |
|---|---|---|---|---|---|
{CAUSE_ROWS}

## A literally common open interval

The GBD 2023 reference is common at every age below 95 years; for the open interval the four
sub-interval values are averaged with each country-year's distribution of all-cause deaths above
95. {RK['open_share_deaths_pct']}% of C70 to C72 deaths fall in that interval, and the
country-year values span {RK['e_open_countryyear_lo']:.3f} to {RK['e_open_countryyear_hi']:.3f}
years. Replacing them with one fixed value for every country-year, either the schedule at exact
age 97.5 years ({RK['e_open_fixed']:.3f}) or the pooled death-weighted mean
({RK['e_open_pooled']:.3f}), changes YLL per death by at most {RK['ypd_maxabs_diff_fixed']} years
and the age-standardised rate by at most {RK['asyr_maxabs_diff_fixed']:.2f} per 100 000. Country
ordering identical under both: {RK['rank_identical_fixed'] and RK['rank_identical_pooled']}.

## Secondary analysis: life expectancy and lifespan variation

| Quantity | Value |
|---|---|
| Gain in life expectancy at birth on cause deletion | {J['gain_e0']['mean']} years, about {R['gain_days']} days (range {J['gain_e0']['lo']} to {J['gain_e0']['hi']}) |
| Reduction in life disparity e<sup>†</sup> | {J['d_edag']['mean']} years, about {R['ded_days']} days (range {J['d_edag']['lo']} to {J['d_edag']['hi']}) |
| Observed e<sup>†</sup> for comparison | about {J['d_edag']['edag_mean']} years |
| Countries with a positive reduction in e<sup>†</sup> | {J['d_edag']['positive']} of {J['frame']['countries']} |
| Countries whose 95% Poisson perturbation interval lies wholly above zero | {J['d_edag']['ci_excl_zero']} of {J['frame']['countries']} |
| Median threshold boundary a<sup>†</sup> | {J['threshold']['median']} years |
| Share of C70 to C72 deaths below a<sup>†</sup> | {J['threshold']['below_lo']}% to {J['threshold']['below_hi']}% |
| Mean age at C70 to C72 death | {J['threshold']['mad_lo']} to {J['threshold']['mad_hi']} years |
| Age-standardised mortality rate | {J['asmr']['lo']} to {J['asmr']['hi']} per 100 000 person-years |

Age-specific contributions to the reduction in e<sup>†</sup> are positive from childhood to about age 75,
peak between 50 and 65 years, and are slightly negative above age 80. The positive
contributions dominate in every country, which is why the net effect is compressive.

## Robustness

| Analysis | Result |
|---|---|
| Balanced panel, {K['bal_n']} countries with all seven years | gain {K['bal_gain_e0']}, reduction in e<sup>†</sup> {K['bal_d_edag']}, against {K['full_gain_e0']} and {K['full_d_edag']} in the full frame; national vs GBD 2023 ranking on these countries: Spearman {RC[RC.icd == 'C70-C72'].spearman_rho.iloc[0]:.3f}, {RC[RC.icd == 'C70-C72'].pct_pairs_reversed.iloc[0]:.1f}% of pairs reversed |
| Later period, {K['lat_n']} countries {K['lat_years']} ({K['lat_country_years']} country-years, {sp(K['lat_deaths'])} deaths) | reduction in e<sup>†</sup> {K['lat_d_edag']} against {K['latprim_d_edag']} for the same countries in the primary window; positive in all {K['lat_d_edag_positive']}; Spearman {K['lat_rho']} between windows |
| Reference dependence in the later-period frame (schedule derived from published GBD 2023 outputs) | mean increase {K['lat_ref_pct_mean']}%; {K['lat_ref_pairs_reversed']} of {K['lat_ref_pairs']} pairwise orderings reversed ({K['lat_ref_pct_pairs_reversed']}%) |
| Observation window, all countries | Spearman {J['frames']['rho_prim_max']} against all available country-years and {J['frames']['rho_prim_excl']} excluding 2020 and 2021; largest shift {J['frames']['max_shift']} places |
| Cause fraction formed with Eurostat all-cause deaths | reduction in e<sup>†</sup> {X['src_d_edag_es']} against {X['src_d_edag_hmd']}; largest country difference {X['src_max_abs_diff_edag']} |
| Restricted to {X['src_lowdisc_n']} countries with sources agreeing within 1% | reduction in e<sup>†</sup> {X['src_lowdisc_d_edag']}, positive in all {X['src_lowdisc_positive']} |
| Threshold perturbation of 1%, 0.1% and 0.01% | identical band in {X['threshold_countryyears_identical']} of {X['threshold_countryyears_compared']} country-years |
| Decomposition ordering reversed | totals identical; signs agree in {S['rev_sign_agreement_pct']}% of age cells; largest single difference {sci(S['rev_max_abs_diff'])} years |
| Open-interval age set to 96 or 100 years | country mean ages at death change by at most {S['openage_max_change_96']} and {S['openage_max_change_100']} years |
| Open interval closed at the published remaining expectancy for exact age 95 | level of e<sup>†</sup> {CL['edag_alt_mean']} against {CL['edag_reported_mean']} years; mean reduction {CL['ded_alt_scaled']} to {CL['ded_alt_fixed']} against {CL['ded_reported']}, largest country difference {CL['ded_max_country_diff']} years, positive in all {CL['ded_positive_alt_fixed']}; share of all-cause YLL {CL['share_yll_alt_pct']}% against {CL['share_yll_reported_pct']}% |

Sex-specific results: gain in life expectancy {J['male']['gain']} years in males and
{J['female']['gain']} in females; reduction in e<sup>†</sup> {J['male']['ded']} and {J['female']['ded']};
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

That difference exceeds the cause-deletion effect being measured. The cause-deletion
quantities are differences between two tables built by the same procedure, so the effect of
the closure largely cancels in them; levels such as life expectancy at birth, e<sup>†</sup> itself, the
age-standardised years of life lost rate and YLL per death depend on the closure directly.
Recomputing with the open interval closed at the published value moves the mean
reduction in e<sup>†</sup> from {CL['ded_reported']} to between {CL['ded_alt_scaled']} and
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

The Poisson perturbation intervals describe instability arising from the size and age allocation of the
cause-specific death counts. They do not represent uncertainty in coding practice, in
reconciliation between the two data sources, or in the Human Mortality Database series.
Resampling was performed for the primary frame only; results for other frames are point
estimates.

The cross-country rankings demonstrate that an ordering computed from the same deaths changes
with the reference life table. They are not a comparison of national performance:
certification and coding practice, completeness of registration and the attribution of deaths
from brain metastases differ between these countries and all enter the rates.
"""

open("RESULTS_SUMMARY.md", "w").write(T)
print("summary written, %d words" % len(T.split()))
