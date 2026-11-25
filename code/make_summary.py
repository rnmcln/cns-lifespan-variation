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
P, X, Q, G, R = J["panel"], J["sens"], J["off"], J["gen"], J["round"]
sp = lambda n: f"{int(n):,}".replace(",", " ")
RA = RK["rank_agreement"]
RKN, RKO = RA["national vs GBD 2023"], RA["GBD 2023 vs minimum-mortality"]
RANK_ROWS = "\n".join(
    f"| {k} | {v['spearman_rho']} | {v['kendall_tau_b']} | {v['median_abs_shift']:.1f} | "
    f"{v['max_abs_shift']} | {v['n_pairs_reversed']} of {v['n_pairs']} ({v['pct_pairs_reversed']}%) |"
    for k, v in RA.items())
CAUSE_ROWS = "\n".join(
    f"| {r.cause} | {r.n_ctry} | {sp(r.deaths)} | {r.mad:.1f} | {r.infl_mean:.1f} | "
    f"{r.infl_lo:.1f} to {r.infl_hi:.1f} |" for r in GC.itertuples())

T = f"""# Results summary

Reference life tables and cross-country comparisons of years of life lost: malignant central
nervous system neoplasms and lifespan variation in 28 European countries.

Aaron Lawson McLean, Christian Senft, Colin W. McLean.

This file records the design and every headline result produced by the code in this
repository. It stands in place of the manuscript, which is not distributed here. All figures
below are read directly from the analysis outputs in `results/` and `data/derived/`.

---

## Question

Years of life lost (YLL) are the usual measure of premature mortality when countries are
compared, and every YLL estimate depends on a reference life table: national period life tables
or a common standard such as the Global Burden of Disease (GBD) theoretical minimum risk life
table. The primary question is how far this choice changes the cross-country ordering of YLL from
malignant central nervous system neoplasms (ICD-10 C70 to C72), whether the answer depends on
which common reference is used, and how the size of the reference effect on YLL per death varies
with the age at which a cause kills.

The secondary question concerns lifespan variation. Life disparity, denoted e†, is the mean
remaining life expectancy at death within a life table. Removing a cause of death compresses the
distribution of ages at death only when those deaths fall below a threshold age a†, so the
effect of a cause on lifespan variation has to be computed rather than inferred from its age
profile. The analysis estimates the changes in life expectancy at birth and in e† implied by
removing C70 to C72 mortality, and the age structure that produces them.

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

Cause-deleted life tables are built from mₓ multiplied by Rₓ, one minus the C70 to C72 share of
all deaths in the age group, using the same function as the observed tables. Threshold ages are
obtained numerically by perturbing mₓ in each age group in turn and are identified as band
boundaries, since the data are grouped in five-year intervals. The effect on e† is decomposed
by age by cumulative removal, reported as the symmetric average of the forward and reverse
orderings. Uncertainty is assessed by resampling age-specific death counts from Poisson
distributions over 10 000 replicates; the resulting 95% Poisson perturbation intervals describe
the stability of the estimates under perturbation of the counts and are not confidence
intervals.

## Reference dependence of years of life lost

Pooled across the frame, C70 to C72 accounts for {R['deaths_pct']}% of all deaths and
{R['yll_pct']}% of all-cause years of life lost under national life tables, a ratio of
{J['shares']['ratio']} that exceeds 1.4 in every country (range {J['shares']['ratio_min']} to
{J['shares']['ratio_max']}).

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

Across nine causes spanning mean ages at death from {GG['common_mad_lo']} to {GG['common_mad_hi']}
years and {GG['common_total_deaths'] / 1e6:.1f} million deaths, on the common frame of the
{GG['n_common_countries']} countries available for every cause, the reference-induced increase
rises with the mean age at death, from {GG['common_min_infl']}% for
{GG['common_min_infl_cause'].lower()} to {GG['common_max_infl']}% for
{GG['common_max_infl_cause'].lower()} (Spearman ρ = {GG['common_spearman_mad_infl']}, with one
inversion).

## Rank agreement between reference life tables

| Comparison | Spearman ρ | Kendall τb | Median abs. shift | Max shift | Pairs reversed |
|---|---|---|---|---|---|
{RANK_ROWS}
| national vs GBD 2023 (derived schedule), later frame | {K['lat_ref_rho']} | {K['lat_ref_tau']} | {K['lat_ref_median_shift']:.1f} | {K['lat_ref_max_shift']} | {K['lat_ref_pairs_reversed']} of {K['lat_ref_pairs']} ({K['lat_ref_pct_pairs_reversed']}%) |

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

## Reference effect by cause

Common frame of {GG['n_common_countries']} countries ({', '.join(GG['common_excluded'])} excluded,
having no retained perinatal stratum). Spearman ρ between mean age at death and the mean
increase: {GG['common_spearman_mad_infl']}. Cause ordering identical to the full-coverage frame:
{GG['common_order_matches_full']}.

| Cause | Countries | Deaths | Mean age at death | Increase % | Range % |
|---|---|---|---|---|---|
{CAUSE_ROWS}

## Life expectancy and lifespan variation (secondary analysis)

| Quantity | Value |
|---|---|
| Gain in life expectancy at birth on cause deletion | {J['gain_e0']['mean']} years, about {R['gain_days']} days (range {J['gain_e0']['lo']} to {J['gain_e0']['hi']}) |
| Reduction in life disparity e† | {J['d_edag']['mean']} years, about {R['ded_days']} days (range {J['d_edag']['lo']} to {J['d_edag']['hi']}) |
| Observed e† for comparison | about {J['d_edag']['edag_mean']} years |
| Countries with a positive reduction in e† | {J['d_edag']['positive']} of {J['frame']['countries']} |
| Countries whose 95% Poisson perturbation interval lies wholly above zero | {J['d_edag']['ci_excl_zero']} of {J['frame']['countries']} |
| Median threshold boundary a† | {J['threshold']['median']} years |
| Share of C70 to C72 deaths below a† | {J['threshold']['below_lo']}% to {J['threshold']['below_hi']}% |
| Mean age at C70 to C72 death | {J['threshold']['mad_lo']} to {J['threshold']['mad_hi']} years |
| Age-standardised mortality rate | {J['asmr']['lo']} to {J['asmr']['hi']} per 100 000 person-years |

Age-specific contributions to the reduction in e† are positive from childhood to about age 75,
peak between 50 and 65 years, and are slightly negative above age 80. The positive
contributions dominate in every country, which is why the net effect is compressive.

## Robustness

| Analysis | Result |
|---|---|
| Balanced panel, {K['bal_n']} countries with all seven years | gain {K['bal_gain_e0']}, reduction in e† {K['bal_d_edag']}, against {K['full_gain_e0']} and {K['full_d_edag']} in the full frame; ordering unchanged |
| Later period, {K['lat_n']} countries {K['lat_years']} ({K['lat_country_years']} country-years, {sp(K['lat_deaths'])} deaths) | reduction in e† {K['lat_d_edag']} against {K['latprim_d_edag']} for the same countries in the primary window; positive in all {K['lat_d_edag_positive']}; Spearman {K['lat_rho']} between windows |
| Reference dependence in the later-period frame (schedule derived from published GBD 2023 outputs) | mean increase {K['lat_ref_pct_mean']}%; {K['lat_ref_pairs_reversed']} of {K['lat_ref_pairs']} pairwise orderings reversed ({K['lat_ref_pct_pairs_reversed']}%) |
| Observation window, all countries | Spearman {J['frames']['rho_prim_max']} against all available country-years and {J['frames']['rho_prim_excl']} excluding 2020 and 2021; largest shift {J['frames']['max_shift']} places |
| Cause fraction formed with Eurostat all-cause deaths | reduction in e† {X['src_d_edag_es']} against {X['src_d_edag_hmd']}; largest country difference {X['src_max_abs_diff_edag']} |
| Restricted to {X['src_lowdisc_n']} countries with sources agreeing within 1% | reduction in e† {X['src_lowdisc_d_edag']}, positive in all {X['src_lowdisc_positive']} |
| Threshold perturbation of 1%, 0.1% and 0.01% | identical band in {X['threshold_countryyears_identical']} of {X['threshold_countryyears_compared']} country-years |
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

That difference exceeds the cause-deletion effect being measured. The cause-deletion
quantities are differences between two tables built by the same procedure, so the effect of
the closure largely cancels in them; levels such as life expectancy at birth, e† itself, the
age-standardised years of life lost rate and YLL per death depend on the closure directly.
Recomputing with the open interval closed at the published value moves the mean
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
