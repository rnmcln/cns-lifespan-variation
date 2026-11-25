# Results summary

Reference life tables and country rankings of years of life lost: malignant central nervous
system neoplasms and eight comparison causes in European countries.

Aaron Lawson McLean, Christian Senft, Colin W. McLean.

This file records the design and every headline result produced by the code in this
repository. It stands in place of the manuscript, which is not distributed here. All figures
below are read directly from the analysis outputs in `results/` and `data/derived/`.

---

## Question

Years of life lost (YLL) depend on a reference life table: national period life tables or a
common standard such as the Global Burden of Disease (GBD) theoretical minimum risk life table.
The question is how far this choice changes the ranking of countries for the same cause of death
when the recorded deaths are held fixed: how often two countries swap places, how far apart their
rates are when they do, and whether the swaps are stable under random variation in the death
counts. The main case is malignant central nervous system neoplasms (ICD-10 C70 to C72) in
28 European countries; eight further causes show where rankings are most
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

Primary frame: the common period 2012 to 2018 for all 28 countries,
comprising 190 of 196 potentially eligible
both-sex country-years, 242 586 C70 to C72 deaths and
34 651 958 all-cause deaths. Six country-years fail the reconciliation
rule, all in Iceland and Luxembourg. Within the frame, Eurostat disclosure control leaves
46 both-sex age cells without a value, all in Iceland and Luxembourg; the
sex-specific extracts have 116 such cells in males and 116 in females.

Countries are ranked on the age-standardised YLL rate (2013 European Standard Population) under
three references: national period life tables, the official GBD 2023 theoretical minimum risk life
table, and an empirical minimum-mortality table built from the lowest age-specific rates in the
frame. For every pair of countries the gap between their rates is recorded under each reference,
as a percentage of the pair mean. Stability is assessed with 2 000 replicates in
which every age-specific death count is resampled from a Poisson distribution and both rankings
are recomputed from the same replicate deaths (seed 20260803); a reversal recurring in at
least 95% of replicates is called stable. These are perturbation summaries, not confidence
intervals.

## Reference dependence of years of life lost

Pooled across the frame, C70 to C72 accounts for 0.70% of all deaths and
1.23% of all-cause years of life lost under national life tables.

| Quantity | Value |
|---|---|
| Age-standardised YLL rate, national reference | 126 to 195 per 100 000 |
| Age-standardised YLL rate, GBD 2023 reference | 167 to 269 per 100 000 |
| Increase in YLL per death, GBD 2023 reference | 35.3% mean (23.1% to 57.8%) |
| Increase in YLL per death, empirical minimum-mortality reference | 13.4% mean |
| Pairwise country orderings reversed, national vs GBD 2023 | 64 of 378 (16.9%) |
| Spearman correlation, national vs GBD 2023 | 0.858 |
| Spearman correlation between the two common references | 0.997 |

Life expectancy at birth is 89.96 years in the official GBD 2023 schedule and
84.76 years in the empirical minimum-mortality schedule. Ranks were computed from
unrounded rates throughout.

## Rank agreement between reference life tables

| Comparison | Spearman ρ | Kendall τ<sub>b</sub> | Median abs. shift | Max shift | Pairs reversed |
|---|---|---|---|---|---|
| national vs GBD 2023 | 0.858 | 0.661 | 2.5 | 12 | 64 of 378 (16.9%) |
| national vs minimum-mortality | 0.857 | 0.661 | 2.5 | 12 | 64 of 378 (16.9%) |
| GBD 2023 vs minimum-mortality | 0.997 | 0.979 | 0.0 | 2 | 4 of 378 (1.1%) |
| national vs GBD 2023 (derived schedule), later frame | 0.835 | 0.626 | 2.0 | 7 | 32 of 171 (18.7%) |

## Size and stability of the reversals

Gaps are differences between the two countries of a pair as a percentage of the pair mean; the
first reference is the national life table.

| Comparison | Spearman ρ | Pairs reversed | Median gap of reversed pairs, first reference % | Median gap, second reference % | Largest first-reference gap reversed % | Reversed with both gaps ≥5% | Reversed with both gaps ≥10% |
|---|---|---|---|---|---|---|---|
| Age-standardised rate, national vs GBD 2023 | 0.858 | 64 of 378 (16.9%) | 3.3 | 5.3 | 15.7 | 14 | 0 |
| Age-standardised rate, national vs minimum-mortality | 0.857 | 64 of 378 (16.9%) | 3.3 | 5.4 | 15.7 | 13 | 0 |
| YLL per death, national vs GBD 2023 | 0.36 | 137 of 378 (36.2%) | 5.3 | 3.6 | 23.3 | 18 | 0 |

For the primary comparison, the median gap is 11.9% over all pairs and
14.8% over pairs whose ordering is kept, against
3.3% (interquartile range 1.6% to
5.9%) over reversed pairs. None of the 150
pairs further apart than 15.7% is reversed.

Across replicates the number of reversed pairs ranges from 47 to
69 (median 57). The national ordering recurs in at
least 95% of replicates for 268 of 378 pairs but for only
19 of the 64 reversed pairs.
14 reversals are stable; in 13 of them a
Central or Eastern European country moves above a country outside that group, and in
0 the reverse. The five largest stable reversals:

| Higher under national tables | Higher under GBD 2023 | Gap national % | Gap GBD 2023 % | Recurrence % |
|---|---|---|---|---|
| Switzerland | Hungary | 8.9 | 10.2 | 99.9 |
| Spain | Hungary | 8.8 | 9.8 | 100.0 |
| France | Hungary | 7.0 | 13.6 | 100.0 |
| Finland | Hungary | 6.9 | 9.1 | 98.7 |
| Switzerland | Czechia | 6.2 | 7.8 | 98.9 |

All pairs, with their gaps and recurrence, are in `data/derived/reversal_pairs.csv`.

## Rank sensitivity across nine causes

Common frame of 26 countries (IS, LU
excluded, having no retained perinatal stratum). Causes ordered by mean age at death. Ratio of
SDs: between-country SD of log(GBD rate / national rate) divided by the between-country SD of
log(national rate).

| Cause | Mean age at death | Spearman ρ | Pairs reversed % | Replicate range % | Stable reversals | SD log national rate | SD log reference effect | Ratio of SDs |
|---|---|---|---|---|---|---|---|---|
| Perinatal conditions | 1.1 | 0.985 | 3.4 | 1.8 to 5.5 | 0 | 0.286 | 0.035 | 0.12 |
| Transport accidents | 49.9 | 0.992 | 2.2 | 1.5 to 4.3 | 2 | 0.397 | 0.057 | 0.14 |
| Intentional self-harm | 54.2 | 0.978 | 4.9 | 3.7 to 6.2 | 4 | 0.385 | 0.068 | 0.18 |
| Chronic liver disease | 63.7 | 0.987 | 2.8 | 1.5 to 3.7 | 1 | 0.582 | 0.063 | 0.11 |
| Malignant CNS neoplasms | 65.0 | 0.845 | 17.8 | 13.2 to 19.1 | 13 | 0.109 | 0.071 | 0.65 |
| Lung cancer | 71.3 | 0.904 | 11.4 | 10.2 to 12.9 | 25 | 0.205 | 0.074 | 0.36 |
| Ischaemic heart disease | 80.2 | 0.995 | 1.8 | 1.5 to 2.5 | 4 | 0.584 | 0.090 | 0.15 |
| Cerebrovascular disease | 82.0 | 0.992 | 2.2 | 1.8 to 2.8 | 4 | 0.500 | 0.071 | 0.14 |
| Dementia | 87.7 | 0.992 | 2.5 | 1.8 to 3.7 | 2 | 1.489 | 0.082 | 0.06 |

Reversed pairs whose rates differ by at least 10% under both references, summed over all nine
causes: 0.

## Reference effect on YLL per death by cause

Spearman ρ between mean age at death and the mean increase in YLL per death:
0.983. Cause ordering identical to the full-coverage frame:
True.

| Cause | Countries | Deaths | Mean age at death | Increase % | Range % |
|---|---|---|---|---|---|
| Perinatal conditions | 26 | 63 553 | 1.1 | 12.9 | 8.2 to 20.8 |
| Transport accidents | 26 | 181 810 | 49.9 | 25.5 | 16.2 to 40.5 |
| Intentional self-harm | 26 | 366 401 | 54.2 | 28.3 | 18.8 to 49.1 |
| Chronic liver disease | 26 | 456 240 | 63.7 | 36.5 | 24.4 to 57.2 |
| Malignant CNS neoplasms | 26 | 242 393 | 65.0 | 35.7 | 23.1 to 57.8 |
| Lung cancer | 26 | 1 839 736 | 71.3 | 44.7 | 27.2 to 65.3 |
| Ischaemic heart disease | 26 | 3 655 522 | 80.2 | 59.0 | 40.1 to 80.9 |
| Cerebrovascular disease | 26 | 2 323 028 | 82.0 | 63.9 | 45.6 to 89.9 |
| Dementia | 26 | 1 298 738 | 87.7 | 85.1 | 59.4 to 111.0 |

## A literally common open interval

The GBD 2023 reference is common at every age below 95 years; for the open interval the four
sub-interval values are averaged with each country-year's distribution of all-cause deaths above
95. 0.28% of C70 to C72 deaths fall in that interval, and the
country-year values span 7.860 to 8.074
years. Replacing them with one fixed value for every country-year, either the schedule at exact
age 97.5 years (8.136) or the pooled death-weighted mean
(7.967), changes YLL per death by at most 0.0009 years
and the age-standardised rate by at most 0.01 per 100 000. Country
ordering identical under both: True.

## Secondary analysis: life expectancy and lifespan variation

| Quantity | Value |
|---|---|
| Gain in life expectancy at birth on cause deletion | 0.117 years, about 43 days (range 0.093 to 0.156) |
| Reduction in life disparity e<sup>†</sup> | 0.0575 years, about 21 days (range 0.0401 to 0.0757) |
| Observed e<sup>†</sup> for comparison | about 10.39 years |
| Countries with a positive reduction in e<sup>†</sup> | 28 of 28 |
| Countries whose 95% Poisson perturbation interval lies wholly above zero | 28 of 28 |
| Median threshold boundary a<sup>†</sup> | 80.0 years |
| Share of C70 to C72 deaths below a<sup>†</sup> | 68.0% to 90.0% |
| Mean age at C70 to C72 death | 62.1 to 67.5 years |
| Age-standardised mortality rate | 5.91 to 10.09 per 100 000 person-years |

Age-specific contributions to the reduction in e<sup>†</sup> are positive from childhood to about age 75,
peak between 50 and 65 years, and are slightly negative above age 80. The positive
contributions dominate in every country, which is why the net effect is compressive.

## Robustness

| Analysis | Result |
|---|---|
| Balanced panel, 26 countries with all seven years | gain 0.117, reduction in e<sup>†</sup> 0.0573, against 0.117 and 0.0575 in the full frame; national vs GBD 2023 ranking on these countries: Spearman 0.845, 17.8% of pairs reversed |
| Later period, 19 countries 2012-2021 (190 country-years, 182 281 deaths) | reduction in e<sup>†</sup> 0.0568 against 0.0576 for the same countries in the primary window; positive in all 19; Spearman 0.9667 between windows |
| Reference dependence in the later-period frame (schedule derived from published GBD 2023 outputs) | mean increase 35.5%; 32 of 171 pairwise orderings reversed (18.7%) |
| Observation window, all countries | Spearman 0.977 against all available country-years and 0.982 excluding 2020 and 2021; largest shift 5 places |
| Cause fraction formed with Eurostat all-cause deaths | reduction in e<sup>†</sup> 0.0588 against 0.0575; largest country difference 0.0056 |
| Restricted to 18 countries with sources agreeing within 1% | reduction in e<sup>†</sup> 0.0573, positive in all 18 |
| Threshold perturbation of 1%, 0.1% and 0.01% | identical band in 190 of 190 country-years |
| Decomposition ordering reversed | totals identical; signs agree in 99.97% of age cells; largest single difference 9.7 × 10<sup>−5</sup> years |
| Open-interval age set to 96 or 100 years | country mean ages at death change by at most 0.0083 and 0.0138 years |
| Open interval closed at the published remaining expectancy for exact age 95 | level of e<sup>†</sup> 10.334 against 10.386 years; mean reduction 0.0579 to 0.0582 against 0.0575, largest country difference 0.00143 years, positive in all 28; share of all-cause YLL 1.231% against 1.227% |

Sex-specific results: gain in life expectancy 0.122 years in males and
0.105 in females; reduction in e<sup>†</sup> 0.0541 and 0.0573;
mean age at death 62.8 and 66.2 years; median threshold boundary
75.0 and 80.0 years. Sex-specific strata are more heavily
affected by Eurostat disclosure control, in the extreme leaving a single retained year for one
country in each sex, so country-level sex-specific estimates should be read with that in mind.

Agreement between the two data sources on all-cause totals across retained country-years:
median absolute difference 0.458%, 95th percentile
2.368%, maximum 5.19%.

## Sensitivity to the closure of the open age interval

The reconstruction collapses ages above 95 into one open interval and closes it with L = l / m.
The published Human Mortality Database series extends to 110 years and closes there, so
remaining expectancy at exact age 95 differs between the two, and because T(x) is cumulated
downward that difference reaches every younger age. Life expectancy at birth in the
reconstruction exceeds the published value by a mean of 0.015
years and at most 0.0535 years, and life disparity by
0.052 years.

That difference exceeds the cause-deletion effect being measured. The cause-deletion
quantities are differences between two tables built by the same procedure, so the effect of
the closure largely cancels in them; levels such as life expectancy at birth, e<sup>†</sup> itself, the
age-standardised years of life lost rate and YLL per death depend on the closure directly.
Recomputing with the open interval closed at the published value moves the mean
reduction in e<sup>†</sup> from 0.0575 to between 0.0579 and
0.0582 years, with a largest country difference of 0.00143
years, and the reduction stays positive in all 28 countries. The
gain in life expectancy moves by at most 0.00048 years. The share of
all-cause years of life lost moves from 1.227% to
1.231%, and the age-standardised years of life lost rate by at most
0.24% in any country, changing the rank of
2 of 28 countries by one place.

0.28% of C70 to C72 deaths fall in the open interval, against
5.4% of all-cause deaths, which is why the cause-specific results
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
