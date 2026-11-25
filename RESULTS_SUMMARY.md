# Results summary

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

Primary frame: the common period 2012 to 2018 for all 28 countries,
comprising 190 of 196 potentially eligible
both-sex country-years, 242 586 C70 to C72 deaths and
34 651 958 all-cause deaths. Six country-years fail the reconciliation
rule, all in Iceland and Luxembourg. Within the frame, Eurostat disclosure control leaves
46 both-sex age cells without a value, all in Iceland and Luxembourg; the
sex-specific extracts have 116 such cells in males and 116 in females.

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

Pooled across the frame, C70 to C72 accounts for 0.70% of all deaths and
1.23% of all-cause years of life lost under national life tables, a ratio of
1.75 that exceeds 1.4 in every country (range 1.46 to
2.04).

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

Across nine causes spanning mean ages at death from 1.1 to 87.7
years and 10.4 million deaths, on the common frame of the
26 countries available for every cause, the reference-induced increase
rises with the mean age at death, from 12.9% for
perinatal conditions to 85.1% for
dementia (Spearman ρ = 0.983, with one
inversion).

## Rank agreement between reference life tables

| Comparison | Spearman ρ | Kendall τb | Median abs. shift | Max shift | Pairs reversed |
|---|---|---|---|---|---|
| national vs GBD 2023 | 0.858 | 0.661 | 2.5 | 12 | 64 of 378 (16.9%) |
| national vs minimum-mortality | 0.857 | 0.661 | 2.5 | 12 | 64 of 378 (16.9%) |
| GBD 2023 vs minimum-mortality | 0.997 | 0.979 | 0.0 | 2 | 4 of 378 (1.1%) |
| national vs GBD 2023 (derived schedule), later frame | 0.835 | 0.626 | 2.0 | 7 | 32 of 171 (18.7%) |

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

## Reference effect by cause

Common frame of 26 countries (IS, LU excluded,
having no retained perinatal stratum). Spearman ρ between mean age at death and the mean
increase: 0.983. Cause ordering identical to the full-coverage frame:
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

## Life expectancy and lifespan variation (secondary analysis)

| Quantity | Value |
|---|---|
| Gain in life expectancy at birth on cause deletion | 0.117 years, about 43 days (range 0.093 to 0.156) |
| Reduction in life disparity e† | 0.0575 years, about 21 days (range 0.0401 to 0.0757) |
| Observed e† for comparison | about 10.39 years |
| Countries with a positive reduction in e† | 28 of 28 |
| Countries whose 95% Poisson perturbation interval lies wholly above zero | 28 of 28 |
| Median threshold boundary a† | 80.0 years |
| Share of C70 to C72 deaths below a† | 68.0% to 90.0% |
| Mean age at C70 to C72 death | 62.1 to 67.5 years |
| Age-standardised mortality rate | 5.91 to 10.09 per 100 000 person-years |

Age-specific contributions to the reduction in e† are positive from childhood to about age 75,
peak between 50 and 65 years, and are slightly negative above age 80. The positive
contributions dominate in every country, which is why the net effect is compressive.

## Robustness

| Analysis | Result |
|---|---|
| Balanced panel, 26 countries with all seven years | gain 0.117, reduction in e† 0.0573, against 0.117 and 0.0575 in the full frame; ordering unchanged |
| Later period, 19 countries 2012-2021 (190 country-years, 182 281 deaths) | reduction in e† 0.0568 against 0.0576 for the same countries in the primary window; positive in all 19; Spearman 0.9667 between windows |
| Reference dependence in the later-period frame (schedule derived from published GBD 2023 outputs) | mean increase 35.5%; 32 of 171 pairwise orderings reversed (18.7%) |
| Observation window, all countries | Spearman 0.977 against all available country-years and 0.982 excluding 2020 and 2021; largest shift 5 places |
| Cause fraction formed with Eurostat all-cause deaths | reduction in e† 0.0588 against 0.0575; largest country difference 0.0056 |
| Restricted to 18 countries with sources agreeing within 1% | reduction in e† 0.0573, positive in all 18 |
| Threshold perturbation of 1%, 0.1% and 0.01% | identical band in 190 of 190 country-years |
| Decomposition ordering reversed | totals identical; signs agree in 99.97% of age cells; largest single difference 9.7e-05 years |
| Open-interval age set to 96 or 100 years | country mean ages at death change by at most 0.0083 and 0.0138 years |
| Open interval closed at the published remaining expectancy for exact age 95 | level of e† 10.334 against 10.386 years; mean reduction 0.0579 to 0.0582 against 0.0575, largest country difference 0.00143 years, positive in all 28; share of all-cause YLL 1.231% against 1.227% |

Sex-specific results: gain in life expectancy 0.122 years in males and
0.105 in females; reduction in e† 0.0541 and 0.0573;
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
the closure largely cancels in them; levels such as life expectancy at birth, e† itself, the
age-standardised years of life lost rate and YLL per death depend on the closure directly.
Recomputing with the open interval closed at the published value moves the mean
reduction in e† from 0.0575 to between 0.0579 and
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
