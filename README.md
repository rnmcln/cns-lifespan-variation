# Lifespan variation and premature mortality from malignant CNS neoplasms in Europe

Analysis code and derived data for a study of how mortality from malignant central nervous
system neoplasms (ICD-10 C70 to C72) affects life expectancy and lifespan variation across 28
European countries, and of how far cross-country comparisons of years of life lost depend on
the reference life table used.

Aaron Lawson McLean, Christian Senft, Colin W. McLean.

`RESULTS_SUMMARY.md` gives the design, the estimands and every headline result. The manuscript
is not included here.

---

## What the analysis does

Eurostat supplies cause-specific deaths by country, sex and five-year age group. The Human
Mortality Database supplies population exposures, all-cause deaths and period life tables. From
these the pipeline computes, for each country-year and sex:

- age-standardised mortality and years of life lost rates, 2013 European Standard Population;
- years of life lost under three reference life tables, national and two common;
- cause-deleted life tables, giving the change in life expectancy at birth and in life
  disparity when the cause is set to zero;
- threshold ages, obtained numerically;
- an exactly additive decomposition of the life-disparity effect by age;
- percentile intervals from 10 000 Poisson perturbation replicates;
- rank agreement between reference life tables, by rank correlation, absolute rank displacement
  and the number of pairwise country orderings that reverse;
- the same comparison with one fixed open-interval reference value applied to every country-year,
  so that the common reference is common at every age.

The primary analysis frame is the common period 2012 to 2018 for all 28 countries. Sensitivity
frames cover a balanced 26-country panel, the 19 countries with complete coverage to 2021, all
available country-years, and all years excluding 2020 and 2021.

## Repository layout

```
code/            pipeline, run order declared in code/gate.json
data/gbd/        Global Burden of Disease extracts and reference schedules
data/raw_hmd/    empty; see the README there for authenticated retrieval
data/derived/    every generated dataset
results/         machine-generated result files consumed by the tables and figures
figures/         figures in vector PDF and 600 dpi PNG
```

## Reproducing the results

Python 3.12.

```bash
pip install -r requirements.txt
printf 'HMD_USER=you@example.org\nHMD_PASS=...\n' > code/.hmdenv
cd code
python3 fetch_hmd.py          # retrieves the 140 Human Mortality Database files
for s in a b c d e f g h i j k m p q; do python3 stage_$s.py; done
python3 stage_n.py
python3 results_json.py
python3 figs2.py
```

`code/gate.json` declares the run order and the full list of expected outputs. A fixed seed
(20260803) is used throughout and figure PDFs are written without a timestamp, so repeated runs
on the same inputs produce byte-identical outputs.

`data/derived/` and `results/` contain the outputs of the run that produced the published
results, including the Eurostat extract, so `stage_a.py` need not be re-run against the live
API; re-running it on 15 August 2026 reproduced the committed extract byte for byte. The Human
Mortality Database series are not redistributable, so `hmd_all.csv` is absent and every stage
that reads it requires a prior `fetch_hmd.py`. Stages `d`, `g`, `k`, `n`, `results_json.py` and
the figure and summary scripts run from the shipped files alone; `b`, `c`, `e`, `f`, `h`, `i`,
`j`, `m`, `p` and `q` need the raw series.

The full pipeline has been run end to end from the raw series and verified: all declared
outputs are produced, are identical between two independent runs, and match the committed
versions byte for byte. `stage_q.py` was added after that gate run and has been checked
separately: it reproduces `hmd_open_interval_deaths.csv` byte for byte from the raw series.

## Data sources and licensing

| Source | Contents | Retrieved | Redistributable |
|---|---|---|---|
| Eurostat `hlth_cd_aro` | deaths by cause, residence, sex and five-year age group | 2 Aug 2026 | yes, open data |
| Human Mortality Database | exposures, period life tables, all-cause deaths, 5x1 | 2 Aug 2026 | **no**, see `data/raw_hmd/README.md` |
| GBD 2023 Demographics | theoretical minimum risk life table, doi:10.6069/6HTR-RB51 | 3 Aug 2026 | yes, IHME free-of-charge non-commercial terms |
| GBD 2023 Results | global CNS deaths and years of life lost by age | 2 Aug 2026 | yes, same terms |

Eurostat and the GBD tools are live services. Re-running retrieval later may return revised
source data and therefore different numbers. The retrieved extracts are included in
`data/derived/` and `data/gbd/` so that results can be reproduced rather than merely re-run.

## Two things worth knowing before modifying the code

**Observed and cause-deleted life tables must be built by the same function.** Taking an
observed value of life disparity from a published table while computing the cause-deleted value
from a reconstruction introduces differences attributable to the method of construction rather
than to removing the cause, and those differences can exceed the effect being measured.
`engine.build_lt` is used for both.

**The age decomposition is path dependent.** It is exactly additive by construction but the
allocation across age groups depends on the order of removal. Reported contributions are the
symmetric average of the forward and reverse orderings; both are retained in
`data/derived/contrib_symmetric.csv`.

## Licence

Code is released under the MIT licence (`LICENSE`). Derived datasets in `data/derived/` and
`results/` are released under CC BY 4.0, except where they are computed from Human Mortality
Database series, which remain subject to that database's terms.

## Citation

See `CITATION.cff`. Please cite the archived release rather than the repository default branch.
