# Reference life tables and country rankings of years of life lost: malignant CNS neoplasms and eight comparison causes in Europe

Analysis code and derived data for a study of how far the choice of reference life table changes
the ranking of countries by years of life lost when the recorded deaths are held fixed. The main
case is malignant central nervous system neoplasms (ICD-10 C70 to C72) in 28 European countries;
eight further causes show where rankings are most and least sensitive to the reference. A
secondary analysis estimates how these deaths affect life expectancy and lifespan variation.

Aaron Lawson McLean, Christian Senft, Colin W. McLean.

`RESULTS_SUMMARY.md` gives the design, the estimands and every headline result. The manuscript
is not included here.

---

## What the analysis does

Eurostat supplies cause-specific deaths by country, sex and five-year age group. The Human
Mortality Database supplies population exposures, all-cause deaths and period life tables. From
these the pipeline computes, for each country-year and sex:

- age-standardised mortality and years of life lost rates, 2013 European Standard Population;
- years of life lost under three reference life tables, national and two common, and rank
  agreement between the resulting country orderings;
- cause-deleted life tables, giving the change in life expectancy at birth and in life
  disparity when the cause is set to zero;
- threshold ages, obtained numerically;
- an exactly additive decomposition of the life-disparity effect by age;
- 95% Poisson perturbation intervals from 10 000 replicates, describing the stability of the
  estimates under perturbation of the cause-specific counts;
- rank agreement between reference life tables, by rank correlation, absolute rank displacement
  and the number of pairwise country orderings that reverse;
- for every pair of countries, the gap between their rates under each reference, and the
  proportion of 2 000 Poisson replicates of the deaths in which each reversal recurs;
- rank agreement, reversal sizes and replicate recurrence for all nine causes on a common country
  frame, with the ratio of the between-country spread of the reference effect to the spread of
  the rates;
- the same comparison with one fixed open-interval reference value applied to every country-year,
  so that the common reference is common at every age;
- the same comparison for eight further causes spanning the range of ages at death;
- counts of Eurostat confidentiality-flagged age cells within the analysed frame.

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

# every script reads and writes in one flat working directory
mkdir work && cd work
cp ../code/*.py .
cp -r ../data/gbd gbd
cp ../data/derived/hmd_vintages.csv .
printf 'HMD_USER=you@example.org\nHMD_PASS=...\n' > .hmdenv

python3 fetch_hmd.py          # retrieves the 140 Human Mortality Database files into raw/
for s in a b c d e f g h i j k m p q r s; do python3 stage_$s.py; done
python3 results_json.py
python3 stage_n.py
python3 figs2.py
python3 make_summary.py
```

The outputs are written into `work/`; their committed copies are in `results/`, `data/derived/`,
`figures/` and the repository root. To re-run only the stages that do not need the raw series,
copy the contents of `data/derived/` and `results/` into the working directory first.

`code/gate.json` declares the run order and the full list of expected outputs. A fixed seed
(20260803) is used throughout and figure PDFs are written without a timestamp, so repeated runs
on the same inputs produce byte-identical outputs.

`data/derived/` and `results/` contain the outputs of the run that produced the reported
results, including the Eurostat extract. The Human Mortality Database series are not
redistributable, so `hmd_all.csv` is absent and every stage that reads it requires a prior
`fetch_hmd.py`. Stages `d`, `k`, `n` and `r` and the figure and summary scripts run from the
shipped files alone; `b`, `c`, `e`, `f`, `g`, `h`, `i`, `j`, `m`, `p`, `q` and `s` need the raw series,
and `results_json.py` reads `decomposition.csv`, which stage `e` writes from them.

On 25 November 2026 the full pipeline was run twice from a clean working directory containing
only the code, the GBD extracts, `hmd_vintages.csv` and the 140 raw series, whose SHA-256 digests
matched the manifest. Stages `a` and `j` retrieved the Eurostat extracts from the live API. Both runs
completed, produced all 49 declared outputs, gave identical outputs, and matched the committed
versions byte for byte. Figure output depends on the fonts installed (Arial, or a
metric-compatible substitute such as Liberation Sans).

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
