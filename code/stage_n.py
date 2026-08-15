"""Stage N: supplementary datasets.

1. Adds the European empirical minimum-mortality schedule to the reference-schedule file,
   so that every reference used in the analysis is supplied.
2. Writes a data dictionary for the country-year results file.
"""
import pandas as pd, json

J = json.load(open("RESULTS.json"))

# ---------------------------------------------- 1. reference schedules
R = pd.read_csv("reference_schedules.csv")
bp = J["sens"]["bp_schedule"]
R["minimum_mortality_at_exact_age"] = R.age_group.map(bp)
assert R.minimum_mortality_at_exact_age.notna().all()
R.to_csv("reference_schedules.csv", index=False)

# ---------------------------------------------- 2. data dictionary
CY = pd.read_csv("results_countryyear.csv")
DEF = {
 "geo": ("Eurostat two-letter country code", "code"),
 "sex": ("T both sexes, M males, F females", "code"),
 "year": ("calendar year", "year"),
 "deaths": ("C70 to C72 deaths in the stratum", "count"),
 "dall": ("all-cause deaths, Human Mortality Database", "count"),
 "expo": ("person-years of exposure, Human Mortality Database", "person-years"),
 "n_suppressed": ("age cells carrying a Eurostat confidentiality flag", "count"),
 "asmr": ("age-standardised mortality rate, 2013 European Standard Population",
          "deaths per 100 000 person-years"),
 "yll_nat": ("years of life lost under the national period life table", "years"),
 "yll_gbd": ("years of life lost under the GBD 2023 reference", "years"),
 "asyr_nat": ("age-standardised years of life lost rate, national reference",
              "years per 100 000 person-years"),
 "asyr_gbd": ("age-standardised years of life lost rate, GBD 2023 reference",
              "years per 100 000 person-years"),
 "e0": ("observed period life expectancy at birth", "years"),
 "e0_del": ("life expectancy at birth in the cause-deleted life table", "years"),
 "gain_e0": ("e0_del minus e0", "years"),
 "edag": ("observed life disparity", "years"),
 "edag_del": ("life disparity in the cause-deleted life table", "years"),
 "d_edag": ("edag minus edag_del; positive values denote compression", "years"),
 "thresh": ("threshold age boundary, lower bound of the identified band", "years"),
 "mad": ("mean age at C70 to C72 death from grouped counts", "years"),
 "frac_below_thresh": ("share of C70 to C72 deaths below the threshold boundary",
                       "proportion"),
 "stepwise_total": ("sum of the age-specific contributions to d_edag; equals d_edag "
                    "to machine precision", "years"),
 "Rx_min": ("smallest value of the cause-deletion multiplier across age groups",
            "ratio"),
}
rows = [{"variable": c,
         "description": DEF.get(c, ("", ""))[0],
         "unit": DEF.get(c, ("", ""))[1]} for c in CY.columns]
D = pd.DataFrame(rows)
assert (D.description != "").all(), D[D.description == ""]
D.to_csv("data_dictionary_countryyear.csv", index=False)

print("reference schedules: %d rows, %d columns" % R.shape)
print("data dictionary: %d variables" % len(D))
