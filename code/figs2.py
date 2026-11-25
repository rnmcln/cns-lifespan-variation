"""Figures 1 to 3 and the two supplementary figures.

Design: minimal ink, sans-serif type, thin rules, no gridlines, no enclosing boxes, no
titles inside the plotting area, greyscale with a single muted accent. Vector PDF is the
primary output; PNG is written at 600 dpi for convenience. Figure widths are set for a
170 mm full-page column with type that remains legible on reduction.
"""
import pandas as pd, numpy as np, matplotlib, json
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from engine import KEYS

INK = "#1A1A1A"      # text and point estimates
RULE = "#8C8C8C"     # axis lines and interval whiskers
FAINT = "#D9D9D9"    # connectors and bands
ACC = "#A03530"      # single muted accent

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
    "font.size": 7, "axes.labelsize": 7, "xtick.labelsize": 6.6, "ytick.labelsize": 6.6,
    "axes.edgecolor": RULE, "axes.linewidth": 0.5,
    "xtick.color": RULE, "ytick.color": RULE,
    "xtick.major.width": 0.5, "ytick.major.width": 0.5,
    "xtick.major.size": 2.5, "ytick.major.size": 2.5,
    "text.color": INK, "axes.labelcolor": INK,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": False, "legend.frameon": False,
    "figure.dpi": 600, "savefig.bbox": "tight", "savefig.pad_inches": 0.02,
    "pdf.fonttype": 42, "ps.fonttype": 42,
})

NAME = {"AT":"Austria","BE":"Belgium","BG":"Bulgaria","HR":"Croatia","CZ":"Czechia",
"DK":"Denmark","EE":"Estonia","FI":"Finland","FR":"France","DE":"Germany","EL":"Greece",
"HU":"Hungary","IS":"Iceland","IE":"Ireland","IT":"Italy","LV":"Latvia","LT":"Lithuania",
"LU":"Luxembourg","NL":"Netherlands","NO":"Norway","PL":"Poland","PT":"Portugal",
"SK":"Slovakia","SI":"Slovenia","ES":"Spain","SE":"Sweden","CH":"Switzerland",
"UK":"United Kingdom"}

S = pd.read_csv("country_summary_all.csv")
P = S[(S.frame == "primary_2012_2018") & (S.sex == "T")].set_index("geo")
P = P.join(pd.read_csv("official_reference_country.csv").set_index("geo")[["asyr_off", "ypd_off"]])
U = pd.read_csv("uncertainty_country.csv").set_index("geo")
CS = pd.read_csv("contrib_symmetric.csv")
J = json.load(open("RESULTS.json"))

MM = 1 / 25.4



def save(fig, stem, extra=None, tight=True):
    # matplotlib stamps a wall-clock /CreationDate into PDF output, which would make the
    # figures differ byte for byte between otherwise identical runs. Fixing it keeps the
    # declared outputs reproducible.
    # bbox_inches="tight" crops to drawn content, and a label placed over the blank part
    # of an axes with hidden spines is not part of that content, so text set outside the
    # plotted data would be cropped unless it is declared here.
    kw = {"bbox_extra_artists": extra} if extra else {}
    if not tight:
        # a figure whose labels are positioned inside the axes is emitted at its declared
        # size; cropping to drawn content would cut a label sitting over the blank part of
        # an axes with hidden spines, which is not part of that content
        kw = {"bbox_inches": None}
    fig.savefig(stem + ".pdf", metadata={"CreationDate": None}, **kw)
    fig.savefig(stem + ".png", **kw)
    plt.close(fig)


# --------------------------------------------------------- Figure 1: life disparity
d = U.join(P[["below"]]).sort_values("ded")
xmax = float(np.ceil(d.dhi.max() * 1000) / 1000) + 0.004
fig, ax = plt.subplots(figsize=(140 * MM, 150 * MM))
y = np.arange(len(d))
ax.axvline(0, color=FAINT, lw=0.5, zorder=1)
ax.hlines(y, d.dlo, d.dhi, color=RULE, lw=0.7, zorder=2)
ax.scatter(d.ded, y, s=9, facecolors="white", edgecolors=INK, linewidths=0.7, zorder=3)
lab_x = -0.0035
for i, g in enumerate(d.index):
    ax.text(lab_x, i, NAME[g], ha="right", va="center", fontsize=6.6, color=INK)
    ax.text(xmax, i, f"{100 * d.below[g]:.0f}", ha="right", va="center", fontsize=6.6, color=INK)
ax.text(xmax, len(d) - 0.1, "% below $a^{\\dagger}$", ha="right", va="center", fontsize=6.6, color=RULE)
ax.set_yticks([])
ax.set_ylim(-1, len(d) + 0.7)
ax.set_xlim(-0.005, xmax)
ax.set_xticks([0, 0.02, 0.04, 0.06, 0.08, 0.10, 0.12])
ax.set_xlabel("reduction in life disparity $e^{\\dagger}$ (years)")
for s in ["left", "top", "right"]:
    ax.spines[s].set_visible(False)
save(fig, "fig1_life_disparity")
assert d.dhi.max() < xmax, "interval clipped"

# ------------------------------------------------- Figure 2: age-specific contributions
piv = CS.pivot_table(index="geo", columns="k", values="sym", aggfunc="mean").reindex(columns=KEYS)
x = np.arange(21)
med, q1, q3 = piv.median(axis=0), piv.quantile(.25), piv.quantile(.75)
fig, ax = plt.subplots(figsize=(140 * MM, 68 * MM))
ax.axhline(0, color=RULE, lw=0.5, zorder=1)
ax.fill_between(x, q1, q3, color=FAINT, lw=0, zorder=2)
ax.plot(x, med, color=INK, lw=0.9, zorder=3)
thr = KEYS.index("80-84") - 0.5   # the boundary at exact age 80
ax.axvline(thr, color=ACC, lw=0.6, ls=(0, (3, 2)), zorder=2)
ax.text(thr - 0.3, q3.max() * 0.95, "median $a^{\\dagger}$", fontsize=6.4,
        color=ACC, ha="right", va="top")
ax.set_xticks(x[::2])
ax.set_xticklabels([KEYS[i] for i in x[::2]], rotation=45, ha="right")
ax.set_xlabel("age group (years)")
ax.set_ylabel("contribution to the reduction\nin $e^{\\dagger}$ (years)")
ax.set_xlim(-0.6, 20.6)
save(fig, "fig2_age_contributions")

# ------------------------------------------------------- Figure 3: reference dependence
d = P.sort_values("asyr_nat")
fig, ax = plt.subplots(figsize=(140 * MM, 140 * MM))
y = np.arange(len(d))
ax.hlines(y, d.asyr_nat, d.asyr_off, color=FAINT, lw=1.4, zorder=1)
ax.scatter(d.asyr_nat, y, s=9, facecolors="white", edgecolors=INK, linewidths=0.7, zorder=3)
ax.scatter(d.asyr_off, y, s=9, facecolors=ACC, edgecolors=ACC, linewidths=0, zorder=3)
lo = float(d.asyr_nat.min())
for i, g in enumerate(d.index):
    ax.text(lo - 6, i, NAME[g], ha="right", va="center", fontsize=6.6, color=INK)
ax.set_yticks([])
ax.set_ylim(-1, len(d) + 1.6)
ax.set_xlim(lo - 42, float(d.asyr_off.max()) + 8)
ax.set_xlabel("age-standardised years of life lost per 100 000 person-years")
for s in ["left", "top", "right"]:
    ax.spines[s].set_visible(False)
# direct labels above the first row in place of a legend: each label sits over the marker
# type it names, so no key has to be decoded
top = d.index[-1]
ax.text(d.asyr_nat[top], len(d) - 0.35, "national\nlife table", ha="center", va="bottom",
        fontsize=6.2, color=INK, linespacing=1.1)
ax.text(d.asyr_off[top], len(d) - 0.35, "GBD 2023\nreference", ha="center", va="bottom",
        fontsize=6.2, color=ACC, linespacing=1.1)
save(fig, "fig3_reference_dependence")

# --------------------------------- Supplementary: contributions under both orderings
CA = pd.read_csv("contrib_symmetric.csv")
geos = sorted(P.index, key=lambda g: NAME[g])   # alphabetical by country name
fig, axes = plt.subplots(7, 4, figsize=(170 * MM, 200 * MM), sharex=True, sharey=True)
for a, g in zip(axes.ravel(), geos):
    s = CA[CA.geo == g]
    f = s.pivot_table(index="k", values="fwd", aggfunc="mean").reindex(KEYS).values.ravel()
    r = s.pivot_table(index="k", values="rev", aggfunc="mean").reindex(KEYS).values.ravel()
    a.axhline(0, color=FAINT, lw=0.4)
    a.plot(np.arange(21), f, color=INK, lw=0.7)
    a.plot(np.arange(21), r, color=ACC, lw=0.7, ls=(0, (2, 1.6)))
    a.set_title(NAME[g], fontsize=6.2, color=INK, pad=2)
    a.tick_params(labelsize=5.4)
    for sp in ["top", "right"]:
        a.spines[sp].set_visible(False)
for a in axes[-1, :]:
    a.set_xticks([0, 5, 10, 15, 20])
    a.set_xticklabels([KEYS[i] for i in [0, 5, 10, 15, 20]], rotation=45, ha="right",
                      fontsize=5.4)
fig.supxlabel("age group (years)", fontsize=7, color=INK)
fig.supylabel("contribution to the reduction in $e^{\\dagger}$ (years)", fontsize=7, color=INK)
fig.tight_layout()
save(fig, "figS1_orderings")

# ------------------- Supplementary: reference choice across nine causes
# An aligned two-panel graphic table. Ordering the rows by mean age at death carries the age
# gradient without spending a second spatial dimension on it, and cause names are the row
# labels, so nothing has to be decoded from a legend.
# The common frame of the countries available for every cause, as in Table 5 of the main text,
# so the figure and the table report the same values.
GC = pd.read_csv("generalisation_common.csv").sort_values("mad").reset_index(drop=True)

fig = plt.figure(figsize=(170 * MM, 84 * MM))
gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.45], wspace=0.50,
                      left=0.175, right=0.880, bottom=0.175, top=0.965)
axA, axB = fig.add_subplot(gs[0]), fig.add_subplot(gs[1])
rows = np.arange(len(GC))
IS_CNS = ["cns" in c.lower() for c in GC["cause"]]

for ax, xmax, ticks in ((axA, 100, [0, 25, 50, 75, 100]),
                        (axB, 120, [0, 25, 50, 75, 100])):
    ax.set_xlim(-2, xmax)
    ax.set_ylim(len(GC) - 0.5, -0.5)
    ax.set_xticks(ticks)
    ax.tick_params(axis="x", labelsize=6.6, length=2.5, width=0.5, colors=RULE)
    ax.tick_params(axis="y", length=0)
    for sp in ["top", "right", "left"]:
        ax.spines[sp].set_visible(False)
    ax.spines["bottom"].set_color(RULE)
    ax.spines["bottom"].set_linewidth(0.5)

# left panel: mean age at death
axA.scatter(GC["mad"], rows, s=13, facecolors=[ACC if c else INK for c in IS_CNS],
            edgecolors="none", zorder=3)
axA.set_yticks(rows)
axA.set_yticklabels(GC["cause"], fontsize=6.6)
for lab, cns in zip(axA.get_yticklabels(), IS_CNS):
    lab.set_color(ACC if cns else INK)

# right panel: country range as a thin rule, mean as the dot
for k, r in GC.iterrows():
    col = ACC if IS_CNS[k] else INK
    # smallest effective difference: the row is distinguished by its name and its dots,
    # so the range rule stays the same neutral grey as every other row
    axB.hlines(k, r["infl_lo"], r["infl_hi"], color=RULE, lw=0.6, zorder=2)
    axB.scatter([r["infl_mean"]], [k], s=13, facecolors=col, edgecolors="none", zorder=3)
axB.set_yticks([])

# exact values, so the figure serves as a lookup table as well as a graphic
for k, r in GC.iterrows():
    col = ACC if IS_CNS[k] else INK
    axA.text(1.10, k, "%.0f" % r["mad"], transform=axA.get_yaxis_transform(),
             ha="right", va="center", fontsize=6.6, color=col, clip_on=False)
    axB.text(1.035, k, "%.1f (%.1f to %.1f)" % (r["infl_mean"], r["infl_lo"], r["infl_hi"]),
             transform=axB.get_yaxis_transform(), ha="left", va="center",
             fontsize=6.6, color=col, clip_on=False)

axA.set_xlabel("mean age at death (years)", fontsize=7, color=INK)
axB.set_xlabel("increase in YLL per death under the common reference (%)", fontsize=7,
               color=INK)
save(fig, "figS2_reference_gradient", tight=False)

print("figures written")
