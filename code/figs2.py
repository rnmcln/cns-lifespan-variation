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
    "font.sans-serif": ["DejaVu Sans", "Helvetica", "Arial"],
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



def _crosses(p1, p2, p3, p4):
    """True if segment p1-p2 properly intersects segment p3-p4."""
    def side(a, b, c):
        return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])
    d1, d2 = side(p3, p4, p1), side(p3, p4, p2)
    d3, d4 = side(p1, p2, p3), side(p1, p2, p4)
    return ((d1 > 0) != (d2 > 0)) and ((d3 > 0) != (d4 > 0))


def place_labels(fig, ax, xs, ys, los, his, labels, colours, fontsize=7.2, audit=False,
                 tagged=False):
    """Place one text label per point, solved against the rendered geometry.

    Each label is anchored to its own marker, never to the end of a whisker: anchoring to a
    whisker end, as an earlier version of this figure did, put a label a distance from its
    own marker that varied with the length of that cause's country range, which is what made
    the association ambiguous for the wide-range causes.

    Candidate positions in a bounded neighbourhood of each marker are scored against the
    area they overlap with markers, whiskers and already-placed labels, plus penalties for
    leader lines that cross each other or another cause's whisker, plus a distance term. The
    greedy pass is repeated under several orderings and the cheapest complete layout wins,
    which avoids the outcome where one label placed early leaves a later one no good option.

    All geometry is in display pixels, converted to figure fractions before being set:
    fractions are invariant to canvas dpi, whereas pixels are not, and the Agg canvas used
    for measurement does not share a dpi with the PDF canvas used for output.
    """
    from matplotlib.transforms import Bbox

    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    PX = fig.dpi / 72.0
    T = ax.transData.transform
    pts = [tuple(T((x, y))) for x, y in zip(xs, ys)]
    n = len(labels)

    sizes = []
    for lab in labels:
        t = ax.annotate(lab, (0, 0), xycoords="figure fraction",
                        textcoords="figure fraction", xytext=(0, 0),
                        ha="center", va="center", fontsize=fontsize)
        w, h = t.get_window_extent(rend).size
        t.remove()
        if tagged:
            # the solver must see the disc, not the numeral inside it, or it will place
            # tags in gaps that only the text would have fitted
            d = max(w, h) + 5.0 * PX
            w = h = d
        sizes.append((w, h))

    marks = [Bbox.from_bounds(px - 5 * PX, py - 5 * PX, 10 * PX, 10 * PX)
             for px, py in pts]
    whisk = []
    for x, lo, hi in zip(xs, los, his):
        (x0, y0), (_, y1) = T((x, lo)), T((x, hi))
        whisk.append(Bbox.from_bounds(x0 - 2.5 * PX, min(y0, y1), 5 * PX, abs(y1 - y0)))
    frame = ax.get_window_extent(rend)

    ANGLES = list(range(0, 360, 10))
    # a label must stay in its own point's neighbourhood: beyond roughly 40 points the
    # association stops being readable however good the leader line is
    RADII = [11, 14, 18, 22, 27, 33, 40, 46]

    def candidates(i):
        w, h = sizes[i]
        out = []
        for rad in RADII:
            for ang in ANGLES:
                dx = rad * PX * np.cos(np.radians(ang))
                dy = rad * PX * np.sin(np.radians(ang))
                sx = 1 if dx > 1e-9 else (-1 if dx < -1e-9 else 0)
                sy = 1 if dy > 1e-9 else (-1 if dy < -1e-9 else 0)
                cx = pts[i][0] + dx + sx * (w / 2 + 2 * PX)
                cy = pts[i][1] + dy + sy * (h / 2 + 1.5 * PX)
                box = Bbox.from_bounds(cx - w / 2 - PX, cy - h / 2 - PX,
                                       w + 2 * PX, h + 2 * PX)
                out.append((box, cx, cy, rad))
        return out

    def cost_of(i, cand, chosen):
        box, cx, cy, rad = cand
        c = rad * 9.0
        for b in marks:
            ov = Bbox.intersection(box, b)
            if ov is not None and ov.width > 0 and ov.height > 0:
                c += 4.0 * ov.width * ov.height / (PX * PX)
        for k, b in enumerate(whisk):
            ov = Bbox.intersection(box, b)
            if ov is not None and ov.width > 0 and ov.height > 0:
                # over a cause's own whisker is unambiguous and only mildly untidy; over
                # another cause's whisker is neither
                c += (2.0 if k == i else 26.0) * ov.width * ov.height / (PX * PX)
        mid = ((box.x0 + box.x1) / 2, (box.y0 + box.y1) / 2)
        for j, jb in chosen.items():
            ov = Bbox.intersection(box, jb[0])
            if ov is not None and ov.width > 0 and ov.height > 0:
                c += 8.0 * ov.width * ov.height / (PX * PX)
            jm = ((jb[0].x0 + jb[0].x1) / 2, (jb[0].y0 + jb[0].y1) / 2)
            if _crosses(pts[i], mid, pts[j], jm):
                c += 1200
        for k, b in enumerate(whisk):
            if k == i:
                continue
            cor = [(b.x0, b.y0), (b.x0, b.y1), (b.x1, b.y1), (b.x1, b.y0)]
            if any(_crosses(pts[i], mid, a, d) for a, d in zip(cor, cor[1:] + cor[:1])):
                c += 500
                break
        if box.x0 < frame.x0 + 3 * PX or box.x1 > frame.x1 - 10 * PX:
            c += 1e5
        if box.y0 < frame.y0 + 3 * PX or box.y1 > frame.y1 - 3 * PX:
            c += 1e5
        return c

    crowd = [sum(1 for q in pts if q != p and np.hypot(p[0] - q[0], p[1] - q[1]) < 95 * PX)
             for p in pts]
    orders = [
        sorted(range(n), key=lambda i: (-crowd[i], -ys[i])),
        sorted(range(n), key=lambda i: (-crowd[i], xs[i])),
        sorted(range(n), key=lambda i: -ys[i]),
        sorted(range(n), key=lambda i: ys[i]),
        sorted(range(n), key=lambda i: xs[i]),
        sorted(range(n), key=lambda i: -xs[i]),
    ]
    best_layout, best_total = None, None
    for order in orders:
        chosen, total = {}, 0.0
        for i in order:
            pick, pc = None, None
            for cand in candidates(i):
                c = cost_of(i, cand, chosen)
                if pc is None or c < pc:
                    pick, pc = cand, c
            chosen[i] = pick
            total += pc
        if best_total is None or total < best_total:
            best_layout, best_total = dict(chosen), total

    texts = []
    for i in range(n):
        box, cx, cy, rad = best_layout[i]
        pos = tuple(fig.transFigure.inverted().transform((cx, cy)))
        if tagged:
            # a filled disc marks the numeral as an identifier rather than a second data
            # annotation, and is visually distinct from the open-circle data markers
            txt = ax.annotate(labels[i], (0, 0), xycoords="figure fraction",
                              textcoords="figure fraction", xytext=pos,
                              ha="center", va="center", fontsize=fontsize,
                              color="white", zorder=6, annotation_clip=False,
                              bbox=dict(boxstyle="circle,pad=0.18", linewidth=0,
                                        facecolor=colours[i], edgecolor="none"))
        else:
            txt = ax.annotate(labels[i], (0, 0), xycoords="figure fraction",
                              textcoords="figure fraction", xytext=pos,
                              ha="center", va="center", fontsize=fontsize,
                              color=colours[i], zorder=5, annotation_clip=False)
        txt.set_clip_on(False)
        texts.append(txt)
        centre = ((box.x0 + box.x1) / 2, (box.y0 + box.y1) / 2)
        d_own = np.hypot(centre[0] - pts[i][0], centre[1] - pts[i][1])
        d_other = min(np.hypot(centre[0] - q[0], centre[1] - q[1])
                      for j, q in enumerate(pts) if j != i)
        # a leader is drawn when the label is materially displaced, or when another marker
        # is comparably close to it, since proximity alone would then not identify the point
        if rad > (7 if len(labels[i]) <= 2 else 12) or d_other <= 1.35 * d_own:
            ex = min(max(pts[i][0], box.x0), box.x1)
            ey = min(max(pts[i][1], box.y0), box.y1)
            inv = ax.transData.inverted()
            ax.annotate("", xy=inv.transform(pts[i]), xytext=inv.transform((ex, ey)),
                        xycoords="data", textcoords="data",
                        arrowprops=dict(arrowstyle="-", color=RULE, lw=0.5,
                                        shrinkA=1.5, shrinkB=4.0),
                        zorder=1, annotation_clip=False)

    if audit:
        bad = 0
        boxes = {i: best_layout[i][0] for i in range(n)}
        for i in range(n):
            for j in range(i + 1, n):
                ov = Bbox.intersection(boxes[i], boxes[j])
                if ov is not None and ov.width > 0 and ov.height > 0:
                    print("  label-label overlap: %s / %s" % (labels[i], labels[j]))
                    bad += 1
            for k, b in enumerate(marks):
                ov = Bbox.intersection(boxes[i], b)
                if ov is not None and ov.width > 0 and ov.height > 0:
                    print("  label-marker overlap: %s over %s" % (labels[i], labels[k]))
                    bad += 1
            for k, b in enumerate(whisk):
                if k == i:
                    continue
                ov = Bbox.intersection(boxes[i], b)
                if ov is not None and ov.width > 0 and ov.height > 0:
                    print("  label-whisker overlap: %s over %s" % (labels[i], labels[k]))
                    bad += 1
            if (boxes[i].x0 < frame.x0 or boxes[i].x1 > frame.x1
                    or boxes[i].y0 < frame.y0 or boxes[i].y1 > frame.y1):
                print("  label outside frame: %s" % labels[i])
                bad += 1
            centre = ((boxes[i].x0 + boxes[i].x1) / 2, (boxes[i].y0 + boxes[i].y1) / 2)
            near = min(range(n), key=lambda j: np.hypot(centre[0] - pts[j][0],
                                                        centre[1] - pts[j][1]))
            if near != i:
                print("  nearest marker to '%s' is '%s' (leader required)"
                      % (labels[i], labels[near]))
        print("  placement audit: %d issue(s)" % bad)
    return texts


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
ax.set_ylim(-1, len(d) + 1.4)
ax.set_xlim(lo - 42, float(d.asyr_off.max()) + 8)
ax.set_xlabel("age-standardised years of life lost per 100 000 person-years")
for s in ["left", "top", "right"]:
    ax.spines[s].set_visible(False)
ax.scatter([], [], s=9, facecolors="white", edgecolors=INK, linewidths=0.7,
           label="national period life table")
ax.scatter([], [], s=9, facecolors=ACC, edgecolors=ACC, linewidths=0,
           label="GBD 2023 reference")
ax.legend(loc="lower right", fontsize=6.4, handletextpad=0.4, borderpad=0.2,
          labelspacing=0.35)
save(fig, "fig3_reference_dependence")

# --------------------------------- Supplementary: contributions under both orderings
CA = pd.read_csv("contrib_symmetric.csv")
geos = sorted(P.index)
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
axes[-1, 0].set_xticks([0, 5, 10, 15, 20])
axes[-1, 0].set_xticklabels([KEYS[i] for i in [0, 5, 10, 15, 20]], rotation=45, ha="right",
                            fontsize=5.4)
fig.supxlabel("age group (years)", fontsize=7, color=INK)
fig.supylabel("contribution to the reduction in $e^{\\dagger}$ (years)", fontsize=7, color=INK)
fig.tight_layout()
save(fig, "figS1_orderings")

# ------------------- Supplementary: reference choice across nine causes
# An aligned two-panel graphic table rather than a scatterplot. The comparison the reader
# needs is where each cause sits relative to the others, and ordering the rows by mean age
# at death carries the age gradient without spending a second spatial dimension on it. Cause
# names are the row labels, so nothing has to be decoded from a legend and no leader lines
# cross the data field. Earlier versions of this figure were a scatterplot whose nine labels,
# and then nine numerals with a key, had to be fitted around the range bars; both spent more
# ink on identification than on the observations.
GC = pd.read_csv("generalisation_causes.csv").sort_values("mad").reset_index(drop=True)

fig = plt.figure(figsize=(170 * MM, 84 * MM))
gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.45], wspace=0.50,
                      left=0.175, right=0.880, bottom=0.155, top=0.855)
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
    for lab in ax.get_xticklabels():
        lab.set_color(INK)

# left panel: mean age at death
axA.scatter(GC["mad"], rows, s=13, facecolors=[ACC if c else INK for c in IS_CNS],
            edgecolors="none", zorder=3)
axA.set_yticks(rows)
axA.set_yticklabels(GC["cause"], fontsize=7.2)
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

axA.text(0.0, 1.055, "Mean age at death (years)", transform=axA.transAxes,
         ha="left", va="bottom", fontsize=7.4, color=INK)
axB.text(0.0, 1.055, "Increase in YLL per death under common reference (%)",
         transform=axB.transAxes, ha="left", va="bottom", fontsize=7.4, color=INK)
save(fig, "figS2_reference_gradient", tight=False)

print("figures written")
