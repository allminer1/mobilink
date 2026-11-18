"""P3 step 3.7: figures 1-6 of the preprint draft, drawn from the result tables (no analysis is run here).

Figure 1  the problem and the measure in public runs: naive share and net linked fraction against read N50
          (native reads and the same reads cut to 3 and 2 kb), per environment      tables/P3/confirm*/partB.tsv
Figure 2  the method (schematic)
Figure 3  calibration in simulation (gate G3)                                       tables/G3/G3a_settings.tsv, G3b_settings.tsv
Figure 4  the registered cut test (H1a, H2), counting runs of both stages           tables/P3/confirm*/partB.tsv
Figure 5  ARG-MGE linkage per environment and its dependence on class definitions   tables/P3/confirm*/per_run.tsv, pooled.tsv,
                                                                                    tables/P3/P3_04_sensitivity.tsv
Figure 6  contigs: the chance model decides the conclusion                          tables/P3/P3_05_*.tsv
Output: manuscript/figures/fig<N>.pdf and .png (300 dpi).
usage: envs/linkstat_p0/bin/python -I scripts/P3_07_figures.py [N ...]
"""
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import FancyBboxPatch, Rectangle  # noqa: E402

L = "/mnt/rmrfc/workflow/iskander/linkstat"
OUT = f"{L}/manuscript/figures"
INK, INK2, MUTED, GRID, AXIS, SURF = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7", "#ffffff"
BLUE, ORANGE, AQUA, YELLOW = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"
ENV = {"human gut": (BLUE, "o"), "livestock": (ORANGE, "s"), "wastewater": (AQUA, "^"), "soil": (YELLOW, "D")}
plt.rcParams.update({
    "font.size": 8, "font.family": "DejaVu Sans", "text.color": INK, "axes.labelcolor": INK2, "axes.edgecolor": AXIS,
    "axes.linewidth": 0.8, "axes.spines.top": False, "axes.spines.right": False, "axes.titlesize": 8.5,
    "axes.titleweight": "bold", "axes.titlelocation": "left", "xtick.color": AXIS, "ytick.color": AXIS,
    "xtick.labelcolor": INK2, "ytick.labelcolor": INK2, "xtick.labelsize": 7.5, "ytick.labelsize": 7.5,
    "grid.color": GRID, "grid.linewidth": 0.6, "axes.axisbelow": True, "legend.frameon": False, "legend.fontsize": 7.5,
    "figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF, "pdf.fonttype": 42,
})


def save(fig, n):
    for ext in ("pdf", "png"):
        fig.savefig(f"{OUT}/fig{n}.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"fig{n} written")


def panel(ax, letter):
    ax.text(-0.12, 1.06, letter, transform=ax.transAxes, fontsize=10, fontweight="bold", ha="right", va="bottom", color=INK)


def part_b():
    a = pd.read_csv(f"{L}/tables/P3/confirm/partB.tsv", sep="\t").assign(stage=1)
    b = pd.read_csv(f"{L}/tables/P3/confirm_ext/partB.tsv", sep="\t").assign(stage=2)
    return pd.concat([a, b], ignore_index=True)


def fig1():
    t = part_b()
    fig, axes = plt.subplots(2, 4, figsize=(7.2, 4.1), sharex=True, sharey="row")
    for j, (env, (col, mk)) in enumerate(ENV.items()):
        e = t[t.env == env]
        for i, (lab, cols) in enumerate((("Naive share", ("r2000_naive", "r3000_naive", "intact_naive")),
                                         ("Net linked fraction within 2 kb", ("r2000_net", "r3000_net", "intact_net")))):
            ax = axes[i, j]
            n = 0
            for r in e.itertuples():
                x = np.array([r.r2000_n50, r.r3000_n50, r.intact_n50]) / 1000
                y = np.array([getattr(r, c) for c in cols], float)
                ok = np.isfinite(y) & (x <= x[-1] + 1e-9)
                if ok.sum() < 2 or not ok[-1]:
                    continue
                n += 1
                ax.plot(x[ok], y[ok], color=col, lw=0.9, alpha=0.75, solid_capstyle="round")
                ax.plot(x[-1], y[-1], mk, color=col, ms=3.6, mec=SURF, mew=0.6)
            ax.set_xscale("log")
            ax.set_xticks([2, 3, 5, 10, 20])
            ax.set_xticklabels(["2", "3", "5", "10", "20"])
            ax.minorticks_off()
            ax.grid(axis="y")
            if i == 0:
                ax.set_title(f"{env.capitalize()} ({n} runs)", fontweight="bold")
            else:
                ax.axhline(0, color=AXIS, lw=0.8)
                ax.set_xlabel("Read N50 (kb)")
            if j == 0:
                ax.set_ylabel(lab)
    axes[0, 0].set_ylim(0, 0.75)
    panel(axes[0, 0], "a")
    panel(axes[1, 0], "b")
    fig.legend(handles=[Line2D([], [], marker="o", color=INK2, lw=0.9, ms=3.6, label="native reads (marker) and the same reads cut to N50 3 kb and 2 kb")],
               loc="lower center", ncol=1, bbox_to_anchor=(0.5, -0.03))
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    save(fig, 1)


def gene(ax, x0, x1, y, color, label=None, h=0.36, text_color="#ffffff"):
    ax.add_patch(FancyBboxPatch((x0, y - h / 2), x1 - x0, h, boxstyle="round,pad=0,rounding_size=0.08", fc=color, ec="none"))
    if label:
        ax.text((x0 + x1) / 2, y, label, ha="center", va="center", fontsize=7, color=text_color, fontweight="bold")


def read(ax, x0, x1, y):
    ax.plot([x0, x1], [y, y], color=AXIS, lw=3.2, solid_capstyle="round", zorder=0)


def fig2():
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.9), gridspec_kw={"width_ratios": [1.25, 1, 1]})
    for ax in axes:
        ax.set_xlim(0, 10)
        ax.set_ylim(0, 6)
        ax.axis("off")
    ax = axes[0]
    ax.set_title("Distance to the nearest B", loc="left")
    read(ax, 0.3, 9.7, 4.6)
    gene(ax, 3.0, 3.9, 4.6, BLUE, "A")
    gene(ax, 4.1, 4.9, 4.6, BLUE, "A")
    gene(ax, 6.0, 7.4, 4.6, ORANGE, "B")
    ax.annotate("", xy=(7.4, 5.25), xytext=(4.9, 5.25), arrowprops=dict(arrowstyle="<->", color=INK2, lw=0.8))
    ax.text(7.55, 5.25, "event: far end of B", ha="left", va="center", fontsize=7, color=INK2)
    ax.annotate("", xy=(0.3, 3.95), xytext=(3.0, 3.95), arrowprops=dict(arrowstyle="<->", color=INK2, lw=0.8))
    ax.text(1.65, 3.8, "censored: read ends", ha="center", va="top", fontsize=7, color=INK2)
    ax.plot([3.0, 4.9], [5.02, 5.02], color=BLUE, lw=1.2)
    ax.text(2.85, 5.02, "one unit", ha="right", va="center", fontsize=7, color=INK2)
    read(ax, 0.3, 6.2, 2.4)
    gene(ax, 2.2, 3.1, 2.4, BLUE, "A")
    gene(ax, 5.3, 6.2, 2.4, ORANGE, "B", h=0.36)
    ax.plot([6.2, 6.2], [2.05, 2.75], color=INK, lw=0.9)
    ax.text(6.35, 2.4, "B cut by the read end:\ncensored, not an event", ha="left", va="center", fontsize=7, color=INK2)
    ax.text(0.3, 1.2, "Kaplan–Meier over all sides gives $F_A(w)$:\nthe share of A units with a whole B within $w$.",
            ha="left", va="center", fontsize=7.5, color=INK)
    ax = axes[1]
    ax.set_title("Chance: anchors, same host", loc="left")
    for y, lab in ((4.9, "host 1"), (4.2, None), (3.5, None)):
        read(ax, 0.5, 9.3, y)
        if lab:
            ax.text(0.5, y + 0.38, "reads of one host", ha="left", va="bottom", fontsize=7, color=INK2)
    gene(ax, 2.0, 2.9, 4.9, BLUE, "A")
    gene(ax, 6.2, 7.5, 4.9, ORANGE, "B")
    gene(ax, 4.6, 5.9, 3.5, ORANGE, "B")
    for x0, y in ((5.0, 4.9), (1.2, 4.2), (6.3, 4.2), (1.5, 3.5), (7.6, 3.5)):
        ax.add_patch(Rectangle((x0, y - 0.18), 0.9, 0.36, fc="none", ec=INK2, lw=0.9, ls=(0, (2, 1.5))))
    ax.text(0.5, 2.55, "Anchors (dashed) give $F_R(w)$.", ha="left", va="center", fontsize=7.5, color=INK)
    ax.text(0.5, 1.35, r"net linked $= 1-\dfrac{1-F_A(w)}{1-F_R(w)}$", ha="left", va="center", fontsize=9, color=INK)
    ax.text(0.5, 0.45, "0 under chance, at any read length", ha="left", va="center", fontsize=7.5, color=INK2)
    ax = axes[2]
    ax.set_title("Test: one locus, one clump", loc="left")
    for k, (x0, x1) in enumerate(((0.6, 6.8), (1.6, 8.4), (2.4, 9.4))):
        y = 4.9 - 0.7 * k
        read(ax, x0, x1, y)
        gene(ax, 3.4, 4.3, y, BLUE, "A")
        gene(ax, 5.2, 6.4, y, ORANGE, "B")
    ax.add_patch(FancyBboxPatch((0.3, 3.05), 9.4, 2.3, boxstyle="round,pad=0,rounding_size=0.25", fc="none", ec=INK2, lw=0.8))
    ax.text(9.6, 2.85, "1 clump, not 3 observations", ha="right", va="top", fontsize=7, color=INK2)
    ax.text(0.3, 1.55, "Null variance over clumps;\np value by Monte Carlo (mid-p).", ha="left", va="center", fontsize=7.5, color=INK)
    ax.text(0.3, 0.45, "Contigs: one clump per contig.", ha="left", va="center", fontsize=7.5, color=INK2)
    for ax, s in zip(axes, "abc"):
        ax.set_title(f"{s}   " + ax.get_title(loc="left"))
    fig.tight_layout(w_pad=2.0)
    save(fig, 2)


def fig3():
    g = pd.read_csv(f"{L}/tables/G3/G3a_settings.tsv", sep="\t")
    f = pd.read_csv(f"{L}/tables/G3/G3b_settings.tsv", sep="\t")
    s = g[g.linked == 0].sort_values("setting")
    lab = [f"{r.profile} · {r.depth:g}× · {'islands' if r.A == 'islands' else 'single genes'}" for r in s.itertuples()]
    y = np.arange(len(s))[::-1]
    ypos = dict(zip(s.setting, y))
    fig, axes = plt.subplots(1, 4, figsize=(7.2, 4.0), sharey=True)
    two = (("net_w", BLUE, "o", "within 2 kb"), ("net_mol", ORANGE, "s", "same molecule"))
    ax = axes[0]
    for k, (c, col, mk, _) in enumerate(two):
        ax.plot(s[f"{c}_reject"], y + (0.16 if k == 0 else -0.16), mk, color=col, ms=4.2, mec=SURF, mew=0.6)
    ax.axvline(0.05, color=INK2, lw=0.8)
    ax.axvline(0.08, color=AXIS, lw=0.8)
    ax.set_xlim(0, 0.1)
    ax.set_xticks([0, 0.05, 0.1])
    ax.set_xticklabels(["0", "5%", "10%"])
    ax.set_title("False positives\n(no linkage)")
    ax.set_yticks(y)
    ax.set_yticklabels(lab)
    ax = axes[1]
    for lk, col in ((0.0, "#86b6ef"), (0.03, "#3987e5"), (0.1, "#104281")):
        e = g[g.linked == lk]
        ax.plot(e.cover_w, [ypos[i] for i in e.setting], "o", color=col, ms=4.0, mec=SURF, mew=0.6, label=f"linked share {lk:g}")
    ax.axvline(0.95, color=INK2, lw=0.8)
    ax.axvline(0.90, color=AXIS, lw=0.8)
    ax.set_xlim(0.88, 1.005)
    ax.set_xticks([0.9, 0.95, 1.0])
    ax.set_xticklabels(["90%", "95%", "100%"])
    ax.set_title("Coverage\n(95% interval)")
    cov_handles = [Line2D([], [], marker="o", color=col, lw=0, ms=4.2, label=f"coverage at linked share {lk:g}")
                   for lk, col in ((0.0, "#86b6ef"), (0.03, "#3987e5"), (0.1, "#104281"))]
    ax = axes[2]
    e = g[g.linked == 0.03]
    for k, (c, col, mk, _) in enumerate(two):
        ax.plot(e[f"{c}_reject"], [ypos[i] + (0.16 if k == 0 else -0.16) for i in e.setting], mk, color=col, ms=4.2, mec=SURF, mew=0.6)
    ax.set_xlim(0, 1.04)
    ax.set_xticks([0, 0.5, 1])
    ax.set_xticklabels(["0", "50%", "100%"])
    ax.set_title("Power\n(linked share 0.03)")
    ax = axes[3]
    for k, (c, col, mk, _) in enumerate(two):
        e = f[f.test == c]
        ax.plot(e.fdr, [ypos[i] + (0.16 if k == 0 else -0.16) for i in e.setting], mk, color=col, ms=4.2, mec=SURF, mew=0.6)
    ax.axvline(0.05, color=INK2, lw=0.8)
    ax.set_xlim(0, 0.1)
    ax.set_xticks([0, 0.05, 0.1])
    ax.set_xticklabels(["0", "5%", "10%"])
    ax.set_title("False discoveries\n(Benjamini–Hochberg)")
    for ax, sx in zip(axes, "abcd"):
        ax.grid(axis="y")
        ax.tick_params(axis="y", length=0)
        ax.set_title(f"{sx}   " + ax.get_title(loc="left").replace("\n", "\n     "))
        for b in (5.5, 11.5):
            ax.axhline(b, color=AXIS, lw=0.8)
    fig.legend(handles=[Line2D([], [], marker=mk, color=col, lw=0, ms=4.6, label=f"test: {nm}") for _, col, mk, nm in two] + cov_handles,
               loc="lower center", ncol=3, bbox_to_anchor=(0.5, -0.03), columnspacing=1.4, handletextpad=0.3)
    fig.tight_layout(rect=(0, 0.08, 1, 1), w_pad=1.6)
    save(fig, 3)


def fig4():
    t = part_b()
    c = t[t.counts].copy()
    c["order"] = c.env.map({k: i for i, k in enumerate(ENV)})
    c = c.sort_values(["order", "intact_net"], ascending=[True, False]).reset_index(drop=True)
    y = np.arange(len(c))[::-1]
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.9), sharey=True)
    for r, yy in zip(c.itertuples(), y):
        col, mk = ENV[r.env]
        ax = axes[0]
        ax.plot([r.r3000_naive, r.intact_naive], [yy, yy], color=col, lw=1.6, alpha=0.55, solid_capstyle="round")
        ax.plot(r.intact_naive, yy, mk, color=col, ms=5.2, mec=SURF, mew=0.7)
        ax.plot(r.r3000_naive, yy, mk, color=SURF, ms=4.6, mec=col, mew=1.2)
        ax.text(r.intact_naive + 0.018, yy, f"{r.naive_ratio:.1f}×", va="center", ha="left", fontsize=7, color=INK2)
        ax = axes[1]
        ax.plot([r.intact_lo, r.intact_hi], [yy + 0.17, yy + 0.17], color=col, lw=1.4, solid_capstyle="round")
        ax.plot(r.intact_net, yy + 0.17, mk, color=col, ms=5.2, mec=SURF, mew=0.7)
        ax.plot([r.r3000_lo, r.r3000_hi], [yy - 0.17, yy - 0.17], color=col, lw=1.4, alpha=0.55, solid_capstyle="round")
        ax.plot(r.r3000_net, yy - 0.17, mk, color=SURF, ms=4.6, mec=col, mew=1.2)
    axes[0].set_yticks(y)
    axes[0].set_yticklabels([f"{r.run}  (stage {r.stage})" for r in c.itertuples()])
    axes[0].set_xlim(0, 0.8)
    axes[0].set_xlabel("Naive share (ARG reads that also carry an MGE)")
    axes[0].set_title("Naive share: native and cut to 3 kb")
    axes[1].axvline(0, color=AXIS, lw=0.8)
    axes[1].set_xlabel("Net linked fraction within 2 kb (95% interval)")
    axes[1].set_title("Net linked fraction: native and cut to 3 kb")
    for ax, sx in zip(axes, "ab"):
        ax.grid(axis="x")
        ax.tick_params(axis="y", length=0)
        ax.set_title(f"{sx}   " + ax.get_title(loc="left"))
    h = [Line2D([], [], marker=mk, color=col, lw=0, ms=5, label=env) for env, (col, mk) in ENV.items() if env in set(c.env)]
    h += [Line2D([], [], marker="o", color=INK2, lw=0, ms=5, label="native reads"),
          Line2D([], [], marker="o", color=SURF, mec=INK2, mew=1.2, lw=0, ms=4.6, label="cut to N50 3 kb")]
    fig.legend(handles=h, loc="lower center", ncol=5, bbox_to_anchor=(0.5, -0.03), columnspacing=1.2, handletextpad=0.3)
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    save(fig, 4)


def fig5():
    pr = pd.concat([pd.read_csv(f"{L}/tables/P3/confirm/per_run.tsv", sep="\t").assign(stage=1),
                    pd.read_csv(f"{L}/tables/P3/confirm_ext/per_run.tsv", sep="\t").assign(stage=2)], ignore_index=True)
    pr = pr[(pr.class_a == "ARG") & (pr.class_b == "MGE")]
    po = pd.read_csv(f"{L}/tables/P3/confirm/pooled.tsv", sep="\t")
    po = po[(po.class_a == "ARG") & (po.class_b == "MGE")].set_index("env")
    se = pd.read_csv(f"{L}/tables/P3/P3_04_sensitivity.tsv", sep="\t")
    se = se[(se.class_a == "ARG") & (se.class_b == "MGE")]
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.3), gridspec_kw={"width_ratios": [1.1, 1]})
    ax = axes[0]
    rng = np.random.default_rng(1)
    for i, (env, (col, mk)) in enumerate(ENV.items()):
        yy = len(ENV) - 1 - i
        for stage, face in ((1, col), (2, SURF)):
            e = pr[(pr.env == env) & (pr.stage == stage) & np.isfinite(pr.net_w)]
            ax.plot(e.net_w, yy - 0.2 + rng.uniform(-0.11, 0.11, len(e)), mk, color=face, mec=col if stage == 2 else SURF,
                    mew=1.0 if stage == 2 else 0.5, ms=4.0, alpha=0.9)
        p = po.loc[env]
        ax.plot([p.net_w_lo, p.net_w_hi], [yy + 0.22, yy + 0.22], color=INK, lw=1.6, solid_capstyle="round")
        ax.plot(p.net_w, yy + 0.22, "o", color=INK, ms=5.4, mec=SURF, mew=0.7)
        ax.text(p.net_w_hi + 0.012, yy + 0.22, f"{p.net_w:.3f}", va="center", ha="left", fontsize=7, color=INK2)
    ax.axvline(0, color=AXIS, lw=0.8)
    ax.set_yticks(range(len(ENV)))
    ax.set_yticklabels([e.capitalize() for e in list(ENV)[::-1]])
    ax.set_ylim(-0.6, len(ENV) - 0.4)
    ax.set_xlabel("ARG–MGE net linked fraction within 2 kb")
    ax.set_title("Per run and pooled per environment")
    h5 = [Line2D([], [], marker="o", color=INK, lw=1.6, ms=5.4, label="pooled, stage 1 (95% interval)"),
          Line2D([], [], marker="o", color=INK2, lw=0, ms=4, label="run, stage 1"),
          Line2D([], [], marker="o", color=SURF, mec=INK2, mew=1.0, lw=0, ms=4, label="run, stage 2")]
    ax = axes[1]
    variants = [("primary", "registered classes"), ("chimera_filtered", "chimera-flagged reads removed"),
                ("genus_strata", "hosts at genus rank"), ("ARG_noint", "without intrinsic ARG"),
                ("MGE_80", "MGE at ≥ 80% identity"), ("ARG_MRG_60", "ARG at ≥ 60% identity"),
                ("MGE_D3", "MGE: broad categories")]
    for k, (v, lab) in enumerate(variants):
        yy = len(variants) - 1 - k
        for env, (col, mk) in ENV.items():
            e = se[(se.variant == v) & (se.env == env)]
            if len(e) and np.isfinite(e.net_w.iloc[0]):
                ax.plot(e.net_w.iloc[0], yy, mk, color=col, ms=5.0, mec=SURF, mew=0.7)
    ax.axvline(0, color=AXIS, lw=0.8)
    ax.set_yticks(range(len(variants)))
    ax.set_yticklabels([lab for _, lab in variants][::-1])
    ax.set_ylim(-0.6, len(variants) - 0.4)
    ax.set_xlabel("ARG–MGE net linked fraction within 2 kb, pooled")
    ax.set_title("Dependence on the analysis choices")
    h5 = [Line2D([], [], marker=mk, color=col, lw=0, ms=5, label=env) for env, (col, mk) in ENV.items()] + h5
    for ax, sx in zip(axes, "ab"):
        ax.grid(axis="x")
        ax.tick_params(axis="y", length=0)
        ax.set_title(f"{sx}   " + ax.get_title(loc="left"))
    fig.legend(handles=h5, loc="lower center", ncol=7, bbox_to_anchor=(0.5, -0.03), columnspacing=1.0, handletextpad=0.3)
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    save(fig, 5)


def fig6():
    s = pd.read_csv(f"{L}/tables/P3/P3_05_null_sim.summary.tsv", sep="\t").set_index("scenario")
    n = pd.read_csv(f"{L}/tables/P3/P3_05_pilot.nulls.tsv", sep="\t")
    n = n[n["sample"] == "Ridder (sum of 4)"]
    ls = pd.read_csv(f"{L}/tables/P3/P3_05_pilot.linkstat.tsv", sep="\t")
    ls = ls[(ls["sample"] == "Ridder (pooled)") & (ls.strata == "MAG") & (ls.class_a == "ARG") & (ls.w == 2000)].iloc[0]
    c = pd.read_csv(f"{L}/tables/P3/P3_05_soil_contigs.linkstat.tsv", sep="\t")
    c = c[c.assembly.str.startswith("mylo_") & (c.class_set == "relaxed60") & (c.strata == "bin") & (c.class_a == "ARG") & (c.class_b == "MRG") & (c.w == 2000)]
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 3.1), gridspec_kw={"width_ratios": [1.0, 1.0, 0.9]})
    ax = axes[0]
    scen = [("S0", "genes placed\nuniformly"), ("S1", "genes in clusters\nwider than 1 kb"), ("S2", "contigs differ in\ngene density")]
    ser = [("perm_loci_depletion_p05", ORANGE, "s", "label permutation"), ("mol_ab_p05", "#86b6ef", "o", "mobilink, no clumps"),
           ("cc_mol_ab_p05", BLUE, "o", "mobilink, contig clumps")]
    for k, (sc, _) in enumerate(scen):
        yy = len(scen) - 1 - k
        for m, (col_name, col, mk, _) in enumerate(ser):
            ax.plot(100 * s.loc[sc, col_name], yy + 0.22 * (1 - m), mk, color=col, ms=5.4, mec=SURF, mew=0.7)
    ax.axvline(5, color=INK2, lw=0.8)
    ax.set_xlim(-4, 104)
    ax.set_xticks([0, 50, 100])
    ax.set_xticklabels(["0", "50%", "100%"])
    ax.set_yticks(range(len(scen)))
    ax.set_yticklabels([lab for _, lab in scen][::-1])
    ax.set_ylim(-0.6, len(scen) - 0.4)
    ax.set_xlabel("No-linkage datasets\ncalled at 0.05")
    ax.set_title("Simulation without linkage")
    ax.legend(handles=[Line2D([], [], marker=mk, color=col, lw=0, ms=5.4, label=nm) for _, col, mk, nm in ser],
              loc="upper center", bbox_to_anchor=(0.3, -0.36), handletextpad=0.3)
    ax = axes[1]
    rows = [("label permutation, loci", n[(n.slots == "loci") & (n.null == "C label permutation")].iloc[0]),
            ("label permutation, units", n[(n.slots == "units") & (n.null == "C label permutation")].iloc[0]),
            ("label permutation, units in MAG", n[(n.slots == "units") & (n.null == "C label permutation within host")].iloc[0]),
            ("length-weighted placement", n[(n.slots == "loci") & (n.null == "A length-weighted placement")].iloc[0])]
    for k, (lab, r) in enumerate(rows):
        yy = len(rows) - k
        ax.plot([r.lo95, r.hi95], [yy, yy], color=ORANGE, lw=1.6, solid_capstyle="round")
        ax.plot(r.expected, yy, "s", color=ORANGE, ms=5.2, mec=SURF, mew=0.7)
        ax.text(r.hi95 + 1.0, yy, f"{r.expected:.1f}", va="center", ha="left", fontsize=7, color=INK2)
    ax.plot(ls.oe_n2u_expected, 0, "o", color=BLUE, ms=5.4, mec=SURF, mew=0.7)
    ax.text(ls.oe_n2u_expected + 1.6, 0, f"{ls.oe_n2u_expected:.1f}", va="center", ha="left", fontsize=7, color=INK2)
    ax.axvline(4, color=INK, lw=0.9)
    ax.text(4.9, len(rows) + 0.8, "observed: 4", ha="left", va="center", fontsize=7, color=INK)
    ax.set_yticks(range(len(rows) + 1))
    ax.set_yticklabels(["mobilink (positions kept)"] + [lab for lab, _ in rows][::-1])
    ax.set_ylim(-0.6, len(rows) + 1.2)
    ax.set_xlim(0, 42)
    ax.set_xlabel("Contigs with ARG and MRG\nexpected by chance")
    ax.set_title("Published loci, four soils")
    ax = axes[2]
    c = c.assign(soil=c.assembly.str.replace("mylo_", "")).sort_values("soil")
    yv = np.arange(len(c))[::-1]
    ax.axvspan(0.8, 1.25, color=GRID, lw=0)
    ax.axvline(1, color=INK2, lw=0.8)
    for r, yy in zip(c.itertuples(), yv):
        if np.isfinite(r.mol_ratio_lo90):
            ax.plot([r.mol_ratio_lo90, r.mol_ratio_hi90], [yy, yy], color=BLUE, lw=1.6, solid_capstyle="round")
        ax.plot(r.mol_ratio, yy, "o", color=BLUE, ms=5.2, mec=SURF, mew=0.7)
    ax.set_xscale("log")
    ax.set_xlim(0.15, 14)
    ax.set_xticks([0.25, 1, 4])
    ax.set_xticklabels(["0.25", "1", "4"])
    ax.minorticks_off()
    ax.set_yticks(yv)
    ax.set_yticklabels([f"soil {x}" for x in c.soil])
    ax.set_xlabel("Share with MRG on the contig:\nARG / chance (90% interval)")
    ax.set_title("Eight soils, ≥ 60% identity")
    for ax, sx in zip(axes, "abc"):
        ax.grid(axis="x")
        ax.tick_params(axis="y", length=0)
        ax.set_title(f"{sx}   " + ax.get_title(loc="left"))
    fig.tight_layout()
    save(fig, 6)


if __name__ == "__main__":
    import os

    os.makedirs(OUT, exist_ok=True)
    todo = [int(x) for x in sys.argv[1:]] or [1, 2, 3, 4, 5, 6]
    for k in todo:
        {1: fig1, 2: fig2, 3: fig3, 4: fig4, 5: fig5, 6: fig6}[k]()
