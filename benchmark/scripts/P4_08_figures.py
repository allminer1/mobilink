"""Phase 5: figures of manuscript draft v2 that come from Phase 4 (no analysis is run here; style of P3_07_figures.py).

Figure 5  ARG-MGE linkage in 241 public runs                tables/P4/confirm/{classes,projects,per_run,partB,h1c}.tsv,
                                                            tables/P4/sens/sens_variants.classes.tsv
Output: manuscript/figures/fig5_v2.pdf and .png (300 dpi). The v1 figure 5 is kept as fig5.pdf.
usage: envs/linkstat_p0/bin/python -I scripts/P4_08_figures.py [5]
"""
import importlib.util
import sys

import numpy as np
import pandas as pd

L = "/mnt/rmrfc/workflow/iskander/linkstat"
_s = importlib.util.spec_from_file_location("p3_07_figures", f"{L}/scripts/P3_07_figures.py")
S = importlib.util.module_from_spec(_s)
_s.loader.exec_module(S)
plt, INK, INK2, MUTED, GRID, AXIS, SURF = S.plt, S.INK, S.INK2, S.MUTED, S.GRID, S.AXIS, S.SURF
BLUE, ORANGE, AQUA = S.BLUE, S.ORANGE, S.AQUA
C = f"{L}/tables/P4/confirm"
H1C = ("human gut", "livestock", "wastewater", "activated sludge")


def save(fig, name):
    for ext in ("pdf", "png"):
        fig.savefig(f"{L}/manuscript/figures/{name}.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(name, "written")


def fig5():
    cl = pd.read_csv(f"{C}/classes.tsv", sep="\t").sort_values("net_w").reset_index(drop=True)
    pj = pd.read_csv(f"{C}/projects.tsv", sep="\t")
    fig = plt.figure(figsize=(7.2, 6.3))
    outer = fig.add_gridspec(2, 1, height_ratios=[1.15, 1], hspace=0.62)
    top = outer[0].subgridspec(1, 2, width_ratios=[1.05, 1], wspace=0.5)
    bot = outer[1].subgridspec(1, 3, width_ratios=[1, 1, 1.75], wspace=0.85)

    def titles(ax, head, sub):
        ax.set_title(head, pad=16)
        ax.text(0, 1.025, sub, transform=ax.transAxes, fontsize=7, color=MUTED, ha="left", va="bottom")

    lab = [e.capitalize() for e in cl.env]

    # a: class estimates with jackknife intervals, and the project estimates behind them
    ax = fig.add_subplot(top[0])
    ax.axvline(0, color=AXIS, lw=0.8)
    for i, r in cl.iterrows():
        p = pj[(pj.env == r.env) & pj.net_w.notna()]
        ax.scatter(p.net_w, np.full(len(p), i) + 0.26, s=6 + 60 * p.units / pj.units.max(), color=BLUE, alpha=0.4, lw=0, zorder=2)
        ax.plot([r.lo95, r.hi95], [i, i], color=INK, lw=2, solid_capstyle="round", zorder=3)
        ax.plot(r.net_w, i, "o", color=INK, ms=6, mec=SURF, mew=1.2, zorder=4)
        ax.text(0.47, i, f"{int(r.projects_with_estimate)} of {int(r.projects)}", ha="right", va="center", fontsize=7, color=MUTED)
    ax.text(0.47, len(cl) - 0.45, "projects with\nan estimate", ha="right", va="bottom", fontsize=6.5, color=MUTED)
    ax.set_yticks(range(len(cl)), lab)
    ax.set_ylim(-0.6, len(cl) + 0.25)
    ax.set_xlim(-0.08, 0.48)
    ax.set_xlabel("ARG–MGE net linked fraction within 2 kb")
    titles(ax, "Class and project estimates", "black: class, 95% interval; blue: projects")
    ax.grid(axis="x")
    ax.tick_params(length=0)
    S.panel(ax, "a")

    # b: rank of the classes by the two measures
    ax = fig.add_subplot(top[1])
    rn = cl.naive_native.rank(ascending=False).astype(int)
    rw = cl.net_w.rank(ascending=False).astype(int)
    for i, r in cl.iterrows():
        moved = abs(rn[i] - rw[i]) >= 2
        col = ORANGE if moved else MUTED
        ax.plot([0, 1], [rn[i], rw[i]], color=col, lw=2 if moved else 1.3, solid_capstyle="round", zorder=3 if moved else 2)
        ax.plot([0, 1], [rn[i], rw[i]], "o", color=col, ms=5, mec=SURF, mew=1, zorder=4)
        ax.text(-0.06, rn[i], f"{lab[i]}  {r.naive_native:.2f}", ha="right", va="center", fontsize=7.5, color=INK if moved else INK2)
        ax.text(1.06, rw[i], f"{r.net_w:.3f}  {lab[i]}", ha="left", va="center", fontsize=7.5, color=INK if moved else INK2)
    ax.set_xlim(-0.95, 1.95)
    ax.set_ylim(7.6, 0.4)
    ax.set_xticks([0, 1], ["Naive share", "Net linked\nfraction"])
    ax.set_yticks([])
    for s in ("left", "bottom"):
        ax.spines[s].set_visible(False)
    ax.tick_params(length=0)
    titles(ax, "Order of the classes by each measure", "orange: moves two or more places")
    S.panel(ax, "b")

    # c: within-project dependence on read N50 (the runs and the centring of H1c)
    pr = pd.read_csv(f"{C}/per_run.tsv", sep="\t")
    pb = pd.read_csv(f"{C}/partB.tsv", sep="\t")
    a = pr[(pr.class_a == "ARG") & (pr.class_b == "MGE")][["run", "env", "project", "net_w"]].merge(
        pb[["run", "intact_naive", "intact_n50", "intact_reads_A"]], on="run")
    a = a[a.env.isin(H1C) & (a.intact_reads_A >= 50)].assign(x=lambda d: np.log2(d.intact_n50), g=lambda d: d.env + "|" + d.project)
    h = pd.read_csv(f"{C}/h1c.tsv", sep="\t").set_index("outcome")
    lim = None
    for j, (y, key, title, col) in enumerate((("intact_naive", "naive share", "Naive share", ORANGE), ("net_w", "net_w", "Net linked fraction", BLUE))):
        d = a[np.isfinite(a[y]) & np.isfinite(a.x)]
        d = d[d.groupby("g")[y].transform("size") >= 2]
        xc, yc = d.x - d.groupby("g").x.transform("mean"), d[y] - d.groupby("g")[y].transform("mean")
        assert len(d) == int(h.loc[key, "runs"]) and abs(float((xc * yc).sum() / (xc ** 2).sum()) - h.loc[key, "slope"]) < 1e-9
        ax = fig.add_subplot(bot[j])
        ax.axhline(0, color=AXIS, lw=0.8)
        ax.scatter(xc, yc, s=11, color=col, alpha=0.55, lw=0)
        xx = np.array([xc.min(), xc.max()])
        ax.plot(xx, h.loc[key, "slope"] * xx, color=INK, lw=1.4)
        lim = lim or (float(np.abs(yc).max()) * 1.12)
        ax.set_ylim(-lim, lim)
        ax.set_xlabel("log2 read N50,\nminus project mean")
        ax.set_ylabel("Value minus project mean" if j == 0 else "")
        titles(ax, title, f"{h.loc[key, 'slope']:+.3f} per doubling of read N50")
        ax.grid(axis="y")
        if j == 0:
            S.panel(ax, "c")

    # d: class estimates under other class definitions
    v = pd.read_csv(f"{L}/tables/P4/sens/sens_variants.classes.tsv", sep="\t")
    ax = fig.add_subplot(bot[2])
    ax.axvline(0, color=AXIS, lw=0.8)
    series = (("primary", "Registered", INK, "o"), ("ARG_MRG_60", "ARG, MRG ≥ 60%", AQUA, "s"), ("MGE_D3", "Broad MGE", ORANGE, "^"))
    for i, r in cl.iterrows():
        w = v[v.env == r.env].set_index("variant").net_w
        ax.plot([w[[s[0] for s in series]].min(), w[[s[0] for s in series]].max()], [i, i], color=GRID, lw=2.2, solid_capstyle="round", zorder=1)
        for name, _, col, mk in series:
            ax.plot(w[name], i, mk, color=col, ms=5.5, mec=SURF, mew=0.9, zorder=3)
    ax.set_yticks(range(len(cl)), lab)
    ax.set_ylim(-0.6, len(cl) - 0.4)
    ax.set_xlabel("ARG–MGE net linked fraction within 2 kb")
    titles(ax, "Other class definitions", "class estimates of panel a under two other definitions")
    ax.grid(axis="x")
    ax.tick_params(length=0)
    ax.legend(handles=[S.Line2D([], [], marker=mk, color=col, ls="", ms=5.5, mec=SURF, label=lb) for _, lb, col, mk in series],
              loc="upper center", bbox_to_anchor=(0.3, -0.27), ncol=3, handletextpad=0.1, columnspacing=0.9)
    ax.set_xlim(-0.015, 0.30)
    S.panel(ax, "d")
    save(fig, "fig5_v2")


def fig6():
    """Reanalysis of a published comparison (a, b: their calls; c: the same runs with this pipeline) and one sample
    with several kits (d). Panels c and d are drawn when their tables exist."""
    import os

    R = f"{L}/tables/P4/reanalysis"
    lc = pd.read_csv(f"{R}/their_length_classes.tsv", sep="\t")
    env = pd.read_csv(f"{R}/their_env.tsv", sep="\t").sort_values("share_pooled", ascending=False).reset_index(drop=True)
    order = ["[0.0, 1000.0)", "[1000.0, 2000.0)", "[2000.0, 4000.0)", "[4000.0, 8000.0)", "[8000.0, 16000.0)", "[16000.0, inf)"]
    ticks = ["<1", "1–2", "2–4", "4–8", "8–16", ">16"]
    short = {"wastewater treatment plants": "Treatment plants", "hospital wastewater": "Hospital wastewater",
             "community sewage": "Community sewage", "municipal wastewater": "Municipal wastewater",
             "swine wastewater": "Swine wastewater", "river": "River", "activated sludge": "Activated sludge", "soil": "Soil"}
    have_c, have_d = os.path.exists(f"{R}/ours.groups.tsv"), os.path.exists(f"{L}/tables/P4/kits/kits.cells.tsv")
    fig = plt.figure(figsize=(7.2, 6.4 if (have_c or have_d) else 3.1))
    outer = fig.add_gridspec(2 if (have_c or have_d) else 1, 1, hspace=0.62)
    top = outer[0].subgridspec(1, 2, width_ratios=[1, 1.2], wspace=0.75)

    def titles(ax, head, sub):
        ax.set_title(head, pad=16)
        ax.text(0, 1.025, sub, transform=ax.transAxes, fontsize=7, color=MUTED, ha="left", va="bottom")

    # a: their share by the length of the ARG read
    ax = fig.add_subplot(top[0])
    for s, g in lc[lc.source != "all"].groupby("source"):
        g = g.set_index("len_class").reindex(order)
        ax.plot(range(6), g.share, color=GRID, lw=1.1, zorder=2)
        if s in ("wastewater treatment plants", "hospital wastewater"):
            nudge = 0.012 if s.startswith("wastewater") else -0.012
            ax.text(5.15, g.share.iloc[-1] + nudge, short[s].split()[0], fontsize=6.5, color=INK2, va="center")
    a = lc[lc.source == "all"].set_index("len_class").reindex(order)
    ax.plot(range(6), a.share, color=ORANGE, lw=2, marker="o", ms=4.5, mec=SURF, mew=0.9, zorder=4)
    ax.text(5.15, a.share.iloc[-1], "All runs", fontsize=7, color=INK, va="center")
    ax.set_xticks(range(6), ticks)
    ax.set_xlim(-0.3, 6.9)
    ax.set_ylim(0, 0.4)
    ax.set_xlabel("Length of the ARG-carrying read (kb)")
    ax.set_ylabel("Share with an MGE call")
    titles(ax, "Published calls, by read length", "grey: the eight environments of panel b")
    ax.grid(axis="y")
    S.panel(ax, "a")

    # b: environments, native reads against reads cut to N50 2 kb
    ax = fig.add_subplot(top[1])
    n = len(env)
    for i, r in env.iterrows():
        y = n - 1 - i
        ax.plot([r.cut_share_pooled, r.share_pooled], [y, y], color=GRID, lw=2.2, solid_capstyle="round", zorder=1)
        ax.plot(r.share_pooled, y, "o", color=ORANGE, ms=5.5, mec=SURF, mew=0.9, zorder=3)
        ax.plot(r.cut_share_pooled, y, "s", color=BLUE, ms=5.2, mec=SURF, mew=0.9, zorder=3)
        ax.text(0.262, y, f"{int(r.rank_native)} → {int(r.rank_cut)}", fontsize=7, color=INK if r.rank_native != r.rank_cut else MUTED,
                ha="right", va="center")
    ax.text(0.262, n - 0.35, "rank", fontsize=6.5, color=MUTED, ha="right", va="bottom")
    lab = [short[s] + (" *" if a_ > 0 else "") for s, a_ in zip(env.source, env.amplified_runs)]
    ax.set_yticks(range(n - 1, -1, -1), lab)
    ax.set_ylim(-0.6, n - 0.1)
    ax.set_xlim(-0.008, 0.27)
    ax.set_xlabel("Share of ARG reads with an MGE call")
    titles(ax, "The same calls after a cut to N50 2 kb", "* amplified libraries")
    ax.grid(axis="x")
    ax.tick_params(length=0)
    ax.legend(handles=[S.Line2D([], [], marker="o", color=ORANGE, ls="", ms=5.5, mec=SURF, label="Native reads"),
                       S.Line2D([], [], marker="s", color=BLUE, ls="", ms=5.2, mec=SURF, label="Cut reads")],
              loc="lower right", bbox_to_anchor=(0.86, 0.02), handletextpad=0.1)
    S.panel(ax, "b")

    bot = outer[1].subgridspec(1, 3, width_ratios=[1.75, 1, 1], wspace=0.85) if (have_c or have_d) else None
    # c: the same runs with this pipeline: naive share and net linked fraction (rows as in panel b)
    if have_c:
        cmp_ = pd.read_csv(f"{R}/compare.env.tsv", sep="\t").set_index("source").reindex(env.source)
        ax = fig.add_subplot(bot[0])
        ax.axvline(0, color=AXIS, lw=0.8)
        for i, (s, r) in enumerate(cmp_.iterrows()):
            y = n - 1 - i
            ax.plot(r.naive, y + 0.17, "o", color=ORANGE, ms=5.5, mec=SURF, mew=0.9, zorder=3)
            ax.plot([r.net_w_lo, r.net_w_hi], [y - 0.17, y - 0.17], color=INK, lw=2, solid_capstyle="round", zorder=2)
            ax.plot(r.net_w, y - 0.17, "o", color=INK, ms=5.5, mec=SURF, mew=1.1, zorder=4)
            ax.text(0.585, y, f"{int(r.rank_their_share)} → {int(r.rank_net_w)}", fontsize=7, ha="right", va="center",
                    color=INK if abs(r.rank_their_share - r.rank_net_w) >= 2 else MUTED)
        ax.text(0.585, n - 0.35, "rank", fontsize=6.5, color=MUTED, ha="right", va="bottom")
        ax.set_yticks(range(n - 1, -1, -1), lab)
        ax.set_ylim(-0.6, n - 0.1)
        ax.set_xlim(-0.08, 0.6)
        ax.set_xlabel("Share or fraction, annotation of this study")
        titles(ax, "The same runs, annotated here", "rank: published share → net linked fraction")
        ax.grid(axis="x")
        ax.tick_params(length=0)
        ax.legend(handles=[S.Line2D([], [], marker="o", color=ORANGE, ls="", ms=5.5, mec=SURF, label="Naive share"),
                           S.Line2D([], [], marker="o", color=INK, ls="-", lw=2, ms=5.5, mec=SURF, label="Net linked fraction")],
                  loc="upper center", bbox_to_anchor=(0.4, -0.27), ncol=2, handletextpad=0.4, columnspacing=1.2)
        S.panel(ax, "c")

    # d: one sample, several kits: values against read N50, both centred on the sample
    if have_d:
        K = f"{L}/tables/P4/kits"
        cells = pd.read_csv(f"{K}/kits.cells.tsv", sep="\t")
        sl = pd.read_csv(f"{K}/kits.slopes.tsv", sep="\t")
        sl = sl[(sl.pair == "ARG:MGE") & (sl.set == "all sets")].set_index("outcome")
        sets = {"stool_kits": ("Stool, 2 kits", BLUE, "o"), "rumen_kits": ("Rumen, kits", ORANGE, "s"),
                "rumen_extraction": ("Rumen, extraction", AQUA, "^")}
        lim = None
        for j, (y, key, title) in enumerate((("ARG_MGE_naive", "naive share", "Naive share"),
                                             ("ARG_MGE_net_w", "net linked fraction within 2 kb", "Net linked fraction"))):
            ax = fig.add_subplot(bot[1 + j])
            ax.axhline(0, color=AXIS, lw=0.8)
            d = cells.dropna(subset=[f"{y}_xc", f"{y}_yc"])
            for s, (name, col, mk) in sets.items():
                g = d[d["set"] == s]
                ax.scatter(g[f"{y}_xc"], g[f"{y}_yc"], s=13, color=col, marker=mk, alpha=0.7, lw=0, label=name)
            if key in sl.index and np.isfinite(sl.loc[key, "slope"]):
                xx = np.array([d[f"{y}_xc"].min(), d[f"{y}_xc"].max()])
                ax.plot(xx, sl.loc[key, "slope"] * xx, color=INK, lw=1.4)
                sub = f"{sl.loc[key, 'slope']:+.3f} per doubling of read N50"
            else:
                sub = "no estimate"
            lim = lim or (float(np.abs(d[f"{y}_yc"]).max()) * 1.12 if len(d) else 0.1)
            ax.set_ylim(-lim, lim)
            ax.set_xlabel("log2 read N50,\nminus sample mean")
            ax.set_ylabel("Value minus sample mean" if j == 0 else "")
            titles(ax, title, sub)
            ax.grid(axis="y")
            if j == 0:
                S.panel(ax, "d")
                ax.legend(loc="upper center", bbox_to_anchor=(1.3, -0.42), ncol=3, handletextpad=0.1, columnspacing=0.9)
    save(fig, "fig6_v2")


if __name__ == "__main__":
    for k in [int(x) for x in sys.argv[1:]] or [5]:
        {5: fig5, 6: fig6}[k]()
