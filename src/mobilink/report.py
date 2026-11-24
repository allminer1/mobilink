"""The HTML report of `mobilink run`: one self-contained file (no script, no external resource)."""

from __future__ import annotations

import html
import math

import pandas as pd

MIN_CURVE_SIDES = 100
SERIES = ("#2a78d6", "#d9642a", "#3d9a6b", "#8a5cc2")  # fixed order: a pair keeps its colour
CSS = """
:root{--bg:#fff;--fg:#1c2330;--mute:#5b6677;--line:#dfe3ea;--card:#f6f7f9;--warn:#fff4d6;--warnline:#e2b93b}
@media (prefers-color-scheme:dark){:root{--bg:#14171c;--fg:#e7eaf0;--mute:#9aa4b4;--line:#2c323c;--card:#1b1f26;
--warn:#3a3115;--warnline:#a8882a}}
*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--fg);
font:15px/1.5 system-ui,-apple-system,"Segoe UI",Roboto,sans-serif}
main{max-width:980px;margin:0 auto;padding:24px 16px 48px}h1{font-size:24px;margin:0 0 4px}
h2{font-size:18px;margin:32px 0 8px}p{margin:8px 0}.mute{color:var(--mute)}
.tiles{display:grid;grid-template-columns:repeat(auto-fit,minmax(150px,1fr));gap:8px;margin:12px 0}
.tile{background:var(--card);border-radius:8px;padding:10px 12px}.tile b{display:block;font-size:20px}
.tile span{color:var(--mute);font-size:13px}
.scroll{overflow-x:auto}table{border-collapse:collapse;width:100%;font-size:14px}
th,td{text-align:left;padding:6px 10px;border-bottom:1px solid var(--line);white-space:nowrap}
th{color:var(--mute);font-weight:600}td.num,th.num{text-align:right;font-variant-numeric:tabular-nums}
.warn{background:var(--warn);border-left:3px solid var(--warnline);padding:8px 12px;margin:8px 0;border-radius:4px}
.verdict{font-size:16px;margin:6px 0}svg text{fill:var(--fg);font-size:12px}svg .mute{fill:var(--mute)}
svg .grid{stroke:var(--line);stroke-width:1}svg .zero{stroke:var(--mute);stroke-width:1;stroke-dasharray:4 3}
.legend{display:flex;gap:16px;flex-wrap:wrap;font-size:13px;margin:4px 0 0}
.legend i{display:inline-block;width:14px;height:3px;border-radius:2px;margin-right:6px;vertical-align:middle}
code{background:var(--card);padding:1px 4px;border-radius:4px}
"""


def _f(x, nd=3, dash="–"):
    try:
        if x is None or (isinstance(x, float) and math.isnan(x)):
            return dash
        return f"{float(x):.{nd}f}"
    except (TypeError, ValueError):
        return dash


def _p(x):
    try:
        x = float(x)
    except (TypeError, ValueError):
        return "–"
    if math.isnan(x):
        return "–"
    return "< 0.001" if x < 0.001 else f"{x:.3f}"


def _peq(x):
    t = _p(x)
    return f"p {t}" if t.startswith("<") else f"p = {t}"


def _n(x):
    try:
        return f"{int(x):,}"
    except (TypeError, ValueError):
        return "–"


def _table(head, rows, num=()):
    th = "".join(f'<th class="{"num" if i in num else ""}">{html.escape(h)}</th>' for i, h in enumerate(head))
    body = "".join(
        "<tr>" + "".join(f'<td class="{"num" if i in num else ""}">{html.escape(str(c))}</td>'
                         for i, c in enumerate(r)) + "</tr>" for r in rows)
    return f'<div class="scroll"><table><thead><tr>{th}</tr></thead><tbody>{body}</tbody></table></div>'


def verdict(row, w: int) -> str:
    """One plain sentence for a row of the test table."""
    a, b = row["class_a"], row["class_b"]
    net, p = row.get("net_w"), row.get("net_w_p")
    if net is None or pd.isna(net):
        why = str(row.get("flags") or "too few observations")
        return f"{a} and {b}: no estimate ({why})."
    lo, hi = row.get("net_w_lo"), row.get("net_w_hi")
    rng = f" (95% interval {_f(lo)} to {_f(hi)})" if not pd.isna(lo) and not pd.isna(hi) else ""
    pct = f"{100 * float(net):.1f}%"
    if pd.isna(p):
        return f"{a} and {b}: net linked fraction {_f(net)}{rng}; no test was possible."
    if float(p) <= 0.05 and float(net) > 0:
        wide = ""
        if not pd.isna(lo) and float(lo) <= 0:
            wide = (" The test and the interval use different methods; the interval is the more cautious one and "
                    "includes 0 here, so treat the size as uncertain.")
        return (f"{a} and {b} are linked: {pct} of the flanks of {a} units have a whole {b} gene within "
                f"{w / 1000:g} kb, beyond what chance gives{rng}; {_peq(p)}.{wide}")
    if float(p) <= 0.05:
        return (f"{a} and {b} sit apart: fewer flanks of {a} units have a {b} gene within {w / 1000:g} kb than "
                f"chance gives (net linked fraction {_f(net)}{rng}; {_peq(p)}).")
    return (f"{a} and {b}: no linkage shown (net linked fraction {_f(net)}{rng}; {_peq(p)}). This is not "
            "evidence of no linkage unless the interval is narrow.")


def curve_svg(cv: pd.DataFrame, w: int) -> str:
    """Net linked fraction against the window, one line per pair; the dashed line at 0 is chance."""
    cv = cv.dropna(subset=["net_w"])
    cv = cv[cv.sides_w >= MIN_CURVE_SIDES]  # windows that few molecules reach are noise; they stay in the table
    if cv.empty:
        return '<p class="mute">Too few observations for a distance curve.</p>'
    W, H, L, R, T, B = 720, 300, 56, 120, 16, 40
    xmax = float(cv.w.max())
    ymin, ymax = min(0.0, float(cv.net_w.min())), max(0.05, float(cv.net_w.max()))
    pad = 0.08 * (ymax - ymin)
    ymin, ymax = ymin - pad, ymax + pad

    def x(v):
        return L + (W - L - R) * v / xmax

    def y(v):
        return T + (H - T - B) * (ymax - v) / (ymax - ymin)

    out = [f'<svg viewBox="0 0 {W} {H}" role="img" aria-label="Net linked fraction against distance" '
           f'style="width:100%;max-width:{W}px;height:auto">']
    step = 10 ** math.floor(math.log10(ymax - ymin))
    step = step / 2 if (ymax - ymin) / step < 4 else step
    t = math.ceil(ymin / step) * step
    while t <= ymax:
        out.append(f'<line class="grid" x1="{L}" x2="{W - R}" y1="{y(t):.1f}" y2="{y(t):.1f}"/>'
                   f'<text class="mute" x="{L - 8}" y="{y(t) + 4:.1f}" text-anchor="end">{t:.2f}</text>')
        t += step
    for v in sorted(set(cv.w)):
        if v in (1000, 2000, 5000, 10000) or v == xmax:
            out.append(f'<text class="mute" x="{x(v):.1f}" y="{H - B + 18}" text-anchor="middle">{v / 1000:g} kb'
                       "</text>")
    out.append(f'<line class="zero" x1="{L}" x2="{W - R}" y1="{y(0):.1f}" y2="{y(0):.1f}"/>'
               f'<text class="mute" x="{W - R + 6}" y="{y(0) + 4:.1f}">chance</text>')
    if w <= xmax:
        out.append(f'<line class="grid" x1="{x(w):.1f}" x2="{x(w):.1f}" y1="{T}" y2="{H - B}"/>')
    out.append(f'<text class="mute" x="{(L + W - R) / 2:.0f}" y="{H - 4}" text-anchor="middle">'
               "distance from the gene (window)</text>")
    legend, ends = [], []
    for i, ((a, b), g) in enumerate(cv.groupby(["class_a", "class_b"], sort=False)):
        col = SERIES[i % len(SERIES)]
        g = g.sort_values("w")
        pts = " ".join(f"{x(r.w):.1f},{y(r.net_w):.1f}" for r in g.itertuples())
        out.append(f'<polyline fill="none" stroke="{col}" stroke-width="2" points="{pts}"/>')
        for r in g.itertuples():
            out.append(f'<circle cx="{x(r.w):.1f}" cy="{y(r.net_w):.1f}" r="4" fill="{col}" stroke="var(--bg)" '
                       f'stroke-width="2"><title>{html.escape(a)}–{html.escape(b)}, within {r.w / 1000:g} kb: '
                       f"{r.net_w:.3f} ({int(r.events_w)} of {int(r.sides_w)} sides)</title></circle>")
        last = g.iloc[-1]
        ends.append([y(last.net_w) + 4, x(last.w) + 8, f"{html.escape(a)}–{html.escape(b)}"])
        legend.append(f'<span><i style="background:{col}"></i>{html.escape(a)}–{html.escape(b)}</span>')
    ends.sort()
    for k in range(1, len(ends)):  # labels at the line ends keep 14 px apart
        ends[k][0] = max(ends[k][0], ends[k - 1][0] + 14)
    out += [f'<text x="{ex:.1f}" y="{ey:.1f}">{label}</text>' for ey, ex, label in ends]
    out.append("</svg>")
    return "".join(out) + f'<div class="legend">{"".join(legend)}</div>'


def write_report(path, name, info, test, sens, cv, genes, fixed, w: int = 2000, contigs: bool = False,
                 clean=None) -> None:
    e = html.escape
    kind = "contigs" if contigs else "reads"
    inp, calls = info.get("input", {}), info.get("calls", {})
    parts = [f"<!doctype html><html lang='en'><head><meta charset='utf-8'><meta name='viewport' "
             f"content='width=device-width,initial-scale=1'><title>mobilink: {e(name)}</title><style>{CSS}</style>"
             f"</head><body><main><h1>{e(name)}</h1><p class='mute'>mobilink {e(str(info.get('mobilink')))} · "
             f"{e(str(info.get('finished', ''))[:10])} · input: {kind}</p>"]

    parts.append("<h2>Result</h2>")
    for _, r in test.iterrows():
        parts.append(f"<p class='verdict'>{e(verdict(r, w))}</p>")
    for fl in info.get("flags", []):
        parts.append(f"<div class='warn'>{e(fl)}</div>")
    rows = []
    for _, r in test.iterrows():
        rows.append([f"{r['class_a']}–{r['class_b']}", _f(r.get("net_w")),
                     f"{_f(r.get('net_w_lo'))} to {_f(r.get('net_w_hi'))}", _p(r.get("net_w_p")),
                     _p(r.get("net_w_q")), _n(r.get("units")), _n(r.get("sides_w")),
                     f"{_f(r.get('net_w_O'), 0)} / {_f(r.get('net_w_E'), 1)}", _f(r.get("net_mol")),
                     str(r.get("flags") if isinstance(r.get("flags"), str) else "")])
    parts.append(_table(["Pair", f"Net linked within {w / 1000:g} kb", "95% interval", "p", "q", "Units",
                         "Sides reaching the window", "Observed / expected", "Net linked, same molecule", "Notes"],
                        rows, num=(1, 2, 3, 4, 5, 6, 7, 8)))
    parts.append("<p class='mute'>Each gene unit of the first class has two flanks, left and right. The net linked "
                 "fraction is the share of flanks with a whole gene of the second class within the window, beyond "
                 "what the gene densities give by chance. A unit with such a gene on one flank only counts one "
                 "half. The value is 0 under chance and does not depend on read length. Compare it between samples; "
                 "do not compare the raw share of reads that carry both.</p>")
    if clean is not None and not clean.empty:
        lp = info.get("laboratory_plasmid", {})
        rows = [[f"{r['class_a']}–{r['class_b']}", _f(r.get("net_w")),
                 f"{_f(r.get('net_w_lo'))} to {_f(r.get('net_w_hi'))}", _p(r.get("net_w_p")),
                 _n(r.get("excluded_molecules"))] for _, r in clean.iterrows()]
        parts.append("<h2>Without laboratory plasmids</h2>" + _table(
            ["Pair", "Net linked", "95% interval", "p", "Reads removed"], rows, num=(1, 2, 3, 4)))
        parts.append(f"<p class='mute'>All reads that carry {e(', '.join(lp.get('genes', [])))} were removed, because "
                     "its carrier reads match cloning-vector sequences. Report these values, not the ones above.</p>")
    if sens is not None and not sens.empty:
        rows = [[f"{r['class_a']}–{r['class_b']}", _f(r.get("net_w")),
                 f"{_f(r.get('net_w_lo'))} to {_f(r.get('net_w_hi'))}", _p(r.get("net_w_p")),
                 _n(r.get("excluded_molecules"))] for _, r in sens.iterrows()]
        parts.append("<h2>Without reads flagged as chimeras</h2>" + _table(
            ["Pair", "Net linked", "95% interval", "p", "Reads removed"], rows, num=(1, 2, 3, 4)))
        parts.append("<p class='mute'>A read that joins two unrelated molecules can create a false link. If the "
                     "values above differ clearly from the main result, report both.</p>")

    parts.append("<h2>How linkage builds up with distance</h2>" + curve_svg(cv, w))
    parts.append("<p class='mute'>Each point is the net linked fraction for one window. A curve that rises and then "
                 "flattens shows the distance over which the two classes travel together. The vertical line is the "
                 f"window of the test. Windows that fewer than {MIN_CURVE_SIDES} flanks reach are not drawn.</p>")

    parts.append("<h2>Input</h2><div class='tiles'>")
    tiles = [(_n(inp.get("molecules")), kind), (f"{inp.get('bases', 0) / 1e9:.2f} Gbp", "bases"),
             (f"{inp.get('n50', 0) / 1000:.1f} kb", "N50")]
    tiles += [(_n(calls.get(c, 0)), f"{c} calls") for c in ("ARG", "MRG", "MGE")]
    if "host_unresolved_share" in info:
        tiles.append((f"{100 * info['host_unresolved_share']:.0f}%", f"{kind} without a host label"))
    if "chimeras" in info:
        tiles.append((_n(info["chimeras"]["flagged"]), f"of {_n(info['chimeras']['judged'])} judged reads flagged "
                                                       "as chimeras"))
    parts += [f"<div class='tile'><b>{e(v)}</b><span>{e(k)}</span></div>" for v, k in tiles]
    parts.append("</div>")
    if info.get("host_unresolved_share", 0) > 0.5:
        parts.append("<div class='warn'>More than half of the molecules have no host label. The host control is then "
                     "weak: genes of one abundant host can look linked. This is usual outside the human gut.</div>")

    if genes is not None and not genes.empty:
        parts.append("<h2>Which genes sit next to which</h2>")
        for (a, b), g in genes.groupby(["class_a", "class_b"], sort=False):
            g = g[g.calls_with_b > 0].head(15)
            if g.empty:
                continue
            extra = [c for c in ("mechanism", "drug_class") if c in g and g[c].notna().any()]
            kind = b == "MGE"  # the label of an MGE call is its mobileOG category; other labels are database ids
            rows = [[r.gene, _n(r.calls), _n(r.calls_with_b), f"{100 * r.share_with_b:.0f}%", _f(r.net_w)]
                    + ([str(r.top_b_label if isinstance(r.top_b_label, str) else "–")] if kind else [])
                    + [str(r.top_b_gene if isinstance(r.top_b_gene, str) else "–")]
                    + [str(getattr(r, c)) if isinstance(getattr(r, c), str) else "–" for c in extra]
                    for r in g.itertuples()]
            parts.append(f"<p><b>{e(a)} next to {e(b)}</b> (the 15 genes with most calls next to {e(b)})</p>" + _table(
                ["Gene", "Calls", f"With {b} within {w / 1000:g} kb", "Share", "Net linked (gene alone)"]
                + ([f"Commonest {b} kind"] if kind else []) + [f"Commonest {b} gene"]
                + [c.replace("_", " ").capitalize() for c in extra],
                rows, num=(1, 2, 3, 4)))
        parts.append("<p class='mute'>Descriptive. The share is not corrected for chance or read length; the "
                     "net linked value is given only for genes with at least 30 calls and has no interval. The full "
                     "table is in the genes file.</p>")

    if fixed is not None and not fixed.empty and fixed.flag.any():
        checked = "laboratory_plasmid" in fixed.columns
        rows = [[r.gene, _n(r.carriers), f"{r.median_length / 1000:.2f} kb", f"{100 * r.share_near_median:.0f}%"]
                + ([f"{100 * r.median_covered:.0f}%", f"{r.median_sources:g}",
                    "laboratory plasmid" if str(r.laboratory_plasmid).lower() == "true" else "not a vector"]
                   if checked else []) for r in fixed[fixed.flag].itertuples()]
        parts.append("<h2>Genes on reads of one fixed length</h2>" + _table(
            ["Gene", "Carrier reads", "Median read length", "Within 3% of the median"]
            + (["Read matching UniVec", "UniVec records", "Verdict"] if checked else []), rows, num=(1, 2, 3, 4, 5)))

    parts.append("<h2>Method and sources</h2><p>Genes were called with DIAMOND blastx (long-read mode) against "
                 "CARD (ARG, identity at least 80%), BacMet2 experimentally confirmed (MRG, at least 80%) and the "
                 "mobility core of mobileOG-db (MGE, at least 60%), each with at least 80% of the reference protein "
                 "covered. " + ("Host labels come from Kraken2 at species rank. " if "host_unresolved_share" in info
                                else "No host labels were used. ")
                 + "The test compares the observed number of linked gene units with random positions in the same "
                 "hosts, and counts reads that cover the same locus once.</p>")
    rows = [[k, str(v.get("version")), str(v.get("cite")), str(v.get("licence"))]
            for k, v in info.get("databases", {}).items()]
    parts.append(_table(["Database", "Version", "Cite", "Licence"], rows))
    files = ", ".join(f"<code>{e(f)}</code>" for f in info.get("files", []))
    parts.append(f"<p class='mute'>Files of this run: {files}</p></main></body></html>")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("".join(parts))
