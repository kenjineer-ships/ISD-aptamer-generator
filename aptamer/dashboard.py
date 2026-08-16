"""Render the pipeline's outputs as one self-contained dashboard.html at the repo root.

Reads whatever the pipeline has written so far (parents.json, switches.csv, *.pdb / *.cif)
and inlines it into a single HTML file as JSON. Nothing is fetched at runtime: the dashboard
is opened over file://, where Chrome blocks local XHR/fetch as cross-origin.

    python aptamer/dashboard.py        -> ./dashboard.html

Missing inputs degrade to a "not yet generated" card; they never raise. No science is
recomputed here -- dg_switch, KD_app, structures and MFEs are read as produced. The only
derived values are display conveniences (residence time 1/koff, fractional occupancy from a
published K_D), and those come from build_parents so the formulas live in one place.

Charts are inline SVG built by ~200 lines of vanilla JS from the inlined JSON. Two panels
prefer a CDN script for something SVG-by-hand cannot do: 3Dmol.js for the 3D structures
(view 8) and fornac for the secondary-structure layout (view 7). View 7 degrades to the
built-in base-pair arc diagram when fornac is unreachable, so only view 8 actually needs
internet; everything else works with no network at all.

fornac ships its stylesheet inside the bundle and style-loader injects it into <head>
unscoped, where its bare `svg {}` / `text {}` rules would resize every other chart. A
scoped copy is inlined from vendor/fornac.css by scoped_fornac_css() and the injected
copy is deleted at runtime -- see dropFornacGlobalCss() in the JS.
"""

import csv
import json
import pathlib
import sys

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parent
OUT = ROOT / "dashboard.html"

# --- pipeline constants: import them, never re-type them ----------------------------------
sys.path.insert(0, str(HERE))

IL6_MW = 23718  # fallback only; overwritten by the import below
WOUND = {
    "post-surgical (median)": 986.6,
    "inflammatory (median)": 4964,
    "infected (median)": 5883,
    "max observed": 135500,
}
SKIN_MODEL_PG_ML = 5000.0  # 5 ng/mL, README "Known constraint: the affinity gap"
NOTES = []


def _occupancy_fallback(kd, pg_ml):
    conc = pg_ml * 1e-9 / IL6_MW
    return conc, conc / (kd + conc)


occupancy = _occupancy_fallback
try:
    from build_parents import IL6_MW as _MW, WOUND as _WOUND, occupancy as _occ

    IL6_MW, WOUND, occupancy = _MW, _WOUND, _occ
except Exception as exc:  # ViennaRNA missing, etc. -- display values still work
    NOTES.append(
        f"build_parents not importable ({exc.__class__.__name__}: {exc}); using local constants")

SWITCH_PARENT = "IL-6-7326"
RANDOMISED = sorted({12, 13, 18, 19, 20, 25, 26, 27, 45, 46, 47, 53, 54, 55, 60, 61})
# The design target as switch_library states it. Older revisions exposed it as a literal
# DG_WINDOW; newer ones state the affinity budget it is derived from. Take whichever exists
# and never re-derive it here -- the shaded band falls back to the passing data range.
DG_WINDOW = None
BUDGET = {}
try:
    import switch_library as _sl

    SWITCH_PARENT = _sl.PARENT
    RANDOMISED = sorted(_sl.randomised_positions())
    if hasattr(_sl, "DG_WINDOW"):
        DG_WINDOW = tuple(_sl.DG_WINDOW)
    for k in ("KD_APP_MAX_NM", "CLOSED_MIN", "MIN_ENGAGEMENT"):
        if hasattr(_sl, k):
            BUDGET[k] = getattr(_sl, k)
except Exception as exc:
    NOTES.append(
        f"switch_library not importable ({exc.__class__.__name__}: {exc}); "
        "parent name and randomised block positions taken from the cached values"
    )


# --- readers -----------------------------------------------------------------------------
def read_parents():
    p = HERE / "parents.json"
    if not p.exists():
        return None
    try:
        rows = json.loads(p.read_text(encoding="utf-8"))
    except Exception as exc:
        NOTES.append(f"parents.json unreadable: {exc}")
        return None
    for r in rows:
        koff = r.get("koff_s") or 0
        r["residence_s"] = (1.0 / koff) if koff else None
        r["kd_nM"] = (r.get("KD_M") or 0) * 1e9
    return rows


NUMERIC = {
    "ds_len": int, "linker_len": int, "rand_covered": int,
    "gc": float, "dg_switch": float, "closed_frac": float,
    "kd_app_nM": float, "engagement": float,
}


def read_switches():
    p = HERE / "switches.csv"
    if not p.exists():
        return None
    try:
        with p.open(newline="", encoding="utf-8") as fh:
            rows = list(csv.DictReader(fh))
    except Exception as exc:
        NOTES.append(f"switches.csv unreadable: {exc}")
        return None
    out = []
    for i, r in enumerate(rows):
        try:
            for k, cast in NUMERIC.items():
                r[k] = cast(r[k])
            start, end = r["window"].split("-")
            r["w0"], r["w1"] = int(start), int(end)
            r["rank"] = i + 1  # file order IS the pipeline ranking
            out.append(r)
        except Exception:
            continue  # a malformed row is dropped, not fatal
    if len(out) != len(rows):
        NOTES.append(f"switches.csv: {len(rows) - len(out)} malformed row(s) skipped")
    return out


def read_structures():
    """1ALU plus any co-folded complexes the structure step has dropped in."""
    out = []
    for p in sorted(HERE.glob("*.pdb")) + sorted(HERE.glob("*.cif")):
        try:
            text = p.read_text(encoding="utf-8", errors="replace")
        except Exception as exc:
            NOTES.append(f"{p.name} unreadable: {exc}")
            continue
        is_ref = p.stem.upper() == "1ALU"
        out.append({
            "name": p.name,
            "format": "pdb" if p.suffix.lower() == ".pdb" else "cif",
            "kind": "receptor" if is_ref else "complex",
            "label": "human IL-6 crystal structure, 1.9 A (PDB 1ALU)" if is_ref else p.stem,
            "text": text,
            "atoms": sum(1 for line in text.splitlines() if line.startswith(("ATOM", "HETATM"))),
        })
    return out


def occupancy_series(parents):
    """Fraction bound vs IL-6 (pg/mL), log-spaced. Display derivation only."""
    if not parents:
        return []
    series = []
    for p in parents:
        kd = p.get("KD_M")
        if not kd:
            continue
        pts = []
        for i in range(121):
            pg = 10.0 ** (1 + 6 * i / 120)  # 1e1 .. 1e7 pg/mL
            _, frac = occupancy(kd, pg)
            pts.append([pg, frac])
        half = kd * IL6_MW / 1e-9  # pg/mL giving 50% occupancy == K_D in mass units
        series.append({"name": p["name"], "kd_nM": p["kd_nM"], "pts": pts, "half_pg": half})
    return series


def occupancy_table(parents):
    rows = []
    landmarks = list(WOUND.items()) + [("skin model (5 ng/mL)", SKIN_MODEL_PG_ML)]
    for label, pg in landmarks:
        row = {"label": label, "pg_ml": pg, "pM": pg * 1e-9 / IL6_MW * 1e12, "frac": {}}
        for p in parents or []:
            if p.get("KD_M"):
                row["frac"][p["name"]] = occupancy(p["KD_M"], pg)[1]
        rows.append(row)
    return rows


def design_window(switches):
    """(lo, hi, caption) for the shaded ΔG_switch band. Stated target if the pipeline
    exposes one, otherwise the range that actually passed -- never re-derived here."""
    if DG_WINDOW:
        return DG_WINDOW[0], DG_WINDOW[1], (
            f"design window {DG_WINDOW[0]} to {DG_WINDOW[1]} kcal/mol")
    if not switches:
        return -2.2, -0.7, "design window (not stated by the pipeline)"
    dgs = [r["dg_switch"] for r in switches]
    bits = []
    if "KD_APP_MAX_NM" in BUDGET:
        bits.append(f"K_D,app ≤ {BUDGET['KD_APP_MAX_NM']:.0f} nM")
    if "CLOSED_MIN" in BUDGET:
        bits.append(f"closed ≥ {BUDGET['CLOSED_MIN']:.0%}")
    cap = "ΔG range that passed" + (" (" + ", ".join(bits) + ")" if bits else "")
    return min(dgs), max(dgs), cap


def build_data():
    parents = read_parents()
    switches = read_switches()
    structures = read_structures()
    lo, hi, cap = design_window(switches)
    return {
        "generated": None,
        "parents": parents,
        "switches": switches,
        "structures": structures,
        "occupancy": occupancy_series(parents),
        "occupancyTable": occupancy_table(parents),
        "wound": WOUND,
        "skinModel": SKIN_MODEL_PG_ML,
        "il6MW": IL6_MW,
        "dgWindow": [lo, hi],
        "dgWindowCaption": cap,
        "budget": BUDGET,
        "randomised": RANDOMISED,
        "switchParent": SWITCH_PARENT,
        "parentSeq": next((p["sequence"] for p in (parents or [])
                           if p["name"] == SWITCH_PARENT), None),
        "notes": NOTES,
    }


# --- the page ----------------------------------------------------------------------------
CSS = r"""
:root { color-scheme: light; }
.viz-root {
  --surface-1:#fcfcfb; --page:#f9f9f7;
  --text-primary:#0b0b0b; --text-secondary:#52514e; --muted:#898781;
  --grid:#e1e0d9; --axis:#c3c2b7; --border:rgba(11,11,11,0.10);
  --series-1:#2a78d6; --series-2:#eb6834; --series-3:#1baf7a; --series-4:#eda100;
  --warning:#fab219; --critical:#d03b3b; --good:#0ca30c;
  --wash:rgba(42,120,214,0.10);
}
@media (prefers-color-scheme: dark) {
  :root:where(:not([data-theme="light"])) .viz-root {
    color-scheme: dark;
    --surface-1:#1a1a19; --page:#0d0d0d;
    --text-primary:#ffffff; --text-secondary:#c3c2b7; --muted:#898781;
    --grid:#2c2c2a; --axis:#383835; --border:rgba(255,255,255,0.10);
    --series-1:#3987e5; --series-2:#d95926; --series-3:#199e70; --series-4:#c98500;
    --wash:rgba(57,135,229,0.14);
  }
}
:root[data-theme="dark"] .viz-root {
  color-scheme: dark;
  --surface-1:#1a1a19; --page:#0d0d0d;
  --text-primary:#ffffff; --text-secondary:#c3c2b7; --muted:#898781;
  --grid:#2c2c2a; --axis:#383835; --border:rgba(255,255,255,0.10);
  --series-1:#3987e5; --series-2:#d95926; --series-3:#199e70; --series-4:#c98500;
  --wash:rgba(57,135,229,0.14);
}

* { box-sizing:border-box; }
html, body { margin:0; }
body { background:var(--page); }
.viz-root {
  font-family: system-ui, -apple-system, "Segoe UI", sans-serif;
  background:var(--page); color:var(--text-primary);
  padding:24px clamp(16px,4vw,48px) 96px; min-height:100vh;
}
h1 { font-size:22px; font-weight:600; margin:0; letter-spacing:-0.01em; }
h2 { font-size:15px; font-weight:600; margin:0 0 2px; }
p.sub, .sub { color:var(--text-secondary); font-size:13px; margin:0; }
.muted { color:var(--muted); }
a { color:var(--series-1); }

header.top { display:flex; align-items:flex-start; gap:16px; flex-wrap:wrap;
  padding-bottom:16px; border-bottom:1px solid var(--border); margin-bottom:20px; }
header.top .grow { flex:1 1 320px; }
button, select, input[type=number], input[type=text] {
  font:inherit; font-size:13px; color:var(--text-primary); background:var(--surface-1);
  border:1px solid var(--border); border-radius:8px; padding:6px 10px;
}
button { cursor:pointer; }
button:hover, select:hover { border-color:var(--axis); }

.card { background:var(--surface-1); border:1px solid var(--border); border-radius:12px;
  padding:16px 18px 18px; margin:0 0 20px; }
.card > header { margin-bottom:12px; }
.grid2 { display:grid; grid-template-columns:repeat(auto-fit,minmax(420px,1fr)); gap:20px; }
.stats { display:grid; grid-template-columns:repeat(auto-fit,minmax(180px,1fr)); gap:12px;
  margin-bottom:20px; }
.stat { background:var(--surface-1); border:1px solid var(--border); border-radius:12px;
  padding:14px 16px; }
.stat .label { font-size:12px; color:var(--text-secondary); }
.stat .value { font-size:26px; font-weight:600; line-height:1.2; margin-top:2px; }
.stat .value.hero { font-size:48px; }
.stat .unit { font-size:13px; color:var(--text-secondary); font-weight:400; margin-left:4px; }

table { border-collapse:separate; border-spacing:0; width:100%; font-size:13px;
  font-variant-numeric:tabular-nums; }
th, td { text-align:right; padding:7px 10px; white-space:nowrap;
  border-bottom:1px solid var(--grid); }
th:first-child, td:first-child, th.l, td.l { text-align:left; }
thead th { position:sticky; top:0; background:var(--surface-1); z-index:2;
  color:var(--text-secondary); font-weight:600; border-bottom:1px solid var(--axis);
  font-size:12px; }
thead th.sortable { cursor:pointer; user-select:none; }
thead th.sortable:hover { color:var(--text-primary); }
thead th .arrow { color:var(--muted); font-size:10px; }
tbody tr:hover { background:var(--wash); }
.scroll { max-height:460px; overflow:auto; border:1px solid var(--border);
  border-radius:10px; }
code, .mono { font-family: ui-monospace, "Cascadia Mono", Consolas, monospace; font-size:12px; }
.seq { font-family: ui-monospace, "Cascadia Mono", Consolas, monospace; font-size:11.5px;
  letter-spacing:0.055em; word-break:break-all; }

.badge { display:inline-flex; align-items:center; gap:5px; font-size:11.5px; font-weight:600;
  border-radius:999px; padding:2px 9px 2px 7px; border:1px solid var(--border);
  color:var(--text-primary); }
.badge.warn { background:rgba(250,178,25,0.18); }
.badge.ok { background:rgba(12,163,12,0.14); color:var(--text-secondary); font-weight:500; }
.badge .ic { font-weight:700; }
.issues { margin:10px 0 0; padding:10px 12px; border-radius:10px;
  background:rgba(250,178,25,0.14); border:1px solid rgba(250,178,25,0.45); font-size:12.5px; }
.issues b { font-weight:600; }
.issues ul { margin:6px 0 0 18px; padding:0; }

.filters { display:flex; flex-wrap:wrap; gap:14px 18px; align-items:flex-end;
  background:var(--surface-1); border:1px solid var(--border); border-radius:12px;
  padding:12px 16px; margin-bottom:16px; }
.filters .f { display:flex; flex-direction:column; gap:4px; }
.filters label { font-size:11.5px; color:var(--text-secondary); }
.filters .row { display:flex; gap:6px; align-items:center; }
.filters input[type=number] { width:88px; }
.filters .count { font-size:12.5px; color:var(--text-secondary); margin-left:auto; }

.legend { display:flex; flex-wrap:wrap; gap:6px 16px; font-size:12px;
  color:var(--text-secondary); margin:8px 0 2px; }
.legend .k { display:inline-flex; align-items:center; gap:6px; }
.legend .sw { width:12px; height:12px; border-radius:3px; display:inline-block; }
.legend .ln { width:16px; height:2px; border-radius:2px; display:inline-block; }
.rampbar { display:flex; align-items:center; gap:8px; font-size:12px;
  color:var(--text-secondary); }
.rampbar .bar { height:10px; width:160px; border-radius:5px; border:1px solid var(--border); }

svg { display:block; max-width:100%; overflow:visible; }
svg text { fill:var(--text-secondary); }
.tip { position:fixed; z-index:50; pointer-events:none; opacity:0; transition:opacity .08s;
  background:var(--surface-1); color:var(--text-primary); border:1px solid var(--border);
  box-shadow:0 6px 24px rgba(0,0,0,.16); border-radius:10px; padding:8px 10px;
  font-size:12.5px; max-width:320px; }
.tip .v { font-weight:600; font-variant-numeric:tabular-nums; }
.tip .n { color:var(--text-secondary); }
.tip .r { display:flex; align-items:center; gap:7px; }
.tip .k { width:14px; height:2px; border-radius:2px; flex:0 0 auto; }
.empty { border:1px dashed var(--axis); border-radius:10px; padding:22px;
  color:var(--text-secondary); font-size:13px; text-align:center; }
details.tv { margin-top:10px; }
details.tv summary { font-size:12.5px; color:var(--text-secondary); cursor:pointer; }
.viewer { position:relative; height:420px; border-radius:10px; overflow:hidden;
  border:1px solid var(--border); background:var(--surface-1); }
footer.foot { color:var(--muted); font-size:12px; border-top:1px solid var(--border);
  padding-top:14px; margin-top:8px; }
@media print { .viewer { display:none; } }
"""

# --- fornac's stylesheet, rewritten so it can only ever touch view 7 ----------------------
FORNA_SCOPE = "#ss-forna"  # the div fornaView() creates inside section 7's #ss-chart


def scoped_fornac_css(scope=FORNA_SCOPE):
    """Read vendor/fornac.css and prefix every selector with `scope`.

    The vendored sheet contains bare `svg { width:100%; min-height:100% }` and
    `text { pointer-events:none }` rules. Every other chart on this page is hand-built
    inline SVG, so an unscoped copy would visually destroy views 2, 4, 5 and 6. The file is
    flat (no at-rules, no nesting), so a split on braces is a sufficient parser -- anything
    unexpected is dropped and reported as a generator note rather than passed through.
    """
    p = HERE / "vendor" / "fornac.css"
    try:
        raw = p.read_text(encoding="utf-8")
    except Exception as exc:
        NOTES.append(f"vendor/fornac.css unreadable ({exc.__class__.__name__}: {exc}); "
                     f"fornac is not loaded and view 7 falls back to the arc diagram")
        return ""
    blocks, dropped = [], []
    for chunk in raw.split("}"):
        if "{" not in chunk:
            continue  # trailing whitespace after the last rule
        sel, body = chunk.split("{", 1)
        sel, body = " ".join(sel.split()), " ".join(body.split())
        if not sel or sel.startswith("@") or "{" in body:
            dropped.append(sel or "(empty selector)")
            continue
        parts = [s.strip() for s in sel.split(",") if s.strip()]
        if not parts:
            continue
        blocks.append(", ".join(f"{scope} {s}" for s in parts) + " { " + body + " }")
    if dropped:
        NOTES.append("vendor/fornac.css: could not scope " + ", ".join(dropped) +
                     " -- rule(s) dropped")
    # our own overrides, after the vendored rules so they win at equal specificity
    blocks += [
        f"{scope} {{ position:relative; overflow:hidden; border-radius:10px;",
        "  border:1px solid var(--border); background:var(--surface-1); }",
        # fornac.css sets `color` on the nucleotide labels, which does nothing to SVG text,
        # so the dashboard's own `svg text` rule would win. fornac tags each label with
        # label_type, which lets the letters and the position numbers be styled apart.
        # Every region colour in ssRegions() is dark enough for white letters (inkOn()).
        f"{scope} svg {{ overflow:hidden; }}",
        f"{scope} text {{ fill:var(--text-secondary); }}",
        f"{scope} text.node-label[label_type=nucleotide] {{ fill:#ffffff; }}",
        f"{scope} text.node-label[label_type=label] {{ fill:var(--text-secondary); }}",
        f"{scope} circle.node.label {{ fill:var(--surface-1); }}",
        f"{scope} circle.node {{ stroke:var(--border); }}",
        # fornac ships raw red / blue / green links; pull them into the page palette.
        # `line.link` is the backbone; base pairs carry both classes, and come later here,
        # so they win. Neither hue collides with the region colours on the nucleotides.
        f"{scope} line.link {{ stroke:var(--muted); }}",
        f"{scope} line.basepair {{ stroke:var(--text-primary); stroke-width:2; }}",
        f"{scope} line.pseudoknot {{ stroke:var(--series-4); }}",
        f"{scope} line.intermolecule {{ stroke:var(--series-3); }}",
    ]
    return "\n".join(blocks)


JS = r"""
'use strict';
const DATA = JSON.parse(document.getElementById('dash-data').textContent);
const NS = 'http://www.w3.org/2000/svg';

/* ---------- palette: both modes are selected, not flipped ---------- */
const PAL = {
  light: { s1:'#2a78d6', s2:'#eb6834', s3:'#1baf7a', s4:'#eda100',
           surface:'#fcfcfb', grid:'#e1e0d9', axis:'#c3c2b7', ink:'#0b0b0b',
           ink2:'#52514e', muted:'#898781', warn:'#fab219', crit:'#d03b3b',
           /* sequential blue 100->700 */
           seq:['#cde2fb','#b7d3f6','#9ec5f4','#86b6ef','#6da7ec','#5598e7','#3987e5',
                '#2a78d6','#256abf','#1c5cab','#184f95','#104281','#0d366b'],
           /* ordinal blue: nothing lighter than step 250 on the light surface */
           ord:['#86b6ef','#6da7ec','#3987e5','#2a78d6','#256abf','#1c5cab','#104281'] },
  dark:  { s1:'#3987e5', s2:'#d95926', s3:'#199e70', s4:'#c98500',
           surface:'#1a1a19', grid:'#2c2c2a', axis:'#383835', ink:'#ffffff',
           ink2:'#c3c2b7', muted:'#898781', warn:'#fab219', crit:'#d03b3b',
           seq:['#0d366b','#104281','#184f95','#1c5cab','#256abf','#2a78d6','#3987e5',
                '#5598e7','#6da7ec','#86b6ef','#9ec5f4','#b7d3f6','#cde2fb'],
           /* ordinal blue on dark: nothing darker than step 600 */
           ord:['#b7d3f6','#86b6ef','#5598e7','#3987e5','#2a78d6','#256abf','#184f95'] }
};
function mode() {
  const t = document.documentElement.dataset.theme;
  if (t === 'dark' || t === 'light') return t;
  return matchMedia('(prefers-color-scheme: dark)').matches ? 'dark' : 'light';
}
let P = PAL[mode()];
const SERIES = m => [m.s1, m.s2, m.s3, m.s4];
function seqColor(t) {          /* t in [0,1] -> sequential blue */
  const r = P.seq, i = Math.max(0, Math.min(r.length - 1, Math.round(t * (r.length - 1))));
  return r[i];
}
function ordColor(i, n) {       /* ordered classes -> ordinal blue */
  const r = P.ord, k = n <= 1 ? 0 : Math.round(i / (n - 1) * (r.length - 1));
  return r[k];
}
function inkOn(hex) {           /* label inside a fill: pick by luminance */
  const h = hex.replace('#',''), n = parseInt(h, 16);
  const l = (0.2126*((n>>16)&255) + 0.7152*((n>>8)&255) + 0.0722*(n&255)) / 255;
  return l > 0.55 ? '#0b0b0b' : '#ffffff';
}

/* ---------- tiny SVG + scale helpers ---------- */
function el(tag, attrs, parent) {
  const n = document.createElementNS(NS, tag);
  for (const k in (attrs || {})) if (attrs[k] !== null && attrs[k] !== undefined)
    n.setAttribute(k, attrs[k]);
  if (parent) parent.appendChild(n);
  return n;
}
function txt(parent, x, y, s, attrs) {
  const t = el('text', Object.assign({x:x, y:y, 'font-size':11}, attrs || {}), parent);
  t.appendChild(document.createTextNode(s));
  return t;
}
const lin = (d0,d1,r0,r1) => v => r0 + (v - d0) / (d1 - d0) * (r1 - r0);
const log = (d0,d1,r0,r1) => {
  const a = Math.log10(d0), b = Math.log10(d1);
  return v => r0 + (Math.log10(Math.max(v, 1e-300)) - a) / (b - a) * (r1 - r0);
};
const fmt = (v, d) => v.toLocaleString(undefined, {minimumFractionDigits:d===undefined?0:d,
                                                  maximumFractionDigits:d===undefined?1:d});
const sup = n => ({'-':'⁻','0':'⁰','1':'¹','2':'²','3':'³','4':'⁴',
                   '5':'⁵','6':'⁶','7':'⁷','8':'⁸','9':'⁹'})[n] || n;
const pow10 = e => '10' + String(e).split('').map(sup).join('');

/* ---------- one tooltip for the whole page ---------- */
const TIP = document.createElement('div');
TIP.className = 'tip'; document.body.appendChild(TIP);
function showTip(ev, rows) {
  TIP.textContent = '';
  rows.forEach(r => {
    const d = document.createElement('div');
    d.className = 'r';
    if (r.color) { const k = document.createElement('span');
      k.className = 'k'; k.style.background = r.color; d.appendChild(k); }
    if (r.value !== undefined) { const v = document.createElement('span');
      v.className = 'v'; v.textContent = r.value; d.appendChild(v); }
    if (r.name) { const n = document.createElement('span');
      n.className = 'n'; n.textContent = r.name; d.appendChild(n); }
    TIP.appendChild(d);
  });
  TIP.style.opacity = 1;
  const pad = 14, w = TIP.offsetWidth, h = TIP.offsetHeight;
  let x = ev.clientX + pad, y = ev.clientY + pad;
  if (x + w > innerWidth - 8) x = ev.clientX - w - pad;
  if (y + h > innerHeight - 8) y = ev.clientY - h - pad;
  TIP.style.left = x + 'px'; TIP.style.top = y + 'px';
}
const hideTip = () => { TIP.style.opacity = 0; };
addEventListener('scroll', hideTip, {passive:true});

function grid(g, x1, x2, y, first) {
  el('line', {x1:x1, x2:x2, y1:y, y2:y, stroke:first ? P.axis : P.grid,
              'stroke-width':1, 'shape-rendering':'crispEdges'}, g);
}
function card(id) { return document.getElementById(id); }
function empty(node, msg) {
  node.textContent = '';
  const d = document.createElement('div');
  d.className = 'empty'; d.textContent = msg; node.appendChild(d);
}
function guard(id, fn) {
  const n = card(id); if (!n) return;
  try { fn(n); } catch (e) { empty(n, 'could not render this view: ' + e.message); }
}

/* =========================================================================
   VIEW 1 - parent aptamer table
   ========================================================================= */
function renderParents() {
  guard('parents-body', node => {
    const rows = DATA.parents;
    if (!rows || !rows.length) return empty(node, 'parents.json not yet generated - run python aptamer/build_parents.py');
    node.textContent = '';
    const wrap = document.createElement('div'); wrap.className = 'scroll';
    const t = document.createElement('table');
    t.innerHTML = '';
    const head = ['Aptamer', 'Variant', 'Length (nt)', 'K_D (nM)', 'k_off (s⁻¹)',
                  'Residence time 1/k_off (s)', 'MFE 22 °C (kcal/mol)',
                  'MFE 37 °C (kcal/mol)', 'Reconstruction'];
    const thead = document.createElement('thead'), htr = document.createElement('tr');
    head.forEach((h, i) => { const th = document.createElement('th');
      th.textContent = h; if (i < 2 || i === 8) th.className = 'l'; htr.appendChild(th); });
    thead.appendChild(htr); t.appendChild(thead);
    const tb = document.createElement('tbody');
    rows.forEach(r => {
      const tr = document.createElement('tr');
      const cells = [
        [r.name, 'l'], [r.variant, 'l'], [String(r.length), ''],
        [fmt(r.kd_nM, 1), ''], [(r.koff_s || 0).toExponential(2), ''],
        [r.residence_s ? fmt(r.residence_s, 0) : '—', ''],
        [(r.mfe_22C).toFixed(2), ''], [(r.mfe_37C).toFixed(2), ''], [null, 'l']
      ];
      cells.forEach(([v, cls]) => {
        const td = document.createElement('td'); if (cls) td.className = cls;
        if (v !== null) td.textContent = v;
        else {
          const b = document.createElement('span');
          const bad = (r.issues || []).length;
          b.className = 'badge ' + (bad ? 'warn' : 'ok');
          const ic = document.createElement('span'); ic.className = 'ic';
          ic.textContent = bad ? '⚠' : '✓'; b.appendChild(ic);
          b.appendChild(document.createTextNode(bad ? bad + ' uncertain nt' : 'clean'));
          td.appendChild(b);
        }
        tr.appendChild(td);
      });
      tb.appendChild(tr);
      const sr = document.createElement('tr');
      const td = document.createElement('td');
      td.colSpan = head.length; td.className = 'l';
      const s = document.createElement('div'); s.className = 'seq';
      s.textContent = r.sequence; td.appendChild(s);
      const d = document.createElement('div');
      d.className = 'seq muted'; d.style.marginTop = '2px';
      d.textContent = r.structure_37C + '   (37 °C)';
      td.appendChild(d);
      sr.appendChild(td); tb.appendChild(sr);
    });
    t.appendChild(tb); wrap.appendChild(t); node.appendChild(wrap);

    const bad = rows.filter(r => (r.issues || []).length);
    if (bad.length) {
      const box = document.createElement('div'); box.className = 'issues';
      const b = document.createElement('b');
      b.textContent = '⚠ Do not order these blind — reconstruction warnings';
      box.appendChild(b);
      const ul = document.createElement('ul');
      bad.forEach(r => (r.issues || []).forEach(i => {
        const li = document.createElement('li');
        li.textContent = r.name + ': ' + i;
        ul.appendChild(li);
      }));
      box.appendChild(ul);
      node.appendChild(box);
    }
  });
}

/* =========================================================================
   VIEW 2 - occupancy curve
   ========================================================================= */
function renderOccupancy() {
  guard('occ-chart', node => {
    const S = DATA.occupancy;
    if (!S || !S.length) return empty(node, 'parents.json not yet generated — no occupancy curve');
    node.textContent = '';
    const W = Math.max(520, Math.min(node.clientWidth || 760, 980)), H = 380;
    const m = {t:14, r:120, b:46, l:52};
    const svg = el('svg', {width:'100%', viewBox:'0 0 ' + W + ' ' + H,
                           role:'img', 'aria-label':'Fraction of aptamer bound vs IL-6 concentration'}, node);
    const x = log(10, 1e7, m.l, W - m.r), y = lin(0, 1, H - m.b, m.t);

    /* wound-fluid bands: context washes, drawn under the data */
    const bands = Object.entries(DATA.wound);
    bands.forEach(([label, pg], i) => {
      const px = x(pg);
      el('line', {x1:px, x2:px, y1:m.t, y2:H - m.b, stroke:P.axis, 'stroke-width':1,
                  'shape-rendering':'crispEdges'}, svg);
      const lab = txt(svg, px + 4, m.t + 12 + (i % 2) * 13, label.replace(' (median)', ''),
                      {'font-size':10.5, fill:P.muted});
      lab.setAttribute('text-anchor', 'start');
    });
    el('rect', {x:x(bands[0][1]), y:m.t, width:Math.max(1, x(bands[2][1]) - x(bands[0][1])),
                height:H - m.b - m.t, fill:P.s1, opacity:0.07}, svg);
    /* skin-model reading */
    const sx = x(DATA.skinModel);
    el('line', {x1:sx, x2:sx, y1:m.t, y2:H - m.b, stroke:P.warn, 'stroke-width':2}, svg);
    txt(svg, sx + 5, H - m.b - 6, 'skin model 5 ng/mL', {'font-size':10.5, fill:P.ink2});

    /* axes */
    for (let f = 0; f <= 1.0001; f += 0.25) {
      grid(svg, m.l, W - m.r, y(f), f === 0);
      txt(svg, m.l - 8, y(f) + 4, Math.round(f * 100) + '%',
          {'text-anchor':'end', 'font-size':10.5, fill:P.muted});
    }
    for (let e = 1; e <= 7; e++) {
      txt(svg, x(Math.pow(10, e)), H - m.b + 16, pow10(e),
          {'text-anchor':'middle', 'font-size':10.5, fill:P.muted});
    }
    txt(svg, (m.l + W - m.r) / 2, H - 8, 'IL-6 (pg/mL, log scale)',
        {'text-anchor':'middle', 'font-size':11.5, fill:P.ink2});
    txt(svg, m.l - 8, m.t - 2, 'fraction bound', {'text-anchor':'end', 'font-size':11.5, fill:P.ink2});

    /* lines + a direct label at each curve's 50% crossing (its K_D) */
    const cols = SERIES(P);
    S.forEach((s, i) => {
      const d = s.pts.map((p, k) => (k ? 'L' : 'M') + x(p[0]).toFixed(1) + ' ' + y(p[1]).toFixed(1)).join(' ');
      el('path', {d:d, fill:'none', stroke:cols[i % 4], 'stroke-width':2,
                  'stroke-linejoin':'round', 'stroke-linecap':'round'}, svg);
      const hx = x(s.half_pg), hy = y(0.5);
      el('circle', {cx:hx, cy:hy, r:4.5, fill:cols[i % 4], stroke:P.surface, 'stroke-width':2}, svg);
      txt(svg, hx + 9, hy + 4 + (i - 1.5) * 13, s.name + ' · K_D ' + fmt(s.kd_nM, 1) + ' nM',
          {'font-size':11, fill:P.ink2});
    });

    /* crosshair: the reader aims at a concentration, never at a 2px line */
    const cross = el('line', {y1:m.t, y2:H - m.b, stroke:P.axis, 'stroke-width':1,
                              opacity:0, 'pointer-events':'none'}, svg);
    const dots = S.map((s, i) => el('circle', {r:4, fill:cols[i % 4], stroke:P.surface,
                                               'stroke-width':2, opacity:0, 'pointer-events':'none'}, svg));
    const hit = el('rect', {x:m.l, y:m.t, width:W - m.r - m.l, height:H - m.b - m.t,
                            fill:'transparent'}, svg);
    const invX = px => Math.pow(10, 1 + (px - m.l) / (W - m.r - m.l) * 6);
    function move(ev) {
      const r = svg.getBoundingClientRect(), sc = W / r.width;
      const px = (ev.clientX - r.left) * sc;
      if (px < m.l || px > W - m.r) return leave();
      const pg = invX(px);
      cross.setAttribute('x1', px); cross.setAttribute('x2', px); cross.setAttribute('opacity', 1);
      const rows = [{value: fmt(pg, 0) + ' pg/mL', name:'IL-6'}];
      S.forEach((s, i) => {
        const k = Math.max(0, Math.min(s.pts.length - 1,
          Math.round((Math.log10(pg) - 1) / 6 * (s.pts.length - 1))));
        const p = s.pts[k];
        dots[i].setAttribute('cx', x(p[0])); dots[i].setAttribute('cy', y(p[1]));
        dots[i].setAttribute('opacity', 1);
        rows.push({color:cols[i % 4], value:(p[1] * 100).toFixed(2) + '%', name:s.name + ' bound'});
      });
      showTip(ev, rows);
    }
    function leave() {
      cross.setAttribute('opacity', 0); dots.forEach(d => d.setAttribute('opacity', 0)); hideTip();
    }
    hit.addEventListener('pointermove', move);
    hit.addEventListener('pointerleave', leave);

    /* legend */
    const lg = document.createElement('div'); lg.className = 'legend';
    S.forEach((s, i) => {
      const k = document.createElement('span'); k.className = 'k';
      const ln = document.createElement('span'); ln.className = 'ln';
      ln.style.background = cols[i % 4]; k.appendChild(ln);
      k.appendChild(document.createTextNode(s.name));
      lg.appendChild(k);
    });
    const wb = document.createElement('span'); wb.className = 'k';
    const ws = document.createElement('span'); ws.className = 'sw';
    ws.style.background = P.s1; ws.style.opacity = 0.25; wb.appendChild(ws);
    wb.appendChild(document.createTextNode('wound-fluid range (post-surgical → infected)'));
    lg.appendChild(wb);
    node.appendChild(lg);
  });

  /* table-view twin */
  guard('occ-table', node => {
    const rows = DATA.occupancyTable, ps = DATA.parents || [];
    if (!rows || !rows.length || !ps.length) return empty(node, 'not yet generated');
    node.textContent = '';
    const t = document.createElement('table');
    const thead = document.createElement('thead'), tr = document.createElement('tr');
    ['IL-6 level', 'pg/mL', 'pM'].concat(ps.map(p => p.name + ' bound (%)')).forEach((h, i) => {
      const th = document.createElement('th'); th.textContent = h;
      if (!i) th.className = 'l'; tr.appendChild(th);
    });
    thead.appendChild(tr); t.appendChild(thead);
    const tb = document.createElement('tbody');
    rows.forEach(r => {
      const trr = document.createElement('tr');
      const c = [[r.label, 'l'], [fmt(r.pg_ml, 0), ''], [fmt(r.pM, 1), '']]
        .concat(ps.map(p => [(100 * (r.frac[p.name] || 0)).toFixed(2), '']));
      c.forEach(([v, cls]) => { const td = document.createElement('td');
        if (cls) td.className = cls; td.textContent = v; trr.appendChild(td); });
      tb.appendChild(trr);
    });
    t.appendChild(tb); node.appendChild(t);
  });
}

/* =========================================================================
   VIEWS 3-6 - everything scoped by the one filter row
   ========================================================================= */
const SW = DATA.switches || [];
const COLS = [
  {k:'rank', h:'Rank', d:0, l:false},
  {k:'ds', h:'DS (5′→3′)', text:true},
  {k:'window', h:'Window (0-idx)', text:true},
  {k:'ds_len', h:'DS length (nt)', d:0},
  {k:'linker_len', h:'Linker (nt)', d:0},
  {k:'gc', h:'GC (fraction)', d:2},
  {k:'dg_switch', h:'ΔG_switch (kcal/mol)', d:2},
  {k:'closed_frac', h:'Closed (fraction)', d:3},
  {k:'kd_app_nM', h:'K_D,app (nM)', d:1},
  {k:'engagement', h:'Engagement (fraction)', d:2},
  {k:'rand_covered', h:'Randomised nt covered', d:0}
];
const FILTERS = [
  {k:'dg_switch', label:'ΔG_switch (kcal/mol)', step:0.05},
  {k:'kd_app_nM', label:'K_D,app (nM)', step:0.1},
  {k:'engagement', label:'Engagement (fraction)', step:0.01}
];
const RANGE = {};
FILTERS.forEach(f => {
  const v = SW.map(r => r[f.k]);
  RANGE[f.k] = SW.length ? [Math.min.apply(null, v), Math.max.apply(null, v)] : [0, 1];
});
let sortKey = 'rank', sortDir = 1, filtered = SW.slice();

function applyFilters() {
  filtered = SW.filter(r => FILTERS.every(f => {
    const a = document.getElementById('lo-' + f.k), b = document.getElementById('hi-' + f.k);
    if (!a || !b) return true;                     /* filter row absent: show everything */
    const lo = parseFloat(a.value), hi = parseFloat(b.value);
    return (isNaN(lo) || r[f.k] >= lo) && (isNaN(hi) || r[f.k] <= hi);
  }));
  const c = document.getElementById('sw-count');
  if (c) c.textContent = filtered.length + ' of ' + SW.length + ' constructs shown';
  renderSwitchTable(); renderScatter(); renderHeatmap(); renderCoverage();
  buildStructurePicker();   /* view 7 sits below the filter row, so it follows the slice */
}
function buildFilterRow() {
  const node = card('sw-filters'); if (!node) return;
  if (!SW.length) { empty(node, 'switches.csv not yet generated — run python aptamer/switch_library.py'); return; }
  node.textContent = '';
  FILTERS.forEach(f => {
    const d = document.createElement('div'); d.className = 'f';
    const lab = document.createElement('label'); lab.textContent = f.label; d.appendChild(lab);
    const row = document.createElement('div'); row.className = 'row';
    ['lo', 'hi'].forEach((side, i) => {
      const inp = document.createElement('input');
      inp.type = 'number'; inp.step = f.step; inp.id = side + '-' + f.k;
      inp.value = RANGE[f.k][i]; inp.setAttribute('aria-label', f.label + ' ' + (i ? 'max' : 'min'));
      inp.addEventListener('change', applyFilters);
      row.appendChild(inp);
      if (!i) row.appendChild(document.createTextNode('–'));
    });
    d.appendChild(row); node.appendChild(d);
  });
  const b = document.createElement('button');
  b.textContent = 'Reset'; b.addEventListener('click', () => {
    FILTERS.forEach(f => { document.getElementById('lo-' + f.k).value = RANGE[f.k][0];
                           document.getElementById('hi-' + f.k).value = RANGE[f.k][1]; });
    sortKey = 'rank'; sortDir = 1; applyFilters();
  });
  node.appendChild(b);
  const c = document.createElement('span'); c.className = 'count'; c.id = 'sw-count';
  node.appendChild(c);
}

/* --- VIEW 3: switch library table --- */
function renderSwitchTable() {
  guard('sw-table', node => {
    if (!SW.length) return empty(node, 'switches.csv not yet generated');
    const rows = filtered.slice().sort((a, b) => {
      const x = a[sortKey], y = b[sortKey];
      if (typeof x === 'string') return sortDir * String(x).localeCompare(String(y));
      return sortDir * (x - y);
    });
    node.textContent = '';
    const wrap = document.createElement('div'); wrap.className = 'scroll';
    const t = document.createElement('table');
    const thead = document.createElement('thead'), tr = document.createElement('tr');
    COLS.forEach(c => {
      const th = document.createElement('th');
      th.className = 'sortable' + (c.text || c.k === 'rank' ? ' l' : '');
      th.textContent = c.h;
      const a = document.createElement('span'); a.className = 'arrow';
      a.textContent = sortKey === c.k ? (sortDir > 0 ? ' ▲' : ' ▼') : '';
      th.appendChild(a);
      th.title = 'sort by ' + c.h;
      th.addEventListener('click', () => {
        if (sortKey === c.k) sortDir = -sortDir;
        else { sortKey = c.k; sortDir = (c.k === 'rank' || c.text) ? 1 : 1; }
        renderSwitchTable();
      });
      tr.appendChild(th);
    });
    thead.appendChild(tr); t.appendChild(thead);
    const tb = document.createElement('tbody');
    rows.forEach(r => {
      const trr = document.createElement('tr');
      COLS.forEach(c => {
        const td = document.createElement('td');
        if (c.text || c.k === 'rank') td.className = 'l';
        if (c.text) { td.textContent = r[c.k]; td.classList.add('mono'); }
        else td.textContent = fmt(r[c.k], c.d);
        trr.appendChild(td);
      });
      trr.title = r.construct;
      tb.appendChild(trr);
    });
    t.appendChild(tb); wrap.appendChild(t); node.appendChild(wrap);
    if (!rows.length) empty(node, 'no constructs match the filters');
  });
}

/* --- VIEW 4: gain-vs-affinity scatter --- */
function renderScatter() {
  guard('sw-scatter', node => {
    if (!SW.length) return empty(node, 'switches.csv not yet generated');
    node.textContent = '';
    const W = Math.max(480, Math.min(node.clientWidth || 720, 900)), H = 400;
    const m = {t:16, r:20, b:48, l:60};
    const svg = el('svg', {width:'100%', viewBox:'0 0 ' + W + ' ' + H, role:'img',
      'aria-label':'Apparent K_D versus switch stabilisation energy'}, node);
    const dg = DATA.dgWindow;
    const x = lin(dg[0] - 0.35, dg[1] + 0.35, m.l, W - m.r);
    const ky = SW.map(r => r.kd_app_nM);
    const y0 = Math.pow(10, Math.floor(Math.log10(Math.min.apply(null, ky))));
    const y1 = Math.pow(10, Math.ceil(Math.log10(Math.max.apply(null, ky))));
    const y = log(y0, y1, H - m.b, m.t);

    /* design window as a shaded region */
    el('rect', {x:x(dg[0]), y:m.t, width:x(dg[1]) - x(dg[0]), height:H - m.b - m.t,
                fill:P.s1, opacity:0.08}, svg);
    txt(svg, (x(dg[0]) + x(dg[1])) / 2, m.t + 13,
        (DATA.dgWindowCaption || 'design window') + ' · ' +
        dg[0].toFixed(2) + ' to ' + dg[1].toFixed(2) + ' kcal/mol',
        {'text-anchor':'middle', 'font-size':10.5, fill:P.ink2});

    for (let e = Math.log10(y0); e <= Math.log10(y1) + 1e-9; e++) {
      const v = Math.pow(10, e);
      grid(svg, m.l, W - m.r, y(v), false);
      txt(svg, m.l - 8, y(v) + 4, fmt(v, 0), {'text-anchor':'end', 'font-size':10.5, fill:P.muted});
      for (let k = 2; k < 10; k++) {
        const vv = v * k; if (vv > y1) break;
        grid(svg, m.l, W - m.r, y(vv), false);
      }
    }
    for (let v = Math.ceil((dg[0] - 0.35) * 4) / 4; v <= dg[1] + 0.35; v += 0.25) {
      txt(svg, x(v), H - m.b + 16, v.toFixed(2), {'text-anchor':'middle', 'font-size':10.5, fill:P.muted});
    }
    el('line', {x1:m.l, x2:W - m.r, y1:H - m.b, y2:H - m.b, stroke:P.axis, 'shape-rendering':'crispEdges'}, svg);
    txt(svg, (m.l + W - m.r) / 2, H - 10, 'ΔG_switch (kcal/mol) — more negative = more strongly closed',
        {'text-anchor':'middle', 'font-size':11.5, fill:P.ink2});
    txt(svg, m.l - 8, m.t - 2, 'K_D,app (nM, log)', {'text-anchor':'end', 'font-size':11.5, fill:P.ink2});

    const rc = SW.map(r => r.rand_covered);
    const rMax = Math.max.apply(null, rc), rMin = Math.min.apply(null, rc);
    const eMin = Math.min.apply(null, SW.map(r => r.engagement));
    const eMax = Math.max.apply(null, SW.map(r => r.engagement));
    const rad = e => 4 + (eMax > eMin ? (e - eMin) / (eMax - eMin) : 0.5) * 4;
    const pts = [];
    filtered.forEach(r => {
      const cx = x(r.dg_switch), cy = y(r.kd_app_nM);
      const col = ordColor(r.rand_covered - rMin, rMax - rMin + 1);
      el('circle', {cx:cx, cy:cy, r:rad(r.engagement), fill:col, 'fill-opacity':0.85,
                    stroke:P.surface, 'stroke-width':2}, svg);
      pts.push({x:cx, y:cy, r:r, col:col});
    });

    /* nearest-point hover: the pointer only has to be closest, not dead centre */
    const halo = el('circle', {r:0, fill:'none', stroke:P.ink, 'stroke-width':1.5, opacity:0}, svg);
    const hit = el('rect', {x:m.l, y:m.t, width:W - m.r - m.l, height:H - m.b - m.t, fill:'transparent'}, svg);
    hit.addEventListener('pointermove', ev => {
      const bb = svg.getBoundingClientRect(), sc = W / bb.width;
      const px = (ev.clientX - bb.left) * sc, py = (ev.clientY - bb.top) * sc;
      let best = null, bd = 1e9;
      pts.forEach(p => { const d = (p.x - px) ** 2 + (p.y - py) ** 2;
                         if (d < bd) { bd = d; best = p; } });
      if (!best || bd > 900) { halo.setAttribute('opacity', 0); return hideTip(); }
      halo.setAttribute('cx', best.x); halo.setAttribute('cy', best.y);
      halo.setAttribute('r', rad(best.r.engagement) + 4); halo.setAttribute('opacity', 0.9);
      const r = best.r;
      showTip(ev, [
        {value:'#' + r.rank + ' ' + r.ds, color:best.col},
        {value:fmt(r.kd_app_nM, 1) + ' nM', name:'K_D,app'},
        {value:r.dg_switch.toFixed(2) + ' kcal/mol', name:'ΔG_switch'},
        {value:r.engagement.toFixed(2), name:'engagement'},
        {value:String(r.rand_covered), name:'randomised nt covered'},
        {value:r.window, name:'window · linker ' + r.linker_len + ' nt'}
      ]);
    });
    hit.addEventListener('pointerleave', () => { halo.setAttribute('opacity', 0); hideTip(); });

    /* legends: colour ramp + size key */
    const lg = document.createElement('div'); lg.className = 'legend';
    const ramp = document.createElement('span'); ramp.className = 'rampbar';
    ramp.appendChild(document.createTextNode('randomised nt covered ' + rMin));
    const bar = document.createElement('span'); bar.className = 'bar';
    const stops = [];
    for (let i = rMin; i <= rMax; i++) stops.push(ordColor(i - rMin, rMax - rMin + 1));
    bar.style.background = 'linear-gradient(90deg,' + stops.join(',') + ')';
    ramp.appendChild(bar);
    ramp.appendChild(document.createTextNode(String(rMax)));
    lg.appendChild(ramp);
    const sz = document.createElement('span'); sz.className = 'k';
    sz.textContent = 'point size = engagement (' + eMin.toFixed(2) + '–' + eMax.toFixed(2) + ')';
    lg.appendChild(sz);
    node.appendChild(lg);
  });
}

/* --- VIEW 5: linker x DS-length heatmap --- */
function renderHeatmap() {
  guard('sw-heatmap', node => {
    if (!SW.length) return empty(node, 'switches.csv not yet generated');
    node.textContent = '';
    const dsL = [...new Set(SW.map(r => r.ds_len))].sort((a, b) => a - b);
    const lnL = [...new Set(SW.map(r => r.linker_len))].sort((a, b) => a - b);
    const best = {}, cnt = {};
    filtered.forEach(r => {
      const k = r.ds_len + '|' + r.linker_len;
      if (best[k] === undefined || r.kd_app_nM < best[k]) best[k] = r.kd_app_nM;
      cnt[k] = (cnt[k] || 0) + 1;
    });
    /* fixed global domain: filtering must not repaint the survivors */
    const all = SW.map(r => r.kd_app_nM);
    const lo = Math.min.apply(null, all), hi = Math.max.apply(null, all);
    const cw = 52, ch = 34, gap = 2, m = {t:24, r:16, b:34, l:96};
    const W = m.l + lnL.length * cw + m.r, H = m.t + dsL.length * ch + m.b;
    const svg = el('svg', {width:'100%', viewBox:'0 0 ' + W + ' ' + H, role:'img',
      'aria-label':'Best apparent K_D for each linker length and DS length'}, node);
    lnL.forEach((ln, j) => txt(svg, m.l + j * cw + cw / 2, m.t - 8, String(ln),
      {'text-anchor':'middle', 'font-size':10.5, fill:P.muted}));
    txt(svg, m.l + lnL.length * cw / 2, H - 8, 'linker length (nt)',
        {'text-anchor':'middle', 'font-size':11.5, fill:P.ink2});
    txt(svg, m.l - 10, m.t - 8, 'DS length (nt)', {'text-anchor':'end', 'font-size':11.5, fill:P.ink2});
    dsL.forEach((ds, i) => {
      txt(svg, m.l - 10, m.t + i * ch + ch / 2 + 4, String(ds),
          {'text-anchor':'end', 'font-size':10.5, fill:P.muted});
      lnL.forEach((ln, j) => {
        const k = ds + '|' + ln, v = best[k];
        const gx = m.l + j * cw + gap / 2, gy = m.t + i * ch + gap / 2;
        if (v === undefined) {
          el('rect', {x:gx, y:gy, width:cw - gap, height:ch - gap, rx:4, fill:'none',
                      stroke:P.grid, 'stroke-width':1}, svg);
          return;
        }
        const t = hi > lo ? (v - lo) / (hi - lo) : 0.5;
        const col = seqColor(t);
        const cell = el('rect', {x:gx, y:gy, width:cw - gap, height:ch - gap, rx:4, fill:col}, svg);
        txt(svg, gx + (cw - gap) / 2, gy + (ch - gap) / 2 + 4, fmt(v, 0),
            {'text-anchor':'middle', 'font-size':10.5, fill:inkOn(col)});
        cell.addEventListener('pointermove', ev => {
          cell.setAttribute('stroke', P.ink); cell.setAttribute('stroke-width', 1.5);
          showTip(ev, [{value:fmt(v, 1) + ' nM', name:'best K_D,app', color:col},
                       {value:'DS ' + ds + ' nt · linker ' + ln + ' nt'},
                       {value:cnt[k] + ' construct(s)', name:'in this cell'}]);
        });
        cell.addEventListener('pointerleave', () => { cell.removeAttribute('stroke'); hideTip(); });
      });
    });
    const lg = document.createElement('div'); lg.className = 'legend';
    const ramp = document.createElement('span'); ramp.className = 'rampbar';
    ramp.appendChild(document.createTextNode('best K_D,app  ' + fmt(lo, 0) + ' nM'));
    const bar = document.createElement('span'); bar.className = 'bar';
    bar.style.background = 'linear-gradient(90deg,' + P.seq.join(',') + ')';
    ramp.appendChild(bar);
    ramp.appendChild(document.createTextNode(fmt(hi, 0) + ' nM'));
    lg.appendChild(ramp);
    const nk = document.createElement('span'); nk.className = 'k';
    const sw = document.createElement('span'); sw.className = 'sw';
    sw.style.background = 'transparent'; sw.style.border = '1px solid ' + P.grid;
    nk.appendChild(sw); nk.appendChild(document.createTextNode('no passing construct'));
    lg.appendChild(nk);
    node.appendChild(lg);
  });
}

/* --- VIEW 6: window coverage track --- */
function renderCoverage() {
  guard('sw-coverage', node => {
    const seq = DATA.parentSeq;
    if (!SW.length || !seq) return empty(node, 'needs both parents.json and switches.csv');
    node.textContent = '';
    const n = seq.length, rand = new Set(DATA.randomised);
    const cov = new Array(n).fill(0);
    filtered.forEach(r => { for (let i = r.w0; i <= r.w1 && i < n; i++) cov[i]++; });
    const maxc = Math.max(1, Math.max.apply(null, cov));
    const cellW = 15, m = {t:18, r:16, b:64, l:44};
    const barH = 150;
    const W = m.l + n * cellW + m.r, H = m.t + barH + m.b;
    const svg = el('svg', {width:'100%', viewBox:'0 0 ' + W + ' ' + H, role:'img',
      'aria-label':'Constructs covering each position of the parent aptamer'}, svg0(node));
    const y = lin(0, maxc, m.t + barH, m.t);
    [0, 0.5, 1].forEach(f => {
      const v = maxc * f;
      grid(svg, m.l, W - m.r, y(v), f === 0);
      txt(svg, m.l - 8, y(v) + 4, fmt(v, 0), {'text-anchor':'end', 'font-size':10.5, fill:P.muted});
    });
    txt(svg, m.l - 8, m.t - 4, 'constructs', {'text-anchor':'end', 'font-size':11, fill:P.ink2});

    /* randomised library blocks: the wash marks them behind everything */
    let i = 0;
    while (i < n) {
      if (rand.has(i)) {
        let j = i; while (j + 1 < n && rand.has(j + 1)) j++;
        el('rect', {x:m.l + i * cellW, y:m.t, width:(j - i + 1) * cellW, height:barH + 34,
                    fill:P.s2, opacity:0.13, rx:3}, svg);
        i = j + 1;
      } else i++;
    }
    /* coverage columns */
    for (let k = 0; k < n; k++) {
      const bw = cellW - 3;
      const gx = m.l + k * cellW + 1.5;
      const h = y(0) - y(cov[k]);
      if (h > 0.5) el('rect', {x:gx, y:y(cov[k]), width:bw, height:h, rx:Math.min(4, bw / 2),
                               fill:P.s1}, svg);
      const t = txt(svg, gx + bw / 2, m.t + barH + 18, seq[k],
                    {'text-anchor':'middle', 'font-size':10.5,
                     fill:rand.has(k) ? P.ink : P.muted,
                     'font-family':'ui-monospace, Consolas, monospace'});
      if (rand.has(k)) t.setAttribute('font-weight', '600');
      if (k % 10 === 0) txt(svg, gx + bw / 2, m.t + barH + 32, String(k),
                            {'text-anchor':'middle', 'font-size':9.5, fill:P.muted});
      const hit = el('rect', {x:m.l + k * cellW, y:m.t, width:cellW, height:barH + 24,
                              fill:'transparent'}, svg);
      hit.addEventListener('pointermove', ev => showTip(ev, [
        {value:cov[k] + ' construct(s)', name:'target position ' + k, color:P.s1},
        {value:seq[k] + ' · ' + (rand.has(k) ? 'randomised block' : 'fixed region')}
      ]));
      hit.addEventListener('pointerleave', hideTip);
    }
    txt(svg, m.l + n * cellW / 2, H - 22, DATA.switchParent + ' — position in the ' + n + '-nt parent (0-indexed)',
        {'text-anchor':'middle', 'font-size':11.5, fill:P.ink2});

    const lg = document.createElement('div'); lg.className = 'legend';
    [[P.s1, 'constructs whose DS window covers this position', 'sw'],
     ['rgba(235,104,52,0.35)', 'randomised library block (bold letters)', 'sw']]
      .forEach(([c, label]) => {
        const k = document.createElement('span'); k.className = 'k';
        const s = document.createElement('span'); s.className = 'sw'; s.style.background = c;
        k.appendChild(s); k.appendChild(document.createTextNode(label)); lg.appendChild(k);
      });
    node.appendChild(lg);
  });
}
/* the coverage track is wider than the card; give it its own horizontal scroller */
function svg0(node) {
  const d = document.createElement('div');
  d.style.overflowX = 'auto'; node.appendChild(d); return d;
}

/* =========================================================================
   VIEW 7 - secondary structure: fornac 2D layout, arc diagram as the fallback
   ========================================================================= */
function pairs(db) {
  const st = [], out = [];
  for (let i = 0; i < db.length; i++) {
    const c = db[i];
    if (c === '(') st.push(i);
    else if (c === ')') { const j = st.pop(); if (j !== undefined) out.push([j, i]); }
  }
  return out;
}
function structureOptions() {
  const opts = [];
  (DATA.parents || []).forEach(p => {
    opts.push({label:'parent ' + p.name + ' @ 22 °C (MFE ' + p.mfe_22C.toFixed(2) + ')',
               seq:p.sequence, db:p.structure_22C, kind:'parent'});
    opts.push({label:'parent ' + p.name + ' @ 37 °C (MFE ' + p.mfe_37C.toFixed(2) + ')',
               seq:p.sequence, db:p.structure_37C, kind:'parent'});
  });
  /* switch constructs follow the filter row, like every other view below it */
  filtered.slice(0, 25).forEach(r => opts.push({
    label:'switch #' + r.rank + ' · DS ' + r.ds + ' · window ' + r.window +
          ' · linker ' + r.linker_len,
    seq:r.construct, db:r.structure, kind:'switch', row:r
  }));
  return opts;
}
/* ---------- fornac (primary renderer) ----------
   fornac 1.1.8 is a webpack UMD bundle, so in a browser it lands on window.fornac; the
   forna site's own build also exposes window.FornaContainer. Accept either, and treat
   "neither is there" as the signal to fall back. d3 3.5.13 is bundled inside fornac.js
   and jQuery is never referenced, so no other CDN script is needed. */
function fornaCtor() {
  return window.FornaContainer ||
         (window.fornac && window.fornac.FornaContainer) || null;
}
/* fornac carries its own stylesheet and style-loader injects it into <head> unscoped. Its
   bare `svg { width:100%; min-height:100% }` and `text {}` rules would resize every other
   chart on this page, so delete that copy -- a copy scoped to #ss-forna is already inlined
   in the page's own <style data-dash>. */
function dropFornacGlobalCss() {
  const styles = document.querySelectorAll('head style:not([data-dash])');
  for (let i = 0; i < styles.length; i++) {
    const t = styles[i].textContent || '';
    if (t.indexOf('circle.outline_node') >= 0 && /(^|\})\s*svg\s*\{/.test(t))
      styles[i].remove();
  }
}
/* [start, end) regions of the construct, in the page palette */
function ssRegions(opt) {
  const n = opt.seq.length;
  if (opt.kind !== 'switch') return [[0, n, 'aptamer', P.s1]];
  const nDs = opt.row.ds_len, nLink = opt.row.linker_len, nAp = n - nDs - nLink;
  return [[0, nAp, 'aptamer', P.s1], [nAp, nAp + nLink, 'poly-T linker', P.s3],
          [nAp + nLink, n, 'displacement strand', P.s2]];
}
function ssLegend(node, keys, shape) {
  const lg = document.createElement('div'); lg.className = 'legend';
  keys.forEach(([c, label]) => {
    const k = document.createElement('span'); k.className = 'k';
    if (c) { const s = document.createElement('span'); s.className = shape || 'sw';
             s.style.background = c; k.appendChild(s); }
    k.appendChild(document.createTextNode(label)); lg.appendChild(k);
  });
  node.appendChild(lg);
}
function fornaView(node, opt) {
  const Forna = fornaCtor(), seq = opt.seq, db = opt.db, n = seq.length;
  const host = document.createElement('div');
  host.id = 'ss-forna';                       /* the id every fornac CSS rule is scoped to */
  const W = Math.max(320, node.clientWidth || 720);
  const H = Math.max(360, Math.min(560, 200 + 3 * n));
  host.style.height = H + 'px';
  node.appendChild(host);
  dropFornacGlobalCss();
  const c = new Forna('#ss-forna', {initialSize:[W, H], labelInterval:10,
    allowPanningAndZooming:true, applyForce:true, transitionDuration:0});
  /* colour every nucleotide by construct region. fornac's customColors is
     {colorValues: {structureName: {1-based position: colour}}}, and '' matches any
     structure, so the RNA does not have to be named. */
  const cv = {}, regions = ssRegions(opt);
  regions.forEach(([a, b, label, col]) => { for (let i = a; i < b; i++) cv[i + 1] = col; });
  c.addCustomColors({colorValues: {'': cv}});
  /* the sequence is DNA and is passed through verbatim -- no T->U conversion, so the
     letters on screen are the letters you would order. fornac only needs the dot-bracket
     for the layout and never validates the alphabet. */
  c.addRNA(db, {sequence: seq, labelInterval: 10});
  c.changeColorScheme('custom');
  ssLegend(node, regions.map(([a, b, label, col]) => [col, label]).concat(
    [[null, pairs(db).length + ' base pairs · ' + n + ' nt'],
     [null, 'drag to pan · scroll to zoom · drag a nucleotide to pull the layout apart']]));
}
/* ---------- base-pair arc diagram (offline fallback) ---------- */
function arcView(node, opt, reason) {
  if (reason) {
    const w = document.createElement('div'); w.className = 'issues';
    w.style.margin = '0 0 12px'; w.textContent = reason; node.appendChild(w);
  }
  const seq = opt.seq, db = opt.db, n = seq.length, pr = pairs(db);
  const cellW = 14, m = {t:118, l:26, r:26, b:60};
  const W = m.l + n * cellW + m.r, H = m.t + m.b;
  const host = document.createElement('div'); host.style.overflowX = 'auto';
  node.appendChild(host);
  const svg = el('svg', {width:'100%', viewBox:'0 0 ' + W + ' ' + H, role:'img',
    'aria-label':'Base-pair arc diagram'}, host);
  const base = m.t, X = i => m.l + i * cellW + cellW / 2;

  /* construct regions: aptamer | poly-T linker | displacement strand */
  let nAp = n, nLink = 0, nDs = 0;
  if (opt.kind === 'switch') { nDs = opt.row.ds_len; nLink = opt.row.linker_len; nAp = n - nDs - nLink; }
  const regions = opt.kind === 'switch'
    ? [[0, nAp, 'aptamer', P.s1], [nAp, nAp + nLink, 'poly-T linker', P.s3],
       [nAp + nLink, n, 'displacement strand', P.s2]]
    : [[0, n, 'aptamer', P.s1]];
  regions.forEach(([a, b, label, col]) => {
    el('rect', {x:X(a) - cellW / 2 + 1, y:base + 20, width:(b - a) * cellW - 2, height:5,
                rx:2.5, fill:col}, svg);
    txt(svg, (X(a) + X(b - 1)) / 2, base + 40, label,
        {'text-anchor':'middle', 'font-size':10.5, fill:P.ink2});
  });
  /* randomised blocks on the parent portion */
  const rand = new Set(DATA.randomised);
  /* arcs above the sequence line */
  const maxSpan = pr.reduce((s, p) => Math.max(s, p[1] - p[0]), 1);
  pr.forEach(([a, b]) => {
    const x1 = X(a), x2 = X(b), r = (b - a) / maxSpan;
    const h = 18 + r * (m.t - 34);
    const cross = opt.kind === 'switch' && ((a < nAp) !== (b < nAp));
    el('path', {d:'M' + x1 + ' ' + base + ' C ' + x1 + ' ' + (base - h) + ' ' +
                   x2 + ' ' + (base - h) + ' ' + x2 + ' ' + base,
                fill:'none', stroke:cross ? P.s2 : P.s1, 'stroke-width':2,
                opacity:0.85, 'stroke-linecap':'round'}, svg);
  });
  el('line', {x1:m.l, x2:W - m.r, y1:base, y2:base, stroke:P.axis, 'shape-rendering':'crispEdges'}, svg);
  for (let i = 0; i < n; i++) {
    const t = txt(svg, X(i), base + 14, seq[i], {'text-anchor':'middle', 'font-size':10.5,
      fill:(i < nAp && rand.has(i)) ? P.ink : P.muted,
      'font-family':'ui-monospace, Consolas, monospace'});
    if (i < nAp && rand.has(i)) t.setAttribute('font-weight', '600');
  }
  /* per-position hover */
  const partner = new Array(n).fill(-1);
  pr.forEach(([a, b]) => { partner[a] = b; partner[b] = a; });
  for (let i = 0; i < n; i++) {
    const hit = el('rect', {x:X(i) - cellW / 2, y:8, width:cellW, height:base + 10,
                            fill:'transparent'}, svg);
    hit.addEventListener('pointermove', ev => showTip(ev, [
      {value:seq[i] + ' at ' + i, color:P.s1},
      {value:partner[i] >= 0 ? 'paired with ' + seq[partner[i]] + ' at ' + partner[i] : 'unpaired'},
      {value:(i < nAp ? (rand.has(i) ? 'aptamer, randomised block' : 'aptamer, fixed region')
                      : (i < nAp + nLink ? 'poly-T linker' : 'displacement strand'))}
    ]));
    hit.addEventListener('pointerleave', hideTip);
  }

  const keys = opt.kind === 'switch'
    ? [[P.s1, 'base pair within the aptamer'], [P.s2, 'DS ↔ aptamer pair (the switch)']]
    : [[P.s1, 'base pair']];
  ssLegend(node, keys.concat([[null, pr.length + ' pairs · ' + n + ' nt']]), 'ln');
}
/* the dot-bracket itself, under whichever renderer drew above it */
function textView(node, opt) {
  const det = document.createElement('details'); det.className = 'tv';
  const sm = document.createElement('summary');
  sm.textContent = 'dot-bracket + sequence (text view)';
  det.appendChild(sm);
  const pre = document.createElement('div'); pre.className = 'seq';
  pre.style.marginTop = '6px';
  pre.textContent = opt.seq + '\n' + opt.db;
  pre.style.whiteSpace = 'pre-wrap';
  det.appendChild(pre); node.appendChild(det);
}
const FORNA_MISSING =
  'fornac is not available — either its CDN script did not load or the generator could ' +
  'not read aptamer/vendor/fornac.css — so this panel is showing the built-in base-pair ' +
  'arc diagram instead. The structure is identical; only the layout differs. Every other ' +
  'view except the 3D panel works offline either way.';
function renderStructure(opt) {
  guard('ss-chart', node => {
    node.textContent = '';
    if (!opt) return empty(node, 'no structures available yet');
    if (fornaCtor()) {
      try {
        fornaView(node, opt);
      } catch (e) {
        node.textContent = '';
        arcView(node, opt, 'fornac loaded but could not draw this structure (' + e.message +
                           '); showing the base-pair arc diagram instead.');
      }
    } else {
      arcView(node, opt, FORNA_MISSING);
    }
    textView(node, opt);
  });
}
let ssOpts = [];
function buildStructurePicker() {
  const sel = document.getElementById('ss-select'); if (!sel) return;
  const keep = sel.selectedIndex >= 0 ? sel.options[sel.selectedIndex].textContent : null;
  ssOpts = structureOptions();
  if (!ssOpts.length) { sel.textContent = ''; renderStructure(null); return; }
  sel.textContent = '';
  ssOpts.forEach((o, i) => { const op = document.createElement('option');
    op.value = i; op.textContent = o.label; sel.appendChild(op); });
  const again = ssOpts.findIndex(o => o.label === keep);
  sel.value = again >= 0 ? again : 0;
  if (!sel.dataset.wired) {
    sel.addEventListener('change', () => renderStructure(ssOpts[+sel.value]));
    sel.dataset.wired = '1';
  }
  window.__ssRedraw = () => renderStructure(ssOpts[+sel.value || 0]);
  renderStructure(ssOpts[+sel.value || 0]);
}

/* =========================================================================
   VIEW 8 - 3D viewer (3Dmol.js, needs internet)
   ========================================================================= */
function render3D() {
  const host = document.getElementById('mol-viewer');
  const sel = document.getElementById('mol-select');
  if (!host || !sel) return;
  const structs = DATA.structures || [];
  if (!structs.length) { empty(host, 'no .pdb / .cif in aptamer/ yet — not yet generated'); return; }
  sel.textContent = '';
  structs.forEach((s, i) => { const o = document.createElement('option');
    o.value = i; o.textContent = s.name + ' — ' + s.label; sel.appendChild(o); });
  function draw() {
    if (!window.$3Dmol) {
      empty(host, '3D viewer could not load. 3Dmol.js is fetched from a CDN, so this one panel needs internet; every other view works offline.');
      return;
    }
    const s = structs[+sel.value || 0];
    host.textContent = '';
    const v = $3Dmol.createViewer(host, {backgroundColor: P.surface});
    v.addModel(s.text, s.format);
    v.setStyle({}, {cartoon:{color:P.s1, opacity:0.95}});
    v.setStyle({resn:['DA','DC','DG','DT','A','C','G','U','T']},
               {stick:{colorscheme:'default', color:P.s2, radius:0.22},
                cartoon:{color:P.s2, style:'trace'}});
    v.setStyle({hetflag:true}, {stick:{color:P.s4, radius:0.18}});
    v.zoomTo(); v.render();
    window.__molRedraw = draw;
  }
  sel.addEventListener('change', draw);
  if (window.$3Dmol) draw();
  else {
    const t = setInterval(() => { if (window.$3Dmol) { clearInterval(t); draw(); } }, 200);
    setTimeout(() => { clearInterval(t); if (!window.$3Dmol) draw(); }, 6000);
  }
}

/* ---------- boot + theme ---------- */
function renderAll() {
  P = PAL[mode()];
  renderParents(); renderOccupancy();
  applyFilters();
  if (window.__ssRedraw) window.__ssRedraw();
  if (window.__molRedraw) window.__molRedraw();
}
function boot() {
  dropFornacGlobalCss();
  buildFilterRow();
  renderParents(); renderOccupancy();
  applyFilters();
  buildStructurePicker();
  /* fornac normally executes before this script, but poll briefly in case the CDN is slow;
     if it never arrives the arc diagram that is already on screen simply stays. */
  if (!fornaCtor()) {
    const ft = setInterval(() => {
      if (fornaCtor()) { clearInterval(ft); if (window.__ssRedraw) window.__ssRedraw(); }
    }, 200);
    setTimeout(() => clearInterval(ft), 6000);
  }
  render3D();
  const btn = document.getElementById('theme-toggle');
  btn.addEventListener('click', () => {
    document.documentElement.dataset.theme = mode() === 'dark' ? 'light' : 'dark';
    btn.textContent = mode() === 'dark' ? 'Light mode' : 'Dark mode';
    renderAll();
  });
  btn.textContent = mode() === 'dark' ? 'Light mode' : 'Dark mode';
  matchMedia('(prefers-color-scheme: dark)').addEventListener('change', renderAll);
  let rt; addEventListener('resize', () => { clearTimeout(rt); rt = setTimeout(renderAll, 200); });
}
boot();
"""


def stat_tiles(data):
    """Small HTML for the stat row; the hero figure is the passing-construct count."""
    parents = data["parents"] or []
    switches = data["switches"] or []
    tiles = []

    def tile(label, value, unit="", hero=False):
        cls = "value hero" if hero else "value"
        u = f'<span class="unit">{unit}</span>' if unit else ""
        return (f'<div class="stat"><div class="label">{label}</div>'
                f'<div class="{cls}">{value}{u}</div></div>')

    if switches:
        tiles.append(tile("Switch constructs passing the design filters", len(switches), "", hero=True))
        best = min(switches, key=lambda r: r["kd_app_nM"])
        tiles.append(tile("Tightest K_D,app in the library", f"{best['kd_app_nM']:.1f}", "nM"))
        top = switches[0]
        tiles.append(tile("Top-ranked DS", top["ds"], f"window {top['window']}"))
    else:
        tiles.append(tile("Switch constructs", "—", "switches.csv not yet generated"))
    if parents:
        b = min(parents, key=lambda r: r["KD_M"])
        tiles.append(tile("Tightest parent K_D", f"{b['kd_nM']:.1f}", f"nM · {b['name']}"))
        flagged = sum(1 for p in parents if p.get("issues"))
        tiles.append(tile("Parents with uncertain nucleotides", flagged, f"of {len(parents)}"))
    return "\n".join(tiles)


def section(title, sub, body_id, extra=""):
    return (f'<section class="card"><header><h2>{title}</h2>'
            f'<p class="sub">{sub}</p></header>{extra}'
            f'<div id="{body_id}"></div></section>')


def build_html(data):
    forna_css = scoped_fornac_css()  # before the payload: it can append a generator note
    # Without the scoped stylesheet fornac would draw strokeless links and its own injected
    # sheet is deleted at runtime, so don't load it at all -- view 7 then uses the arcs.
    forna_script = (
        '<!-- fornac 1.1.8 draws view 7. Note the path: dist/scripts/fornac.js, not\n'
        '     dist/fornac.js (that one 404s). d3 3.5.13 is bundled inside it and jQuery is\n'
        '     never referenced, so it needs no companion script. If this fails to load,\n'
        '     view 7 falls back to the built-in base-pair arc diagram. -->\n'
        '<script src="https://cdn.jsdelivr.net/npm/fornac@1.1.8/dist/scripts/fornac.js">'
        "</script>"
    ) if forna_css else "<!-- fornac not loaded: vendor/fornac.css was unreadable -->"
    payload = json.dumps(data, allow_nan=False).replace("</", "<\\/")
    n_sw = len(data["switches"] or [])
    n_par = len(data["parents"] or [])
    complexes = [s for s in data["structures"] if s["kind"] == "complex"]
    notes = data["notes"]
    note_html = ""
    if notes:
        items = "".join(f"<li>{n}</li>" for n in notes)
        note_html = (f'<div class="issues" style="margin-bottom:20px"><b>Generator notes</b>'
                     f'<ul>{items}</ul></div>')

    complex_note = (
        f"{len(complexes)} co-folded complex(es) found." if complexes else
        "No co-folded aptamer–IL-6 complexes have been generated yet, so only the "
        "IL-6 receptor structure is available here."
    )

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>IL-6 aptamer switch design — dashboard</title>
<style data-dash>{CSS}
/* vendor/fornac.css, every selector rewritten to sit under {FORNA_SCOPE} by
   dashboard.py:scoped_fornac_css(). Unscoped it would resize every chart on the page. */
{forna_css}</style>
</head>
<body data-palette="#2a78d6,#eb6834,#1baf7a,#eda100">
<div class="viz-root">

<header class="top">
  <div class="grow">
    <h1>Structure-switching IL-6 aptamers</h1>
    <p class="sub">{n_par} parent aptamers · {n_sw} ranked ISD switch constructs ·
       generated by <code>aptamer/dashboard.py</code>. All data is inlined — works offline
       except the 3D panel.</p>
  </div>
  <button id="theme-toggle" type="button">Dark mode</button>
</header>

{note_html}

<div class="stats">{stat_tiles(data)}</div>

{section("1 · Parent aptamers",
         "Published Neomer candidates, reconstructed and folded. Reconstruction warnings flag a "
         "genuinely uncertain nucleotide — those constructs must not be ordered blind.",
         "parents-body")}

<section class="card">
  <header><h2>2 · Fractional occupancy in wound fluid</h2>
  <p class="sub">Fraction of aptamer bound vs IL-6 concentration, from each published K_D.
  Vertical rules mark the wound-fluid medians; the wash spans post-surgical to infected;
  the amber rule is the skin-model reading (5 ng/mL).</p></header>
  <div id="occ-chart"></div>
  <details class="tv" open><summary>Table view — occupancy at each landmark</summary>
    <div id="occ-table" style="margin-top:8px"></div></details>
</section>

<h2 style="margin:28px 0 10px">Switch library</h2>
<p class="sub" style="margin-bottom:12px">These filters scope everything below them — the
table, the scatter, the heatmap, the coverage track and the switch structures all re-render
against the same slice.</p>
<div class="filters" id="sw-filters"></div>

{section("3 · Ranked switch constructs",
         "Default order is the pipeline's ranking: randomised nucleotides covered, then "
         "engagement, then short DS. Click any header to sort.",
         "sw-table")}

<div class="grid2">
  {section("4 · Gain vs affinity",
           "Every construct: stabilisation energy against the affinity it costs. Colour = "
           "randomised nucleotides the DS covers; size = engagement.",
           "sw-scatter")}
  {section("5 · Linker × DS length",
           "Best (lowest) K_D,app reachable at each combination — the two knobs the "
           "design turns.",
           "sw-heatmap")}
</div>

{section("6 · Where the DS windows land",
         "The parent sequence as an axis. Bars count the passing constructs whose displacement "
         "strand targets each position; the amber wash marks the randomised library blocks.",
         "sw-coverage")}

<section class="card">
  <header><h2>7 · Secondary structure</h2>
  <p class="sub">The real 2D layout, drawn by <a href="https://github.com/ViennaRNA/fornac"
  target="_blank" rel="noopener">fornac</a> from the dot-bracket the pipeline folded —
  stems, loops and junctions where they actually sit. Nothing is re-folded here and the DNA
  sequence is shown verbatim (no T→U conversion). Nucleotides are coloured by construct
  region: aptamer, poly-T linker, displacement strand. Drag to pan, scroll to zoom.
  fornac comes from a CDN; without it this panel falls back to the built-in base-pair arc
  diagram and says so, so it still works offline.</p></header>
  <div style="margin-bottom:10px">
    <label class="sub" for="ss-select">Structure&nbsp;</label>
    <select id="ss-select" style="max-width:min(560px,100%)"></select>
  </div>
  <div id="ss-chart"></div>
</section>

<section class="card">
  <header><h2>8 · 3D structures</h2>
  <p class="sub">{complex_note} This panel loads 3Dmol.js from a CDN, so it is the one view
  that needs internet.</p></header>
  <div style="margin-bottom:10px">
    <label class="sub" for="mol-select">Structure&nbsp;</label>
    <select id="mol-select"></select>
  </div>
  <div class="viewer" id="mol-viewer"></div>
</section>

<footer class="foot">
Regenerate with <code>python aptamer/dashboard.py</code>. Never hand-edit dashboard.html.
Values are read as the pipeline produced them; only residence time (1/k_off) and fractional
occupancy are derived for display.
</footer>
</div>

<script id="dash-data" type="application/json">{payload}</script>
<script src="https://cdn.jsdelivr.net/npm/3dmol@2.4.2/build/3Dmol-min.js"></script>
{forna_script}
<script>{JS}</script>
</body>
</html>
"""


def main():
    data = build_data()
    html = build_html(data)
    OUT.write_text(html, encoding="utf-8")
    kb = OUT.stat().st_size / 1024
    print(f"parents.json   : {len(data['parents'] or [])} aptamers")
    print(f"switches.csv   : {len(data['switches'] or [])} constructs")
    for s in data["structures"]:
        print(f"{s['name']:<15}: {s['atoms']} atoms ({s['kind']})")
    if not any(s["kind"] == "complex" for s in data["structures"]):
        print("co-folded      : none yet -> 3D panel shows the receptor only")
    for n in data["notes"]:
        print(f"note           : {n}")
    print(f"\nwrote {OUT} ({kb:.0f} KB)")


if __name__ == "__main__":
    main()
