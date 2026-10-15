#!/usr/bin/env python3
# =====================================================================================================================
#  digitise_master.py -- MASTER FILE of the figure-digitisation protocol (Runs/R/digitisation/, 14 Oct 2026)
# =====================================================================================================================
#
#  WHAT THIS FILE IS
#  -----------------
#  1. The single shared LIBRARY of every function used to turn published figures (and raw data files) into the
#     time-series format read by the change-point pipelines (PELT+CROPS and segmented regression). Every per-source
#     script (in source_papers/ or utilities/) adds this folder to the Python path and imports `from digitise_master import *`,
#     so the functions below are the only copy.
#  2. The REGISTRY of all digitisation scripts (section 10): source paper, figure, method, output series, run order.
#  3. A command-line RUNNER (section 11):
#        python3 R/digitisation/digitise_master.py --list               table of all scripts (run order, sources, outputs)
#        python3 R/digitisation/digitise_master.py --which <series>     which script makes a data_digitized/ series
#        python3 R/digitisation/digitise_master.py --run all            re-run every script in the correct order
#        python3 R/digitisation/digitise_master.py --run guo2024 ...    re-run selected scripts (or --run group:<group>)
#        python3 R/digitisation/digitise_master.py --check              re-run everything in a scratch copy and verify
#                                                                       that data_digitized/ is reproduced exactly
#     Always run from the Runs root (the runner changes to it automatically).
#
#  FOLDER LAYOUT
#  -------------
#     R/digitisation/digitise_master.py      this file (library + registry + runner)
#     R/digitisation/README.md               overview, run order, inputs/outputs, old-to-new file names
#     R/digitisation/source_papers/*.py      one script per source paper (37 scripts)
#     R/digitisation/utilities/*.py          steps not tied to one paper: raster extraction (first), QC panels,
#                                            PI age cutoffs (always last)
#
#  PROTOCOL IN BRIEF
#  -----------------
#  * One script per source paper (file named after the paper, e.g. foth2021.py); it reads the figure from
#    data_raw/articles_PDF/ and writes one CSV per series to data_digitized/<Clade>_<Parameter>[_<data_type>].csv.
#  * Output format (one row per time point; section 2):
#        Clade, Parameter, Age_Ma, Value, Value_err, Err_lo, Err_hi, Err_type, Source, Age_cal_sd
#    Err_type = "95CI" (Err_lo/Err_hi = published 95% interval), "SE" (+/-1 SE), "Poisson" (counts without published
#    error: Value +/- sqrt(Value)) or "none" (no published error; only the digitisation error remains).
#    Age_cal_sd = SD of the age calibration (Myr; linear-fit RMSE combined with half a pixel), used by both pipelines
#    to jitter ages in the Monte-Carlo replicates, so digitisation error propagates into the onset/offset ranges.
#  * Two kinds of source figure:
#      VECTOR figures (lines stored as PDF drawing paths): values are read from the path coordinates with PyMuPDF
#        (section 4). Exact apart from the axis calibration.
#      RASTER figures (embedded images): the image is extracted from the PDF (extract_raster_figures.py, or saved
#        in data_digitized/digitisation/*_img/), the curve is isolated by colour / darkness masks and traced per pixel
#        column (section 5).
#    In both cases the axes are calibrated by a least-squares LINEAR fit of pixel (or PDF point) position against
#    the tick values (section 3); using all ticks (not just two clicked points, as in metaDigitise) gives a residual
#    that measures the calibration error.
#  * Dense traced curves are averaged onto a regular 1-Myr grid (section 6); node-based raw rate files are summarised
#    in 5-Myr bins (section 7).
#  * Every series gets a QC panel (source crop | digitised series) in data_raw/<parameter folder>/ (section 8).
#  * PI age cutoffs (data_digitized/digitisation/age_cutoffs.csv) are applied LAST by age_cutoffs.py.
#  * Choices and caveats per series: data_digitized/digitisation/notes_{vert,invert,massext}.csv (appended to the
#    pipeline Flags).
#
#  OUTPUT "DIALECTS"
#  -----------------
#  The scripts were written in three batches (vertebrates, invertebrates, mass-extinction sources). Their output
#  writers differ only in rounding and in the order of their arguments; all three are kept (save_vert, save_invert,
#  save_massext) so that every existing series is reproduced byte for byte. New scripts should use save_series().
#
#  HOW TO ADD A NEW FIGURE
#  -----------------------
#  1. Copy the closest existing script in source_papers/ (vector: guo2024.py / navalon2022.py; raster line:
#     whalen2020.py / rawson2022.py; raster points with error bars: lamsdell_selden2016.py / wright2017.py; raw node
#     data: raw_data_convert.py), rename it after the paper and keep its two import lines (sys.path + digitise_master).
#  2. Calibrate the axes with lin() on as many ticks as possible; extract the curve; build a DataFrame with
#     Age_Ma, Value (+ Err_lo, Err_hi, Value_err when the figure shows an interval).
#  3. Write it with save_series(df, clade, parameter, err_type, source_text, age_sd).
#  4. Add a QC panel (qc_panel) to the matching utilities/qc_panels_*.py and a line to REGISTRY (section 10); add a
#     note to the notes_*.csv file if a choice was made (line colour, axis, excluded points ...).
#  5. Run `digitise_master.py --run <script> qc_panels_<group> age_cutoffs` and then `--check`.
#
#  Requirements: Python 3 with numpy, pandas, scipy, pillow, matplotlib, pymupdf (and scikit-image for wright2017.py).
# =====================================================================================================================
import os, re, sys, json, shutil, filecmp, hashlib, tempfile, subprocess, argparse
import numpy as np, pandas as pd

class _Missing:
    """Placeholder for an optional package that is not installed. Importing this file (e.g. from
    R/prep_mass_extinction.py, which only needs node_bins/COLS/REGISTRY) then works without it; the first call that
    actually needs the package stops with an instruction to install it."""
    def __init__(self, pkg): self._pkg = pkg
    def __getattr__(self, name):
        raise ImportError(f"this step needs the Python package '{self._pkg}': pip install {self._pkg}")
    def __call__(self, *a, **k): self.__getattr__("__call__")

try: import pymupdf                       # reading PDF vector paths, extracting embedded images, rendering crops
except ImportError: pymupdf = _Missing("pymupdf")
try: from PIL import Image                # reading/writing raster images
except ImportError: Image = _Missing("pillow")
try: from scipy import ndimage as ndi     # connected components / morphological opening of raster masks
except ImportError: ndi = _Missing("scipy")

# ---------------------------------------------------------------------------------------------------------------------
# 1. PATHS AND CONSTANTS (relative to the Runs root; every script is run from there)
# ---------------------------------------------------------------------------------------------------------------------
PDF = "data_raw/articles_PDF/"            # source PDFs and raw data files listed in the keys
A = PDF                                   # alias used by the vertebrate scripts
OUT = "data_digitized/"                   # pipeline-format series (vertebrate / mass-extinction dialect, with "/")
OUT_INVERT = "data_digitized"             # same folder, invertebrate dialect (paths written as f"{OUT}/...")
IMG_INVERT = "data_digitized/digitisation/invert_img/"    # rasters of the invertebrate figures (pdfimages extraction)
IMG_VERT = "data_digitized/digitisation/vert_img/"        # rasters of the vertebrate figures (extract_raster_figures.py)
IMG_MASSEXT = "data_digitized/digitisation/massext_img/"  # rasters of the mass-extinction figures (extracted by each script)
IMG = IMG_MASSEXT                         # default image folder (mass-extinction scripts); others re-bind IMG
QC_VERT = "data_digitized/digitisation/vert_qc/"          # trace-over-source QC images of the vertebrate raster scripts
COLS = ["Clade", "Parameter", "Age_Ma", "Value", "Value_err", "Err_lo", "Err_hi", "Err_type", "Source", "Age_cal_sd"]
HERE = os.path.dirname(os.path.abspath(__file__))         # .../Runs/R/digitisation
SCRIPTS = os.path.join(HERE, "source_papers")             # .../Runs/R/digitisation/source_papers (per-source scripts)
UTILS = os.path.join(HERE, "utilities")                   # .../Runs/R/digitisation/utilities (extraction, QC, cutoffs)
UTIL_GROUPS = ("utility", "qc", "cutoff")                 # registry groups whose scripts live in utilities/
RUNS = os.path.dirname(os.path.dirname(HERE))             # .../Runs

def _ensure_dirs():
    """Output folders the scripts expect (created only when the working directory is the Runs root)."""
    if os.path.isdir("data_raw"):
        for d in (OUT, IMG_INVERT, IMG_MASSEXT): os.makedirs(d, exist_ok=True)
_ensure_dirs()


# ---------------------------------------------------------------------------------------------------------------------
# 2. OUTPUT: writing pipeline-format series
# ---------------------------------------------------------------------------------------------------------------------
def save_series(df, clade, param, err_type, source, age_sd, suffix="", digits=6):
    """Write one series to data_digitized/<clade>_<param>[_<suffix>].csv in the pipeline format (COLS).
    df      : DataFrame with at least Age_Ma and Value; Value_err / Err_lo / Err_hi are optional (if absent:
              Value_err = 0 and Err_lo = Err_hi = Value, i.e. no published error).
    err_type: "95CI", "SE", "Poisson" or "none" (see the header).
    source  : free text recorded in every row (reference, figure, line, method); shown in the pipeline tables.
    age_sd  : age-calibration SD in Myr (RMSE of the x-axis fit combined with half a pixel), stored as Age_cal_sd.
    suffix  : data_type slug when a clade x parameter has several series (e.g. "FBDD", "continuous_data_GMM").
    Rows are sorted from old to young. This is the writer to use for NEW scripts (= save_vert)."""
    df = df.copy(); df["Clade"] = clade; df["Parameter"] = param; df["Err_type"] = err_type; df["Source"] = source; df["Age_cal_sd"] = round(age_sd, 3)
    for c in ["Value_err", "Err_lo", "Err_hi"]:
        if c not in df: df[c] = 0.0 if c == "Value_err" else df["Value"]
    name = f"{clade}_{param}{('_' + suffix) if suffix else ''}.csv"
    df = df[COLS].sort_values("Age_Ma", ascending=False); df.round(digits).to_csv(OUT + name, index=False)
    print(name, len(df), "rows; ages %.1f-%.1f; value %.4g-%.4g" % (df.Age_Ma.max(), df.Age_Ma.min(), df.Value.min(), df.Value.max()))
    return df

def save_vert(df, clade, param, err_type, source, age_sd, suffix=""):
    """Vertebrate dialect (formerly vcommon.save): identical to save_series (values rounded to 6 decimals)."""
    return save_series(df, clade, param, err_type, source, age_sd, suffix, digits=6)

def save_invert(df, clade, param, age_sd, err_type, source):
    """Invertebrate dialect (formerly common.save). NOTE the argument order (age_sd BEFORE err_type) and the
    rounding to 4 decimals; Clade/Parameter are inserted as the first columns."""
    df = df.copy(); df.insert(0, "Parameter", param); df.insert(0, "Clade", clade)
    for c in ["Value_err", "Err_lo", "Err_hi"]:
        if c not in df: df[c] = 0.0 if c == "Value_err" else df["Value"]
    df["Err_type"] = err_type; df["Source"] = source; df["Age_cal_sd"] = round(age_sd, 3)
    df = df.sort_values("Age_Ma", ascending=False)
    df[COLS].round(4).to_csv(f"{OUT_INVERT}/{clade}_{param}.csv", index=False); print(clade, param, len(df)); return df[COLS]

def save_massext(df, clade, param, err_type, source, age_sd):
    """Mass-extinction dialect (formerly mcommon.save): as save_series but rounded to 5 decimals, no suffix."""
    df = df.copy(); df["Clade"] = clade; df["Parameter"] = param; df["Err_type"] = err_type
    df["Source"] = source; df["Age_cal_sd"] = round(float(age_sd), 3)
    for c in ["Value_err", "Err_lo", "Err_hi"]:
        if c not in df: df[c] = 0.0 if c == "Value_err" else df["Value"]
    df = df[COLS].sort_values("Age_Ma", ascending=False)
    df.round(5).to_csv(f"{OUT}{clade}_{param}.csv", index=False)
    print(f"{clade}_{param}.csv  n={len(df)}  ages {df.Age_Ma.max():.1f}-{df.Age_Ma.min():.1f}  values {df.Value.min():.4g}-{df.Value.max():.4g}")
    return df

def save_table(df, name):
    """Writer of the raw-data conversions (formerly raw_convert.save): df already has all COLS; file name given."""
    df = df[COLS].sort_values("Age_Ma", ascending=False); df.round(5).to_csv(OUT + name, index=False)
    print(name, len(df), "rows; ages", df.Age_Ma.max().round(1), "-", df.Age_Ma.min().round(1))

def poisson(d):
    """Counts without published error bars (genera, species): Value_err = sqrt(Value), Err_lo/hi = Value -/+ sqrt(Value).
    Negative digitised counts are clipped to 0 first. Used with err_type = "Poisson"."""
    d["Value"] = d["Value"].clip(lower=0)
    d["Value_err"] = np.sqrt(d["Value"]); d["Err_lo"] = d["Value"] - d["Value_err"]; d["Err_hi"] = d["Value"] + d["Value_err"]
    return d


# ---------------------------------------------------------------------------------------------------------------------
# 3. AXIS CALIBRATION
# ---------------------------------------------------------------------------------------------------------------------
def lin(pix, val):
    """Linear axis calibration: least-squares fit value = a + b * pixel over ALL ticks given.
    pix: tick positions (PDF points for vector figures, pixel rows/columns for rasters); val: their axis values.
    Returns (f, rmse): f converts positions to axis values; rmse = residual SD in axis units, i.e. the calibration
    error (for the age axis it becomes Age_cal_sd, usually combined with half a pixel: np.hypot(rmse, 0.5 / px_per_Myr)).
    Typical use:  X, sx = lin([94.3, 136.7, 178.3], [160, 140, 120]);  ages = X(xs)"""
    pix = np.asarray(pix, float); val = np.asarray(val, float); b, a = np.polyfit(pix, val, 1); res = val - (a + b * pix)
    return (lambda p: a + b * np.asarray(p, float)), float(np.sqrt(np.mean(res ** 2)))


# ---------------------------------------------------------------------------------------------------------------------
# 4. VECTOR FIGURES: reading PDF drawing paths (PyMuPDF page.get_drawings())
#    A drawing is a dict with 'type' ('s' stroke, 'f' fill, 'fs' both), 'color' / 'fill' (RGB 0-1), 'width',
#    'dashes', 'rect' (bounding box) and 'items' (path operators: 'l' line, 'c' Bezier curve, 're' rectangle,
#    'qu' quad, each followed by points). Lines are identified by colour, width and dash pattern; ticks are short
#    black strokes next to the axis; markers are small filled circles (all items 'c').
# ---------------------------------------------------------------------------------------------------------------------
def path_pts(g):
    """All points (x, y) of a drawing path, in drawing order (PDF points; y grows downwards)."""
    out = []
    for it in g["items"]:
        for q in it[1:]:
            if hasattr(q, "x"): out.append((q.x, q.y))
    return out
pts = path_pts                            # alias used by the vertebrate scripts

def dedupe(P, tol=1e-3):
    """Drop consecutive duplicate points of a path (segment ends repeat the next segment's start)."""
    o = []
    for x, y in P:
        if not o or abs(o[-1][0] - x) > tol or abs(o[-1][1] - y) > tol: o.append((x, y))
    return o

def segs(g):
    """Path as a list of straight segments (x0, y0, x1, y1); curves are approximated by their chords, rectangles and
    quads by their outlines. Used to intersect lines / filled bands with vertical lines (crossings)."""
    S = []
    for it in g["items"]:
        P = [q for q in it[1:] if hasattr(q, "x")]
        if it[0] == "re":
            r = it[1]; P = [pymupdf.Point(r.x0, r.y0), pymupdf.Point(r.x1, r.y0), pymupdf.Point(r.x1, r.y1), pymupdf.Point(r.x0, r.y1), pymupdf.Point(r.x0, r.y0)]
        if it[0] == "qu":
            q = it[1]; P = [q.ul, q.ur, q.lr, q.ll, q.ul]
        for a, b in zip(P[:-1], P[1:]): S.append((a.x, a.y, b.x, b.y))
    return S

def crossings(S, x):
    """y values where the segments S cross the vertical line at x (linear interpolation within each segment).
    For a line this gives one y; for a filled band (closed polygon) the upper and lower edges (min / max y)."""
    ys = []
    for x0, y0, x1, y1 in S:
        if (x0 - x) * (x1 - x) <= 0 and x0 != x1: ys.append(y0 + (y1 - y0) * (x - x0) / (x1 - x0))
    return ys

def ys_at(S, x, tol=0.3):
    """crossings() that also works exactly at the ends of a path: if no segment strictly spans x, the end vertices
    lying within tol of x are used (needed when sampling a line at its first / last vertex)."""
    ys = crossings(S, x)
    if not ys:
        ys = [y0 for x0, y0, x1, y1 in S if abs(x0 - x) <= tol] + [y1 for x0, y0, x1, y1 in S if abs(x1 - x) <= tol]
    return ys

def poly_profile(groups, xs, mode="center"):
    """Vertical profile of filled polygons (list of drawings) at each x: returns rows (x, centre, min_y, max_y).
    Centre = mid-point of the outline (thick lines drawn as filled outlines); min/max = band edges."""
    S = [s for g in groups for s in segs(g)]; out = []
    for x in xs:
        ys = crossings(S, x)
        if len(ys) < 2: out.append((x, np.nan, np.nan, np.nan)); continue
        out.append((x, (min(ys) + max(ys)) / 2, min(ys), max(ys)))
    return np.array(out)

def circle_centres(drawings, color, width=None, box=None):
    """Centres of circular markers (paths made only of Bezier curves) of a given stroke/fill colour (RGB 0-1,
    tolerance 0.02), optional stroke width and bounding box (x0, y0, x1, y1). Sorted by x."""
    out = []
    for g in drawings:
        if g["type"] != "s" and g["type"] != "fs" and g["type"] != "f": continue
        col = g.get("color") if g["type"] != "f" else g.get("fill")
        if col is None or np.max(np.abs(np.array(col) - np.array(color))) > 0.02: continue
        if width is not None and abs((g.get("width") or 0) - width) > 0.02: continue
        if not all(i[0] == "c" for i in g["items"]): continue
        r = g["rect"]
        if box and not (box[0] <= r.x0 and r.x1 <= box[2] and box[1] <= r.y0 and r.y1 <= box[3]): continue
        out.append(((r.x0 + r.x1) / 2, (r.y0 + r.y1) / 2))
    return sorted(out)

def line_and_band(dr, col, linebox, bandcol=None, width=1.87):
    """A coloured line plus its same-coloured confidence band (filled polygon) inside linebox (x0, y0, x1, y1).
    Returns rows (x, y_line, y_upper, y_lower) at every line vertex (band edges read at the same x; if no band
    point is found the line value is used). Used for Novack-Gottshall 2022 (line width 1.87)."""
    line = None; band = None
    for g in dr:
        c = g.get("color") if g["type"] == "s" else g.get("fill")
        if c is None or max(abs(a - b) for a, b in zip(c, col)) > 0.01: continue
        r = g["rect"]
        if not (r.x1 > linebox[0] and r.x0 < linebox[2] and r.y0 >= linebox[1] - 15 and r.y1 <= linebox[3] + (45 if g["type"] == "f" else 5)): continue
        if g["type"] == "s" and abs(g["width"] - width) < 0.05: line = g
        if g["type"] == "f": band = g
    L = []
    for x, y in path_pts(line):
        if not L or abs(L[-1][0] - x) > 1e-3 or abs(L[-1][1] - y) > 1e-3: L.append((x, y))
    B = path_pts(band) if band else []
    rows = []
    for x, y in L:
        if x < linebox[0] - 0.01 or x > linebox[2] + 0.01: continue
        ys = [by for bx, by in B if abs(bx - x) < 0.06]
        up = min(ys) if ys else y; lo = max(ys) if ys else y
        rows.append((x, y, up, lo))
    d = {}                                   # one row per x
    for r in rows: d.setdefault(round(r[0], 2), r)
    return [d[k] for k in sorted(d)]


# ---------------------------------------------------------------------------------------------------------------------
# 5. RASTER FIGURES: masks, tick detection and per-column tracing (arrays are rows x columns x RGB, ints 0-255)
# ---------------------------------------------------------------------------------------------------------------------
def load(f):
    """Image file -> int array (rows, columns, 3)."""
    return np.asarray(Image.open(f).convert("RGB")).astype(int)

def dark(a, t=90):
    """Dark pixels (all channels < t): black lines, axes, ticks, text."""
    return (a.max(2) < t)

def axis_rows(a, minlen=0.5, t=90):
    """Rows that are dark over more than minlen of the image width (horizontal axes / frames)."""
    d = dark(a, t); w = a.shape[1]; return [r for r in range(a.shape[0]) if d[r].sum() > minlen * w]

def axis_cols(a, minlen=0.5, t=90):
    """Columns that are dark over more than minlen of the image height (vertical axes / frames)."""
    d = dark(a, t); h = a.shape[0]; return [c for c in range(a.shape[1]) if d[:, c].sum() > minlen * h]

def runs(idx):
    """Group sorted integer indices into consecutive runs; returns (centre, first, last) per run.
    E.g. the columns of a tick mark, or the rows covered by a line in one column."""
    idx = sorted(idx); out = []
    for i in idx:
        if out and i == out[-1][1] + 1: out[-1][1] = i
        else: out.append([i, i])
    return [((a + b) / 2, a, b) for a, b in out]

def xticks(a, row, below=(2, 6), t=110):
    """Tick marks under a horizontal axis at image row `row`: columns dark over rows row+below[0] .. row+below[1]."""
    d = dark(a, t); cols = [c for c in range(a.shape[1]) if d[row + below[0]:row + below[1], c].all()]
    return runs(cols)

def yticks(a, col, left=(2, 6), t=110):
    """Tick marks left of a vertical axis at image column `col`."""
    d = dark(a, t); rows = [r for r in range(a.shape[0]) if d[r, col - left[1]:col - left[0]].all()]
    return runs(rows)

def colmask(a, ref, tol=60):
    """Pixels whose RGB is within `tol` (sum of absolute channel differences) of the reference colour ref."""
    return (np.abs(a - np.array(ref)).sum(2) < tol)

def trace(mask, c0=0, c1=None, minpix=2):
    """Per-column trace of a line mask: for every column with >= minpix line pixels, keep the LARGEST contiguous
    run and return rows (column, centre_row, first_row, last_row)."""
    c1 = c1 or mask.shape[1]; out = []
    for c in range(c0, c1):
        r = np.where(mask[:, c])[0]
        if len(r) >= minpix:
            rr = runs(r); big = max(rr, key=lambda x: x[2] - x[1]); out.append((c, big[0], big[1], big[2]))
    return np.array(out, float)

def line_mask(a, ref, tol, box, min_comp=40, keep="largest", excl=()):
    """Colour mask of one line restricted to the plot box (x0, y0, x1, y1), with optional excluded rectangles
    (legend, labels). Keeps the largest connected component (keep="largest") or all components >= min_comp pixels."""
    m = colmask(a, ref, tol); x0, y0, x1, y1 = box; z = np.zeros_like(m); z[y0:y1, x0:x1] = True; m &= z
    for ex0, ey0, ex1, ey1 in excl: m[ey0:ey1, ex0:ex1] = False
    lab, n = ndi.label(m)
    if n == 0: return m
    sz = ndi.sum(m, lab, range(1, n + 1))
    ids = [np.argmax(sz) + 1] if keep == "largest" else [i + 1 for i in range(n) if sz[i] >= min_comp]
    return np.isin(lab, ids)

def centre_trace(m, trim=True):
    """Per-column centre of the line mask (largest run per column); the two ends are trimmed by half the line width
    (the rounded line caps would otherwise bias the first / last values). Returns (trace, line_width_px)."""
    t = trace(m, minpix=1)
    if len(t) == 0: return t
    w = np.median(t[:, 3] - t[:, 2] + 1)
    if trim:
        k = int(round(w / 2)); t = t[k:len(t) - k] if len(t) > 2 * k else t
    return t, w

def interp_gaps(cols, vals):
    """Fill missing columns (line hidden behind other lines or labels) by linear interpolation."""
    full = np.arange(cols.min(), cols.max() + 1); return full, np.interp(full, cols, vals)


# ---------------------------------------------------------------------------------------------------------------------
# 6. RESAMPLING DENSE CURVES
# ---------------------------------------------------------------------------------------------------------------------
def to_series(ages, vals, step=1.0):
    """Average a dense polyline / pixel trace onto a regular age grid (bins of `step` Myr centred on multiples of
    step): one value per bin = mean of the points in it, at their mean age. Keeps traced curves to a manageable,
    evenly spaced series (the pipelines resample to 1 Myr anyway)."""
    d = pd.DataFrame({"a": ages, "v": vals}).sort_values("a")
    d["b"] = np.round(d.a / step) * step
    return d.groupby("b").agg(Age_Ma=("a", "mean"), Value=("v", "mean")).reset_index(drop=True)


# ---------------------------------------------------------------------------------------------------------------------
# 7. RAW NODE-BASED RATE FILES (one row per branch / node)
# ---------------------------------------------------------------------------------------------------------------------
def node_bins(f, age, rate, clade, param, src, bw=5.0, drop_zero=False, max_age=None, anchor=0.0):
    """Branch/node rate file (data_raw/articles_PDF/<f>) -> binned series.
    Nodes are grouped in bins of bw Myr by node age (column `age`); per bin: Value = MEDIAN rate (robust to the long
    right tail of relative rates), Age_Ma = mean node age, Value_err = SE = 1.2533 * pooled robust SD / sqrt(n)
    (robust SD = 1.4826 * MAD of each bin with >= 3 nodes, pooled as the root mean square; 1.2533 = SE factor of a
    median). Err_type = "SE".
    drop_zero: drop extant tips (age 0).  max_age: drop nodes older than max_age BEFORE binning (PI age cutoffs).
    anchor: bin edges at anchor - k*bw. 0 (default) = edges at multiples of bw (clade-origin series); the
    mass-extinction series use anchor = extinction age, so the first bin covers the first bw Myr after the event."""
    d = pd.read_csv(PDF + f)[[age, rate]].dropna(); d.columns = ["age", "rate"]
    if drop_zero: d = d[d.age > 0]
    if max_age is not None: d = d[d.age <= max_age + 1e-9]
    d["bin"] = np.floor((d.age - anchor) / bw)
    mad_sd = d.groupby("bin").rate.apply(lambda x: 1.4826 * np.median(np.abs(x - np.median(x))) if len(x) >= 3 else np.nan)
    pooled = float(np.sqrt(np.nanmean(mad_sd.values ** 2)))
    g = d.groupby("bin").agg(Age_Ma=("age", "mean"), Value=("rate", "median"), n=("rate", "size")).reset_index()
    g["Value_err"] = 1.2533 * pooled / np.sqrt(g.n); g["Err_lo"] = g.Value - g.Value_err; g["Err_hi"] = g.Value + g.Value_err
    g["Clade"] = clade; g["Parameter"] = param; g["Err_type"] = "SE"; g["Age_cal_sd"] = 0.0
    g["Source"] = (f"{src} (raw data {f}: {len(d)} nodes" + (f" <= {max_age:g} Ma" if max_age is not None else "") +
                   f", ages={age}, rates={rate}; {bw:g}-Myr bin medians" + (f" (bins anchored at {anchor:g} Ma)" if anchor else "") +
                   f", SE from pooled robust SD={pooled:.3f})")
    return g[COLS + ["n"]]


# ---------------------------------------------------------------------------------------------------------------------
# 8. QC PANELS: source crop | digitised series, saved next to the crop in data_raw/<parameter folder>/
# ---------------------------------------------------------------------------------------------------------------------
FOLD = {"phenotypic_rate": "phenotypic_rates", "phenotypic_disparity": "phenotypic_disparity", "taxonomic_diversity": "taxonomic_diversity"}

def pdf_crop(pdf, page, clip, path=None, dpi=400):
    """Render a rectangle (PDF points: x0, y0, x1, y1) of page `page` (0-based) of data_raw/articles_PDF/<pdf> at
    `dpi`; saved to `path` if given (QC crops of vector figures)."""
    pix = pymupdf.open(PDF + pdf)[page].get_pixmap(dpi=dpi, clip=pymupdf.Rect(*clip))
    if path: pix.save(path)
    return pix

def raster_crop(src, box, path):
    """Crop an extracted raster (pixel box x0, y0, x1, y1; None = whole image) and save it (QC crops)."""
    im = Image.open(src).convert("RGB"); (im.crop(box) if box else im).save(path)

def qc_panel(clade, param, crop_path, series_csv, title_chars=80):
    """Two-panel QC figure: left = source crop, right = digitised series (Value with its Err_lo-Err_hi band) on a
    reversed age axis. Written as <crop>_QC_digitised.png. Used by qc_panels_invert.py and qc_panels_massext.py
    (qc_panels_vert.py draws a variant: red line, grey band)."""
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    d = pd.read_csv(series_csv); im = Image.open(crop_path)
    fig, ax = plt.subplots(1, 2, figsize=(13, 4.6), gridspec_kw=dict(width_ratios=[1, 1.1]))
    ax[0].imshow(im); ax[0].axis("off"); ax[0].set_title("Source: " + d.Source[0][:title_chars], fontsize=7)
    ax[1].fill_between(d.Age_Ma, d.Err_lo, d.Err_hi, color="tab:blue", alpha=.2, lw=0)
    ax[1].plot(d.Age_Ma, d.Value, "o-", ms=3, color="tab:blue", lw=1)
    ax[1].invert_xaxis(); ax[1].set_xlabel("Age (Ma)"); ax[1].set_ylabel(param)
    ax[1].set_title(f"Digitised: {clade} {param} (n={len(d)}; error = {d.Err_type[0]})", fontsize=9)
    fig.tight_layout(); fig.savefig(crop_path.replace(".png", "_QC_digitised.png"), dpi=110); plt.close(fig)


# ---------------------------------------------------------------------------------------------------------------------
# 9. FILES THAT ARE INPUTS (not regenerated by the scripts)
# ---------------------------------------------------------------------------------------------------------------------
INPUTS = [
    "data_raw/articles_PDF/",                                   # source PDFs and raw data CSVs
    "data_raw/phenotypic_rates/Crinoidea_Eucladida_rates_Wright2017_FigS6.png",      # Wright 2017 SI crops (wright2017.py
    "data_raw/phenotypic_disparity/Crinoidea_Eucladida_disparity_Wright2017_FigS8.png",  #  reads these images)
    "data_raw/taxonomic_diversity/Crinoidea_Eucladida_diversity_Wright2017_FigS10.png",
    "data_digitized/digitisation/invert_img/",                 # fan-000.png, wills-001.png, foote-000.png (pdfimages), foote_A_full.png
    "data_digitized/digitisation/close2019_work/f5/",          # img-000.jpg (pdfimages of close_2019 p. 4), qc_tops_zoom.png
    "data_digitized/digitisation/age_cutoffs.csv",             # PI age cutoffs
]


# ---------------------------------------------------------------------------------------------------------------------
# 10. REGISTRY: every script, in run order
#     method: vector | raster | raw (raw/converted data files) | derived (from other series) | utility | qc
# ---------------------------------------------------------------------------------------------------------------------
def _r(script, group, method, reference, figure, outputs, note=""):
    return dict(script=script, group=group, method=method, reference=reference, figure=figure, outputs=outputs, note=note)

REGISTRY = [
 _r("extract_raster_figures", "utility", "utility", "foth_2021, sookias_2026, quintero_2025_SI, rawson_2022, crouch_tobias_2022, felice_goswami_2017, close_2015, goswami_2022, caron_pie_2024, allvarez-carretero_2022",
    "embedded rasters", [], "Run first: extracts the raster figures of the vertebrate sources to data_digitized/digitisation/vert_img/"),
 _r("close2019_1_extract_localities", "close2019", "vector", "close_2019", "Fig. 1a", [], "locality richness points -> digitisation/close2019_fig1a_localities.csv"),
 _r("close2019_2_stage_richness", "close2019", "derived", "close_2019", "Fig. 1a", ["Tetrapoda_nonflying_taxonomic_diversity.csv"], "summed richness per ICS stage"),
 _r("close2019_3_max_richness", "close2019", "derived", "close_2019", "Fig. 1a", ["Tetrapoda_nonflying_stage_max_local_richness.csv"], "maximum local richness per stage + bootstrap CI"),
 _r("close2019_4_squamata_extract", "close2019", "raster", "close_2019", "Fig. 5c", [], "squamate dot envelope (green mask)"),
 _r("close2019_5_squamata_stages", "close2019", "derived", "close_2019", "Fig. 5c", ["Squamata_stage_max_local_richness.csv", "Squamata_max_local_richness.csv"], "max local richness of squamates per stage"),
 _r("close2019_6_squamata_plot", "close2019", "qc", "close_2019", "Fig. 5c", [], "summary plots + QC panel"),
 _r("ghezelayagh2022_fig2", "vertebrates", "vector", "ghezelayagh_2022", "Fig. 2", ["Acanthomorpha_phenotypic_disparity.csv"]),
 _r("ghezelayagh2022_edfig5e", "vertebrates", "vector", "ghezelayagh_2022", "Extended Data Fig. 5e", ["Acanthomorpha_taxonomic_diversity.csv"], "ED Fig. 5e used instead of 5b (PI)"),
 _r("brownstein_griffin2026", "vertebrates", "vector", "brownstein_griffin_2026", "Fig. 2A, 2B", ["Theropoda_non-avian_phenotypic_rate.csv", "Dinosauria_early_phenotypic_rate.csv"]),
 _r("yu2021", "vertebrates", "vector", "yu_2021_SI", "Fig. S40", ["Aves_stem_taxonomic_diversity.csv", "Aves_crown_taxonomic_diversity.csv"], "split at 100 Ma (PI)"),
 _r("navalon2022", "vertebrates", "vector", "navalon_2022", "Fig. 1a", ["Aves_crown_phenotypic_disparity.csv"]),
 _r("quintero2024", "vertebrates", "vector", "quintero_2024", "Fig. 3A", ["Mammalia_early_taxonomic_diversity_BDD.csv"]),
 _r("alvarez_carretero2022", "vertebrates", "raster", "allvarez-carretero_2022", "Fig. 3b", ["Placentalia_taxonomic_diversity.csv"]),
 _r("foth2021", "vertebrates", "raster", "foth_2021", "Fig. 9A rows 1-3", ["Archosauromorpha_stem_phenotypic_disparity.csv", "Dinosauria_early_phenotypic_disparity.csv", "Theropoda_non-avian_phenotypic_disparity.csv"]),
 _r("sookias2026", "vertebrates", "raster", "sookias_2026", "Fig. 1b, 1f, 1h", ["Archosauromorpha_stem_taxonomic_diversity.csv", "Dinosauria_early_taxonomic_diversity_PyRate.csv", "Theropoda_non-avian_taxonomic_diversity.csv"]),
 _r("quintero2025", "vertebrates", "raster", "quintero_2025_SI", "Suppl. Fig. 9", ["Dinosauria_early_taxonomic_diversity_FBDD.csv"]),
 _r("rawson2022", "vertebrates", "raster", "rawson_2022", "Fig. 5A", ["Tetrapoda_stem_phenotypic_disparity.csv"]),
 _r("crouch_tobias2022", "vertebrates", "raster", "crouch_tobias_2022", "Fig. 3A", ["Aves_crown_phenotypic_rate_continuous_data_LMM.csv"]),
 _r("felice_goswami2017", "vertebrates", "raster", "felice_goswami_2017", "Fig. 1C", ["Aves_crown_phenotypic_rate_continuous_data_GMM.csv"]),
 _r("close2015", "vertebrates", "raster", "close_2015", "Fig. 1B-D", ["Mammalia_early_phenotypic_rate.csv", "Mammalia_early_phenotypic_disparity.csv", "Mammalia_early_taxonomic_diversity.csv"]),
 _r("goswami2022", "vertebrates", "raster", "goswami_2022", "Fig. 3A", ["Placentalia_phenotypic_rate.csv"], "6 root-age lines pooled"),
 _r("raw_data_convert", "vertebrates", "raw", "simoes_pierce_2021, simoes_2022, simoes_2026, zhang_2022, simoes_2020, close_2019 (converted)", "raw data",
    ["Tetrapoda_stem_phenotypic_rate.csv", "Archosauromorpha_stem_phenotypic_rate.csv", "Lepidosauria_phenotypic_rate_continuous_data_GMM.csv", "Aves_stem_phenotypic_rate.csv",
     "Lepidosauria_phenotypic_rate_discrete_data.csv", "Lepidosauria_phenotypic_disparity.csv", "Tetrapoda_stem_taxonomic_diversity.csv", "Lepidosauria_taxonomic_diversity.csv"],
    "needs the Close 2019 outputs (run after close2019_*)"),
 _r("zhou2023", "invertebrates", "vector", "zhou_et_al_2023", "Fig. 1a", ["Bivalvia_phenotypic_disparity.csv", "Bivalvia_taxonomic_diversity.csv"]),
 _r("paterson2019", "invertebrates", "vector", "paterson_et_al_2019_withSI", "Fig. 3", ["Trilobita_phenotypic_rate.csv", "Trilobita_phenotypic_rate_rawpoints.csv"]),
 _r("novack_gottshall2022", "invertebrates", "vector", "novack-gottshall_et_al_2022", "Fig. 1b, 1c", ["Echinodermata_phenotypic_rate.csv", "Echinodermata_phenotypic_disparity.csv"]),
 _r("hopkins_smith2015", "invertebrates", "vector", "hopkins_smith_2015", "Fig. 1, 2A", ["Echinoidea_Irregularia_phenotypic_rate.csv", "Echinoidea_Irregularia_phenotypic_rate_Fig1_allEchinoidea.csv"], "Fig. 2A (Irregularia) used; Fig. 1 written for reference"),
 _r("boivin2018", "invertebrates", "vector", "boivin-et_al_2018", "Fig. 3", ["Echinoidea_Irregularia_phenotypic_disparity.csv", "Echinoidea_Irregularia_taxonomic_diversity.csv"]),
 _r("fan2020", "invertebrates", "raster", "fan_et_al_2020_SI", "Fig. S10", ["Trilobita_taxonomic_diversity.csv", "Pancrustacea_taxonomic_diversity.csv"]),
 _r("wills1998", "invertebrates", "raster", "wills_1998", "Fig. 7a", ["Pancrustacea_phenotypic_disparity.csv"]),
 _r("foote1993", "invertebrates", "raster", "foote_1993", "Fig. 5A", ["Trilobita_phenotypic_disparity.csv"]),
 _r("wright2017", "invertebrates", "raster", "Wright_2017_SI", "Figs S6, S8, S10", ["Crinoidea_Eucladida_phenotypic_rate.csv", "Crinoidea_Eucladida_phenotypic_disparity.csv", "Crinoidea_Eucladida_taxonomic_diversity.csv"],
    "reads the source crops in data_raw/; writes calibration_ticks.csv and metaDigitise_style_output.csv"),
 _r("guo2024", "mass_extinction", "vector", "guo_2024", "Fig. 2a-d", ["Brachiopoda_" + c + "_" + p + ".csv" for c in ("Terebratuliida", "Rhynchonellida", "Spiniferida", "Athyridida") for p in ("phenotypic_disparity", "taxonomic_diversity")]),
 _r("whalen2020", "mass_extinction", "raster", "whalen_2020", "Fig. 2A, 2E", ["Ammonoidea_taxonomic_diversity.csv", "Ammonoidea_phenotypic_disparity.csv"]),
 _r("cole_hopkins2021", "mass_extinction", "raster", "cole_hopkins_2021", "Fig. 1, 2B", ["Crinoidea_Diplobathrida_taxonomic_diversity.csv", "Crinoidea_Diplobathrida_phenotypic_disparity.csv"]),
 _r("lamsdell_selden2016", "mass_extinction", "raster", "lamsdell_selden_2016", "Fig. 2, 4", ["Eurypterida_taxonomic_diversity.csv", "Eurypterida_phenotypic_disparity.csv"]),
 _r("halliday_goswami2016", "mass_extinction", "raster", "halliday_goswami_2016", "Fig. 4", ["Placentalia_phenotypic_disparity.csv"]),
 _r("raw_data_squamata_crown", "mass_extinction", "raw", "simoes_2026, simoes_2020, close_2019 (converted)", "raw data",
    ["Squamata_crown_phenotypic_rate_continuous_data_GMM.csv", "Squamata_crown_phenotypic_rate_discrete_data.csv", "Squamata_crown_phenotypic_disparity.csv", "Squamata_crown_taxonomic_diversity.csv"]),
 _r("qc_panels_vert", "qc", "qc", "", "", [], "source crops + QC panels of the vertebrate figure series"),
 _r("qc_panels_invert", "qc", "qc", "", "", [], "source crops + QC panels of the invertebrate figure series"),
 _r("qc_panels_massext", "qc", "qc", "", "", [], "source crops + QC panels of the mass-extinction figure series"),
 _r("age_cutoffs", "cutoff", "utility", "", "", [], "PI age cutoffs (age_cutoffs.csv); ALWAYS LAST"),
]
# one-line description of how each script extracts its series (shown by --list --long and in each script's header)
HOW = {
 "extract_raster_figures": "saves the embedded raster images of the vertebrate source figures (page + image index) as PNGs for the tracing scripts",
 "close2019_1_extract_localities": "reads the locality points of Fig. 1a from the vector PDF and calibrates age and richness on the axis ticks",
 "close2019_2_stage_richness": "summed locality richness per ICS stage (ICS chart v2023/09)",
 "close2019_3_max_richness": "maximum local species richness per stage with bootstrap 95% CI (1000 replicates) + summary plot",
 "close2019_4_squamata_extract": "isolates the green squamate dots of Fig. 5c (raster) with a colour mask and records their upper envelope",
 "close2019_5_squamata_stages": "highest dot top per ICS stage = maximum local richness of squamates (+/-1 species reading error)",
 "close2019_6_squamata_plot": "summary plot of squamates vs all non-flying tetrapods; source crop + QC panel",
 "ghezelayagh2022_fig2": "blue DTT line read from its vector path (colour + width), averaged to a 1-Myr grid",
 "ghezelayagh2022_edfig5e": "posterior-mean net-diversification line and 95% credible band read from the vector paths",
 "brownstein_griffin2026": "coral Theropoda line (2A) and black pooled line with grey 95% band (2B) read from vector paths",
 "yu2021": "net-diversification rate bars (vector) per time bin with SE; split at 100 Ma into Aves_stem / Aves_crown",
 "navalon2022": "grey empirical disparity line read from its vector path",
 "quintero2024": "Theria lineages-through-time median + 50% CI from vector paths (log10; CI expanded to a 95% equivalent)",
 "alvarez_carretero2022": "orange eutherian-genera line traced with a colour mask in the raster inset",
 "foth2021": "coloured disparity lines of three panels traced with colour masks (largest component, per-column centre)",
 "sookias2026": "black lineage step lines of three sub-panels traced with a darkness mask (dashed zero line removed; relative units)",
 "quintero2025": "FBDD median diversity and 95% CI traced in the raster",
 "rawson2022": "black skull DTT line traced with a darkness mask; short hidden segments interpolated",
 "crouch_tobias2022": "mean log-rate line and 95% band traced in the raster (grey levels)",
 "felice_goswami2017": "black whole-skull rate line traced; segments hidden beneath the module lines interpolated",
 "close2015": "rate line (1B), disparity points + 95% CI (1C) and lineage step line (1D) traced; ages on GTS2012 epoch boundaries",
 "goswami2022": "six root-age rate lines traced by colour and pooled per 1-Myr bin (median; min-max as 95% interval)",
 "raw_data_convert": "node-based rate files binned with node_bins(); Lepidosauria disparity per time slice; Close 2019 tables reformatted",
 "zhou2023": "blue (functional groups) and black (genera) circle markers read from the vector PDF; Poisson error",
 "paterson2019": "blue per-branch rate markers (vector) binned in 2-Myr bins (mean +/- SE, n >= 3)",
 "novack_gottshall2022": "red dashed morphology lines and their bands read from vector paths (line_and_band)",
 "hopkins_smith2015": "black mean markers and vertical CI bars read from vector paths",
 "boivin2018": "top-most envelope of the stacked areas (sum of all clades) read from vector paths",
 "fan2020": "species-level curves traced per pixel column (median of line pixels), averaged in 1-Myr bins",
 "wills1998": "light-grey upper envelope traced in the raster; ages calibrated on Harland et al. 1990 period boundaries",
 "foote1993": "total-variance points and SE bars traced in the raster; time axis from the Fig. 5B ticks",
 "wright2017": "markers and error-bar ends detected in the source crops (morphological opening, peak detection); calibration on all ticks",
 "guo2024": "median SOV line, 95% envelope and dashed genera line read from vector paths at the time-bin midpoints",
 "whalen2020": "blue richness line and black MPD line traced per pixel column (colour masks), averaged in 1-Myr bins",
 "cole_hopkins2021": "navy diversity line and black SOV line + blue CI band read at the GTS2012 stage midpoints",
 "lamsdell_selden2016": "top edge of the black TDE area at stage-bin midpoints; SOV squares + error bars found by morphological opening",
 "halliday_goswami2016": "red/black markers of six MPT sets isolated by morphological opening and pooled per bin (median, min-max)",
 "raw_data_squamata_crown": "Squamata rates binned with node_bins(); disparity and diversity reformatted from the same files as Lepidosauria",
 "qc_panels_vert": "source crops + QC panels (red line, grey band) in data_raw/<parameter folder>/",
 "qc_panels_invert": "source crops + QC panels (qc_panel) in data_raw/<parameter folder>/",
 "qc_panels_massext": "source crops + QC panels (qc_panel) in data_raw/<parameter folder>/",
 "age_cutoffs": "drops rows older than the PI cutoffs (age_cutoffs.csv); rewrites a file only when its content changes",
}
for _e in REGISTRY:
    _e["how"] = HOW.get(_e["script"], "")
    _e["folder"] = "utilities" if _e["group"] in UTIL_GROUPS else "source_papers"     # sub-folder of R/digitisation/
    _e["path"] = f"R/digitisation/{_e['folder']}/{_e['script']}.py"                   # path from the Runs root
RUN_ORDER = [r["script"] for r in REGISTRY]
_BY_NAME = {r["script"]: r for r in REGISTRY}


# ---------------------------------------------------------------------------------------------------------------------
# 11. RUNNER (command line)
# ---------------------------------------------------------------------------------------------------------------------
def _select(names):
    if not names or names == ["all"]: return RUN_ORDER
    out = []
    for n in names:
        n = n[:-3] if n.endswith(".py") else n
        if n.startswith("group:"):
            out += [r["script"] for r in REGISTRY if r["group"] == n[6:]]
        elif n in RUN_ORDER: out.append(n)
        else: sys.exit(f"unknown script '{n}'; see --list")
    return [s for s in RUN_ORDER if s in out]          # always in registry order

def run_scripts(names, root=RUNS, quiet=False):
    """Run the selected scripts (registry order) as separate processes with the working directory `root`."""
    failed = []
    for s in names:
        path = os.path.join(UTILS if _BY_NAME[s]["group"] in UTIL_GROUPS else SCRIPTS, s + ".py")
        p = subprocess.run([sys.executable, path], cwd=root, capture_output=True, text=True)
        status = "ok" if p.returncode == 0 else "FAILED"
        print(f"{status:6s} {s}")
        if p.returncode != 0: failed.append(s); print(p.stderr[-2000:])
        elif not quiet and p.stdout.strip(): print("       " + "\n       ".join(p.stdout.strip().splitlines()[-3:]))
    return failed

def check(keep=False):
    """Reproducibility check: copy the INPUTS into a scratch Runs folder, run every script there in registry order,
    and compare each series in data_digitized/ (and the digitisation/ tables) with the current files (byte for byte;
    if bytes differ, the values are compared too)."""
    tmp = tempfile.mkdtemp(prefix="digitise_check_"); root = os.path.join(tmp, "Runs")
    os.makedirs(os.path.join(root, "data_digitized", "digitisation"), exist_ok=True)
    for p in INPUTS:
        src = os.path.join(RUNS, p); dst = os.path.join(root, p)
        if p.endswith("/"): shutil.copytree(src, dst, dirs_exist_ok=True)
        else: os.makedirs(os.path.dirname(dst), exist_ok=True); shutil.copy2(src, dst)
    print("scratch copy:", root)
    failed = run_scripts(RUN_ORDER, root, quiet=True)
    cur = os.path.join(RUNS, "data_digitized"); new = os.path.join(root, "data_digitized")
    rows = []
    for sub in ("", "digitisation"):
        for f in sorted(os.listdir(os.path.join(new, sub))):
            if not f.endswith(".csv"): continue
            a, b = os.path.join(cur, sub, f), os.path.join(new, sub, f)
            if not os.path.exists(a): rows.append((os.path.join(sub, f), "new (not in Runs)")); continue
            if filecmp.cmp(a, b, shallow=False): rows.append((os.path.join(sub, f), "identical")); continue
            try:
                x, y = pd.read_csv(a), pd.read_csv(b)
                same = list(x.columns) == list(y.columns) and x.shape == y.shape and x.equals(y)
                sub_cols = set(x.columns) < set(y.columns) and x.equals(y[list(x.columns)])
                rows.append((os.path.join(sub, f), "same values (format differs)" if same else
                             "same values in the current columns (current file has fewer columns)" if sub_cols else "DIFFERENT"))
            except Exception as e: rows.append((os.path.join(sub, f), f"DIFFERENT ({e})"))
    for f in sorted(x for x in os.listdir(cur) if x.endswith(".csv")):
        if not os.path.exists(os.path.join(new, f)): rows.append((f, "not produced by the scripts"))
    R = pd.DataFrame(rows, columns=["file", "status"])
    print(R.to_string(index=False)); print(R.status.value_counts().to_string())
    if failed: print("FAILED scripts:", failed)
    if not keep: shutil.rmtree(tmp, ignore_errors=True)
    return R

def main():
    ap = argparse.ArgumentParser(description="Figure-digitisation protocol: list, run and check the per-source scripts.")
    ap.add_argument("--list", action="store_true", help="table of all scripts in run order")
    ap.add_argument("--long", action="store_true", help="with --list: also print how each script extracts its series")
    ap.add_argument("--which", metavar="SERIES", help="which script writes data_digitized/<SERIES>(.csv)")
    ap.add_argument("--run", nargs="+", metavar="SCRIPT", help="'all', script names, or group:<close2019|vertebrates|invertebrates|mass_extinction|qc|cutoff>")
    ap.add_argument("--check", action="store_true", help="re-run everything in a scratch copy and compare with data_digitized/")
    ap.add_argument("--keep", action="store_true", help="with --check: keep the scratch copy")
    a = ap.parse_args()
    if a.list:
        print(f"{'step':>4}  {'script':<32} {'group':<16} {'method':<8} {'source (key reference)':<42} figure")
        for k, r in enumerate(REGISTRY, 1):
            ref = r["reference"] if len(r["reference"]) <= 42 else r["reference"][:39] + "..."
            print(f"{k:>4}  {r['script']:<32} {r['group']:<16} {r['method']:<8} {ref:<42} {r['figure']}")
            if a.long:
                print(f"{'':6}how    : {r['how']}")
                if r["outputs"]: print(f"{'':6}outputs: " + ", ".join(r["outputs"]))
                if r["note"]: print(f"{'':6}note   : {r['note']}")
        print("\nscripts of the groups utility, qc and cutoff are in R/digitisation/utilities/; all others in R/digitisation/source_papers/")
    if a.which:
        f = a.which if a.which.endswith(".csv") else a.which + ".csv"
        hits = [r for r in REGISTRY if os.path.basename(f) in r["outputs"]]
        print("\n".join(f"{r['path']}  ({r['reference']}, {r['figure']}; {r['method']})" for r in hits) or "not found in the registry")
    if a.run:
        os.chdir(RUNS); failed = run_scripts(_select(a.run))
        if failed: sys.exit(1)
    if a.check:
        os.chdir(RUNS); R = check(a.keep); sys.exit(1 if (R.status.str.startswith("DIFFERENT")).any() else 0)
    if not (a.list or a.which or a.run or a.check): ap.print_help()

if __name__ == "__main__":
    main()
