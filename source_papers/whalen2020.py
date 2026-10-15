"""whalen2020.py -- step 34 of the digitisation protocol: whalen_2020, Fig. 2A, 2E [raster]
Method  : blue richness line and black MPD line traced per pixel column (colour masks), averaged in 1-Myr bins
Outputs : data_digitized/Ammonoidea_taxonomic_diversity.csv, data_digitized/Ammonoidea_phenotypic_disparity.csv
Run     : python3 R/digitisation/source_papers/whalen2020.py from the Runs root, or python3 R/digitisation/digitise_master.py --run whalen2020
Helpers : shared functions imported from R/digitisation/digitise_master.py (documented there)
(formerly R/digitise_massext/whalen.py)

Whalen et al. 2020 (Sci. Adv. 6:eabc2365) Fig. 2, Paleozoic ammonoids, ranged-through per-million-year data.
Raster figure (embedded PNG, 1884 x 2341 px; extracted to data_digitized/digitisation/massext_img/whalen_fig2.png).
  A  richness: BLUE line, right axis = total species in the AMMON database (key: right axis, blue line)
  E  disparity: mean pairwise distance, solid BLACK line (left axis)
Calibration: x = centres of the age labels 400 ... 260 Ma (linear fit); y = label centres (A right axis 600/400/200)
and grid lines (E: 2, 1.5, 1, 0.5). Curves traced per pixel column (median of the line pixels), then averaged in
1-Myr bins centred on whole Ma (the data are per-million-year); bins without visible line pixels (hidden behind the
black line in A, or behind period boundary lines in E) are linearly interpolated from neighbouring bins.
"""
import os, sys; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))   # R/digitisation/, where digitise_master.py lives
from digitise_master import *                                   # shared helpers (see digitise_master.py)
from digitise_master import save_massext as save, IMG_MASSEXT as IMG                     # mass-extinction dialect (5 decimals)
from PIL import Image
src = IMG + "whalen_fig2.png"
if not os.path.exists(src):
    d = pymupdf.open(PDF + "whalen_2020.pdf"); open(src, "wb").write(d.extract_image(92)["image"])
im = np.asarray(Image.open(src).convert("RGB")).astype(int); R, G, B = im[..., 0], im[..., 1], im[..., 2]
gray = im.mean(2)
X, sx = lin([298.0, 491.5, 686.0, 879.0, 1073.5, 1266.5, 1461.0, 1654.5], [400, 380, 360, 340, 320, 300, 280, 260])
px_per_myr = abs(1 / (X(1) - X(0)))
age_sd = np.sqrt(sx ** 2 + (0.5 / px_per_myr) ** 2)
x0, x1 = 232, 1742                                   # inner plot frame (left/right axis lines excluded)

def trace(mask, y0, y1, drop_cols=()):
    cols = np.arange(x0, x1); ys = np.full(len(cols), np.nan)
    for i, c in enumerate(cols):
        if c in drop_cols: continue
        yy = np.where(mask[y0:y1, c])[0]
        if len(yy): ys[i] = y0 + np.median(yy)
    return cols, ys

def to_myr(cols, ys, Y):
    age = X(cols); v = Y(ys)
    d = pd.DataFrame({"a": age, "v": v}).dropna(); d["b"] = np.round(d.a)
    g = d.groupby("b").v.mean()
    full = np.arange(int(np.ceil(X(x1))), int(np.floor(X(x0))) + 1)
    s = g.reindex(full).interpolate(limit_area="inside")
    n_interp = int(g.reindex(full).isna().sum() - s.isna().sum())
    s = s.dropna()
    return pd.DataFrame({"Age_Ma": s.index.astype(float), "Value": s.values}), n_interp

# --- A: blue line, right axis (AMMON species)
YA, sya = lin([181.0, 245.0, 309.0], [600, 400, 200])
blue = (B > 150) & (R < 130) & (G > 70) & (G < 170) & (B - R > 70)
c, y = trace(blue, 183, 373)
rich, ni = to_myr(c, y, YA)
rich["Value"] = rich.Value.clip(lower=0).round(1)
rich = poisson(rich)
print(f"A richness: x RMSE {sx:.3f} Myr, y RMSE {sya:.2f}; {len(rich)} 1-Myr bins ({ni} interpolated)")
save(rich, "Ammonoidea", "taxonomic_diversity", "Poisson",
     "whalen_2020 Fig. 2A (blue line, right axis: total species in the AMMON database; ranged-through, per-million-year; raster traced per pixel column, 1-Myr bin means)", age_sd)

# --- E: black line, mean pairwise distance
YE, sye = lin([1035.5, 1153.0, 1211.5], [2.0, 1.0, 0.5])
black = gray < 70
ey0, ey1 = 1005, 1238
vcols = [cc for cc in range(x0 - 3, x1 + 3) if black[ey0:ey1, cc].sum() > 0.8 * (ey1 - ey0)]   # period boundary lines
drop = set(v + k for v in vcols for k in range(-3, 4))
c, y = trace(black, ey0, ey1, drop)
mpd, ni = to_myr(c, y, YE)
print(f"E MPD: y RMSE {sye:.4f}; boundary-line columns removed: {len(vcols)}; {len(mpd)} bins ({ni} interpolated)")
save(mpd, "Ammonoidea", "phenotypic_disparity", "none",
     "whalen_2020 Fig. 2E (solid black line: mean pairwise distance; ranged-through, per-million-year; no error band; raster traced per pixel column, 1-Myr bin means)", age_sd)
