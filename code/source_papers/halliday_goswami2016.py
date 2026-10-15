"""halliday_goswami2016.py -- step 37 of the digitisation protocol: halliday_goswami_2016, Fig. 4 [raster]
Method  : red/black markers of six MPT sets isolated by morphological opening and pooled per bin (median, min-max)
Outputs : data_digitized/Placentalia_phenotypic_disparity.csv
Run     : python3 R/digitisation/source_papers/halliday_goswami2016.py from the Runs root, or python3 R/digitisation/digitise_master.py --run halliday_goswami2016
Helpers : shared functions imported from R/digitisation/digitise_master.py (documented there)
(formerly R/digitise_massext/halliday.py)

Halliday & Goswami 2016 (Biol. J. Linn. Soc. 118:152-168) Fig. 4: mean pairwise distance (MPD) of eutherian
morphologies per time bin, for six sets of MPTs (red = discrete-continuous, black = discrete only; circles /
triangles / squares = minimum / full / Purgatorius constraint). Raster (JPEG 623 x 379 px).
Key: "all lines; pooling disparity uncertainty among all lines" -> per time bin, all plotted markers are pooled:
Value = median across lines, Err_lo / Err_hi = min / max across lines (treated as a 95% interval), as for the
Placentalia rates (Goswami 2022). Markers are isolated by a 3x3 morphological opening (removes the 1-px lines);
overlapping markers are split into their top / centre / bottom marker positions.
x calibration: the figure has no age ticks ("Time before present"; geoscale). Ages are calibrated on the epoch
boundaries of the time-scale bar (GTS2012: 145.0, 100.5, 66.0, 33.9, 23.03, 5.33 Ma); the six lines share the time-bin x positions,
so marker blobs are clustered in x (3-px tolerance) and each cluster is one time bin. y calibration: tick marks 0.6-1.2.
"""
import os, sys; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))   # R/digitisation/, where digitise_master.py lives
from digitise_master import *                                   # shared helpers (see digitise_master.py)
from digitise_master import save_massext as save, IMG_MASSEXT as IMG                     # mass-extinction dialect (5 decimals)
from PIL import Image
from scipy import ndimage as ndi
src = IMG + "halliday_fig4.jpg"
if not os.path.exists(src): open(src, "wb").write(pymupdf.open(PDF + "halliday_goswami_2016.pdf").extract_image(67)["image"])
im = np.asarray(Image.open(src).convert("RGB")).astype(int); R, G, B = im[..., 0], im[..., 1], im[..., 2]; a = im.mean(2)
X, sx = lin([57, 230, 365, 490, 532, 601], [145.0, 100.5, 66.0, 33.9, 23.03, 5.33])
Y, sy = lin([46.0, 111.0, 176.0, 241.0], [1.2, 1.0, 0.8, 0.6])
x0, x1, ytop, ybot = 57, 622, 0, 262
red = (R > 150) & (G < 110) & (B < 110) & (R - G > 70)
blk = (R < 75) & (G < 75) & (B < 75)
mk = np.zeros_like(red)
for m in (red, blk): mk |= ndi.binary_opening(m, structure=np.ones((3, 3)))
mk[ybot:, :] = False; mk[:, :x0] = False
# markers of all six lines share the time-bin x positions: blob centres clustered in x (tolerance 3 px)
lb, k = ndi.label(mk); blobs = []
for j in range(1, k + 1):
    yy, xx = np.where(lb == j)
    if len(yy) >= 6: blobs.append((xx.mean(), yy.min(), yy.max()))
blobs.sort(); clusters = []
for b in blobs:
    if clusters and b[0] - np.mean([q[0] for q in clusters[-1]]) <= 3: clusters[-1].append(b)
    else: clusters.append([b])
rows = []
for cl in clusters:
    ys = []
    for xc, t, b in cl:
        ys += [(t + b) / 2] if b - t <= 8 else [t + 3, (t + b) / 2, b - 3]   # merged markers: top / centre / bottom
    v = Y(np.array(ys)); xm = np.mean([q[0] for q in cl])
    rows.append((float(X(xm)), float(np.median(v)), float(v.min()), float(v.max()), len(ys)))
d = pd.DataFrame(rows, columns=["Age_Ma", "Value", "Err_lo", "Err_hi", "n_markers"])
d["Value_err"] = (d.Err_hi - d.Err_lo) / 3.92
print(d.round(3).to_string(index=False))
age_sd = np.sqrt(sx ** 2 + (0.5 / 3.89) ** 2)
print(f"x RMSE {sx:.3f} Myr (epoch boundaries), y RMSE {sy:.4f}")
save(d[["Age_Ma", "Value", "Value_err", "Err_lo", "Err_hi"]].round(4), "Placentalia", "phenotypic_disparity", "95CI",
     "halliday_goswami_2016 Fig. 4 (mean pairwise distance, all six MPT sets pooled per time bin: median across lines, min-max treated as 95% interval; low-resolution raster; ages calibrated on GTS2012 epoch boundaries of the time-scale bar)", age_sd)
