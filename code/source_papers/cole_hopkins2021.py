"""cole_hopkins2021.py -- step 35 of the digitisation pipeline: cole_hopkins_2021, Fig. 1, 2B [raster]
Method  : navy diversity line and black SOV line + blue CI band read at the GTS2012 stage midpoints
Outputs : outputs/data_tables/clade_origin/Crinoidea_Diplobathrida_taxonomic_diversity.csv, outputs/data_tables/clade_origin/Crinoidea_Diplobathrida_phenotypic_disparity.csv
Run     : python code/source_papers/cole_hopkins2021.py (from any folder), or python code/digitise_master.py --run cole_hopkins2021
Helpers : shared functions imported from code/digitise_master.py (documented there)

Cole & Hopkins 2021 (Sci. Adv. 7:eabf4072), Crinoidea Diplobathrida.
  Fig. 1  range-through genus diversity (navy line; raster 1501 x 1188 px)  -> taxonomic_diversity
  Fig. 2B sum of variances, morphology (black line = mean of 1000 bootstraps; blue band = 95% CI;
          raster 2169 x 1506 px)                                           -> phenotypic_disparity
Both curves are polylines with one vertex per stage (Ordovician-Carboniferous stages, GTS2012). Vertex ages = stage
midpoints; inner stage boundaries are read from the stage boxes of each figure (calibrated on the x tick marks,
470 ... 330 Ma); the two stages clipped by the plot frame use GTS2012 ages (Floian base 477.7 Ma, Serpukhovian top
323.2 Ma). Each curve is read at the vertex x positions (median of line pixels in a 3-px window); the CI band is the
outer edge of the blue fill in that column.
"""
import os, sys; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))   # code/, where digitise_master.py lives
from digitise_master import *                                   # shared helpers (see digitise_master.py)
from digitise_master import save_massext as save, IMG_MASSEXT as IMG                     # mass-extinction-source writer (5 decimals)
from PIL import Image
pdf = pymupdf.open(PDF + "cole_hopkins_2021.pdf")
for x, nm in [(43, "cole_hopkins_fig1.png"), (77, "cole_hopkins_fig2.png")]:
    if not os.path.exists(IMG + nm): open(IMG + nm, "wb").write(pdf.extract_image(x)["image"])
STAGES = ["Floian", "Dapingian", "Darriwilian", "Sandbian", "Katian", "Hirnantian", "Llandovery", "Wenlock", "Ludlow",
          "Pridoli", "Lochkovian", "Pragian", "Emsian", "Eifelian", "Givetian", "Frasnian", "Famennian", "Tournaisian",
          "Visean", "Serpukhovian"]

def stage_mids(lines_px, X):
    b = list(X(np.array(lines_px)))
    b[0], b[-1] = 477.7, 323.2                      # clipped first/last stage: GTS2012 ages
    return np.array([(b[i] + b[i + 1]) / 2 for i in range(len(b) - 1)])

def col_median(mask, x, y0, y1, w=1):
    yy = np.where(mask[y0:y1, int(round(x)) - w:int(round(x)) + w + 1].any(1))[0]
    return y0 + np.median(yy) if len(yy) else np.nan

# ---------------- Fig. 1: diversity
im = np.asarray(Image.open(IMG + "cole_hopkins_fig1.png").convert("RGB")).astype(int); R, G, B = im[..., 0], im[..., 1], im[..., 2]
X1, sx1 = lin([134.5, 228.5, 322.0, 415.5, 509.0, 602.5, 696.5, 789.5, 883.5, 977.0, 1070.5, 1164.5, 1257.5, 1351.5, 1445.0],
              np.arange(470, 329, -10))
Y1, sy1 = lin([34.5, 116.5, 198.0, 280.0, 362.0, 444.0, 525.5, 607.5, 689.5, 771.5], np.arange(20, 1, -2))
lines1 = [82.5, 134.5, 160.5, 243.5, 293.5, 366.5, 380.0, 477.5, 533.5, 574.5, 610.5, 688.5, 718.5, 852.5, 905.0, 951.5, 1050.0, 1174.5, 1288.5, 1436.5, 1498.5]
mid1 = stage_mids(lines1, X1)
inv1 = lambda a: np.interp(a, X1(np.array([0, 1500]))[::-1], [1500, 0])
navy = (B > 90) & (R < 80) & (G < 80) & (B - R > 60)
v = []
for a in mid1:
    v.append(col_median(navy, inv1(a), 0, 850))
div = pd.DataFrame({"Stage": STAGES, "Age_Ma": mid1.round(2), "Value": np.round(Y1(np.array(v)), 2)})
print(div.to_string(index=False))
div["Value"] = div.Value.round(0)
div = poisson(div[["Age_Ma", "Value"]].dropna())
age_sd1 = np.sqrt(sx1 ** 2 + (0.5 / 9.35) ** 2)
print(f"Fig 1: x RMSE {sx1:.3f} Myr, y RMSE {sy1:.3f}")
save(div, "Crinoidea_Diplobathrida", "taxonomic_diversity", "Poisson",
     "cole_hopkins_2021 Fig. 1 (range-through genus diversity, one vertex per stage at the stage midpoint, GTS2012; raster; values rounded to integer genera)", age_sd1)

# ---------------- Fig. 2B: SOV morphology
im = np.asarray(Image.open(IMG + "cole_hopkins_fig2.png").convert("RGB")).astype(int); R, G, B = im[..., 0], im[..., 1], im[..., 2]
X2, sx2 = lin([1201.5, 1267.5, 1334.0, 1400.0, 1466.0, 1532.5, 1599.0, 1665.0, 1731.0, 1797.5, 1863.5, 1930.0, 1996.0, 2062.0, 2128.5],
              np.arange(470, 329, -10))
Y2, sy2 = lin([70.0, 159.0, 247.5, 337.0, 425.5, 515.0, 604.0], [0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.0])
lines2 = [1164.5, 1201.5, 1219.5, 1278.0, 1314.0, 1365.5, 1375.0, 1444.0, 1483.5, 1512.5, 1538.0, 1593.5, 1614.5, 1709.5, 1746.5, 1779.0, 1849.0, 1937.0, 2018.0, 2122.5, 2165.0]
mid2 = stage_mids(lines2, X2)
inv2 = lambda a: np.interp(a, X2(np.array([0, 2169]))[::-1], [2169, 0])
blk = (R < 70) & (G < 70) & (B < 70)
band = (B - R > 30) & (B > 180)
yt, yb = 30, 622                                     # panel B inner frame (rows)
rows = []
for st, a in zip(STAGES, mid2):
    x = int(round(inv2(a)))
    ym = col_median(blk, x, yt, yb)
    bb = np.where(band[yt:yb, x - 1:x + 2].any(1))[0]
    hi = yt + bb.min() if len(bb) else np.nan; lo = yt + bb.max() if len(bb) else np.nan
    rows.append((st, a, ym, hi, lo))
r = pd.DataFrame(rows, columns=["Stage", "Age_Ma", "ym", "yhi", "ylo"])
r["Value"] = Y2(r.ym); r["Err_hi"] = Y2(r.yhi); r["Err_lo"] = Y2(r.ylo)
print(r[["Stage", "Age_Ma", "Value", "Err_lo", "Err_hi"]].round(3).to_string(index=False))
d = r.dropna(subset=["Value", "Err_lo", "Err_hi"]).copy()      # stage plotted = line AND band present (drops the panel letter "B" above the Floian)
d["Err_lo"] = d.Err_lo.fillna(d.Value).clip(lower=0); d["Err_hi"] = d.Err_hi.fillna(d.Value)
d["Value_err"] = (d.Err_hi - d.Err_lo) / 3.92
age_sd2 = np.sqrt(sx2 ** 2 + (0.5 / 6.62) ** 2)
print(f"Fig 2B: x RMSE {sx2:.3f} Myr, y RMSE {sy2:.4f}")
save(d[["Age_Ma", "Value", "Value_err", "Err_lo", "Err_hi"]], "Crinoidea_Diplobathrida", "phenotypic_disparity", "95CI",
     "cole_hopkins_2021 Fig. 2B (sum of variances, morphology; black line = mean of 1000 bootstraps, blue band = 95% CI; one vertex per stage at the stage midpoint, GTS2012; raster)", age_sd2)
