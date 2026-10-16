"""lamsdell_selden2016.py -- step 36 of the digitisation pipeline: lamsdell_selden_2016, Fig. 2, 4 [raster]
Method  : top edge of the black TDE area at stage-bin midpoints; SOV squares + error bars found by morphological opening
Outputs : outputs/data_tables/clade_origin/Eurypterida_taxonomic_diversity.csv, outputs/data_tables/clade_origin/Eurypterida_phenotypic_disparity.csv
Run     : python code/source_papers/lamsdell_selden2016.py (from any folder), or python code/digitise_master.py --run lamsdell_selden2016
Helpers : shared functions imported from code/digitise_master.py (documented there)

Lamsdell & Selden 2016 (Evolution 71:95-110), Eurypterida. Raster figures (low-resolution JPEGs).
  Fig. 2  species-level diversity: solid black line = taxonomic diversity estimate (TDE, black area; the grey
          area is the phylogenetic diversity estimate PDE, not used). Top edge of the black area traced per column
          (contiguous dark region above the zero line) and read at the stage-bin midpoints (bin edges = the alternating
          grey/white stage bands of the figure).                                      -> taxonomic_diversity
  Fig. 4  sum of variances (black squares + 95% CI error bars, epoch-level bins).     -> phenotypic_disparity
x calibration: tick marks every 5 Myr (485 ... 255 Ma); y: axis tick marks.
"""
import os, sys; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))   # code/, where digitise_master.py lives
from digitise_master import *                                   # shared helpers (see digitise_master.py)
from digitise_master import save_massext as save, IMG_MASSEXT as IMG                     # mass-extinction-source writer (5 decimals)
from PIL import Image
from scipy import ndimage as ndi
pdf = pymupdf.open(PDF + "lamsdell_selden_2016.pdf")
for x, nm in [(15, "lamsdell_fig2.jpg"), (25, "lamsdell_fig4.jpg")]:
    if not os.path.exists(IMG + nm): open(IMG + nm, "wb").write(pdf.extract_image(x)["image"])

# ---------------- Fig. 2: TDE (black area)
a = np.asarray(Image.open(IMG + "lamsdell_fig2.jpg").convert("L")).astype(int)
tick2 = [57.0, 72.5, 87.5, 103.0, 118.5, 133.5, 149.0, 164.5, 180.0, 195.0, 210.5, 226.0, 241.5, 257.0, 272.0, 287.5, 303.0, 318.5,
         333.5, 349.0, 364.5, 380.0, 395.5, 411.0, 426.0, 441.5, 457.0, 472.5, 488.0, 503.0, 518.5, 534.0, 549.5, 564.5, 580.0, 595.5,
         611.0, 626.0, 641.5, 657.0, 672.5, 688.0, 703.5, 718.5, 734.0, 749.5, 765.0]
X2, sx2 = lin(tick2, np.arange(485, 254, -5))
Y2, sy2 = lin([110.5, 195.5, 280.5, 365.5, 450.0, 535.0, 620.0], [120, 100, 80, 60, 40, 20, 0])
base = 618                                                   # just above the zero line (row 620)
# stage bands: per column, among background pixels (> 205: band ~232, white 255; curves/areas excluded) take the majority class
sub = a[112:600, 57:772]; bg = (sub > 225)                  # band ~232, white 255 (curve/area edges excluded)
nb = bg.sum(0); gb = ((sub < 245) & bg).sum(0)
g = np.where(nb >= 25, (gb > 0.5 * nb).astype(float), np.nan)  # columns mostly covered by the curves: unknown
g = pd.Series(g).ffill().bfill().values.astype(int)
lab, n = ndi.label(g)
edges = []
for i in range(1, n + 1):
    w = np.where(lab == i)[0]; edges += [57 + w.min() - 0.5, 57 + w.max() + 0.5]
edges = sorted(set(e for e in edges if 60 < e < 769))        # inner stage-bin edges
mids = [(edges[i] + edges[i + 1]) / 2 for i in range(len(edges) - 1)]
dark = a < 72
def top_edge(x):
    vals = []
    for c in range(int(round(x)) - 1, int(round(x)) + 2):
        y = base
        while y > base - 45 and not dark[y, c]: y -= 1         # the dashed line may cover the area at the zero line
        if not dark[y, c]: vals.append(np.nan); continue
        while y > 100:                                       # climb the black area; jump over the dashed line drawn on top of it (gaps up to 45 px; only the TDE area is this dark)
            if dark[y - 1, c]: y -= 1; continue
            nxt = [k for k in range(2, 46) if dark[y - k, c]]
            if nxt: y -= nxt[0]; continue
            break
        vals.append(y)
    return np.nanmedian(vals) if not all(np.isnan(vals)) else 620.0
rows = [(float(X2(x)), float(Y2(top_edge(x)))) for x in mids]
d = pd.DataFrame(rows, columns=["Age_Ma", "Value"]).round(3)
d = d[d.Age_Ma <= 472]                                       # TDE starts in the Middle Ordovician
d["Value"] = d.Value.clip(lower=0).round(0)
print(d.T.to_string(header=False)); print('bin edges (Ma):', np.round(X2(np.array(edges)), 1).tolist())
age_sd2 = np.sqrt(sx2 ** 2 + (0.5 / 3.1) ** 2)
print(f"Fig 2: x RMSE {sx2:.3f} Myr, y RMSE {sy2:.3f}; {len(d)} stage bins")
save(poisson(d), "Eurypterida", "taxonomic_diversity", "Poisson",
     "lamsdell_selden_2016 Fig. 2 (black area = taxonomic diversity estimate TDE, species; read at stage-bin midpoints; low-resolution raster, values rounded to integers)", age_sd2)

# ---------------- Fig. 4: SOV (black squares + 95% CI)
a = np.asarray(Image.open(IMG + "lamsdell_fig4.jpg").convert("L")).astype(int)
tick4 = [64.0, 79.0, 94.0, 109.5, 124.5, 139.5, 154.5, 170.0, 185.0, 200.0, 215.5, 230.0, 245.5, 260.5, 275.5, 291.0, 306.0, 321.0,
         336.0, 351.0, 366.5, 381.5, 396.5, 412.0, 427.0, 442.0, 457.0, 472.5, 487.5, 502.5, 518.0, 533.0, 548.0, 563.0, 578.5, 593.5,
         608.5, 623.5, 639.0, 654.0, 669.0, 684.0, 699.5, 714.5, 729.5, 745.0, 759.5]
X4, sx4 = lin(tick4, np.arange(485, 254, -5))
Y4, sy4 = lin([115.0, 162.0, 209.0, 257.5, 301.0, 348.5, 393.5], [7000, 6500, 6000, 5500, 5000, 4500, 4000])
blk = a < 90
op = ndi.binary_opening(blk[110:420, 64:768], structure=np.ones((5, 5)))
lab, n = ndi.label(op)
rows = []
for i in range(1, n + 1):
    yy, xx = np.where(lab == i); cy, cx = 110 + yy.mean(), 64 + xx.mean(); top, bot = 110 + yy.min(), 110 + yy.max()
    c = int(round(cx)); up = top; dn = bot
    bar = a < 150                                            # error bar (black) may cross the grey SOR line/markers
    while up > 100 and bar[up - 1, c - 1:c + 2].any() and blk[up - 1, c - 3:c + 4].sum() <= 5: up -= 1
    while dn < 445 and bar[dn + 1, c - 1:c + 2].any() and blk[dn + 1, c - 3:c + 4].sum() <= 5: dn += 1
    rows.append((X4(cx), Y4(cy), Y4(dn), Y4(up), dn - bot, top - up))
d = pd.DataFrame(rows, columns=["Age_Ma", "Value", "Err_lo", "Err_hi", "bar_dn_px", "bar_up_px"]).sort_values("Age_Ma", ascending=False)
d["Err_lo"] = np.where(d.bar_dn_px <= 1, d.Value, d.Err_lo); d["Err_hi"] = np.where(d.bar_up_px <= 1, d.Value, d.Err_hi)   # bar hidden by the square
d["Value_err"] = (d.Err_hi - d.Err_lo) / 3.92
print(d.round(1).to_string(index=False))
age_sd4 = np.sqrt(sx4 ** 2 + (0.5 / 3.02) ** 2)
print(f"Fig 4: x RMSE {sx4:.3f} Myr, y RMSE {sy4:.2f}")
save(d[["Age_Ma", "Value", "Value_err", "Err_lo", "Err_hi"]].round(3), "Eurypterida", "phenotypic_disparity", "95CI",
     "lamsdell_selden_2016 Fig. 4 (sum of variances, black squares = epoch-level bins, error bars = 95% CI; rarefied; low-resolution raster)", age_sd4)
