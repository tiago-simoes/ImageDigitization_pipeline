"""guo2024.py -- step 33 of the digitisation pipeline: guo_2024, Fig. 2a-d [vector]
Method  : median SOV line, 95% envelope and dashed genera line read from vector paths at the time-bin midpoints
Outputs : outputs/data_tables/clade_origin/Brachiopoda_Terebratuliida_phenotypic_disparity.csv, outputs/data_tables/clade_origin/Brachiopoda_Terebratuliida_taxonomic_diversity.csv, outputs/data_tables/clade_origin/Brachiopoda_Rhynchonellida_phenotypic_disparity.csv, outputs/data_tables/clade_origin/Brachiopoda_Rhynchonellida_taxonomic_diversity.csv, outputs/data_tables/clade_origin/Brachiopoda_Spiniferida_phenotypic_disparity.csv, outputs/data_tables/clade_origin/Brachiopoda_Spiniferida_taxonomic_diversity.csv, outputs/data_tables/clade_origin/Brachiopoda_Athyridida_phenotypic_disparity.csv, outputs/data_tables/clade_origin/Brachiopoda_Athyridida_taxonomic_diversity.csv
Run     : python code/source_papers/guo2024.py (from any folder), or python code/digitise_master.py --run guo2024
Helpers : shared functions imported from code/digitise_master.py (documented there)

Guo et al. 2024 (Nat Ecol Evol 8:1948) Fig. 2: disparity (SOV, median + 95% bootstrap envelope) and
diversity (number of genera, dashed grey line; vertex values rounded to integers, max deviation 0.04) of four brachiopod orders. Vector figure (PDF page 3):
values are read from the PDF path coordinates (exact), calibrated on the axis tick marks.
  a Terebratulida (green), b Rhynchonellida (blue), c Spiriferinida (red), d Athyridida (orange).
Clade names follow the spelling used in the keys (Brachiopoda_Terebratuliida, _Rhynchonellida, _Spiniferida, _Athyridida).
The diversity line vertices are the time-bin midpoints; the disparity median and envelope are read at the
same x positions (one value per bin).
"""
import os, sys; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))   # code/, where digitise_master.py lives
from digitise_master import *                                   # shared helpers (see digitise_master.py)
from digitise_master import save_massext as save, IMG_MASSEXT as IMG                     # mass-extinction-source writer (5 decimals)
pg = pymupdf.open(PDF + "guo_2024.pdf")[2]; dr = pg.get_drawings()
close = lambda c, ref: c is not None and max(abs(a - b) for a, b in zip(c, ref)) < 0.01
# x axes (tick marks): left column (a, b) and right column (c, d)
XL, sxl = lin([270.0, 236.7, 203.4, 170.1, 136.7, 103.4], [0, 50, 100, 150, 200, 250])
XR, sxr = lin([510.3, 478.2, 446.1, 414.1, 382.0, 349.9], [190, 210, 230, 250, 270, 290])
P = {  # panel: clade, colour, x calib, left-y ticks (px, value), right-y ticks, frame (x0, y0, x1, y1)
 "a": ("Brachiopoda_Terebratuliida", (0.09, 0.63, 0.55), XL, sxl, ([188.3, 165.5, 142.7, 119.8, 97.0, 74.1], [1.05, 1.10, 1.15, 1.20, 1.25, 1.30]),
       ([180.3, 149.6, 119.0, 88.3], [20, 40, 60, 80]), (70.8, 60.4, 270.0, 202.0), "Terebratulida"),
 "b": ("Brachiopoda_Rhynchonellida", (0.09, 0.45, 0.8), XL, sxl, ([342.1, 310.6, 279.1, 247.7], [2.2, 2.4, 2.6, 2.8]),
       ([337.9, 304.8, 271.8, 238.8], [20, 40, 60, 80]), (70.8, 222.5, 270.0, 364.1), "Rhynchonellida"),
 "c": ("Brachiopoda_Spiniferida", (1.0, 0.39, 0.28), XR, sxr, ([202.0, 167.9, 133.8, 99.7, 65.6], [0, .2, .4, .6, .8]),
       ([183.1, 162.4, 141.6, 120.9, 100.1, 79.3], [5, 10, 15, 20, 25, 30]), (335.7, 60.4, 534.8, 202.0), "Spiriferinida"),
 "d": ("Brachiopoda_Athyridida", (1.0, 0.5, 0.0), XR, sxr, ([364.1, 332.6, 301.2, 269.7, 238.2], [0, .1, .2, .3, .4]),
       ([336.2, 300.5, 264.7, 228.9], [5, 10, 15, 20]), (335.7, 222.5, 534.8, 364.1), "Athyridida")}
for k, (clade, col, X, sx, yl, yr, fr, order) in P.items():
    inside = lambda g: g["rect"].x0 >= fr[0] - 1 and g["rect"].x1 <= fr[2] + 1 and g["rect"].y0 >= fr[1] - 1 and g["rect"].y1 <= fr[3] + 1
    line = [g for g in dr if g["type"] == "s" and close(g.get("color"), col) and inside(g)]
    band = [g for g in dr if g["type"] == "f" and close(g.get("fill"), col) and inside(g)]
    dash = [g for g in dr if g["type"] == "s" and close(g.get("color"), (0.57, 0.58, 0.59)) and g.get("dashes") not in (None, "[] 0") and inside(g)]
    assert len(line) == 1 and len(band) == 1 and len(dash) == 1, (k, len(line), len(band), len(dash))
    YL, syl = lin(*yl); YR, syr = lin(*yr)
    dv = dedupe(path_pts(dash[0]))                       # diversity vertices = bin midpoints
    xs = np.array([p[0] for p in dv]); ydv = np.array([p[1] for p in dv])
    LS = segs(line[0]); BS = segs(band[0])
    med, lo, hi = [], [], []
    for x in xs:
        ym = ys_at(LS, x); yb = ys_at(BS, x)
        med.append(np.mean(ym) if ym else np.nan)
        hi.append(min(yb) if yb else np.nan); lo.append(max(yb) if yb else np.nan)   # smaller y = larger value
    age = X(xs)
    disp = pd.DataFrame({"Age_Ma": age.round(3), "Value": YL(med).clip(min=0), "Err_lo": YL(lo).clip(min=0), "Err_hi": YL(hi)})
    disp["Value_err"] = (disp.Err_hi - disp.Err_lo) / 3.92
    div = poisson(pd.DataFrame({"Age_Ma": age.round(3), "Value": YR(ydv).round(0)}))
    print(f"panel {k} {clade}: x RMSE {sx:.3f} Myr, yL RMSE {syl:.4f}, yR RMSE {syr:.3f}; {len(xs)} bins")
    save(disp, clade, "phenotypic_disparity", "95CI",
         f"guo_2024 Fig. 2{k} ({order}; solid line = median SOV of bootstrapped disparity, envelope = 95% CI; vector PDF, read at the time-bin midpoints)", sx)
    save(div, clade, "taxonomic_diversity", "Poisson",
         f"guo_2024 Fig. 2{k} ({order}; dashed grey line = number of genera per time bin; vector PDF vertices)", sx)
