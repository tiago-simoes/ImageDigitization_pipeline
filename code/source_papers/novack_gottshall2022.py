"""novack_gottshall2022.py -- step 26 of the digitisation pipeline: novack-gottshall_et_al_2022, Fig. 1b, 1c; Extended Data Fig. 1a [vector + raster]
Method  : Fig. 1b, 1c: red dashed morphology lines and their bands read from vector paths (line_and_band);
          Extended Data Fig. 1a (raster): top edge of the stacked lineage-richness area (sum of all classes) traced per pixel column, 1-Myr means
Outputs : outputs/data_tables/clade_origin/Echinodermata_phenotypic_rate.csv, outputs/data_tables/clade_origin/Echinodermata_phenotypic_disparity.csv,
          outputs/data_tables/clade_origin/Echinodermata_taxonomic_diversity.csv
Run     : python code/source_papers/novack_gottshall2022.py (from any folder), or python code/digitise_master.py --run novack_gottshall2022
Helpers : shared functions imported from code/digitise_master.py (documented there)
"""
import os, sys; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))   # code/, where digitise_master.py lives
from digitise_master import *                                   # shared helpers (see digitise_master.py)
from digitise_master import save_invert as save, OUT_INVERT as OUT, IMG_INVERT as IMG   # invertebrate writer (4 decimals)
p=pymupdf.open(PDF+"novack-gottshall_et_al_2022.pdf")[1]; dr=p.get_drawings()
BLUE=(0.8,0.27,0.47)   # pink/red dashed = Morphology
panels={"phenotypic_rate":dict(x=([344.8,437.48,506.82],[541,485.4,443.8]),y=([236.95,192.25,147.55,102.85],[0,.01,.02,.03]),box=(344.8,60,506.9,290),fig="Fig. 1b"),
        "phenotypic_disparity":dict(x=([114.14,206.82,276.15],[541,485.4,443.8]),y=([438.65,405.4,372.15,338.85,305.55,272.35],[0,.2,.4,.6,.8,1.0]),box=(114.14,250,276.2,445),fig="Fig. 1c")}
for par,c in panels.items():
    X,sx=lin(*c["x"]); Y,sy=lin(*c["y"])
    rows=line_and_band(dr,BLUE,c["box"])
    d=pd.DataFrame(rows,columns=["px","py","up","lo"])
    d["Age_Ma"]=X(d.px).round(2); d["Value"]=Y(d.py); d["Err_hi"]=Y(d.up); d["Err_lo"]=Y(d.lo)
    d["Value_err"]=(d.Err_hi-d.Err_lo)/2
    print(par,"cal RMSE",round(sx,3),round(sy,5)); print(d[["Age_Ma","Value","Err_lo","Err_hi"]].round(4).T.to_string())
    save(d,"Echinodermata",par,sx,"SE","novack-gottshall_et_al_2022 "+c["fig"]+" (red/pink dashed: Morphology; band = SD across 50 trees / SE)")

# ---------------- Extended Data Fig. 1a: Echinodermata taxonomic diversity
# Key: "stacked plot with multiple groups: consider sum of all groups (top-most line)". Panel a = tree-based lineage
# richness (mean of 366 tip and 365 node ranges across the time-scaled trees), stacked by class; the sum of all classes
# is the top edge of the stacked area. Raster figure (JPEG 1772 x 1968 px, PDF page 11), panel a = x 1001-1721, y 41-738.
from PIL import Image
src = IMG + "novack_edfig1.jpeg"                  # outputs/intermediates/invert_img/ (extracted here if missing)
if not os.path.exists(src):
    _d = pymupdf.open(PDF + "novack-gottshall_et_al_2022.pdf")
    open(src, "wb").write(_d.extract_image(_d[10].get_images()[0][0])["image"])
im = np.asarray(Image.open(src).convert("RGB")).astype(int)
# ages: labelled ticks 541, 485.4, 443.8 Ma (tick centres x = 1001, 1412.5, 1721) + the GTS2012 boundaries of the series
# boxes of the time scale (Terreneuvian/Series 2 521, Series 2/Miaolingian 509, Miaolingian/Furongian 497, Lower/Middle
# Ordovician 470.0, Middle/Upper Ordovician 458.4 Ma); values: left-axis ticks 50-250 lineages (0 falls on the baseline)
X, sx = lin([1001, 1148.5, 1237.5, 1326.5, 1412.5, 1526, 1612.5, 1721], [541, 521, 509, 497, 485.4, 470.0, 458.4, 443.8])
Y, sy = lin([100.5, 229, 355.5, 482.5, 610], [250, 200, 150, 100, 50])
px_per_myr = abs(1 / (X(1) - X(0))); age_sd = float(np.sqrt(sx ** 2 + (0.5 / px_per_myr) ** 2))
rows = []
for x in range(1004, 1719):                       # inside the plot frame
    c = im[44:737, x].min(axis=1) < 200           # coloured (any class) vs white background
    run = c[:-1] & c[1:]                          # two coloured pixels in a row (ignores isolated JPEG noise)
    if run.any(): rows.append((x, 44 + int(np.argmax(run)) - 0.5))   # edge between the last white and first coloured pixel
    # columns without coloured pixels (white Cambrian/Ordovician boundary line at 485.4 Ma) are skipped
d = pd.DataFrame(rows, columns=["px", "py"]); d["Age_Ma"] = X(d.px); d["Value"] = Y(d.py)
d["bin"] = np.round(d.Age_Ma)                     # 1-Myr bins
g = d.groupby("bin").agg(Age_Ma=("Age_Ma", "mean"), Value=("Value", "mean")).reset_index(drop=True)
g = poisson(g)                                    # lineage counts without published error: +/- sqrt(n)
print("ED Fig. 1a cal RMSE", round(sx, 3), "Myr,", round(sy, 3), "lineages; columns at the frame top:", int((d.py < 45).sum()))
save(g, "Echinodermata", "taxonomic_diversity", age_sd, "Poisson",
     "novack-gottshall_et_al_2022 Extended Data Fig. 1a: lineage richness by class, top edge of the stacked area "
     "(sum of all classes); traced per pixel column, 1-Myr means")
