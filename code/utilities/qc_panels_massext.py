"""qc_panels_massext.py -- step 41 of the digitisation pipeline [qc]
Method  : source crops + QC panels (qc_panel) in outputs/figures/<parameter folder>/
Note    : source crops + QC panels of the mass-extinction figure series
Run     : python code/utilities/qc_panels_massext.py (from any folder), or python code/digitise_master.py --run qc_panels_massext
Helpers : shared functions imported from code/digitise_master.py (documented there)

Source crops + QC panels (source | digitised series) for the sources of the mass_extinction group.
Writes outputs/figures/<parameter folder>/<Clade>_<parameter>_<Source>.png and ..._QC_digitised.png. Run after the
digitisation scripts (guo2024.py, whalen2020.py, cole_hopkins2021.py, lamsdell_selden2016.py, halliday_goswami2016.py).
"""
import os, sys; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))   # code/, where digitise_master.py lives
from digitise_master import *                                   # shared helpers (see digitise_master.py)
from digitise_master import save_massext as save, IMG_MASSEXT as IMG                     # mass-extinction-source writer (5 decimals)
from PIL import Image
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
vec = pdf_crop                # digitise_master.pdf_crop(pdf, page, clip, path, dpi=400)
ras = raster_crop             # digitise_master.raster_crop(src, box, path)
W = IMG + "whalen_fig2.png"; C1 = IMG + "cole_hopkins_fig1.png"; C2 = IMG + "cole_hopkins_fig2.png"
L2 = IMG + "lamsdell_fig2.jpg"; L4 = IMG + "lamsdell_fig4.jpg"; H4 = IMG + "halliday_fig4.jpg"
S = [("Ammonoidea", "taxonomic_diversity", "Ammonoidea_diversity_Whalen2020_Fig2A.png", lambda p: ras(W, (0, 0, 1884, 380), p)),
     ("Ammonoidea", "phenotypic_disparity", "Ammonoidea_disparity_Whalen2020_Fig2E.png", lambda p: ras(W, (0, 990, 1884, 1245), p))]
for k, cl, box in [("a", "Brachiopoda_Terebratuliida", (35, 45, 300, 210)), ("b", "Brachiopoda_Rhynchonellida", (35, 207, 300, 395)),
                   ("c", "Brachiopoda_Spiniferida", (300, 45, 565, 210)), ("d", "Brachiopoda_Athyridida", (300, 207, 565, 395))]:
    S += [(cl, "phenotypic_disparity", f"{cl}_disparity_Guo2024_Fig2{k}.png", lambda p, b=box: vec("guo_2024.pdf", 2, b, p)),
          (cl, "taxonomic_diversity", f"{cl}_diversity_Guo2024_Fig2{k}.png", lambda p, b=box: vec("guo_2024.pdf", 2, b, p))]
S += [("Crinoidea_Diplobathrida", "taxonomic_diversity", "Crinoidea_Diplobathrida_diversity_ColeHopkins2021_Fig1.png", lambda p: ras(C1, None, p)),
      ("Crinoidea_Diplobathrida", "phenotypic_disparity", "Crinoidea_Diplobathrida_disparity_ColeHopkins2021_Fig2B.png", lambda p: ras(C2, (1080, 0, 2169, 640), p)),
      ("Eurypterida", "taxonomic_diversity", "Eurypterida_diversity_LamsdellSelden2016_Fig2.png", lambda p: ras(L2, None, p)),
      ("Eurypterida", "phenotypic_disparity", "Eurypterida_disparity_LamsdellSelden2016_Fig4.png", lambda p: ras(L4, None, p)),
      ("Placentalia", "phenotypic_disparity", "Placentalia_disparity_HallidayGoswami2016_Fig4.png", lambda p: ras(H4, None, p))]
for cl, par, name, fn in S:
    F = os.path.join(FIGS, FOLD[par]); os.makedirs(F, exist_ok=True); path = os.path.join(F, name); fn(path)
    qc_panel(cl, par, path, f"{OUT}{cl}_{par}.csv", title_chars=95)   # source | digitised (digitise_master.qc_panel)
    print(os.path.relpath(path, ROOT))
