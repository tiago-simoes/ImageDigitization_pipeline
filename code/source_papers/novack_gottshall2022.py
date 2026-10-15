"""novack_gottshall2022.py -- step 26 of the digitisation protocol: novack-gottshall_et_al_2022, Fig. 1b, 1c [vector]
Method  : red dashed morphology lines and their bands read from vector paths (line_and_band)
Outputs : data_digitized/Echinodermata_phenotypic_rate.csv, data_digitized/Echinodermata_phenotypic_disparity.csv
Run     : python3 R/digitisation/source_papers/novack_gottshall2022.py from the Runs root, or python3 R/digitisation/digitise_master.py --run novack_gottshall2022
Helpers : shared functions imported from R/digitisation/digitise_master.py (documented there)
(formerly R/digitise_invert/novack.py)
"""
import os, sys; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))   # R/digitisation/, where digitise_master.py lives
from digitise_master import *                                   # shared helpers (see digitise_master.py)
from digitise_master import save_invert as save, OUT_INVERT as OUT, IMG_INVERT as IMG   # invertebrate dialect (4 decimals)
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
