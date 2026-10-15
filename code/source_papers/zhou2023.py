"""zhou2023.py -- step 24 of the digitisation protocol: zhou_et_al_2023, Fig. 1a [vector]
Method  : blue (functional groups) and black (genera) circle markers read from the vector PDF; Poisson error
Outputs : data_digitized/Bivalvia_phenotypic_disparity.csv, data_digitized/Bivalvia_taxonomic_diversity.csv
Run     : python3 R/digitisation/source_papers/zhou2023.py from the Runs root, or python3 R/digitisation/digitise_master.py --run zhou2023
Helpers : shared functions imported from R/digitisation/digitise_master.py (documented there)
(formerly R/digitise_invert/zhou.py)
"""
import os, sys; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))   # R/digitisation/, where digitise_master.py lives
from digitise_master import *                                   # shared helpers (see digitise_master.py)
from digitise_master import save_invert as save, OUT_INVERT as OUT, IMG_INVERT as IMG   # invertebrate dialect (4 decimals)
p=pymupdf.open(PDF+"zhou_et_al_2023.pdf")[1]; dr=p.get_drawings()
# x calibration: age labels (centre), y: dotted gridlines 0..70
labs=[(w[0]+w[2])/2 for w in p.get_text("words") if 175<w[1]<180 and w[0]<215 and w[4].isdigit()]
X,sx=lin(sorted(labs),[520,510,500,490,480,470,460,450,440])
Y,sy=lin([148.86,136.71,124.55,112.4,100.25,88.09,75.93,63.79],[0,10,20,30,40,50,60,70])
box=(95,50,215,150.5)
for col,param in [((0,0,1),"phenotypic_disparity"),((0,0,0),"taxonomic_diversity")]:
    c=[q for q in circle_centres(dr,col,0.52,box) if q[1]>62.5 or q[0]>110]  # drop legend
    c=[q for q in c if not (q[0]<110 and q[1]<75)]
    d=pd.DataFrame({"Age_Ma":X([q[0] for q in c]).round(2),"Value":Y([q[1] for q in c]).round(2)})
    d["Value"]=d["Value"].clip(lower=0)
    d["Value_err"]=np.sqrt(d["Value"]); d["Err_lo"]=d["Value"]-d["Value_err"]; d["Err_hi"]=d["Value"]+d["Value_err"]
    print(param,"x RMSE",round(sx,3),"y RMSE",round(sy,3)); print(d.T.to_string())
    save(d,"Bivalvia",param,sx,"Poisson","zhou_et_al_2023 Fig. 1a ("+("blue: no. FGs" if col==(0,0,1) else "black: no. genera")+")")
