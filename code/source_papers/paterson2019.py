"""paterson2019.py -- step 25 of the digitisation pipeline: paterson_et_al_2019_withSI, Fig. 3 [vector]
Method  : blue per-branch rate markers (vector) binned in 2-Myr bins (mean +/- SE, n >= 3)
Outputs : outputs/data_tables/clade_origin/Trilobita_phenotypic_rate.csv, outputs/data_tables/clade_origin/Trilobita_phenotypic_rate_rawpoints.csv
Run     : python code/source_papers/paterson2019.py (from any folder), or python code/digitise_master.py --run paterson2019
Helpers : shared functions imported from code/digitise_master.py (documented there)
"""
import os, sys; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))   # code/, where digitise_master.py lives
from digitise_master import *                                   # shared helpers (see digitise_master.py)
from digitise_master import save_invert as save, OUT_INVERT as OUT, IMG_INVERT as IMG   # invertebrate writer (4 decimals)
p=pymupdf.open(PDF+"paterson_et_al_2019_withSI.pdf")[3]
X,sx=lin([166.73,195.04,223.34,251.64,279.8,308.1,336.41,364.71,393.01,421.17,449.47],list(range(535,484,-5)))
Y,sy=lin([221.85,139.47,56.91],[0,.5,1])
P=[]
for g in p.get_drawings():
    r=g['rect']
    if g['type']=='f' and g['fill'] and max(abs(a-b) for a,b in zip(g['fill'],(0.31,0.5,0.75)))<0.02 and r.width<3 and r.y1<230:
        P.append(((r.x0+r.x1)/2,(r.y0+r.y1)/2))
pts=pd.DataFrame(P,columns=["px","py"]); pts["Age_Ma"]=X(pts.px); pts["Value"]=Y(pts.py)
pts=pts[(pts.px>166.5)&(pts.px<449.5)]
pts.round(4).to_csv(OUT+"/Trilobita_phenotypic_rate_rawpoints.csv",index=False)
print(len(pts),"points; age range",pts.Age_Ma.min().round(1),pts.Age_Ma.max().round(1), "cal",sx,sy)
BW=2.0
edges=np.arange(np.floor(pts.Age_Ma.min()/BW)*BW, pts.Age_Ma.max()+BW, BW)
pts["bin"]=pd.cut(pts.Age_Ma,edges)
g=pts.groupby("bin",observed=True).agg(Age_Ma=("Age_Ma","mean"),Value=("Value","mean"),sd=("Value","std"),n=("Value","size")).reset_index(drop=True)
g=g[g.n>=3].copy(); g["Value_err"]=g.sd/np.sqrt(g.n); g["Err_lo"]=g.Value-g.Value_err; g["Err_hi"]=g.Value+g.Value_err
print(g.round(3).T.to_string())
save(g,"Trilobita","phenotypic_rate",sx,"SE","paterson_et_al_2019 Fig. 3 (blue = discrete characters; 211 branch rates binned in 2-Myr bins, mean +/- SE, n>=3)")
