"""fan2020.py -- step 29 of the digitisation pipeline: fan_et_al_2020_SI, Fig. S10 [raster]
Method  : species-level curves traced per pixel column (median of line pixels), averaged in 1-Myr bins
Outputs : outputs/data_tables/clade_origin/Trilobita_taxonomic_diversity.csv, outputs/data_tables/clade_origin/Pancrustacea_taxonomic_diversity.csv
Run     : python code/source_papers/fan2020.py (from any folder), or python code/digitise_master.py --run fan2020
Helpers : shared functions imported from code/digitise_master.py (documented there)
"""
import os, sys; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))   # code/, where digitise_master.py lives
from digitise_master import *                                   # shared helpers (see digitise_master.py)
from digitise_master import save_invert as save, OUT_INVERT as OUT, IMG_INVERT as IMG   # invertebrate writer (4 decimals)
from PIL import Image
a=np.asarray(Image.open(FIGIMG+'fan_et_al_2020_SI_FigS10.png').convert('L')).astype(int)   # PDF p. 17 image (supplied: see README)
xt=[98.5,215.5,332.5,449.0,566.0,682.5,799.5,916.5]; xv=[540,500,460,420,380,340,300,260]
P={"Trilobita":dict(x0=99.5,x1=973.5,y0=525.5,y1=1040,off=0,yt=[582.5,639.5,696.5,754.5,811.5,868.5,925.5,982.5],yv=[160,140,120,100,80,60,40,20],
                    excl=(800,973,526,600),lab='"Trilobite" (row 2, left)'),
   "Pancrustacea":dict(x0=1134.5,x1=2008,y0=1040,y1=1554.5,off=1035,yt=[1125.5,1211.5,1297.5,1383.5,1469.5],yv=[50,40,30,20,10],
                    excl=(1134,1325,1041,1125),lab='"Crustacea" (row 3, right)')}
for clade,c in P.items():
    X,sx=lin(np.array(xt)+c["off"],xv); Y,sy=lin(c["yt"],c["yv"])
    rows=[]
    for x in range(int(c["x0"])+4,int(c["x1"])-3):
        col=a[int(c["y0"])+4:int(c["y1"])-4,x]
        ys=np.where(col<110)[0]+int(c["y0"])+4
        if c["excl"][0]<=x<=c["excl"][1]: ys=ys[(ys<c["excl"][2])|(ys>c["excl"][3])]
        v=float(np.median(ys)) if len(ys) else c["y1"]
        rows.append((x,v,len(ys)))
    d=pd.DataFrame(rows,columns=["px","py","n"]); d["Age_Ma"]=X(d.px); d["Value"]=np.clip(Y(d.py),0,None)
    d["bin"]=np.round(d.Age_Ma)   # 1-Myr bins
    g=d.groupby("bin").agg(Age_Ma=("Age_Ma","mean"),Value=("Value","mean")).reset_index(drop=True)
    g["Value_err"]=np.sqrt(g.Value); g["Err_lo"]=g.Value-g.Value_err; g["Err_hi"]=g.Value+g.Value_err
    print(clade,"cal",round(sx,3),round(sy,3),len(g)); print(g.iloc[::10].round(1).T.to_string())
    save(g,clade,"taxonomic_diversity",sx,"Poisson","fan_et_al_2020_SI Fig. S10 "+c["lab"]+"; traced per pixel column, 1-Myr means")
