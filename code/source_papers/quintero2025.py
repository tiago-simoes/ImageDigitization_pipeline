"""quintero2025.py -- step 17 of the digitisation pipeline: quintero_2025_SI, Suppl. Fig. 9 [raster]
Method  : FBDD median diversity and 95% CI traced in the raster
Outputs : outputs/data_tables/clade_origin/Dinosauria_early_taxonomic_diversity_FBDD.csv
Run     : python code/source_papers/quintero2025.py (from any folder), or python code/digitise_master.py --run quintero2025
Helpers : shared functions imported from code/digitise_master.py (documented there)
"""
import os, sys; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))   # code/, where digitise_master.py lives
from digitise_master import *                                   # shared helpers (see digitise_master.py)
from digitise_master import save_vert as save                     # vertebrate writer (rounded to 6 decimals)
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
IMG=IMG_VERT
FIG=QC_VERT
a=load(IMG+"quintero25_S9.png")
X,sx=lin([875,975.5,1076,1177],[216,166,116,66]); Y,sy=lin([1119,1079,1039,999,959.5],[0,200,400,600,800])
x0,x1,y0,y1=796,1185,945,1123
sub=a[y0:y1,x0:x1]; r,gg,b=sub[...,0],sub[...,1],sub[...,2]
blue=(b-r>18)&(sub.max(2)<235)                       # any band / line pixel
med=(np.abs(sub-np.array([52,88,170])).sum(2)<75)     # median line (dark blue)
rows=[]
for c in range(sub.shape[1]):
    yb=np.where(blue[:,c])[0]; ym=np.where(med[:,c])[0]
    if len(yb)<1: continue
    rr=runs(ym) if len(ym) else []
    ymid=max(rr,key=lambda q:q[2]-q[1])[0] if rr else np.nan
    rows.append((c+x0,ymid+y0,yb.min()+y0,yb.max()+y0))
t=np.array(rows,float); ok=~np.isnan(t[:,1])
cm,ym=interp_gaps(t[ok,0],t[ok,1])
hi=np.interp(cm,t[:,0],t[:,2]); lo=np.interp(cm,t[:,0],t[:,3])
d=pd.DataFrame(dict(a=X(cm),v=np.clip(Y(ym),0,None),lo=np.clip(Y(lo),0,None),hi=Y(hi)))
d["k"]=np.round(d.a)
s=d.groupby("k").agg(Age_Ma=("a","mean"),Value=("v","mean"),Err_lo=("lo","mean"),Err_hi=("hi","mean")).reset_index(drop=True)
s["Err_lo"]=np.minimum(s.Err_lo,s.Value); s["Err_hi"]=np.maximum(s.Err_hi,s.Value); s["Value_err"]=(s.Err_hi-s.Err_lo)/3.92
src="quintero_2025_SI Supplementary Fig. 9 (FBDD, Dinosauria), bottom-right 'Diversity' panel: dark-blue median diversity trajectory; outer light band = 95% CI (raster, colour-traced)"
asd=float(np.hypot(sx,0.5/2.01)); print("cal",sx,sy,"age_sd",asd)
save(s,"Dinosauria_early","taxonomic_diversity","95CI",src,asd,"FBDD")
fig,ax=plt.subplots(figsize=(8,5)); ax.imshow(a[930:1140,780:1190],extent=(780,1190,1140,930))
ax.plot(cm,ym,'y-',lw=.8); ax.plot(cm,hi,'r-',lw=.5); ax.plot(cm,lo,'r-',lw=.5); ax.set_title("Dinosauria_early diversity FBDD (Quintero 2025 SI Fig. 9) - QC")
fig.savefig(FIG+"Dinosauria_early_taxonomic_diversity_FBDD_QC_digitised.png",dpi=120,bbox_inches="tight")
