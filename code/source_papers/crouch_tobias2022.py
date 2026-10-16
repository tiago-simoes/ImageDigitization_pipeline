"""crouch_tobias2022.py -- step 19 of the digitisation pipeline: crouch_tobias_2022, Fig. 3A [raster]
Method  : mean log-rate line and 95% band traced in the raster (grey levels)
Outputs : outputs/data_tables/clade_origin/Aves_crown_phenotypic_rate_continuous_data_LMM.csv
Run     : python code/source_papers/crouch_tobias2022.py (from any folder), or python code/digitise_master.py --run crouch_tobias2022
Helpers : shared functions imported from code/digitise_master.py (documented there)
"""
import os, sys; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))   # code/, where digitise_master.py lives
from digitise_master import *                                   # shared helpers (see digitise_master.py)
from digitise_master import save_vert as save                     # vertebrate writer (rounded to 6 decimals)
from PIL import Image
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
IMG=IMG_VERT
FIG=QC_VERT
a=np.asarray(Image.open(IMG+"crouch_fig3.png").convert("L")).astype(int)
X,sx=lin([127,334,540.5,747,953.5],[80,60,40,20,0]); Y,sy=lin([90.5,311,531,751],[15,10,5,0])
x0,x1,y0,y1=107,991,55,896
m=np.zeros(a.shape,bool); m[y0:y1,x0:x1]=a[y0:y1,x0:x1]<55
lab,n=ndi.label(m,structure=np.ones((3,3))); sz=ndi.sum(m,lab,range(1,n+1))
m=np.isin(lab,[i+1 for i in range(n) if sz[i]>=25])
band=np.zeros(a.shape,bool); band[y0:y1,x0:x1]=a[y0:y1,x0:x1]<160
rows=[]
for c in range(x0,x1):
    ym=np.where(m[:,c])[0]; yb=np.where(band[:,c])[0]
    if len(ym)==0: continue
    rr=runs(ym); big=max(rr,key=lambda q:q[2]-q[1])
    # band: contiguous run containing the line
    br=[q for q in runs(yb) if q[1]<=big[0]<=q[2]]; br=br[0] if br else big
    rows.append((c,big[0],br[1],br[2]))
t=np.array(rows,float)
d=pd.DataFrame(dict(a=X(t[:,0]),v=Y(t[:,1]),hi=Y(t[:,2]),lo=Y(t[:,3]))); d["k"]=np.round(d.a)
s=d.groupby("k").agg(Age_Ma=("a","mean"),Value=("v","mean"),Err_lo=("lo","mean"),Err_hi=("hi","mean")).reset_index(drop=True)
s["Err_lo"]=np.minimum(s.Err_lo,s.Value); s["Err_hi"]=np.maximum(s.Err_hi,s.Value); s["Value_err"]=(s.Err_hi-s.Err_lo)/3.92
src="crouch_tobias_2022 Fig. 3a (all morphological traits, n=10): black solid line = mean log rate of morphological evolution of all branches (BayesTraits variable-rates; continuous data, LMM); grey band = 95% CI - raster, traced"
asd=float(np.hypot(sx,0.5/10.33)); print("cal",sx,sy)
save(s,"Aves_crown","phenotypic_rate","95CI",src,asd,"continuous_data_LMM")
fig,ax=plt.subplots(figsize=(8,7)); ax.imshow(a[:1010,:1010],cmap="gray",extent=(0,1010,1010,0))
ax.plot(t[:,0],t[:,1],'r-',lw=.6); ax.plot(t[:,0],t[:,2],'c-',lw=.4); ax.plot(t[:,0],t[:,3],'c-',lw=.4)
ax.set_title("Aves_crown rates LMM (Crouch & Tobias 2022 Fig. 3a) - QC"); fig.savefig(FIG+"Aves_crown_phenotypic_rate_continuous_data_LMM_QC_digitised.png",dpi=110,bbox_inches="tight")
