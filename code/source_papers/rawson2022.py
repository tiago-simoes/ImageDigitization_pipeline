"""rawson2022.py -- step 18 of the digitisation pipeline: rawson_2022, Fig. 5A [raster]
Method  : black skull DTT line traced with a darkness mask; short hidden segments interpolated
Outputs : outputs/data_tables/clade_origin/Tetrapoda_stem_phenotypic_disparity.csv
Run     : python code/source_papers/rawson2022.py (from any folder), or python code/digitise_master.py --run rawson2022
Helpers : shared functions imported from code/digitise_master.py (documented there)
"""
import os, sys; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))   # code/, where digitise_master.py lives
from digitise_master import *                                   # shared helpers (see digitise_master.py)
from digitise_master import save_vert as save                     # vertebrate writer (rounded to 6 decimals)
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
IMG=IMG_VERT
FIG=QC_VERT
a=load(IMG+"rawson_fig5.png"); g=a.max(2)
X,sx=lin([458.5,581,703.5,826.5],[400,350,300,250]); Y,sy=lin([521,396,271.5,147],[0,1,2,3])
m=np.zeros(g.shape,bool); sub=a[90:537,388:850]; m[90:537,388:850]=(g[90:537,388:850]<75)|((sub[...,0]<160)&(sub[...,1]<45)&(sub[...,2]<45))   # black, or black under the translucent red Romer's-gap bar
lab,n=ndi.label(m,structure=np.ones((3,3))); sz=ndi.sum(m,lab,range(1,n+1))
order=np.argsort(sz)[::-1]; print("largest comps",sz[order[:6]], [ndi.find_objects(lab)[i] for i in order[:3]])
mm=np.isin(lab,[i+1 for i in order if sz[i]>=500])   # line is split by the two event stars -> keep both main pieces, gap interpolated
t=trace(mm,minpix=1); w=np.median(t[:,3]-t[:,2]+1); k=int(round(w/2)); t=t[k:len(t)-k]
c,r=interp_gaps(t[:,0],t[:,1])
s=to_series(X(c),np.clip(Y(r),0,None),1.0); s["Value_err"]=0.0; s["Err_lo"]=s.Value; s["Err_hi"]=s.Value
src="rawson_2022 Fig. 5A (Skull): black solid line = mean relative subclade disparity-through-time (topological DTT, tree A); grey band = BM simulation 95% CI (not data error) - raster, traced"
asd=float(np.hypot(sx,0.5/2.453)); print("w",w,"cal",sx,sy)
save(s,"Tetrapoda_stem","phenotypic_disparity","none",src,asd)
fig,ax=plt.subplots(figsize=(8,6)); ax.imshow(a[60:580,360:870],extent=(360,870,580,60)); ax.plot(c,r,'c-',lw=.8)
ax.set_title("Tetrapoda_stem disparity (Rawson 2022 Fig. 5A) - QC (cyan = extracted)")
fig.savefig(FIG+"Tetrapoda_stem_phenotypic_disparity_QC_digitised.png",dpi=110,bbox_inches="tight")
