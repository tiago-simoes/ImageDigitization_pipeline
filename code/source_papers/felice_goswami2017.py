"""felice_goswami2017.py -- step 20 of the digitisation pipeline: felice_goswami_2017, Fig. 1C [raster]
Method  : black whole-skull rate line traced; segments hidden beneath the module lines interpolated
Outputs : outputs/data_tables/clade_origin/Aves_crown_phenotypic_rate_continuous_data_GMM.csv
Run     : python code/source_papers/felice_goswami2017.py (from any folder), or python code/digitise_master.py --run felice_goswami2017
Helpers : shared functions imported from code/digitise_master.py (documented there)
"""
import os, sys; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))   # code/, where digitise_master.py lives
from digitise_master import *                                   # shared helpers (see digitise_master.py)
from digitise_master import save_vert as save                     # vertebrate writer (rounded to 6 decimals)
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
IMG=IMG_VERT
FIG=QC_VERT
a=load(IMG+"felice_fig1.png"); g=a.max(2); sat=a.max(2)-a.min(2)
xt=[98,138,177.5,217.5,257,297,337,377,416.5,456,496,536,576,615.5,655]
X,sx=lin(xt,list(range(70,-1,-5))); Y,sy=lin([761,888.5,1016,1143],[4,3,2,1])
x0,x1,y0,y1=55,674,648,1220
m=np.zeros(g.shape,bool); m[y0:y1,x0:x1]=(g[y0:y1,x0:x1]<75)&(sat[y0:y1,x0:x1]<40)
m[1141:1146,:]=False; m[:,127:132]=False          # reference lines (rate = 1; K-Pg at 66 Ma)
m[660:905,460:674]=False                           # legend box
lab,n=ndi.label(m,structure=np.ones((3,3))); sz=ndi.sum(m,lab,range(1,n+1))
m=np.isin(lab,[i+1 for i in range(n) if sz[i]>=15])
t=trace(m,minpix=1); k=2; t=t[k:len(t)-k]
c,r=interp_gaps(t[:,0],t[:,1])
s=to_series(X(c),Y(r),1.0); s["Value_err"]=0.0; s["Err_lo"]=s.Value; s["Err_hi"]=s.Value
src=("felice_goswami_2017 Fig. 1C: black line 'Whole Skull' = relative evolutionary rate of whole-skull shape (3D GMM, continuous data) per 1-Myr bin; "
     "black segments hidden under coloured module lines interpolated - raster, traced")
asd=float(np.hypot(sx,0.5/7.96)); print("cal",sx,sy,"cols traced",len(t),"of",int(t[-1,0]-t[0,0]+1))
save(s,"Aves_crown","phenotypic_rate","none",src,asd,"continuous_data_GMM")
fig,ax=plt.subplots(figsize=(8,7)); ax.imshow(a[640:1240,40:679],extent=(40,679,1240,640)); ax.plot(c,r,'m-',lw=.7)
ax.plot(t[:,0],t[:,1],'y.',ms=1); ax.set_title("Aves_crown rates GMM (Felice & Goswami 2017 Fig. 1C) - QC (magenta = extracted)")
fig.savefig(FIG+"Aves_crown_phenotypic_rate_continuous_data_GMM_QC_digitised.png",dpi=110,bbox_inches="tight")
