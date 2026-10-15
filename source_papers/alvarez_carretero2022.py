"""alvarez_carretero2022.py -- step 14 of the digitisation protocol: allvarez-carretero_2022, Fig. 3b [raster]
Method  : orange eutherian-genera line traced with a colour mask in the raster inset
Outputs : data_digitized/Placentalia_taxonomic_diversity.csv
Run     : python3 R/digitisation/source_papers/alvarez_carretero2022.py from the Runs root, or python3 R/digitisation/digitise_master.py --run alvarez_carretero2022
Helpers : shared functions imported from R/digitisation/digitise_master.py (documented there)
(formerly R/digitise_vert/alvarez.py)

Alvarez-Carretero et al. 2022 (Nature) Fig. 3b inset: orange line = number of EUTHERIAN GENERA through time
(PaleoDB), log10 y axis -> Placentalia taxonomic diversity (second, independent series). Raster trace.
"""
import os, sys; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))   # R/digitisation/, where digitise_master.py lives
from digitise_master import *                                   # shared helpers (see digitise_master.py)
from digitise_master import save_vert as save                     # output writer, vertebrate dialect (rounded to 6 decimals)
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
IMG="data_digitized/digitisation/vert_img/"; FIG="data_digitized/digitisation/vert_qc/"; os.makedirs(FIG,exist_ok=True)
a=load(IMG+"alvarez_fig3b.png"); R,G,B=a[...,0],a[...,1],a[...,2]
X,sx=lin([32.5,257.5,482.5,706,931,1155.5],[250,200,150,100,50,0])                       # tick marks inside the image
Y,sy=lin([30.5,152.5,274.5,396.5,482],np.log10([5000,500,50,5,1]))                        # log10(no. of lineages)
m=(R>200)&(G>120)&(G<210)&(B<90); m[:,:34]=False; m[481:,:]=False                         # orange; drop axis frame
lb,n=ndi.label(m,structure=np.ones((3,3))); sz=ndi.sum(m,lb,range(1,n+1)); m=np.isin(lb,[i+1 for i in range(n) if sz[i]>=6])
t=trace(m,minpix=1); t=t[1:-1]
c,r=interp_gaps(t[:,0],t[:,1])
gap=np.ones(len(c),bool); gap[np.isin(c,t[:,0])]=False                                    # stretches hidden under the blue/black lines:
lbl,ng=ndi.label(gap)                                                                      # <= 15 px (~3.3 Myr) interpolated, longer left missing
for j in range(1,ng+1):
    if (lbl==j).sum()<=15: gap[lbl==j]=False
c,r=c[~gap],r[~gap]
s=to_series(X(c),np.clip(Y(r),0,None),1.0); s["Value_err"]=0.0; s["Err_lo"]=s.Value; s["Err_hi"]=s.Value
src=("allvarez-carretero_2022 Fig. 3b: orange line = number of eutherian genera through time (PaleoDB); "
     "Value = log10(no. of genera) as plotted on the log axis; raster trace, 1-Myr means")
asd=float(np.hypot(sx,0.5/4.49)); print("cal",round(sx,3),round(sy,4),"cols",len(t))
save(s,"Placentalia","taxonomic_diversity","none",src,asd)
fig,ax=plt.subplots(figsize=(9,4)); ax.imshow(a); ax.plot(c,r,'k-',lw=.7); ax.set_title("Placentalia diversity (Alvarez-Carretero 2022 Fig. 3b, orange) - QC")
fig.savefig(FIG+"Placentalia_taxonomic_diversity_QC_trace.png",dpi=110,bbox_inches="tight")
