"""sookias2026.py -- step 16 of the digitisation pipeline: sookias_2026, Fig. 1b, 1f, 1h [raster]
Method  : black lineage step lines of three sub-panels traced with a darkness mask (dashed zero line removed; relative units)
Outputs : outputs/data_tables/clade_origin/Archosauromorpha_stem_taxonomic_diversity.csv, outputs/data_tables/clade_origin/Dinosauria_early_taxonomic_diversity_PyRate.csv, outputs/data_tables/clade_origin/Theropoda_non-avian_taxonomic_diversity.csv
Run     : python code/source_papers/sookias2026.py (from any folder), or python code/digitise_master.py --run sookias2026
Helpers : shared functions imported from code/digitise_master.py (documented there)
"""
import os, sys; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))   # code/, where digitise_master.py lives
from digitise_master import *                                   # shared helpers (see digitise_master.py)
from digitise_master import save_vert as save                     # vertebrate writer (rounded to 6 decimals)
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
IMG=IMG_VERT
FIG=QC_VERT
a=load(IMG+"sookias_fig1.png"); g=a.max(2)
X,sx=lin([274.5,429.25,584.0,738.5],[240,220,200,180])     # sub-panel ticks (same time axis as Fig. 1a: 240-180 Ma)
src=("sookias_2026 Fig. 1{p} left plot ({lab}): black solid line = number of lineages through time (PyRate); "
     "lineage axis unlabelled in sub-panels -> Value = relative lineage count (fraction of panel height, 0 = panel floor)")
P=[("b","Stem-Archosauria","Archosauromorpha_stem",858.5,1310.5,(1107,1112),""),
   ("f","Dinosauria","Dinosauria_early",1942.5,2394.5,(2190,2196),"PyRate"),
   ("h","Theropoda","Theropoda_non-avian",2486.0,2935.0,(2733,2738),"")]
for p,lab,cl,top,bot,dash,suf in P:
    y0,y1=int(top)+4,int(bot)-2
    m=np.zeros(g.shape,bool); m[y0:y1,148:779]=g[y0:y1,148:779]<80
    m[dash[0]:dash[1]+1,:]=False                               # drop dashed zero-rate line
    lab_,n=ndi.label(m,structure=np.ones((3,3))); sz=ndi.sum(m,lab_,range(1,n+1))
    keep=[i+1 for i in range(n) if sz[i]>=0.15*sz.max()]        # main step line (may be split where it crosses the dashed line)
    m=np.isin(lab_,keep)
    t=trace(m,minpix=1); w=np.median(t[:,3]-t[:,2]+1)
    k=int(round(w/2)); t=t[k:len(t)-k]
    c,r=interp_gaps(t[:,0],t[:,1])
    rel=(bot-r)/(bot-top)
    s=to_series(X(c),rel,1.0); s["Value_err"]=0.0; s["Err_lo"]=s.Value; s["Err_hi"]=s.Value
    asd=float(np.hypot(sx,0.5/7.73))
    print(cl,"components kept",len(keep),"width",w,"cal",round(sx,4))
    save(s,cl,"taxonomic_diversity","none",src.format(p=p,lab=lab),asd,suf)
    fig,ax=plt.subplots(figsize=(8,6)); ax.imshow(a[int(top)-10:int(bot)+40,130:800],extent=(130,800,int(bot)+40,int(top)-10))
    ax.plot(c,r,'y-',lw=.8); ax.set_title(f"{cl} diversity (Sookias 2026 Fig.1{p}) - QC (yellow = extracted)")
    fig.savefig(FIG+f"{cl}_taxonomic_diversity{('_'+suf) if suf else ''}_QC_digitised.png",dpi=110,bbox_inches="tight"); plt.close(fig)
