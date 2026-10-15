"""close2015.py -- step 21 of the digitisation protocol: close_2015, Fig. 1B-D [raster]
Method  : rate line (1B), disparity points + 95% CI (1C) and lineage step line (1D) traced; ages on GTS2012 epoch boundaries
Outputs : data_digitized/Mammalia_early_phenotypic_rate.csv, data_digitized/Mammalia_early_phenotypic_disparity.csv, data_digitized/Mammalia_early_taxonomic_diversity.csv
Run     : python3 R/digitisation/source_papers/close2015.py from the Runs root, or python3 R/digitisation/digitise_master.py --run close2015
Helpers : shared functions imported from R/digitisation/digitise_master.py (documented there)
(formerly R/digitise_vert/close15.py)
"""
import os, sys; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))   # R/digitisation/, where digitise_master.py lives
from digitise_master import *                                   # shared helpers (see digitise_master.py)
from digitise_master import save_vert as save                     # output writer, vertebrate dialect (rounded to 6 decimals)
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
IMG="data_digitized/digitisation/vert_img/"
FIG="data_digitized/digitisation/vert_qc/"; os.makedirs(FIG,exist_ok=True)
a=load(IMG+"close15_fig1.png"); g=a.max(2); sat=a.max(2)-a.min(2)
# time axis: epoch boundaries of the timescale bar (GTS2012)
bx=[43.5,104,144,225.5,510.5,728,813,960.5,1316,1591.5,1671]; ba=[259.8,252.17,247.2,237.0,201.3,174.1,163.5,145.0,100.5,66.0,56.0]
X,sx=lin(bx,ba); PX=1/abs(np.polyfit(bx,ba,1)[0]); print("time cal rmse",round(sx,3),"px/Myr",round(PX,2))
asd=float(np.hypot(sx,0.5/PX))
blk=(g<50)&(sat<30)
def qc(name,img_box,curves):
    x0,y0,x1,y1=img_box; fig,ax=plt.subplots(figsize=(10,4)); ax.imshow(a[y0:y1,x0:x1],extent=(x0,x1,y1,y0))
    for c,r,st in curves: ax.plot(c,r,st,lw=.8,ms=4)
    ax.set_title(name+" - QC"); fig.savefig(FIG+name+"_QC_digitised.png",dpi=110,bbox_inches="tight"); plt.close(fig)
# ---- B: mean rate line
Y,sy=lin([981.5,1302.5,1623.5,1944.5,2265],[8,6,4,2,0])
m=np.zeros(g.shape,bool); m[968:2290,86:1680]=blk[968:2290,86:1680]
lab,n=ndi.label(m,structure=np.ones((3,3))); sz=ndi.sum(m,lab,range(1,n+1)); m=lab==(np.argmax(sz)+1)
t=trace(m,minpix=1); t=t[2:-2]; c,r=interp_gaps(t[:,0],t[:,1])
s=to_series(X(c),Y(r),1.0)
save(s,"Mammalia_early","phenotypic_rate","none","close_2015 Fig. 1B: black line = mean rate of phenotypic evolution (character changes per lineage Myr) across 100 posterior trees, age-level bins (raster, traced); spread of the 100 grey lines not used as error",asd)
qc("Mammalia_early_phenotypic_rate",(40,950,1700,2300),[(c,r,'y-')])
# ---- C: disparity points + 95% CI bars
Yc,syc=lin([2082,2265.5,2449.5],[1.0,0.8,0.6])
reg=np.zeros(g.shape,bool); reg[2040:2340,280:1500]=g[2040:2340,280:1500]<70
yy,xx=np.mgrid[-4:5,-4:5]; op=ndi.binary_opening(reg,structure=(xx**2+yy**2)<=16)
lab,n=ndi.label(op); P=[]
for i,(cy,cx) in enumerate(ndi.center_of_mass(op,lab,range(1,n+1))):
    cx=int(round(cx)); cy=int(round(cy)); col=g[:,cx-1:cx+2].min(1)
    up=cy
    while up>2000 and col[up-1]<170: up-=1
    dn=cy
    while dn<2470 and col[dn+1]<170: dn+=1
    P.append((cx,cy,up,dn))
P=sorted(P); print("C points",P)
dc=pd.DataFrame([dict(Age_Ma=float(X(cx)),Value=float(Yc(cy)),Err_hi=float(Yc(up)),Err_lo=float(Yc(dn))) for cx,cy,up,dn in P])
dc["Value_err"]=(dc.Err_hi-dc.Err_lo)/3.92
save(dc,"Mammalia_early","phenotypic_disparity","95CI","close_2015 Fig. 1C: weighted mean pairwise disparity (WMPD) of epoch-level bins, black points with 95% CI bars (raster, traced)",asd)
qc("Mammalia_early_phenotypic_disparity",(250,2020,1550,2340),[([p[0] for p in P],[p[1] for p in P],'ro'),([p[0] for p in P],[p[2] for p in P],'c_'),([p[0] for p in P],[p[3] for p in P],'c_')])
# ---- D: lineage diversity step line
Yd,syd=lin([2318,2368.5,2419,2469.5,2520],[20,15,10,5,0])
m=np.zeros(g.shape,bool); m[2290:2512,172:1690]=blk[2290:2512,172:1690]
lab,n=ndi.label(m,structure=np.ones((3,3))); sz=ndi.sum(m,lab,range(1,n+1))
m=np.isin(lab,[i+1 for i in range(n) if sz[i]>=0.2*sz.max()])
t=trace(m,minpix=1); t=t[2:-2]; c,r=interp_gaps(t[:,0],t[:,1])
s=to_series(X(c),np.clip(Yd(r),0,None),1.0); s["Value_err"]=np.sqrt(s.Value); s["Err_lo"]=s.Value-s.Value_err; s["Err_hi"]=s.Value+s.Value_err
save(s,"Mammalia_early","taxonomic_diversity","Poisson","close_2015 Fig. 1D: phylogenetic lineage diversity (observed taxa + ghost lineages) step line (raster, traced)",asd)
qc("Mammalia_early_taxonomic_diversity",(150,2280,1700,2530),[(c,r,'y-')])
