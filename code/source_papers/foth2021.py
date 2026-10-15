"""foth2021.py -- step 15 of the digitisation protocol: foth_2021, Fig. 9A rows 1-3 [raster]
Method  : coloured disparity lines of three panels traced with colour masks (largest component, per-column centre)
Outputs : data_digitized/Archosauromorpha_stem_phenotypic_disparity.csv, data_digitized/Dinosauria_early_phenotypic_disparity.csv, data_digitized/Theropoda_non-avian_phenotypic_disparity.csv
Run     : python3 R/digitisation/source_papers/foth2021.py from the Runs root, or python3 R/digitisation/digitise_master.py --run foth2021
Helpers : shared functions imported from R/digitisation/digitise_master.py (documented there)
(formerly R/digitise_vert/foth.py)
"""
import os, sys; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))   # R/digitisation/, where digitise_master.py lives
from digitise_master import *                                   # shared helpers (see digitise_master.py)
from digitise_master import save_vert as save                     # output writer, vertebrate dialect (rounded to 6 decimals)
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
IMG="data_digitized/digitisation/vert_img/"
FIG="data_digitized/digitisation/vert_qc/"; os.makedirs(FIG,exist_ok=True)
src="foth_2021 Fig. 9A (cranial disparity, sum of variance; 'gradualistic' model; raster panel, colour-traced)"
P=[ # file, clade, line colour, plot box, x label centres->ages, y tick rows->values
 ("foth9A_r1.png","Archosauromorpha_stem","purple line (stem-Archosauria), row 1 left",(104,56,152),(85,48,590,345),
   [130.5,222.5,314.5,406.5,498.5],[245,225,205,185,165],[345.5,260.0,174.75,89.25],[0,.01,.02,.03]),
 ("foth9A_r2.png","Dinosauria_early","blue line (Dinosauria), row 2 left",(0,0,245),(85,0,590,293),
   [110.5,164.0,217.5,270.5,324.5,377.5,431.0,484.5,538.0],[235,225,215,205,195,185,175,165,155],[293.0,173.75,54.5],[0,.01,.02]),
 ("foth9A_r3.png","Theropoda_non-avian","magenta line (Theropoda), row 3 left",(176,82,148),(85,0,590,241),
   [83.5,140.0,196.5,252.5,309.5,365.5,422.0,478.5,534.5],[235,225,215,205,195,185,175,165,155],[241.5,167.5,92.75,18.25],[0,.01,.02,.03]),
]
for f,cl,lab,ref,box,xp,xa,yp,yv in P:
    a=load(IMG+f); X,sx=lin(xp,xa); Y,sy=lin(yp,yv)
    m=line_mask(a,ref,90,box); t,w=centre_trace(m)
    c,r=interp_gaps(t[:,0],t[:,1])
    ages=X(c); vals=np.clip(Y(r),0,None)
    s=to_series(ages,vals,1.0)
    s["Value_err"]=0.0; s["Err_lo"]=s.Value; s["Err_hi"]=s.Value
    pxmy=abs(xp[1]-xp[0])/abs(xa[1]-xa[0])
    asd=float(np.hypot(sx,0.5/pxmy))      # label-centre calibration RMSE + half-pixel
    print(cl,"line width px",w,"cal rmse",round(sx,3),round(sy,6),"px/Myr",round(pxmy,2))
    save(s,cl,"phenotypic_disparity","none",src+"; "+lab,asd)
    fig,ax=plt.subplots(figsize=(8,6)); ax.imshow(a); ax.plot(c,r,'k-',lw=.8); ax.plot(t[:,0],t[:,1],'y.',ms=1)
    ax.set_title(f"{cl} disparity - QC (black = extracted centreline)"); fig.savefig(FIG+f"{cl}_phenotypic_disparity_QC_digitised.png",dpi=110,bbox_inches="tight"); plt.close(fig)
