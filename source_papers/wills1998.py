"""wills1998.py -- step 30 of the digitisation protocol: wills_1998, Fig. 7a [raster]
Method  : light-grey upper envelope traced in the raster; ages calibrated on Harland et al. 1990 period boundaries
Outputs : data_digitized/Pancrustacea_phenotypic_disparity.csv
Run     : python3 R/digitisation/source_papers/wills1998.py from the Runs root, or python3 R/digitisation/digitise_master.py --run wills1998
Helpers : shared functions imported from R/digitisation/digitise_master.py (documented there)
(formerly R/digitise_invert/wills.py)
"""
import os, sys; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))   # R/digitisation/, where digitise_master.py lives
from digitise_master import *                                   # shared helpers (see digitise_master.py)
from digitise_master import save_invert as save, OUT_INVERT as OUT, IMG_INVERT as IMG   # invertebrate dialect (4 decimals)
from PIL import Image
a=np.asarray(Image.open(IMG+'wills-001.png')).astype(float)
print("grey levels sample", np.percentile(a[100:120,95:130],[5,50,95]), np.percentile(a[300:400,300:320],[5,50,95]))
lines=[86.4,140.6,201.6,227.6,272.6,330.5,368.4,401.4,456.6,468.3,529.0,581.6,586.6]
X,sx=lin([86.4,140.6,201.6,227.6,272.6,330.5,368.4,401.4,456.6,529.0,586.6],[570,510,439,409,363,290,245,208,146,65,0])
Y,sy=lin([439.0,228.0,17.5],[0,100,200])
rows=[]
for x in range(89,585):
    if min(abs(x-l) for l in lines)<2.0: continue
    col=a[8:437,x]
    nonw=(col<232)&(col>150)
    run=np.convolve(nonw.astype(int),np.ones(4,int),"valid")>=4
    idx=np.where(run)[0]
    if len(idx)==0: continue
    y=8+idx[0]
    # sub-pixel: edge at half-intensity between white and grey
    rows.append((x,y))
d=pd.DataFrame(rows,columns=["px","py"]); d["Age_Ma"]=X(d.px).round(2); d["Value"]=Y(d.py)
print("cal",sx,sy, len(d)); print(d.iloc[::12].round(1).T.to_string())
save(d[["Age_Ma","Value"]],"Pancrustacea","phenotypic_disparity",sx,"none","wills_1998 Fig. 7A (light grey upper envelope = observed + ghost-lineage disparity; Harland et al. 1990 timescale calibration)")
