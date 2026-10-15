"""foote1993.py -- step 31 of the digitisation protocol: foote_1993, Fig. 5A [raster]
Method  : total-variance points and SE bars traced in the raster; time axis from the Fig. 5B ticks
Outputs : data_digitized/Trilobita_phenotypic_disparity.csv
Run     : python3 R/digitisation/source_papers/foote1993.py from the Runs root, or python3 R/digitisation/digitise_master.py --run foote1993
Helpers : shared functions imported from R/digitisation/digitise_master.py (documented there)
(formerly R/digitise_invert/foote.py)
"""
import os, sys; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))   # R/digitisation/, where digitise_master.py lives
from digitise_master import *                                   # shared helpers (see digitise_master.py)
from digitise_master import save_invert as save, OUT_INVERT as OUT, IMG_INVERT as IMG   # invertebrate dialect (4 decimals)
from PIL import Image
a=np.asarray(Image.open(IMG+'foote_A_full.png'))<128
X,sx=lin([377.0,638.0,898.0,1157.5,1417.5,1678.5],[560,497,434,371,308,245])
Y,sy=lin([142.5,321.0,498.5,675.5,850.5,1026.0],[0.071,0.0568,0.0426,0.0284,0.0142,0.0])
exp={"C1":420,"C2":529,"O1":650,"O2":790,"S":930,"D1":1044,"D2":1142,"EC":1268,"LC":1425,"P":1585}
rows=[]
for k,xc in exp.items():
    best=None
    for x in range(xc-14,xc+15):
        col=a[125:1015,x]; mx=cur=0; end=0
        for i,v in enumerate(col):
            cur=cur+1 if v else 0
            if cur>mx: mx=cur; end=i
        if best is None or mx>best[1]: best=(x,mx,125+end-mx+1,125+end)
    x,L,top,bot=best
    # stem centre: average columns within +-4 having run >= 0.8L
    top_c=top+3; bot_c=bot-3            # cap line half-thickness (~6px caps)
    rows.append((k,x,(top_c+bot_c)/2,top_c,bot_c))
d=pd.DataFrame(rows,columns=["interval","px","py","top","bot"])
d["Age_Ma"]=X(d.px).round(2); d["Value"]=Y(d.py); d["Err_hi"]=Y(d.top); d["Err_lo"]=Y(d.bot); d["Value_err"]=(d.Err_hi-d.Err_lo)/2
print("cal",sx,sy); print(d.round(4).to_string())
save(d,"Trilobita","phenotypic_disparity",sx,"SE","foote_1993 Fig. 5A (total variance, +/-1 SE bootstrap; time axis from Fig. 5B ticks; zero-anchors at 560/245 Ma excluded)")
