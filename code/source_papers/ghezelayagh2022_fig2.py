"""ghezelayagh2022_fig2.py -- step 8 of the digitisation pipeline: ghezelayagh_2022, Fig. 2 [vector]
Method  : blue DTT line read from its vector path (colour + width), averaged to a 1-Myr grid
Outputs : outputs/data_tables/clade_origin/Acanthomorpha_phenotypic_disparity.csv
Run     : python code/source_papers/ghezelayagh2022_fig2.py (from any folder), or python code/digitise_master.py --run ghezelayagh2022_fig2
Helpers : shared functions imported from code/digitise_master.py (documented there)
"""
import os, sys; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))   # code/, where digitise_master.py lives
from digitise_master import *                                   # shared helpers (see digitise_master.py)
from digitise_master import save_vert as save                     # vertebrate writer (rounded to 6 decimals)
p=pymupdf.open(A+"ghezelayagh_2022.pdf")[3]
X,sx=lin([94.31,136.71,178.3,220.29,261.88,303.46,345.45,387.04,428.62],[160,140,120,100,80,60,40,20,0])
Y,sy=lin([501.57,460.24,418.92,377.59,336.27,294.94],[0,.2,.4,.6,.8,1.0])
for g in p.get_drawings():
    c=g.get('color')
    if g['type']=='s' and c and abs(c[0]-0.262)<0.01 and abs(c[2]-0.68)<0.01 and g['width']>1.4: L=dedupe(pts(g))
L=np.array(L); i=np.where(L[:,0]>=114.2)[0][0]; x0,y0=L[i-1]; x1,y1=L[i]; yb=y0+(y1-y0)*(114.2-x0)/(x1-x0); L=np.vstack([[114.2,yb],L[i:]]); print("blue pts",len(L),"x",L[:,0].min(),L[:,0].max(), "cal",sx,sy)
s=to_series(X(L[:,0]),Y(L[:,1]),step=1.0)
save(s,"Acanthomorpha","phenotypic_disparity","none","ghezelayagh_2022 Fig. 2 (solid blue line: observed mean relative subclade disparity, DTT; vector-extracted, averaged to 1-Myr grid; line hidden left of the crown-age band at ~150.7 Ma excluded)",sx)
