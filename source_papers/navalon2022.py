"""navalon2022.py -- step 12 of the digitisation protocol: navalon_2022, Fig. 1a [vector]
Method  : grey empirical disparity line read from its vector path
Outputs : data_digitized/Aves_crown_phenotypic_disparity.csv
Run     : python3 R/digitisation/source_papers/navalon2022.py from the Runs root, or python3 R/digitisation/digitise_master.py --run navalon2022
Helpers : shared functions imported from R/digitisation/digitise_master.py (documented there)
(formerly R/digitise_vert/navalon.py)
"""
import os, sys; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))   # R/digitisation/, where digitise_master.py lives
from digitise_master import *                                   # shared helpers (see digitise_master.py)
from digitise_master import save_vert as save                     # output writer, vertebrate dialect (rounded to 6 decimals)
d=pymupdf.open(A+'navalon_2022.pdf'); p=d[1]; dr=p.get_drawings()
X,sx=lin([341.95+12.9614*i for i in range(8)],[70,60,50,40,30,20,10,0]); Y,sy=lin([99.26,114.81,130.37,145.92,161.47,177.03],[1,.8,.6,.4,.2,0])
G=[g for g in dr if g['type']=='f' and g.get('fill') and abs(g['fill'][0]-0.48)<0.03 and abs(g['fill'][2]-0.49)<0.03 and g['rect'].x1<440]
print(len(G),[g['rect'] for g in G])
x0=min(g['rect'].x0 for g in G)+0.35; x1=max(g['rect'].x1 for g in G)-0.35
xs=np.arange(x0,x1,0.05); pr=poly_profile(G,xs); ok=~np.isnan(pr[:,1])
s=to_series(X(pr[ok,0]),Y(pr[ok,1]),1.0)
src="navalon_2022 Fig. 1a (skeletal form): grey continuous line = empirical disparity (mean pairwise distance) time series for Neornithes from ancestral character estimation (vector, exact); blue band = BM simulation CI (not data error)"
print("cal",sx,sy)
save(s,"Aves_crown","phenotypic_disparity","none",src,sx)
