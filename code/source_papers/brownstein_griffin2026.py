"""brownstein_griffin2026.py -- step 10 of the digitisation pipeline: brownstein_griffin_2026, Fig. 2A, 2B [vector]
Method  : coral Theropoda line (2A) and black pooled line with grey 95% band (2B) read from vector paths
Outputs : outputs/data_tables/clade_origin/Theropoda_non-avian_phenotypic_rate.csv, outputs/data_tables/clade_origin/Dinosauria_early_phenotypic_rate.csv
Run     : python code/source_papers/brownstein_griffin2026.py (from any folder), or python code/digitise_master.py --run brownstein_griffin2026
Helpers : shared functions imported from code/digitise_master.py (documented there)
"""
import os, sys; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))   # code/, where digitise_master.py lives
from digitise_master import *                                   # shared helpers (see digitise_master.py)
from digitise_master import save_vert as save                     # vertebrate writer (rounded to 6 decimals)
p=pymupdf.open(A+"brownstein_griffin_2026.pdf")[5]; dr=p.get_drawings()
xt=[g['rect'] for g in dr if g['type']=='f' and g['rect'].width<0.5 and 1<g['rect'].height<2 and max(g['fill'])<0.1]
xtA=sorted([(r.x0+r.x1)/2 for r in xt if 200<r.y0<206]); xtB=sorted([(r.x0+r.x1)/2 for r in xt if 345<r.y0<380])
print("A xticks",len(xtA),"B xticks",len(xtB), xtB[:3])
XA,sxa=lin(xtA,[260,250,240,230,220,210,200,190,180,170]); YA,sya=lin([71.35,106.54,141.73,176.92],[-1.0,-1.5,-2.0,-2.5])
XB,sxb=(lin(xtB,[260,250,240,230,220,210,200,190,180,170]) if len(xtB)==10 else (XA,sxa)); YB,syb=lin([233.74,269.65,305.68,341.72],[-1.6,-1.8,-2.0,-2.2])
def col(g,c): return g['type']=='f' and max(abs(a-b) for a,b in zip(g['fill'],c))<0.01
thero=[g for g in dr if col(g,(1.0,0.498,0.314)) and len(g['items'])>200]
xs=np.linspace(thero[0]['rect'].x0+0.8,thero[0]['rect'].x1-0.8,300); L=poly_profile(thero,xs)
d=pd.DataFrame({"Age_Ma":XA(L[:,0]),"Value":YA(L[:,1])}).dropna(); d=d.groupby(np.round(d.Age_Ma)).mean().reset_index(drop=True)
save(d,"Theropoda_non-avian","phenotypic_rate","none","brownstein_griffin_2026 Fig. 2A (coral Theropoda line, log-transformed morphological rate; clade 95% band not separable from overlapping bands, no error used)",sxa)
blk=[g for g in dr if col(g,(0,0,0)) and len(g['items'])>200 and 220<g['rect'].y0 and g['rect'].y1<360]
band=[g for g in dr if col(g,(0.789,0.793,0.801)) and 220<g['rect'].y0 and g['rect'].y1<362]
xs=np.linspace(blk[0]['rect'].x0+0.8,blk[0]['rect'].x1-0.8,300); L=poly_profile(blk,xs); B=poly_profile(band,xs)
d=pd.DataFrame({"Age_Ma":XB(L[:,0]),"Value":YB(L[:,1]),"Err_hi":YB(B[:,2]),"Err_lo":YB(B[:,3])}).dropna()
d=d.groupby(np.round(d.Age_Ma)).mean().reset_index(drop=True); d["Value_err"]=(d.Err_hi-d.Err_lo)/3.92
print(d.iloc[::15].round(3).to_string())
save(d,"Dinosauria_early","phenotypic_rate","95CI","brownstein_griffin_2026 Fig. 2B (pooled log-transformed morphological rate, black line; grey 95% CI)",sxb)
