"""ghezelayagh2022_edfig5e.py -- step 9 of the digitisation protocol: ghezelayagh_2022, Extended Data Fig. 5e [vector]
Method  : posterior-mean net-diversification line and 95% credible band read from the vector paths
Outputs : data_digitized/Acanthomorpha_taxonomic_diversity.csv
Note    : ED Fig. 5e used instead of 5b (PI)
Run     : python3 R/digitisation/source_papers/ghezelayagh2022_edfig5e.py from the Runs root, or python3 R/digitisation/digitise_master.py --run ghezelayagh2022_edfig5e
Helpers : shared functions imported from R/digitisation/digitise_master.py (documented there)
(formerly R/digitise_vert/ghez5e.py)
"""
import os, sys; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))   # R/digitisation/, where digitise_master.py lives
from digitise_master import *                                   # shared helpers (see digitise_master.py)
from digitise_master import save_vert as save                     # output writer, vertebrate dialect (rounded to 6 decimals)
p=pymupdf.open(A+"ghezelayagh_2022.pdf")[14]
X,sx=lin([90.15,116.18,141.63,168.06,193.42,218.84,246.09,272.19],[140,120,100,80,60,40,20,0])
Y,sy=lin([464.42,449.52,434.52,419.6,404.64,389.74],[0,.05,.10,.15,.20,.25])
dr=p.get_drawings()
line=[g for g in dr if g['type']=='f' and abs(g['fill'][0]-0.278)<0.01 and abs(g['fill'][2]-0.639)<0.01 and g['rect'].y0>365 and g['rect'].y1<480]
band=[g for g in dr if g['type']=='f' and abs(g['fill'][0]-0.775)<0.01 and abs(g['fill'][2]-0.882)<0.01 and g['rect'].y0>365 and g['rect'].y1<480]
x0=max(min(g['rect'].x0 for g in line),75.3)+0.6; x1=max(g['rect'].x1 for g in line)-0.6
xs=np.linspace(x0,x1,400)
L=poly_profile(line,xs); B=poly_profile(band+line,xs)
ok=~np.isnan(L[:,1])
d=pd.DataFrame({"Age_Ma":X(L[ok,0]),"Value":Y(L[ok,1]),"Err_hi":Y(B[ok,2]),"Err_lo":Y(B[ok,3])})
d=d.groupby(np.round(d.Age_Ma)).mean().reset_index(drop=True)
d["Value_err"]=(d.Err_hi-d.Err_lo)/3.92
print("cal",sx,sy); print(d.iloc[::10].round(4).to_string())
save(d,"Acanthomorpha","taxonomic_diversity","95CI","ghezelayagh_2022 Extended Data Fig. 5e (posterior mean net-diversification rate, TESS-CoMET; 95% credible interval; used instead of ED Fig. 5b as confirmed by PI)",sx)
