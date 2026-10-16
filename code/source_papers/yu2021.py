"""yu2021.py -- step 11 of the digitisation pipeline: yu_2021_SI, Fig. S40 [vector]
Method  : net-diversification rate bars (vector) per time bin with SE; split at 100 Ma into Aves_stem / Aves_crown
Outputs : outputs/data_tables/clade_origin/Aves_stem_taxonomic_diversity.csv, outputs/data_tables/clade_origin/Aves_crown_taxonomic_diversity.csv
Note    : split at 100 Ma
Run     : python code/source_papers/yu2021.py (from any folder), or python code/digitise_master.py --run yu2021
Helpers : shared functions imported from code/digitise_master.py (documented there)
"""
import os, sys; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))   # code/, where digitise_master.py lives
from digitise_master import *                                   # shared helpers (see digitise_master.py)
from digitise_master import save_vert as save                     # vertebrate writer (rounded to 6 decimals)
p=pymupdf.open(A+"yu_2021_SI.pdf")[62]; dr=p.get_drawings()
xt=[479.46,462.3,445.14,427.98,410.81,393.65,376.49,359.33,342.16,325.0,307.84,290.67,273.51,256.34,239.18,222.02,204.94,187.78,170.62,153.45,136.29,119.13,101.97]
X,sx=lin(xt,list(range(0,221,10))); Y,sy=lin([456.84,427.83,398.82,369.82,340.81,311.8,282.78],[0,.05,.10,.15,.20,.25,.30])
red=[g for g in dr if g['type']=='s' and g.get('color') and abs(g['color'][0]-0.78)<0.01 and abs(g['color'][1]-0.145)<0.01 and g['rect'].y0>280]
H=[]
for g in red:
    for it in g['items']:
        P=[q for q in it[1:] if hasattr(q,'x')]
        if len(P)==2 and abs(P[0].y-P[1].y)<0.05 and abs(P[0].x-P[1].x)>0.5: H.append((min(P[0].x,P[1].x),max(P[0].x,P[1].x),P[0].y))
H=sorted(set((round(a,2),round(b,2),round(c,2)) for a,b,c in H))
# red CI bars: reddish fills (opaque, colour varies where blended with tree)
def reddish(c): return c and c[0]>0.6 and c[0]-c[1]>0.25 and c[0]-c[2]>0.25
F=[g['rect'] for g in dr if g['type']=='f' and reddish(g.get('fill')) and g["rect"].y0>270 and g["rect"].y1<=458 and g["rect"].height>0.5]
rows=[]
for x0,x1,y in H:
    xm=(x0+x1)/2; bars=[r for r in F if r.x0<=xm<=r.x1]
    top=min([r.y0 for r in bars],default=np.nan); bot=max([r.y1 for r in bars],default=np.nan)
    occl = 300<xm<363 and 443.5<bot<448.5   # lower bar edge hidden under the opaque Tyrannosauroidea triangle (y 447.6-458.8)
    if occl: bot=456.84                       # visible bars in the same region reach the zero line -> lower bound set to 0
    rows.append(dict(occluded=occl,Age_old=X(x0),Age_young=X(x1),Age_Ma=X(xm),Value=Y(y),Err_hi=Y(top),Err_lo=Y(bot),nbars=len(bars)))
d=pd.DataFrame(rows); d["Value_err"]=(d.Err_hi-d.Err_lo)/2
print("cal",sx,sy); print(d.round(3).to_string())
d.to_csv(WORK+"yu_S40_bins.csv",index=False)
src="yu_2021_SI Fig. S40 (solid red line = net diversification rate of birds, Myr^-1, step line per time bin; shaded bars = interval, treated as +/-1 error); split at 100 Ma"
save(d[d.Age_Ma>100],"Aves_stem","taxonomic_diversity","SE",src+" - part older than 100 Ma",sx)
save(d[d.Age_Ma<=100],"Aves_crown","taxonomic_diversity","SE",src+" - part 100 Ma and younger",sx)
