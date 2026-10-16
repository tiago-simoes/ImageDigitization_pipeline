"""quintero2024.py -- step 13 of the digitisation pipeline: quintero_2024, Fig. 3A [vector]
Method  : Theria lineages-through-time median + 50% CI from vector paths (log10; CI expanded to a 95% equivalent)
Outputs : outputs/data_tables/clade_origin/Mammalia_early_taxonomic_diversity_BDD.csv
Run     : python code/source_papers/quintero2024.py (from any folder), or python code/digitise_master.py --run quintero2024
Helpers : shared functions imported from code/digitise_master.py (documented there)
"""
import os, sys; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))   # code/, where digitise_master.py lives
from digitise_master import *                                   # shared helpers (see digitise_master.py)
from digitise_master import save_vert as save                     # vertebrate writer (rounded to 6 decimals)
d=pymupdf.open(A+'quintero_2024.pdf'); p=d[2]; dr=p.get_drawings()
X,sx=lin([66.33,99.40,132.47,165.54,198.95],[200,150,100,50,0]); Y,sy=lin([383.44,363.64,343.84,324.03,304.23],[0,1,2,3,4])
line=[g for g in dr if g['type']=='f' and g.get('fill') and abs(g['fill'][0]-0.39)<0.01 and abs(g['fill'][2]-0.39)<0.01 and g['rect'].x0<70 and g['rect'].width>100][0]
band=[g for g in dr if g['type']=='f' and g.get('fill') and abs(g['fill'][0]-0.76)<0.01 and g['rect'].x0<70 and g['rect'].width>100][0]
x0=line['rect'].x0+0.3; x1=line['rect'].x1-0.3; xs=np.arange(x0,x1,0.05)
L=poly_profile([line],xs); Bd=poly_profile([band],xs)
df=pd.DataFrame(dict(a=X(xs),c=Y(L[:,1]),bhi=Y(Bd[:,2]),blo=Y(Bd[:,3])))      # log10 units (y axis inverted)
df["up"]=np.clip(df.bhi-df.c,0,None); df["dn"]=np.clip(df.c-df.blo,0,None)
inb=df.up.notna(); first,last=np.where(inb)[0][[0,-1]]
for col in ("up","dn"):
    v=df[col].copy(); v.iloc[first:last+1]=v.iloc[first:last+1].interpolate(); df[col]=v.fillna(0)   # band gaps (vertical K-Pg drop) interpolated
df["k"]=np.round(df.a)
s=df.groupby("k").agg(Age_Ma=("a","mean"),lc=("c","mean"),up=("up","mean"),dn=("dn","mean")).reset_index(drop=True)
# 50% CI (interquartile) -> 95%-equivalent half-widths in log10 space (x 1.96/0.674); no band -> no error
k=1.96/0.6745
out=pd.DataFrame(dict(Age_Ma=s.Age_Ma,Value=s.lc,Err_lo=s.lc-k*s.dn,Err_hi=s.lc+k*s.up))   # kept in plotted units: log10(number of lineages)
out["Value_err"]=(out.Err_hi-out.Err_lo)/3.92
src=("quintero_2024 Fig. 3A: solid line = median number of therian mammal lineages through time (birth-death diffusion, data-augmented trees); Value = log10(number of lineages), as plotted on the log axis; "
     "shaded band = 50% CI, expanded to a 95%-equivalent interval in log10 space (x1.96/0.674) - vector, exact; clade in figure = Theria (attributed to Mammalia_early)")
print("cal",sx,sy)
save(out,"Mammalia_early","taxonomic_diversity","95CI",src,sx,"BDD")   # Theria, attributed to Mammalia_early
