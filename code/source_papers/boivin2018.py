"""boivin2018.py -- step 28 of the digitisation protocol: boivin-et_al_2018, Fig. 3 [vector]
Method  : top-most envelope of the stacked areas (sum of all clades) read from vector paths
Outputs : data_digitized/Echinoidea_Irregularia_phenotypic_disparity.csv, data_digitized/Echinoidea_Irregularia_taxonomic_diversity.csv
Run     : python3 R/digitisation/source_papers/boivin2018.py from the Runs root, or python3 R/digitisation/digitise_master.py --run boivin2018
Helpers : shared functions imported from R/digitisation/digitise_master.py (documented there)
(formerly R/digitise_invert/boivin.py)
"""
import os, sys; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))   # R/digitisation/, where digitise_master.py lives
from digitise_master import *                                   # shared helpers (see digitise_master.py)
from digitise_master import save_invert as save, OUT_INVERT as OUT, IMG_INVERT as IMG   # invertebrate dialect (4 decimals)
p=pymupdf.open(PDF+"boivin-et_al_2018.pdf")[7]; dr=p.get_drawings()
X,sx=lin([225.95,297.63,388.24,545.26],[182.7,174.2,163.5,145])
for par,yr,(yp,yv) in [("taxonomic_diversity",(70,281.5),([280.53,239.11,197.34,155.56,113.78],[0,50,100,150,200])),
                       ("phenotypic_disparity",(282,479.6),([478.6,414.03,349.37],[0,.005,.01]))]:
    Y,sy=lin(yp,yv)
    V=[]
    for g in dr:
        r=g['rect']
        if g['type'] in ('f','fs') and r.x0>224 and r.y0>=yr[0]-2 and r.y1<=yr[1]+1 and r.width>50:
            if g.get('fill')==(0.0,0.0,0.0): continue
            if r.width>315 and r.height>0.9*(yr[1]-yr[0]-10) and abs(r.y0-yr[0])<3: continue   # panel background
            V+= [q for q in path_pts(g)]
    xs=sorted(set(round(x,1) for x,y in V))
    cl=[]
    for x in xs:
        if not cl or x-cl[-1][-1]>1.5: cl.append([x])
        else: cl[-1].append(x)
    rows=[]
    for c in cl:
        xm=np.mean(c); ys=[y for x,y in V if c[0]-0.05<=round(x,1)<=c[-1]+0.05]
        rows.append((xm,min(ys),len(ys)))
    rows=[r for r in rows if r[2]>=20]
    d=pd.DataFrame(rows,columns=["px","py","n"]); d["Age_Ma"]=X(d.px).round(2); d["Value"]=Y(d.py)
    print(par,sx,sy); print(d.round(4).T.to_string())
    d=d[(d.px>227)&(d.px<544)]   # drop frame-edge padding vertices
    et="none"
    if par=="taxonomic_diversity":
        d["Value_err"]=np.sqrt(d.Value); d["Err_lo"]=d.Value-d.Value_err; d["Err_hi"]=d.Value+d.Value_err; et="Poisson"
    save(d,"Echinoidea_Irregularia",par,sx,et,"boivin_et_al_2018 Fig. 3 "+("top" if par=="taxonomic_diversity" else "bottom")+" (top-most envelope = sum of all clades)")
