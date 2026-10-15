"""close2019_1_extract_localities.py -- step 2 of the digitisation protocol: close_2019, Fig. 1a [vector]
Method  : reads the locality points of Fig. 1a from the vector PDF and calibrates age and richness on the axis ticks
Note    : locality richness points -> digitisation/close2019_fig1a_localities.csv
Run     : python3 R/digitisation/source_papers/close2019_1_extract_localities.py from the Runs root, or python3 R/digitisation/digitise_master.py --run close2019_1_extract_localities
(formerly R/digitise_close2019/extract.py)
"""
import pymupdf, numpy as np, pandas as pd
PDF="data_raw/articles_PDF/close_2019.pdf"   # run from the Runs root: python3 R/digitisation/source_papers/close2019_<n>_*.py
DIG="data_digitized/digitisation/"; DAT="data_digitized/"; QC="data_raw/taxonomic_diversity/"
WORK="data_digitized/digitisation/close2019_work/"   # intermediates; f5/img-000.jpg = embedded Fig. 5 raster of close_2019.pdf p. 4 (pdfimages -j -f 4 -l 4 ... f5/img)
F=PDF
p=pymupdf.open(F)[1]
bx,ax_=np.polyfit([336.44,388.62,440.81,493.0,545.19],[400,300,200,100,0],1)
by,ay=np.polyfit([170.38,144.45,118.51,92.58,66.63],[0,20,40,60,80],1)
rows=[]
for g in p.get_drawings():
    if g['type']!='f' or not g.get('fill'): continue
    r=g['rect']
    if not (333<r.x0 and r.x1<550 and 55<r.y0 and r.y1<175 and 2.6<r.width<2.9): continue
    col=g['fill']; red=col[0]>0.9 and col[1]<0.1
    grey=abs(col[0]-0.77)<0.02
    if not (red or grey): continue
    cx,cy=(r.x0+r.x1)/2,(r.y0+r.y1)/2
    rows.append((cx,cy,"exceptional" if red else "standard"))
d=pd.DataFrame(rows,columns=["px","py","type"])
# drop the two legend markers (top-left legend box)
leg=(d.px<350)&(d.py<80); print("legend dots dropped:",leg.sum()); d=d[~leg]
d["Age_Ma"]=ax_+bx*d.px; d["Richness_raw"]=ay+by*d.py
d["Richness"]=d.Richness_raw.round().astype(int)
res=d.Richness_raw-d.Richness
print("n dots",len(d), d.type.value_counts().to_dict()); print("integer residual (species) mean abs",round(res.abs().mean(),3),"max",round(res.abs().max(),3))
print("age range",d.Age_Ma.min().round(2),d.Age_Ma.max().round(2),"richness range",d.Richness.min(),d.Richness.max())
d.round(4).to_csv(DIG+"close2019_fig1a_localities.csv",index=False)
