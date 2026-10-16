"""close2019_1_extract_localities.py -- step 2 of the digitisation pipeline: close_2019, Fig. 1a [vector]
Method  : reads the locality points of Fig. 1a from the vector PDF and calibrates age and richness on the axis ticks
Note    : locality richness points -> outputs/intermediates/close2019/close2019_fig1a_localities.csv; source crop + QC panel
Run     : python code/source_papers/close2019_1_extract_localities.py (from any folder), or python code/digitise_master.py --run close2019_1_extract_localities
"""
import pymupdf, numpy as np, pandas as pd
import os, sys; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))   # code/, where digitise_master.py lives
from digitise_master import PDF as PDF_DIR, OUT as DAT, FIGS, FIGIMG, CLOSE2019   # paths (digitise_master.py, section 1)
PDF=PDF_DIR+"close_2019.pdf"; DIG=CLOSE2019; WORK=CLOSE2019; QC=FIGS+"taxonomic_diversity/"   # intermediates in outputs/intermediates/close2019/
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
# source crop + QC panel (source | re-plot of the extracted dots)
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from PIL import Image
crop=QC+"Tetrapoda_nonflying_diversity_Close2019_Fig1a.png"
p.get_pixmap(dpi=400,clip=pymupdf.Rect(300,50,552,320)).save(crop)
fig,q=plt.subplots(1,2,figsize=(13,4.6),gridspec_kw=dict(width_ratios=[1,1.55]))
q[0].imshow(Image.open(crop)); q[0].axis("off"); q[0].set_title("Source: Close et al. 2019 Fig. 1a",fontsize=8)
st_,ex_=d[d.type=="standard"],d[d.type=="exceptional"]
q[1].scatter(st_.Age_Ma,st_.Richness,s=3,color="#b5b8bf",lw=0); q[1].scatter(ex_.Age_Ma,ex_.Richness,s=12,color="red",lw=0)
q[1].set_xlim(410,-3); q[1].set_xlabel("Age (Ma)"); q[1].set_ylabel("Species richness")
q[1].set_title(f"Re-plot of extracted vector dots (n={len(d)}; {len(ex_)} exceptional)",fontsize=9)
fig.tight_layout(); fig.savefig(crop.replace(".png","_QC_digitised.png"),dpi=110); plt.close(fig)
