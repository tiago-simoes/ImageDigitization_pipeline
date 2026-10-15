"""close2019_4_squamata_extract.py -- step 5 of the digitisation protocol: close_2019, Fig. 5c [raster]
Method  : isolates the green squamate dots of Fig. 5c (raster) with a colour mask and records their upper envelope
Note    : squamate dot envelope (green mask)
Run     : python3 R/digitisation/source_papers/close2019_4_squamata_extract.py from the Runs root, or python3 R/digitisation/digitise_master.py --run close2019_4_squamata_extract
(formerly R/digitise_close2019/sq_extract.py)
"""
import numpy as np, pandas as pd
from PIL import Image
from scipy import ndimage as ndi
PDF="data_raw/articles_PDF/close_2019.pdf"   # run from the Runs root: python3 R/digitisation/source_papers/close2019_<n>_*.py
DIG="data_digitized/digitisation/"; DAT="data_digitized/"; QC="data_raw/taxonomic_diversity/"
WORK="data_digitized/digitisation/close2019_work/"   # intermediates; f5/img-000.jpg = embedded Fig. 5 raster of close_2019.pdf p. 4 (pdfimages -j -f 4 -l 4 ... f5/img)
import os; os.makedirs(WORK+"f5",exist_ok=True)
a=np.asarray(Image.open(WORK+"f5/img-000.jpg").convert("RGB")).astype(int)
X0,Y0=330.90655517578125,64.2503662109375; SX=895/(545.5865478515625-X0); SY=1297/(375.5303649902344-Y0)
px=lambda xpt:(xpt-X0)*SX; py=lambda ypt:(ypt-Y0)*SY
# calibration (pt) -> pixels
xt=[px(v) for v in [333.95,386.47,439.03,491.59,544.2]]; yt=[py(v) for v in [360.41,337.63,314.86,292.09,269.32]]
bx,ax_=np.polyfit(xt,[400,300,200,100,0],1); by,ay=np.polyfit(yt,[0,20,40,60,80],1)
top,bot=int(py(262)),int(py(361.5))          # panel c plot region (below panel-b frame, above time bar)
R,G,B=a[...,0],a[...,1],a[...,2]
green=(G-np.maximum(R,B)>35)&(G>110)
green[:top]=0; green[bot:]=0
# exclude the legend box (top-left of panel c)
green[top:int(py(287)),:int(px(400))]=0
lab,n=ndi.label(green); sizes=ndi.sum(green,lab,range(1,n+1))
# dot diameter from small isolated components
iso=[ndi.find_objects(lab)[i] for i in range(n) if 30<sizes[i]<120]
dia=np.median([ (s[0].stop-s[0].start + s[1].stop-s[1].start)/2 for s in iso]); r=dia/2
print("panel rows",top,bot,"components",n,"dot diameter px",round(dia,1),"species per px",round(-by,3),"Myr per px",round(-bx,3))
# upper envelope: topmost green pixel per column -> local minima = dot centres on the envelope
H,Wd=green.shape
ytop=np.full(Wd,np.nan)
for x in range(Wd):
    ys=np.where(green[:,x])[0]
    if len(ys): ytop[x]=ys.min()
cent=[]
k=int(round(r))
for x in range(Wd):
    if np.isnan(ytop[x]): continue
    win=ytop[max(0,x-k):x+k+1]
    if ytop[x]<=np.nanmin(win)+0.01:
        # merge plateaus: keep middle of run of equal minima
        cent.append((x,ytop[x]))
c=pd.DataFrame(cent,columns=["px","ytop"])
c["grp"]=(c.px.diff()>1).cumsum(); c=c.groupby("grp").agg(px=("px","mean"),ytop=("ytop","min")).reset_index(drop=True)
c["py"]=c.ytop+r
c["Age_Ma"]=ax_+bx*c.px; c["Richness_raw"]=ay+by*c.py; c["Richness"]=c.Richness_raw.round().clip(lower=1).astype(int)
print("envelope dots",len(c), "age range",c.Age_Ma.max().round(1),c.Age_Ma.min().round(1),"max",c.Richness.max())
c.round(3).to_csv(WORK+"close2019_fig5c_squamata_envelope_dots.csv",index=False)
np.save(WORK+"f5/green_mask.npy",green)
