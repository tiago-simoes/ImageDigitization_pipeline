"""goswami2022.py -- step 22 of the digitisation protocol: goswami_2022, Fig. 3A [raster]
Method  : six root-age rate lines traced by colour and pooled per 1-Myr bin (median; min-max as 95% interval)
Outputs : data_digitized/Placentalia_phenotypic_rate.csv
Note    : 6 root-age lines pooled
Run     : python3 R/digitisation/source_papers/goswami2022.py from the Runs root, or python3 R/digitisation/digitise_master.py --run goswami2022
Helpers : shared functions imported from R/digitisation/digitise_master.py (documented there)
(formerly R/digitise_vert/goswami.py)
"""
import os, sys; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))   # R/digitisation/, where digitise_master.py lives
from digitise_master import *                                   # shared helpers (see digitise_master.py)
from digitise_master import save_vert as save                     # output writer, vertebrate dialect (rounded to 6 decimals)
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
IMG="data_digitized/digitisation/vert_img/"
FIG="data_digitized/digitisation/vert_qc/"; os.makedirs(FIG,exist_ok=True)
a=load(IMG+"goswami_fig3.png")
X,sx=lin([97.5,239.5,382,524.5,666.5],[100,75,50,25,0]); Y,sy=lin([146,253,360.75,468.25,575.75],[.5,.4,.3,.2,.1])
L={"70-75 Ma":((175,30,35),70),"75-80 Ma":((245,225,80),80),"80-85 Ma":((50,50,50),60),
   "85-90 Ma":((192,202,230),40),"90-95 Ma":((75,180,225),60),"95-100 Ma":((68,93,160),55)}
box=(100,66,664,637); excl=[(503,88,666,306),(300,64,480,96),(345,105,425,135)]   # legend, K/Pg + PETM labels
fig,ax=plt.subplots(figsize=(9,8)); ax.imshow(a[60:645,90:670],extent=(90,670,645,60))
ser=[]
for lab,(ref,tol) in L.items():
    R,G,B=a[...,0],a[...,1],a[...,2]; mx=a.max(2); sat=mx-a.min(2)
    rule={"70-75 Ma":(R-G>55)&(R>120)&(R-B>45),
          "75-80 Ma":(R>190)&(G>170)&(B<170)&(R-B>60),
          "80-85 Ma":(mx<110)&(sat<25),
          "85-90 Ma":(B-R>18)&(R>150)&(B>200)&(G<225),
          "90-95 Ma":(B>170)&(R<140)&(G>130)&(B-R>60),
          "95-100 Ma":(B-R>22)&(mx<185)&(R<130)&(G<150)}[lab]
    m=np.zeros(rule.shape,bool); m[box[1]:box[3],box[0]:box[2]]=rule[box[1]:box[3],box[0]:box[2]]
    for ex0,ey0,ex1,ey1 in excl: m[ey0:ey1,ex0:ex1]=False
    lb,n=ndi.label(m,structure=np.ones((3,3))); sz=ndi.sum(m,lb,range(1,n+1)); m=np.isin(lb,[i+1 for i in range(n) if sz[i]>=8])
    if lab=="80-85 Ma":      # drop the dashed K/Pg line (short near-vertical dashes)
        lb,n=ndi.label(m,structure=np.ones((3,3)))
        for i,o in enumerate(ndi.find_objects(lb)):
            h=o[0].stop-o[0].start; w=o[1].stop-o[1].start
            if w<=6 and 280<=o[1].start<=310: m[lb==i+1]=False
    root_max=float(lab.split("-")[1].split()[0]); m[:,:int(np.ceil(97.5+(100-root_max)*5.69))-2]=False   # a line cannot start before its root age
    t=trace(m,minpix=1)
    c,r=interp_gaps(t[:,0],t[:,1])
    gap=np.ones(len(c),bool); gap[np.isin(c,t[:,0])]=False      # hidden stretches > 15 px (~2.6 Myr) left missing, shorter ones interpolated
    lbl,ng=ndi.label(gap)
    for j in range(1,ng+1):
        if (lbl==j).sum()<=15: gap[lbl==j]=False
    c,r=c[~gap],r[~gap]
    s=to_series(X(c),Y(r),1.0); s["line"]=lab; ser.append(s)
    ax.plot(c,r,'-',lw=.6,label=lab)
    print(lab,"cols",len(t),"span %.1f-%.1f"%(s.Age_Ma.max(),s.Age_Ma.min()))
ax.legend(fontsize=7); ax.set_title("Placentalia rates (Goswami 2022 Fig. 3A) - QC: traced lines"); fig.savefig(FIG+"Placentalia_phenotypic_rate_QC_digitised.png",dpi=110,bbox_inches="tight")
D=pd.concat(ser); D["k"]=np.round(D.Age_Ma)
D.to_csv("data_digitized/digitisation/goswami_lines.csv",index=False)
s=D.groupby("k").agg(Age_Ma=("k","first"),Value=("Value","median"),Err_lo=("Value","min"),Err_hi=("Value","max"),n=("Value","size")).reset_index(drop=True)
s["Value_err"]=(s.Err_hi-s.Err_lo)/3.92; print(s.n.value_counts().sort_index().to_dict())
src=("goswami_2022 Fig. 3A: rates of cranial evolution through time, all 6 root-age lines (70-75 ... 95-100 Ma) traced (raster) and pooled per 1-Myr bin: "
     "Value = median across lines, Err_lo/Err_hi = min/max across lines (treated as 95% interval)")
save(s,"Placentalia","phenotypic_rate","95CI",src,float(np.hypot(sx,0.5/5.69)))
