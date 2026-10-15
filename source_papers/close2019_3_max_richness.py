"""close2019_3_max_richness.py -- step 4 of the digitisation protocol: close_2019, Fig. 1a [derived]
Method  : maximum local species richness per stage with bootstrap 95% CI (1000 replicates) + summary plot
Outputs : data_digitized/Tetrapoda_nonflying_stage_max_local_richness.csv
Note    : maximum local richness per stage + bootstrap CI
Run     : python3 R/digitisation/source_papers/close2019_3_max_richness.py from the Runs root, or python3 R/digitisation/digitise_master.py --run close2019_3_max_richness
(formerly R/digitise_close2019/maxrich.py)
"""
import numpy as np, pandas as pd, matplotlib, re
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
PDF="data_raw/articles_PDF/close_2019.pdf"   # run from the Runs root: python3 R/digitisation/source_papers/close2019_<n>_*.py
DIG="data_digitized/digitisation/"; DAT="data_digitized/"; QC="data_raw/taxonomic_diversity/"
WORK="data_digitized/digitisation/close2019_work/"   # intermediates; f5/img-000.jpg = embedded Fig. 5 raster of close_2019.pdf p. 4 (pdfimages -j -f 4 -l 4 ... f5/img)
import os; os.makedirs(WORK,exist_ok=True)
st=pd.read_csv(DIG+"Tetrapoda_nonflying_stage_summed_richness.csv")[["Stage","Period","Base_Ma","Top_Ma","Mid_Ma"]]
d=pd.read_csv(DIG+"close2019_fig1a_localities.csv"); d["Age_Ma"]=d.Age_Ma.clip(lower=0)
idx=np.searchsorted(-st.Base_Ma.values,-d.Age_Ma.values,side="right")-1
d["Stage"]=st.Stage.values[np.clip(idx,0,len(st)-1)]
rng=np.random.default_rng(1); N=20; B=1000
rows=[]
for s in st.Stage:
    r=d.loc[d.Stage==s,"Richness"].values; n=len(r)
    row=dict(Stage=s,N_localities=n)
    if n:
        row.update(Max_local_richness=int(r.max()),Q95_local_richness=float(np.quantile(r,.95)))
        # bootstrap 95% CI of the maximum (resample localities)
        bm=np.array([rng.choice(r,n,replace=True).max() for _ in range(B)])
        row.update(Max_boot_lo=float(np.quantile(bm,.025)),Max_boot_hi=float(np.quantile(bm,.975)))
        if n>=N:   # rarefied expected maximum at N localities (subsampling without replacement)
            rm=np.array([rng.choice(r,N,replace=False).max() for _ in range(B)])
            row.update(**{f"RarefiedMax_n{N}":float(rm.mean()),f"RarefiedMax_n{N}_lo":float(np.quantile(rm,.025)),f"RarefiedMax_n{N}_hi":float(np.quantile(rm,.975))})
    rows.append(row)
o=st.merge(pd.DataFrame(rows),on="Stage")
o.round(2).to_csv(DAT+"Tetrapoda_nonflying_stage_max_local_richness.csv",index=False)
print(o[["Stage","N_localities","Max_local_richness","Max_boot_lo","Max_boot_hi",f"RarefiedMax_n{N}"]].round(1).to_string(index=False))
# tidy series for pipeline (raw maximum; bootstrap CI as error)
t=o[o.N_localities>0]
pd.DataFrame({"Clade":"Tetrapoda_nonflying","Parameter":"max_local_richness","Age_Ma":t.Mid_Ma.round(3),"Value":t.Max_local_richness,
  "Value_err":((t.Max_boot_hi-t.Max_boot_lo)/3.92).round(3),"Err_lo":t.Max_boot_lo,"Err_hi":t.Max_boot_hi,"Err_type":"95CI",
  "Source":"close_2019 Fig. 1a (maximum species richness per locality per ICS stage; bootstrap 95% CI, 1000 reps)","Age_cal_sd":0.0}).to_csv(WORK+"Tetrapoda_nonflying_max_local_richness.csv",index=False)
# plot
PER=[("D",382.7,358.9,"#CB8C37"),("C",358.9,298.9,"#67A599"),("P",298.9,251.902,"#F04028"),("Tr",251.902,201.4,"#812B92"),("J",201.4,145.0,"#34B2C9"),
     ("K",145.0,66.0,"#7FC64E"),("Pg",66.0,23.03,"#FD9A52"),("Ng",23.03,2.58,"#FFE619"),("Q",2.58,0,"#F9F97F")]
INK="#2b2b2b"; MUTED="#6b6b6b"; GRID="#e9e9e9"; C1="#2a78d6"; C2="#eb6834"
plt.rcParams.update({"font.size":9,"axes.edgecolor":MUTED,"xtick.color":MUTED,"ytick.color":MUTED,"axes.labelcolor":INK})
fig,(ax,axb)=plt.subplots(2,1,figsize=(8,4.8),sharex=True,gridspec_kw=dict(height_ratios=[1,0.07],hspace=0.05))
for _,a,b,_c in PER: ax.axvline(b,color=GRID,lw=0.8,zorder=0)
for x in [252.2,66.0]: ax.axvline(x,color=MUTED,lw=0.8,ls=(0,(3,3)),zorder=1)
ax.vlines(t.Mid_Ma,t.Max_boot_lo,t.Max_boot_hi,color=C1,alpha=.35,lw=2,zorder=2)
ax.plot(t.Mid_Ma,t.Max_local_richness,"-o",color=C1,ms=3.5,lw=1.4,mec="white",mew=.5,zorder=4,label="Maximum local richness (bootstrap 95% CI)")
rr=o   # NaN rows break the line at stages with < N localities
ax.fill_between(rr.Mid_Ma,rr[f"RarefiedMax_n{N}_lo"],rr[f"RarefiedMax_n{N}_hi"],color=C2,alpha=.15,lw=0,zorder=1)
ax.plot(rr.Mid_Ma,rr[f"RarefiedMax_n{N}"],"-",color=C2,lw=1.6,zorder=3,label=f"Rarefied maximum at {N} localities (stages with ≥{N})")
ax.set_ylabel("Species per locality"); ax.set_ylim(0,90); ax.spines[["top","right"]].set_visible(False)
ax.legend(loc="upper left",frameon=False,fontsize=7.5)
ax.set_title("Non-flying terrestrial tetrapods – maximum local species richness per stage (Close et al. 2019, Fig. 1a)",loc="left",fontsize=9.5,color=INK)
for n,a,b,c in PER:
    axb.add_patch(Rectangle((b,0),a-b,1,facecolor=c,edgecolor="white",lw=1))
    if a-b>8: axb.text((a+b)/2,.5,n,ha="center",va="center",fontsize=7.5,color=INK)
axb.set_ylim(0,1); axb.set_yticks([]); [axb.spines[s].set_visible(False) for s in axb.spines]
axb.set_xlim(390,-3); axb.set_xticks(range(400,-1,-50)); axb.set_xlabel("Age (Ma)")
fig.savefig(WORK+"Tetrapoda_nonflying_max_local_richness_Close2019.png",dpi=300,bbox_inches="tight"); fig.savefig(WORK+"Tetrapoda_nonflying_max_local_richness_Close2019.pdf",bbox_inches="tight")
