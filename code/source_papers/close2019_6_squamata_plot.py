"""close2019_6_squamata_plot.py -- step 7 of the digitisation pipeline: close_2019, Fig. 5c [qc]
Method  : summary plot of squamates vs all non-flying tetrapods; source crop + QC panel
Note    : summary plots + QC panel
Run     : python code/source_papers/close2019_6_squamata_plot.py (from any folder), or python code/digitise_master.py --run close2019_6_squamata_plot
"""
import numpy as np, pandas as pd, matplotlib, pymupdf
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from PIL import Image
import os, sys; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))   # code/, where digitise_master.py lives
from digitise_master import PDF as PDF_DIR, OUT as DAT, FIGS, FIGIMG, CLOSE2019   # paths (digitise_master.py, section 1)
PDF=PDF_DIR+"close_2019.pdf"; DIG=CLOSE2019; WORK=CLOSE2019; QC=FIGS+"taxonomic_diversity/"   # intermediates in outputs/intermediates/close2019/
s=pd.read_csv(DAT+"Squamata_stage_max_local_richness.csv"); a=pd.read_csv(DAT+"Tetrapoda_nonflying_stage_max_local_richness.csv")
PER=[("D",382.7,358.9,"#CB8C37"),("C",358.9,298.9,"#67A599"),("P",298.9,251.902,"#F04028"),("Tr",251.902,201.4,"#812B92"),("J",201.4,145.0,"#34B2C9"),
     ("K",145.0,66.0,"#7FC64E"),("Pg",66.0,23.03,"#FD9A52"),("Ng",23.03,2.58,"#FFE619"),("Q",2.58,0,"#F9F97F")]
INK="#2b2b2b"; MUTED="#6b6b6b"; GRID="#e9e9e9"; C3="#1baf7a"
plt.rcParams.update({"font.size":9,"axes.edgecolor":MUTED,"xtick.color":MUTED,"ytick.color":MUTED,"axes.labelcolor":INK})
fig,(ax,axb)=plt.subplots(2,1,figsize=(8,4.6),sharex=True,gridspec_kw=dict(height_ratios=[1,0.07],hspace=0.05))
for _,x0,x1,_c in PER: ax.axvline(x1,color=GRID,lw=0.8,zorder=0)
ax.axvline(66.0,color=MUTED,lw=0.8,ls=(0,(3,3)))
ax.plot(a.Mid_Ma,a.Max_local_richness,"-",color="#b5b5b5",lw=1.2,zorder=2,label="All non-flying tetrapods (Fig. 1a)")
unc=s.Note.fillna("").str.contains("narrower|covered")
ax.plot(s.Mid_Ma,s.Max_local_richness,"-",color=C3,lw=1.8,zorder=3,label="Squamata (Fig. 5c)")
ax.plot(s.Mid_Ma[~unc],s.Max_local_richness[~unc],"o",color=C3,ms=4,mec="white",mew=.5,zorder=4,label="dot top inside stage")
ax.plot(s.Mid_Ma[unc],s.Max_local_richness[unc],"o",mfc="white",mec=C3,ms=4,mew=1.2,zorder=4,label="stage narrower than a dot / value from envelope")
ax.set_ylabel("Maximum species per locality"); ax.set_ylim(0,90); ax.spines[["top","right"]].set_visible(False)
ax.legend(loc="upper left",frameon=False,fontsize=7.5)
ax.set_title("Squamata – maximum local species richness per stage (Close et al. 2019, Fig. 5c)",loc="left",fontsize=9.5,color=INK)
for n,x0,x1,c in PER:
    axb.add_patch(Rectangle((x1,0),x0-x1,1,facecolor=c,edgecolor="white",lw=1))
    if x0-x1>8: axb.text((x0+x1)/2,.5,n,ha="center",va="center",fontsize=7.5,color=INK)
axb.set_ylim(0,1); axb.set_yticks([]); [axb.spines[k].set_visible(False) for k in axb.spines]
axb.set_xlim(390,-3); axb.set_xticks(range(400,-1,-50)); axb.set_xlabel("Age (Ma)")
fig.savefig(WORK+"Squamata_max_local_richness_Close2019.png",dpi=300,bbox_inches="tight"); fig.savefig(WORK+"Squamata_max_local_richness_Close2019.pdf",bbox_inches="tight")
# source crop + QC
F=PDF
pymupdf.open(F)[3].get_pixmap(dpi=400,clip=pymupdf.Rect(318,262,548,392)).save(QC+"Squamata_diversity_Close2019_Fig5c.png")
fig,q=plt.subplots(1,2,figsize=(13,4.4),gridspec_kw=dict(width_ratios=[1.1,1]))
if os.path.exists(FIGIMG+"close_2019_Fig5c_dot_tops_zoom.png"): q[0].imshow(Image.open(FIGIMG+"close_2019_Fig5c_dot_tops_zoom.png"))
q[0].axis("off"); q[0].set_title("Fig. 5c raster (zoom): detected green dot tops (red rings)",fontsize=8)
q[1].plot(s.Mid_Ma,s.Max_local_richness,"o-",color=C3,ms=3); q[1].invert_xaxis(); q[1].set_xlabel("Age (Ma)"); q[1].set_ylabel("Max species per locality"); q[1].set_title("Stage maxima from detected tops",fontsize=8)
fig.tight_layout(); fig.savefig(QC+"Squamata_diversity_Close2019_Fig5c_QC_digitised.png",dpi=110)
