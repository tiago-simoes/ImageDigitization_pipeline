"""qc_panels_vert.py -- step 39 of the digitisation protocol [qc]
Method  : source crops + QC panels (red line, grey band) in data_raw/<parameter folder>/
Note    : source crops + QC panels of the vertebrate figure series
Run     : python3 R/digitisation/utilities/qc_panels_vert.py from the Runs root, or python3 R/digitisation/digitise_master.py --run qc_panels_vert
Helpers : shared functions imported from R/digitisation/digitise_master.py (documented there)
(formerly R/digitise_vert/assets.py)

Source-figure crops + QC panels (source | digitised series) for the vertebrate key rows (figure-derived series only).
"""
import os, sys; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))   # R/digitisation/, where digitise_master.py lives
from digitise_master import *                                   # shared helpers (see digitise_master.py)
from digitise_master import save_vert as save                     # output writer, vertebrate dialect (rounded to 6 decimals)
import re
from PIL import Image
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
V=""
IMG="data_digitized/digitisation/vert_img/"; DAT="data_digitized/"; OUTD="data_raw/"
FOLD={"phenotypic_rate":"phenotypic_rates","phenotypic_disparity":"phenotypic_disparity","taxonomic_diversity":"taxonomic_diversity"}
SHORT={"phenotypic_rate":"rates","phenotypic_disparity":"disparity","taxonomic_diversity":"diversity"}
def pdfclip(f,pg,rect,dpi=300):
    pix=pymupdf.open(A+f+".pdf")[pg].get_pixmap(dpi=dpi,clip=pymupdf.Rect(*rect),annots=False); return Image.frombytes("RGB",(pix.width,pix.height),pix.samples)
def imgclip(f,box=None):
    im=Image.open(IMG+f).convert("RGB"); return im.crop(box) if box else im
S=[ # clade, parameter, data_type slug, tag, source image
 ("Acanthomorpha","phenotypic_disparity","","Ghezelayagh2022_Fig2",lambda:pdfclip("ghezelayagh_2022",3,(60,270,445,525))),
 ("Acanthomorpha","taxonomic_diversity","","Ghezelayagh2022_EDFig5e",lambda:pdfclip("ghezelayagh_2022",14,(55,370,285,485))),
 ("Tetrapoda_stem","phenotypic_disparity","","Rawson2022_Fig5A",lambda:imgclip("rawson_fig5.png",(0,0,882,700))),
 ("Archosauromorpha_stem","phenotypic_disparity","","Foth2021_Fig9A_row1",lambda:imgclip("foth9A_r1.png")),
 ("Archosauromorpha_stem","taxonomic_diversity","","Sookias2026_Fig1b",lambda:imgclip("sookias_fig1.png",(0,850,1550,1400))),
 ("Dinosauria_early","phenotypic_rate","","BrownsteinGriffin2026_Fig2B",lambda:pdfclip("brownstein_griffin_2026",5,(150,222,485,390))),
 ("Dinosauria_early","phenotypic_disparity","","Foth2021_Fig9A_row2",lambda:imgclip("foth9A_r2.png")),
 ("Dinosauria_early","taxonomic_diversity","PyRate","Sookias2026_Fig1f",lambda:imgclip("sookias_fig1.png",(0,1930,1550,2480))),
 ("Dinosauria_early","taxonomic_diversity","FBDD","Quintero2025SI_FigS9",lambda:imgclip("quintero25_S9.png",(690,930,1200,1200))),
 ("Theropoda_non-avian","phenotypic_rate","","BrownsteinGriffin2026_Fig2A",lambda:pdfclip("brownstein_griffin_2026",5,(150,40,485,228))),
 ("Theropoda_non-avian","phenotypic_disparity","","Foth2021_Fig9A_row3",lambda:imgclip("foth9A_r3.png")),
 ("Theropoda_non-avian","taxonomic_diversity","","Sookias2026_Fig1h",lambda:imgclip("sookias_fig1.png",(0,2470,1550,2974))),
 ("Aves_stem","taxonomic_diversity","","Yu2021SI_FigS40",lambda:pdfclip("yu_2021_SI",62,(40,40,500,485))),
 ("Aves_crown","phenotypic_rate","continuous_data_LMM","CrouchTobias2022_Fig3a",lambda:imgclip("crouch_fig3.png",(0,0,1010,1010))),
 ("Aves_crown","phenotypic_rate","continuous_data_GMM","FeliceGoswami2017_Fig1C",lambda:imgclip("felice_fig1.png",(0,630,679,1278))),
 ("Aves_crown","phenotypic_disparity","","Navalon2022_Fig1a",lambda:pdfclip("navalon_2022",1,(300,90,440,205))),
 ("Aves_crown","taxonomic_diversity","","Yu2021SI_FigS40",lambda:pdfclip("yu_2021_SI",62,(40,40,500,485))),
 ("Mammalia_early","phenotypic_rate","","Close2015_Fig1B",lambda:imgclip("close15_fig1.png",(0,940,1876,2300))),
 ("Mammalia_early","phenotypic_disparity","","Close2015_Fig1C",lambda:imgclip("close15_fig1.png",(0,2000,1876,2591))),
 ("Mammalia_early","taxonomic_diversity","","Close2015_Fig1D",lambda:imgclip("close15_fig1.png",(0,2280,1876,2591))),
 ("Placentalia","phenotypic_rate","","Goswami2022_Fig3A",lambda:imgclip("goswami_fig3.png",(0,0,700,700))),
 ("Mammalia_early","taxonomic_diversity","BDD","Quintero2024_Fig3A",lambda:pdfclip("quintero_2024",2,(40,290,210,400))),
 ("Placentalia","taxonomic_diversity","","AlvarezCarretero2022_Fig3b",lambda:pdfclip("allvarez-carretero_2022",3,(195,280,420,388))),
]
for cl,pr,dt,tag,src in S:
    sid=f"{cl}_{pr}"+(f"_{dt}" if dt else ""); d=pd.read_csv(DAT+sid+".csv")
    dts=dt if dt and not re.search(r"_\d{4}$",dt) else ""      # reference-derived suffix (e.g. quintero_2024) is already in the tag
    stem=f"{cl}_{SHORT[pr]}"+(f"_{dts}" if dts else "")+f"_{tag}"
    im=src(); fold=OUTD+FOLD[pr]+"/"; os.makedirs(fold,exist_ok=True); im.save(fold+stem+".png")
    fig,ax=plt.subplots(1,2,figsize=(13,4.8),gridspec_kw=dict(width_ratios=[1.15,1]))
    ax[0].imshow(im); ax[0].axis("off"); ax[0].set_title("source: "+tag,fontsize=9)
    has=(d.Err_hi-d.Err_lo).abs().sum()>0
    if has: ax[1].fill_between(d.Age_Ma,d.Err_lo,d.Err_hi,color="0.8",label=d.Err_type.iloc[0])
    ax[1].plot(d.Age_Ma,d.Value,"-o",ms=2.5,lw=1,color="#c0392b",label="digitised")
    ax[1].invert_xaxis(); ax[1].set_xlabel("Age (Ma)"); ax[1].set_ylabel(sid.replace(cl+"_",""),fontsize=8)
    ax[1].set_title(f"{sid}: n={len(d)}, {d.Age_Ma.max():.1f}-{d.Age_Ma.min():.1f} Ma",fontsize=9); ax[1].legend(fontsize=7)
    fig.tight_layout(); fig.savefig(fold+stem+"_QC_digitised.png",dpi=110); plt.close(fig)
    print(fold.split("/")[-2],stem)
