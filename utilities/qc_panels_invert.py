"""qc_panels_invert.py -- step 40 of the digitisation protocol [qc]
Method  : source crops + QC panels (qc_panel) in data_raw/<parameter folder>/
Note    : source crops + QC panels of the invertebrate figure series
Run     : python3 R/digitisation/utilities/qc_panels_invert.py from the Runs root, or python3 R/digitisation/digitise_master.py --run qc_panels_invert
Helpers : shared functions imported from R/digitisation/digitise_master.py (documented there)
(formerly R/digitise_invert/crops.py)
"""
import os, sys; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))   # R/digitisation/, where digitise_master.py lives
from digitise_master import *                                   # shared helpers (see digitise_master.py)
from digitise_master import save_invert as save, OUT_INVERT as OUT, IMG_INVERT as IMG   # invertebrate dialect (4 decimals)
from PIL import Image
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
F="data_raw"   # crops + QC panels go to data_raw/<parameter folder>/ (set per series below)
def vec(pdf,page,clip,name,dpi=400):
    pg=pymupdf.open(PDF+pdf)[page]; pix=pg.get_pixmap(dpi=dpi,clip=pymupdf.Rect(*clip)); pix.save(f"{F}/{name}")
def ras(src,box,name):
    Image.open(src).convert("RGB").crop(box).save(f"{F}/{name}")
S=[("phenotypic_disparity","Bivalvia","Bivalvia_disparity_Zhou2023_Fig1a.png",lambda n:vec("zhou_et_al_2023.pdf",1,(60,45,225,205),n)),
   ("taxonomic_diversity","Bivalvia","Bivalvia_diversity_Zhou2023_Fig1a.png",lambda n:vec("zhou_et_al_2023.pdf",1,(60,45,225,205),n)),
   ("phenotypic_disparity","Pancrustacea","Pancrustacea_disparity_Wills1998_Fig7A.png",lambda n:ras(IMG+"wills-001.png",(0,0,595,475),n)),
   ("taxonomic_diversity","Pancrustacea","Pancrustacea_diversity_Fan2020SI_FigS10.png",lambda n:ras(IMG+"fan-000.png",(1060,1030,2011,2740),n)),
   ("phenotypic_rate","Trilobita","Trilobita_rates_Paterson2019_Fig3.png",lambda n:vec("paterson_et_al_2019_withSI.pdf",3,(125,45,465,240),n)),
   ("phenotypic_disparity","Trilobita","Trilobita_disparity_Foote1993_Fig5A.png",lambda n:ras(IMG+"foote-000.png",(2100,480,3850,2950),n)),
   ("taxonomic_diversity","Trilobita","Trilobita_diversity_Fan2020SI_FigS10.png",lambda n:ras(IMG+"fan-000.png",(30,515,980,2740),n)),
   ("phenotypic_rate","Echinodermata","Echinodermata_rates_NovackGottshall2022_Fig1b.png",lambda n:vec("novack-gottshall_et_al_2022.pdf",1,(300,55,525,300),n)),
   ("phenotypic_disparity","Echinodermata","Echinodermata_disparity_NovackGottshall2022_Fig1c.png",lambda n:vec("novack-gottshall_et_al_2022.pdf",1,(75,255,300,500),n)),
   ("phenotypic_rate","Echinoidea_Irregularia","Echinoidea_Irregularia_rates_HopkinsSmith2015_Fig2A.png",lambda n:vec("hopkins_smith_2015.pdf",2,(100,48,215,145),n,600)),
   ("phenotypic_disparity","Echinoidea_Irregularia","Echinoidea_Irregularia_disparity_Boivin2018_Fig3bottom.png",lambda n:vec("boivin-et_al_2018.pdf",7,(200,280,560,555),n)),
   ("taxonomic_diversity","Echinoidea_Irregularia","Echinoidea_Irregularia_diversity_Boivin2018_Fig3top.png",lambda n:vec("boivin-et_al_2018.pdf",7,(200,65,560,290),n))]
import json; manifest=[]
for par,cl,name,fn in S:
    folder={"phenotypic_disparity":"phenotypic_disparity","taxonomic_diversity":"taxonomic_diversity"}.get(par,"phenotypic_rates")
    F=os.path.join("data_raw",folder); os.makedirs(F,exist_ok=True)
    fn(name)
    qc_panel(cl, par, f"{F}/{name}", f"{OUT}/{cl}_{par}.csv", title_chars=80)   # source | digitised (digitise_master.qc_panel)
    folder={"phenotypic_disparity":"phenotypic_disparity","taxonomic_diversity":"taxonomic_diversity"}.get(par,"phenotypic_rates")
    manifest.append((folder,name)); manifest.append((folder,name.replace(".png","_QC_digitised.png")))
json.dump(manifest,open("data_digitized/digitisation/invert_crops_manifest.json","w")); print(len(manifest))
