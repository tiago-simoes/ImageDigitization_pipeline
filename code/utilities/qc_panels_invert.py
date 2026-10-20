"""qc_panels_invert.py -- step 40 of the digitisation pipeline [qc]
Method  : source crops + QC panels (qc_panel) in outputs/figures/<parameter folder>/
Note    : source crops + QC panels of the invertebrate figure series; a panel whose source raster has not been
          extracted yet (its script was not part of the run) is skipped with a note
Run     : python code/utilities/qc_panels_invert.py (from any folder), or python code/digitise_master.py --run qc_panels_invert
Helpers : shared functions imported from code/digitise_master.py (documented there)
"""
import os, sys; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))   # code/, where digitise_master.py lives
from digitise_master import *                                   # shared helpers (see digitise_master.py)
from digitise_master import save_invert as save, OUT_INVERT as OUT, IMG_INVERT as IMG   # invertebrate writer (4 decimals)
from PIL import Image
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
F=FIGS   # crops + QC panels go to outputs/figures/<parameter folder>/ (set per series below)
def vec(pdf,page,clip,name,dpi=400):
    pg=pymupdf.open(PDF+pdf)[page]; pix=pg.get_pixmap(dpi=dpi,clip=pymupdf.Rect(*clip)); pix.save(f"{F}/{name}")
def ras(src,box,name):
    Image.open(src).convert("RGB").crop(box).save(f"{F}/{name}")
S=[("phenotypic_disparity","Bivalvia","Bivalvia_disparity_Zhou2023_Fig1a.png",lambda n:vec("zhou_et_al_2023.pdf",1,(60,45,225,205),n)),
   ("taxonomic_diversity","Bivalvia","Bivalvia_diversity_Zhou2023_Fig1a.png",lambda n:vec("zhou_et_al_2023.pdf",1,(60,45,225,205),n)),
   ("phenotypic_disparity","Pancrustacea","Pancrustacea_disparity_Wills1998_Fig7A.png",lambda n:ras(IMG+"wills_1998_p26.png",(0,0,595,475),n)),
   ("taxonomic_diversity","Pancrustacea","Pancrustacea_diversity_Fan2020SI_FigS10.png",lambda n:ras(FIGIMG+"fan_et_al_2020_SI_FigS10.png",(1060,1030,2011,2740),n)),
   ("phenotypic_rate","Trilobita","Trilobita_rates_Paterson2019_Fig3.png",lambda n:vec("paterson_et_al_2019_withSI.pdf",3,(125,45,465,240),n)),
   ("phenotypic_disparity","Trilobita","Trilobita_disparity_Foote1993_Fig5A.png",lambda n:ras(IMG+"foote_1993_p11.png",(2100,480,3850,2950),n)),
   ("taxonomic_diversity","Trilobita","Trilobita_diversity_Fan2020SI_FigS10.png",lambda n:ras(FIGIMG+"fan_et_al_2020_SI_FigS10.png",(30,515,980,2740),n)),
   ("phenotypic_rate","Echinodermata","Echinodermata_rates_NovackGottshall2022_Fig1b.png",lambda n:vec("novack-gottshall_et_al_2022.pdf",1,(300,55,525,300),n)),
   ("phenotypic_disparity","Echinodermata","Echinodermata_disparity_NovackGottshall2022_Fig1c.png",lambda n:vec("novack-gottshall_et_al_2022.pdf",1,(75,255,300,500),n)),
   ("taxonomic_diversity","Echinodermata","Echinodermata_diversity_NovackGottshall2022_EDFig1a.png",lambda n:ras(IMG+"novack_edfig1.jpeg",(950,0,1772,960),n)),
   ("phenotypic_rate","Echinoidea_Irregularia","Echinoidea_Irregularia_rates_HopkinsSmith2015_Fig2A.png",lambda n:vec("hopkins_smith_2015.pdf",2,(100,48,215,145),n,600)),
   ("phenotypic_disparity","Echinoidea_Irregularia","Echinoidea_Irregularia_disparity_Boivin2018_Fig3bottom.png",lambda n:vec("boivin-et_al_2018.pdf",7,(200,280,560,555),n)),
   ("taxonomic_diversity","Echinoidea_Irregularia","Echinoidea_Irregularia_diversity_Boivin2018_Fig3top.png",lambda n:vec("boivin-et_al_2018.pdf",7,(200,65,560,290),n))]
for par,cl,name,fn in S:
    folder={"phenotypic_disparity":"phenotypic_disparity","taxonomic_diversity":"taxonomic_diversity"}.get(par,"phenotypic_rates")
    F=os.path.join(FIGS,folder); os.makedirs(F,exist_ok=True)
    try: fn(name)
    except FileNotFoundError as e:   # raster not extracted yet (its script was not part of this run, e.g. --key with a subset)
        print(f"skipped {folder}/{name}: {os.path.relpath(e.filename, ROOT) if e.filename else 'source image'} not found (run the script that extracts it)"); continue
    qc_panel(cl, par, f"{F}/{name}", f"{OUT}/{cl}_{par}.csv", title_chars=80)   # source | digitised (digitise_master.qc_panel)
    print(f"{folder}/{name}")
