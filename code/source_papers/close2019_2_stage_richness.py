"""close2019_2_stage_richness.py -- step 3 of the digitisation protocol: close_2019, Fig. 1a [derived]
Method  : summed locality richness per ICS stage (ICS chart v2023/09)
Outputs : data_digitized/Tetrapoda_nonflying_taxonomic_diversity.csv
Note    : summed richness per ICS stage
Run     : python3 R/digitisation/source_papers/close2019_2_stage_richness.py from the Runs root, or python3 R/digitisation/digitise_master.py --run close2019_2_stage_richness
(formerly R/digitise_close2019/stages.py)
"""
PDF="data_raw/articles_PDF/close_2019.pdf"   # run from the Runs root: python3 R/digitisation/source_papers/close2019_<n>_*.py
DIG="data_digitized/digitisation/"; DAT="data_digitized/"; QC="data_raw/taxonomic_diversity/"
WORK="data_digitized/digitisation/close2019_work/"   # intermediates; f5/img-000.jpg = embedded Fig. 5 raster of close_2019.pdf p. 4 (pdfimages -j -f 4 -l 4 ... f5/img)
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
# ICS International Chronostratigraphic Chart v2023/09: (stage, period, base Ma)
S=[("Frasnian","Devonian",382.7),("Famennian","Devonian",372.2),
("Tournaisian","Carboniferous",358.9),("Visean","Carboniferous",346.7),("Serpukhovian","Carboniferous",330.9),("Bashkirian","Carboniferous",323.2),("Moscovian","Carboniferous",315.2),("Kasimovian","Carboniferous",307.0),("Gzhelian","Carboniferous",303.7),
("Asselian","Permian",298.9),("Sakmarian","Permian",293.52),("Artinskian","Permian",290.1),("Kungurian","Permian",283.5),("Roadian","Permian",273.01),("Wordian","Permian",266.9),("Capitanian","Permian",264.28),("Wuchiapingian","Permian",259.51),("Changhsingian","Permian",254.14),
("Induan","Triassic",251.902),("Olenekian","Triassic",251.2),("Anisian","Triassic",247.2),("Ladinian","Triassic",242.0),("Carnian","Triassic",237.0),("Norian","Triassic",227.0),("Rhaetian","Triassic",208.5),
("Hettangian","Jurassic",201.4),("Sinemurian","Jurassic",199.5),("Pliensbachian","Jurassic",192.9),("Toarcian","Jurassic",184.2),("Aalenian","Jurassic",174.7),("Bajocian","Jurassic",170.9),("Bathonian","Jurassic",168.2),("Callovian","Jurassic",165.3),("Oxfordian","Jurassic",161.5),("Kimmeridgian","Jurassic",154.8),("Tithonian","Jurassic",149.2),
("Berriasian","Cretaceous",145.0),("Valanginian","Cretaceous",139.8),("Hauterivian","Cretaceous",132.6),("Barremian","Cretaceous",125.77),("Aptian","Cretaceous",121.4),("Albian","Cretaceous",113.0),("Cenomanian","Cretaceous",100.5),("Turonian","Cretaceous",93.9),("Coniacian","Cretaceous",89.8),("Santonian","Cretaceous",86.3),("Campanian","Cretaceous",83.6),("Maastrichtian","Cretaceous",72.1),
("Danian","Paleogene",66.0),("Selandian","Paleogene",61.6),("Thanetian","Paleogene",59.2),("Ypresian","Paleogene",56.0),("Lutetian","Paleogene",48.07),("Bartonian","Paleogene",41.2),("Priabonian","Paleogene",37.71),("Rupelian","Paleogene",33.9),("Chattian","Paleogene",27.82),
("Aquitanian","Neogene",23.03),("Burdigalian","Neogene",20.44),("Langhian","Neogene",15.98),("Serravallian","Neogene",13.82),("Tortonian","Neogene",11.63),("Messinian","Neogene",7.246),("Zanclean","Neogene",5.333),("Piacenzian","Neogene",3.6),
("Gelasian","Quaternary",2.58),("Calabrian","Quaternary",1.8),("Chibanian","Quaternary",0.774),("Upper Pleistocene","Quaternary",0.129),("Holocene","Quaternary",0.0117)]
st=pd.DataFrame(S,columns=["Stage","Period","Base_Ma"]); st["Top_Ma"]=list(st.Base_Ma[1:])+[0.0]; st["Mid_Ma"]=(st.Base_Ma+st.Top_Ma)/2
d=pd.read_csv(DIG+"close2019_fig1a_localities.csv")
d["Age_Ma"]=d.Age_Ma.clip(lower=0)
idx=np.searchsorted(-st.Base_Ma.values,-d.Age_Ma.values,side="right")-1   # stage whose [top,base) contains age
d["Stage"]=st.Stage.values[np.clip(idx,0,len(st)-1)]
g=d.groupby("Stage").agg(N_localities=("Richness","size"),Summed_richness=("Richness","sum"),Max_local_richness=("Richness","max"),Median_local_richness=("Richness","median"))
out=st.merge(g,on="Stage",how="left").fillna({"N_localities":0,"Summed_richness":0})
out=out[out.Base_Ma<=d.Age_Ma.max()+30]   # from the first sampled stage onward
for c in ["N_localities","Summed_richness"]: out[c]=out[c].astype(int)
out.round(3).to_csv(DIG+"Tetrapoda_nonflying_stage_summed_richness.csv",index=False)
print(out[["Stage","Period","Base_Ma","N_localities","Summed_richness","Max_local_richness"]].to_string(index=False))
# tidy series for the pipeline (stage midpoints; Poisson error on counts)
t=out[out.N_localities>0].copy()
tidy=pd.DataFrame({"Clade":"Tetrapoda_nonflying","Parameter":"taxonomic_diversity","Age_Ma":t.Mid_Ma.round(3),"Value":t.Summed_richness,
  "Value_err":np.sqrt(t.Summed_richness).round(3),"Err_lo":(t.Summed_richness-np.sqrt(t.Summed_richness)).round(3),"Err_hi":(t.Summed_richness+np.sqrt(t.Summed_richness)).round(3),
  "Err_type":"Poisson","Source":"close_2019 Fig. 1a (sum of per-locality species richness per ICS stage; vector-extracted, n=24360 localities)","Age_cal_sd":0.0})
tidy.to_csv(DAT+"Tetrapoda_nonflying_taxonomic_diversity.csv",index=False)
