"""close2019_5_squamata_stages.py -- step 6 of the digitisation protocol: close_2019, Fig. 5c [derived]
Method  : highest dot top per ICS stage = maximum local richness of squamates (+/-1 species reading error)
Outputs : data_digitized/Squamata_stage_max_local_richness.csv, data_digitized/Squamata_max_local_richness.csv
Note    : max local richness of squamates per stage
Run     : python3 R/digitisation/source_papers/close2019_5_squamata_stages.py from the Runs root, or python3 R/digitisation/digitise_master.py --run close2019_5_squamata_stages
(formerly R/digitise_close2019/sq_stage.py)
"""
import numpy as np, pandas as pd
from PIL import Image
PDF="data_raw/articles_PDF/close_2019.pdf"   # run from the Runs root: python3 R/digitisation/source_papers/close2019_<n>_*.py
DIG="data_digitized/digitisation/"; DAT="data_digitized/"; QC="data_raw/taxonomic_diversity/"
WORK="data_digitized/digitisation/close2019_work/"   # intermediates; f5/img-000.jpg = embedded Fig. 5 raster of close_2019.pdf p. 4 (pdfimages -j -f 4 -l 4 ... f5/img)
import os; os.makedirs(WORK+"f5",exist_ok=True)
green=np.load(WORK+"f5/green_mask.npy")
X0,Y0=330.90655517578125,64.2503662109375; SX=895/(545.5865478515625-X0); SY=1297/(375.5303649902344-Y0)
px=lambda v:(v-X0)*SX; py=lambda v:(v-Y0)*SY
bx,ax_=np.polyfit([px(v) for v in [333.95,386.47,439.03,491.59,544.2]],[400,300,200,100,0],1)
by,ay=np.polyfit([py(v) for v in [360.41,337.63,314.86,292.09,269.32]],[0,20,40,60,80],1)
r=10.5/2
H,Wd=green.shape; ytop=np.full(Wd,np.inf)
for x in range(Wd):
    ys=np.where(green[:,x])[0]
    if len(ys): ytop[x]=ys.min()
# dot tops = local minima of the top profile (window +-2 px); plateaus collapsed to their middle
cand=[x for x in range(2,Wd-2) if np.isfinite(ytop[x]) and ytop[x]<=ytop[x-2:x+3].min()]
runs=[];
for x in cand:
    if runs and x-runs[-1][-1]<=1 and ytop[x]==ytop[runs[-1][-1]]: runs[-1].append(x)
    else: runs.append([x])
tops=pd.DataFrame([(np.mean(rn),ytop[rn[0]]) for rn in runs],columns=["px","ytop"])
tops["Age_Ma"]=ax_+bx*tops.px; tops["Richness"]=np.clip(np.round(ay+by*(tops.ytop+r)),1,None).astype(int)
tops.round(3).to_csv(DIG+"close2019_fig5c_squamata_dot_tops.csv",index=False)
st=pd.read_csv(DAT+"Tetrapoda_nonflying_stage_max_local_richness.csv")[["Stage","Period","Base_Ma","Top_Ma","Mid_Ma"]]
age_col=ax_+bx*np.arange(Wd)
rows=[]
for _,s in st.iterrows():
    t=tops[(tops.Age_Ma<=s.Base_Ma)&(tops.Age_Ma>s.Top_Ma)]
    cols=(age_col<=s.Base_Ma)&(age_col>s.Top_Ma)
    has_green=bool(np.isfinite(ytop[cols]).any())
    width_px=cols.sum()
    note=""
    if len(t): mx=int(t.Richness.max())
    elif has_green:
        mx=int(np.clip(np.round(ay+by*(ytop[cols].min()+r)),1,None)); note="no dot top inside stage (covered by neighbouring dots); value read from envelope"
    else: mx=np.nan
    if width_px< 10.5 and has_green: note=(note+"; " if note else "")+f"stage narrower than one dot ({width_px} px = {s.Base_Ma-s.Top_Ma:.2f} Myr): assignment uncertain"
    rows.append(dict(Stage=s.Stage,N_localities=np.nan,Max_local_richness=mx,Q95_local_richness=np.nan,Max_boot_lo=np.nan,Max_boot_hi=np.nan,
                     RarefiedMax_n20=np.nan,RarefiedMax_n20_lo=np.nan,RarefiedMax_n20_hi=np.nan,Envelope_dot_tops_in_stage=len(t),Note=note))
o=st.merge(pd.DataFrame(rows),on="Stage")
first=o.Max_local_richness.first_valid_index(); o=o.loc[first:]   # from first stage with squamates
o.to_csv(DAT+"Squamata_stage_max_local_richness.csv",index=False)
print(o[["Stage","Base_Ma","Max_local_richness","Envelope_dot_tops_in_stage","Note"]].to_string(index=False))
t=o.dropna(subset=["Max_local_richness"])
pd.DataFrame({"Clade":"Squamata","Parameter":"max_local_richness","Age_Ma":t.Mid_Ma.round(3),"Value":t.Max_local_richness,
  "Value_err":0.5,"Err_lo":t.Max_local_richness-1,"Err_hi":t.Max_local_richness+1,"Err_type":"none",
  "Source":"close_2019 Fig. 5c (squamates; max species richness per locality per ICS stage; raster upper-envelope extraction, +/-1 species reading error)","Age_cal_sd":2.4}).to_csv(DAT+"Squamata_max_local_richness.csv",index=False)
