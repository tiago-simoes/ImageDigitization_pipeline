"""hopkins_smith2015.py -- step 27 of the digitisation pipeline: hopkins_smith_2015, Fig. 1, 2A [vector]
Method  : black mean markers and vertical CI bars read from vector paths
Outputs : outputs/data_tables/clade_origin/Echinoidea_Irregularia_phenotypic_rate.csv
Note    : Fig. 2A (Irregularia) used; Fig. 1 (all crown-group echinoids) written to outputs/intermediates/ for reference
Run     : python code/source_papers/hopkins_smith2015.py (from any folder), or python code/digitise_master.py --run hopkins_smith2015
Helpers : shared functions imported from code/digitise_master.py (documented there)
"""
import os, sys; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))   # code/, where digitise_master.py lives
from digitise_master import *                                   # shared helpers (see digitise_master.py)
from digitise_master import save_invert as save, OUT_INVERT as OUT, IMG_INVERT as IMG   # invertebrate writer (4 decimals)
D=pymupdf.open(PDF+"hopkins_smith_2015.pdf")
def ext(page,box,xt,xv,lw,ywidth_tick):
    dr=D[page].get_drawings()
    yt=sorted([ (g['rect'].y0) for g in dr if g['type']=='s' and abs(g['width']-lw)<0.02 and abs(g['rect'].width-ywidth_tick)<0.15 and g['rect'].height<0.05 and g['rect'].x1<box[0]+3 and box[1]<g['rect'].y0<box[3]],reverse=True)
    return dr,yt
def points(dr,box,lw):
    C=[]
    for g in dr:
        r=g['rect']
        if g['type']=='fs' and g['fill']==(0,0,0) and all(i[0]=='c' for i in g['items']) and box[0]<=r.x0 and r.x1<=box[2] and box[1]<=r.y0 and r.y1<=box[3]:
            C.append(((r.x0+r.x1)/2,(r.y0+r.y1)/2))
    C=sorted(C); U=[]
    for q in C:
        if not U or abs(U[-1][0]-q[0])>0.3: U.append(q)
    C=U
    V=[g['rect'] for g in dr if g['type']=='s' and abs(g['width']-lw)<0.02 and g['rect'].width<0.05 and g['rect'].height>0.3 and box[0]<=g['rect'].x0<=box[2] and box[1]<=g['rect'].y0<=box[3]]
    rows=[]
    for x,y in C:
        segs=[v for v in V if abs(v.x0-x)<0.3 and v.y0-0.5<=y<=v.y1+0.5]
        top=min([v.y0 for v in segs],default=y); bot=max([v.y1 for v in segs],default=y)
        rows.append((x,y,top,bot))
    return rows
# Fig 1 (crown-group echinoids), page 2
dr=D[1].get_drawings()
X,sx=lin([94.63,120.43,146.23,172.02,197.82,223.62,249.41],[300,250,200,150,100,50,0])
Y,sy=lin([614.95,595.74,576.54,557.34,538.13,518.93,499.73],[0,.1,.2,.3,.4,.5,.6])
rows=points(dr,(92,495,252,615),0.88)
d=pd.DataFrame(rows,columns=["px","py","top","bot"]); d["Age_Ma"]=X(d.px).round(2); d["Value"]=Y(d.py); d["Err_hi"]=Y(d.top); d["Err_lo"]=Y(d.bot); d["Value_err"]=(d.Err_hi-d.Err_lo)/3.92
print("Fig1",len(d),sx,sy); print(d[["Age_Ma","Value","Err_lo","Err_hi"]].round(3).T.to_string())
save(d,"Echinoidea_Irregularia","phenotypic_rate_Fig1_allEchinoidea",sx,"95CI","hopkins_smith_2015 Fig. 1 (crown-group Echinoidea; black = mean)",out=WORK)   # reference only
# Fig 2A (Irregularia), page 3
dr=D[2].get_drawings()
X,sx=lin([132.18,143.71,155.23,166.76,178.29,189.81,201.34],[300,250,200,150,100,50,0])
ylab=[(119.8+126.4)/2,(107.6+114.3)/2,(95.4+102.1)/2,(82.7+89.3)/2,(69.9+76.6)/2]
yt=sorted([g['rect'].y0 for g in dr if g['type']=='s' and g['rect'].height<0.05 and 0.8<g['rect'].width<2 and 124<g['rect'].x0<130 and 60<g['rect'].y0<130],reverse=True)
print("ytick candidates",[round(v,2) for v in yt])
Y,sy=lin(yt if len(yt)==5 else ylab,[0,.5,1,1.5,2])
rows=points(dr,(128,60,206,128),0.35)
d=pd.DataFrame(rows,columns=["px","py","top","bot"]); d["Age_Ma"]=X(d.px).round(2); d["Value"]=Y(d.py); d["Err_hi"]=Y(d.top); d["Err_lo"]=Y(d.bot); d["Value_err"]=(d.Err_hi-d.Err_lo)/3.92
print("Fig2A",len(d),sx,sy); print(d[["Age_Ma","Value","Err_lo","Err_hi"]].round(3).T.to_string())
save(d,"Echinoidea_Irregularia","phenotypic_rate",sx,"95CI","hopkins_smith_2015 Fig. 2A (Irregularia; black = mean)")
