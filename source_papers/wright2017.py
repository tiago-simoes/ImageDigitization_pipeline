"""wright2017.py -- step 32 of the digitisation protocol: Wright_2017_SI, Figs S6, S8, S10 [raster]
Method  : markers and error-bar ends detected in the source crops (morphological opening, peak detection); calibration on all ticks
Outputs : data_digitized/Crinoidea_Eucladida_phenotypic_rate.csv, data_digitized/Crinoidea_Eucladida_phenotypic_disparity.csv, data_digitized/Crinoidea_Eucladida_taxonomic_diversity.csv
Note    : reads the source crops in data_raw/; writes calibration_ticks.csv and metaDigitise_style_output.csv
Run     : python3 R/digitisation/source_papers/wright2017.py from the Runs root, or python3 R/digitisation/digitise_master.py --run wright2017
(formerly R/digitise_wright2017_SI.py)

Automated metaDigitise-equivalent extraction for Wright (2017) SI Figs S6, S8, S10.
Calibration: least-squares linear fit of pixel->value on ALL labelled axis ticks
(metaDigitise uses 2 clicked calibration points; using all ticks is more robust and
gives a calibration residual). Points: marker centroids; error bars: bar-end pixels.
"""
import numpy as np, pandas as pd, json, os
from PIL import Image, ImageDraw
from scipy import ndimage as ndi
from skimage.morphology import disk, opening
from skimage.feature import peak_local_max

def runs(idx,gap=3):
    idx=list(idx); out=[]
    s=p=idx[0]
    for i in idx[1:]:
        if i-p>gap: out.append((s+p)/2); s=i
        p=i
    out.append((s+p)/2); return np.array(out)

def ticks(blk,xr,yc):
    xs=np.where(blk[xr[0]:xr[1]].sum(0)>=(xr[1]-xr[0])*0.6)[0]
    ys=np.where(blk[:,yc[0]:yc[1]].sum(1)>=(yc[1]-yc[0])*0.6)[0]
    return runs(xs),runs(ys)

def calib(pix,val):
    b,a=np.polyfit(pix,val,1); res=val-(a+b*pix)
    return a,b,res

def bar_end(blk,cx,cy,direction,r0,halfw=2,gap=3,stop=None):
    """walk from marker centre along column until black run ends"""
    y=cy; last=cy; miss=0; H=blk.shape[0]
    # skip a non-black marker fill (S6 colored circles) until the outline is hit
    k=0
    while k<25 and not blk[y,int(cx)-halfw:int(cx)+halfw+1].any(): y+=direction; k+=1
    last=y
    while 0<y<H-1:
        y+=direction
        if stop is not None and ((direction<0 and y<stop) or (direction>0 and y>stop)): break
        if blk[y,int(cx)-halfw:int(cx)+halfw+1].any(): last=y; miss=0
        else:
            miss+=1
            if miss>gap: break
    return last

FIG={
 "rates":dict(file="Crinoidea_Eucladida_rates_Wright2017_FigS6.png",folder="phenotypic_rates",
   fig="Figure S6",param="phenotypic_rate",ylab="Rate (# changes / Lmy)",
   xr=(1424,1445),yc=(225,246),xv=[450,400,350,300,250],yv=[1.2,1.0,0.8,0.6,0.4,0.2],
   err="95CI"),
 "disparity":dict(file="Crinoidea_Eucladida_disparity_Wright2017_FigS8.png",folder="phenotypic_disparity",
   fig="Figure S8",param="phenotypic_disparity",ylab="Disparity (mean pairwise Gower)",
   xr=(1268,1300),yc=(335,372),xv=[450,400,350,300,250],yv=[.7,.6,.5,.4,.3,.2,.1,0],err="SE"),
 "diversity":dict(file="Crinoidea_Eucladida_diversity_Wright2017_FigS10.png",folder="taxonomic_diversity",
   fig="Figure S10",param="taxonomic_diversity",ylab="Genus richness",
   xr=(1422,1450),yc=(340,375),xv=[450,400,350,300,250],yv=[50,40,30,20,10,0],err="SE"),
}
calrows=[]; tidy=[]; md=[]
for k,c in FIG.items():
    im=Image.open(os.path.join("data_raw",c["folder"],c["file"])).convert("RGB"); a=np.asarray(im).astype(int)
    blk=a.max(2)<110
    tx,ty=ticks(blk,c["xr"],c["yc"])
    assert len(tx)==len(c["xv"]) and len(ty)==len(c["yv"]),(k,tx,ty)
    ax_,bx,rx=calib(tx,np.array(c["xv"],float)); ay,by,ry=calib(ty,np.array(c["yv"],float))
    for p,v,r in zip(tx,c["xv"],rx): calrows.append([c["fig"],"x",p,v,r])
    for p,v,r in zip(ty,c["yv"],ry): calrows.append([c["fig"],"y",p,v,r])
    X=lambda p: ax_+bx*p; Y=lambda p: ay+by*p
    cal_sd=(float(np.sqrt(np.mean(rx**2))),float(np.sqrt(np.mean(ry**2))))
    # --- markers
    if k=="rates":
        r,g,b=a[...,0],a[...,1],a[...,2]
        fill=((a.max(2)-a.min(2))>60)|((abs(r-g)<12)&(abs(g-b)<12)&(r>140)&(r<215))
        fill[:, :252]=0; fill[1415:]=0
        lab,n=ndi.label(opening(fill,disk(5)))
        cents=ndi.center_of_mass(lab>0,lab,range(1,n+1))
        pts=[(cx,cy) for cy,cx in cents]
        cls=[ "high" if a[int(cy),int(cx),0]>200 and a[int(cy),int(cx),1]<100 else
              ("low" if a[int(cy),int(cx),2]>200 and a[int(cy),int(cx),0]<100 else "background") for cx,cy in pts]
        stopY=(25,1415)
    elif k=="disparity":
        m=a.max(2)<100; m[:, :390]=0; m[1200:]=0
        lab,n=ndi.label(opening(m,np.ones((31,31),bool)))
        pts=[(cx,cy) for cy,cx in ndi.center_of_mass(lab>0,lab,range(1,n+1))]; cls=[""]*len(pts)
        stopY=(5,1205)
    else:
        m=a.max(2)<100; m[:, :395]=0; m[1355:]=0
        dt=ndi.distance_transform_edt(m)
        pk=peak_local_max(dt,min_distance=12,threshold_abs=14); pts=[]
        for y,x in pk:
            yy,xx=np.mgrid[y-10:y+11,x-10:x+11]; w=dt[yy,xx]>=dt[y,x]-1.5
            pts.append((xx[w].mean(),yy[w].mean()))
        cls=[""]*len(pts); stopY=(100,1355)
    order=np.argsort([p[0] for p in pts]); pts=[pts[i] for i in order]; cls=[cls[i] for i in order]
    dr=ImageDraw.Draw(im)
    for i,((cx,cy),cl) in enumerate(zip(pts,cls)):
        top=bar_end(blk,cx,int(round(cy)),-1,0,stop=stopY[0]); bot=bar_end(blk,cx,int(round(cy)),+1,0,stop=stopY[1])
        lw=1.5 if k=="rates" else 2.5     # half line thickness correction
        top+=lw; bot-=lw
        v=Y(cy); vhi=Y(top); vlo=Y(bot)
        if k=="diversity": v=max(v,0.0)   # genus counts cannot be negative (calibration noise)
        if c["err"]=="95CI":
            se=(vhi-vlo)/(2*1.96); lo,hi=vlo,vhi
        else:
            half=min(vhi-v, v-vlo)   # symmetric 1 SE; take uncontaminated side
            se=half; lo,hi=v-half,v+half
        tidy.append(dict(Clade="Crinoidea_Eucladida",Parameter=c["param"],Age_Ma=round(X(cx),2),
                         Value=round(v,4),Value_err=round(se,4),Err_lo=round(lo,4),Err_hi=round(hi,4),
                         Err_type=c["err"],Rate_class=cl,Source=f"Wright_2017_SI {c['fig']}",
                         px_x=round(cx,1),px_y=round(cy,1),
                         Age_cal_sd=round(cal_sd[0],2)))   # x-axis calibration RMSE (Ma), used for age jitter in the MC
        md.append(dict(filename=c["file"],group_id=f"pt{i+1:02d}",x=round(X(cx),3),y=round(v,4),
                       y_upper=round(hi,4),y_lower=round(lo,4),plot_type="scatterplot",
                       x_variable="Age_Ma",y_variable=c["ylab"],error_type=c["err"]))
        dr.ellipse([cx-8,cy-8,cx+8,cy+8],outline=(0,200,255),width=4)
        dr.line([cx+14,top,cx+14,bot],fill=(255,0,200),width=3)
    for p in tx: dr.line([p,c["xr"][0]-30,p,c["xr"][1]+30],fill=(0,170,0),width=2)
    for p in ty: dr.line([c["yc"][0]-30,p,c["yc"][1]+30,p],fill=(0,170,0),width=2)
    im.save(os.path.join("data_raw",c["folder"],c["file"].replace(".png","_QC_digitised.png")))
    print(k,len(pts),"cal RMSE x(Ma)=%.2f y=%.4f"%cal_sd)
pd.DataFrame(calrows,columns=["figure","axis","pixel","value","residual"]).to_csv("data_digitized/digitisation/calibration_ticks.csv",index=False)
T=pd.DataFrame(tidy); T.to_csv("data_digitized/digitisation/wright2017_all_tidy.csv",index=False)
pd.DataFrame(md).to_csv("data_digitized/digitisation/metaDigitise_style_output.csv",index=False)
for p,g in T.groupby("Parameter"):
    g.to_csv(f"data_digitized/Crinoidea_Eucladida_{p}.csv",index=False)
pd.set_option("display.width",200); print(T.drop(columns=["Source","Clade"]).to_string())
