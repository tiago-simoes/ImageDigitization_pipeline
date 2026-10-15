"""extract_raster_figures.py -- step 1 of the digitisation protocol: foth_2021, sookias_2026, quintero_2025_SI, rawson_2022, crouch_tobias_2022, felice_goswami_2017, close_2015, goswami_2022, caron_pie_2024, allvarez-carretero_2022, embedded rasters [utility]
Method  : saves the embedded raster images of the vertebrate source figures (page + image index) as PNGs for the tracing scripts
Note    : Run first: extracts the raster figures of the vertebrate sources to data_digitized/digitisation/vert_img/
Run     : python3 R/digitisation/utilities/extract_raster_figures.py from the Runs root, or python3 R/digitisation/digitise_master.py --run extract_raster_figures
(formerly R/digitise_vert/extract_images.py)

Run first: pull the embedded raster figures out of the source PDFs (run from the Runs root).
Output: data_digitized/digitisation/vert_img/*.png, used by the raster tracing scripts.
"""
import pymupdf, os
A="data_raw/articles_PDF/"; O="data_digitized/digitisation/vert_img/"; os.makedirs(O,exist_ok=True)
def save_xref(pdf,xref,name):
    d=pymupdf.open(A+pdf+".pdf"); pix=pymupdf.Pixmap(d,xref)
    if pix.alpha: pix=pymupdf.Pixmap(pix,0)
    if pix.n>3: pix=pymupdf.Pixmap(pymupdf.csRGB,pix)
    pix.save(O+name+".png"); print(name,pix.width,pix.height)
def save_page_image(pdf,page,index,name):        # n-th embedded image on a page (xref looked up, robust to re-saves)
    p=pymupdf.open(A+pdf+".pdf")[page]; save_xref(pdf,p.get_image_info(xrefs=True)[index]["xref"],name)
save_page_image("foth_2021",11,0,"foth9A_r1"); save_page_image("foth_2021",11,3,"foth9A_r2"); save_page_image("foth_2021",11,6,"foth9A_r3")
save_page_image("sookias_2026",3,0,"sookias_fig1")
save_page_image("quintero_2025_SI",20,0,"quintero25_S9")
save_page_image("rawson_2022",5,0,"rawson_fig5")
save_page_image("crouch_tobias_2022",5,0,"crouch_fig3")
save_page_image("felice_goswami_2017",1,1,"felice_fig1")
save_page_image("close_2015",2,0,"close15_fig1")
save_page_image("goswami_2022",3,2,"goswami_fig3")
save_page_image("caron_pie_2024",6,0,"caron_fig3")
save_page_image("allvarez-carretero_2022",3,39,"alvarez_fig3b")   # Fig. 3b inset (index of the inset raster on page 4)
