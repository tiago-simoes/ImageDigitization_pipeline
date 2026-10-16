"""extract_raster_figures.py -- step 1 of the digitisation pipeline: embedded raster figures [utility]
Method  : saves the embedded raster images of the source figures (page + image index) for the tracing scripts
Outputs : outputs/intermediates/vert_img/*.png (vertebrate figures), outputs/intermediates/invert_img/*.png
          (Wills 1998, Foote 1993), outputs/intermediates/close2019/close2019_fig5.jpg (Close et al. 2019 Fig. 5)
Note    : run first. The rasters of the mass-extinction sources are extracted by their own scripts; the Fan et al. 2020
          SI figure is supplied in articles_PDF/figure_images/ (see README.md)
Run     : python code/utilities/extract_raster_figures.py (from any folder), or python code/digitise_master.py --run extract_raster_figures
"""
import os, sys; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))   # code/, where digitise_master.py lives
from digitise_master import PDF as A, IMG_VERT, IMG_INVERT, CLOSE2019, pymupdf      # paths (digitise_master.py, section 1)

def save_xref(pdf, xref, path):
    """Decode an embedded image (by xref) and save it as PNG (alpha dropped, CMYK converted to RGB)."""
    d = pymupdf.open(A + pdf + ".pdf"); pix = pymupdf.Pixmap(d, xref)
    if pix.alpha: pix = pymupdf.Pixmap(pix, 0)
    if pix.n > 3: pix = pymupdf.Pixmap(pymupdf.csRGB, pix)
    pix.save(path); print(os.path.basename(path), pix.width, pix.height)

def xref_of(pdf, page, index):
    """xref of the n-th embedded image on a page (0-based page and index; looked up, robust to re-saved PDFs)."""
    return pymupdf.open(A + pdf + ".pdf")[page].get_image_info(xrefs=True)[index]["xref"]

def save_page_image(pdf, page, index, name, folder=IMG_VERT):
    save_xref(pdf, xref_of(pdf, page, index), folder + name + ".png")

def save_page_image_raw(pdf, page, index, path):
    """Save the image stream as stored in the PDF (e.g. the original JPEG), without decoding."""
    d = pymupdf.open(A + pdf + ".pdf"); b = d.extract_image(xref_of(pdf, page, index))
    open(path, "wb").write(b["image"]); print(os.path.basename(path), b["width"], b["height"])

# vertebrate sources
save_page_image("foth_2021", 11, 0, "foth9A_r1"); save_page_image("foth_2021", 11, 3, "foth9A_r2"); save_page_image("foth_2021", 11, 6, "foth9A_r3")
save_page_image("sookias_2026", 3, 0, "sookias_fig1")
save_page_image("quintero_2025_SI", 20, 0, "quintero25_S9")
save_page_image("rawson_2022", 5, 0, "rawson_fig5")
save_page_image("crouch_tobias_2022", 5, 0, "crouch_fig3")
save_page_image("felice_goswami_2017", 1, 1, "felice_fig1")
save_page_image("close_2015", 2, 0, "close15_fig1")
save_page_image("goswami_2022", 3, 2, "goswami_fig3")
save_page_image("allvarez-carretero_2022", 3, 39, "alvarez_fig3b")   # Fig. 3b inset (index of the inset raster on page 4)
# invertebrate sources (scanned pages)
save_page_image("wills_1998", 25, 1, "wills_1998_p26", IMG_INVERT)                            # Fig. 7 (second image of PDF p. 26)
save_page_image("foote_1993", 10, 0, "foote_1993_p11", IMG_INVERT)                            # Fig. 5 (PDF p. 11)
# Close et al. 2019 Fig. 5 (PDF p. 4), kept as the stored JPEG
save_page_image_raw("close_2019", 3, 0, CLOSE2019 + "close2019_fig5.jpg")
