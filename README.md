# Figure digitisation 

All Python code that extracts data points directly from figures in published papers using PDF files (and raw data files)

```
code/
├── digitise_master.py     master file: shared functions + registry of scripts + runner
├── README.md              this file
├── source_papers/         one script per source paper (37 scripts)
└── utilities/             steps not tied to one paper (5 scripts): raster extraction, QC panels, age cutoffs
```

* **`digitise_master.py`**: the master file. It is
  1. the **shared library**: every helper used by more than one script (axis calibration, PDF vector-path reading,
     raster masks and tracing, regridding, node-rate binning, output writers, QC panels), each with a commented
     explanation. 
  2. the **registry**: one entry per script with source paper, figure, method (vector / raster / raw / derived),
     output series and run order;
  3. the **runner**: list, run and check the scripts from the command line (below).
* **`source_papers/`: one script per source paper**, named after the paper (`guo2024.py`, `foth2021.py`, ...). Each
  script starts with a header giving its step number, source figure, method, outputs and former file name. It then
  adds `code/digitization/` to the Python path and imports the shared functions with `from digitise_master import *`
  (these two lines must stay at the top of any new script). Figure-specific values (tick positions, line colours,
  crop boxes) stay in the script, so each figure's extraction can be read and audited on its own.
* **`utilities/`: the steps that serve all sources**: `extract_raster_figures.py` (runs first: saves the embedded
  raster figures of the vertebrate PDFs), `qc_panels_vert.py`, `qc_panels_invert.py`, `qc_panels_massext.py` (source
  crops + QC panels) and `age_cutoffs.py` (PI age cutoffs; always last). The runner knows which folder each script
  is in (groups `utility`, `qc`, `cutoff` = `utilities/`).

## Quick start (from the `Runs` root)

```
python3 R/digitisation/digitise_master.py --list            # all scripts in run order (add --long for method + outputs)
python3 R/digitisation/digitise_master.py --which Eurypterida_phenotypic_disparity   # which script makes a series
python3 R/digitisation/digitise_master.py --run all         # re-run everything in the correct order
python3 R/digitisation/digitise_master.py --run whalen2020 age_cutoffs              # selected scripts
python3 R/digitisation/digitise_master.py --run group:mass_extinction               # one group
python3 R/digitisation/digitise_master.py --check           # re-run in a scratch copy and compare with data_digitized/
python3 R/digitisation/source_papers/whalen2020.py          # any script can also be run on its own
python3 R/digitisation/utilities/age_cutoffs.py
```

`age_cutoffs.py` applies the PI age cutoffs (`data_digitized/digitisation/age_cutoffs.csv`) and must run **after**
any re-digitised series (`--run all` does this). The QC panels are drawn before the cutoffs, so they show the full
traced series (a `qc_panels_*.py` script re-run on its own after `age_cutoffs.py` draws the cut series instead).
After re-digitising, re-run `R/prep_mass_extinction.py` and the R pipelines (see `R/README.md`).

Requirements: Python 3 with numpy, pandas, scipy, pillow, matplotlib and pymupdf (`pip install pymupdf`);
scikit-image for `wright2017.py`. R is not needed for this step.

## Inputs and outputs

Inputs (the only files `--check` copies into its scratch folder):

* `data_raw/articles_PDF/`: source PDFs and the raw data CSVs listed in the keys;
* `data_raw/<parameter folder>/Crinoidea_Eucladida_*_Wright2017_FigS*.png`: Wright 2017 SI figure crops;
* `data_digitized/digitisation/invert_img/` and `close2019_work/f5/`: rasters extracted earlier with `pdfimages`;
* `data_digitized/digitisation/age_cutoffs.csv`.

Outputs:

* `data_digitized/<Clade>_<Parameter>[_<data_type>].csv`, one row per time point:
  `Clade, Parameter, Age_Ma, Value, Value_err, Err_lo, Err_hi, Err_type, Source, Age_cal_sd`.
  `Err_type` is `95CI`, `SE`, `Poisson` (counts without published error) or `none`; `Age_cal_sd` is the SD of the
  age calibration (Myr), used for the Monte-Carlo age jitter in both pipelines;
* `data_digitized/digitisation/`: intermediates (extracted rasters `*_img/`, tick tables, Close 2019 working files)
  and the per-series notes `notes_{vert,invert,massext}.csv` that are appended to the pipeline Flags;
* `data_raw/<parameter folder>/`: a source crop and a `_QC_digitised.png` panel (source | digitised series) per series.

## Method in brief

* **Vector figures** (curves stored as PDF drawing paths): coordinates read directly with PyMuPDF; exact apart from
  the axis calibration.
* **Raster figures** (embedded images): the image is extracted from the PDF, the curve isolated with colour or
  darkness masks and traced per pixel column (median of the line pixels); markers and error bars are isolated by
  morphological opening.
* **Calibration**: least-squares linear fit of pixel (or PDF point) position on all available tick marks or labels;
  the fit RMSE combined with half a pixel is stored as `Age_cal_sd`.
* Dense traces are averaged on a 1-Myr grid; node-based raw rate files are summarised in 5-Myr bins
  (`node_bins`, also used by `R/prep_mass_extinction.py` with bins anchored at the extinction age).

Each script header and `--list --long` describe the specific extraction; `digitise_master.py` documents every helper.

## Scripts in run order (groups `utility`, `qc`, `cutoff` in `utilities/`; all others in `source_papers/`)

| step | script | group | method | source (key reference) | figure | series written to `data_digitized/` |
|---:|---|---|---|---|---|---|
| 1 | `extract_raster_figures.py` | utility | utility | foth_2021, sookias_2026, quintero_2025_SI, rawson_2022, crouch_tobias_2022, felice_goswami_2017, close_2015, goswami_2022, caron_pie_2024, allvarez-carretero_2022 | embedded rasters | Run first: extracts the raster figures of the vertebrate sources to data_digitized/digitisation/vert_img/ |
| 2 | `close2019_1_extract_localities.py` | close2019 | vector | close_2019 | Fig. 1a | locality richness points -> digitisation/close2019_fig1a_localities.csv |
| 3 | `close2019_2_stage_richness.py` | close2019 | derived | close_2019 | Fig. 1a | Tetrapoda_nonflying_taxonomic_diversity |
| 4 | `close2019_3_max_richness.py` | close2019 | derived | close_2019 | Fig. 1a | Tetrapoda_nonflying_stage_max_local_richness |
| 5 | `close2019_4_squamata_extract.py` | close2019 | raster | close_2019 | Fig. 5c | squamate dot envelope (green mask) |
| 6 | `close2019_5_squamata_stages.py` | close2019 | derived | close_2019 | Fig. 5c | Squamata_stage_max_local_richness, Squamata_max_local_richness |
| 7 | `close2019_6_squamata_plot.py` | close2019 | qc | close_2019 | Fig. 5c | summary plots + QC panel |
| 8 | `ghezelayagh2022_fig2.py` | vertebrates | vector | ghezelayagh_2022 | Fig. 2 | Acanthomorpha_phenotypic_disparity |
| 9 | `ghezelayagh2022_edfig5e.py` | vertebrates | vector | ghezelayagh_2022 | Extended Data Fig. 5e | Acanthomorpha_taxonomic_diversity |
| 10 | `brownstein_griffin2026.py` | vertebrates | vector | brownstein_griffin_2026 | Fig. 2A, 2B | Theropoda_non-avian_phenotypic_rate, Dinosauria_early_phenotypic_rate |
| 11 | `yu2021.py` | vertebrates | vector | yu_2021_SI | Fig. S40 | Aves_stem_taxonomic_diversity, Aves_crown_taxonomic_diversity |
| 12 | `navalon2022.py` | vertebrates | vector | navalon_2022 | Fig. 1a | Aves_crown_phenotypic_disparity |
| 13 | `quintero2024.py` | vertebrates | vector | quintero_2024 | Fig. 3A | Mammalia_early_taxonomic_diversity_BDD |
| 14 | `alvarez_carretero2022.py` | vertebrates | raster | allvarez-carretero_2022 | Fig. 3b | Placentalia_taxonomic_diversity |
| 15 | `foth2021.py` | vertebrates | raster | foth_2021 | Fig. 9A rows 1-3 | Archosauromorpha_stem_phenotypic_disparity, Dinosauria_early_phenotypic_disparity, Theropoda_non-avian_phenotypic_disparity |
| 16 | `sookias2026.py` | vertebrates | raster | sookias_2026 | Fig. 1b, 1f, 1h | Archosauromorpha_stem_taxonomic_diversity, Dinosauria_early_taxonomic_diversity_PyRate, Theropoda_non-avian_taxonomic_diversity |
| 17 | `quintero2025.py` | vertebrates | raster | quintero_2025_SI | Suppl. Fig. 9 | Dinosauria_early_taxonomic_diversity_FBDD |
| 18 | `rawson2022.py` | vertebrates | raster | rawson_2022 | Fig. 5A | Tetrapoda_stem_phenotypic_disparity |
| 19 | `crouch_tobias2022.py` | vertebrates | raster | crouch_tobias_2022 | Fig. 3A | Aves_crown_phenotypic_rate_continuous_data_LMM |
| 20 | `felice_goswami2017.py` | vertebrates | raster | felice_goswami_2017 | Fig. 1C | Aves_crown_phenotypic_rate_continuous_data_GMM |
| 21 | `close2015.py` | vertebrates | raster | close_2015 | Fig. 1B-D | Mammalia_early_phenotypic_rate, Mammalia_early_phenotypic_disparity, Mammalia_early_taxonomic_diversity |
| 22 | `goswami2022.py` | vertebrates | raster | goswami_2022 | Fig. 3A | Placentalia_phenotypic_rate |
| 23 | `raw_data_convert.py` | vertebrates | raw | simoes_pierce_2021, simoes_2022, simoes_2026, zhang_2022, simoes_2020, close_2019 (converted) | raw data | Tetrapoda_stem_phenotypic_rate, Archosauromorpha_stem_phenotypic_rate, Lepidosauria_phenotypic_rate_continuous_data_GMM, Aves_stem_phenotypic_rate, Lepidosauria_phenotypic_rate_discrete_data, Lepidosauria_phenotypic_disparity, Tetrapoda_stem_taxonomic_diversity, Lepidosauria_taxonomic_diversity |
| 24 | `zhou2023.py` | invertebrates | vector | zhou_et_al_2023 | Fig. 1a | Bivalvia_phenotypic_disparity, Bivalvia_taxonomic_diversity |
| 25 | `paterson2019.py` | invertebrates | vector | paterson_et_al_2019_withSI | Fig. 3 | Trilobita_phenotypic_rate, Trilobita_phenotypic_rate_rawpoints |
| 26 | `novack_gottshall2022.py` | invertebrates | vector | novack-gottshall_et_al_2022 | Fig. 1b, 1c | Echinodermata_phenotypic_rate, Echinodermata_phenotypic_disparity |
| 27 | `hopkins_smith2015.py` | invertebrates | vector | hopkins_smith_2015 | Fig. 1, 2A | Echinoidea_Irregularia_phenotypic_rate, Echinoidea_Irregularia_phenotypic_rate_Fig1_allEchinoidea |
| 28 | `boivin2018.py` | invertebrates | vector | boivin-et_al_2018 | Fig. 3 | Echinoidea_Irregularia_phenotypic_disparity, Echinoidea_Irregularia_taxonomic_diversity |
| 29 | `fan2020.py` | invertebrates | raster | fan_et_al_2020_SI | Fig. S10 | Trilobita_taxonomic_diversity, Pancrustacea_taxonomic_diversity |
| 30 | `wills1998.py` | invertebrates | raster | wills_1998 | Fig. 7a | Pancrustacea_phenotypic_disparity |
| 31 | `foote1993.py` | invertebrates | raster | foote_1993 | Fig. 5A | Trilobita_phenotypic_disparity |
| 32 | `wright2017.py` | invertebrates | raster | Wright_2017_SI | Figs S6, S8, S10 | Crinoidea_Eucladida_phenotypic_rate, Crinoidea_Eucladida_phenotypic_disparity, Crinoidea_Eucladida_taxonomic_diversity |
| 33 | `guo2024.py` | mass_extinction | vector | guo_2024 | Fig. 2a-d | Brachiopoda_Terebratuliida_phenotypic_disparity, Brachiopoda_Terebratuliida_taxonomic_diversity, Brachiopoda_Rhynchonellida_phenotypic_disparity, Brachiopoda_Rhynchonellida_taxonomic_diversity, Brachiopoda_Spiniferida_phenotypic_disparity, Brachiopoda_Spiniferida_taxonomic_diversity, Brachiopoda_Athyridida_phenotypic_disparity, Brachiopoda_Athyridida_taxonomic_diversity |
| 34 | `whalen2020.py` | mass_extinction | raster | whalen_2020 | Fig. 2A, 2E | Ammonoidea_taxonomic_diversity, Ammonoidea_phenotypic_disparity |
| 35 | `cole_hopkins2021.py` | mass_extinction | raster | cole_hopkins_2021 | Fig. 1, 2B | Crinoidea_Diplobathrida_taxonomic_diversity, Crinoidea_Diplobathrida_phenotypic_disparity |
| 36 | `lamsdell_selden2016.py` | mass_extinction | raster | lamsdell_selden_2016 | Fig. 2, 4 | Eurypterida_taxonomic_diversity, Eurypterida_phenotypic_disparity |
| 37 | `halliday_goswami2016.py` | mass_extinction | raster | halliday_goswami_2016 | Fig. 4 | Placentalia_phenotypic_disparity |
| 38 | `raw_data_squamata_crown.py` | mass_extinction | raw | simoes_2026, simoes_2020, close_2019 (converted) | raw data | Squamata_crown_phenotypic_rate_continuous_data_GMM, Squamata_crown_phenotypic_rate_discrete_data, Squamata_crown_phenotypic_disparity, Squamata_crown_taxonomic_diversity |
| 39 | `qc_panels_vert.py` | qc | qc |  |  | source crops + QC panels of the vertebrate figure series |
| 40 | `qc_panels_invert.py` | qc | qc |  |  | source crops + QC panels of the invertebrate figure series |
| 41 | `qc_panels_massext.py` | qc | qc |  |  | source crops + QC panels of the mass-extinction figure series |
| 42 | `age_cutoffs.py` | cutoff | utility |  |  | PI age cutoffs (age_cutoffs.csv); ALWAYS LAST |

Groups: `close2019` (Close et al. 2019 locality data, steps 2-7, in this order), `vertebrates` and `invertebrates`
(clade-origin key), `mass_extinction` (new sources of the mass-extinction key, 14 Oct 2026), `qc` (source crops and
QC panels), `cutoff` (always last). `hopkins_smith2015.py` also writes
`Echinoidea_Irregularia_phenotypic_rate_Fig1_allEchinoidea.csv` (all echinoids, Fig. 1), which is not used by the keys.


## Adding a new figure

1. Copy the closest script in `source_papers/` and name it after the paper, keeping its two import lines (vector: `guo2024.py`, `navalon2022.py`; raster line:
   `whalen2020.py`, `rawson2022.py`; raster points with error bars: `lamsdell_selden2016.py`, `wright2017.py`;
   raw node data: `raw_data_convert.py`).
2. Calibrate both axes with `lin()` on as many ticks as possible, extract the curve, and build a table with
   `Age_Ma`, `Value` (plus `Err_lo`, `Err_hi`, `Value_err` when the figure shows an interval).
3. Write it with `save_series(df, clade, parameter, err_type, source_text, age_sd)`.
4. Add a QC panel to the matching `utilities/qc_panels_*.py` (`pdf_crop` / `raster_crop` + `qc_panel`), an entry in
   `REGISTRY` (section 10 of the master) and a note in `notes_*.csv` for any judgement call (line choice, axis,
   excluded points).
5. Run `digitise_master.py --run <script> qc_panels_massext age_cutoffs` (with the matching QC script), then `--check`.

The three batches of scripts write with slightly different rounding (`save_vert` 6 decimals, `save_invert` 4,
`save_massext` 5); these writers are kept so that existing series are reproduced exactly. New scripts should use
`save_series()`.


