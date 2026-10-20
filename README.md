# Figure digitisation pipeline

Python code that extracts time series of taxonomic diversity, phenotypic disparity and phenotypic evolutionary rates
directly from figures in published papers (PDF files), and converts raw data files into the same table format. Every
series gets a QC figure that shows the source figure next to the digitised values.

## Repository layout

```
articles_PDF/                     source PDFs and raw data files (CSV) named in the keys
  figure_images/                  figure images read directly by three scripts (see "Inputs")
input_keys/                       keys of the analysis, example keys and settings (age_cutoffs.csv,
                                  massext_copy_from_origin.csv)
code/
  digitise_master.py              shared functions + registry of scripts + command-line runner
  source_papers/                  one script per source paper (37 scripts)
  utilities/                      steps not tied to one paper (6 scripts)
outputs/
  data_tables/clade_origin/       one table per digitised series (full time range)
  data_tables/mass_extinction/    the same series cropped to the time after each mass extinction (+ crop_log.csv)
  figures/<parameter>/            source crop + QC panel (source | digitised series) per series
  intermediates/                  extracted rasters, tick tables, trace plots (created when the code runs; not tracked)
```

* **`code/digitise_master.py`**: the master file. It is
  1. the **shared library**: every helper used by more than one script (axis calibration, PDF vector-path reading,
     raster masks and tracing, regridding, node-rate binning, output writers, QC panels), each with a commented
     explanation;
  2. the **registry**: one entry per script with source paper, figure, method (vector / raster / raw / derived),
     output series, run order and the steps it depends on;
  3. the **runner**: list, run and check the scripts from the command line (below).
* **`code/source_papers/`: one script per source paper**, named after the paper (`guo2024.py`, `foth2021.py`, ...).
  Each script starts with a header giving its step number, source figure, method and outputs. It then adds `code/` to
  the Python path and imports the shared functions with `from digitise_master import *` (these two lines must stay at
  the top of any new script). Figure-specific values (tick positions, line colours, crop boxes) stay in the script, so
  each figure's extraction can be read and audited on its own.
* **`code/utilities/`: the steps that serve all sources**: `extract_raster_figures.py` (runs first: saves the embedded
  raster figures that the tracing scripts read), `qc_panels_vert.py`, `qc_panels_invert.py`, `qc_panels_massext.py`
  (source crops + QC panels), `age_cutoffs.py` (age cutoffs, after all series are written) and
  `mass_extinction_crop.py` (mass-extinction keys only; last).

All paths are built from the location of `code/`, so the commands below work from any folder (shown here from the
repository root).

## Installation

Python 3.9 or later with numpy, pandas, scipy, pillow, matplotlib, pymupdf and scikit-image (`wright2017.py` only):

```
pip install -r requirements.txt
```

The tested versions are listed in `requirements.txt`. Use `python3` instead of `python` where that is the name of the
Python 3 interpreter (macOS, Linux).

## Quick start: run the example keys

```
python code/digitise_master.py --key input_keys/key_test_origin.csv
python code/digitise_master.py --key input_keys/key_test_MassExct.csv
```

A key has one row per clade x parameter (`clade`, `parameter`, `reference`, `figure/table`; optional `data_type`).
For each row the runner finds the script that writes that series (here `wright2017.py`, for Wright 2017 SI Figs S6,
S8 and S10), runs it together with the steps it depends on and its QC step, and then applies the age cutoffs. A
mass-extinction key has two more columns, `extinction_Age` (Ma) and `extinction_name` (`extinction_age(Ma)` and
`extinction_event` are also accepted); each series is then also cropped
to the time after that extinction. The example keys write:

* `outputs/data_tables/clade_origin/Crinoidea_Eucladida_{phenotypic_rate,phenotypic_disparity,taxonomic_diversity}.csv`
* `outputs/figures/<parameter>/Crinoidea_Eucladida_*_Wright2017_FigS*.png` and `..._QC_digitised.png`
* `outputs/data_tables/mass_extinction/Crinoidea_Eucladida_*__end_Ordovician.csv`, `..._end_Devoninan_Famennian.csv`
  and their rows in `crop_log.csv`

The same command runs any key in this format, as long as its rows refer to series in the registry
(`python code/digitise_master.py --list --long` shows every series and its source).

## Reproduce all tables: the keys of the analysis

```
python code/digitise_master.py --key "input_keys/key_origin_invert&vert.csv"
python code/digitise_master.py --key "input_keys/key_massext_invert&vert.csv"
```

(The quotes are needed because the file names contain `&`.)

* `key_origin_invert&vert.csv` (clade origins): 52 rows (clade x parameter); 48 have a published source and are
  written to `outputs/data_tables/clade_origin/`, 4 are `na` (no published series).
* `key_massext_invert&vert.csv` (mass extinctions): 63 rows (17 clades, six extinctions). 37 series are cropped to
  the time after their extinction (`outputs/data_tables/mass_extinction/`; one of them, Crinoidea_Eucladida disparity
  after the Famennian, keeps only 3 points); 15 rows of the five clades in `massext_copy_from_origin.csv` use the
  full clade-origin series; 11 rows are `na`.

Re-running the two keys reproduces the committed tables exactly. `outputs/data_tables/clade_origin/` also holds a few
tables that no key row names directly: the Close et al. 2019 locality summaries (`Tetrapoda_nonflying_*`,
`Squamata_*max_local_richness.csv`, read by `raw_data_convert.py`), `Trilobita_phenotypic_rate_rawpoints.csv`
(per-branch points behind the binned Trilobita rates) and the four `Squamata_crown_*` series written by
`raw_data_squamata_crown.py` (Squamata_crown is not part of the final mass-extinction key).

## Other commands

```
python code/digitise_master.py --list            # all scripts in run order (add --long for method, outputs, dependencies)
python code/digitise_master.py --which Eurypterida_phenotypic_disparity   # which script makes a series
python code/digitise_master.py --run all         # re-run every digitisation script in the correct order
python code/digitise_master.py --run whalen2020 age_cutoffs              # selected scripts
python code/digitise_master.py --run group:mass_extinction               # one group
python code/digitise_master.py --check           # re-run everything in a scratch copy and compare with outputs/data_tables/
python code/source_papers/whalen2020.py          # any script can also be run on its own
python code/utilities/mass_extinction_crop.py input_keys/key_test_MassExct.csv
```

Run order: `extract_raster_figures.py` first; the Close et al. 2019 scripts in their numbered order, before
`raw_data_convert.py` and `raw_data_squamata_crown.py` (which read their tables); the QC panels; `age_cutoffs.py`;
`mass_extinction_crop.py` last. `--run all` and `--key` keep this order. The QC panels are drawn before the cutoffs,
so they show the full traced series (a `qc_panels_*.py` script re-run on its own after `age_cutoffs.py` draws the cut
series instead).

`--check` copies `articles_PDF/`, `input_keys/` and `code/` to a temporary folder, runs every script there and compares
each table with `outputs/data_tables/` (byte for byte; if the bytes differ, the values are compared). Add
`--key <mass-extinction key>` to check the crops of that key too (`--check --key "input_keys/key_massext_invert&vert.csv"`
checks every committed table). Figures are compared pixel by pixel and reported
separately, because PNG encoding and anti-aliasing can differ slightly between library versions. On Windows the tables
are written with CRLF line endings; they are then reported as "same values (format differs)".

## Inputs

* `articles_PDF/`: the source PDFs and the raw data CSVs named in the keys (`simoes_*.csv`, `zhang_2022_rates.csv`).
* `articles_PDF/figure_images/`: figure images read directly rather than extracted by the scripts:
  * `wright_2017_SI_FigS6.png`, `wright_2017_SI_FigS8.png`, `wright_2017_SI_FigS10.png`: crops of the Wright 2017 SI
    figures read by `wright2017.py`;
  * `fan_et_al_2020_SI_FigS10.png`: the embedded image of `fan_et_al_2020_SI.pdf` p. 17 (Fig. S10), extracted with
    poppler (`pdfimages -png -f 17 -l 17`). It is supplied because this image is stored in CMYK, and PyMuPDF converts
    it to RGB slightly differently;
  * `close_2019_Fig5c_dot_tops_zoom.png`: zoomed view of the detected dot tops, shown in the Close et al. 2019 Fig. 5c
    QC panel.
* `input_keys/`:
  * `key_origin_invert&vert.csv`, `key_massext_invert&vert.csv`: the keys of the analysis (above). Columns: `group`,
    `clade`, (mass extinctions: `extinction_event`, `extinction_age(Ma)`), `parameter`, `reference`, `figure/table`,
    `legend/line`, `note` (extraction instructions), `data_file` (raw data rows), `data_type` (tells apart two
    series of one clade x parameter, e.g. `PyRate` / `FBDD`);
  * `key_test_origin.csv`, `key_test_MassExct.csv`: example keys (format above);
  * `age_cutoffs.csv`: series cut at a maximum age (`Clade, Parameter, Max_age_Ma, Applied_to, Reason`), applied by
    `age_cutoffs.py` (node-based rates are cut on node ages inside `raw_data_convert.py`);
  * `massext_copy_from_origin.csv`: clades whose full clade-origin series are used for their extinction (clades that
    originate right after it, or whose full series show the parameters across it); they are not cropped by
    `mass_extinction_crop.py`.

## Outputs

* `outputs/data_tables/clade_origin/<Clade>_<Parameter>[_<data_type>].csv`, one row per time point:
  `Clade, Parameter, Age_Ma, Value, Value_err, Err_lo, Err_hi, Err_type, Source, Age_cal_sd`.
  `Err_type` is `95CI`, `SE`, `Poisson` (counts without published error) or `none`; `Age_cal_sd` is the SD of the
  age calibration (Myr), so that the digitisation error can be propagated into downstream analyses (e.g. by jittering
  ages in Monte-Carlo replicates). `Source` records the reference, figure, line and method (and any age cutoff).
  `wright2017.py` also writes `Rate_class` and the marker pixel positions. A few Close et al. 2019 tables are stage
  summaries (`*_stage_max_local_richness.csv`).
* `outputs/data_tables/mass_extinction/<series>__<event>.csv`: rows with `Age_Ma` <= the key's extinction age (node-based
  raw rates are re-binned from the node ages, in 5-Myr bins anchored at the extinction age); `crop_log.csv` has one
  row per cropped series (points before and after the crop). Running a smaller key updates only its rows of
  `crop_log.csv`; rows of a key that are no longer cropped (clade added to `massext_copy_from_origin.csv`) are
  removed from the log, and their old crop files are reported.
* `outputs/figures/<parameter folder>/`: a source crop and a `_QC_digitised.png` panel (source | digitised series) per
  series.
* `outputs/intermediates/`: extracted rasters (`*_img/`), trace-over-source images (`vert_qc/`), tick and calibration
  tables, the Close et al. 2019 locality tables and summary plots. Regenerated on every run.

## Method in brief

* **Vector figures** (curves stored as PDF drawing paths): coordinates read directly with PyMuPDF; exact apart from
  the axis calibration.
* **Raster figures** (embedded images): the image is extracted from the PDF, the curve isolated with colour or
  darkness masks and traced per pixel column (median of the line pixels); markers and error bars are isolated by
  morphological opening.
* **Calibration**: least-squares linear fit of pixel (or PDF point) position on all available tick marks or labels;
  the fit RMSE combined with half a pixel is stored as `Age_cal_sd`.
* Dense traces are averaged on a 1-Myr grid; node-based raw rate files are summarised in 5-Myr bins (`node_bins`).

Each script header and `--list --long` describe the specific extraction; `digitise_master.py` documents every helper.

## Scripts in run order (groups `utility`, `qc`, `cutoff`, `crop` in `code/utilities/`; all others in `code/source_papers/`)

| step | script | group | method | source (key reference) | figure | series written to `outputs/data_tables/clade_origin/` |
|---:|---|---|---|---|---|---|
| 1 | `extract_raster_figures.py` | utility | utility | foth_2021, sookias_2026, quintero_2025_SI, rawson_2022, crouch_tobias_2022, felice_goswami_2017, close_2015, goswami_2022, allvarez-carretero_2022, close_2019, wills_1998, foote_1993 | embedded rasters | Run first: extracts the raster figures read by the tracing scripts to outputs/intermediates/ |
| 2 | `close2019_1_extract_localities.py` | close2019 | vector | close_2019 | Fig. 1a | locality richness points -> intermediates/close2019/close2019_fig1a_localities.csv; source crop + QC panel |
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
| 26 | `novack_gottshall2022.py` | invertebrates | vector (ED Fig. 1a: raster) | novack-gottshall_et_al_2022 | Fig. 1b, 1c; ED Fig. 1a | Echinodermata_phenotypic_rate, Echinodermata_phenotypic_disparity, Echinodermata_taxonomic_diversity |
| 27 | `hopkins_smith2015.py` | invertebrates | vector | hopkins_smith_2015 | Fig. 2A (Fig. 1) | Echinoidea_Irregularia_phenotypic_rate |
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
| 42 | `age_cutoffs.py` | cutoff | utility |  |  | age cutoffs (input_keys/age_cutoffs.csv); after all series are written |
| 43 | `mass_extinction_crop.py` | crop | utility |  |  | mass-extinction keys only: crops each series at its extinction age -> outputs/data_tables/mass_extinction/ |

Groups: `close2019` (Close et al. 2019 locality data, steps 2-7, in this order), `vertebrates` and `invertebrates`
(clade-origin sources), `mass_extinction` (sources added for the mass-extinction series; `halliday_goswami2016.py`
also serves the clade-origin key, Placentalia disparity cut at 83.3 Ma), `qc` (source crops and QC panels), `cutoff`
and `crop`. `hopkins_smith2015.py` also writes
`Echinoidea_Irregularia_phenotypic_rate_Fig1_allEchinoidea.csv` (all crown-group echinoids, Fig. 1) to
`outputs/intermediates/` for reference.

## Adding a new figure

1. Copy the closest script in `code/source_papers/` and name it after the paper, keeping its two import lines
   (vector: `guo2024.py`, `navalon2022.py`; raster line: `whalen2020.py`, `rawson2022.py`; raster points with error
   bars: `lamsdell_selden2016.py`, `wright2017.py`; raw node data: `raw_data_convert.py`).
2. Calibrate both axes with `lin()` on as many ticks as possible, extract the curve, and build a table with
   `Age_Ma`, `Value` (plus `Err_lo`, `Err_hi`, `Value_err` when the figure shows an interval).
3. Write it with `save_series(df, clade, parameter, err_type, source_text, age_sd)`.
4. Add a QC panel to the matching `code/utilities/qc_panels_*.py` (`pdf_crop` / `raster_crop` + `qc_panel`) and an
   entry in `REGISTRY` (section 10 of the master; `requires` lists steps whose outputs the script reads, `qc` the
   script that draws its QC panel).
5. Add a row for the new series to a key and run `python code/digitise_master.py --key <your key>`, then `--check`.

The three groups of scripts write with slightly different rounding (`save_vert` 6 decimals, `save_invert` 4,
`save_massext` 5); these writers are kept so that the published tables are reproduced exactly. New scripts should use
`save_series()`.
