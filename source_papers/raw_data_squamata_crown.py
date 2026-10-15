"""raw_data_squamata_crown.py -- step 38 of the digitisation protocol: simoes_2026, simoes_2020, close_2019 (converted), raw data [raw]
Method  : Squamata rates binned with node_bins(); disparity and diversity reformatted from the same files as Lepidosauria
Outputs : data_digitized/Squamata_crown_phenotypic_rate_continuous_data_GMM.csv, data_digitized/Squamata_crown_phenotypic_rate_discrete_data.csv, data_digitized/Squamata_crown_phenotypic_disparity.csv, data_digitized/Squamata_crown_taxonomic_diversity.csv
Run     : python3 R/digitisation/source_papers/raw_data_squamata_crown.py from the Runs root, or python3 R/digitisation/digitise_master.py --run raw_data_squamata_crown
Helpers : shared functions imported from R/digitisation/digitise_master.py (documented there)
(formerly R/digitise_massext/squamata_raw.py)

Squamata_crown series of the mass-extinction key (raw / converted data; no digitisation).
  phenotypic_rate continuous data (GMM): simoes_2026_squam_rates_cont.csv (age_median, rate_BT)
  phenotypic_rate discrete data:          simoes_2026_squam_rates_disc.csv (node_age_median, relrates; the key's
                                          "ages = age; median; ci95" describes the disparity file, not this one)
  phenotypic_disparity:                   simoes_2020_lepido_disparity.csv (same source as Lepidosauria, as in the key)
  taxonomic_diversity:                    Squamata_max_local_richness.csv (Close 2019 Fig. 5c, squamates only)
Node-based rates: 5-Myr bins of node ages (bin median, SE from pooled robust SD), exactly as the Lepidosauria rates
in R/digitisation/source_papers/raw_data_convert.py. For the mass-extinction runs, R/prep_mass_extinction.py re-bins these raw node
files with bins anchored at the extinction age.
"""
import os, sys; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))   # R/digitisation/, where digitise_master.py lives
from digitise_master import *                                   # shared helpers (see digitise_master.py)
out = save_table             # writer for complete tables (digitise_master.save_table)
out(node_bins("simoes_2026_squam_rates_cont.csv", "age_median", "rate_BT", "Squamata_crown", "phenotypic_rate", "simoes_2026", drop_zero=True),
    "Squamata_crown_phenotypic_rate_continuous_data_GMM.csv")
out(node_bins("simoes_2026_squam_rates_disc.csv", "node_age_median", "relrates", "Squamata_crown", "phenotypic_rate", "simoes_2026"),
    "Squamata_crown_phenotypic_rate_discrete_data.csv")
u = pd.read_csv(PDF + "simoes_2020_lepido_disparity.csv").drop_duplicates("age")
dd = pd.DataFrame({"Clade": "Squamata_crown", "Parameter": "phenotypic_disparity", "Age_Ma": u.age, "Value": u["median"], "Value_err": u.ci95 / 1.96,
                   "Err_lo": u["median"] - u.ci95, "Err_hi": u["median"] + u.ci95, "Err_type": "95CI", "Age_cal_sd": 0.0,
                   "Source": "simoes_2020 (raw data simoes_2020_lepido_disparity.csv: age, median, ci95; 8 time bins, 100 bootstrap replicates each; same file as Lepidosauria, as listed in the key)"})
out(dd, "Squamata_crown_phenotypic_disparity.csv")
q = pd.read_csv(OUT + "Squamata_max_local_richness.csv"); q["Clade"] = "Squamata_crown"; q["Parameter"] = "taxonomic_diversity"
q["Source"] = "close_2019 Fig. 5c converted data (Squamata_max_local_richness.csv): maximum local species richness per ICS stage (squamates only), +/-1 species reading error"
out(q, "Squamata_crown_taxonomic_diversity.csv")
