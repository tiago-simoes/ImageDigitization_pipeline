"""raw_data_convert.py -- step 23 of the digitisation protocol: simoes_pierce_2021, simoes_2022, simoes_2026, zhang_2022, simoes_2020, close_2019 (converted), raw data [raw]
Method  : node-based rate files binned with node_bins(); Lepidosauria disparity per time slice; Close 2019 tables reformatted
Outputs : data_digitized/Tetrapoda_stem_phenotypic_rate.csv, data_digitized/Archosauromorpha_stem_phenotypic_rate.csv, data_digitized/Lepidosauria_phenotypic_rate_continuous_data_GMM.csv, data_digitized/Aves_stem_phenotypic_rate.csv, data_digitized/Lepidosauria_phenotypic_rate_discrete_data.csv, data_digitized/Lepidosauria_phenotypic_disparity.csv, data_digitized/Tetrapoda_stem_taxonomic_diversity.csv, data_digitized/Lepidosauria_taxonomic_diversity.csv
Note    : needs the Close 2019 outputs (run after close2019_*)
Run     : python3 R/digitisation/source_papers/raw_data_convert.py from the Runs root, or python3 R/digitisation/digitise_master.py --run raw_data_convert
Helpers : shared functions imported from R/digitisation/digitise_master.py (documented there)
(formerly R/digitise_vert/raw_convert.py)

Convert raw / converted data files listed in the key into the pipeline series format
(Clade, Parameter, Age_Ma, Value, Value_err, Err_lo, Err_hi, Err_type, Source, Age_cal_sd).
Node-based rate files (one row per branch/node): binned in 5-Myr bins by node age;
value = bin median (robust to the long right tail of relative rates), error = SE approximated
as 1.2533 * pooled robust SD (MAD-based) / sqrt(n); age = mean node age in the bin.
"""
import os, sys; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))   # R/digitisation/, where digitise_master.py lives
from digitise_master import *                                   # shared helpers (see digitise_master.py)
D = OUT                      # data_digitized/ (converted Close 2019 tables are read from there)
save = save_table            # writer for complete tables (digitise_master.save_table); node_bins() is also from digitise_master
save(node_bins("simoes_pierce_2021_tetrapods.csv","node_age_median","relrates","Tetrapoda_stem","phenotypic_rate","simoes_pierce_2021",max_age=380.0),"Tetrapoda_stem_phenotypic_rate.csv")   # PI 9 Oct 2026: nodes > 380 Ma are not Tetrapoda (8 Oct: 390 Ma); see data_digitized/digitisation/age_cutoffs.csv
save(node_bins("simoes_2022_archo_stem.csv","node_age_median","rates","Archosauromorpha_stem","phenotypic_rate","simoes_2022"),"Archosauromorpha_stem_phenotypic_rate.csv")
save(node_bins("simoes_2026_lepido_rates_cont.csv","age_median","rate_BT","Lepidosauria","phenotypic_rate","simoes_2026",drop_zero=True),"Lepidosauria_phenotypic_rate_continuous_data_GMM.csv")
save(node_bins("zhang_2022_rates.csv","node_age_median","rates","Aves_stem","phenotypic_rate","zhang_2022",drop_zero=True),"Aves_stem_phenotypic_rate.csv")   # 2 extant tips (age 0) excluded
save(node_bins("simoes_2026_lepido_rates_disc.csv","node_age_median","relrates","Lepidosauria","phenotypic_rate","simoes_2026"),"Lepidosauria_phenotypic_rate_discrete_data.csv")
# Lepidosauria disparity: one row per bootstrap replicate -> one row per age (median, ci95)
u=pd.read_csv(A+"simoes_2020_lepido_disparity.csv").drop_duplicates("age")
dd=pd.DataFrame({"Clade":"Lepidosauria","Parameter":"phenotypic_disparity","Age_Ma":u.age,"Value":u["median"],"Value_err":u.ci95/1.96,
                 "Err_lo":u["median"]-u.ci95,"Err_hi":u["median"]+u.ci95,"Err_type":"95CI","Age_cal_sd":0.0,
                 "Source":"simoes_2020 (raw data simoes_2020_lepido_disparity.csv: age, median, ci95; 8 time bins, 100 bootstrap replicates each)"})
save(dd,"Lepidosauria_phenotypic_disparity.csv")
# converted Close 2019 data
s=pd.read_csv(D+"Tetrapoda_nonflying_stage_max_local_richness.csv"); s=s[s.N_localities>0]
tt=pd.DataFrame({"Clade":"Tetrapoda_stem","Parameter":"taxonomic_diversity","Age_Ma":s.Mid_Ma,"Value":s.Max_local_richness,
   "Value_err":(s.Max_boot_hi-s.Max_boot_lo)/3.92,"Err_lo":s.Max_boot_lo,"Err_hi":s.Max_boot_hi,"Err_type":"95CI","Age_cal_sd":0.0,
   "Source":"close_2019 Fig. 1a converted data (Tetrapoda_nonflying_stage_max_local_richness.csv): maximum local species richness per ICS stage, bootstrap 95% CI"})
save(tt,"Tetrapoda_stem_taxonomic_diversity.csv")
q=pd.read_csv(D+"Squamata_max_local_richness.csv"); q["Clade"]="Lepidosauria"; q["Parameter"]="taxonomic_diversity"
q["Source"]="close_2019 Fig. 5c converted data (Squamata_max_local_richness.csv): maximum local species richness per ICS stage (squamates only), +/-1 species reading error"
save(q,"Lepidosauria_taxonomic_diversity.csv")
