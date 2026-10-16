"""age_cutoffs.py -- step 42 of the digitisation pipeline [utility]
Method  : drops rows older than the age cutoffs (input_keys/age_cutoffs.csv); rewrites a file only when its content changes
Note    : run after all series are written (the digitisation scripts write full series)
Run     : python code/utilities/age_cutoffs.py (from any folder), or python code/digitise_master.py --run age_cutoffs

For every row of input_keys/age_cutoffs.csv (Clade, Parameter, Max_age_Ma, Applied_to, Reason), rows with
Age_Ma > Max_age_Ma are dropped from outputs/data_tables/clade_origin/<Clade>_<Parameter>.csv and the cutoff is
recorded in the Source column. Idempotent: a file that is already cut is not rewritten; series not present are
skipped. Node-based rate series are cut on node ages inside raw_data_convert.py (before binning); for them this
script only re-checks the bin ages.
"""
import os, sys, re; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))   # code/, where digitise_master.py lives
from digitise_master import AGE_CUTOFFS, OUT                                         # paths (digitise_master.py, section 1)
import pandas as pd
cut = pd.read_csv(AGE_CUTOFFS)
for _, r in cut.iterrows():
    f = f"{OUT}{r.Clade}_{r.Parameter}.csv"; name = os.path.basename(f)
    if not os.path.exists(f): print(f"{name}: not present, skipped"); continue
    d = pd.read_csv(f); n0 = len(d)
    d = d[d.Age_Ma <= r.Max_age_Ma + 1e-9].copy()
    tag = f"; age cutoff: data <= {r.Max_age_Ma:g} Ma ({r.Reason})"
    d["Source"] = [re.sub(r"; age cutoff: .*$", "", s) + tag for s in d.Source.astype(str)]
    old = open(f, encoding="utf-8", newline="").read()
    new = d.to_csv(index=False, lineterminator="\r\n" if "\r\n" in old else "\n")
    if new == old:                       # already cut: leave the file (and its date) untouched
        print(f"{name}: unchanged ({n0} rows)"); continue
    with open(f, "w", encoding="utf-8", newline="") as g: g.write(new)
    print(f"{name}: {n0} -> {len(d)} rows; oldest {d.Age_Ma.max():.1f} Ma")
