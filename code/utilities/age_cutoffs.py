"""age_cutoffs.py -- step 42 of the digitisation protocol [utility]
Method  : drops rows older than the PI cutoffs (age_cutoffs.csv); rewrites a file only when its content changes
Note    : PI age cutoffs (age_cutoffs.csv); ALWAYS LAST
Run     : python3 R/digitisation/utilities/age_cutoffs.py from the Runs root, or python3 R/digitisation/digitise_master.py --run age_cutoffs
(formerly R/digitise_vert/age_cutoffs.py)

Apply PI age cutoffs (data_digitized/digitisation/age_cutoffs.csv) to pipeline series in data_digitized/.
Run from the Runs root AFTER the digitisation / raw_data_convert.py scripts (they write full series):
    python3 R/digitisation/utilities/age_cutoffs.py
Rows with Age_Ma > Max_age_Ma are dropped from data_digitized/<Clade>_<Parameter>.csv and the cutoff is
recorded in the Source column. Idempotent: a file that is already cut is not rewritten. Node-based rate series are cut on node ages inside
raw_data_convert.py (before binning); for them this script only re-checks the bin ages.
"""
import pandas as pd, re
cut = pd.read_csv("data_digitized/digitisation/age_cutoffs.csv")
for _, r in cut.iterrows():
    f = f"data_digitized/{r.Clade}_{r.Parameter}.csv"
    d = pd.read_csv(f); n0 = len(d)
    d = d[d.Age_Ma <= r.Max_age_Ma + 1e-9].copy()
    tag = f"; PI age cutoff: data <= {r.Max_age_Ma:g} Ma ({r.Reason})"
    d["Source"] = [re.sub(r"; PI age cutoff: .*$", "", s) + tag for s in d.Source.astype(str)]
    old = open(f, encoding="utf-8", newline="").read()
    new = d.to_csv(index=False, lineterminator="\r\n" if "\r\n" in old else "\n")
    if new == old:                       # already cut: leave the file (and its date) untouched
        print(f"{f}: unchanged ({n0} rows)"); continue
    with open(f, "w", encoding="utf-8", newline="") as g: g.write(new)
    print(f"{f}: {n0} -> {len(d)} rows; oldest {d.Age_Ma.max():.1f} Ma")
