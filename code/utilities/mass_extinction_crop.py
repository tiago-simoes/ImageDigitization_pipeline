"""mass_extinction_crop.py -- step 43 of the digitisation pipeline: mass-extinction keys only [utility]
Method  : keeps the rows with Age_Ma <= extinction age of each key row (node-based rates re-binned from the raw nodes)
Outputs : outputs/data_tables/mass_extinction/<series>__<event>.csv and crop_log.csv
Run     : python code/utilities/mass_extinction_crop.py input_keys/key_test_MassExct.csv
          (or python code/digitise_master.py --key input_keys/key_test_MassExct.csv, which first re-runs the scripts
          behind the key's rows)

For every row of a mass-extinction key (columns clade, parameter, reference, figure/table, extinction age
[extinction_Age or extinction_age(Ma)] and extinction name [extinction_name or extinction_event]; data_type optional):
  * the full series outputs/data_tables/clade_origin/<series>.csv is cut to the time AFTER the extinction, i.e. the rows
    with Age_Ma <= extinction age, and written to outputs/data_tables/mass_extinction/<series>__<event>.csv; the Source
    column records the crop. The full series are not modified.
  * node-based raw rate files (Tetrapoda_stem, Squamata_crown) are re-binned from the raw node ages instead (nodes <=
    extinction age, 5-Myr bins anchored at the extinction age, so the first bin is the first 5 Myr after the event).
  * clades listed in input_keys/massext_copy_from_origin.csv are not cropped (their clade-origin series are used).
Series ids: <clade>_<parameter>[_<data_type>] (digitise_master.key_series); the event slug is appended after "__".
crop_log.csv gets one row per cropped series; rows of series not in this key are kept, so running a small key does
not remove the log entries of other series. A series file is rewritten only when its content changes.
"""
import os, sys; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))   # code/, where digitise_master.py lives
import numpy as np, pandas as pd
from digitise_master import (node_bins, COLS, OUT, OUT_MASSEXT, KEYS, ROOT, read_key, key_series, is_massext_key, slug)

KEY = sys.argv[1] if len(sys.argv) > 1 else os.path.join(KEYS, "key_test_MassExct.csv")
MIN_N = 4                                   # fewer points after the extinction: flagged in the printed summary
key = read_key(KEY)
if not is_massext_key(key): sys.exit(f"{KEY}: not a mass-extinction key (needs extinction age and extinction name columns)")
key = key_series(key)
key["event_slug"] = key.extinction_event.map(slug); key["sid"] = key.base_sid + "__" + key.event_slug

cf = os.path.join(KEYS, "massext_copy_from_origin.csv")
COPY = dict(pd.read_csv(cf)[["Clade", "Reason"]].values) if os.path.exists(cf) else {}
# node-based raw rate files: re-binned per event (bins anchored at the extinction age)
NODE = {"Tetrapoda_stem_phenotypic_rate": ("simoes_pierce_2021_tetrapods.csv", "node_age_median", "relrates", "simoes_pierce_2021", dict(max_age=380.0)),
        "Squamata_crown_phenotypic_rate_continuous_data_GMM": ("simoes_2026_squam_rates_cont.csv", "age_median", "rate_BT", "simoes_2026", dict(drop_zero=True)),
        "Squamata_crown_phenotypic_rate_discrete_data": ("simoes_2026_squam_rates_disc.csv", "node_age_median", "relrates", "simoes_2026", {})}

summary, log = [], []
for _, r in key.iterrows():
    ea = float(r.extinction_age); ref = r.reference.strip(); base = r.base_sid
    if ref.lower() in ("na", ""):
        summary.append((r.sid, "no published source", "")); continue
    if r.clade in COPY:
        summary.append((r.sid, "not cropped (input_keys/massext_copy_from_origin.csv)", "")); continue
    f = OUT + base + ".csv"
    if not os.path.exists(f):
        summary.append((r.sid, "MISSING: " + os.path.relpath(f, ROOT), "")); continue
    full = pd.read_csv(f)
    tag = f"; mass-extinction crop: data <= {ea:g} Ma ({r.extinction_event})"
    if base in NODE:
        nf, ac, rc, src, kw = NODE[base]
        kw = dict(kw); mx = kw.pop("max_age", None)
        mx = ea if mx is None else min(mx, ea)
        g = node_bins(nf, ac, rc, r.clade, r.pr_std, src, max_age=mx, anchor=ea, **kw)
        g = g[COLS].copy(); g["Source"] = g.Source + tag
        cut = g.sort_values("Age_Ma", ascending=False); how = "raw nodes <= extinction age, 5-Myr bins anchored at the extinction age"
    else:
        cut = full[full.Age_Ma <= ea + 1e-9].copy()
        cut["Source"] = cut.Source.astype(str) + tag; how = "rows with Age_Ma <= extinction age"
    out = OUT_MASSEXT + r.sid + ".csv"
    new_txt = cut.to_csv(index=False)
    if not (os.path.exists(out) and open(out).read() == new_txt):     # rewrite only on change
        open(out, "w").write(new_txt)
    n = len(cut)
    first = round(cut.Age_Ma.max(), 2) if n else np.nan; last = round(cut.Age_Ma.min(), 2) if n else np.nan
    summary.append((r.sid, f"{n} points after the extinction" + (f" (fewer than {MIN_N})" if n < MIN_N else ""), os.path.relpath(out, ROOT)))
    log.append(dict(Series=r.sid, Extinction_event=r.extinction_event, Extinction_age_Ma=ea, Crop=how,
                    N_before=len(full), N_after=n, Oldest_after=first, Youngest=last))

L = pd.DataFrame(log, columns=["Series", "Extinction_event", "Extinction_age_Ma", "Crop", "N_before", "N_after", "Oldest_after", "Youngest"])
lf = OUT_MASSEXT + "crop_log.csv"
if os.path.exists(lf):                       # keep the rows of series not in this key; replace the others in place
    old = pd.read_csv(lf); new = L.set_index("Series")
    rows = [new.loc[s].to_dict() | {"Series": s} if s in new.index else old.iloc[i].to_dict() for i, s in enumerate(old.Series)]
    rows += [L.iloc[i].to_dict() for i, s in enumerate(L.Series) if s not in set(old.Series)]
    L = pd.DataFrame(rows)[L.columns]
L.to_csv(lf, index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 90)
print(pd.DataFrame(summary, columns=["series", "status", "file"]).to_string(index=False))
