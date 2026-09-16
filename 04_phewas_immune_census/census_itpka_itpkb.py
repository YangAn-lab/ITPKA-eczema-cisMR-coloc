import cellxgene_census, pandas as pd, numpy as np

GENES = {"ITPKA": "ENSG00000137825", "ITPKB": "ENSG00000143772"}
results = []
with cellxgene_census.open_soma(census_version="2025-01-30") as census:
    for tissue, tlabel in [("skin of body", "skin"), ("blood", "blood")]:
        for gname, gid in GENES.items():
            ad = cellxgene_census.get_anndata(
                census, organism="homo_sapiens",
                obs_value_filter=f"tissue_general == '{tissue}' and disease == 'normal' and is_primary_data == True",
                var_value_filter=f"feature_id == '{gid}'",
                obs_column_names=["cell_type", "donor_id"],
                X_name="normalized")
            X = np.asarray(ad.X.todense()).ravel()
            df = pd.DataFrame({"cell_type": ad.obs["cell_type"].values, "expr": X})
            df["is_t"] = df["cell_type"].str.contains("T cell|thymocyte", case=False, regex=True)
            for subset, sub in [("all", df), ("T_cells", df[df["is_t"]])]:
                if len(sub) == 0: continue
                results.append({
                    "tissue": tlabel, "gene": gname, "subset": subset,
                    "n_cells": len(sub), "n_pos": int((sub["expr"]>0).sum()),
                    "pct_pos": round(100*(sub["expr"]>0).mean(), 3),
                    "mean_expr": round(float(sub["expr"].mean()), 4),
                    "mean_pos": round(float(sub.loc[sub["expr"]>0,"expr"].mean()), 4) if (sub["expr"]>0).any() else 0.0})
            print(tlabel, gname, "done", flush=True)
out = pd.DataFrame(results)
out.to_csv("census_itpka_vs_itpkb.csv", index=False)
print(out.to_string(index=False))
