ext = pd.read_csv('/workspace/coloc_panel/coloc_panel_extended_results.csv')
cols = ['gene','origin','endpoint','nsnps','PP.H3','PP.H4','PP.H4_p12_1e6','PP.H4_p12_1e4','PP.H4_w250','PP.H4_w125','min_p_eqtl']
v = ext[cols].copy()
for c in ['PP.H3','PP.H4','PP.H4_p12_1e6','PP.H4_p12_1e4','PP.H4_w250','PP.H4_w125']:
    v[c] = v[c].round(4)
print(v.to_string(index=False))