#!/usr/bin/env python3
"""Build P2-02 supplementary tables S9-S17 (new) and update S1, S2, S3, S8.
All sources verified against manuscript v4 legends. Output: /workspace/supplementary_v4/
"""
import pandas as pd
import numpy as np
import subprocess, io, os

SRC = "/workspace/fable_fix/pqtl"
OUT = "/workspace/supplementary_v4"
os.makedirs(OUT, exist_ok=True)

def w(df, name):
    df.to_csv(f"{OUT}/{name}", index=False)
    print(f"{name}: {len(df)} rows")

# ---------- S9: plasma pQTL colocalization ----------
abf = pd.read_csv(f"{SRC}/pqtl_coloc_abf_results.csv")
abf.insert(0, "analysis", "coloc.abf (single-signal)")
# susie CS-level from RDS
r = subprocess.run(["Rscript", "-e", '''
x <- readRDS("/workspace/fable_fix/pqtl/pqtl_coloc_susie_results.rds")
for (nm in names(x)) {
  s <- x[[nm]]$summary
  if (!is.null(s) && nrow(as.data.frame(s)) > 0) {
    d <- as.data.frame(s)
    d$pair <- nm
    write.csv(d, stdout(), row.names=FALSE)
  } else {
    cat(sprintf("EMPTY,%s\\n", nm))
  }
}'''], capture_output=True, text=True)
blocks, empties = [], []
cur = []
for line in r.stdout.splitlines():
    if line.startswith("EMPTY,"):
        empties.append(line.split(",",1)[1]); continue
    if line.startswith('"nsnps"') or line.startswith("nsnps"):
        if cur: blocks.append("\n".join(cur)); cur = []
    cur.append(line)
if cur: blocks.append("\n".join(cur))
sus = pd.concat([pd.read_csv(io.StringIO(b)) for b in blocks], ignore_index=True)
sus = sus.rename(columns={"hit1":"lead_pqtl_cs","hit2":"lead_gwas_cs","idx1":"pqtl_cs_index","idx2":"gwas_cs_index"})
sus.insert(0, "analysis", "coloc.susie (credible-set level)")
for c in ["pqtl","protein","cohort","N_pqtl"]:
    sus[c] = ""
for i, row in sus.iterrows():
    pq, ep = row["pair"].split(" x ")
    sus.loc[i, "pqtl"] = pq; sus.loc[i, "endpoint"] = ep
    sus.loc[i, "protein"] = "TYRO3"
    sus.loc[i, "cohort"] = "Icelandic (Gudjonsson 2022)" if pq.startswith("ICE") else "HELIC-MANOLIS (Gilly 2020)"
sus["note"] = ""
for e in empties:
    pq, ep = e.split(" x ")
    sus.loc[len(sus)] = {"analysis":"coloc.susie (credible-set level)","pqtl":pq,"protein":"TYRO3",
        "cohort":"Icelandic (Gudjonsson 2022)" if pq.startswith("ICE") else "HELIC-MANOLIS (Gilly 2020)",
        "endpoint":ep,"note":"no credible set constructed on at least one side; CS-level colocalization not computable"}
abf["note"] = ""
common_cols = ["analysis","pqtl","protein","cohort","endpoint","nsnps","PP.H0","PP.H1","PP.H2","PP.H3","PP.H4",
               "H4_p12_1e6","H4_p12_5e6","H4_p12_5e5","H4_p12_1e4","min_p_pqtl","min_p_gwas","lead_pqtl","lead_gwas","note"]
for c in common_cols:
    if c not in sus.columns: sus[c] = np.nan
s9 = pd.concat([abf[[c for c in common_cols if c in abf.columns]], sus[common_cols]], ignore_index=True)
w(s9, "P2-02_SupplementaryTableS9_plasma_pQTL_coloc.csv")

# ---------- S10: exhaustive +/-1 Mb locus panel ----------
genes = pd.read_csv(f"{SRC}/locus_genes_1mb.csv")
eq = pd.read_csv(f"{SRC}/rs11635906_all_locus_genes_eqtl.csv")
ext = pd.read_csv(f"{SRC}/coloc_panel_extension_results.csv")
panel12 = ['CHP1','EHD4','ITPKA','JMJD7','LINC_260814','NDUFAF1','NUSAP1','OIP5','OIP5-AS1','PLA2G4B','RPAP1','RTF1']
ensg_panel = {"ENSG00000260814":"LINC_260814"}  # unnamed lncRNA in panel under this placeholder
eq_sub = eq[["ensg","beta","se","p","eaf"]].rename(columns={"beta":"rs11635906_eqtl_beta","se":"rs11635906_eqtl_se","p":"rs11635906_eqtl_p","eaf":"rs11635906_eaf"})
s10 = genes.merge(eq_sub, on="ensg", how="left")
s10["in_12gene_panel"] = s10.apply(lambda r: r["symbol"] in panel12 or r["ensg"] in ensg_panel, axis=1)
# coloc for 6 out-of-panel genes x 3 endpoints -> wide
ext_w = ext.pivot_table(index="gene", columns="endpoint", values=["PP.H3","PP.H4"], aggfunc="first")
ext_w.columns = [f"coloc_{a}_{b.replace('_FinnGen','')}" for a,b in ext_w.columns]
ext_w = ext_w.reset_index().rename(columns={"gene":"symbol"})
s10 = s10.merge(ext_w, on="symbol", how="left")
s10 = s10.sort_values("start")
w(s10, "P2-02_SupplementaryTableS10_locus_1Mb_gene_panel.csv")

# ---------- S11: rs11635906 pQTL replication ----------
s11 = pd.read_csv(f"{SRC}/rs11635906_pqtl_replication.csv")
w(s11, "P2-02_SupplementaryTableS11_pQTL_replication.csv")

# ---------- S12: FDR-level scan hits coloc ----------
s12 = pd.read_csv(f"{SRC}/fdr_hit_coloc_results.csv")
w(s12, "P2-02_SupplementaryTableS12_FDR_hit_coloc.csv")

# ---------- S13: Finnish-matched LD sensitivity ----------
r = subprocess.run(["Rscript", "-e", '''
x <- readRDS("/workspace/fable_fix/pqtl/fin_ld_sensitivity.rds")
for (panel in c("EUR","FIN")) {
  d <- as.data.frame(x[[panel]]$summary); d$LD_panel <- panel
  write.csv(d, stdout(), row.names=FALSE)
}'''], capture_output=True, text=True)
blocks, cur = [], []
for line in r.stdout.splitlines():
    if line.startswith('"nsnps"') or line.startswith("nsnps"):
        if cur: blocks.append("\n".join(cur)); cur = []
    cur.append(line)
if cur: blocks.append("\n".join(cur))
s13 = pd.concat([pd.read_csv(io.StringIO(b)) for b in blocks], ignore_index=True)
s13 = s13.rename(columns={"hit1":"lead_eqtl_cs_variant","hit2":"lead_gwas_cs_variant","idx1":"eqtl_cs_index","idx2":"gwas_cs_index"})
s13["pair_note"] = np.where((s13.LD_panel=="EUR"), "single eQTL credible set x single GWAS credible set (primary result)",
                     np.where(s13["PP.H4.abf"]>0.5, "leading FIN eQTL credible set x GWAS credible set (primary result)",
                              "additional FIN-resolved eQTL substructure; supports distinct signal (PP.H3>=0.997)"))
meta = pd.DataFrame({"pair_note":["Per-variant LD correlation between FIN (n=99) and EUR (n=525) panels: r=0.925 over 2,013 shared variants"]})
s13 = pd.concat([s13, meta], ignore_index=True)
w(s13, "P2-02_SupplementaryTableS13_FIN_LD_sensitivity.csv")

# ---------- S14: extended FinnGen skin endpoints ----------
s14 = pd.read_csv(f"{SRC}/rs11635906_extended_skin_endpoints.csv")
s14["multiple_testing"] = "nominal only; none survives Bonferroni across the 7 extended endpoints (0.05/7=0.0071) except none; hypothesis-generating"
s14.loc[s14.phenocode=="L12_PSORIASIS","multiple_testing"] = "nominal (p=0.0026 passes 0.05/7=0.0071); hypothesis-generating, not scan-wide"
w(s14, "P2-02_SupplementaryTableS14_extended_skin_endpoints.csv")

# ---------- S15: regulatory annotation ----------
rows = [
 ("VEP consequence","Open Targets VEP (GRCh38)","rs11635906 most severe consequence","upstream_gene_variant (MODIFIER); no coding consequence","6,812 bp from nearest protein-coding ITPKA TSS; 72,150 bp from TYRO3 TSS; 57,195 bp from RPAP1"),
 ("L2G score","Open Targets locus-to-gene","Asthma association locus containing rs11635906","ITPKA 0.34 (rank 1); NDUFAF1 0.23 (rank 2)","only atopic-disease-relevant GWAS credible set available in OT for this variant"),
 ("L2G score","Open Targets locus-to-gene","Blood-cell trait loci (eosinophil/granulocyte/neutrophil/MCHC)","NDUFAF1 0.36-0.58 and RPAP1 dominate","consistent with blood-cell associations being non-ITPKA"),
 ("L2G score","Open Targets locus-to-gene","Strongest molQTL credible set at locus","NDUFAF1 (eQTL lead p=4.5e-119)","NDUFAF1 x eczema coloc PP.H3=0.9995 (distinct signal), concordant with Table 3"),
 ("Enhancer-gene prediction","ENCODE-rE2G / scE2G (1,458 cell/tissue types)","rs11635906","No predictions found","variant does not overlap any predicted enhancer-gene link; LD proxies not individually queryable"),
 ("Promoter-capture Hi-C","Javierre 2016 (17 blood cell types, CHiCAGO>=5)","PIR coverage of rs11635906","None","no promoter-interacting region covers the index variant"),
 ("Promoter-capture Hi-C","Javierre 2016","ITPKA promoter interactions","None in any of 17 cell types","consistent with low, T-cell-restricted ITPKA expression"),
 ("Promoter-capture Hi-C","Javierre 2016","TYRO3 promoter interactions","5 interactions, predominantly macrophage (Mac0/1/2)","PIRs at 41.99-42.11 Mb (hg19), 200-330 kb downstream of index; no overlap with GWAS signal region"),
 ("Promoter-capture Hi-C","Javierre 2016","Locus-wide interactions within +/-1 Mb","276 interactions","mainly MAPKBP1 / JMJD7-PLA2G4B / SPTBN5 / CHP1-EXD1 / NUSAP1-OIP5 baits; enriched in erythroid/megakaryocyte (EP/Ery/MK); full list below"),
]
s15a = pd.DataFrame(rows, columns=["component","resource","item","result","detail"])
pc = pd.read_csv(f"{SRC}/pchic_locus_interactions.csv")
cellcols = ["Mon","Mac0","Mac1","Mac2","Neu","MK","EP","Ery","FoeT","nCD4","tCD4","aCD4","naCD4","nCD8","tCD8","nB","tB"]
pc["celltypes_CHiCAGO_ge5"] = pc[cellcols].apply(lambda r: ";".join([c for c in cellcols if pd.notna(r[c]) and r[c]>=5]), axis=1)
pc["max_CHiCAGO"] = pc[cellcols].max(axis=1)
s15b = pd.DataFrame({
    "component":"PCHi-C interaction (Javierre 2016)",
    "resource":"PCHiC_peak_matrix_cutoff5",
    "item":pc["baitName"].astype(str)+" bait (chr15:"+pc["baitStart"].astype(str)+"-"+pc["baitEnd"].astype(str)+", hg19)",
    "result":"OE chr15:"+pc["oeStart"].astype(str)+"-"+pc["oeEnd"].astype(str)+"; max CHiCAGO="+pc["max_CHiCAGO"].round(2).astype(str),
    "detail":"cell types with CHiCAGO>=5: "+pc["celltypes_CHiCAGO_ge5"].replace("","none")})
s15 = pd.concat([s15a, s15b], ignore_index=True)
w(s15, "P2-02_SupplementaryTableS15_regulatory_annotation.csv")

# ---------- S16: ITPKA vs ITPKB ----------
s16 = pd.read_csv(f"{SRC}/census_itpka_vs_itpkb.csv")
s16["note"] = "CZ CELLxGENE Census detection fraction; scRNA-seq dropout makes low-expression detection a lower-bound estimate"
s16.loc[len(s16)] = {"tissue":"genetic","gene":"ITPKB","subset":"1q42 locus, +/-500 kb","n_cells":np.nan,"n_pos":np.nan,
    "pct_pos":np.nan,"mean_expr":np.nan,"mean_pos":np.nan,
    "note":"no genome-wide significant dermatitis/eczema or AD association at the ITPKB locus (FinnGen R12 minimum p=1.3e-6)"}
w(s16, "P2-02_SupplementaryTableS16_ITPKA_vs_ITPKB.csv")

# ---------- S17: IL6ST/ANKRD55 discrimination ----------
s17 = pd.read_csv(f"{SRC}/il6st_ankrd55_discrimination.csv")
s17["interpretation"] = "identical posteriors for the two genes: locus-level colocalization is gene-indistinguishable at 5q11.2"
w(s17, "P2-02_SupplementaryTableS17_IL6ST_ANKRD55_discrimination.csv")

# ---------- S8 rebuild: 9 traits x 3 endpoints + cystatin C attribution ----------
s8m = pd.read_csv(f"{SRC}/s8_coloc_abf_results.csv")
s8m.insert(0, "analysis", "coloc.abf")
cond = pd.read_csv(f"{SRC}/cystatinc_cond_coloc_results.csv")
label = {"unconditioned":"cystatin C x dermatitis/eczema GWAS, unconditioned (reference)",
         "cond_rs13329240":"cystatin C x dermatitis/eczema GWAS, both conditioned on rs13329240 (r=0.61 with rs11635906)",
         "cond_rs7165675":"cystatin C x dermatitis/eczema GWAS, both conditioned on rs7165675 (r=0.76 with rs11635906; over-conditioned - variant tags signal A itself)",
         "drop_region":"ITPKA eQTL x dermatitis/eczema GWAS after excluding chr15:41.15-41.35 Mb (cystatin C peak region)"}
cond["analysis"] = cond["analysis"].map(label)
cond = cond.rename(columns={"H3":"PP.H3","H4":"PP.H4"})
s8 = pd.concat([s8m, cond], ignore_index=True)
s8["note"] = ""
s8.loc[s8.analysis.str.contains("conditioned on rs7165675", na=False),"note"] = "abolishes colocalization by construction; uninformative for signal separation"
s8.loc[len(s8)] = {"analysis":"eGFR z-profile correlation across locus","trait":"cystatinC","note":"Pearson r=-0.76 between cystatin C and eGFR z-score profiles across the locus; supports kidney-function-related third signal"}
s8.loc[len(s8)] = {"analysis":"cluster LD with index variant","trait":"cystatinC","note":"rs13329240 r=0.612, rs7165675 r=0.757 with rs11635906 (1000 Genomes EUR); full pairwise LD in cystatinc_cluster_ld.csv (code repository)"}
w(s8, "P2-02_SupplementaryTableS8_PheWAS_trait_coloc.csv")

# ---------- S2 update ----------
s2 = pd.read_csv("/mnt/results/P2-02_SupplementaryTableS2_PheWAS_rs11635906.csv")
# 1. childhood asthma beta fix (OR 0.936 was stored as beta)
m = s2.id=="ebi-a-GCST007800"
assert abs(s2.loc[m,"beta"].iloc[0] - (-0.935996)) < 1e-6
s2.loc[m,"beta"] = round(np.log(0.936), 4)   # -0.0661
# 2. harmonization status column
def harm(r):
    if pd.isna(r["eaf"]): return "EAF not reported in source"
    if abs(r["eaf"]-0.262) < 0.02: return "concordant (EAF within 0.02 of eQTLGen 0.262)"
    return "EAF discordant (possible allele-coding or ancestry difference); direction interpreted with caution"
s2["harmonization_status"] = s2.apply(harm, axis=1)
# 3. notes
s2["note"] = ""
s2.loc[m,"note"] = "beta corrected from OR scale (OR=0.936) to log-OR; SE and p unchanged (log-scale)"
s2.loc[s2.id=="ebi-a-GCST003045","note"] = "de Lange 2017 UC; same statistics re-deposited as ieu-a-970 (duplicate, counted once in text)"
s2.loc[s2.id=="ieu-a-970","note"] = "duplicate deposition of de Lange 2017 (ebi-a-GCST003045) statistics"
s2.loc[s2.id=="ieu-a-968","note"] = "Liu 2015 UC (independent of de Lange 2017)"
sakaue = ["ebi-a-GCST90018959","ebi-a-GCST90018979","ebi-a-GCST90018953","ebi-a-GCST90018968"]
s2.loc[s2.id.isin(sakaue),"note"] = "Sakaue 2021; beta/SE as deposited (rounded to 2-4 significant digits in source)"
w(s2, "P2-02_SupplementaryTableS2_PheWAS_rs11635906.csv")

# ---------- S3 update: NGFRAP1 -> BEX3 ----------
s3 = pd.read_csv("/mnt/results/P2-02_SupplementaryTableS3_scan498.csv")
n = (s3.gene=="NGFRAP1").sum()
s3["pair_id"] = s3["pair_id"].str.replace("NGFRAP1","BEX3")
s3["gene"] = s3["gene"].replace("NGFRAP1","BEX3")
print(f"S3: renamed NGFRAP1 -> BEX3 in {n} rows (current HGNC symbol; ENSG00000166681 unchanged)")
w(s3, "P2-02_SupplementaryTableS3_scan498.csv")

# ---------- S1 update: units ----------
s1 = pd.read_csv("/mnt/results/P2-02_SupplementaryTableS1_dilution_calculations.csv")
s1 = s1[~s1.quantity.str.contains("Genetic effect diluted")]
new_rows = pd.DataFrame([
 {"quantity":"Expected bulk-skin per-allele expression shift (Census basis)","value":"~0.7% of the within-cell allelic effect",
  "note":"Mixture model: bulk expression = sum_c(f_c x E_c); per-allele bulk shift = f_ITPKA+ x dE where dE is the per-cell eQTL effect on the same (linear, normalized) expression scale. With f_ITPKA+ = 0.66%, the bulk shift is ~0.0066 x dE - of order 10^-2 to 10^-3 of the within-cell effect and below bulk-tissue detection limits. eQTLGen beta (SD units of normalized expression) is NOT multiplied by a cell fraction to yield log2FC; the two scales are not commensurable."},
 {"quantity":"Expected bulk-skin per-allele expression shift (He 2020 basis)","value":"~0.07% of the within-cell allelic effect",
  "note":"f_ITPKA+ = 0.067% in lesional AD skin; bulk shift ~0.0007 x dE - undetectable in bulk tissue"},
])
s1 = pd.concat([s1, new_rows], ignore_index=True)
w(s1, "P2-02_SupplementaryTableS1_dilution_calculations.csv")

print("\nAll tables written to", OUT)
