# 2) ITPKA versioned gencodeId
gene = gtex_get("/reference/gene", geneId="ITPKA", datasetId="gtex_v8")
g = gene["data"][0]
print("gene:", g["geneSymbol"], g["gencodeId"], g["chromosome"], g["start"], "-", g["end"])
gencode_id = g["gencodeId"]

# 3) rs11635906 本身在两个皮肤组织是否为显著eQTL（对任何基因）
tissues = ["Skin_Sun_Exposed_Lower_leg", "Skin_Not_Sun_Exposed_Suprapubic"]
rows = []
for t in tissues:
    res = gtex_get("/association/singleTissueEqtl", variantId="chr15_41487062_A_G_b38",
                   tissueSiteDetailId=t, datasetId="gtex_v8")
    for d in res["data"]:
        rows.append({"tissue": t, "gene": d.get("geneSymbol"), "gencodeId": d.get("gencodeId"),
                     "nes": d.get("nes"), "p": d.get("pValue"), "variantId": d.get("variantId")})
df_snp = pd.DataFrame(rows)
print("\nrs11635906 作为显著eQTL的记录:")
print(df_snp.to_string(index=False) if len(df_snp) else "  （两个皮肤组织中均非显著eQTL）")