VID = "chr15_41487062_A_G_b38"
skin = ["Skin_Sun_Exposed_Lower_leg", "Skin_Not_Sun_Exposed_Suprapubic"]

# 4) ITPKA 在两个皮肤组织中的全部显著eQTL（v8 与 v10 都查）+ rs11635906 在血液的交叉验证
for ds in ["gtex_v8", "gtex_v10"]:
    print("="*70); print("dataset:", ds)
    try:
        gid = gtex_get("/reference/gene", geneId="ITPKA", datasetId=ds)["data"][0]["gencodeId"]
    except Exception as e:
        print("  gene resolve failed:", e); continue
    for t in skin:
        try:
            res = gtex_get("/association/singleTissueEqtl", gencodeId=gid, tissueSiteDetailId=t, datasetId=ds)
            d = res["data"]
            if d:
                df = pd.DataFrame(d)[["variantId","snpId","nes","pValue"]].sort_values("pValue")
                print(f"  ITPKA显著eQTL @ {t}: n={len(df)}, top={df.iloc[0]['snpId']} nes={df.iloc[0]['nes']:.3f} p={df.iloc[0]['pValue']:.2e}")
                if VID in set(df['variantId']): print("    -> rs11635906 在显著列表中!")
            else:
                print(f"  ITPKA显著eQTL @ {t}: 无（该组织无ITPKA显著eQTL）")
        except Exception as e:
            print(f"  {t} query failed: {e}")
    # rs11635906 在全血（对照eQTLGen血液发现）
    try:
        res = gtex_get("/association/singleTissueEqtl", variantId=VID, tissueSiteDetailId="Whole_Blood", datasetId=ds)
        d = res["data"]
        genes = {x.get("geneSymbol"): (round(x.get("nes"),3), x.get("pValue")) for x in d}
        print("  rs11635906 全血显著eQTL基因:", genes if genes else "无")
    except Exception as e:
        print("  Whole_Blood query failed:", e)

# 5) ITPKA 在皮肤/血液的中位表达（判断皮肤eQTL缺失是否为低表达/功效问题）
expr = gtex_get("/expression/medianGeneExpression", gencodeId="ENSG00000137825.10",
                tissueSiteDetailId=skin+["Whole_Blood"], datasetId="gtex_v8")
print("\nITPKA 中位表达 (TPM, gtex_v8):")
for d in expr["data"]:
    print(f"  {d['tissueSiteDetailId']:35s} {d['median']:.3f}")