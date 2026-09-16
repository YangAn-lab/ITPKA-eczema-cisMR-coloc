# 9) METASOFT 多组织 meta 分析：rs11635906-ITPKA 跨组织证据（含皮肤m-value）
print("--- multitissueEqtl (METASOFT, gtex_v8) ---")
try:
    r = S.get(f"{BASE}/association/multitissueEqtl", timeout=60,
              params={"variantId": VID, "gencodeId": "ENSG00000137825.10", "datasetId": "gtex_v8"})
    print(r.status_code)
    mt = r.json()
    d = mt["data"]
    if isinstance(d, list) and d:
        d0 = d[0]
        print("keys:", list(d0.keys()))
        print(json.dumps(d0, indent=1)[:1200])
    else:
        print(str(mt)[:600])
except Exception as e:
    print("failed:", e)

# 10) 组织样本量（看原始字段名）
print("\n--- tissueSiteDetail 原始字段 ---")
td = gtex_get("/dataset/tissueSiteDetail", datasetId="gtex_v8")
for d in td["data"]:
    if d.get("tissueSiteDetailId") in skin + ["Whole_Blood"]:
        print({k: v for k, v in d.items() if k != "tissueSiteDetail"})