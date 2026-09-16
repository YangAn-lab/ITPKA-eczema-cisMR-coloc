import os

# 6) 尝试 dyneqtl 动态计算 rs11635906->ITPKA 在皮肤的名义关联（显著对列表只含过阈值对）
print("--- dyneqtl 名义关联（gtex_v8）---")
for t in skin:
    try:
        r = S.get(f"{BASE}/association/dyneqtl", timeout=60,
                  params={"variantId": VID, "gencodeId": "ENSG00000137825.10",
                          "tissueSiteDetailId": t, "datasetId": "gtex_v8"})
        print(t, r.status_code, r.text[:400])
    except Exception as e:
        print(t, "failed:", e)

# 7) 皮肤/血液 eQTL 样本量（功效背景）
print("\n--- 组织样本量 ---")
try:
    td = gtex_get("/dataset/tissueSiteDetail", datasetId="gtex_v8")
    for d in td["data"]:
        if d.get("tissueSiteDetailId") in skin + ["Whole_Blood"]:
            print(f"  {d['tissueSiteDetailId']:35s} n={d.get('rnaSeqAndGenotypeSampleCount', d.get('rnaSeqSampleCount'))}")
except Exception as e:
    print("tissueSiteDetail failed:", e)

# 8) OpenGWAS JWT 是否可用（只查存在性，不打印值）
jwt = os.environ.get("OPENGWAS_JWT", "")
print("\nOPENGWAS_JWT 已配置:", bool(jwt))
if jwt:
    for hdr in [{"X-API-TOKEN": jwt}, {"Authorization": f"Bearer {jwt}"}]:
        try:
            r = requests.get("https://api.opengwas.io/api/associations",
                             params={"variant": "rs11635906", "id": "eqtl-a-ENSG00000137825"},
                             headers=hdr, timeout=60)
            print(list(hdr)[0], "->", r.status_code)
            if r.status_code == 200:
                dat = r.json()
                print(json.dumps(dat, indent=1)[:800]); break
        except Exception as e:
            print("OpenGWAS query failed:", e)