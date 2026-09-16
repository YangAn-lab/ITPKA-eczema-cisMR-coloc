import requests, json, pandas as pd

BASE = "https://gtexportal.org/api/v2"
S = requests.Session()

def gtex_get(path, **params):
    r = S.get(f"{BASE}{path}", params=params, timeout=60)
    r.raise_for_status()
    return r.json()

# 1) rs11635906 -> GTEx variant ID (hg38)
var = gtex_get("/dataset/variant", snpId="rs11635906", datasetId="gtex_v8")
print(json.dumps(var, indent=1)[:1500])