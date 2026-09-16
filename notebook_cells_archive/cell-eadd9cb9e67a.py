# 11) ITPKA 在 hg19/GRCh37 的坐标（Ensembl GRCh37 REST）——核实移交文档“基因内”表述
r = requests.get("https://grch37.rest.ensembl.org/lookup/symbol/homo_sapiens/ITPKA",
                 headers={"Content-Type": "application/json"}, timeout=60)
g37 = r.json()
print("ITPKA GRCh37:", g37["seq_region_name"], g37["start"], "-", g37["end"], "strand:", g37["strand"])
var_b37 = 41779260
inside37 = g37["start"] <= var_b37 <= g37["end"]
print(f"rs11635906 (b37 pos {var_b37}) 在ITPKA基因内(hg19注释): {inside37}")
print(f"hg38 (GTEx v8): ITPKA 41,493,393-41,503,551; 变异 41,487,062 -> 距基因起点 {41493393-41487062} bp 上游")

# 12) 日晒皮肤 eQTL 样本量精确值
td = gtex_get("/dataset/tissueSiteDetail", datasetId="gtex_v8")
n_map = {d["tissueSiteDetailId"]: d["eqtlSampleSummary"]["totalCount"] for d in td["data"]
         if d.get("tissueSiteDetailId") in skin + ["Whole_Blood"]}
print("\neQTL样本量:", n_map)