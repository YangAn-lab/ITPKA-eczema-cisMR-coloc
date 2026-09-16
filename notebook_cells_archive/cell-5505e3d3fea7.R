JWT=$(cat /workspace/.opengwas_jwt)
echo "=== 1. GWAS Catalog: rs11635906 associations ==="
curl -sL --max-time 40 "https://www.ebi.ac.uk/gwas/rest/api/variants/rs11635906/associations?size=20" -H "Accept: application/json" | python3 -c "
import json,sys
raw=sys.stdin.read()
try:
    d=json.loads(raw)
    emb=d.get('_embedded',{}).get('associations',[])
    print('associations found:', len(emb) if emb else 0, '| totalElements:', d.get('page',{}).get('totalElements','?'))
except Exception as e:
    print('parse fail:', e, raw[:200])
"
echo "=== 2. OpenGWAS phewas (POST) ==="
curl -s --max-time 120 -X POST "https://api.opengwas.io/api/phewas" -H "Authorization: Bearer $JWT" -H "Content-Type: application/json" -d '{\"rsid\":[\"rs11635906\"],\"pval\":1e-5}' -o /workspace/coloc_oliva/phewas_opengwas.json -w "http=%{http_code} size=%{size_download}\n"
head -c 300 /workspace/coloc_oliva/phewas_opengwas.json; echo
echo "=== 3. FinnGen PheWeb variant API ==="
for v in "15-41487062-A-G" "15-41487062-G-A"; do
  curl -sL --max-time 40 "https://r12.finngen.fi/api/v1/variant/$v" -o /workspace/coloc_oliva/pheweb_$v.json -w "$v http=%{http_code} size=%{size_download}\n"
done