JWT=$(cat /workspace/.opengwas_jwt)
echo "=== OpenGWAS phewas (variant 参数) ==="
curl -s --max-time 180 -X POST "https://api.opengwas.io/api/phewas" -H "Authorization: Bearer $JWT" -H "Content-Type: application/json" -d '{\"variant\":[\"rs11635906\"],\"pval\":1e-5}' -o /workspace/coloc_oliva/phewas_opengwas.json -w "http=%{http_code} size=%{size_download}\n"
python3 -c "
import json
try:
    d=json.load(open('/workspace/coloc_oliva/phewas_opengwas.json'))
    print('records:', len(d) if isinstance(d,list) else d)
except Exception as e:
    print('parse:', e); print(open('/workspace/coloc_oliva/phewas_opengwas.json').read()[:300])
"