JWT=$(cat /workspace/.opengwas_jwt)
echo "=== POST form-encoded ==="
curl -s --max-time 180 -X POST "https://api.opengwas.io/api/phewas" -H "Authorization: Bearer $JWT" -H "Content-Type: application/x-www-form-urlencoded" -d "variant=rs11635906&pval=0.00001" -o /workspace/coloc_oliva/phewas_opengwas.json -w "http=%{http_code} size=%{size_download}\n"
head -c 300 /workspace/coloc_oliva/phewas_opengwas.json; echo
echo "=== POST JSON 字符串非数组 ==="
curl -s --max-time 60 -X POST "https://api.opengwas.io/api/phewas" -H "Authorization: Bearer $JWT" -H "Content-Type: application/json" -d '{\"variant\":\"rs11635906\",\"pval\":0.00001}' -w "http=%{http_code}\n" | head -c 300