JWT=$(cat /workspace/.opengwas_jwt)
echo "=== OpenGWAS phewas GET query ==="
curl -s --max-time 180 "https://api.opengwas.io/api/phewas?variant=rs11635906&pval=0.00001" -H "Authorization: Bearer $JWT" -o /workspace/coloc_oliva/phewas_opengwas.json -w "http=%{http_code} size=%{size_download}\n"
head -c 200 /workspace/coloc_oliva/phewas_opengwas.json; echo
echo "=== FinnGen autocomplete ==="
curl -sL --max-time 40 "https://r12.finngen.fi/api/v1/autocomplete?query=rs11635906" | head -c 500; echo
echo "=== FinnGen variant API v2 ==="
curl -sL --max-time 40 "https://r12.finngen.fi/api/v1/variant/rs11635906" -w "\nhttp=%{http_code}\n" | head -c 300