#!/bin/bash
BASE="https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics"
declare -A FILES=(
  [INTERVAL_ITPKA]="GCST90241001-GCST90242000/GCST90241504/harmonised/GCST90241504.h.tsv.gz"
  [ICE_ITPKA]="GCST90087001-GCST90088000/GCST90087462/harmonised/35078996-GCST90087462-EFO_0007937.h.tsv.gz"
  [ICE_TYRO3]="GCST90087001-GCST90088000/GCST90087596/harmonised/35078996-GCST90087596-EFO_0007937.h.tsv.gz"
  [HELIC_TYRO3]="GCST90010001-GCST90011000/GCST90010356/harmonised/33303764-GCST90010356-EFO_0007937.h.tsv.gz"
)
for K in INTERVAL_ITPKA ICE_ITPKA ICE_TYRO3 HELIC_TYRO3; do
  URL="$BASE/${FILES[$K]}"
  EXPECTED=$(curl -sI "$URL" | grep -i content-length | awk '{print $2}' | tr -d '\r')
  echo "[$(date +%H:%M:%S)] $K 预期 $EXPECTED"
  aria2c -x 4 -s 4 -k 20M --max-tries=8 --retry-wait=10 --console-log-level=error --summary-interval=0 "$URL" -o "${K}.tsv.gz"
  ACTUAL=$(stat -c%s "${K}.tsv.gz" 2>/dev/null || echo 0)
  echo "[$(date +%H:%M:%S)] $K: $ACTUAL / $EXPECTED $([ "$ACTUAL" == "$EXPECTED" ] && echo OK || echo MISMATCH)"
done
echo "ALL DONE"
