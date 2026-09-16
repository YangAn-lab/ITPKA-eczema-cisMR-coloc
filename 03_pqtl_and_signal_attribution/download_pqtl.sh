#!/bin/bash
set -e
BASE="https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics"
declare -A FILES=(
  [INTERVAL_ITPKA]="GCST90241001-GCST90242000/GCST90241504/harmonised/GCST90241504.h.tsv.gz"
  [INTERVAL_TYRO3]="GCST90243001-GCST90244000/GCST90243208/harmonised/GCST90243208.h.tsv.gz"
  [ICE_ITPKA]="GCST90087001-GCST90088000/GCST90087462/harmonised/35078996-GCST90087462-EFO_0007937.h.tsv.gz"
  [ICE_TYRO3]="GCST90087001-GCST90088000/GCST90087596/harmonised/35078996-GCST90087596-EFO_0007937.h.tsv.gz"
  [HELIC_TYRO3]="GCST90010001-GCST90011000/GCST90010356/harmonised/33303764-GCST90010356-EFO_0007937.h.tsv.gz"
)
# 按大小升序下载，逐个校验字节数
for K in INTERVAL_TYRO3 ICE_ITPKA INTERVAL_ITPKA ICE_TYRO3 HELIC_TYRO3; do
  URL="$BASE/${FILES[$K]}"
  OUT="${K}.tsv.gz"
  EXPECTED=$(curl -sI "$URL" | grep -i content-length | awk '{print $2}' | tr -d '\r')
  echo "[$(date +%H:%M:%S)] 下载 $K (预期 ${EXPECTED} 字节)..."
  curl -sL --retry 3 --retry-delay 10 "$URL" -o "$OUT"
  ACTUAL=$(stat -c%s "$OUT")
  if [ "$ACTUAL" != "$EXPECTED" ]; then
    echo "  !! 字节不符: $ACTUAL != $EXPECTED，重试一次"
    curl -sL --retry 3 "$URL" -o "$OUT"
    ACTUAL=$(stat -c%s "$OUT")
  fi
  echo "[$(date +%H:%M:%S)] $K 完成: $ACTUAL 字节 $([ "$ACTUAL" == "$EXPECTED" ] && echo OK || echo MISMATCH)"
done
echo "全部下载完成"
ls -la *.tsv.gz
