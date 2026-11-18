#!/bin/bash
# Descarga las 12 bibliotecas de hoja de PRJNA486291 con reanudacion y verificacion de tamano.
W="${NAC_ROOT:-.}/rnaseq_udec"
cd "$W" || exit 1
mkdir -p fastq logs

dl() {
  run="${1%%|*}"; url="${1##*|}"
  want=$(awk -v r="$run" '$1==r{print $2}' "$W/tamanos.tsv")
  f="$W/fastq/$run.fastq.gz"
  for i in 1 2 3 4 5 6 7 8 9 10; do
    have=$(stat -c%s "$f" 2>/dev/null || echo 0)
    if [ -n "$want" ] && [ "$have" = "$want" ]; then
      echo "OK   $run $have" >> "$W/logs/descarga_hoja.log"; return 0
    fi
    curl -sS -L -C - --speed-limit 200000 --speed-time 90 --retry 3 --retry-delay 5 \
         -o "$f" "$url" 2>> "$W/logs/curl.err"
  done
  echo "FALLO $run $(stat -c%s "$f" 2>/dev/null || echo 0)/$want" >> "$W/logs/descarga_hoja.log"
}
export -f dl; export W

awk '{print $1"|"$2}' "$W/urls_hoja.txt" | xargs -P 4 -I{} bash -c 'dl "$@"' _ {}
