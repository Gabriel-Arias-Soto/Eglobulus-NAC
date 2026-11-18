#!/bin/bash
# Descarga verificada por MD5 de las bibliotecas de hoja de PRJNA486291.
# Reanuda lo parcial, reinicia desde cero lo que exceda el tamano esperado o falle el MD5.
W="${NAC_ROOT:-.}/rnaseq_udec"
cd "$W" || exit 1
mkdir -p fastq logs
LOG="$W/logs/descarga_v3.log"

dl() {
  run="$1"
  want=$(awk -v r="$run" '$1==r{print $2}' "$W/ena_md5.tsv")
  md5w=$(awk -v r="$run" '$1==r{print $3}' "$W/ena_md5.tsv")
  url="https://$(awk -v r="$run" '$1==r{print $4}' "$W/ena_md5.tsv")"
  f="$W/fastq/$run.fastq.gz"
  [ -z "$want" ] && { echo "SINMETA $run" >> "$LOG"; return 1; }

  for intento in 1 2 3; do
    have=$(stat -c%s "$f" 2>/dev/null || echo 0)
    # exceso de tamano = fichero corrupto: vaciar y empezar de nuevo
    if [ "$have" -gt "$want" ]; then
      echo "EXCESO $run $have>$want, reinicio" >> "$LOG"; : > "$f"; have=0
    fi
    if [ "$have" -lt "$want" ]; then
      curl -sS -f -L -C - --speed-limit 150000 --speed-time 120 --retry 5 --retry-delay 10 \
           -o "$f" "$url" 2>> "$W/logs/curl_v3.err"
    fi
    have=$(stat -c%s "$f" 2>/dev/null || echo 0)
    [ "$have" != "$want" ] && continue
    md5g=$(md5sum "$f" | cut -d' ' -f1)
    if [ "$md5g" = "$md5w" ]; then
      echo "OK $run $have md5ok" >> "$LOG"; return 0
    fi
    echo "MD5MAL $run intento $intento, reinicio" >> "$LOG"; : > "$f"
  done
  echo "FALLO $run" >> "$LOG"; return 1
}
export -f dl; export W LOG

cut -f1 "$W/urls_hoja.txt" | xargs -P 4 -I{} bash -c 'dl "$@"' _ {}
echo "=== fin tanda $(date +%H:%M) ===" >> "$LOG"
