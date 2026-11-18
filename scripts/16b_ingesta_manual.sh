#!/usr/bin/env bash
# 16b_ingesta_manual.sh — recoge los FASTQ bajados a mano, los verifica y los
# deja listos para 16_quant_hoja.sh.
#
# Por que existe: la descarga automatica desde este entorno iba a ~1,1 MB/s por
# conexion y las 42,4 GB restantes se bajaron desde el equipo del usuario con un
# gestor. Este guion hace la unica parte que no se puede delegar: comprobar que
# lo recibido es exactamente lo que ENA declara, antes de que entre en la
# cuantificacion. Un fichero con el tamano correcto puede llevar un tramo
# duplicado dentro y dar valores de expresion plausibles y falsos.
#
# Que hace, por fichero de rnaseq_udec/fastq_manual/:
#   1. contrasta el tamano con fastq_bytes de ENA
#   2. contrasta el md5 con fastq_md5 de ENA
#   3. si ambos casan, lo mueve a rnaseq_udec/fastq/ y escribe el testigo
#      .md5ok, que es lo que 16_quant_hoja.sh mira para no volver a bajarlo
#   4. si falla, lo DEJA donde esta y lo dice: no borra nada bajado a mano
#
# Uso:  bash scripts/16b_ingesta_manual.sh
set -uo pipefail

RAIZ="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BASE="$RAIZ/rnaseq_udec"
ORIGEN="$BASE/fastq_manual"
DESTINO="$BASE/fastq"
REF="$BASE/ena_md5.tsv"

[[ -s "$REF" ]] || { echo "falta $REF" >&2; exit 1; }
mkdir -p "$DESTINO" "$BASE/logs"

tam_ena() { awk -F'\t' -v s="$1" '$1==s{print $2}' "$REF"; }
md5_ena() { awk -F'\t' -v s="$1" '$1==s{print $3}' "$REF"; }

ok=0; mal=0; ausente=0
printf '%-12s %-10s %s\n' RUN ESTADO DETALLE

while IFS=$'\t' read -r s md5esp tamesp _resto; do
  [[ "$s" == "run" ]] && continue
  f="$ORIGEN/$s.fastq.gz"

  if [[ ! -f "$f" ]]; then
    printf '%-12s %-10s %s\n' "$s" ausente "no esta en fastq_manual/"
    ausente=$((ausente+1)); continue
  fi

  hay=$(stat -c%s "$f")
  if [[ "$hay" -ne "$tamesp" ]]; then
    printf '%-12s %-10s %s\n' "$s" TAMANO "tiene $hay, espera $tamesp (faltan $((tamesp-hay)))"
    mal=$((mal+1)); continue
  fi

  obt=$(md5sum "$f" | cut -d' ' -f1)
  if [[ "$obt" != "$md5esp" ]]; then
    printf '%-12s %-10s %s\n' "$s" MD5 "obtenido $obt"
    mal=$((mal+1)); continue
  fi

  mv -f "$f" "$DESTINO/$s.fastq.gz"
  echo "verificado $(date -Iseconds) (descarga manual)" > "$DESTINO/$s.fastq.gz.md5ok"
  printf '%-12s %-10s %s\n' "$s" ok "movido y verificado"
  ok=$((ok+1))
done < "$BASE/PENDIENTES_descarga.tsv"

echo
echo "verificados $ok · fallidos $mal · sin recibir $ausente"
if [[ "$mal" -gt 0 ]]; then
  echo "Los fallidos siguen en fastq_manual/: hay que volver a bajarlos enteros." >&2
fi
if [[ "$ok" -gt 0 && "$mal" -eq 0 && "$ausente" -eq 0 ]]; then
  echo "Todo correcto. Siguiente paso:  bash scripts/16_quant_hoja.sh"
fi
