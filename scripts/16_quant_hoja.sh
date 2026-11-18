#!/usr/bin/env bash
# 16_quant_hoja.sh -> cuantificacion de las 12 librerias de HOJA de PRJNA486291
# contra el genoma conspecifico GCA_046226855.1 (E. globulus), no contra
# E. grandis. El motivo esta en el memo: con la referencia de E. grandis la
# tabla de DE queda en loci Eucgr y hay que cruzarla con el censo por el
# NOMBRE del simbolo, y el mismo texto nombra loci distintos en cada
# nomenclatura. Cuantificando contra E. globulus la tabla sale en ACJRO7_*,
# que es el identificador del censo Y el del transcrito que amplifican los
# cebadores de qPCR (t2g.tsv), asi que no queda ningun cruce por nombre.
#
# Se ejecuta desde la RAIZ del proyecto, en el disco Z:
#   bash scripts/16_quant_hoja.sh
# Es reanudable: salta las librerias que ya tienen quant.sf.
set -euo pipefail

BASE="rnaseq_udec"
IDX="$BASE/salmon_idx"
# El registro lo abre el script, no quien lo lanza: asi queda constancia
# aunque se ejecute en segundo plano o con la salida redirigida a otro sitio.
# Se anade, no se trunca, porque el guion es reanudable y el historial de
# reanudaciones es justamente lo que hay que poder releer.
mkdir -p "$BASE/logs"
LOG="$BASE/logs/pipeline_hoja.log"
exec > >(tee -a "$LOG") 2>&1
echo "=== $(date '+%F %T') : arranca 16_quant_hoja.sh (PID $$) ==="
HILOS=6          # los mismos que las dos librerias ya cuantificadas
EN_PARALELO=2    # 2 x 6 = 12 nucleos

# Las 12 de hoja segun los metadatos de ENA (sample_alias), 4 condiciones x 3
# replicas. Las de raiz quedan fuera: el mapa de calor y la correlacion con
# qPCR son solo de hoja.
HOJA=(SRR7726470 SRR7726471 SRR7726474      # NA   rep 1,2,3
      SRR7726478 SRR7726481 SRR7726475      # CABF rep 1,2,3
      SRR7726479 SRR7726483 SRR7726476      # CAAF rep 1,2,3
      SRR9887586 SRR9888129 SRR9888508)     # DS   rep 1,2,3

mkdir -p "$BASE/quant" "$BASE/logs" "$BASE/fastq"

# --- fase 1: integridad de los fastq ----------------------------------------
# De las 12 descargas previas 7 estaban INCOMPLETAS (de un 7 % a un 85 % del
# fichero), y salmon no siempre avisa: con el flujo gzip cortado limpio
# cuantifica un subconjunto de lecturas y devuelve exito. Las tres replicas del
# control NA fallaron, pero no con el mismo mensaje: SRR7726471 y SRR7726474 con
# "incomplete deflate stream" y SRR7726470 con "Invalid header: (C): expected (@)".
# El diagnostico es el mismo en las tres -descarga truncada, confirmado por el
# cotejo de tamano contra ENA-; lo que cambia es donde cayo el corte. Asi que aqui se compara
# el tamano contra fastq_bytes de ENA, se reanuda lo que falte (curl -C -, que
# ahorra lo ya bajado) y se comprueba el md5. Sin este filtro la tabla de DE
# sale de librerias truncadas de distinto modo en cada grupo.
tam_ena() { awk -F'\t' -v s="$1" '$1==s{print $2}' "$BASE/ena_md5.tsv"; }
md5_ena() { awk -F'\t' -v s="$1" '$1==s{print $3}' "$BASE/ena_md5.tsv"; }
url_ena() { awk -F'\t' -v s="$1" '$1==s{print $4}' "$BASE/ena_md5.tsv"; }

descargar() {
  local s="$1" fq="$BASE/fastq/$1.fastq.gz"
  local esp; esp=$(tam_ena "$s")
  # Tres reglas que este bucle aprendio a golpes, cada una por un fallo real:
  #
  # 1. --speed-limit/--speed-time: las conexiones de ENA se quedan colgadas sin
  #    cerrar el socket. Sin corte, curl espera indefinidamente y la descarga se
  #    para sin que nadie la reintente.
  #
  # 2. NO se usa el --retry de curl junto con -C -. curl fija el desplazamiento
  #    de reanudacion UNA vez, al arrancar, y en su reintento interno vuelve a
  #    pedir desde ese mismo desplazamiento mientras sigue ANADIENDO al fichero:
  #    duplica un tramo y el fichero acaba pesando de mas. Quien reintenta es el
  #    bucle de aqui abajo, que invoca curl de nuevo y por tanto recalcula el
  #    desplazamiento contra el tamano real del fichero.
  #
  # 3. El md5 solo se comprueba cuando el fichero esta COMPLETO. Comprobarlo
  #    sobre un fichero a medias siempre falla, y la rama de fallo borra lo ya
  #    bajado: un corte de red a mitad costaba varios GB (le paso a SRR7726474).
  #    Incompleto = seguir reanudando; pasado de tamano = corrupto, se reinicia.
  for intento in 1 2 3 4 5 6; do
    local hay=0; [[ -f "$fq" ]] && hay=$(stat -c%s "$fq")
    if [[ "$hay" -gt "$esp" ]]; then
      echo "  $s pesa $hay > $esp esperado: corrupto, se reinicia" >&2
      rm -f "$fq"; hay=0
    fi
    if [[ "$hay" -lt "$esp" ]]; then
      echo "  bajando $s (tiene $hay de $esp, intento $intento)"
      # -f es obligatorio: sin el, curl escribe el CUERPO de las respuestas HTTP
      # de error dentro del .fastq.gz, y como venimos reanudando con -C - ese
      # cuerpo queda anexado al final. Ya ocurrio en este mismo proyecto
      # (SRR7726475 acabo pesando 1.405.195.701 frente a 1.356.084.661).
      curl -fsS -C - --speed-limit 50000 --speed-time 60 \
           --connect-timeout 30 -o "$fq" "https://$(url_ena "$s")" || true
      hay=0; [[ -f "$fq" ]] && hay=$(stat -c%s "$fq")
    fi
    if [[ "$hay" -eq "$esp" ]]; then
      if [[ "$(md5sum "$fq" | cut -d' ' -f1)" == "$(md5_ena "$s")" ]]; then
        # El testigo evita re-calcular el md5 de 8 GB en cada relanzamiento, y
        # sobre todo deja constancia EN DISCO de que este fichero concreto se
        # verifico: sin el, un fichero corrupto que da la talla se cuela.
        echo "verificado $(date -Iseconds)" > "$fq.md5ok"
        echo "  md5 ok  $s"; return 0
      fi
      echo "  MD5 MAL $s con tamano correcto: se descarta y se reinicia" >&2
      rm -f "$fq" "$fq.md5ok"
    fi
  done
  echo "NO SE PUDO DESCARGAR $s" >&2; return 1
}
export -f descargar tam_ena md5_ena url_ena; export BASE

# Un fichero entra en la cola si le falta tamano O si nunca se le comprobo el
# md5. Lo segundo no es redundante: el fallo real de este proyecto fue un
# fichero que alcanzo el tamano EXACTO y aun asi tenia un tramo duplicado
# dentro. Seleccionar solo por tamano lo habria mandado a cuantificar.
# Las ya cuantificadas se saltan porque salmon dejo num_processed en su
# meta_info.json y ese recuento se contrasto contra read_count de ENA, que es
# una comprobacion de integridad tan buena como el md5.
bajar=()
for s in "${HOJA[@]}"; do
  [[ -s "$BASE/quant/$s/quant.sf" ]] && continue
  fq="$BASE/fastq/$s.fastq.gz"
  esp=$(tam_ena "$s"); hay=0; [[ -f "$fq" ]] && hay=$(stat -c%s "$fq")
  if [[ "$hay" -ne "$esp" || ! -f "$fq.md5ok" ]]; then bajar+=("$s"); fi
done
if [[ ${#bajar[@]} -gt 0 ]]; then
  echo "== fase 1: ${#bajar[@]} ficheros incompletos, reanudando =="
  # Una conexion por biblioteca. Medido en este equipo: ENA entrega ~1,16 MB/s
  # POR CONEXION y una conexion nueva abierta con otras seis ya en curso
  # consigue su porcion integra sin quitarsela a las demas, luego el limite es
  # del servidor por flujo y no del cliente. El disco (Z:, 183 MB/s de
  # escritura) no interviene. Con 12 el techo son ~14 MB/s.
  printf '%s\n' "${bajar[@]}" | xargs -P 12 -I{} bash -c 'descargar {}'
else
  echo "== fase 1: los 12 fastq estan completos =="
fi

# --- fase 2: cuantificacion --------------------------------------------------
pendientes=()
for s in "${HOJA[@]}"; do
  if [[ -s "$BASE/quant/$s/quant.sf" ]]; then
    echo "ya esta  $s"
  elif [[ ! -s "$BASE/fastq/$s.fastq.gz" ]]; then
    echo "SIN FASTQ $s" >&2; exit 1
  else
    pendientes+=("$s")
  fi
done
echo "== fase 2: por cuantificar ${#pendientes[@]} =="
[[ ${#pendientes[@]} -eq 0 ]] && { echo "== nada que hacer =="; exit 0; }

cuantificar() {
  local s="$1"
  # libType A: la deteccion automatica dio SR (single-end, stranded reverse)
  # en las dos primeras; se deja automatica para que quede registrado por
  # libreria en lib_format_counts.json en vez de imponerlo.
  salmon quant -i "$IDX" -l A -r "$BASE/fastq/$s.fastq.gz" \
    -o "$BASE/quant/$s" -p "$HILOS" > "$BASE/logs/quant_$s.log" 2>&1
  echo "listo $s"
}
export -f cuantificar; export BASE IDX HILOS

printf '%s\n' "${pendientes[@]}" | xargs -P "$EN_PARALELO" -I{} bash -c 'cuantificar {}'
echo "== cuantificacion terminada =="
