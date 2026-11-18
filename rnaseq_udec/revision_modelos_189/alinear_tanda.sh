#!/usr/bin/env bash
# Una tanda de alineamiento contra el GENOMA COMPLETO (no contra los loci), en primer plano.
# Uso: bash alinear_tanda.sh SRR saltar_lecturas n_lecturas etiqueta [hilos]
# (copia de alineamiento_parciales; aquí regiones_5kb.bed = 189 NAC ±5 kb)
# Guarda TODOS los alineamientos que caen en regiones_5kb.bed (únicos y múltiples; el filtro NH:i:1 se aplica
# en el análisis y los múltiples quedan para diagnóstico) y el resumen de HISAT2 del genoma completo.
set -euo pipefail
SRR="$1"; SALTA="$2"; N="$3"; ET="$4"; HILOS="${5:-8}"
BIN="${CONDA_ENVS:-$HOME/conda/envs}/rnaseq/bin"; export PATH="$BIN:$PATH"
cd "$(dirname "$0")"
IDX=../hisat2_idx/genoma; SS=../hisat2_idx/uniones_anotadas.txt; FQ=../fastq/$SRR.fastq.gz
OUT=tandas/${SRR}_${ET}
t0=$(date +%s)
# Lectura: pigz descomprime, tail salta 4*SALTA líneas y head corta 4*N (más rápido que hisat2 -s).
# El total leído se comprueba después contra "Total reads" del resumen de HISAT2.
{ pigz -dc -p 2 "$FQ" | tail -n +$((4*SALTA+1)) | head -n $((4*N)) || true; } \
 | hisat2 -p "$HILOS" -x "$IDX" --known-splicesite-infile "$SS" -U - \
       --new-summary --summary-file "$OUT.resumen.txt" --rg-id "$SRR" --rg "SM:$SRR" 2> "$OUT.hisat2.log" \
 | samtools view -u -L regiones_5kb.bed - \
 | samtools sort -@ 2 -m 512M -o "$OUT.bam" -
samtools index "$OUT.bam"
t1=$(date +%s)
nreg=$(samtools view -c "$OUT.bam"); nuni=$(samtools view -c -e '[NH]==1' -F 0x904 "$OUT.bam")
echo -e "$SRR\t$ET\t$SALTA\t$N\t$((t1-t0))\t$nreg\t$nuni" | tee -a tandas/registro_tandas.tsv
