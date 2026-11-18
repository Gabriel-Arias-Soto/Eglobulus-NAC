#!/usr/bin/env bash
# 19_versiones.sh -> manuscrito/versiones_software.csv
#
# Metodos no puede llevar versiones tecleadas: se sondean los binarios y las
# librerias que efectivamente corrieron. Lo que no se encuentra se escribe como
# "no disponible" en vez de omitirse, para que la ausencia se vea en la tabla y
# no se confunda con un programa que no se uso.
#
# Uso, desde la raiz de "Revision NAC":   bash scripts/19_versiones.sh
set -uo pipefail

ENVS="${CONDA_ENVS:-${CONDA_ENVS:-$HOME/conda/envs}}"

# Version de un paquete de Python instalado en un entorno dado. Se lee de los
# metadatos de la distribucion, no de un atributo __version__: no todos los
# paquetes lo exponen y los metadatos son lo que registra el instalador.
vpy() {  # vpy <entorno> <paquete>
  "$ENVS/$1/bin/python" -c "import importlib.metadata as m;print(m.version('$2'))" \
    2>/dev/null || echo "no disponible"
}
SALIDA="manuscrito/versiones_software.csv"
mkdir -p manuscrito

# Primer token con forma de version de lo que imprima el binario. Un patron
# unico para todos: los programas anteponen su nombre de formas distintas
# ("HMMER 3.4", "trimAl v1.5.rev1", "IQ-TREE version 3.1.3") y ninguno lleva
# digitos en el nombre, asi que el primer numero de la linea es la version.
v() {  # v <entorno> <binario> <argumento de version>
  local p="$ENVS/$1/bin/$2" s
  [ -x "$p" ] || { echo "no disponible"; return; }
  s=$("$p" $3 2>&1 | tr -d '\r' | grep -oE '[0-9]+\.[0-9][0-9a-z.]*' | head -1)
  echo "${s:-no disponible}"
}

{
  echo "programa,version,paso,entorno"
  echo "HMMER,$(v hmmer hmmsearch -h),identificacion del dominio NAM,hmmer"
  echo "MAFFT,$(v filo mafft --version),alineamiento de proteinas,filo"
  echo "trimAl,$(v filo trimal --version),depuracion del alineamiento,filo"
  echo "IQ-TREE,$(v filo iqtree --version),inferencia filogenetica,filo"
  echo "MEME,$(v rnaseq meme -version),motivos C-terminales,rnaseq"
  echo "salmon,$(v rnaseq salmon --version),cuantificacion de RNA-seq,rnaseq"
  echo "PyMOL,$("$ENVS/struct/bin/python" -c 'import pymol;print(pymol.cmd.get_version()[0])' 2>/dev/null || echo "no disponible"),render de estructura,struct"

  # Capa estructural y de composicion. Faltaba entera aunque el Metodos la
  # cita: Boltz-2 predijo el complejo proteina-ADN sobre PyTorch en CPU, gemmi
  # lee las coordenadas y matplotlib compone F10.
  echo "Boltz-2,$(vpy boltz-cpu boltz),prediccion del complejo proteina-ADN,boltz-cpu"
  echo "PyTorch,$(vpy boltz-cpu torch),motor de calculo de Boltz-2,boltz-cpu"
  echo "gemmi,$(vpy struct gemmi),lectura de coordenadas y mmCIF,struct"
  echo "matplotlib,$(vpy struct matplotlib),composicion de la figura estructural,struct"
  echo "NumPy,$(vpy struct numpy),calculo numerico de la capa estructural,struct"
  echo "pandas,$(vpy struct pandas),tablas de la capa estructural,struct"
  # ColabFold no se instalo: se uso su SERVIDOR de alineamientos multiples, por
  # la opcion --use_msa_server de Boltz-2. No hay binario local que sondear y
  # declarar una version seria inventarla, asi que se declara el recurso.
  echo "ColabFold (servidor MSA),servicio remoto,alineamiento multiple para el plegado,ninguno"
  "$ENVS/r-nac/bin/Rscript" --vanilla -e '
    ps <- c(edgeR = "expresion diferencial", tximport = "lectura de la cuantificacion",
            ggtree = "figura de filogenia", treeio = "lectura de arboles",
            ComplexHeatmap = "mapas de calor", gggenes = "arquitectura y vecindades",
            ape = "manejo de arboles", ggplot2 = "figuras", patchwork = "paneles")
    for (p in names(ps)) cat(sprintf("%s,%s,%s,r-nac\n", p,
      tryCatch(as.character(packageVersion(p)), error = function(e) "no disponible"), ps[[p]]))
    cat(sprintf("R,%s,entorno de analisis,r-nac\n",
                paste(R.version$major, R.version$minor, sep = ".")))' 2>/dev/null
} > "$SALIDA"

printf "  %s: %d programas\n" "$SALIDA" "$(($(wc -l < "$SALIDA") - 1))"
FALTAN=$(grep -c "no disponible" "$SALIDA" || true)
[ "${FALTAN:-0}" -gt 0 ] && echo "  AVISO: $FALTAN sin version detectada"
exit 0
