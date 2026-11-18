#!/usr/bin/env bash
# Árbol de máxima verosimilitud de los dominios NAC con 185 E. globulus + 156 anclas de E. grandis (18-11-2026).
# Entrada: NAC_dom.trim_341.fa (revision_arbol_185.py). Mismos parámetros que el árbol de 183 (08-09-2026):
# ModelFinder restringido a modelos estándar (-m TEST), 1000 réplicas UFBoot y 1000 SH-aLRT. IQ-TREE 3.1.3.
# KMP_AFFINITY=disabled evita que OpenMP aborte al fijar afinidad de hilos en WSL2.
# La búsqueda se puede cortar y retomar: IQ-TREE continúa desde arbol_NAC_189_v2.ckp.gz con el mismo comando.
set -euo pipefail
cd "$(dirname "$0")"
export KMP_AFFINITY=disabled OMP_NUM_THREADS="${HILOS:-10}"
iqtree3 -s NAC_dom.trim_341.fa -m TEST -bb 1000 -alrt 1000 -nt "${HILOS:-10}" -seed 20261118 --prefix arbol_NAC_189_v2 --quiet
grep -E "Best-fit model|Log-likelihood of consensus" arbol_NAC_189_v2.log
# después: cp arbol_NAC_189_v2.treefile ../arbol_NAC_189.treefile y python scripts/20_subfamilias_Hussey.py --escribir
