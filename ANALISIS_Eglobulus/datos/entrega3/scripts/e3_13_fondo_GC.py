"""Entrega 3, paso «Contraste con fondo emparejado por GC» (preparación): conjuntos de fondo.

Se construyen N_CONJ conjuntos aleatorios de 189 promotores del genoma (no NAC) con la misma
distribución de GC que los 189 promotores NAC. En cada conjunto, a cada promotor NAC se le asigna
un promotor del fondo del mismo intervalo de GC (anchura ANCHO_GC), sin reemplazo dentro del
conjunto. Si un intervalo tiene menos de MIN_RAZON promotores de fondo por promotor NAC, se fusiona
con el intervalo vecino hacia la mediana del fondo. Semilla fija.

Este paso no depende de la salida de PlantCARE: los sitios se cuentan después en estos conjuntos,
con las secuencias de sitio que devuelva PlantCARE.

Controles: GC recalculado del FASTA = columna GC de la tabla; promotores del FASTA de los 189 =
los del FASTA genómico (por locus); fondo sin NAC; cada promotor elegido cae en el intervalo de
su NAC; sin repetidos dentro de un conjunto; KS del GC (NAC frente a fondo elegido y frente al fondo
completo); reutilización de promotores entre conjuntos.

Entradas:
  datos/entrega3/promotores/promotores_2kb_todos.{fa,tsv}, NAC189_promotores_2kb.fa
  tablas/T40_promotores_NAC_2kb.csv (símbolo <-> locus)
Salidas (datos/entrega3/cis_plantcare/fondo/):
  fondo_GC_conjuntos.npz   nac (símbolos, orden del FASTA de 189), locus_fondo (reserva de fondo),
                           idx int32 [N_CONJ, 189]: idx[s, j] = posición en locus_fondo del
                           promotor emparejado con el NAC j en el conjunto s; gc_nac, gc_fondo
  fondo_GC_intervalos.tsv  intervalos finales: GC inicial y final, n NAC, n fondo, razón
  resumen_fondo_GC.json    parámetros, versiones y controles
Uso (desde la raíz «Revisión NAC», con el python de nac3):  python e3_13_fondo_GC.py [--n N]
  N = 1000 por defecto (colección original). Con otro N las salidas llevan el sufijo _N
  (fondo_GC_conjuntos_N.npz, resumen_fondo_GC_N.json), el archivo de intervalos no se reescribe (se
  comprueba que coincide) y se verifica que los primeros 1000 conjuntos son los de la colección original
  (misma semilla y mismo orden de muestreo: la colección de N conjuntos amplía la de 1000).
"""
import argparse, json, os, platform, time
import numpy as np
import pandas as pd
import scipy
from scipy.stats import ks_2samp

BASE = (os.environ.get("NAC_ROOT", ".") + "")
A = os.path.join(BASE, "ANALISIS_Eglobulus")
E3 = os.path.join(A, "datos", "entrega3")
PROM = os.path.join(E3, "promotores")
SAL = os.path.join(E3, "cis_plantcare", "fondo")
_ap = argparse.ArgumentParser()
_ap.add_argument("--n", type=int, default=1000, help="número de conjuntos (1000 = colección original)")
N_CONJ = _ap.parse_args().n
SUF = "" if N_CONJ == 1000 else f"_{N_CONJ}"
ANCHO_GC = 0.01
MIN_RAZON = 20
SEMILLA = 20261018
NB = round(1 / ANCHO_GC)


def leer_fasta(p):
    seqs, h, buf = {}, None, []
    with open(p) as fh:
        for l in fh:
            if l.startswith(">"):
                if h is not None:
                    seqs[h] = "".join(buf).upper()
                h, buf = l[1:].split()[0], []
            else:
                buf.append(l.strip())
    if h is not None:
        seqs[h] = "".join(buf).upper()
    return seqs


def n_gc(s):
    return s.count("G") + s.count("C")


t0 = time.time()
os.makedirs(SAL, exist_ok=True)
tab = pd.read_csv(os.path.join(PROM, "promotores_2kb_todos.tsv"), sep="\t")
gen = leer_fasta(os.path.join(PROM, "promotores_2kb_todos.fa"))
nacf = leer_fasta(os.path.join(PROM, "NAC189_promotores_2kb.fa"))
t40 = pd.read_csv(os.path.join(A, "tablas", "T40_promotores_NAC_2kb.csv"))
print(f"tabla {len(tab)} filas; FASTA genómico {len(gen)}; FASTA NAC {len(nacf)}; T40 {len(t40)}")

# Control 1: tabla y FASTA genómico describen las mismas secuencias, con el mismo GC
assert tab.locus.is_unique and set(tab.locus) == set(gen), "loci de la tabla != FASTA genómico"
assert not any(set(s) - set("ACGT") for s in gen.values()), "hay caracteres distintos de ACGT"
tab["largo_fa"] = tab.locus.map(lambda x: len(gen[x]))
tab["ngc"] = tab.locus.map(lambda x: n_gc(gen[x]))
dif_gc = (tab.ngc / tab.largo_fa - tab.GC).abs().max()
assert (tab.largo_fa == tab.largo).all() and dif_gc < 6e-5, f"GC/largo no coinciden (máx. dif. {dif_gc})"
print(f"control 1: largo y GC del FASTA = tabla (máx. |dif. GC| {dif_gc:.1e}; GC de la tabla con 4 decimales)")

# Control 2: los 189 NAC (símbolo -> locus) y su secuencia
sim2loc = dict(zip(t40.simbolo, t40.locus))
nac = list(nacf)
assert len(nac) == 189 and set(nac) == set(t40.simbolo), "símbolos del FASTA NAC != T40"
assert set(t40.locus) == set(tab.locus[tab.es_NAC]), "loci NAC de T40 != es_NAC de la tabla"
iguales = sum(nacf[s] == gen[sim2loc[s]] for s in nac)
assert iguales == 189, f"solo {iguales}/189 promotores NAC idénticos al FASTA genómico"
print("control 2: 189/189 promotores NAC idénticos a su locus en el FASTA genómico")

# Intervalos de GC con aritmética entera (sin errores de coma flotante en los bordes)
tab["bin"] = (tab.ngc * NB) // tab.largo_fa
fondo = tab[~tab.es_NAC].reset_index(drop=True)
assert len(fondo) == len(tab) - 189 and not set(fondo.locus) & set(t40.locus)
loc2fila = dict(zip(tab.locus, range(len(tab))))
bin_nac = np.array([tab.bin.iat[loc2fila[sim2loc[s]]] for s in nac])
gc_nac = np.array([tab.ngc.iat[loc2fila[sim2loc[s]]] / tab.largo_fa.iat[loc2fila[sim2loc[s]]] for s in nac])
bin_fondo = fondo.bin.to_numpy()
gc_fondo = (fondo.ngc / fondo.largo_fa).to_numpy()
cn = np.bincount(bin_nac, minlength=NB + 1)
cf = np.bincount(bin_fondo, minlength=NB + 1)
med = int(np.median(bin_fondo))

# Fusión de intervalos con poco fondo, hacia la mediana del fondo
grupos = [[b, b] for b in range(NB + 1)]
def cuenta(g):
    return cn[g[0]:g[1] + 1].sum(), cf[g[0]:g[1] + 1].sum()
fusiones = 0
while True:
    malos = [(cuenta(g)[1] / cuenta(g)[0], i) for i, g in enumerate(grupos) if cuenta(g)[0] > 0 and cuenta(g)[1] < MIN_RAZON * cuenta(g)[0]]
    if not malos:
        break
    _, i = min(malos)
    if grupos[i][1] < med:
        j = i + 1
    elif grupos[i][0] > med:
        j = i - 1
    else:
        j = max((k for k in (i - 1, i + 1) if 0 <= k < len(grupos)), key=lambda k: cuenta(grupos[k])[1])
    a, b = sorted((i, j))
    grupos[a] = [grupos[a][0], grupos[b][1]]
    del grupos[b]
    fusiones += 1
bin2grupo = np.empty(NB + 1, dtype=int)
for k, g in enumerate(grupos):
    bin2grupo[g[0]:g[1] + 1] = k
g_nac, g_fondo = bin2grupo[bin_nac], bin2grupo[bin_fondo]
filas = []
for k, g in enumerate(grupos):
    n_n, n_f = cuenta(g)
    if n_n:
        filas.append({"intervalo": len(filas) + 1, "GC_desde": g[0] / NB, "GC_hasta": (g[1] + 1) / NB,
                      "n_NAC": int(n_n), "n_fondo": int(n_f), "razon_fondo_NAC": round(n_f / n_n, 1)})
iv = pd.DataFrame(filas)
F_IV = os.path.join(SAL, "fondo_GC_intervalos.tsv")
if SUF:   # los intervalos no dependen de N: se comprueba que coinciden con los de la colección original
    iv0 = pd.read_csv(F_IV, sep="\t")
    assert iv0.shape == iv.shape and np.allclose(iv0.to_numpy(float), iv.to_numpy(float)), "intervalos distintos del original"
else:
    iv.to_csv(F_IV, sep="\t", index=False)
print(f"intervalos con NAC: {len(iv)} (fusiones {fusiones}); razón fondo/NAC mínima {iv.razon_fondo_NAC.min()}")

# Muestreo: N_CONJ conjuntos, un promotor de fondo por NAC, del mismo intervalo, sin reemplazo
rng = np.random.default_rng(SEMILLA)
idx = np.empty((N_CONJ, 189), dtype=np.int32)
reservas = {k: np.flatnonzero(g_fondo == k) for k in np.unique(g_nac)}
posiciones = {k: np.flatnonzero(g_nac == k) for k in np.unique(g_nac)}
for s in range(N_CONJ):
    for k, P in reservas.items():
        J = posiciones[k]
        idx[s, J] = rng.choice(P, size=len(J), replace=False)

# Controles 3-5
assert (g_fondo[idx] == g_nac[None, :]).all(), "un promotor elegido no cae en el intervalo de su NAC"
assert all(len(np.unique(r)) == 189 for r in idx), "promotor repetido dentro de un conjunto"
gc_sel = gc_fondo[idx]
ks_emp = ks_2samp(gc_nac, gc_sel.ravel())
ks_bruto = ks_2samp(gc_nac, gc_fondo)
dmed = gc_sel.mean(axis=1) - gc_nac.mean()
usos = np.bincount(idx.ravel(), minlength=len(fondo))
print(f"control 3: {N_CONJ} conjuntos x 189, cada promotor en el intervalo de su NAC; sin repetidos")
print(f"control 4: GC medio NAC {gc_nac.mean():.4f}; fondo completo {gc_fondo.mean():.4f}; "
      f"conjuntos {gc_sel.mean():.4f} (dif. por conjunto {dmed.min():+.4f} a {dmed.max():+.4f})")
print(f"          KS NAC vs fondo elegido D={ks_emp.statistic:.4f} p={ks_emp.pvalue:.3g}; "
      f"NAC vs fondo completo D={ks_bruto.statistic:.4f} p={ks_bruto.pvalue:.3g}")
print(f"control 5: promotores de fondo usados {np.count_nonzero(usos)}/{len(fondo)}; "
      f"usos por promotor: mediana {int(np.median(usos[usos > 0]))}, máx. {usos.max()}")

igual_original = None
if SUF:
    z0 = np.load(os.path.join(SAL, "fondo_GC_conjuntos.npz"))
    n0 = z0["idx"].shape[0]
    igual_original = bool(N_CONJ >= n0 and (idx[:n0] == z0["idx"]).all() and list(z0["nac"]) == list(nac)
                          and list(z0["locus_fondo"]) == fondo.locus.astype(str).tolist())
    print(f"control 6: los primeros {n0} conjuntos son idénticos a los de la colección original: {igual_original}")
np.savez_compressed(os.path.join(SAL, f"fondo_GC_conjuntos{SUF}.npz"), nac=np.array(nac), locus_fondo=fondo.locus.to_numpy().astype(str),
                    idx=idx, gc_nac=gc_nac, gc_fondo=gc_fondo)
resumen = {
    "fecha": time.strftime("%Y-%m-%d %H:%M"), "semilla": SEMILLA, "n_conjuntos": N_CONJ, "ancho_GC": ANCHO_GC,
    "min_razon_fondo_NAC": MIN_RAZON, "n_NAC": 189, "n_fondo": int(len(fondo)), "intervalos_con_NAC": int(len(iv)),
    "fusiones": fusiones, "razon_minima": float(iv.razon_fondo_NAC.min()),
    "GC_medio": {"NAC": round(float(gc_nac.mean()), 5), "fondo_completo": round(float(gc_fondo.mean()), 5),
                 "conjuntos": round(float(gc_sel.mean()), 5)},
    "dif_GC_medio_por_conjunto": [round(float(dmed.min()), 5), round(float(dmed.max()), 5)],
    "KS_NAC_vs_fondo_elegido": {"D": round(float(ks_emp.statistic), 5), "p": float(ks_emp.pvalue)},
    "KS_NAC_vs_fondo_completo": {"D": round(float(ks_bruto.statistic), 5), "p": float(ks_bruto.pvalue)},
    "promotores_fondo_usados": int(np.count_nonzero(usos)), "usos_max": int(usos.max()),
    "primeros_conjuntos_iguales_a_la_coleccion_original": igual_original,
    "controles": {"GC_tabla_igual_FASTA": True, "NAC_identicos_al_genomico": 189, "intervalo_correcto": True,
                  "sin_repetidos_en_conjunto": True},
    "versiones": {"python": platform.python_version(), "numpy": np.__version__, "pandas": pd.__version__, "scipy": scipy.__version__},
}
with open(os.path.join(SAL, f"resumen_fondo_GC{SUF}.json"), "w", encoding="utf-8") as fh:
    json.dump(resumen, fh, ensure_ascii=False, indent=1)
print(f"salidas en {os.path.relpath(SAL, BASE)} ({time.time() - t0:.0f} s)")
