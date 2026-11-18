#!/usr/bin/env python3
# 20_subfamilias_Hussey.py -> reclasifica las subfamilias NAC del censo con las etiquetas
# PUBLICADAS de Hussey et al. (2015), no con las del manuscrito anterior.
# Uso, desde la raiz de "Revision NAC":
#     python scripts/20_subfamilias_Hussey.py              # solo calcula y compara
#     python scripts/20_subfamilias_Hussey.py --escribir   # respalda y reescribe las tablas
#
# Por que existe: hasta el 12-oct-2026 las anclas de E. grandis del arbol NAC llevaban las
# etiquetas del manuscrito previo proyectadas sobre 85 co-ortologos. Eso hacia circular la
# concordancia con el manuscrito y dejaba fuera VIIIa (no habia anclas VIII), de modo que
# EglNAC36 y EglNAC204 caian en IIIc.
#
# 1. Fig. S1 de Hussey (suppinfohussey/nph13139-sup-0001-figs1.pdf) es un PDF vectorial: los
#    nombres de las puntas son texto (x ~ 500) y cada subfamilia es un corchete (objeto path
#    de 5,6 pt de ancho, x ~ 541-547) con su rotulo centrado (x ~ 549-556). Cada punta EgrNAC
#    se asigna al corchete que contiene su coordenada y. EgrNAC151 queda fuera de todo
#    corchete (Unassigned). Faltan en la figura EgrNAC29 (sin etiqueta publicada) y EgrNAC91,
#    que toma la de la Tabla S10 (IVa, bloque en tandem 13). La Tabla S10 contradice a la
#    figura en EgrNAC137-139 (S10: IVc; Fig. S1: VIc, clado de EgrNAC152 y ANAC104); se usa
#    la figura, que es la evidencia filogenetica.  -> Hussey_FigS1_subfamilias_189.csv
# 2. Las puntas GR_EUGRSUZ_<id> del arbol se traducen a EgrNAC con la Tabla S4
#    (EUGRSUZ_ -> Eucgr.). Cada gen de E. globulus (puntas GL_<locus>) hereda la subfamilia
#    del ancla etiquetada mas proxima por distancia patristica. margen = distancia al ancla
#    mas proxima de OTRA subfamilia menos la distancia al ancla mas proxima;
#    fiable = dist_ancla <= 1,0 y margen >= 0,10 (misma regla que antes; con las etiquetas
#    antiguas reproduce T01 183/183 en subfamilia, ancla, distancia y fiabilidad).
# 3. Con --escribir: copia las tablas a ANALISIS_Eglobulus/_respaldo_pre_subfamilias_Hussey/
#    (solo la primera vez) y reescribe, campo a campo y sin tocar las demas columnas,
#    subfamilia_final/ancla/dist_ancla/fiable en T01 y la copia de subfamilia_final en T02,
#    T06, T13, T16, T14 (subfam_1, subfam_2, misma_subfamilia) y T11 (subfamilia).
#    T36 lo regenera 17_de_globulus.R; las figuras, scripts/correr_todo.R.
import csv, os, re, shutil, sys
import numpy as np
import pandas as pd
import pypdfium2 as pdfium
import pypdfium2.raw as pdfium_c
from Bio import Phylo

RAIZ = os.getcwd()
TABLAS = os.path.join(RAIZ, "ANALISIS_Eglobulus", "tablas")
RESPALDO = os.path.join(RAIZ, "ANALISIS_Eglobulus", "_respaldo_pre_subfamilias_Hussey")
FIG_S1 = os.path.join(RAIZ, "suppinfohussey", "nph13139-sup-0001-figs1.pdf")
S4 = os.path.join(RAIZ, "Hussey_TablaS4_189simbolos.csv")
ARBOL = os.path.join(RAIZ, "arbol_NAC_189.treefile")
SALIDA_HUSSEY = os.path.join(RAIZ, "Hussey_FigS1_subfamilias_189.csv")
NO_EVAL = "no evaluada (dominio parcial)"
ROTULOS = r"(I[abc]|II|IIIa/b|IIIc|IV[abcd]|Va\([12]\)|Vb|VI[abc]|VII|VIII[ab]|IX|X|XI)"
FUERA_DE_CORCHETE = {"EgrNAC151": "Unassigned"}   # en la figura, sin corchete
DE_S10 = {"EgrNAC91": "IVa"}                      # ausente de la figura; Tabla S10, bloque 13
CONFLICTO_S10 = {"EgrNAC137": "IVc", "EgrNAC138": "IVc", "EgrNAC139": "IVc"}  # S10, bloque 18
DIST_MAX, MARGEN_MIN = 1.0, 0.10
TABLAS_COPIA = {"T02_equivalencia_nomenclaturas_189.csv": "subfamilia_final",
                "T06_subdominios_189.csv": "subfamilia_final",
                "T13_vecindades_NAC.csv": "subfamilia_final",
                "T16_arquitectura_dominios.csv": "subfamilia_final"}


def etiquetas_fig_s1(ruta=FIG_S1):
    """{EgrNAC: subfamilia} leido de los corchetes de la Fig. S1."""
    pdf = pdfium.PdfDocument(ruta); pag = pdf[0]; tp = pag.get_textpage()
    runs = []
    for i in range(tp.count_rects()):
        l, b, r, t = tp.get_rect(i)
        runs.append((tp.get_text_bounded(l, b, r, t).strip(), l, (b + t) / 2))
    corchetes = [o.get_bounds() for o in pag.get_objects(max_depth=4)
                 if o.type == pdfium_c.FPDF_PAGEOBJ_PATH and o.get_bounds()[0] > 505]
    rotulos = [(txt, y) for txt, x, y in runs if x > 540 and re.fullmatch(ROTULOS, txt)]
    if len(corchetes) != len(rotulos):
        raise ValueError("Fig. S1: %d corchetes y %d rotulos" % (len(corchetes), len(rotulos)))
    nombre = {}
    for txt, y in rotulos:
        c = [k for k, (x0, y0, x1, y1) in enumerate(corchetes) if y0 <= y <= y1]
        if len(c) != 1:
            raise ValueError("rotulo %s en y=%.1f cae en %d corchetes" % (txt, y, len(c)))
        nombre[c[0]] = txt
    out, fuera = {}, []
    for txt, x, y in runs:
        m = re.search(r"EgrNAC\d+", txt)
        if not m or not 480 < x < 520:
            continue
        c = [k for k, (x0, y0, x1, y1) in enumerate(corchetes) if y0 <= y <= y1]
        if len(c) == 1:
            out[m.group(0)] = nombre[c[0]]
        else:
            fuera.append(m.group(0))
    if sorted(fuera) != sorted(FUERA_DE_CORCHETE):
        raise ValueError("puntas fuera de corchete inesperadas: %s" % fuera)
    out.update(FUERA_DE_CORCHETE)
    return out


def tabla_hussey():
    fig = etiquetas_fig_s1()
    s4 = pd.read_csv(S4)
    h = s4[["simbolo_EgrNAC", "locus_v1_0"]].copy()
    h["subfamilia_FigS1"] = h.simbolo_EgrNAC.map(fig)
    h["nota_S10"] = h.simbolo_EgrNAC.map(
        {**{g: "S10 %s; ausente de la Fig. S1" % s for g, s in DE_S10.items()},
         **{g: "S10 %s; contradice la Fig. S1, se usa la figura" % s for g, s in CONFLICTO_S10.items()}})
    h["subfamilia_Hussey"] = h.subfamilia_FigS1.fillna(h.simbolo_EgrNAC.map(DE_S10))
    h["fuente"] = np.where(h.subfamilia_FigS1.notna(), "Fig. S1",
                           np.where(h.subfamilia_Hussey.notna(), "Tabla S10", "sin etiqueta"))
    return h


class Patristica:
    def __init__(self, arbol):
        self.dep = arbol.depths()
        self.ruta = {t.name: [arbol.root] + arbol.get_path(t) for t in arbol.get_terminals()}

    def __call__(self, a, b):
        pa, pb = self.ruta[a], self.ruta[b]; k = 0
        while k < min(len(pa), len(pb)) and pa[k] is pb[k]:
            k += 1
        return self.dep[pa[-1]] + self.dep[pb[-1]] - 2 * self.dep[pa[k - 1]]


def asignar(h):
    arbol = Phylo.read(ARBOL, "newick"); d = Patristica(arbol)
    puntas = [t.name for t in arbol.get_terminals()]
    egr = dict(zip(h.locus_v1_0, h.simbolo_EgrNAC)); lab = dict(zip(h.simbolo_EgrNAC, h.subfamilia_Hussey))
    anclas = {}
    for p in puntas:
        if p.startswith("GR_"):
            s = egr.get(p[3:].replace("EUGRSUZ_", "Eucgr."))
            if s and isinstance(lab.get(s), str):
                anclas[p] = (lab[s], s)
    filas = []
    for g in (p for p in puntas if p.startswith("GL_")):
        ds = sorted((d(g, a), s, sym) for a, (s, sym) in anclas.items())
        d1, s1, sym1 = ds[0]
        d2 = next((x for x, s, _ in ds if s != s1), np.inf)
        m = round(d2 - d1, 4) if np.isfinite(d2) else np.inf
        filas.append(dict(locus=g[3:], subfamilia_final=s1, ancla=sym1, dist_ancla=round(d1, 4),
                          margen=m, fiable=bool(round(d1, 4) <= DIST_MAX and m >= MARGEN_MIN)))
    n_gr = sum(p.startswith("GR_") for p in puntas)
    return pd.DataFrame(filas), len(anclas), n_gr


def _fmt(v):
    if isinstance(v, (bool, np.bool_)):
        return "True" if v else "False"
    if isinstance(v, float):
        return repr(round(v, 4)) if np.isfinite(v) else ""
    return "" if v is None else str(v)


def reescribir(ruta, clave, valores, cambios_extra=None):
    """Sustituye in situ las columnas de `valores` ({col: {clave: valor}}) sin tocar el resto."""
    with open(ruta, encoding="utf-8", newline="") as f:
        filas = list(csv.reader(f))
    cab = filas[0]; jk = cab.index(clave); n = 0
    for fila in filas[1:]:
        for col, mapa in valores.items():
            if fila[jk] in mapa:
                nuevo = _fmt(mapa[fila[jk]]); j = cab.index(col)
                if fila[j] != nuevo:
                    fila[j] = nuevo; n += 1
        if cambios_extra:
            n += cambios_extra(cab, fila)
    with open(ruta, "w", encoding="utf-8", newline="") as f:
        csv.writer(f, lineterminator="\n").writerows(filas)
    return n


def main(escribir):
    h = tabla_hussey()
    nueva, n_anclas, n_gr = asignar(h)
    t01 = pd.read_csv(os.path.join(TABLAS, "T01_censo_maestro_189.csv"))
    ev = t01[t01.subfamilia_final != NO_EVAL]
    if set(ev.locus) != set(nueva.locus):
        raise ValueError("las puntas GL del arbol no son los genes evaluables de T01")
    cmp = nueva.merge(t01[["locus", "simbolo", "subfamilia_final", "fiable"]], on="locus", suffixes=("", "_antes"))
    camb = cmp[cmp.subfamilia_final != cmp.subfamilia_final_antes]
    print("anclas etiquetadas: %d de %d puntas de E. grandis" % (n_anclas, n_gr))
    print("genes evaluados: %d | fiables: %d (antes %d)" % (len(cmp), cmp.fiable.sum(), cmp.fiable_antes.sum()))
    print("cambios de subfamilia: %d" % len(camb))
    print(camb.groupby(["subfamilia_final_antes", "subfamilia_final"]).simbolo.apply(lambda s: ", ".join(s)).to_string())
    print("subfamilias con genes: %d" % nueva.subfamilia_final.nunique())
    print(nueva.subfamilia_final.value_counts().to_string())
    if not escribir:
        return
    os.makedirs(RESPALDO, exist_ok=True)
    afect = ["T01_censo_maestro_189.csv", "T11_especificidad_primers_qpcr.csv",
             "T14_pares_NAC_adyacentes.csv"] + list(TABLAS_COPIA)
    for f in afect:
        dst = os.path.join(RESPALDO, f)
        if not os.path.exists(dst):
            shutil.copy2(os.path.join(TABLAS, f), dst)
    h.to_csv(SALIDA_HUSSEY, index=False)
    nv = nueva.set_index("locus")
    sub = nv.subfamilia_final.to_dict()
    n = reescribir(os.path.join(TABLAS, "T01_censo_maestro_189.csv"), "locus",
                   {"subfamilia_final": sub, "ancla": nv.ancla.to_dict(),
                    "dist_ancla": nv.dist_ancla.to_dict(), "fiable": nv.fiable.to_dict()})
    print("T01: %d campos" % n)
    for f, col in TABLAS_COPIA.items():
        print("%s: %d campos" % (f[:3], reescribir(os.path.join(TABLAS, f), "locus", {col: sub})))
    t01n = pd.read_csv(os.path.join(TABLAS, "T01_censo_maestro_189.csv"))
    por_simbolo = t01n.set_index("simbolo").subfamilia_final.to_dict()
    print("T11: %d campos" % reescribir(os.path.join(TABLAS, "T11_especificidad_primers_qpcr.csv"),
                                        "simbolo_censo_189", {"subfamilia": por_simbolo}))
    def misma(cab, fila):
        j = cab.index("misma_subfamilia")
        v = "True" if fila[cab.index("subfam_1")] == fila[cab.index("subfam_2")] else "False"
        if fila[j] != v:
            fila[j] = v; return 1
        return 0
    por_locus = t01n.set_index("locus").subfamilia_final.to_dict()
    with open(os.path.join(TABLAS, "T14_pares_NAC_adyacentes.csv"), encoding="utf-8") as f:
        cab14 = next(csv.reader(f))
    n14 = reescribir(os.path.join(TABLAS, "T14_pares_NAC_adyacentes.csv"), "locus_1", {"subfam_1": por_locus})
    n14 += reescribir(os.path.join(TABLAS, "T14_pares_NAC_adyacentes.csv"), "locus_2", {"subfam_2": por_locus}, misma)
    print("T14: %d campos" % n14)


if __name__ == "__main__":
    main("--escribir" in sys.argv[1:])
