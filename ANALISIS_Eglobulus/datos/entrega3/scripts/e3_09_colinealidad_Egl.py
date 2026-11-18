"""Entrega 3, paso 4.1: colinealidad y modos de duplicación dentro de E. globulus.

Entradas
  datos/entrega3/isoformas/Egl.bed, Egl.pep   un gen por locus en los 11 cromosomas (paso 1.3):
                                              36 278 genes; en los 189 NAC las coordenadas son las de T01
  tablas/T01_censo_maestro_189.csv
  tablas/T15_bloques_tandem.csv

Proceso
  1. DIAMOND blastp de todos contra todos (--ultra-sensitive, e <= 1e-10, 5 dianas por consulta
     contando la propia): la receta BLASTP del manual de MCScanX (-e 1e-10 -b 5 -v 5).
  2. MCScanX con sus parámetros por defecto (-k 50 -g -1 -s 5 -e 1e-5 -m 25 -w 5): bloques
     colineales (.collinearity) y pares en tándem (.tandem).
  3. duplicate_gene_classifier con los mismos parámetros y -n 10. Según el código de MCScanX
     (cls.cc), cada gen queda como singleton (0, sin acierto), disperso (1), proximal (2, homólogo
     a menos de 10 genes en el mismo cromosoma), tándem (3, homólogo adyacente) o WGD/segmental
     (4, ancla de un bloque colineal), y manda la clase más alta. Los aciertos consigo mismo no cuentan.
  4. Regla del artículo 3 (Zhang et al. 2026, BMC Genomics): pares de NAC en el mismo cromosoma a
     <= 200 kb, medidos entre los bordes de los genes; los grupos salen por enlace simple. Contraste con
     T01.bloque_tandem (<= 5 genes intercalados y <= 200 kb entre NAC consecutivos) y con T15.
  5. Proporción de cada clase en los 189 NAC frente al resto del genoma (Fisher bilateral, BH).

Salidas
  datos/entrega3/colinealidad/mcscanx/Egl.{gff,blast,collinearity,tandem,gene_type}, Egl.html/
  datos/entrega3/colinealidad/Egl_pares_colineales.tsv   todos los pares ancla del genoma (para Ks)
  tablas/T45_modos_duplicacion_NAC189.csv                 un NAC por fila
  tablas/T46_pares_duplicados_NAC.csv                     pares NAC-NAC segmentales, proximales, en tándem o a <= 200 kb
  datos/entrega3/colinealidad/resumen_colinealidad_Egl.json

Uso, desde la raíz de Revisión NAC:
  <nac3>/bin/python -u ANALISIS_Eglobulus/datos/entrega3/scripts/e3_09_colinealidad_Egl.py [--forzar]
  --forzar repite DIAMOND aunque ya exista Egl.blast.
"""
import argparse, json, os, re, subprocess, sys, time
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

RAIZ = (os.environ.get("NAC_ROOT", ".") + "")
A = f"{RAIZ}/ANALISIS_Eglobulus"
E3 = f"{A}/datos/entrega3"
ISO = f"{E3}/isoformas"
COL = f"{E3}/colinealidad"
MC = f"{COL}/mcscanx"
PRE = f"{MC}/Egl"
BIN = (os.environ.get("CONDA_ENVS", os.path.expanduser("~/conda/envs")) + "/nac3/bin")
T01 = f"{A}/tablas/T01_censo_maestro_189.csv"
T15 = f"{A}/tablas/T15_bloques_tandem.csv"
T45 = f"{A}/tablas/T45_modos_duplicacion_NAC189.csv"
T46 = f"{A}/tablas/T46_pares_duplicados_NAC.csv"

EVALUE, DIANAS, HILOS = 1e-10, 5, 12
KB_ART3 = 200.0
N_PROXIMAL = 10
CLASES = {0: "singleton", 1: "dispersed", 2: "proximal", 3: "tandem", 4: "WGD/segmental"}
COLS = "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore".split()
N_GENES, N_NAC = 36278, 189


def corto(c):
    return c.replace(BIN + "/", "").replace(RAIZ, ".")


def correr(cmd, **kw):
    print("$", " ".join(corto(c) for c in cmd), flush=True)
    t = time.time()
    r = subprocess.run(cmd, check=True, **kw)
    print(f"  ({time.time() - t:.0f} s)", flush=True)
    return r


def genes():
    """Egl.bed -> genes en el orden de MCScanX (struct.h: cromosoma, inicio, nombre)."""
    g = pd.read_csv(f"{ISO}/Egl.bed", sep="\t", header=None, names=["crom", "ini0", "fin", "gen", "puntaje", "hebra"])
    g["inicio"] = g.ini0 + 1
    g = g.sort_values(["crom", "inicio", "gen"], kind="mergesort").reset_index(drop=True)
    assert len(g) == N_GENES and g.gen.is_unique and g.crom.nunique() == 11, g.shape
    g["rango"] = g.groupby("crom").cumcount()
    return g


def diamond(forzar):
    out = f"{PRE}.blast"
    if os.path.exists(out) and not forzar:
        print(f"DIAMOND: se reutiliza {corto(out)} (--forzar para repetir)", flush=True)
        return out
    correr([f"{BIN}/diamond", "makedb", "--in", f"{ISO}/Egl.pep", "-d", f"{MC}/Egl_db", "--threads", str(HILOS)])
    correr([f"{BIN}/diamond", "blastp", "-q", f"{ISO}/Egl.pep", "-d", f"{MC}/Egl_db", "-o", out + ".tmp",
            "--ultra-sensitive", "--evalue", str(EVALUE), "--max-target-seqs", str(DIANAS),
            "--threads", str(HILOS), "--outfmt", "6", *COLS])
    os.replace(out + ".tmp", out)
    return out


def leer_colinealidad(p):
    est, bloques, anclas = {}, [], []
    rb = re.compile(r"## Alignment (\d+): score=(\S+) e_value=(\S+) N=(\d+) (\S+)&(\S+) (plus|minus)")
    ra = re.compile(r"\s*(\d+)-\s*(\d+):\s+(\S+)\s+(\S+)\s+(\S+)")
    for l in open(p):
        if l.startswith("# Number of collinear genes"):
            m = re.search(r"genes: (\d+), Percentage: ([\d.]+)", l)
            est["genes_colineales"], est["porcentaje"] = int(m[1]), float(m[2])
        elif l.startswith("# Number of all genes"):
            est["genes_totales"] = int(re.search(r"(\d+)", l)[1])
        elif l.startswith("## Alignment"):
            m = rb.match(l)
            assert m, l
            bloques.append(dict(bloque=int(m[1]), puntaje=float(m[2]), evalue_bloque=float(m[3]), n_anclas=int(m[4]),
                                crom1=m[5], crom2=m[6], orientacion=m[7]))
        elif l.strip() and not l.startswith("#"):
            m = ra.match(l)
            assert m, l
            anclas.append(dict(bloque=int(m[1]), gen1=m[3], gen2=m[4], evalue_par=float(m[5])))
    b, an = pd.DataFrame(bloques), pd.DataFrame(anclas)
    assert (an.groupby("bloque").size().reindex(b.bloque).values == b.n_anclas.values).all()
    return est, b, an


def grupos_enlace_simple(pares, orden):
    """pares [(x, y)] -> {gen: id de grupo}, grupos numerados por posición en el genoma."""
    padre = {}
    def raiz(x):
        while padre.setdefault(x, x) != x:
            padre[x] = padre[padre[x]]
            x = padre[x]
        return x
    for x, y in pares:
        padre[raiz(x)] = raiz(y)
    comp = {}
    for x in padre:
        comp.setdefault(raiz(x), []).append(x)
    comp = sorted(comp.values(), key=lambda m: min(orden[v] for v in m))
    return {v: f"A{i + 1:02d}" for i, m in enumerate(comp) for v in m}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--forzar", action="store_true")
    arg = ap.parse_args()
    os.makedirs(MC, exist_ok=True)
    print(time.strftime("%Y-%m-%d %H:%M"), "| python", sys.version.split()[0], "| pandas", pd.__version__, flush=True)
    print(subprocess.run([f"{BIN}/diamond", "version"], capture_output=True, text=True).stdout.strip(), flush=True)

    # --- 1. genes y DIAMOND ------------------------------------------------------------------------
    g = genes()
    g[["crom", "gen", "inicio", "fin"]].to_csv(f"{PRE}.gff", sep="\t", header=False, index=False)
    gi = g.set_index("gen")
    t01 = pd.read_csv(T01)
    nac = t01[["simbolo", "locus", "cromosoma", "inicio", "fin", "hebra", "subfamilia_final", "bloque_tandem"]]
    nac = nac.rename(columns={"subfamilia_final": "subfamilia", "bloque_tandem": "bloque_tandem_T01"})
    assert len(nac) == N_NAC and nac.locus.isin(gi.index).all()
    assert (gi.loc[nac.locus, "inicio"].values == nac.inicio.values).all()
    assert (gi.loc[nac.locus, "fin"].values == nac.fin.values).all()
    nac = nac.assign(crom=gi.loc[nac.locus, "crom"].values, rango=gi.loc[nac.locus, "rango"].values)
    nac = nac.sort_values(["crom", "inicio", "locus"], kind="mergesort").reset_index(drop=True)
    nac["orden_NAC"] = nac.groupby("crom").cumcount()
    NACS = set(nac.locus)
    N = nac.set_index("locus")
    orden = {v: i for i, v in enumerate(g.gen)}

    bl = pd.read_csv(diamond(arg.forzar), sep="\t", names=COLS)
    propio = bl.qseqid == bl.sseqid
    sin_propio = len(set(g.gen) - set(bl.qseqid[propio]))
    otros = bl[~propio].copy()
    otros["a"] = np.where(otros.qseqid < otros.sseqid, otros.qseqid, otros.sseqid)
    otros["b"] = np.where(otros.qseqid < otros.sseqid, otros.sseqid, otros.qseqid)
    otros = otros.sort_values(["bitscore", "evalue"], ascending=[False, True]).drop_duplicates(["a", "b"])
    hn = otros[otros.a.isin(NACS) & otros.b.isin(NACS)]
    hits_nac = {(r.a, r.b): r for r in hn.itertuples()}
    print(f"DIAMOND: {len(bl)} aciertos | {bl.qseqid.nunique()} consultas | {len(otros)} pares no propios"
          f" | {sin_propio} genes sin acierto propio | {len(hits_nac)} pares NAC-NAC", flush=True)

    # --- 2-3. MCScanX y duplicate_gene_classifier --------------------------------------------------
    correr([f"{BIN}/MCScanX", PRE], stdout=subprocess.DEVNULL)
    r = correr([f"{BIN}/duplicate_gene_classifier", PRE], capture_output=True, text=True)
    tabla_dgc = [l for l in r.stdout.splitlines() if re.match(r"^(Type of dup|Singleton|Dispersed|Proximal|Tandem|WGD)", l)]
    print("\n".join("  " + l for l in tabla_dgc), flush=True)
    est, bloques, an = leer_colinealidad(f"{PRE}.collinearity")
    # msa.cc (mark_tandem): par con acierto y genes contiguos en el mismo cromosoma; sin pares no hay archivo
    tan = (pd.read_csv(f"{PRE}.tandem", header=None, names=["gen1", "gen2"]) if os.path.exists(f"{PRE}.tandem")
           else pd.DataFrame(columns=["gen1", "gen2"]))
    gt = pd.read_csv(f"{PRE}.gene_type", sep="\t", header=None, names=["gen", "clase"])
    assert set(gt.gen) == set(g.gen) and gt.gen.is_unique
    clase = dict(zip(gt.gen, gt.clase))
    anclas_genes = set(an.gen1) | set(an.gen2)
    clase4 = {v for v, c in clase.items() if c == 4}
    print(f"MCScanX: {len(bloques)} bloques | {len(an)} pares ancla | {len(anclas_genes)} genes ancla"
          f" ({est['genes_colineales']} y {est['porcentaje']} % según la cabecera) | {len(tan)} pares en tándem"
          f" | clase 4 = anclas: {clase4 == anclas_genes}", flush=True)

    an = an.merge(bloques, on="bloque")
    an["NAC1"], an["NAC2"] = an.gen1.isin(NACS), an.gen2.isin(NACS)
    an[["bloque", "crom1", "crom2", "orientacion", "puntaje", "evalue_bloque", "n_anclas", "gen1", "gen2",
        "evalue_par", "NAC1", "NAC2"]].to_csv(f"{COL}/Egl_pares_colineales.tsv", sep="\t", index=False)

    # --- pares NAC-NAC ------------------------------------------------------------------------------
    def ordenar(x, y):
        return (x, y) if orden[x] <= orden[y] else (y, x)

    seg = {}
    for q in an[an.NAC1 & an.NAC2].itertuples():
        seg.setdefault(ordenar(q.gen1, q.gen2), []).append(q)
    tanN = {ordenar(x, y) for x, y in zip(tan.gen1, tan.gen2) if x in NACS and y in NACS}
    claves = set(seg) | tanN
    for _, d in nac.groupby("crom"):
        loc = d.locus.tolist()
        for i in range(len(loc)):
            for j in range(i + 1, len(loc)):
                x, y = loc[i], loc[j]
                gap = max(0, max(N.at[x, "inicio"], N.at[y, "inicio"]) - min(N.at[x, "fin"], N.at[y, "fin"])) / 1e3
                dr = abs(int(N.at[x, "rango"]) - int(N.at[y, "rango"]))
                if gap <= KB_ART3 or ((min(x, y), max(x, y)) in hits_nac and dr < N_PROXIMAL):
                    claves.add(ordenar(x, y))

    t15 = pd.read_csv(T15)
    def bloque15(c, i, f):
        m = t15[(t15.cromosoma == c) & (t15.inicio <= f) & (t15.fin >= i)]
        assert len(m) <= 1
        return m.bloque.iloc[0] if len(m) else ""
    nac["bloque_T15"] = [bloque15(r.cromosoma, r.inicio, r.fin) for r in nac.itertuples()]
    N = nac.set_index("locus")

    filas = []
    for x, y in sorted(claves, key=lambda k: (orden[k[0]], orden[k[1]])):
        a, b = N.loc[x], N.loc[y]
        mismo = a.crom == b.crom
        h = hits_nac.get((min(x, y), max(x, y)))
        dr = abs(int(a.rango) - int(b.rango)) if mismo else np.nan
        gap = max(0, max(a.inicio, b.inicio) - min(a.fin, b.fin)) / 1e3 if mismo else np.nan
        s = seg.get((x, y), [])
        if s:
            modo = "WGD/segmental"
        elif mismo and h is not None and dr == 1:
            modo = "tandem"
        elif mismo and h is not None and dr < N_PROXIMAL:
            modo = "proximal"
        else:
            modo = ""
        filas.append(dict(
            simbolo1=a.simbolo, simbolo2=b.simbolo, locus1=x, locus2=y, cromosoma1=a.cromosoma, cromosoma2=b.cromosoma,
            subfamilia1=a.subfamilia, subfamilia2=b.subfamilia, misma_subfamilia=a.subfamilia == b.subfamilia,
            modo_par_MCScanX=modo, tandem_MCScanX=(x, y) in tanN,
            bloque_MCScanX=";".join(str(q.bloque) for q in s), orientacion_bloque=";".join(q.orientacion for q in s),
            n_anclas_bloque=";".join(str(q.n_anclas) for q in s), evalue_bloque=";".join(f"{q.evalue_bloque:.2g}" for q in s),
            evalue_par_MCScanX=";".join(f"{q.evalue_par:.2g}" for q in s),
            mismo_cromosoma=mismo, distancia_kb=round(gap, 1) if mismo else np.nan,
            distancia_inicios_kb=round(abs(a.inicio - b.inicio) / 1e3, 1) if mismo else np.nan,
            genes_intercalados=dr - 1 if mismo else np.nan,
            NAC_intercalados=abs(int(a.orden_NAC) - int(b.orden_NAC)) - 1 if mismo else np.nan,
            tandem_art3=bool(mismo and gap <= KB_ART3),
            acierto_DIAMOND=h is not None, identidad_pct=h.pident if h is not None else np.nan,
            bitscore=h.bitscore if h is not None else np.nan, evalue_DIAMOND=h.evalue if h is not None else np.nan,
            mismo_bloque_T01=bool(pd.notna(a.bloque_tandem_T01) and a.bloque_tandem_T01 == b.bloque_tandem_T01),
            mismo_bloque_T15=bool(a.bloque_T15 != "" and a.bloque_T15 == b.bloque_T15)))
    t46 = pd.DataFrame(filas)
    t46 = t46.astype({"genes_intercalados": "Int64", "NAC_intercalados": "Int64"})
    p3 = t46[t46.tandem_art3]
    grupo = grupos_enlace_simple(zip(p3.locus1, p3.locus2), orden)
    t46.insert(t46.columns.get_loc("tandem_art3") + 1, "grupo_art3",
               [grupo.get(x, "") if t else "" for x, t in zip(t46.locus1, t46.tandem_art3)])
    t46.to_csv(T46, index=False)

    # --- T45: un NAC por fila --------------------------------------------------------------------------
    def socios(x, conj, nac_si):
        out = []
        for u, v in conj:
            if x in (u, v):
                o = v if u == x else u
                if (o in NACS) == nac_si:
                    out.append(N.at[o, "simbolo"] if nac_si else o)
        return ";".join(sorted(set(out)))
    seg_todos = {ordenar(q.gen1, q.gen2) for q in an[an.NAC1 | an.NAC2].itertuples()}
    tan_todos = {ordenar(u, v) for u, v in zip(tan.gen1, tan.gen2) if u in NACS or v in NACS}
    bloques_de = {}
    for q in an[an.NAC1 | an.NAC2].itertuples():
        for v in (q.gen1, q.gen2):
            if v in NACS:
                bloques_de.setdefault(v, set()).add(q.bloque)
    tam_grupo = pd.Series(grupo).value_counts()
    t45 = nac[["simbolo", "locus", "cromosoma", "inicio", "fin", "hebra", "subfamilia"]].copy()
    t45["rango_en_cromosoma"] = nac.rango + 1
    t45["clase_MCScanX"] = [clase[v] for v in nac.locus]
    t45["modo_duplicacion"] = t45.clase_MCScanX.map(CLASES)
    t45["bloques_MCScanX"] = [";".join(map(str, sorted(bloques_de.get(v, [])))) for v in nac.locus]
    t45["socios_segmentales"] = [socios(v, seg_todos, True) for v in nac.locus]
    t45["socios_segmentales_no_NAC"] = [socios(v, seg_todos, False) for v in nac.locus]
    t45["socios_tandem"] = [socios(v, tan_todos, True) for v in nac.locus]
    t45["socios_tandem_no_NAC"] = [socios(v, tan_todos, False) for v in nac.locus]
    t45["grupo_art3"] = [grupo.get(v, "") for v in nac.locus]
    t45["n_genes_grupo_art3"] = [int(tam_grupo.get(grupo[v], 0)) if v in grupo else 0 for v in nac.locus]
    t45["bloque_tandem_T01"] = nac.bloque_tandem_T01.fillna("")
    t45["bloque_T15"] = nac.bloque_T15
    t45.to_csv(T45, index=False)

    # --- 5. NAC frente al genoma ----------------------------------------------------------------------
    cn = t45.clase_MCScanX
    resto = gt[~gt.gen.isin(NACS)].clase
    fis = []
    for c, nombre in CLASES.items():
        a1, c1 = int((cn == c).sum()), int((resto == c).sum())
        OR, p = fisher_exact([[a1, N_NAC - a1], [c1, len(resto) - c1]])
        fis.append(dict(clase=nombre, NAC=a1, pct_NAC=round(100 * a1 / N_NAC, 1), resto=c1,
                        pct_resto=round(100 * c1 / len(resto), 1), OR=round(float(OR), 2) if np.isfinite(OR) else None,
                        p=float(p)))
    q = multipletests([f["p"] for f in fis], method="fdr_bh")[1]
    for f, qq in zip(fis, q):
        f["q_BH"] = qq

    # --- contrastes de tándem -----------------------------------------------------------------------------
    def conjuntos(col):
        s = t45[t45[col] != ""]
        return {k: frozenset(v) for k, v in s.groupby(col).locus}
    G3, G01, G15 = conjuntos("grupo_art3"), conjuntos("bloque_tandem_T01"), conjuntos("bloque_T15")
    gen3, gen01, gen15 = set().union(*G3.values()), set().union(*G01.values()), set().union(*G15.values())
    s15_listados = {s.strip() for x in t15.simbolos for s in x.split(",")}
    ady = (t46.mismo_cromosoma & t46.acierto_DIAMOND & (t46.genes_intercalados == 0)).fillna(False).astype(bool)
    ady_acierto = set(zip(t46.locus1[ady], t46.locus2[ady]))
    contraste = dict(
        art3=dict(pares=int(t46.tandem_art3.sum()), pares_NAC_consecutivos=int((t46.tandem_art3 & (t46.NAC_intercalados == 0)).sum()),
                  pares_con_acierto_DIAMOND=int((t46.tandem_art3 & t46.acierto_DIAMOND).sum()),
                  pares_inicios_200kb=int((t46.distancia_inicios_kb <= KB_ART3).sum()),
                  grupos=len(G3), genes=len(gen3),
                  tamanos={int(k): int(v) for k, v in pd.Series([len(v) for v in G3.values()]).value_counts().sort_index().items()}),
        T01=dict(bloques=len(G01), genes=len(gen01), genes_en_art3=len(gen01 & gen3), art3_fuera_de_T01=len(gen3 - gen01),
                 T01_fuera_de_art3=len(gen01 - gen3), bloques_identicos_a_un_grupo_art3=len(set(G01.values()) & set(G3.values()))),
        T15=dict(bloques=len(G15), genes=len(gen15), coincide_con_simbolos_listados=set(N.loc[list(gen15), "simbolo"]) == s15_listados,
                 genes_en_art3=len(gen15 & gen3), genes_en_T01=len(gen15 & gen01),
                 bloques_identicos_a_un_grupo_art3=len(set(G15.values()) & set(G3.values()))),
        MCScanX=dict(genes_clase_tandem=int((cn == 3).sum()), genes_clase_proximal=int((cn == 2).sum()),
                     pares_tandem=len(tanN), pares_tandem_igual_a_adyacentes_con_acierto=ady_acierto == tanN,
                     clase_de_genes_en_art3=t45[t45.grupo_art3 != ""].modo_duplicacion.value_counts().to_dict(),
                     clase_de_genes_en_T01=t45[t45.bloque_tandem_T01 != ""].modo_duplicacion.value_counts().to_dict()))

    segp = t46[t46.modo_par_MCScanX == "WGD/segmental"]
    res = dict(
        parametros=dict(DIAMOND=f"blastp --ultra-sensitive --evalue {EVALUE:g} --max-target-seqs {DIANAS}",
                        MCScanX="por defecto: -k 50 -g -1 -s 5 -e 1e-5 -m 25 -w 5", duplicate_gene_classifier=f"por defecto, -n {N_PROXIMAL}",
                        regla_art3=f"pares NAC en el mismo cromosoma a <= {KB_ART3:g} kb entre bordes de gen; grupos por enlace simple"),
        genoma=dict(genes=N_GENES, cromosomas=int(g.crom.nunique()), aciertos_DIAMOND=len(bl), pares_no_propios=len(otros),
                    genes_sin_acierto_propio=sin_propio,
                    MCScanX=dict(bloques=len(bloques), pares_ancla=len(an), genes_ancla=len(anclas_genes), **est,
                                 bloques_intracromosomicos=int((bloques.crom1 == bloques.crom2).sum()),
                                 bloques_invertidos=int((bloques.orientacion == "minus").sum()),
                                 anclas_por_bloque_mediana=float(bloques.n_anclas.median()), pares_tandem=len(tan),
                                 clase4_igual_a_anclas=clase4 == anclas_genes),
                    clases={CLASES[c]: int((gt.clase == c).sum()) for c in CLASES}, salida_duplicate_gene_classifier=tabla_dgc),
        NAC189=dict(clases={CLASES[c]: int((cn == c).sum()) for c in CLASES}, fisher_NAC_vs_resto=fis,
                    pares=dict(T46=len(t46), segmentales=len(segp), bloques_con_pares_NAC=int(segp.bloque_MCScanX.str.split(";").explode().nunique()) if len(segp) else 0,
                               segmentales_intercromosomicos=int((~segp.mismo_cromosoma).sum()),
                               segmentales_misma_subfamilia=int(segp.misma_subfamilia.sum()),
                               tandem_MCScanX=len(tanN), proximales=int((t46.modo_par_MCScanX == "proximal").sum())),
                    NAC_ancla_con_socio_no_NAC=int((t45.socios_segmentales_no_NAC != "").sum()),
                    NAC_tandem_con_socio_no_NAC=int((t45.socios_tandem_no_NAC != "").sum()),
                    contraste_tandem=contraste))
    json.dump(res, open(f"{COL}/resumen_colinealidad_Egl.json", "w"), ensure_ascii=False, indent=1, default=str)

    print("\nClases en el genoma:", res["genoma"]["clases"])
    print("Clases en los 189 NAC:", res["NAC189"]["clases"])
    print(pd.DataFrame(fis).to_string(index=False, float_format=lambda v: f"{v:.3g}"))
    print("Pares NAC:", json.dumps(res["NAC189"]["pares"], ensure_ascii=False))
    for k, v in contraste.items():
        print(f"{k}:", json.dumps(v, ensure_ascii=False, default=str))
    print(f"\nEscritos: {corto(T45)} ({len(t45)} filas), {corto(T46)} ({len(t46)} filas),"
          f" {corto(COL)}/Egl_pares_colineales.tsv ({len(an)} filas), resumen_colinealidad_Egl.json", flush=True)


if __name__ == "__main__":
    main()
