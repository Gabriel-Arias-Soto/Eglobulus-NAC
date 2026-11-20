"""Entrega 3, paso «Leer y verificar la salida de PlantCARE».

Los 189 promotores NAC de 2 kb se enviaron a PlantCARE en 4 lotes (cis_plantcare/envios/lote_{1..4}.fa). Los
resultados (un .tab por lote y un HTML por gen) están en cis_plantcare/resultados/PlantCARE_<trabajo>_<lote>_plantCARE.tar.gz.

1. Entrada: los lotes son los registrados (md5); un HTML por gen y los mismos genes en el HTML, en el .tab y en el
   FASTA de cada lote; las secuencias de los lotes son las del FASTA de los 189 (2000 pb, solo ACGT).
2. Sitios: cada fila con nombre se comprueba en la secuencia enviada, en la posición y la hebra declaradas.
   Puntaje (5.º campo del .tab) = bases fijas que coinciden + 0,5 x posiciones degeneradas que coinciden. Se exige
   puntaje calculado = declarado y >= 0,9 del máximo del registro. Posición = base más a la izquierda en la hebra
   enviada (1 = -2000 del ATG), también para los sitios de la hebra -.
3. Regla de búsqueda: con los registros (elemento, secuencia) devueltos se recorren las dos hebras de los 189
   promotores y se marcan las posiciones con puntaje >= 0,9 del máximo. PlantCARE informa una fila por gen, elemento
   y posición (un solo registro y una sola hebra aunque coincidan varios). Se exige que el conjunto
   (gen, elemento, posición) del escaneo sea idéntico al de PlantCARE y se acota el umbral (mayor fracción no
   informada < 0,9 <= menor informada). Esta es la regla que aplica e3_17 a todos los promotores del genoma.
4. Filas sin nombre: PlantCARE devuelve filas con los campos de marcador (motif_sequence, organism, short_function)
   y 4 pb. No tienen registro ni función: se excluyen y se informan aparte.
5. Clases: las 5 del artículo 1 (luz, desarrollo, fitohormonas, estrés, otros), según la leyenda de su Figura 6
   (65 elementos). Primero por nombre exacto y después por nombre normalizado (sin espacios, guiones ni guiones bajos,
   en minúsculas). Los elementos con función en PlantCARE que no están en el artículo 1 se asignan según el texto de
   su función (CLASE_EXTRA). Los que no tienen función en PlantCARE (o tienen «?») quedan como «Sin función anotada»,
   fuera de las clases, como en el artículo 1. TATA-box y CAAT-box se informan aparte («Núcleo del promotor»).
   Unidad de análisis = nombre de PlantCARE, con dos excepciones: G-box y G-Box se funden en «G-box» (el mismo
   elemento con dos nombres; se cuentan las posiciones distintas) y los nombres Unnamed_* se sustituyen por su
   registro, p. ej. «Unnamed (GAATTTAATTAA)»: PlantCARE usa un mismo nombre para motivos distintos (Unnamed__1 tiene
   3 registros) y nombres distintos para un mismo motivo (Unnamed__8, __10, __12 y __14 = TCCACGTAGA).
   El texto de función se limpia (comillas sueltas y espacios dobles de la salida de PlantCARE).

Entradas: datos/entrega3/cis_plantcare/{resultados/*.tar.gz, envios/lote_{1..4}.fa, envios/trabajos_plantcare.tsv};
          datos/entrega3/promotores/NAC189_promotores_2kb.fa; tablas/T40_promotores_NAC_2kb.csv,
          T45_modos_duplicacion_NAC189.csv (subfamilia)
Salidas:  tablas/T54_elementos_cis_PlantCARE_clases.csv      unidad -> clase -> función; sitios y genes
          tablas/T55_sitios_PlantCARE_NAC189.csv             un sitio por fila (gen, unidad, posición)
          tablas/T56_recuentos_PlantCARE_por_gen_NAC189.csv  189 genes x (clases, núcleo, unidades)
          datos/entrega3/cis_plantcare/registros_PlantCARE.tsv          registros para el escaneo del genoma (e3_17)
          datos/entrega3/cis_plantcare/resumen_lectura_PlantCARE.json   controles y cifras
Uso: python e3_16_lectura_plantcare.py   (python de nac3; desde cualquier carpeta). Las funciones de la regla
     (tokens, matriz, codificar, mascara) las importa e3_17.
"""
import glob
import hashlib
import json
import os
import re
import tarfile
import time

import numpy as np
import pandas as pd

BASE = (os.environ.get("NAC_ROOT", ".") + "")
A = os.path.join(BASE, "ANALISIS_Eglobulus")
E3 = os.path.join(A, "datos", "entrega3")
CP = os.path.join(E3, "cis_plantcare")
PROM = os.path.join(E3, "promotores")
TAB = os.path.join(A, "tablas")
LARGO = 2000
NUM, DEN = 9, 10                     # umbral de PlantCARE: puntaje >= 9/10 del máximo del registro

# Artículo 1, leyenda de la Figura 6: elementos por clase
ART1 = {
    "Luz": ["3-AF1 binding site", "AAAC-motif", "ACE", "AE-box", "AT1-motif", "ATC-motif", "ATCT-motif", "Box 4",
            "Box II", "C-box", "chs-CMA1a", "chs-CMA2a", "chs-CMA2b", "chs-Unit 1 m1", "GA-motif", "Gap-box",
            "GATA-motif", "G-Box", "GGA-motif", "GT1-motif", "GTGGC-motif", "I-box", "LAMP-element", "L-box", "LS7",
            "MRE", "Sp1", "TCCC-motif", "TCT-motif"],
    "Fitohormonas": ["ABRE", "AuxRE", "AuxRR-core", "CGTCA-motif", "GARE-motif", "P-box", "TATC-box", "TCA-element",
                     "TGA-box", "TGACG-motif", "TGA-element"],
    "Desarrollo": ["3-AF3 binding site", "5UTR Py-rich stretch", "AACA_motif", "ACA-motif", "CAT-box", "circadian",
                   "GCN4_motif", "HD-Zip 1", "HD-Zip 3", "MSA-like", "O2-site", "RY-element"],
    "Estrés": ["ARE", "AT-rich sequence", "GC-motif", "LTR", "MBS", "MBSI", "TC-rich repeats", "WUN-motif"],
    "Otros": ["A-box", "AT-rich element", "Box II-like sequence", "Box III", "CCAAT-box"],
}
# Elementos con función en PlantCARE que no están en el artículo 1 (clase según el texto de la función)
CLASE_EXTRA = {
    "CAG-motif": "Luz",                     # part of a light response element
    "GATT-motif": "Luz",                    # part of a light responsive element
    "Pc-CMA2c": "Luz",                      # part of a light responsive element
    "sbp-CMA1c": "Luz",                     # part of a light responsive element
    "Unnamed (GGATTTTACAGT)": "Luz",        # Unnamed_1: cis-acting element involved in phytochrome down-regulation expression
    "DRE": "Estrés",                        # cis-acting element involved in dehydration, low-temp, salt stresses
    "NON-box": "Desarrollo",                # cis-acting regulatory element related to meristem specific activation
    "motif I": "Desarrollo",                # cis-acting regulatory element root specific
    "Unnamed (GAATTTAATTAA)": "Otros",      # Unnamed__1: 60K protein binding site
}
NUCLEO = ["TATA-box", "CAAT-box"]
FUSION = {"G-Box": "G-box"}
SIN_FUNCION = {"", "?"}
CLASES = ["Luz", "Desarrollo", "Fitohormonas", "Estrés", "Otros"]
NUC, SINF = "Núcleo del promotor", "Sin función anotada"
CLASE_EN = {"Luz": "Light-responsive", "Desarrollo": "Development-related", "Fitohormonas": "Phytohormone-responsive",
            "Estrés": "Stress-responsive", "Otros": "Other functions", NUC: "Core promoter", SINF: "No annotated function"}
GRAM = re.compile(r"(?:\([ACGTacgt](?:/[ACGTacgt])+\)|[ACGTacgt])+")
BASE4 = {"A": 0, "C": 1, "G": 2, "T": 3}
LUT = np.full(256, 4, np.uint8)
for _b, _i in BASE4.items():
    LUT[ord(_b)] = LUT[ord(_b.lower())] = _i


def leer_fasta(p):
    d, k = {}, None
    for l in open(p, encoding="ascii"):
        l = l.rstrip("\n")
        if l.startswith(">"):
            k = l[1:].split()[0]
            assert k not in d, f"{p}: cabecera repetida {k}"
            d[k] = []
        elif k is not None:
            d[k].append(l.strip())
    return {k: "".join(v) for k, v in d.items()}


def norm(s):
    return re.sub(r"[\s_\-]+", "", s).lower()


def tokens(m):
    """Posiciones del registro: base fija (mayúscula o minúscula) o grupo degenerado «(A/G)»."""
    t = re.findall(r"\(([ACGTacgt](?:/[ACGTacgt])+)\)|([ACGTacgt])", m)
    assert "".join(f"({d})" if d else f for d, f in t) == m, f"registro fuera de la gramática: {m}"
    return t


def matriz(m):
    """Pesos x2 (enteros): base fija 2, cada base de un grupo degenerado 1 (PlantCARE: 1 y 0,5). Columna 4 = N."""
    t = tokens(m)
    W = np.zeros((len(t), 5), np.int16)
    for i, (deg, fijo) in enumerate(t):
        if deg:
            for b in deg.upper().split("/"):
                W[i, BASE4[b]] = 1
        else:
            W[i, BASE4[fijo.upper()]] = 2
    return W


def rc_matriz(W):
    return W[::-1][:, [3, 2, 1, 0, 4]]


def codificar(seqs):
    """Lista de secuencias de igual largo -> matriz uint8 (A0 C1 G2 T3, otra cosa 4)."""
    L = len(seqs[0])
    assert all(len(s) == L for s in seqs)
    return LUT[np.frombuffer("".join(seqs).encode("ascii"), np.uint8)].reshape(len(seqs), L)


def puntajes(X, W):
    """Puntaje x2 en cada posición de inicio (n x (L - ancho + 1))."""
    n = X.shape[1] - len(W) + 1
    S = W[0][X[:, :n]].astype(np.int16)
    for k in range(1, len(W)):
        S += W[k][X[:, k:k + n]]
    return S


def mascara(X, motivos, exacto=False):
    """Posiciones de inicio (base izquierda de la hebra enviada) donde algún registro de la unidad alcanza el umbral
    de PlantCARE en alguna de las dos hebras (exacto=True: solo el puntaje máximo)."""
    M = np.zeros(X.shape, bool)
    for m in motivos:
        W = matriz(m)
        mx = int(W.max(1).sum())
        for w in (W, rc_matriz(W)):
            S = puntajes(X, w)
            M[:, :S.shape[1]] |= (S == mx) if exacto else (DEN * S >= NUM * mx)
    return M


def unidad(elemento, registro):
    e = FUSION.get(elemento, elemento)
    return f"Unnamed ({registro})" if e.startswith("Unnamed") else e


def limpiar_funcion(f):
    return re.sub(r"\s+", " ", f.strip().strip('"').strip())


def rc(s):
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def orden_simbolo(s):
    return int(re.sub(r"\D", "", s))


def main():
    t0 = time.time()
    res = {"fecha": time.strftime("%Y-%m-%d %H:%M:%S")}
    reg = pd.read_csv(os.path.join(CP, "envios", "trabajos_plantcare.tsv"), sep="\t")
    nac = leer_fasta(os.path.join(PROM, "NAC189_promotores_2kb.fa"))
    t40 = pd.read_csv(os.path.join(TAB, "T40_promotores_NAC_2kb.csv"))
    t45 = pd.read_csv(os.path.join(TAB, "T45_modos_duplicacion_NAC189.csv"))
    s2l, s2sub = dict(zip(t40.simbolo, t40.locus)), dict(zip(t40.simbolo, t40.subfamilia_final))
    assert len(nac) == 189 and set(nac) == set(s2l) == set(s2sub) == set(t45.simbolo)
    assert all(s2sub[s] == x for s, x in zip(t45.simbolo, t45.subfamilia)), "subfamilia de T40 != T45"

    # 1. Lotes, HTML y .tab
    tars = sorted(glob.glob(os.path.join(CP, "resultados", "PlantCARE_*_*_plantCARE.tar.gz")))
    assert len(tars) == 4, tars
    filas, env, lotes = [], {}, {}
    for f in tars:
        job, k = re.search(r"PlantCARE_(\d+)_(\d)_plantCARE\.tar\.gz$", f).groups()
        k = int(k)
        fa = os.path.join(CP, "envios", f"lote_{k}.fa")
        lote = leer_fasta(fa)
        md5 = hashlib.md5(open(fa, "rb").read()).hexdigest()
        r = reg[reg.lote == k].iloc[0]
        with tarfile.open(f) as tf:
            nombres = tf.getnames()
            html = sorted(os.path.basename(n)[:-len("_plantcare.html")] for n in nombres if n.endswith("_plantcare.html"))
            tabn = [n for n in nombres if n.endswith(".tab")]
            assert len(tabn) == 1, f"lote {k}: {len(tabn)} .tab"
            raw = tf.extractfile(tabn[0]).read().decode("utf-8")
        lin = [x.split("\t") for x in raw.splitlines() if x.strip()]
        assert all(len(x) == 8 for x in lin), f"lote {k}: filas sin 8 campos"
        assert md5 == r.md5 and int(r.secuencias_enviadas) == int(r.secuencias_aceptadas) == len(lote), f"lote {k}: md5 o número de secuencias"
        assert html == sorted(lote), f"lote {k}: genes del HTML != FASTA del lote"
        assert {x[0] for x in lin} == set(lote), f"lote {k}: genes del .tab != FASTA del lote"
        assert not set(lote) & set(env), f"lote {k}: genes repetidos entre lotes"
        env.update(lote)
        filas += [[k, f"PlantCARE_{job}"] + x for x in lin]
        lotes[k] = {"archivo": os.path.basename(f), "trabajo_resultado": f"PlantCARE_{job}",
                    "trabajo_registro": r.trabajo, "md5_lote": md5, "md5_igual_al_registro": True,
                    "secuencias": len(lote), "html": len(html), "filas_tab": len(lin)}
    assert set(env) == set(nac) and all(env[g].upper() == nac[g].upper() for g in nac)
    seqs = {g: s.upper() for g, s in env.items()}
    assert all(len(s) == LARGO and set(s) <= set("ACGT") for s in seqs.values())
    res["lotes"] = lotes
    res["entrada"] = {"genes": len(seqs), "secuencias_lotes_iguales_a_NAC189_fa": True, "largo": LARGO, "solo_ACGT": True,
                      "nota": "los resultados de los trabajos del envío automático (PlantCARE_10418, _11195, _11962, "
                              "_12688) no llegaron por correo; los mismos 4 lotes se volvieron a subir desde la web "
                              "(PlantCARE_27152, _27588, _27843, _28051). La identidad se comprueba por el contenido: md5 "
                              "de los lotes, genes de cada lote y puntaje de cada sitio en la secuencia enviada"}
    print(f"lotes: {[(k, v['trabajo_resultado'], v['secuencias'], v['filas_tab']) for k, v in lotes.items()]}")

    pc = pd.DataFrame(filas, columns=["lote", "trabajo", "gen", "elemento", "registro", "posicion", "puntaje", "hebra",
                                      "organismo", "funcion"])
    marc = (pc.elemento == "") & (pc.registro == "motif_sequence")
    assert (pc.loc[marc, "organismo"] == "organism").all() and (pc.loc[marc, "funcion"] == "short_function").all()
    ms, sit = pc[marc].copy(), pc[~marc].copy()
    assert (sit.elemento != "").all() and sit.registro.map(lambda m: bool(GRAM.fullmatch(m))).all()
    sit["funcion"] = sit.funcion.map(limpiar_funcion)
    assert (sit.groupby(["elemento", "registro"])[["funcion", "organismo"]].nunique() == 1).all().all()
    for d in (ms, sit):
        d["posicion"] = d.posicion.astype(int)
        d["puntaje"] = d.puntaje.astype(float)
        assert d.hebra.isin(["+", "-"]).all() and d.posicion.between(1, LARGO).all()

    # 2. Cada sitio en su posición y hebra
    genes = sorted(seqs, key=orden_simbolo)
    gi = {g: i for i, g in enumerate(genes)}
    X = codificar([seqs[g] for g in genes])
    W = {m: matriz(m) for m in sit.registro.unique()}
    mx = {m: int(w.max(1).sum()) for m, w in W.items()}
    sc, frag = [], []
    for g, m, p, h in zip(sit.gen, sit.registro, sit.posicion, sit.hebra):
        w = W[m] if h == "+" else rc_matriz(W[m])
        x = X[gi[g], p - 1:p - 1 + len(w)]
        sc.append(int(w[np.arange(len(w)), x].sum()) if len(x) == len(w) else -1)
        s = seqs[g][p - 1:p - 1 + len(w)]
        frag.append(s if h == "+" else rc(s))
    sit["puntaje_x2"], sit["sitio"] = sc, frag
    sit["max_x2"] = sit.registro.map(mx)
    dentro = sit.puntaje_x2 >= 0
    igual = sit.puntaje_x2 == (2 * sit.puntaje).round().astype(int)
    umbral = DEN * sit.puntaje_x2 >= NUM * sit.max_x2
    exacto = sit.puntaje_x2 == sit.max_x2
    res["sitios"] = {"filas_con_nombre": len(sit), "dentro_de_la_secuencia": int(dentro.sum()),
                     "puntaje_calculado_igual_al_declarado": int(igual.sum()), "puntaje_>=0.9_del_maximo": int(umbral.sum()),
                     "coincidencia_exacta": int(exacto.sum()), "fraccion_exacta": round(float(exacto.mean()), 4),
                     "con_un_fallo_o_degenerada_incompleta": int((~exacto).sum())}
    assert dentro.all() and igual.all() and umbral.all(), res["sitios"]
    k_pc = set(zip(sit.gen, sit.elemento, sit.posicion))
    assert len(k_pc) == len(sit), "más de una fila por gen, elemento y posición"
    print(f"sitios: {len(sit)} filas con nombre; en posición y hebra declaradas con el puntaje declarado y >= 0,9 del "
          f"máximo: {int((dentro & igual & umbral).sum())} ({(dentro & igual & umbral).mean() * 100:.1f} %); exactos {exacto.mean() * 100:.1f} %")

    # 3. Regla de búsqueda: escaneo de los 189 con todos los registros
    k_sc, max_bajo, varias, dos_hebras = set(), 0.0, 0, 0
    for e, mots in sit.groupby("elemento").registro.unique().items():
        M = np.zeros(X.shape, bool)
        C = np.zeros(X.shape, np.int16)
        for m in mots:
            w0 = W[m]
            Sf, Sr = puntajes(X, w0), puntajes(X, rc_matriz(w0))
            hf, hr = DEN * Sf >= NUM * mx[m], DEN * Sr >= NUM * mx[m]
            n = Sf.shape[1]
            M[:, :n] |= hf | hr
            C[:, :n] += hf.astype(np.int16) + hr
            dos_hebras += int((hf & hr).sum())
            S = np.maximum(Sf, Sr)
            b = S[~(hf | hr)]
            if b.size:
                max_bajo = max(max_bajo, float(b.max()) / mx[m])
        varias += int((C > 1).sum())
        gg, jj = np.nonzero(M)
        k_sc |= {(genes[a], e, int(b) + 1) for a, b in zip(gg, jj)}
    min_inf = float((sit.puntaje_x2 / sit.max_x2).min())
    res["regla"] = {"descripcion": "puntaje >= 0,9 del máximo del registro (fija 1, degenerada 0,5), ambas hebras; una fila "
                                   "por gen, elemento y posición",
                    "registros": int(len(W)), "elementos": int(sit.elemento.nunique()),
                    "posiciones_PlantCARE": len(k_pc), "posiciones_escaneo": len(k_sc), "comunes": len(k_pc & k_sc),
                    "solo_PlantCARE": len(k_pc - k_sc), "solo_escaneo": len(k_sc - k_pc),
                    "mayor_fraccion_no_informada": round(max_bajo, 4), "menor_fraccion_informada": round(min_inf, 4),
                    "posiciones_con_varias_coincidencias_registro_hebra": varias,
                    "coincidencias_en_las_dos_hebras_misma_posicion_y_registro": dos_hebras}
    assert k_sc == k_pc and max_bajo < NUM / DEN <= min_inf + 1e-12, res["regla"]
    print(f"regla: escaneo {len(k_sc)} = PlantCARE {len(k_pc)} posiciones (gen, elemento, posición), sin diferencias; "
          f"umbral acotado: {max_bajo:.4f} < 0,9 <= {min_inf:.4f}")

    # 4. Filas sin nombre
    ms["tetramero"] = [seqs[g][p - 1:p + 3] if h == "+" else rc(seqs[g][p - 1:p + 3])
                       for g, p, h in zip(ms.gen, ms.posicion, ms.hebra)]
    res["filas_sin_nombre"] = {"filas": len(ms), "genes": int(ms.gen.nunique()), "puntaje": sorted(ms.puntaje.unique().tolist()),
                               "tetrameros_distintos": int(ms.tetramero.nunique()),
                               "tetrameros": ms.tetramero.value_counts().to_dict(),
                               "decision": "excluidas: sin nombre, registro ni función en la salida de PlantCARE"}

    # 5. Unidades y clases
    sit["unidad"] = [unidad(e, m) for e, m in zip(sit.elemento, sit.registro)]
    e2c = {e: c for c, l in ART1.items() for e in l}
    nart = {norm(e): e for e in e2c}
    assert len(e2c) == 65 and len(nart) == 65
    fun_u = sit.groupby("unidad").funcion.agg(lambda x: sorted(set(v for v in x if v not in SIN_FUNCION)))
    falta_extra = set(CLASE_EXTRA) - set(fun_u.index)
    assert not falta_extra, f"CLASE_EXTRA con unidades que no están en la salida: {falta_extra}"
    clase, fuente, nom1 = {}, {}, {}
    for u, fs in fun_u.items():
        if u in NUCLEO:
            clase[u], fuente[u] = NUC, "núcleo del promotor (se informa aparte)"
        elif u in e2c:
            clase[u], fuente[u], nom1[u] = e2c[u], "artículo 1, nombre exacto", u
        elif norm(u) in nart:
            nom1[u] = nart[norm(u)]
            clase[u], fuente[u] = e2c[nom1[u]], "artículo 1, nombre normalizado"
        elif u in CLASE_EXTRA:
            assert fs, f"{u}: en CLASE_EXTRA pero sin función en PlantCARE"
            clase[u], fuente[u] = CLASE_EXTRA[u], "texto de función de PlantCARE"
        elif not fs:
            clase[u], fuente[u] = SINF, "sin función en PlantCARE"
        else:
            raise SystemExit(f"unidad con función en PlantCARE y sin clase: {u}: {fs}")
    sit["clase"] = sit.unidad.map(clase)
    ausentes_art1 = sorted(set(e2c) - set(nom1.values()))

    # Registros con más de un nombre (además de los fundidos en una unidad)
    rn = sit.groupby("registro").unidad.agg(lambda x: sorted(set(x)))
    compartidos = {m: {"unidades": u, "clases": sorted({clase[x] for x in u})} for m, u in rn.items() if len(u) > 1}

    # Sitios por (gen, unidad, posición): si varias filas caen en la misma posición (G-box y G-Box, o varios nombres
    # Unnamed con el mismo registro), queda la de mayor fracción del máximo y, a igualdad, la primera del .tab
    sit["simbolo"] = sit.gen
    sit["_o"] = np.arange(len(sit))
    sit["_f"] = sit.puntaje_x2 / sit.max_x2
    nombres_pos = sit.groupby(["gen", "unidad", "posicion"]).elemento.agg(lambda x: ";".join(sorted(set(x))))
    t55 = sit.sort_values(["_f", "_o"], ascending=[False, True]).drop_duplicates(["gen", "unidad", "posicion"]).copy()
    t55["elementos_PlantCARE"] = [nombres_pos[(g, u, p)] for g, u, p in zip(t55.gen, t55.unidad, t55.posicion)]
    fusion = int(len(sit) - len(t55))
    fusion_por_unidad = (sit.groupby("unidad").size() - t55.groupby("unidad").size()).loc[lambda x: x > 0].to_dict()
    t55["locus"] = t55.simbolo.map(s2l)
    t55["subfamilia"] = t55.simbolo.map(s2sub)
    t55["clase_en"] = t55.clase.map(CLASE_EN)
    t55["ancho"] = t55.registro.map(lambda m: len(tokens(m)))
    t55["inicio_rel_ATG"] = t55.posicion - (LARGO + 1)
    t55["fin_rel_ATG"] = t55.posicion + t55.ancho - (LARGO + 2)
    t55["puntaje_max"] = t55.max_x2 / 2
    t55["exacto"] = t55.puntaje_x2 == t55.max_x2
    t55["_g"] = t55.simbolo.map(orden_simbolo)
    t55 = t55.sort_values(["_g", "posicion", "unidad"])
    t55 = t55[["simbolo", "locus", "subfamilia", "clase", "clase_en", "unidad", "elementos_PlantCARE", "registro", "sitio",
               "posicion", "inicio_rel_ATG", "fin_rel_ATG", "hebra", "ancho", "puntaje", "puntaje_max", "exacto",
               "organismo", "funcion", "lote", "trabajo"]]
    assert len(set(zip(t55.simbolo, t55.unidad, t55.posicion))) == len(t55)

    # Recuentos por gen
    orden_c = {c: i for i, c in enumerate(CLASES + [NUC, SINF])}
    unidades = sorted(clase, key=lambda u: (orden_c[clase[u]], u.lower()))
    cnt = pd.crosstab(t55.simbolo, t55.unidad).reindex(index=genes, columns=unidades, fill_value=0)
    t56 = pd.DataFrame({"simbolo": genes, "locus": [s2l[g] for g in genes], "subfamilia": [s2sub[g] for g in genes]})
    for c in CLASES:
        t56[c] = cnt[[u for u in unidades if clase[u] == c]].sum(1).to_numpy()
    t56["total_5_clases"] = t56[CLASES].sum(1)
    for u in NUCLEO:
        t56[u] = cnt[u].to_numpy()
    t56["sin_funcion_anotada"] = cnt[[u for u in unidades if clase[u] == SINF]].sum(1).to_numpy()
    for u in unidades:
        if u not in NUCLEO:
            t56[u] = cnt[u].to_numpy()
    assert (t56[CLASES + NUCLEO + ["sin_funcion_anotada"]].sum(1).to_numpy() == cnt.sum(1).to_numpy()).all()

    # Tabla de elementos
    ag = sit.groupby("unidad").agg(elementos_PlantCARE=("elemento", lambda x: ";".join(sorted(set(x)))),
                                   registros=("registro", lambda x: ";".join(sorted(set(x)))),
                                   n_registros=("registro", "nunique"),
                                   organismos=("organismo", lambda x: ";".join(sorted(set(x)))))
    t54 = pd.DataFrame({"unidad": unidades})
    t54["clase"] = t54.unidad.map(clase)
    t54["clase_en"] = t54.clase.map(CLASE_EN)
    t54["fuente_clase"] = t54.unidad.map(fuente)
    t54["nombre_articulo1"] = t54.unidad.map(lambda u: nom1.get(u, ""))
    t54 = t54.join(ag, on="unidad")
    t54["funcion_PlantCARE"] = t54.unidad.map(lambda u: " | ".join(fun_u[u]))
    t54["sitios_NAC189"] = t54.unidad.map(cnt.sum(0))
    t54["genes_con_sitio"] = t54.unidad.map((cnt > 0).sum(0))
    t54["media_por_gen"] = t54.unidad.map(cnt.mean(0)).round(3)
    t54["sitios_exactos"] = t54.unidad.map(t55.groupby("unidad").exacto.sum()).astype(int)

    # Registros para el escaneo del genoma
    rg = sit.groupby(["unidad", "elemento", "registro"]).agg(sitios_informados=("gen", "size"),
                                                           funcion=("funcion", "first"),
                                                           organismo=("organismo", "first")).reset_index()
    assert (sit.groupby(["elemento", "registro"]).funcion.nunique() == 1).all()
    rg["clase"] = rg.unidad.map(clase)
    rg["ancho"] = rg.registro.map(lambda m: len(tokens(m)))
    rg["n_degeneradas"] = rg.registro.map(lambda m: sum(1 for d, _ in tokens(m) if d))
    rg["puntaje_max"] = rg.registro.map(mx) / 2
    rg = rg[["unidad", "clase", "elemento", "registro", "ancho", "n_degeneradas", "puntaje_max", "sitios_informados",
             "funcion", "organismo"]].sort_values(["unidad", "registro"])

    t54.to_csv(os.path.join(TAB, "T54_elementos_cis_PlantCARE_clases.csv"), index=False)
    t55.to_csv(os.path.join(TAB, "T55_sitios_PlantCARE_NAC189.csv"), index=False)
    t56.to_csv(os.path.join(TAB, "T56_recuentos_PlantCARE_por_gen_NAC189.csv"), index=False)
    rg.to_csv(os.path.join(CP, "registros_PlantCARE.tsv"), sep="\t", index=False)

    por_clase = t54.groupby("clase").agg(unidades=("unidad", "size"), sitios=("sitios_NAC189", "sum")).reindex(CLASES + [NUC, SINF])
    res["clases"] = {"unidades": len(unidades), "por_clase": {c: {"unidades": int(r.unidades), "sitios": int(r.sitios)}
                                                              for c, r in por_clase.iterrows()},
                     "por_fuente": t54.fuente_clase.value_counts().to_dict(),
                     "emparejados_por_nombre_normalizado": {u: nom1[u] for u in nom1 if u != nom1[u]},
                     "articulo1_ausentes_aqui": ausentes_art1, "clase_extra": CLASE_EXTRA,
                     "sin_funcion": {u: int(cnt[u].sum()) for u in unidades if clase[u] == SINF},
                     "filas_fundidas_misma_posicion": fusion, "filas_fundidas_por_unidad": fusion_por_unidad,
                     "sitios_tras_fusion": len(t55),
                     "registros_en_varias_unidades": compartidos,
                     "nota_registros_compartidos": "un mismo registro bajo nombres distintos se cuenta en cada unidad "
                                                   "(como en la salida de PlantCARE y en el artículo 1)"}
    res["segundos"] = round(time.time() - t0, 1)
    json.dump(res, open(os.path.join(CP, "resumen_lectura_PlantCARE.json"), "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    print(f"filas sin nombre: {len(ms)} en {ms.gen.nunique()} genes ({ms.tetramero.nunique()} tetrámeros); excluidas")
    print(f"unidades: {len(unidades)}; fuentes: {res['clases']['por_fuente']}")
    print(f"por nombre normalizado: {res['clases']['emparejados_por_nombre_normalizado']}; del artículo 1 ausentes aquí: {ausentes_art1}")
    print(f"filas en la misma posición fundidas: {fusion} {fusion_por_unidad}; sitios {len(sit)} -> {len(t55)}")
    print(f"registros con más de una unidad: {len(compartidos)}; con unidades de clases distintas: "
          f"{sum(len(v['clases']) > 1 for v in compartidos.values())}")
    print(por_clase.to_string())
    print(f"listo en {time.time() - t0:.1f} s")


if __name__ == "__main__":
    main()
