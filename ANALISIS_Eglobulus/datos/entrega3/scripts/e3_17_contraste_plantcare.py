"""Entrega 3, elementos cis de PlantCARE (artículo 1), paso 2: contraste con el fondo emparejado por GC.

Escaneo del genoma. Los registros (secuencias consenso) que PlantCARE informó en los 189 promotores NAC
(registros_PlantCARE.tsv de e3_16) se buscan en ambas hebras de los 35.943 promotores de 2 kb con la regla de PlantCARE
que e3_16 reprodujo sin diferencias en los 189: puntaje >= 0,9 del máximo del registro (base fija 1, cada base de un
grupo degenerado 0,5), con las funciones de e3_16. Sitio = posición de inicio (base izquierda de la hebra enviada) donde
algún registro de la unidad alcanza el umbral en alguna hebra: una posición por gen y unidad, como en T55 (G-box y G-Box,
y los Unnamed con el mismo registro, cuentan una vez). Se guardan también los recuentos por registro y, para las unidades
con varios registros, qué registros dan cada sitio (para la calibración con selección).
Controles: en los 189 NAC el escaneo da los mismos sitios que T55 (gen, unidad, posición) y los mismos recuentos que T56
(por unidad y por clase) y T54; los sitios exactos de T55 están entre los exactos del escaneo, y los que el escaneo añade
son sitios de T55 donde PlantCARE informó un registro no exacto y otro registro de la unidad es exacto (T55 y T54 marcan
el registro informado, uno por elemento y posición; el escaneo, cualquier registro de la unidad, igual en NAC y fondo);
todos los registros tienen sitio en los NAC; recuentos por registro, por unidad y marcas de registro concuerdan entre sí.

Contraste. Los 189 NAC se comparan con N conjuntos de 189 promotores de fondo de e3_13 (cada promotor de fondo en el
mismo intervalo de GC que su NAC; fondo = 35.754 promotores no NAC). Principal: N = 10.000, como en e3_15 (los primeros
1000 son la colección original; con ellos se repite el principal como control de resolución).
Medidas por promotor:
  unidad: número de sitios (principal); presencia (>= 1 sitio); sitios exactos (algún registro con el puntaje máximo);
  clase:  suma de los sitios de sus unidades (definición de T56 y del artículo 1). Grupos: las 5 clases del artículo 1,
          núcleo del promotor (TATA-box, CAAT-box) y sin función anotada.
Estadístico de grupo = suma en los 189; razón, z, p empírica bilateral, q de BH (dentro de cada nivel y de cada
  análisis) y dirección con las funciones de e3_15.
Sensibilidad (unidades): presencia; sitios exactos; un gen por bloque en tándem (T45.bloque_T15, como e3_15); sin
  redundancia, un promotor por grupo de promotores con secuencia compartida (20-mers canónicos sin baja complejidad
  [<= 2 bases distintas] ni presentes en más de 10 promotores NAC; dos promotores se unen si comparten más del 10 % de
  los 20-mers del que tiene menos; enlace simple; queda el primero por símbolo; capta también duplicados que no están en
  un mismo bloque en tándem); prueba condicionada a presencia (referencia = conjuntos de fondo con >= 1 sitio: una
  unidad solo se prueba si PlantCARE la halló en los NAC). Los promotores que quedan conservan su columna de fondo
  emparejada. Clases: exactos, sin tándem y sin redundancia.
Calibración (N_CALIB conjuntos de fondo como falso «NAC» frente a los N - 1 restantes, como e3_15):
  a) sin selección: todas las unidades con todos sus registros;
  b) con selección, como el análisis real: solo los registros con sitio en el falso NAC y solo las unidades con algún
     registro así (los recuentos de las unidades con registros no vistos se recalculan en todo el genoma); unidades con
     la prueba principal y la condicionada, por estratos de sitios en el falso NAC; clases con esa misma selección.
Por gen: para cada NAC y grupo (y cada unidad, en un .npz), sitios frente a la media de sus N promotores emparejados
  (columna j de idx): log2((obs + 1) / (esperado + 1)) y percentil en su fondo (como T53).
Limitación: los registros de PlantCARE sin sitio en los 189 no se prueban (la base no se puede descargar): no se detecta
  el empobrecimiento de elementos ausentes y la selección de registros y unidades por presencia favorece el
  enriquecimiento de los raros; la calibración b) lo mide dentro del universo de registros conocidos.

Uso: python e3_17_contraste_plantcare.py [--n N] [--procesos P] [--reusar_escaneo]   (python de nac3)
Entradas: datos/entrega3/promotores/{promotores_2kb_todos.fa, NAC189_promotores_2kb.fa};
          datos/entrega3/cis_plantcare/{registros_PlantCARE.tsv, fondo/fondo_GC_conjuntos[_N].npz,
          fondo/resumen_fondo_GC_10000.json}; datos/entrega3/cis_jaspar/resumen_contraste_JASPAR2026.json;
          tablas/{T40, T45, T54, T55, T56}
Salidas (con N distinto de 10.000 llevan el sufijo _N<N>): tablas/T57_enriquecimiento_PlantCARE_elementos.csv,
         T58_enriquecimiento_PlantCARE_clases.csv, T59_clases_PlantCARE_por_gen_NAC189.csv;
         datos/entrega3/cis_plantcare/{escaneo_PlantCARE_genoma.npz, fondo_sumas_PlantCARE.npz,
         esperado_por_gen_PlantCARE.npz, resumen_contraste_PlantCARE.json}
"""
import argparse, json, os, platform, sys, time
import multiprocessing as mp
import numpy as np
import pandas as pd
import scipy
from collections import Counter
from scipy.stats import spearmanr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from e3_16_lectura_plantcare import (leer_fasta, matriz, rc_matriz, codificar, puntajes, orden_simbolo,  # noqa: E402
                                     NUM, DEN, LARGO, CLASES, NUC, SINF)
from e3_15_contraste_jaspar import bh, sumas_fondo, contraste, concordantes, direccion  # noqa: E402

BASE = (os.environ.get("NAC_ROOT", ".") + "")
A = os.path.join(BASE, "ANALISIS_Eglobulus")
E3 = os.path.join(A, "datos", "entrega3")
CP = os.path.join(E3, "cis_plantcare")
PROM = os.path.join(E3, "promotores")
TAB = os.path.join(A, "tablas")
N_NAC = 189
N_TROZOS = 64
N_CALIB = 100
SEMILLA = 20261018
COL_TANDEM = "bloque_T15"          # definición autoritativa (T15), la misma de e3_15
GRUPOS = CLASES + [NUC, SINF]
ESTRATOS = [(1, 2), (3, 9), (10, 49), (50, 10 ** 9)]
K_RED = 20                         # redundancia entre promotores NAC: 20-mers canónicos compartidos
UMBRAL_RED = 0.1                   # se unen dos promotores si comparten más del 10 % de los 20-mers del que tiene menos
COMUN_RED = 10                     # se quitan los 20-mers presentes en más de 10 promotores NAC (repeticiones dispersas)

_X, _REG, _RSEQ = None, None, None   # promotores codificados, registros por unidad y secuencias (a los procesos por fork)


def mascara_registro(X, m):
    """Inicios (base izquierda de la hebra enviada) donde el registro alcanza el umbral de PlantCARE en alguna hebra,
    y donde alcanza el máximo (exacto)."""
    W = matriz(m)
    mx = int(W.max(1).sum())
    Mt = np.zeros(X.shape, bool)
    Me = np.zeros(X.shape, bool)
    for w in (W, rc_matriz(W)):
        S = puntajes(X, w)
        n = S.shape[1]
        Mt[:, :n] |= DEN * S >= NUM * mx
        Me[:, :n] |= S == mx
    return Mt, Me


def _trozo(rango):
    a, b = rango
    Xc = _X[a:b]
    st = np.zeros((b - a, len(_REG)), np.int16)
    ex = np.zeros((b - a, len(_REG)), np.int16)
    cr = np.zeros((b - a, len(_RSEQ)), np.int16)
    ef, ep, ej, eb = [], [], [], []
    for j, regs in enumerate(_REG):
        Ut = np.zeros(Xc.shape, bool)
        Ue = np.zeros(Xc.shape, bool)
        B = np.zeros(Xc.shape, np.int32) if len(regs) > 1 else None
        for k, r in enumerate(regs):
            mt, me = mascara_registro(Xc, _RSEQ[r])
            cr[:, r] = mt.sum(axis=1)
            Ut |= mt
            Ue |= me
            if B is not None:
                B[mt] |= np.int32(1 << k)
        st[:, j], ex[:, j] = Ut.sum(axis=1), Ue.sum(axis=1)
        if B is not None:
            f, p = np.nonzero(B)
            ef.append((f + a).astype(np.int32)); ep.append(p.astype(np.int16))
            ej.append(np.full(len(f), j, np.int16)); eb.append(B[f, p])
    return a, st, ex, cr, tuple(np.concatenate(x) for x in (ef, ep, ej, eb))


def contraste_condicional(obs, S):
    """Como contraste() de e3_15, condicionado a >= 1 sitio: la referencia son los conjuntos de fondo con algún sitio
    (la prueba de una unidad solo existe si PlantCARE la halló en el grupo). No evaluable sin conjuntos así."""
    V = S >= 1
    n1 = V.sum(0)
    p_enr = (1 + ((S >= obs) & V).sum(0)) / (1 + n1)
    p_emp = (1 + ((S <= obs) & V).sum(0)) / (1 + n1)
    ev = (obs >= 1) & (n1 >= 1)
    p = np.where(ev, np.minimum(1, 2 * np.minimum(p_enr, p_emp)), np.nan)
    q = np.full(len(p), np.nan)
    q[ev] = bh(p[ev])
    with np.errstate(divide="ignore", invalid="ignore"):
        mu = np.where(n1 > 0, (S * V).sum(0) / n1, np.nan)
        de = np.array([S[V[:, k], k].std(ddof=1) if n1[k] > 1 else np.nan for k in range(S.shape[1])])
        razon = np.where(mu > 0, obs / mu, np.nan)
        z = np.where(de > 0, (obs - mu) / de, np.nan)
        l2 = np.log2(razon)
    return {"obs": obs, "media_fondo": mu, "de_fondo": de, "razon": razon, "log2_razon": l2, "z": z,
            "p_enr": p_enr, "p_emp": p_emp, "p": p, "q": q, "evaluable": ev, "conjuntos_con_sitio": n1}


def tabla(base, res, claves, prefijos):
    t = base.copy()
    for clave, pre in zip(claves, prefijos):
        r = res[clave]
        cols_r = ["obs", "media_fondo", "de_fondo", "razon", "log2_razon", "z", "p", "q"] if pre == "" else ["obs", "razon", "p", "q"]
        for c in cols_r:
            t[(pre + c) if pre else c] = r[c]
    return t


def por_gen(M, fila_nac, filas):
    """Observado de cada NAC frente a sus N promotores emparejados (columna j de idx): esperado y percentil (e3_15)."""
    obs = M[fila_nac]
    suma, men, igu = (np.zeros(obs.shape) for _ in range(3))
    N = filas.shape[0]
    for i in range(0, N, 1000):
        fg = M[filas[i:i + 1000]]                         # b x 189 x K
        suma += fg.sum(axis=0)
        men += (fg < obs[None]).sum(axis=0)
        igu += (fg == obs[None]).sum(axis=0)
    return obs, suma / N, (men + 0.5 * igu) / N * 100


def estrato(o):
    return next(f"{a}-{b}" if b < 10 ** 9 else f">={a}" for a, b in ESTRATOS if a <= o <= b)


def grupos_redundantes(seqs, k=K_RED, umbral=UMBRAL_RED, comun=COMUN_RED):
    """Grupos de promotores con secuencia compartida. k-mers canónicos (solo A, C, G, T; sin baja complejidad: <= 2 bases
    distintas, lo que incluye las repeticiones de mono- y dinucleótidos), sin los presentes en más de `comun` promotores;
    dos promotores se unen si comparten más de `umbral` de los k-mers del que tiene menos; grupos por enlace simple.
    Devuelve etiquetas de grupo, matriz de fracciones compartidas y número de k-mers comunes quitados."""
    comp = str.maketrans("ACGT", "TGCA")
    acgt = set("ACGT")
    ks = []
    for s in seqs:
        s = s.upper()
        out = set()
        for i in range(len(s) - k + 1):
            x = s[i:i + k]
            b = set(x)
            if b <= acgt and len(b) > 2:
                out.add(min(x, x.translate(comp)[::-1]))
        ks.append(out)
    cnt = Counter(x for s in ks for x in s)
    comunes = {x for x, c in cnt.items() if c > comun}
    ks = [s - comunes for s in ks]
    n = len(ks)
    F = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            if ks[i] and ks[j]:
                F[i, j] = F[j, i] = len(ks[i] & ks[j]) / min(len(ks[i]), len(ks[j]))
    lab = list(range(n))

    def raiz(x):
        while lab[x] != x:
            lab[x] = lab[lab[x]]
            x = lab[x]
        return x
    for i, j in zip(*np.nonzero(np.triu(F > umbral, 1))):
        lab[raiz(i)] = raiz(j)
    return np.array([raiz(i) for i in range(n)]), F, len(comunes)


def main(n_conj=10000, procesos=None, reusar=False):
    global _X, _REG, _RSEQ
    t0 = time.time()
    procesos = procesos or max(1, min(10, (os.cpu_count() or 2) - 2))
    suf = "" if n_conj == 10000 else f"_N{n_conj}"
    FONDO = os.path.join(CP, "fondo", "fondo_GC_conjuntos.npz" if n_conj == 1000 else f"fondo_GC_conjuntos_{n_conj}.npz")
    ESC = os.path.join(CP, "escaneo_PlantCARE_genoma.npz")

    # 1. Promotores
    t40 = pd.read_csv(os.path.join(TAB, "T40_promotores_NAC_2kb.csv"))
    s2l = dict(zip(t40.simbolo, t40.locus))
    seqs = leer_fasta(os.path.join(PROM, "promotores_2kb_todos.fa"))
    nac_fa = leer_fasta(os.path.join(PROM, "NAC189_promotores_2kb.fa"))
    loci = list(seqs)
    l2i = {x: i for i, x in enumerate(loci)}
    nP = len(loci)
    assert len(nac_fa) == N_NAC and set(nac_fa) == set(s2l) and len(l2i) == nP
    iguales = sum(seqs[s2l[g]].upper() == nac_fa[g].upper() for g in nac_fa)
    assert iguales == N_NAC, f"NAC189_promotores_2kb.fa y promotores_2kb_todos.fa difieren en {N_NAC - iguales} genes"
    X = codificar([seqs[x] for x in loci])
    assert X.shape == (nP, LARGO) and (X < 4).all()
    del seqs

    # 2. Unidades (orden y clase de T54) y sus registros
    t54 = pd.read_csv(os.path.join(TAB, "T54_elementos_cis_PlantCARE_clases.csv"), keep_default_na=False)
    rg = pd.read_csv(os.path.join(CP, "registros_PlantCARE.tsv"), sep="\t", keep_default_na=False)
    unidades = t54.unidad.tolist()
    clase = dict(zip(t54.unidad, t54.clase))
    c2en = dict(zip(t54.clase, t54.clase_en))
    registros = sorted(set(rg.registro))
    r2i = {m: i for i, m in enumerate(registros)}
    reg_idx = [[r2i[m] for m in sorted(set(rg.registro[rg.unidad == u]))] for u in unidades]
    assert len(set(unidades)) == len(unidades) and set(rg.unidad) == set(unidades) and all(reg_idx)
    assert all(clase[u] == c for u, c in zip(rg.unidad, rg.clase)) and set(clase.values()) == set(GRUPOS)
    assert max(len(r) for r in reg_idx) <= 30, "más de 30 registros en una unidad: no caben en las marcas int32"
    nU, nR = len(unidades), len(registros)
    n_pares = int(rg[["unidad", "registro"]].drop_duplicates().shape[0])
    firma = np.array([";".join(registros[r] for r in regs) for regs in reg_idx])
    print(f"promotores {nP}; unidades {nU} ({sum(len(r) > 1 for r in reg_idx)} con varios registros); registros {nR} "
          f"({n_pares} pares unidad-registro)", flush=True)

    # 3. Escaneo del genoma (o el guardado, si coincide con las entradas)
    seg_escaneo = None
    if reusar and os.path.exists(ESC):
        z = np.load(ESC)
        assert (list(z["loci"]) == loci and list(z["unidades"]) == unidades and list(z["registros"]) == registros
                and list(z["registros_por_unidad"]) == list(firma)), f"{ESC} no corresponde a las entradas actuales"
        CT, CE, CR = z["sitios"], z["exactos"], z["por_registro"]
        EF, EP, EJ, EB = z["ent_promotor"], z["ent_posicion"], z["ent_unidad"], z["ent_marcas"]
        print(f"escaneo reutilizado: {os.path.relpath(ESC, BASE)}", flush=True)
    else:
        _X, _REG, _RSEQ = X, reg_idx, registros
        bordes = np.linspace(0, nP, N_TROZOS + 1).astype(int)
        CT = np.full((nP, nU), -1, np.int16)
        CE = np.full((nP, nU), -1, np.int16)
        CR = np.full((nP, nR), -1, np.int16)
        ents = {}
        t1 = time.time()
        print(f"escaneo: {N_TROZOS} trozos en {procesos} procesos", flush=True)
        with mp.get_context("fork").Pool(procesos) as pool:
            for k, (a, st, ex, cr, e) in enumerate(pool.imap_unordered(_trozo, list(zip(bordes[:-1], bordes[1:]))), 1):
                CT[a:a + len(st)], CE[a:a + len(st)], CR[a:a + len(st)] = st, ex, cr
                ents[a] = e
                if k % 8 == 0:
                    print(f"  {k}/{N_TROZOS} trozos ({time.time() - t1:.0f} s)", flush=True)
        seg_escaneo = round(time.time() - t1)
        EF, EP, EJ, EB = (np.concatenate([ents[a][i] for a in sorted(ents)]) for i in range(4))
        np.savez_compressed(ESC, sitios=CT, exactos=CE, por_registro=CR, ent_promotor=EF, ent_posicion=EP, ent_unidad=EJ,
                            ent_marcas=EB, loci=np.array(loci), unidades=np.array(unidades), registros=np.array(registros),
                            registros_por_unidad=firma)
        print(f"escaneo listo en {seg_escaneo} s; guardado {os.path.relpath(ESC, BASE)}", flush=True)
    assert (CT >= 0).all() and (CE >= 0).all() and (CR >= 0).all() and (CE <= CT).all()
    CT32, CE32, CR32 = CT.astype(np.int32), CE.astype(np.int32), CR.astype(np.int32)

    # Coherencia interna: unidad de un registro = ese registro; marcas = sitios de la unidad y de cada registro
    um = [j for j, regs in enumerate(reg_idx) if len(regs) > 1]
    jl = np.full(nU, -1, np.int64)
    jl[um] = np.arange(len(um))
    EJL = jl[EJ.astype(np.int64)]
    assert (EJL >= 0).all() and (EB > 0).all()
    for j, regs in enumerate(reg_idx):
        if len(regs) == 1:
            assert (CT32[:, j] == CR32[:, regs[0]]).all(), unidades[j]
    for l, j in enumerate(um):
        e = EJL == l
        f_, b_ = EF[e], EB[e]
        assert (np.bincount(f_, minlength=nP) == CT32[:, j]).all(), unidades[j]
        for k, r in enumerate(reg_idx[j]):
            assert (np.bincount(f_[((b_ >> k) & 1) == 1], minlength=nP) == CR32[:, r]).all(), (unidades[j], registros[r])
    G = np.zeros((nU, len(GRUPOS)), np.int32)
    for j, u in enumerate(unidades):
        G[j, GRUPOS.index(clase[u])] = 1
    CTc, CEc = CT32 @ G, CE32 @ G

    # 4. Controles en los 189 NAC: mismos sitios que T55 y mismos recuentos que T56 y T54
    t55 = pd.read_csv(os.path.join(TAB, "T55_sitios_PlantCARE_NAC189.csv"), keep_default_na=False)
    t56 = pd.read_csv(os.path.join(TAB, "T56_recuentos_PlantCARE_por_gen_NAC189.csv"), keep_default_na=False)
    genes = t56.simbolo.tolist()
    fg = np.array([l2i[s2l[g]] for g in genes])
    Xn = X[fg]
    kt, ke = set(), set()
    for u, regs in zip(unidades, reg_idx):
        Ut = np.zeros(Xn.shape, bool)
        Ue = np.zeros(Xn.shape, bool)
        for r in regs:
            mt, me = mascara_registro(Xn, registros[r])
            Ut |= mt
            Ue |= me
        kt |= {(genes[a], u, int(b) + 1) for a, b in zip(*np.nonzero(Ut))}
        ke |= {(genes[a], u, int(b) + 1) for a, b in zip(*np.nonzero(Ue))}
    ex55 = t55.exacto.astype(str).str.lower() == "true"
    k55 = set(zip(t55.simbolo, t55.unidad, t55.posicion.astype(int)))
    k55e = set(zip(t55.simbolo[ex55], t55.unidad[ex55], t55.posicion[ex55].astype(int)))
    # Exacto en el escaneo = algún registro de la unidad alcanza el máximo en alguna hebra (la misma regla en NAC y fondo).
    # Exacto en T55 = el registro que informó PlantCARE en esa posición (uno por elemento y posición) alcanza el máximo.
    # Los exactos del escaneo que T55 no marca tienen que ser sitios de T55 con registro informado no exacto en los que
    # otro registro de la unidad sí es exacto; se comprueba uno a uno y se guardan.
    informado = {(g, u, int(p)): (r, h, bool(x)) for g, u, p, r, h, x in
                 zip(t55.simbolo, t55.unidad, t55.posicion, t55.registro, t55.hebra, ex55)}
    g2f = {g: i for i, g in enumerate(genes)}
    u2r = dict(zip(unidades, reg_idx))
    difs_exactos = []
    for g, u, p in sorted(ke - k55e, key=lambda k: (orden_simbolo(k[0]), k[2], k[1])):
        exactos_aqui = []
        for r in u2r[u]:
            W = matriz(registros[r])
            mx = int(W.max(1).sum())
            for h, w in (("+", W), ("-", rc_matriz(W))):
                S = puntajes(Xn[g2f[g]:g2f[g] + 1], w)[0]
                if p - 1 < len(S) and S[p - 1] == mx:
                    exactos_aqui.append(f"{registros[r]} ({h})")
        r_inf, h_inf, x_inf = informado.get((g, u, p), (None, None, None))
        difs_exactos.append({"simbolo": g, "unidad": u, "posicion": p, "registro_informado": r_inf,
                             "hebra_informada": h_inf, "registro_informado_exacto": x_inf,
                             "registros_exactos_en_la_posicion": exactos_aqui})
    difs_u = pd.Series([d["unidad"] for d in difs_exactos], dtype=object).value_counts().to_dict()
    d54 = {u: int(v) for u, v in zip(t54.unidad, CE[fg].sum(axis=0) - t54.sitios_exactos.to_numpy()) if v != 0}
    nuc = t56[["TATA-box", "CAAT-box"]].sum(axis=1).to_numpy()
    controles = {
        "secuencias_NAC_iguales_en_los_dos_FASTA": int(iguales),
        "recuentos_NAC_iguales_a_T56_unidades": bool((CT[fg] == t56[unidades].to_numpy()).all()),
        "recuentos_NAC_iguales_a_T56_clases": bool((CTc[fg][:, :5] == t56[CLASES].to_numpy()).all()
                                                   and (CTc[fg][:, 5] == nuc).all()
                                                   and (CTc[fg][:, 6] == t56.sin_funcion_anotada.to_numpy()).all()),
        "T54_sitios_y_genes_iguales": bool((CT[fg].sum(axis=0) == t54.sitios_NAC189.to_numpy()).all()
                                           and ((CT[fg] > 0).sum(axis=0) == t54.genes_con_sitio.to_numpy()).all()),
        "todos_los_registros_con_sitio_en_los_NAC": bool((CR32[fg].sum(axis=0) > 0).all()),
        "sitios_NAC_escaneo": len(kt), "sitios_T55": len(k55), "sitios_iguales_a_T55": kt == k55,
        "exactos_T55": len(k55e), "exactos_escaneo": len(ke), "exactos_T55_incluidos_en_el_escaneo": k55e <= ke,
        "exactos_escaneo_no_marcados_en_T55": len(ke - k55e),
        "exactos_no_marcados_explicados_por_otro_registro_de_la_unidad": all(
            d["registro_informado"] is not None and d["registro_informado_exacto"] is False
            and len(d["registros_exactos_en_la_posicion"]) > 0 for d in difs_exactos),
        "exactos_matriz_igual_al_escaneo": int(CE[fg].sum()) == len(ke),
        "T54_sitios_exactos_iguales_salvo_esos_sitios": d54 == difs_u,
        "exactos_no_marcados_por_unidad": difs_u,
        "exactos_no_marcados_en_T55": difs_exactos}
    assert all(v for v in controles.values() if isinstance(v, bool)), {k: v for k, v in controles.items()
                                                                      if k != "exactos_no_marcados_en_T55"}
    print(f"controles NAC: sitios = T55 ({len(kt)}); recuentos = T56 (189 x {nU} y clases) y T54; exactos {len(ke)} "
          f"(T55 {len(k55e)}, incluidos; {len(difs_exactos)} sitios más, donde PlantCARE informó un registro no exacto y "
          f"otro registro de la unidad es exacto: {difs_u}); los {nR} registros tienen sitio en los NAC", flush=True)
    print(f"genoma: {int(CT32.sum()):,} sitios ({int(CE32.sum()):,} exactos); {len(EF):,} sitios de unidades con varios "
          f"registros, con sus marcas", flush=True)

    # 5. Fondo emparejado por GC
    fz = np.load(FONDO)
    nac_sim, loc_fondo, idx = list(fz["nac"]), list(fz["locus_fondo"]), fz["idx"]
    assert len(nac_sim) == N_NAC and set(nac_sim) == set(genes)
    fila_nac = np.array([l2i[s2l[s]] for s in nac_sim])
    fila_fondo = np.array([l2i[x] for x in loc_fondo])
    assert len(set(fila_fondo) & set(fila_nac)) == 0 and len(fila_fondo) + len(fila_nac) == nP
    filas = fila_fondo[idx]                      # N x 189, filas globales
    N = filas.shape[0]
    assert N == n_conj, f"{FONDO}: {N} conjuntos, se esperaban {n_conj}"
    rf = json.load(open(os.path.join(CP, "fondo", "resumen_fondo_GC_10000.json"), encoding="utf-8"))
    control_1000 = n_conj == 10000 and rf["primeros_conjuntos_iguales_a_la_coleccion_original"]

    # 6. Contrastes: principal y sensibilidad
    medidas = {"u_sitios": CT32, "u_presencia": (CT32 > 0).astype(np.int32), "u_exactos": CE32,
               "c_sitios": CTc, "c_exactos": CEc}
    res, sumas = {}, {}
    for nom, M in medidas.items():
        sumas[nom] = sumas_fondo(M, filas)
        res[nom] = contraste(M[fila_nac].sum(axis=0), sumas[nom])
    res["u_sitios_condicional"] = contraste_condicional(res["u_sitios"]["obs"], sumas["u_sitios"])

    t45 = pd.read_csv(os.path.join(TAB, "T45_modos_duplicacion_NAC189.csv")).sort_values(["cromosoma", "inicio"])
    bloq = t45.dropna(subset=[COL_TANDEM])
    quitar = set(bloq.simbolo) - set(bloq.groupby(COL_TANDEM).simbolo.first())
    rj = json.load(open(os.path.join(E3, "cis_jaspar", "resumen_contraste_JASPAR2026.json"), encoding="utf-8"))
    assert quitar == set(rj["tandem_quitados"]), "los genes quitados por tándem no son los de e3_15"
    cols = np.array([j for j, s in enumerate(nac_sim) if s not in quitar])
    for nom in ("u_sitios", "c_sitios"):
        M = medidas[nom]
        res[nom + "_sin_tandem"] = contraste(M[fila_nac[cols]].sum(axis=0), sumas_fondo(M, filas[:, cols]))
    print(f"tándem ({COL_TANDEM}): {bloq[COL_TANDEM].nunique()} bloques con {len(bloq)} genes; se quitan {len(quitar)} "
          f"(los mismos de e3_15); quedan {len(cols)} NAC", flush=True)

    # Redundancia de secuencia entre promotores (también fuera de los bloques en tándem): un promotor por grupo, el
    # primero por símbolo; cada promotor que queda conserva su columna de fondo emparejada por GC
    grp, Fr, n_comunes = grupos_redundantes([nac_fa[s] for s in nac_sim])
    miembros = {g: sorted(np.nonzero(grp == g)[0].tolist(), key=lambda k: orden_simbolo(nac_sim[k]))
                for g in set(grp.tolist())}
    cols_r = np.array(sorted(v[0] for v in miembros.values()))
    for nom in ("u_sitios", "c_sitios"):
        M = medidas[nom]
        res[nom + "_sin_redundancia"] = contraste(M[fila_nac[cols_r]].sum(axis=0), sumas_fondo(M, filas[:, cols_r]))
    b15 = dict(zip(bloq.simbolo, bloq[COL_TANDEM]))
    multi = sorted(([nac_sim[k] for k in v] for v in miembros.values() if len(v) > 1), key=lambda x: orden_simbolo(x[0]))
    quitar_r = {s for x in multi for s in x[1:]}
    tam = np.array([len(v) for v in miembros.values()])
    redund = {
        "k": K_RED, "umbral_fraccion_compartida": UMBRAL_RED,
        "kmers_comunes_quitados": n_comunes, "kmers_comunes_definicion": f"presentes en más de {COMUN_RED} promotores NAC",
        "pares_sobre_umbral": int((Fr[np.triu_indices(len(nac_sim), 1)] > UMBRAL_RED).sum()),
        "grupos": len(miembros), "grupos_con_mas_de_1": len(multi), "genes_en_el_mayor": int(tam.max()),
        "quedan": int(len(cols_r)), "quitados": sorted(quitar_r, key=orden_simbolo),
        "quitados_tambien_por_tandem": len(quitar_r & quitar),
        "grupos_de_mas_de_1": [", ".join(x) + " [" + (", ".join(sorted({str(b15[s]) for s in x if s in b15})) or "sin bloque")
                               + "]" for x in multi],
        "genes_agrupados_fuera_de_bloques_en_tandem": sorted((s for x in multi for s in x if s not in b15), key=orden_simbolo)}
    print(f"redundancia ({K_RED}-mers, > {UMBRAL_RED:.0%} compartido; {n_comunes} k-mers comunes quitados): "
          f"{redund['pares_sobre_umbral']} pares, {len(miembros)} grupos ({len(multi)} con > 1, máx {int(tam.max())}); "
          f"quedan {len(cols_r)} NAC; de los {len(quitar_r)} quitados, {len(quitar_r & quitar)} también por tándem", flush=True)

    # Control de resolución: el principal con los primeros 1000 conjuntos (colección original)
    ctrl1000 = {}
    if control_1000:
        for nom in ("u_sitios", "c_sitios"):
            r0, r1 = res[nom], contraste(res[nom]["obs"], sumas[nom][:1000])
            s10, s1 = r0["q"] < 0.05, r1["q"] < 0.05
            mismo = np.sign(r0["obs"] - r0["media_fondo"]) == np.sign(r1["obs"] - r1["media_fondo"])
            ctrl1000[nom] = {"q<0.05_N10000": int(s10.sum()), "q<0.05_N1000": int(s1.sum()),
                             "en_ambos_mismo_signo": int((s10 & s1 & mismo).sum()),
                             "spearman_p": round(float(spearmanr(r0["p"], r1["p"], nan_policy="omit")[0]), 4),
                             "en_p_minima_N1000": int((r1["p"] <= 2 / 1001 * (1 + 1e-9)).sum())}
            print(f"control N=1000 {nom}: {ctrl1000[nom]}", flush=True)

    # Calibración: conjuntos de fondo como falso NAC; a) sin selección, b) con la selección del análisis real
    nM = len(um)
    lleno = np.array([(1 << len(reg_idx[j])) - 1 for j in um], np.int64)
    rng = np.random.default_rng(SEMILLA)
    bs = rng.choice(N, N_CALIB, replace=False)
    calib = {}
    for nom in ("u_sitios", "c_sitios"):
        S = sumas[nom]
        fr05, nq = [], []
        for b in bs:
            r = contraste(S[b], np.delete(S, b, axis=0))
            fr05.append(float((r["p"][r["evaluable"]] < 0.05).mean()))
            nq.append(int((r["q"] < 0.05).sum()))
        calib[nom + "_sin_seleccion"] = {"fraccion_p<0.05_media": float(np.mean(fr05)),
                                         "fraccion_p<0.05_rango": [min(fr05), max(fr05)], "n_q<0.05_por_conjunto": nq}
    t1 = time.time()
    sel_obs, sel_p, sel_pc = [], [], []
    nq_u, nq_uc, npr, fr_c, nq_c, n_camb = [], [], [], [], [], []
    for b in bs:
        visto = CR32[filas[b]].sum(axis=0) > 0
        sb = np.array([sum(1 << k for k, r in enumerate(reg_idx[j]) if visto[r]) for j in um], np.int64)
        uv = np.array([bool(visto[regs].any()) for regs in reg_idx])
        cl = np.nonzero((sb != lleno) & (sb != 0))[0]
        M = CT32.copy()
        if len(cl):
            es = np.zeros(nM, bool)
            es[cl] = True
            ok = es[EJL] & ((EB.astype(np.int64) & sb[EJL]) != 0)
            cnt = np.bincount(EF[ok].astype(np.int64) * nM + EJL[ok], minlength=nP * nM).reshape(nP, nM)
            for l in cl:
                M[:, um[l]] = cnt[:, l]
        M[:, ~uv] = 0
        Su = sumas["u_sitios"].copy()
        camb = [um[l] for l in cl]
        if camb:
            Su[:, camb] = sumas_fondo(M[:, camb], filas)
        So = np.delete(Su[:, uv], b, axis=0)
        r, rc_ = contraste(Su[b, uv], So), contraste_condicional(Su[b, uv], So)
        sel_obs.append(Su[b, uv]); sel_p.append(r["p"]); sel_pc.append(rc_["p"])
        nq_u.append(int((r["q"] < 0.05).sum())); nq_uc.append(int((rc_["q"] < 0.05).sum()))
        npr.append(int(uv.sum())); n_camb.append(len(camb))
        Sc = sumas_fondo(M @ G, filas)
        rcl = contraste(Sc[b], np.delete(Sc, b, axis=0))
        fr_c.append(float((rcl["p"][rcl["evaluable"]] < 0.05).mean())); nq_c.append(int((rcl["q"] < 0.05).sum()))
    so, sp, spc = np.concatenate(sel_obs), np.concatenate(sel_p), np.concatenate(sel_pc)
    est = np.array([estrato(o) for o in so])
    por_estrato = {}
    for a_, b_ in ESTRATOS:
        e = f"{a_}-{b_}" if b_ < 10 ** 9 else f">={a_}"
        m = est == e
        por_estrato[e] = {"pruebas": int(m.sum()), "fraccion_p<0.05": float((sp[m] < 0.05).mean()) if m.any() else None,
                          "fraccion_p_condicional<0.05": float((spc[m][np.isfinite(spc[m])] < 0.05).mean()) if m.any() else None,
                          "no_evaluables_condicional": int((~np.isfinite(spc[m])).sum())}
    calib["u_sitios_con_seleccion"] = {
        "unidades_probadas_media": float(np.mean(npr)), "unidades_con_registros_no_vistos_media": float(np.mean(n_camb)),
        "fraccion_p<0.05": float((sp < 0.05).mean()),
        "fraccion_p_condicional<0.05": float((spc[np.isfinite(spc)] < 0.05).mean()),
        "n_q<0.05_por_conjunto": nq_u, "n_q_condicional<0.05_por_conjunto": nq_uc, "por_estrato_de_sitios": por_estrato}
    calib["c_sitios_con_seleccion"] = {"fraccion_p<0.05_media": float(np.mean(fr_c)), "n_q<0.05_por_conjunto": nq_c}
    for k, v in calib.items():
        resumen_k = {kk: (vv if not isinstance(vv, list) or len(vv) <= 2 else f"máx {max(vv)}, media {np.mean(vv):.2f}")
                     for kk, vv in v.items()}
        print(f"calibración {k}: {resumen_k}", flush=True)
    print(f"calibración con selección en {time.time() - t1:.0f} s", flush=True)

    # 7. Tablas
    tu = t54[["unidad", "clase", "clase_en", "nombre_articulo1", "elementos_PlantCARE", "n_registros"]].copy()
    tu = tabla(tu, res, ["u_sitios", "u_presencia", "u_exactos", "u_sitios_sin_tandem", "u_sitios_sin_redundancia",
                         "u_sitios_condicional"],
               ["", "presencia_", "exactos_", "sin_tandem_", "sin_redundancia_", "condicional_"])
    tu = tu.drop(columns="condicional_obs").rename(columns={
        "obs": "sitios_NAC", "presencia_obs": "promotores_NAC_con_sitio", "exactos_obs": "exactos_sitios_NAC",
        "sin_tandem_obs": "sin_tandem_sitios_NAC", "sin_redundancia_obs": "sin_redundancia_sitios_NAC"})
    tu.insert(tu.columns.get_loc("sitios_NAC"), "promotores_NAC_con_sitio_pct",
              np.round(res["u_presencia"]["obs"] / N_NAC * 100, 1))
    tu["conjuntos_fondo_con_sitio_pct"] = np.round(res["u_sitios_condicional"]["conjuntos_con_sitio"] / N * 100, 2)
    with np.errstate(divide="ignore", invalid="ignore"):
        tu["razon_fondo_sin_emparejar"] = res["u_sitios"]["obs"] / (N_NAC * CT32[fila_fondo].mean(axis=0))
    tu["estrato_sitios_NAC"] = [estrato(o) for o in res["u_sitios"]["obs"]]
    tu["direccion"] = direccion(tu, "sitios_NAC")
    variantes_u = ["u_presencia", "u_exactos", "u_sitios_sin_tandem", "u_sitios_sin_redundancia", "u_sitios_condicional"]
    nombres_u = ["presencia", "exactos", "sin_tandem", "sin_redundancia", "condicional"]
    tu["variantes_concordantes"], tu["variantes_concordantes_lista"] = concordantes(res, "u_sitios", variantes_u, nombres_u)
    tu = tu.sort_values(["q", "p", "z"], ascending=[True, True, False])

    sit_nac = CT32[fila_nac].sum(axis=0)
    top = [", ".join(f"{unidades[j]} ({sit_nac[j]})" for j in sorted(np.nonzero(G[:, i])[0], key=lambda j: -sit_nac[j])[:3])
           for i in range(len(GRUPOS))]
    tc = pd.DataFrame({"clase": GRUPOS, "clase_en": [c2en[c] for c in GRUPOS], "clase_del_articulo1": [c in CLASES for c in GRUPOS],
                       "n_unidades": G.sum(axis=0), "unidades_con_mas_sitios": top})
    variantes_c = ["c_exactos", "c_sitios_sin_tandem", "c_sitios_sin_redundancia"]
    nombres_c = ["exactos", "sin_tandem", "sin_redundancia"]
    tc = tabla(tc, res, ["c_sitios"] + variantes_c, [""] + [x + "_" for x in nombres_c])
    tc = tc.rename(columns={"obs": "sitios_NAC", "exactos_obs": "exactos_sitios_NAC", "sin_tandem_obs": "sin_tandem_sitios_NAC",
                            "sin_redundancia_obs": "sin_redundancia_sitios_NAC"})
    tc.insert(tc.columns.get_loc("sitios_NAC"), "promotores_NAC_con_sitio_pct", np.round((CTc[fila_nac] > 0).mean(axis=0) * 100, 1))
    tc["razon_fondo_sin_emparejar"] = res["c_sitios"]["obs"] / (N_NAC * CTc[fila_fondo].mean(axis=0))
    tc["direccion"] = direccion(tc, "sitios_NAC")
    tc["variantes_concordantes"], tc["variantes_concordantes_lista"] = concordantes(res, "c_sitios", variantes_c, nombres_c)

    obs_c, esp_c, pct_c = por_gen(CTc, fila_nac, filas)
    obs_u, esp_u, pct_u = por_gen(CT32, fila_nac, filas)
    sub = dict(zip(t45.simbolo, t45.subfamilia))
    filas_g = []
    for j, s in sorted(enumerate(nac_sim), key=lambda x: orden_simbolo(x[1])):
        for i, c in enumerate(GRUPOS):
            filas_g.append((s, s2l[s], sub.get(s, ""), c, c2en[c], int(obs_c[j, i]), round(float(esp_c[j, i]), 3),
                            round(float(np.log2((obs_c[j, i] + 1) / (esp_c[j, i] + 1))), 3), round(float(pct_c[j, i]), 1)))
    tg = pd.DataFrame(filas_g, columns=["simbolo", "locus", "subfamilia", "clase", "clase_en", "sitios", "esperado_fondo_GC",
                                        "log2_obs_esp", "percentil_en_fondo"])
    assert len(tg) == N_NAC * len(GRUPOS)

    for t, n in ((tu, "T57_enriquecimiento_PlantCARE_elementos"), (tc, "T58_enriquecimiento_PlantCARE_clases"),
                 (tg, "T59_clases_PlantCARE_por_gen_NAC189")):
        t.to_csv(os.path.join(TAB, f"{n}{suf}.csv"), index=False, float_format="%.6g")
    np.savez_compressed(os.path.join(CP, f"fondo_sumas_PlantCARE{suf}.npz"), **{k: v.astype(np.int32) for k, v in sumas.items()},
                        unidades=np.array(unidades), grupos=np.array(GRUPOS))
    np.savez_compressed(os.path.join(CP, f"esperado_por_gen_PlantCARE{suf}.npz"), genes=np.array(nac_sim),
                        unidades=np.array(unidades), observado=obs_u.astype(np.int16), esperado=esp_u.astype(np.float32),
                        percentil=pct_u.astype(np.float32))

    # 8. Concordancia, resolución y resumen
    conc = {}
    for base, otras in (("u_sitios", variantes_u), ("c_sitios", variantes_c)):
        r0 = res[base]
        sig = r0["q"] < 0.05
        for o in otras:
            r1 = res[o]
            ok = np.isfinite(r0["log2_razon"]) & np.isfinite(r1["log2_razon"])
            misma = sig & (r1["q"] < 0.05) & (np.sign(r1["obs"] - r1["media_fondo"]) == np.sign(r0["obs"] - r0["media_fondo"]))
            conc[o] = {"spearman_log2_razon": round(float(spearmanr(r0["log2_razon"][ok], r1["log2_razon"][ok])[0]), 3),
                       "sig_principal": int(sig.sum()), "sig_tambien_aqui_mismo_signo": int(misma.sum()),
                       "sig_aqui": int((r1["q"] < 0.05).sum())}
            print(f"sensibilidad {o}: rho {conc[o]['spearman_log2_razon']}; de {int(sig.sum())} significativos en el "
                  f"principal, {int(misma.sum())} siguen con q<0,05 y el mismo signo; aquí {int((r1['q'] < 0.05).sum())}", flush=True)
    piso = 2 / (N + 1)
    resol = {k: {"evaluables": int(r["evaluable"].sum()), "en_p_minima": int((r["p"] <= piso * (1 + 1e-9)).sum()),
                 "evaluables_sin_sitios_en_fondo": int((r["evaluable"] & ~(r["media_fondo"] > 0)).sum()),
                 "q_minima": float(np.nanmin(r["q"])), "q<0.05": int((r["q"] < 0.05).sum())} for k, r in res.items()}
    sig_est = tu[tu.q < 0.05].groupby("estrato_sitios_NAC").size().to_dict()
    resumen = {
        "fecha": time.strftime("%Y-%m-%d %H:%M:%S"), "archivo_fondo": os.path.relpath(FONDO, BASE),
        "n_conjuntos_fondo": int(N), "n_unidades": nU, "n_registros": nR, "pares_unidad_registro": n_pares,
        "regla": "puntaje >= 0,9 del máximo del registro (fija 1, degenerada 0,5), ambas hebras; una posición por gen y "
                 "unidad (funciones de e3_16)",
        "escaneo": {"archivo": os.path.relpath(ESC, BASE), "reutilizado": seg_escaneo is None, "promotores": nP,
                    "trozos": N_TROZOS, "procesos": procesos, "segundos": seg_escaneo, "sitios_genoma": int(CT32.sum()),
                    "exactos_genoma": int(CE32.sum()), "sitios_con_marcas_de_registro": int(len(EF))},
        "controles": controles, "p_minima_posible_bilateral": piso, "resolucion": resol, "control_N1000": ctrl1000,
        "n_calibracion": N_CALIB, "semilla": SEMILLA, "calibracion": calib, "concordancia": conc,
        "tandem_columna": COL_TANDEM, "tandem_quitados": sorted(quitar), "redundancia": redund,
        "variantes_unidades": nombres_u, "variantes_clases": nombres_c,
        "unidades_por_estrato_de_sitios_NAC": tu.estrato_sitios_NAC.value_counts().to_dict(),
        "unidades_q<0.05_por_estrato": {k: int(v) for k, v in sig_est.items()},
        "unidades_q<0.05": {"enriquecidas": int((tu.direccion == "enriquecido").sum()),
                            "empobrecidas": int((tu.direccion == "empobrecido").sum())},
        "unidades_q<0.05_por_variantes_concordantes": {
            str(k): int(v) for k, v in tu.variantes_concordantes.dropna().astype(int).value_counts().sort_index().items()},
        "unidades_concordantes_en_todas_las_variantes": [
            f"{u} ({c}; razón {r:.2f}, q {q:.2g})" for u, c, r, q in
            tu.loc[tu.variantes_concordantes == len(variantes_u), ["unidad", "clase", "razon", "q"]].itertuples(index=False)],
        "clases": {c: {"sitios_NAC": int(o), "media_fondo": round(float(m), 1), "razon": round(float(r), 3),
                       "p": float(p), "q": float(q), "direccion": d,
                       "variantes_concordantes": (None if pd.isna(v) else int(v))}
                   for c, o, m, r, p, q, d, v in tc[["clase", "sitios_NAC", "media_fondo", "razon", "p", "q", "direccion",
                                                     "variantes_concordantes"]].itertuples(index=False)},
        "limitacion": "solo se prueban los registros que PlantCARE informó en los 189 (la base no se puede descargar)",
        "versiones": {"python": platform.python_version(), "numpy": np.__version__, "pandas": pd.__version__,
                      "scipy": scipy.__version__},
        "segundos": round(time.time() - t0)}
    json.dump(resumen, open(os.path.join(CP, f"resumen_contraste_PlantCARE{suf}.json"), "w", encoding="utf-8"),
              ensure_ascii=False, indent=1)
    for k, v in resol.items():
        print(f"resolución {k}: {v}")
    print(f"unidades por estrato de sitios NAC: {resumen['unidades_por_estrato_de_sitios_NAC']}; q<0,05 por estrato: "
          f"{resumen['unidades_q<0.05_por_estrato']}")
    print(f"unidades q<0,05: {resumen['unidades_q<0.05']}; por variantes concordantes de {len(variantes_u)}: "
          f"{resumen['unidades_q<0.05_por_variantes_concordantes']}")
    pd.set_option("display.width", 250)
    print(tc[["clase", "n_unidades", "sitios_NAC", "media_fondo", "razon", "z", "p", "q", "exactos_q", "sin_tandem_q",
              "sin_redundancia_q", "razon_fondo_sin_emparejar", "direccion"]].to_string(
        index=False, float_format=lambda v: f"{v:.3g}"))
    print(tu[["unidad", "clase", "sitios_NAC", "promotores_NAC_con_sitio_pct", "media_fondo", "razon", "z", "p", "q",
              "sin_redundancia_razon", "sin_redundancia_q", "condicional_q", "variantes_concordantes",
              "razon_fondo_sin_emparejar"]].head(30).to_string(index=False, float_format=lambda v: f"{v:.3g}"))
    print(f"listo en {time.time() - t0:.0f} s")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=10000, help="conjuntos de fondo (10000 principal; 1000 = colección original)")
    ap.add_argument("--procesos", type=int, default=None, help="procesos para el escaneo (por defecto, núcleos - 2, máx 10)")
    ap.add_argument("--reusar_escaneo", action="store_true", help="usar escaneo_PlantCARE_genoma.npz si coincide con las entradas")
    a = ap.parse_args()
    main(a.n, a.procesos, a.reusar_escaneo)
