"""Entrega 3, análisis cis local (JASPAR 2026 + FIMO), paso 2: contraste con el fondo emparejado por GC.

Los 189 promotores NAC se comparan con N conjuntos de 189 promotores de fondo de e3_13 (cada promotor de
fondo cae en el mismo intervalo de GC que su NAC; fondo = 35.754 promotores no NAC). Análisis principal:
N = 10.000 (fondo_GC_conjuntos_10000.npz, cuyos primeros 1000 conjuntos son la colección original). Con
1000 conjuntos la p empírica bilateral no baja de 2/1001 y, con cerca de 900 pruebas, BH no resuelve:
todos los motivos significativos quedaban en ese mínimo.

Medidas por promotor (sitios FIMO con p < 1e-4, ancho del motivo, cualquier hebra):
  motivo:   número de sitios del motivo (principal) y presencia (>= 1 sitio, secundaria);
  familia:  número de regiones de la familia = sitios de cualquier motivo de la familia, en cualquier
            hebra, fusionados cuando se solapan (los motivos redundantes de una familia cuentan una vez).
Estadístico de grupo = suma en los 189 promotores. Para cada conjunto de fondo se calcula la misma suma.
  razón = observado / media del fondo; z = (observado - media) / DE del fondo;
  p empírica bilateral = min(1, 2 * min(P_enr, P_emp)), con P_enr = (1 + #{fondo >= obs}) / (N + 1) y
  P_emp = (1 + #{fondo <= obs}) / (N + 1); q de Benjamini-Hochberg dentro de cada nivel (motivos, familias)
  y de cada análisis, solo con las pruebas evaluables. No evaluable = ningún sitio ni en los NAC ni en
  ningún conjunto de fondo (p = 1 por construcción; son motivos cortos que no llegan al umbral de FIMO):
  p y q quedan vacías y la dirección es «no evaluable».
Sensibilidad: umbral FIMO p < 1e-5; presencia por motivo; un solo gen por bloque en tándem (T45.bloque_T15 =
  T15_bloques_tandem.csv, la definición del proyecto: 21 bloques, 119 genes, por adyacencia en el orden génico;
  se conserva el primero por posición y se quitan las columnas de fondo de los demás).
Calibración: N_CALIB conjuntos de fondo hacen de falso «NAC» frente a los N - 1 restantes; se espera ~5 %
  de p < 0,05 entre las pruebas evaluables y casi ningún q < 0,05.
Por gen: para cada NAC y familia, regiones observadas frente a la media de sus N promotores de fondo
  emparejados (columna j de idx): log2((obs + 1) / (esperado + 1)) y percentil en su fondo.

Uso: python e3_15_contraste_jaspar.py [--n N]   (N = 10000 por defecto, 1000 = colección original;
     desde cualquier carpeta; python de nac3)
Entradas: datos/entrega3/cis_jaspar/{conteos_JASPAR2026.npz, sitios/trozo_KK.npz, motivos_JASPAR2026.tsv,
          loci_orden.txt}; datos/entrega3/cis_plantcare/fondo/fondo_GC_conjuntos[_N].npz;
          tablas/T40_promotores_NAC_2kb.csv, T45_modos_duplicacion_NAC189.csv
Salidas: tablas/T51_enriquecimiento_JASPAR_motivos.csv, T52_enriquecimiento_JASPAR_familias.csv,
         T53_familias_JASPAR_por_gen_NAC189.csv; datos/entrega3/cis_jaspar/{regiones_familia_JASPAR2026.npz,
         fondo_sumas_JASPAR2026.npz, resumen_contraste_JASPAR2026.json}
"""
import argparse, json, os, time
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

BASE = (os.environ.get("NAC_ROOT", ".") + "")
A = os.path.join(BASE, "ANALISIS_Eglobulus")
E3 = os.path.join(A, "datos", "entrega3")
CJ = os.path.join(E3, "cis_jaspar")
TAB = os.path.join(A, "tablas")
N_TROZOS = 64
UMBRAL_ESTRICTO = 1e-5
N_CALIB = 100
SEMILLA = 20261018
DESP = 4096  # > largo del promotor (2000): separa los grupos en la fusión de intervalos

# Familias armonizadas (sinónimos de JASPAR y familias tomadas de la clase cuando falta la familia)
ARMON = {"Group S": "bZIP grupo S", "S": "bZIP grupo S", "Group A": "bZIP grupo A", "group A": "bZIP grupo A",
         "Group D": "bZIP grupo D", "D": "bZIP grupo D", "Group G": "bZIP grupo G", "Group H": "bZIP grupo H",
         "Group C": "bZIP grupo C", "Group I": "bZIP grupo I", "Group B": "bZIP grupo B", "Group K": "bZIP grupo K",
         "bZIP": "bZIP sin grupo", "Basic leucine zipper factors (bZIP)": "bZIP sin grupo",
         "Basic helix-loop-helix factors (bHLH)": "bHLH", "C2H2 zinc finger factors": "C2H2",
         "HC3 zinc ribbon factors": "ALOG", "Type II": "MADS (MIKC)", "MIKC": "MADS (MIKC)",
         "TALE-type HD": "TALE", "PLINC": "ZF-HD", "AP2/EREBP": "AP2/EREBP sin subfamilia", "": "sin anotación"}


def bh(p):
    p = np.asarray(p, float)
    n = len(p)
    o = np.argsort(p)
    q = p[o] * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n)
    out[o] = np.minimum(q, 1)
    return out


def sumas_fondo(X, filas, bloque=50):
    """X (promotores x K), filas (conjuntos x columnas) -> sumas (conjuntos x K)."""
    S = np.empty((filas.shape[0], X.shape[1]), dtype=np.int64)
    for i in range(0, filas.shape[0], bloque):
        S[i:i + bloque] = X[filas[i:i + bloque]].sum(axis=1, dtype=np.int64)
    return S


def contraste(obs, S):
    mu, de = S.mean(0), S.std(0, ddof=1)
    n = S.shape[0]
    p_enr = (1 + (S >= obs).sum(0)) / (n + 1)
    p_emp = (1 + (S <= obs).sum(0)) / (n + 1)
    # Evaluable: algún sitio en el grupo o en algún conjunto de fondo; si todo es cero, p = 1 por construcción
    ev = (obs > 0) | (S.max(0) > 0)
    p = np.where(ev, np.minimum(1, 2 * np.minimum(p_enr, p_emp)), np.nan)
    q = np.full(len(p), np.nan)
    q[ev] = bh(p[ev])
    with np.errstate(divide="ignore", invalid="ignore"):
        razon = np.where(mu > 0, obs / mu, np.nan)
        z = np.where(de > 0, (obs - mu) / de, np.nan)
        l2 = np.log2(razon)
    return {"obs": obs, "media_fondo": mu, "de_fondo": de, "razon": razon, "log2_razon": l2,
            "z": z, "p_enr": p_enr, "p_emp": p_emp, "p": p, "q": q, "evaluable": ev}


def concordantes(res, base, otras, nombres):
    """Para las pruebas con q < 0,05 en el análisis principal: cuántas variantes de sensibilidad (y cuáles) también
    tienen q < 0,05 con el mismo signo (observado frente a la media del fondo). Vacío si no es significativa."""
    r0 = res[base]
    sig = r0["q"] < 0.05
    s0 = np.sign(r0["obs"] - r0["media_fondo"])
    m = np.array([(res[o]["q"] < 0.05) & (np.sign(res[o]["obs"] - res[o]["media_fondo"]) == s0) for o in otras])
    n = np.where(sig, m.sum(0), np.nan)
    lista = [";".join(nm for nm, x in zip(nombres, m[:, i]) if x) if sig[i] else "" for i in range(len(sig))]
    return n, lista


def direccion(t, col_obs):
    # signo por observado frente a la media del fondo (la razón no existe si la media del fondo es 0)
    d = np.where(t[col_obs] > t.media_fondo, "enriquecido", "empobrecido")
    return np.where(t.p.isna(), "no evaluable", np.where(t.q < 0.05, d, "ns"))


def main(n_conj=10000):
    t0 = time.time()
    FONDO = os.path.join(E3, "cis_plantcare", "fondo",
                         "fondo_GC_conjuntos.npz" if n_conj == 1000 else f"fondo_GC_conjuntos_{n_conj}.npz")
    mot = pd.read_csv(os.path.join(CJ, "motivos_JASPAR2026.tsv"), sep="\t", keep_default_na=False)
    loci = open(os.path.join(CJ, "loci_orden.txt")).read().split()
    l2i = {x: i for i, x in enumerate(loci)}
    cz = np.load(os.path.join(CJ, "conteos_JASPAR2026.npz"))
    C4, C5 = cz["p1e4"].astype(np.int32), cz["p1e5"].astype(np.int32)
    assert list(cz["loci"]) == loci and list(cz["motivos"]) == mot.motivo_id.tolist()
    nS, nM = C4.shape

    # Familias armonizadas
    mot["familia_armonizada"] = mot.familia.map(lambda f: ARMON.get(f, f))
    fams = sorted(mot.familia_armonizada.unique())
    f2i = {f: i for i, f in enumerate(fams)}
    m2f = mot.familia_armonizada.map(f2i).to_numpy()
    nF = len(fams)
    ancho = mot.ancho.to_numpy()

    # Fondo emparejado por GC y NAC
    fz = np.load(FONDO)
    nac_sim, loc_fondo, idx = list(fz["nac"]), list(fz["locus_fondo"]), fz["idx"]
    t40 = pd.read_csv(os.path.join(TAB, "T40_promotores_NAC_2kb.csv"))
    s2l = dict(zip(t40.simbolo, t40.locus))
    assert set(nac_sim) == set(s2l) and len(nac_sim) == 189
    fila_nac = np.array([l2i[s2l[s]] for s in nac_sim])
    fila_fondo = np.array([l2i[x] for x in loc_fondo])
    assert len(set(fila_fondo) & set(fila_nac)) == 0 and len(fila_fondo) + len(fila_nac) == nS
    filas = fila_fondo[idx]                      # N x 189, filas globales
    N = filas.shape[0]
    assert N == n_conj, f"{FONDO}: {N} conjuntos, se esperaban {n_conj}"
    print(f"promotores {nS}; motivos {nM}; familias armonizadas {nF}; conjuntos de fondo {N} x {filas.shape[1]}")

    # Regiones por familia (fusión de solapes, ambas hebras), p<1e-4 y p<1e-5
    F4 = np.zeros(nS * nF, dtype=np.int64)
    F5 = np.zeros(nS * nF, dtype=np.int64)
    S4 = np.zeros(nS * nF, dtype=np.int64)
    S5 = np.zeros(nS * nF, dtype=np.int64)
    tot = 0
    for k in range(1, N_TROZOS + 1):
        z = np.load(os.path.join(CJ, "sitios", f"trozo_{k:02d}.npz"))
        s, m, ini, p = z["seq"].astype(np.int64), z["motivo"].astype(np.int64), z["inicio"].astype(np.int64), z["p"]
        tot += len(s)
        f = m2f[m]
        fin = ini + ancho[m] - 1
        S4 += np.bincount(s * nF + f, minlength=nS * nF)
        S5 += np.bincount((s * nF + f)[p <= UMBRAL_ESTRICTO], minlength=nS * nF)   # mismo criterio que e3_14
        for F, sel in ((F4, slice(None)), (F5, p <= UMBRAL_ESTRICTO)):
            g, a, b = (s * nF + f)[sel], ini[sel], fin[sel]
            o = np.lexsort((a, g))
            g, a, b = g[o], a[o], b[o]
            cm = np.maximum.accumulate(g * DESP + b)
            nuevo = np.ones(len(g), bool)
            nuevo[1:] = (g[1:] != g[:-1]) | (g[1:] * DESP + a[1:] > cm[:-1])
            F += np.bincount(g[nuevo], minlength=nS * nF)
    F4, F5, S4, S5 = (x.reshape(nS, nF) for x in (F4, F5, S4, S5))
    assert tot == C4.sum(), f"sitios en trozos {tot} != matriz {C4.sum()}"
    # Sitios por promotor y familia leídos de los trozos = matriz de recuentos de e3_14 agregada por familia
    C4f = np.stack([C4[:, m2f == i].sum(1) for i in range(nF)], axis=1)
    C5f = np.stack([C5[:, m2f == i].sum(1) for i in range(nF)], axis=1)
    assert (S4 == C4f).all() and (S5 == C5f).all(), "sitios por promotor y familia (trozos) != matriz de recuentos"
    assert (F4 <= S4).all() and ((F4 > 0) == (S4 > 0)).all()
    assert (F5 <= S5).all() and ((F5 > 0) == (S5 > 0)).all()
    # F5 puede superar a F4: si un sitio con 1e-5 < p < 1e-4 unía dos sitios estrictos en una sola región,
    # al quitarlo la región se parte en dos. No es un error; se cuenta.
    n_partidas = int((F5 > F4).sum())
    print(f"control: {tot:,} sitios = suma de la matriz; sitios por promotor y familia = matriz agregada (p<1e-4 y "
          f"p<=1e-5); regiones <= sitios; fusión: {S4.sum():,} sitios -> {F4.sum():,} regiones (p<1e-4), "
          f"{S5.sum():,} -> {F5.sum():,} (p<=1e-5); celdas con más regiones estrictas que a p<1e-4: {n_partidas:,}")
    np.savez_compressed(os.path.join(CJ, "regiones_familia_JASPAR2026.npz"), p1e4=F4.astype(np.int16),
                        p1e5=F5.astype(np.int16), loci=np.array(loci), familias=np.array(fams))

    # Contrastes: principal y sensibilidad
    P4 = (C4 > 0).astype(np.int32)
    medidas = {"mot_sitios": C4, "mot_presencia": P4, "mot_sitios_p1e5": C5, "fam_regiones": F4, "fam_regiones_p1e5": F5}
    res, sumas = {}, {}
    for nom, X in medidas.items():
        sumas[nom] = sumas_fondo(X, filas)
        res[nom] = contraste(X[fila_nac].sum(0), sumas[nom])

    # Un gen por bloque en tándem
    t45 = pd.read_csv(os.path.join(TAB, "T45_modos_duplicacion_NAC189.csv"))
    t45 = t45.sort_values(["cromosoma", "inicio"])
    COL_TANDEM = "bloque_T15"          # definición autoritativa (T15); bloque_tandem_T01 es la de T07 (distancia física)
    bloq = t45.dropna(subset=[COL_TANDEM])
    quitar = set(bloq.simbolo) - set(bloq.groupby(COL_TANDEM).simbolo.first())
    cols = np.array([j for j, s in enumerate(nac_sim) if s not in quitar])
    for nom in ("mot_sitios", "fam_regiones"):
        X = medidas[nom]
        res[nom + "_sin_tandem"] = contraste(X[fila_nac[cols]].sum(0), sumas_fondo(X, filas[:, cols]))
    print(f"tándem ({COL_TANDEM}): {bloq[COL_TANDEM].nunique()} bloques con {len(bloq)} genes; se quitan {len(quitar)}; "
          f"quedan {len(cols)} NAC")

    # Calibración: conjuntos de fondo como falso NAC
    rng = np.random.default_rng(SEMILLA)
    calib = {}
    for nom in ("mot_sitios", "fam_regiones"):
        S = sumas[nom]
        fr05, nq = [], []
        for b in rng.choice(N, N_CALIB, replace=False):
            r = contraste(S[b], np.delete(S, b, axis=0))
            fr05.append(float((r["p"][r["evaluable"]] < 0.05).mean()))
            nq.append(int((r["q"] < 0.05).sum()))
        calib[nom] = {"fraccion_p<0.05_media": float(np.mean(fr05)), "fraccion_p<0.05_rango": [min(fr05), max(fr05)],
                      "n_q<0.05_por_conjunto": nq}
        print(f"calibración {nom}: p<0,05 en {np.mean(fr05) * 100:.1f} % de las pruebas (rango "
              f"{min(fr05) * 100:.1f}-{max(fr05) * 100:.1f} %); q<0,05 por conjunto: máx {max(nq)}, media {np.mean(nq):.2f}")

    # Tablas
    def tabla(base, claves, prefijos):
        t = base.copy()
        for clave, pre in zip(claves, prefijos):
            r = res[clave]
            cols_r = ["obs", "media_fondo", "de_fondo", "razon", "log2_razon", "z", "p", "q"] if pre == "" else ["obs", "razon", "p", "q"]
            for c in cols_r:
                t[(pre + c) if pre else c] = r[c]
        return t

    tm = mot[["motivo_id", "nombre", "familia_armonizada", "familia_JASPAR", "clase_JASPAR", "ancho", "tipo_dato", "uniprot"]].copy()
    tm = tabla(tm, ["mot_sitios", "mot_presencia", "mot_sitios_p1e5", "mot_sitios_sin_tandem"],
               ["", "presencia_", "p1e5_", "sin_tandem_"])
    tm = tm.rename(columns={"obs": "sitios_NAC", "presencia_obs": "promotores_NAC_con_sitio", "p1e5_obs": "p1e5_sitios_NAC",
                            "sin_tandem_obs": "sin_tandem_sitios_NAC"})
    tm.insert(tm.columns.get_loc("sitios_NAC"), "promotores_NAC_con_sitio_pct", np.round(res["mot_presencia"]["obs"] / 189 * 100, 1))
    tm["direccion"] = direccion(tm, "sitios_NAC")
    tm["variantes_concordantes"], tm["variantes_concordantes_lista"] = concordantes(
        res, "mot_sitios", ["mot_presencia", "mot_sitios_p1e5", "mot_sitios_sin_tandem"], ["presencia", "p1e5", "sin_tandem"])
    tm = tm.sort_values(["q", "p", "z"], ascending=[True, True, False])

    nmot = mot.groupby("familia_armonizada").motivo_id.size()
    ej = mot.groupby("familia_armonizada").nombre.apply(lambda s: ", ".join(sorted(s)[:6]))
    tf = pd.DataFrame({"familia": fams, "n_motivos": nmot.reindex(fams).to_numpy(), "ejemplos": ej.reindex(fams).to_numpy()})
    tf = tabla(tf, ["fam_regiones", "fam_regiones_p1e5", "fam_regiones_sin_tandem"], ["", "p1e5_", "sin_tandem_"])
    tf = tf.rename(columns={"obs": "regiones_NAC", "p1e5_obs": "p1e5_regiones_NAC", "sin_tandem_obs": "sin_tandem_regiones_NAC"})
    tf.insert(3, "promotores_NAC_con_region_pct", np.round((F4[fila_nac] > 0).mean(0) * 100, 1))
    tf["direccion"] = direccion(tf, "regiones_NAC")
    tf["variantes_concordantes"], tf["variantes_concordantes_lista"] = concordantes(
        res, "fam_regiones", ["fam_regiones_p1e5", "fam_regiones_sin_tandem"], ["p1e5", "sin_tandem"])
    tf = tf.sort_values(["q", "p", "z"], ascending=[True, True, False])

    # Por gen: observado frente a sus N promotores de fondo emparejados (por bloques de 1000 conjuntos)
    obs_g = F4[fila_nac]                                  # 189 x nF
    suma_g, men_g, igu_g = (np.zeros(obs_g.shape) for _ in range(3))
    for i in range(0, N, 1000):
        fg = F4[filas[i:i + 1000]]                        # b x 189 x nF
        suma_g += fg.sum(0); men_g += (fg < obs_g[None]).sum(0); igu_g += (fg == obs_g[None]).sum(0)
    esp_g = suma_g / N
    pct_g = (men_g + 0.5 * igu_g) / N * 100
    sub = dict(zip(t45.simbolo, t45.subfamilia))
    filas_g = []
    for j, s in enumerate(nac_sim):
        for i, fa in enumerate(fams):
            filas_g.append((s, s2l[s], sub.get(s, ""), fa, int(obs_g[j, i]), round(float(esp_g[j, i]), 3),
                            round(float(np.log2((obs_g[j, i] + 1) / (esp_g[j, i] + 1))), 3), round(float(pct_g[j, i]), 1)))
    tg = pd.DataFrame(filas_g, columns=["simbolo", "locus", "subfamilia", "familia", "regiones", "esperado_fondo_GC",
                                        "log2_obs_esp", "percentil_en_fondo"])

    for t, n in ((tm, "T51_enriquecimiento_JASPAR_motivos.csv"), (tf, "T52_enriquecimiento_JASPAR_familias.csv"),
                 (tg, "T53_familias_JASPAR_por_gen_NAC189.csv")):
        t.to_csv(os.path.join(TAB, n), index=False, float_format="%.6g")
    np.savez_compressed(os.path.join(CJ, "fondo_sumas_JASPAR2026.npz"), **{k: v.astype(np.int32) for k, v in sumas.items()},
                        motivos=mot.motivo_id.to_numpy().astype(str), familias=np.array(fams))

    # Concordancia de las sensibilidades con el análisis principal
    conc = {}
    for nivel, base, otras in (("motivos", "mot_sitios", ["mot_presencia", "mot_sitios_p1e5", "mot_sitios_sin_tandem"]),
                               ("familias", "fam_regiones", ["fam_regiones_p1e5", "fam_regiones_sin_tandem"])):
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
                  f"principal, {int(misma.sum())} siguen con q<0,05 y el mismo signo; aquí {int((r1['q'] < 0.05).sum())}")

    piso = 2 / (N + 1)
    resol = {k: {"evaluables": int(r["evaluable"].sum()),
                 "en_p_minima": int((r["p"] <= piso * (1 + 1e-9)).sum()),
                 "evaluables_sin_sitios_en_fondo": int((r["evaluable"] & (r["media_fondo"] == 0)).sum()),
                 "q_minima": float(np.nanmin(r["q"])), "q<0.05": int((r["q"] < 0.05).sum())} for k, r in res.items()}
    for k, v in resol.items():
        print(f"resolución {k}: {v}")
    resumen = {"fecha": time.strftime("%Y-%m-%d %H:%M:%S"), "archivo_fondo": os.path.relpath(FONDO, BASE),
               "n_conjuntos_fondo": int(N), "n_motivos": int(nM), "n_familias": int(nF), "familias_armonizadas": ARMON,
               "p_minima_posible_bilateral": piso, "p_minima_posible_unilateral": 1 / (N + 1), "resolucion": resol,
               "n_calibracion": N_CALIB,
               "motivos_q<0.05": {"enriquecidos": int((tm.direccion == "enriquecido").sum()),
                                  "empobrecidos": int((tm.direccion == "empobrecido").sum())},
               "motivos_no_evaluables": int((tm.direccion == "no evaluable").sum()),
               "motivos_q<0.05_por_variantes_concordantes_de_3": {str(k): int(v) for k, v in
                                                                   tm.variantes_concordantes.dropna().astype(int).value_counts().sort_index().items()},
               "motivos_concordantes_en_las_3_variantes": [f"{a} {b} ({c})" for a, b, c in
                                                           tm.loc[tm.variantes_concordantes == 3, ["motivo_id", "nombre", "familia_armonizada"]].itertuples(index=False)],
               "familias_q<0.05": {"enriquecidas": int((tf.direccion == "enriquecido").sum()),
                                   "empobrecidas": int((tf.direccion == "empobrecido").sum())},
               "familias_no_evaluables": int((tf.direccion == "no evaluable").sum()),
               "familias_q<0.05_por_variantes_concordantes_de_2": {str(k): int(v) for k, v in
                                                                    tf.variantes_concordantes.dropna().astype(int).value_counts().sort_index().items()},
               "tandem_columna": COL_TANDEM, "tandem_quitados": sorted(quitar), "calibracion": calib, "concordancia": conc,
               "umbral_estricto": "p <= 1e-5 (float32, mismo criterio que la matriz p1e5 de e3_14)",
               "regiones_familia": {"sitios_p1e4": int(S4.sum()), "regiones_p1e4": int(F4.sum()),
                                    "sitios_p1e5": int(S5.sum()), "regiones_p1e5": int(F5.sum()),
                                    "celdas_promotor_familia_con_F5_mayor_que_F4": n_partidas},
               "segundos": round(time.time() - t0)}
    json.dump(resumen, open(os.path.join(CJ, "resumen_contraste_JASPAR2026.json"), "w"), ensure_ascii=False, indent=1)
    print(f"motivos q<0,05: {resumen['motivos_q<0.05']} (no evaluables {resumen['motivos_no_evaluables']}); por variantes "
          f"concordantes de 3: {resumen['motivos_q<0.05_por_variantes_concordantes_de_3']}; en las 3: "
          f"{resumen['motivos_concordantes_en_las_3_variantes']}")
    print(f"familias q<0,05: {resumen['familias_q<0.05']} (no evaluables {resumen['familias_no_evaluables']}); por variantes "
          f"concordantes de 2: {resumen['familias_q<0.05_por_variantes_concordantes_de_2']}")
    pd.set_option("display.width", 220)
    print(tf[["familia", "n_motivos", "regiones_NAC", "media_fondo", "razon", "z", "p", "q", "sin_tandem_q", "p1e5_q"]]
          .head(25).to_string(index=False, float_format=lambda v: f"{v:.3g}"))
    print(f"listo en {time.time() - t0:.0f} s")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=10000, help="conjuntos de fondo (10000 principal; 1000 = colección original)")
    main(ap.parse_args().n)
