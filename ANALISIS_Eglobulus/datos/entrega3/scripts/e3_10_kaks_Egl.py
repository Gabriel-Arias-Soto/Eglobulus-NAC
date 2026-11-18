"""Entrega 3, paso 4.2 (artículo 3): Ka/Ks y datación de los pares NAC duplicados de E. globulus.

Entradas
  tablas/T46_pares_duplicados_NAC.csv                    pares NAC-NAC (e3_09_colinealidad_Egl.py)
  datos/entrega3/colinealidad/Egl_pares_colineales.tsv   pares ancla de MCScanX en todo el genoma
  datos/entrega3/isoformas/Egl.pep, Egl.cds              isoforma representativa por gen

Pares analizados
  Los 615 pares NAC-NAC de T46 (columna tipo_par): WGD/segmentales, tándem y proximales de MCScanX;
  pares a <= 200 kb (regla del artículo 3) con acierto DIAMOND pero sin modo MCScanX; y pares a
  <= 200 kb sin acierto DIAMOND (ninguna de las 5 dianas por consulta con e <= 1e-10; se alinean igual
  porque los dos genes tienen el dominio NAM). Contraste: todos los pares
  ancla del genoma.

Método
  1. Conciliación CDS-proteína: se busca el marco (0, 1 o 2) cuya traducción contiene la proteína;
     sin coincidencia exacta se acepta el marco con <= 2 % de diferencias; si no, el gen queda fuera.
  2. Alineamiento de proteínas por par con MAFFT L-INS-i (--localpair --maxiterate 1000) y
     retrotraducción a codones (como PAL2NAL); se quitan columnas con hueco, codones ambiguos o de parada.
  3. Ka y Ks con yn00 de PAML (Yang y Nielsen 2000; estimación principal, columnas *_YN00) y con
     NG86 y YN00 de Biopython (Bio.Align.analysis.calculate_dn_ds; columnas *_NG86 y *_YN00_Biopython)
     sobre el mismo alineamiento de codones. El YN00 de Biopython 1.88 devuelve NaN en muchos pares
     (división por cero en la ponderación de caminos), por eso PAML es la estimación principal.
  4. Prueba de razón de verosimilitudes por par NAC con codeml de PAML (runmode -2, F3x4): omega libre
     frente a omega = 1 (chi2, 1 g.l.), FDR de Benjamini-Hochberg sobre los pares NAC.
  5. T = Ks / (2 lambda), lambda = 6,5e-9 sustituciones/sitio/año (artículo 1). Ks > 2 se marca como
     saturado y no se data.
  6. Picos de la distribución de Ks de los pares ancla por densidad de núcleo gaussiana.

Salidas
  tablas/T47_KaKs_pares_NAC.csv
  datos/entrega3/colinealidad/Ks_anclas_genoma.tsv
  datos/entrega3/colinealidad/alineamientos_codones_NAC.fa
  datos/entrega3/colinealidad/resumen_kaks_Egl.json

Uso: python e3_10_kaks_Egl.py [--prueba N] [--cache]
  --prueba: N pares NAC con acierto DIAMOND, N sin acierto y N anclas, sin escribir salidas; --cache: reusa los alineamientos y estimaciones guardados en colinealidad/cache_e3_10_kaks.pkl
"""
import argparse
import json
import os
import re
import subprocess
import tempfile
import time
import warnings
from multiprocessing import Pool

import numpy as np
import pandas as pd
from Bio.Align import Alignment
from Bio.Align.analysis import calculate_dn_ds
from Bio.Seq import Seq
from scipy.stats import chi2, gaussian_kde, mannwhitneyu, spearmanr
from statsmodels.stats.multitest import multipletests

warnings.filterwarnings("ignore", category=RuntimeWarning)

RAIZ = (os.environ.get("NAC_ROOT", ".") + "/ANALISIS_Eglobulus")
ISO = f"{RAIZ}/datos/entrega3/isoformas"
COL = f"{RAIZ}/datos/entrega3/colinealidad"
BIN = (os.environ.get("CONDA_ENVS", os.path.expanduser("~/conda/envs")) + "/nac3/bin")
T46 = f"{RAIZ}/tablas/T46_pares_duplicados_NAC.csv"
T47 = f"{RAIZ}/tablas/T47_KaKs_pares_NAC.csv"
ANCLAS = f"{COL}/Egl_pares_colineales.tsv"
KS_GEN = f"{COL}/Ks_anclas_genoma.tsv"
ALN = f"{COL}/alineamientos_codones_NAC.fa"
RESUMEN = f"{COL}/resumen_kaks_Egl.json"
CACHE = f"{COL}/cache_e3_10_kaks.pkl"   # resultados crudos por par, para reanudar (--cache)
LAMBDA = 6.5e-9          # sustituciones/sitio/año (artículo 1)
KS_SAT = 2.0             # Ks > 2: saturado, no se data
MAX_DIF = 0.02           # diferencias aceptadas entre traducción del CDS y proteína
MIN_COD = 30             # codones alineados mínimos para estimar
HILOS = 10
METODOS = ("YN00", "NG86", "YN00_Biopython")   # YN00 = yn00 de PAML (principal)
ACGT = set("ACGT")
PARADA = {"TAA", "TAG", "TGA"}
PEP, COD = {}, {}


def corto(p):
    return p.replace(RAIZ + "/", "")


def leer_fasta(fuente, es_texto=False):
    seqs, k = {}, None
    for l in (fuente.splitlines() if es_texto else open(fuente)):
        l = l.strip()
        if not l:
            continue
        if l.startswith(">"):
            k = l[1:].split()[0]
            seqs[k] = []
        else:
            seqs[k].append(l)
    return {k: "".join(v) for k, v in seqs.items()}


def conciliar(cds, pep):
    """(codones, diferencias, marco) o None."""
    cds, pep = cds.upper(), pep.rstrip("*").upper()
    mejor = None
    for marco in (0, 1, 2):
        sub = cds[marco:]
        tr = str(Seq(sub[: len(sub) // 3 * 3]).translate())
        i = tr.find(pep)
        if i >= 0:
            ini = marco + 3 * i
            return [cds[ini + 3 * j: ini + 3 * j + 3] for j in range(len(pep))], 0, marco
        if len(tr) >= len(pep):
            dif = sum(a != b for a, b in zip(tr, pep))
            if mejor is None or dif < mejor[1]:
                mejor = ([cds[marco + 3 * j: marco + 3 * j + 3] for j in range(len(pep))], dif, marco)
    return mejor if mejor is not None and mejor[1] <= MAX_DIF * len(pep) else None


def yn00_paml(s1, s2, d):
    open(f"{d}/par.nuc", "w").write(f"  2  {len(s1)}\nsecA  {s1}\nsecB  {s2}\n")
    open(f"{d}/yn00.ctl", "w").write("seqfile = par.nuc\noutfile = yn.out\nverbose = 0\nicode = 0\n"
                                     "weighting = 0\ncommonf3x4 = 0\n")
    subprocess.run([f"{BIN}/yn00", "yn00.ctl"], cwd=d, capture_output=True, text=True, timeout=300)
    txt = open(f"{d}/yn.out").read() if os.path.exists(f"{d}/yn.out") else ""
    out = dict(Ka_YN00=np.nan, Ks_YN00=np.nan, kappa_YN00=np.nan)
    i = txt.find("Yang & Nielsen (2000)")
    if i >= 0:
        for l in txt[i:].splitlines():
            if re.match(r"^\s*2\s+1\s", l):
                v = [float(x) + 0.0 if re.match(r"^-?\d", x) else np.nan for x in re.findall(r"-?nan|-?inf|-?\d+\.\d+", l)]
                if len(v) >= 9:      # S N t kappa omega dN SE dS SE
                    out = dict(Ka_YN00=v[5], Ks_YN00=v[7], kappa_YN00=v[3])
                break
    return out, txt


def procesar(tarea):
    clave, g1, g2 = tarea
    r = dict(clave=clave, locus1=g1, locus2=g2, estado="ok")
    if g1 not in COD or g2 not in COD:
        r["estado"] = "sin_conciliar"
        return r, None
    with tempfile.TemporaryDirectory(prefix="kaks_") as d:
        open(f"{d}/p.fa", "w").write(f">a\n{PEP[g1]}\n>b\n{PEP[g2]}\n")
        m = subprocess.run([f"{BIN}/mafft", "--quiet", "--localpair", "--maxiterate", "1000", "--thread", "1", "p.fa"],
                           cwd=d, capture_output=True, text=True)
        al = leer_fasta(m.stdout, es_texto=True)
        a1, a2 = al.get("a", "").upper(), al.get("b", "").upper()
        if len(a1.replace("-", "")) != len(PEP[g1]) or len(a2.replace("-", "")) != len(PEP[g2]):
            r["estado"] = "error_retrotraduccion"
            return r, None
        c1, c2 = iter(COD[g1]), iter(COD[g2])
        k1, k2, iguales = [], [], 0
        for x, y in zip(a1, a2):
            cx = next(c1) if x != "-" else None
            cy = next(c2) if y != "-" else None
            if cx is None or cy is None or set(cx) - ACGT or set(cy) - ACGT or cx in PARADA or cy in PARADA:
                continue
            k1.append(cx)
            k2.append(cy)
            iguales += x == y
        L = len(k1)
        s1, s2 = "".join(k1), "".join(k2)
        r.update(codones_alineados=L, identidad_aa_pct=round(100 * iguales / L, 1) if L else np.nan,
                 cobertura_pct=round(100 * L / min(len(PEP[g1]), len(PEP[g2])), 1))
        if L < MIN_COD:
            r["estado"] = "alineamiento_corto"
            return r, (s1, s2)
        aln = Alignment([s1, s2], coordinates=np.array([[0, 3 * L], [0, 3 * L]]))
        for met, col in (("NG86", "NG86"), ("YN00", "YN00_Biopython")):
            try:
                dn, ds = calculate_dn_ds(aln, method=met)
            except Exception:
                dn, ds = np.nan, np.nan
            r[f"Ka_{col}"] = dn if np.isfinite(dn) and dn >= 0 else np.nan     # NG86 devuelve -1 si satura
            r[f"Ks_{col}"] = ds if np.isfinite(ds) and ds >= 0 else np.nan
        r.update(yn00_paml(s1, s2, d)[0])
    return r, (s1, s2)


def codeml_par(tarea, devolver_texto=False):
    """lnL con omega libre y con omega = 1 (codeml, runmode -2, F3x4)."""
    clave, s1, s2 = tarea
    out = dict(clave=clave, lnL_omega_libre=np.nan, lnL_omega_1=np.nan, omega_ML=np.nan, dN_ML=np.nan, dS_ML=np.nan)
    txt_libre = ""
    with tempfile.TemporaryDirectory(prefix="cml_") as d:
        open(f"{d}/par.nuc", "w").write(f"  2  {len(s1)}\nsecA  {s1}\nsecB  {s2}\n")
        for nombre, fijo in (("libre", 0), ("1", 1)):
            open(f"{d}/w{nombre}.ctl", "w").write(
                f"seqfile = par.nuc\noutfile = w{nombre}.out\nrunmode = -2\nseqtype = 1\nCodonFreq = 2\nmodel = 0\n"
                f"NSsites = 0\nicode = 0\nfix_kappa = 0\nkappa = 2\nfix_omega = {fijo}\nomega = {1 if fijo else 0.5}\n"
                "cleandata = 1\nverbose = 0\nnoisy = 0\ngetSE = 0\nRateAncestor = 0\n")
            subprocess.run([f"{BIN}/codeml", f"w{nombre}.ctl"], cwd=d, capture_output=True, text=True, timeout=900)
            txt = open(f"{d}/w{nombre}.out").read() if os.path.exists(f"{d}/w{nombre}.out") else ""
            m = re.findall(r"lnL\s*=\s*(-?\d+\.\d+)", txt)
            out[f"lnL_omega_{nombre}"] = float(m[-1]) if m else np.nan
            if not fijo:
                txt_libre = txt
                m = re.search(r"dN/dS\s*=\s*(-?[\d.]+|nan)\s+dN\s*=\s*(-?[\d.]+|nan)\s+dS\s*=\s*(-?[\d.]+|nan)", txt)
                if m:
                    out.update(omega_ML=float(m.group(1)), dN_ML=float(m.group(2)), dS_ML=float(m.group(3)))
    return (out, txt_libre) if devolver_texto else out


def derivar(df):
    for met in METODOS:
        ka, ks = df[f"Ka_{met}"], df[f"Ks_{met}"]
        df[f"KaKs_{met}"] = np.where((ks > 0) & np.isfinite(ka), ka / ks.where(ks > 0), np.nan)
        df[f"T_Mya_{met}"] = np.where(ks <= KS_SAT, ks / (2 * LAMBDA) / 1e6, np.nan)
    return df


def picos(x, hasta=3.0):
    x = np.asarray(x, float)
    x = x[np.isfinite(x) & (x > 0) & (x <= hasta)]
    kde = gaussian_kde(x)
    g = np.linspace(0, hasta, 601)
    dens = kde(g)
    out = []
    for i in range(1, len(g) - 1):
        if dens[i] > dens[i - 1] and dens[i] >= dens[i + 1] and dens[i] >= 0.1 * dens.max():
            mitad = dens[i] / 2
            a = i
            while a > 0 and dens[a] > mitad:
                a -= 1
            b = i
            while b < len(g) - 1 and dens[b] > mitad:
                b += 1
            out.append(dict(Ks=round(float(g[i]), 3), densidad_rel=round(float(dens[i] / dens.max()), 3),
                            mitad_altura=[round(float(g[a]), 3), round(float(g[b]), 3)]))
    return out, len(x), float(kde.factor)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--prueba", type=int, default=0)
    ap.add_argument("--cache", action="store_true")
    arg = ap.parse_args()
    t0 = time.time()
    print("== Ka/Ks de pares NAC duplicados de E. globulus (e3_10_kaks_Egl.py)", flush=True)
    t46 = pd.read_csv(T46)
    t46["acierto_DIAMOND"] = t46.acierto_DIAMOND.astype(bool)
    t46["tipo_par"] = np.where(t46.modo_par_MCScanX.notna(), t46.modo_par_MCScanX,
                               np.where(t46.acierto_DIAMOND, "a_200kb_sin_modo", "a_200kb_sin_acierto"))
    nac = t46.reset_index(drop=True)
    an = pd.read_csv(ANCLAS, sep="\t")
    if arg.prueba:
        nac = pd.concat([nac[nac.acierto_DIAMOND].head(arg.prueba), nac[~nac.acierto_DIAMOND].head(arg.prueba)]).reset_index(drop=True)
        an = an.head(arg.prueba)
    print(f"pares NAC analizados: {len(nac)} | con acierto DIAMOND: {int(nac.acierto_DIAMOND.sum())} | pares ancla: {len(an)}", flush=True)
    print("tipos de par NAC:", nac.tipo_par.value_counts().to_dict(), flush=True)

    pep, cds = leer_fasta(f"{ISO}/Egl.pep"), leer_fasta(f"{ISO}/Egl.cds")
    genes = sorted(set(nac.locus1) | set(nac.locus2) | set(an.gen1) | set(an.gen2))
    con = dict(genes=len(genes), sin_cds=0, exactos=0, con_diferencias=0, sin_conciliar=0, marco_1_o_2=0, genes_sin_conciliar=[])
    for g in genes:
        if g not in cds or g not in pep:
            con["sin_cds"] += 1
            continue
        c = conciliar(cds[g], pep[g])
        if c is None:
            con["sin_conciliar"] += 1
            con["genes_sin_conciliar"].append(g)
            continue
        COD[g], PEP[g] = c[0], pep[g].rstrip("*").upper()
        con["exactos" if c[1] == 0 else "con_diferencias"] += 1
        con["marco_1_o_2"] += c[2] != 0
    print("conciliación CDS-proteína:", con, flush=True)

    tareas = [(f"NAC:{i}", r.locus1, r.locus2) for i, r in nac.iterrows()] + \
             [(f"ancla:{i}", r.gen1, r.gen2) for i, r in an.iterrows()]
    if arg.prueba:
        r, (s1, s2) = procesar(tareas[0])
        with tempfile.TemporaryDirectory(prefix="kaks_") as d:
            txt = yn00_paml(s1, s2, d)[1]
        i = txt.find("Yang & Nielsen (2000)")
        print("--- salida yn00 (PAML), primer par:\n" + txt[i: i + 400], flush=True)
        c, txt = codeml_par(("NAC:0", s1, s2), devolver_texto=True)
        print("--- salida codeml (omega libre), primer par (final):\n" + txt[-700:], "\n", c, flush=True)
    t1 = time.time()
    cache = pd.read_pickle(CACHE) if arg.cache and os.path.exists(CACHE) else {}
    pend = [t for t in tareas if (t[1], t[2]) not in cache]
    print(f"pares en caché: {len(tareas) - len(pend)} | por alinear: {len(pend)}", flush=True)
    if pend:
        with Pool(HILOS) as pool:
            for k, (r, sq) in enumerate(pool.imap(procesar, pend, chunksize=4), 1):
                cache[(r["locus1"], r["locus2"])] = (r, sq)
                if k % 250 == 0 or k == len(pend):
                    print(f"  {k}/{len(pend)} pares ({time.time() - t1:.0f} s)", flush=True)
    res, alins = [], {}
    for clave, g1, g2 in tareas:
        r, sq = cache[(g1, g2)]
        res.append(dict(r, clave=clave))
        if sq is not None:
            alins[clave] = sq
    assert len(res) == len(tareas), (len(res), len(tareas))
    gpar = {r["clave"]: (r["locus1"], r["locus2"]) for r in res if r["clave"].startswith("NAC:") and r["estado"] == "ok"}
    pc = [c for c in gpar if ("codeml",) + gpar[c] not in cache]
    print(f"codeml (omega libre frente a omega = 1): {len(gpar)} pares NAC, {len(pc)} por calcular", flush=True)
    if pc:
        t2 = time.time()
        with Pool(HILOS) as pool:
            for k, out in enumerate(pool.imap(codeml_par, [(c,) + alins[c] for c in pc], chunksize=2), 1):
                cache[("codeml",) + gpar[out["clave"]]] = out
                if k % 100 == 0 or k == len(pc):
                    print(f"  codeml {k}/{len(pc)} ({time.time() - t2:.0f} s)", flush=True)
    if (pend or pc) and not arg.prueba:
        pd.to_pickle(cache, CACHE)
    cml = pd.DataFrame([dict(cache[("codeml",) + gpar[c]], clave=c) for c in gpar])
    R = derivar(pd.DataFrame(res))
    if len(cml):
        R = R.merge(cml, on="clave", how="left")
        R["LRT_2dlnL"] = (2 * (R.lnL_omega_libre - R.lnL_omega_1)).clip(lower=0)
        R["p_LRT"] = chi2.sf(R.LRT_2dlnL, 1)
        m = R.p_LRT.notna()
        R.loc[m, "q_LRT"] = multipletests(R.loc[m, "p_LRT"], method="fdr_bh")[1]
    print("estados:", R.estado.value_counts().to_dict(), flush=True)
    if arg.prueba:
        cols = ["clave", "codones_alineados", "identidad_aa_pct"] + [f"{a}_{m}" for m in METODOS for a in ("Ka", "Ks")] + \
               ["omega_ML", "dS_ML", "LRT_2dlnL", "p_LRT"]
        print(R[cols].round(4).to_string(index=False))
        print(f"tiempo por par: {(time.time() - t1) / len(tareas):.2f} s con {HILOS} procesos | total {time.time() - t0:.0f} s")
        return

    # ---- T47: pares NAC ----
    rn = R[R.clave.str.startswith("NAC:")].copy()
    rn.index = rn.clave.str.split(":").str[1].astype(int)
    base = nac[["simbolo1", "simbolo2", "locus1", "locus2", "cromosoma1", "cromosoma2", "subfamilia1", "subfamilia2",
                "misma_subfamilia", "tipo_par", "modo_par_MCScanX", "tandem_art3", "grupo_art3", "acierto_DIAMOND",
                "mismo_cromosoma", "distancia_kb", "genes_intercalados", "NAC_intercalados", "identidad_pct",
                "mismo_bloque_T15"]].rename(columns={"identidad_pct": "identidad_DIAMOND_pct"})
    t47 = base.join(rn.drop(columns=["clave", "locus1", "locus2"]))
    t47["Ks_saturado"] = ~(t47.Ks_YN00 <= KS_SAT)
    t47["KaKs_puntual"] = np.select([t47.KaKs_YN00 < 1, t47.KaKs_YN00 > 1], ["<1", ">1"], "no estimable")
    t47["seleccion_LRT"] = np.select(
        [t47.q_LRT.isna(), (t47.q_LRT < 0.05) & (t47.omega_ML < 1), (t47.q_LRT < 0.05) & (t47.omega_ML > 1)],
        ["no estimable", "purificadora", "positiva"], "no distinta de neutral")
    t47["Ks_menor_0.01"] = t47.Ks_YN00 < 0.01
    n209 = "modelo génico de EglNAC209 probablemente partido (ACJRO7_019604 contiene el C-terminal)"
    nks = "Ks < 0,01: duplicación muy reciente o posible artefacto de ensamblaje"
    t47["nota"] = ["; ".join(x for x, f in ((n209, a), (nks, b)) if f) for a, b in
                   zip(t47.simbolo1.eq("EglNAC209") | t47.simbolo2.eq("EglNAC209"), t47["Ks_menor_0.01"])]
    orden = {"WGD/segmental": 0, "tandem": 1, "proximal": 2, "a_200kb_sin_modo": 3, "a_200kb_sin_acierto": 4}
    t47 = t47.sort_values(["tipo_par", "Ks_YN00"], key=lambda s: s.map(orden) if s.name == "tipo_par" else s)
    num = [c for c in t47.columns if c.startswith(("Ka_", "Ks_", "KaKs_", "kappa", "omega", "dN_", "dS_", "lnL", "LRT"))
           and pd.api.types.is_float_dtype(t47[c])]
    t47[num] = t47[num].astype(float).round(4)
    t47[[c for c in t47.columns if c.startswith("T_Mya")]] = t47[[c for c in t47.columns if c.startswith("T_Mya")]].round(2)
    t47.to_csv(T47, index=False)
    with open(ALN, "w") as fh:
        for i, r in nac.iterrows():
            sq = alins.get(f"NAC:{i}")
            if sq:
                fh.write(f">{r.simbolo1}|{r.locus1}|par_{r.simbolo1}_{r.simbolo2}\n{sq[0]}\n>{r.simbolo2}|{r.locus2}|par_{r.simbolo1}_{r.simbolo2}\n{sq[1]}\n")

    # ---- pares ancla del genoma ----
    ra = R[R.clave.str.startswith("ancla:")].copy()
    ra.index = ra.clave.str.split(":").str[1].astype(int)
    ga = an[["bloque", "crom1", "crom2", "n_anclas", "gen1", "gen2", "NAC1", "NAC2"]].join(
        ra.drop(columns=["clave", "locus1", "locus2"]))
    ga.round(4).to_csv(KS_GEN, sep="\t", index=False)

    # ---- controles y contrastes ----
    ok = t47[t47.estado == "ok"]
    rng = np.random.default_rng(1)
    muestra = [tareas[i] for i in sorted(rng.choice(len(nac), size=min(5, len(nac)), replace=False))]
    recal = pd.DataFrame([procesar(t)[0] for t in muestra]).set_index("clave")
    orig = R.set_index("clave").loc[recal.index]
    dif_recalculo = float(np.nanmax(np.abs(recal[["Ka_YN00", "Ks_YN00"]].values - orig[["Ka_YN00", "Ks_YN00"]].values)))
    concord = {}
    for nombre, df in (("NAC", ok), ("anclas", ga[ga.estado == "ok"])):
        m = df[df.Ks_YN00 <= KS_SAT]
        b = m.dropna(subset=["Ks_YN00_Biopython"])
        n_ = m.dropna(subset=["Ks_NG86"])
        concord[nombre] = dict(
            pares_PAML_no_saturados=int(len(m)),
            YN00_Biopython_no_estimable=int(df.Ks_YN00_Biopython.isna().sum()), de=int(len(df)),
            rho_Ks_YN00_PAML_vs_Biopython=round(float(spearmanr(b.Ks_YN00, b.Ks_YN00_Biopython)[0]), 4),
            mediana_dif_rel_Ks_PAML_vs_Biopython_pct=round(float(100 * ((b.Ks_YN00 - b.Ks_YN00_Biopython).abs() / b.Ks_YN00).median()), 2),
            NG86_no_estimable=int(df.Ks_NG86.isna().sum()),
            rho_Ks_YN00_PAML_vs_NG86=round(float(spearmanr(n_.Ks_YN00, n_.Ks_NG86)[0]), 4),
            rho_Ka_YN00_PAML_vs_NG86=round(float(spearmanr(n_.Ka_YN00, n_.Ka_NG86)[0]), 4),
            mediana_dif_rel_Ks_PAML_vs_NG86_pct=round(float(100 * ((n_.Ks_YN00 - n_.Ks_NG86).abs() / n_.Ks_YN00).median()), 2))
    por_modo = {}
    for modo, df in ok.groupby("tipo_par"):
        ns = df[~df.Ks_saturado]
        por_modo[modo] = dict(
            pares=int(len(df)), Ks_saturado=int(df.Ks_saturado.sum()), seleccion_LRT=df.seleccion_LRT.value_counts().to_dict(),
            KaKs_YN00_mediana=round(float(df.KaKs_YN00.median()), 3),
            KaKs_YN00_rango=[round(float(df.KaKs_YN00.min()), 3), round(float(df.KaKs_YN00.max()), 3)],
            KaKs_mayor_1=int((df.KaKs_YN00 > 1).sum()),
            Ks_YN00_mediana=round(float(df.Ks_YN00.median()), 3),
            T_Mya_mediana_no_saturados=round(float(ns.T_Mya_YN00.median()), 1) if len(ns) else None,
            T_Mya_rango_no_saturados=[round(float(ns.T_Mya_YN00.min()), 1), round(float(ns.T_Mya_YN00.max()), 1)] if len(ns) else None)
    p_gen, n_gen, bw = picos(ga.Ks_YN00)
    bloq = ga[ga.estado == "ok"].groupby("bloque").agg(n=("Ks_YN00", "size"), Ks=("Ks_YN00", "median"))
    p_blo, n_blo, _ = picos(bloq.Ks)
    p_ng86, _, _ = picos(ga.Ks_NG86)
    seg = ok[ok.tipo_par == "WGD/segmental"]
    pico = max(p_gen, key=lambda p: p["densidad_rel"]) if p_gen else None
    en_pico = seg[(seg.Ks_YN00 >= pico["mitad_altura"][0]) & (seg.Ks_YN00 <= pico["mitad_altura"][1])] if pico else seg.iloc[:0]
    def mw(a, b):
        a, b = a.Ks_YN00.dropna(), b.Ks_YN00.dropna()
        if not (len(a) and len(b)):
            return None
        r = mannwhitneyu(a, b)
        return dict(n=[int(len(a)), int(len(b))], U=float(r.statistic), AUC=round(float(r.statistic) / (len(a) * len(b)), 3),
                    p=float(r.pvalue))
    locales = ok[(ok.tipo_par != "WGD/segmental") & ok.acierto_DIAMOND]
    sin_ac = ok[~ok.acierto_DIAMOND]
    anc_ok = ga[(ga.estado == "ok")]
    resumen = dict(
        script="e3_10_kaks_Egl.py", segundos=round(time.time() - t0),
        parametros=dict(MAFFT="--localpair --maxiterate 1000 (L-INS-i), retrotraducción a codones",
                        YN00_principal="PAML yn00, icode 0, weighting 0, commonf3x4 0",
                        Biopython="Bio.Align.analysis.calculate_dn_ds, NG86 (k=1) y YN00 (secundarios)", lambda_sitio_anio=LAMBDA, Ks_saturado=KS_SAT,
                        codones_minimos=MIN_COD, diferencias_CDS_proteina_max=MAX_DIF),
        conciliacion_CDS_proteina=con,
        pares_NAC=dict(en_T46=int(len(t46)), sin_acierto_DIAMOND=int((~t46.acierto_DIAMOND).sum()),
                       analizados=int(len(t47)), estados=t47.estado.value_counts().to_dict(),
                       seleccion_LRT=ok.seleccion_LRT.value_counts().to_dict(),
                       seleccion_LRT_con_acierto_DIAMOND=ok[ok.acierto_DIAMOND].seleccion_LRT.value_counts().to_dict(),
                       KaKs_YN00_menor_1=int((ok.KaKs_YN00 < 1).sum()), KaKs_YN00_mayor_1=int((ok.KaKs_YN00 > 1).sum()),
                       KaKs_no_estimable=int(ok.KaKs_YN00.isna().sum()), Ks_saturado=int(ok.Ks_saturado.sum()),
                       Ks_menor_0_01=int(ok["Ks_menor_0.01"].sum()), por_modo=por_modo),
        segmentales_NAC=dict(
            pares=int(len(seg)), Ks_YN00_mediana=round(float(seg.Ks_YN00.median()), 3),
            Ks_YN00_rango=[round(float(seg.Ks_YN00.min()), 3), round(float(seg.Ks_YN00.max()), 3)],
            en_pico_principal_del_genoma=int(len(en_pico)),
            mann_whitney_Ks_segmental_vs_locales_con_acierto=mw(seg, locales),
            mann_whitney_Ks_segmental_vs_200kb_sin_acierto=mw(seg, sin_ac),
            mann_whitney_Ks_locales_con_acierto_vs_200kb_sin_acierto=mw(locales, sin_ac)),
        anclas_genoma=dict(
            pares=int(len(ga)), estimados=int(len(anc_ok)), estados=ga.estado.value_counts().to_dict(),
            Ks_YN00_mediana=round(float(anc_ok.Ks_YN00.median()), 3),
            Ks_YN00_saturado=int((~(anc_ok.Ks_YN00 <= KS_SAT)).sum()),
            KaKs_YN00_mediana=round(float(anc_ok.KaKs_YN00.median()), 3),
            picos_KDE_anclas_YN00_PAML=dict(n=n_gen, factor_banda=round(bw, 4), picos=p_gen),
            picos_KDE_anclas_NG86=p_ng86,
            picos_KDE_medianas_por_bloque=dict(bloques=n_blo, picos=p_blo),
            T_Mya_del_pico_principal=round(pico["Ks"] / (2 * LAMBDA) / 1e6, 1) if pico else None),
        concordancia=concord,
        control_recalculo_5_pares=dict(claves=list(recal.index), dif_abs_max_Ka_Ks_YN00=dif_recalculo))
    json.dump(resumen, open(RESUMEN, "w"), ensure_ascii=False, indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o))
    print(json.dumps(resumen, ensure_ascii=False, indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o)), flush=True)
    print(f"escritos: {corto(T47)} ({len(t47)} filas), {corto(KS_GEN)} ({len(ga)} filas), {corto(ALN)}, {corto(RESUMEN)}"
          f" | {time.time() - t0:.0f} s", flush=True)


if __name__ == "__main__":
    main()
