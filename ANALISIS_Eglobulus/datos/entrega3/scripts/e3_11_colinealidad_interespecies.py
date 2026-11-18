"""Entrega 3, paso 4.3 (artículo 3): colinealidad de E. globulus con E. grandis v2 y P. trichocarpa v4.1.

Entradas
  datos/entrega3/isoformas/{Egl,Egr,Ptr}.{bed,cds,pep}   isoforma representativa por gen, solo cromosomas
  referencia/PF02365.hmm                                  Pfam NAM (PF02365.21, InterPro), umbral GA 27 bits
  tablas/T01_censo_maestro_189.csv, tablas/T45_modos_duplicacion_NAC189.csv
  ../Hussey_TablaS4_189simbolos.csv                       símbolo EgrNAC -> locus Eucgr. (v1.0)
  ../_respaldo_app/colinealidad_grandis_globulus.csv      mapeo posicional anterior (98 símbolos del manuscrito)
  referencia/Eglobulus_ASM4622685v1/sequence_report.json  largos de cromosoma del ensamblado en uso
  referencia/Eglobulus_GCA_014182545.1/sequence_report.json  largos del ensamblado declarado en el manuscrito

Método
  1. Familia NAC de cada especie: hmmsearch --cut_ga con PF02365 sobre la isoforma representativa. Control: en
     E. globulus devuelve exactamente los 189 loci de T01.
  2. JCVI catalog ortholog (MCscan; el programa que cita el artículo 3): LAST sobre CDS (lastal -u0 -i3G), filtro
     C-score >= 0,7, tándems fusionados (<= 10 genes), bloques de >= 4 anclas con huecos <= 20 genes y rescate
     (liftover) de pares dentro de los bloques, marcados con 'L'. Se usa el archivo .lifted.anchors.
  3. Ortólogo colineal de cada NAC de E. globulus: su mejor socio dentro de los bloques (mayor puntaje LAST). El
     rescate también trae a los parálogos en tándem vecinos, que quedan como socios colineales. Recíproco si el gen
     de E. globulus es a su vez el mejor socio colineal de su ortólogo.
  4. Contraste (a): el ortólogo colineal en E. grandis frente al co-ortólogo EgrNAC de T01 (símbolo de Hussey ->
     locus Eucgr. -> EUGRSUZ_ de la v2).
     Contraste (b): el mapeo posicional anterior. Sus coordenadas en E. globulus no son del ensamblado en uso, así que
     se compara el cromosoma y la identidad del gen: el locus Eucgr. de cada símbolo del manuscrito -> su par colineal
     en E. globulus -> ¿es el locus que T01 asigna a ese símbolo por secuencia?

Ejecución (en primer plano, por lotes)
  python e3_11_colinealidad_interespecies.py --solo Egr    # JCVI E. globulus - E. grandis
  python e3_11_colinealidad_interespecies.py --solo Ptr    # JCVI E. globulus - P. trichocarpa
  python e3_11_colinealidad_interespecies.py               # tablas y resumen (reutiliza las anclas)

Salidas
  tablas/T48_ortologos_colineales_NAC189.csv          un renglón por NAC de E. globulus
  tablas/T49_pares_colineales_NAC_interespecies.csv   pares de anclas con un NAC en al menos un lado
  tablas/T50_contraste_mapeo_posicional.csv           los 98 símbolos del mapeo anterior frente a JCVI
  datos/entrega3/colinealidad/jcvi/                   anclas, LAST filtrado, gráficos de puntos, salidas de hmmsearch
  datos/entrega3/colinealidad/resumen_colinealidad_interespecies.json
"""

import argparse
import json
import os
import shutil
import subprocess
import tempfile
import time

import numpy as np
import pandas as pd

RAIZ = (os.environ.get("NAC_ROOT", ".") + "/ANALISIS_Eglobulus")
REV = os.path.dirname(RAIZ)
ISO = f"{RAIZ}/datos/entrega3/isoformas"
COL = f"{RAIZ}/datos/entrega3/colinealidad"
JC = f"{COL}/jcvi"
BIN = (os.environ.get("CONDA_ENVS", os.path.expanduser("~/conda/envs")) + "/nac3/bin")
HMMSEARCH = (os.environ.get("CONDA_ENVS", os.path.expanduser("~/conda/envs")) + "/hmmer/bin/hmmsearch")
HMM = f"{RAIZ}/referencia/PF02365.hmm"
T01 = f"{RAIZ}/tablas/T01_censo_maestro_189.csv"
T45 = f"{RAIZ}/tablas/T45_modos_duplicacion_NAC189.csv"
HUSSEY = f"{REV}/Hussey_TablaS4_189simbolos.csv"
ANTIGUO = f"{REV}/_respaldo_app/colinealidad_grandis_globulus.csv"
SR_ACT = f"{RAIZ}/referencia/Eglobulus_ASM4622685v1/sequence_report.json"
SR_ANT = f"{RAIZ}/referencia/Eglobulus_GCA_014182545.1/sequence_report.json"
T48 = f"{RAIZ}/tablas/T48_ortologos_colineales_NAC189.csv"
T49 = f"{RAIZ}/tablas/T49_pares_colineales_NAC_interespecies.csv"
T50 = f"{RAIZ}/tablas/T50_contraste_mapeo_posicional.csv"
RESUMEN = f"{COL}/resumen_colinealidad_interespecies.json"
OTRAS = {"Egr": "E. grandis v2 (GCA_000612305.2)", "Ptr": "P. trichocarpa v4.1 (GCF_000002775.5)"}
CSCORE = 0.7


def corto(p):
    return p.replace(RAIZ + "/", "")


def primera_linea(cmd, clave=""):
    r = subprocess.run(cmd, capture_output=True, text=True)
    ls = [x for x in (r.stdout + r.stderr).splitlines() if clave in x]
    return ls[0].strip("# ").strip() if ls else ""


def correr_jcvi(b, trabajo, cpus):
    """JCVI ortholog Egl-b en un directorio local (LAST trabaja mal sobre la unidad de red) y copia a Z."""
    os.makedirs(trabajo, exist_ok=True)
    os.makedirs(JC, exist_ok=True)
    for s in ("Egl", b):
        for ext in ("bed", "cds"):
            o, d = f"{ISO}/{s}.{ext}", f"{trabajo}/{s}.{ext}"
            if not (os.path.exists(d) and os.path.getsize(d) == os.path.getsize(o)):
                shutil.copyfile(o, d)
    env = dict(os.environ, PATH=BIN + os.pathsep + os.environ.get("PATH", ""), MPLBACKEND="Agg")
    cmd = [f"{BIN}/python", "-m", "jcvi.compara.catalog", "ortholog", "Egl", b,
           f"--cscore={CSCORE}", "--no_strip_names", f"--cpus={cpus}", "--notex"]
    print(f"$ (en {trabajo}) jcvi.compara.catalog ortholog Egl {b} --cscore={CSCORE} --no_strip_names "
          f"--cpus={cpus} --notex", flush=True)
    t = time.time()
    r = subprocess.run(cmd, cwd=trabajo, env=env, capture_output=True, text=True)
    print(f"  JCVI Egl-{b}: código {r.returncode}, {time.time() - t:.0f} s", flush=True)
    for ln in (r.stderr + r.stdout).splitlines()[-12:]:
        print("   |", ln[:220])
    pre = f"Egl.{b}"
    for s in (f"{pre}.last", f"{pre}.last.filtered", f"{pre}.anchors", f"{pre}.lifted.anchors"):
        assert os.path.exists(f"{trabajo}/{s}"), f"falta {s}"
    with open(f"{trabajo}/{pre}.last") as fh:
        n_last = sum(1 for _ in fh)
    for s in (f"{pre}.last.filtered", f"{pre}.anchors", f"{pre}.lifted.anchors", f"{pre}.pdf"):
        if os.path.exists(f"{trabajo}/{s}"):
            shutil.copyfile(f"{trabajo}/{s}", f"{JC}/{s}")
    json.dump(dict(especie=b, codigo=r.returncode, segundos=round(time.time() - t), lineas_last=n_last,
                   comando=" ".join(cmd[1:])), open(f"{JC}/{pre}.ejecucion.json", "w"), indent=1)
    print(f"  copiado a {corto(JC)}: {pre}.last.filtered, .anchors, .lifted.anchors, .pdf "
          f"(el .last completo, {n_last} líneas, queda en {trabajo})", flush=True)


def nac_hmm(sp):
    os.makedirs(f"{JC}/hmm", exist_ok=True)
    dom = f"{JC}/hmm/{sp}_PF02365_cutga.domtbl"
    if not os.path.exists(dom):
        subprocess.run([HMMSEARCH, "-o", os.devnull, "--domtblout", dom, "--cut_ga", "--cpu", "8", HMM,
                        f"{ISO}/{sp}.pep"], check=True)
    return {ln.split()[0] for ln in open(dom) if not ln.startswith("#")}


def leer_bed(sp):
    b = pd.read_csv(f"{ISO}/{sp}.bed", sep="\t", header=None, names=["chr", "ini", "fin", "gen", "s", "hebra"])
    return b.set_index("gen")


def leer_anclas(p):
    filas, blo = [], 0
    for ln in open(p):
        if ln.startswith("#"):
            blo += 1
            continue
        c = ln.split()
        if len(c) >= 3:
            filas.append((blo, c[0], c[1], c[2].endswith("L"), float(c[2].rstrip("L"))))
    d = pd.DataFrame(filas, columns=["bloque", "gen_Egl", "gen_otro", "rescatado", "puntaje"])
    d["anclas_en_bloque"] = d.groupby("bloque").gen_Egl.transform("size")
    return d


def largos(p):
    rep = json.load(open(p))["reports"]
    return {int(r["chr_name"]): int(r["length"]) for r in rep
            if r.get("role") == "assembled-molecule" and str(r.get("chr_name", "")).isdigit()}


def registros(df):
    return df.astype(object).where(df.notna(), None).to_dict("records")


def unir(xs):
    xs = sorted(set(x for x in xs if isinstance(x, str) and x))
    return ";".join(xs) if xs else np.nan


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--trabajo", default=os.path.join(tempfile.gettempdir(), "jcvi_Egl"))
    ap.add_argument("--cpus", type=int, default=12)
    ap.add_argument("--solo", choices=sorted(OTRAS), help="solo correr JCVI frente a esta especie")
    arg = ap.parse_args()
    t0 = time.time()
    if arg.solo:
        correr_jcvi(arg.solo, arg.trabajo, arg.cpus)
        return

    versiones = dict(
        jcvi=primera_linea([f"{BIN}/python", "-c", "import jcvi; print('jcvi', jcvi.__version__)"], "jcvi"),
        last=primera_linea([f"{BIN}/lastal", "--version"], "lastal"),
        hmmer=primera_linea([HMMSEARCH, "-h"], "HMMER"),
        perfil=" ".join(primera_linea(["grep", "-m1", f"^{k}", HMM], k) for k in ("ACC", "GA")))
    print("versiones:", versiones)

    # 1. Familia NAC por especie
    nac = {sp: nac_hmm(sp) for sp in ("Egl", "Egr", "Ptr")}
    bed = {sp: leer_bed(sp) for sp in ("Egl", "Egr", "Ptr")}
    t01 = pd.read_csv(T01)
    t45 = pd.read_csv(T45)
    assert nac["Egl"] == set(t01.locus), "el perfil no reproduce los 189 de T01"
    for sp in ("Egr", "Ptr"):
        assert nac[sp] <= set(bed[sp].index)
    print("NAC por hmmsearch --cut_ga:", {sp: len(v) for sp, v in nac.items()})

    hus = pd.read_csv(HUSSEY)
    hus["locus_v2"] = hus.locus_v1_0.str.strip().str.replace("Eucgr.", "EUGRSUZ_", regex=False)
    sim_de_egr = dict(zip(hus.locus_v2, hus.simbolo_EgrNAC.str.strip()))
    egr_de_sim = dict(zip(hus.simbolo_EgrNAC.str.strip(), hus.locus_v2))
    hus_v2 = dict(simbolos=len(hus), en_cromosomas_v2=int(hus.locus_v2.isin(bed["Egr"].index).sum()),
                  con_NAM_cutga=int(hus.locus_v2.isin(nac["Egr"]).sum()),
                  NAC_v2_fuera_de_Hussey=int(len(nac["Egr"] - set(hus.locus_v2))))
    print("Hussey S4 en la v2:", hus_v2)

    base = t45[["simbolo", "locus", "subfamilia", "modo_duplicacion"]].merge(
        t01[["locus", "coortologo_EgrNAC", "regla", "simbolo_manuscrito_2026"]], on="locus", how="left",
        validate="1:1")
    assert len(base) == 189
    base.insert(2, "cromosoma", base.locus.map(bed["Egl"].chr))
    base.insert(3, "inicio", base.locus.map(bed["Egl"].ini))
    base["coortologo_v2"] = base.coortologo_EgrNAC.map(egr_de_sim)
    sim_egl = dict(zip(base.locus, base.simbolo))

    # 2. JCVI
    anc, genoma = {}, {}
    for b in OTRAS:
        p = f"{JC}/Egl.{b}.lifted.anchors"
        if not os.path.exists(p):
            correr_jcvi(b, arg.trabajo, arg.cpus)
        d = leer_anclas(p)
        d["chr_Egl"] = d.gen_Egl.map(bed["Egl"].chr)
        d["chr_otro"] = d.gen_otro.map(bed[b].chr)
        assert d.chr_Egl.notna().all() and d.chr_otro.notna().all(), "genes de anclas fuera de los BED"
        d["ini_Egl"] = d.gen_Egl.map(bed["Egl"].ini)
        d["ini_otro"] = d.gen_otro.map(bed[b].ini)
        d["Egl_es_NAC"] = d.gen_Egl.isin(nac["Egl"])
        d["otro_es_NAC"] = d.gen_otro.isin(nac[b])
        anc[b] = d
        ej = json.load(open(f"{JC}/Egl.{b}.ejecucion.json")) if os.path.exists(f"{JC}/Egl.{b}.ejecucion.json") else {}
        n_filt = sum(1 for _ in open(f"{JC}/Egl.{b}.last.filtered"))
        ct = pd.crosstab(d.chr_Egl, d.chr_otro)
        principal = {c: {k: round(v / f.sum() * 100, 1) for k, v in f.sort_values(ascending=False).head(2).items()}
                     for c, f in ct.iterrows()}
        prof = d.groupby("gen_Egl").gen_otro.nunique().value_counts().sort_index()
        prof_i = d.groupby("gen_otro").gen_Egl.nunique().value_counts().sort_index()
        genoma[b] = dict(
            especie=OTRAS[b], lineas_salida_LAST=ej.get("lineas_last"), aciertos_filtrados=n_filt,
            segundos_JCVI=ej.get("segundos"), pares_ancla=len(d), pares_rescatados_L=int(d.rescatado.sum()),
            bloques=int(d.bloque.nunique()), anclas_por_bloque_mediana=float(d.groupby("bloque").size().median()),
            bloques_con_20_o_mas=int((d.groupby("bloque").size() >= 20).sum()),
            genes_Egl=int(d.gen_Egl.nunique()), pct_genes_Egl=round(d.gen_Egl.nunique() / len(bed["Egl"]) * 100, 1),
            genes_otro=int(d.gen_otro.nunique()), pct_genes_otro=round(d.gen_otro.nunique() / len(bed[b]) * 100, 1),
            socios_por_gen_Egl={int(k): int(v) for k, v in prof.items()},
            socios_por_gen_otro={int(k): int(v) for k, v in prof_i.items()},
            cromosoma_principal_pct=principal)
        print(f"Egl-{b}:", {k: v for k, v in genoma[b].items() if k != "cromosoma_principal_pct"})

    # 3. Ortólogos colineales de los 189 NAC. El rescate (liftover) añade todos los aciertos LAST que caen dentro de un
    #    bloque, también los parálogos en tándem del vecindario. El ortólogo de cada gen es su mejor socio colineal
    #    (mayor puntaje LAST en bits; a igual puntaje, el ancla original antes que la rescatada); los demás quedan
    #    como socios colineales.
    mejor, inv = {}, {}
    for b in OTRAS:
        d = anc[b].sort_values(["puntaje", "rescatado"], ascending=[False, True])
        mejor[b] = d.drop_duplicates("gen_Egl").set_index("gen_Egl")
        inv[b] = d.drop_duplicates("gen_otro").set_index("gen_otro")
        d["es_mejor_socio_Egl"] = d.gen_otro.values == d.gen_Egl.map(mejor[b].gen_otro).values
        d["es_mejor_socio_otro"] = d.gen_Egl.values == d.gen_otro.map(inv[b].gen_Egl).values
        anc[b] = d
        m, dn = mejor[b], d[d.Egl_es_NAC]
        o = base.locus.map(m.gen_otro)
        base[f"ortologo_{b}"] = o
        base[f"ortologo_{b}_es_NAC"] = o.map(lambda x, s=nac[b]: (x in s) if isinstance(x, str) else np.nan)
        base[f"ortologo_{b}_puntaje"] = base.locus.map(m.puntaje)
        base[f"ortologo_{b}_rescatado"] = base.locus.map(m.rescatado)
        base[f"ortologo_{b}_reciproco"] = [(inv[b].gen_Egl.get(x) == l) if isinstance(x, str) else np.nan
                                           for l, x in zip(base.locus, o)]
        base[f"cromosoma_{b}"] = o.map(bed[b].chr)
        base[f"bloque_{b}"] = base.locus.map(m.bloque).astype("Int64")
        base[f"anclas_bloque_{b}"] = base.locus.map(m.anclas_en_bloque).astype("Int64")
        g = dn.groupby("gen_Egl")
        base[f"n_socios_{b}"] = base.locus.map(g.gen_otro.nunique()).fillna(0).astype(int)
        base[f"n_socios_{b}_NAC"] = base.locus.map(
            dn[dn.otro_es_NAC].groupby("gen_Egl").gen_otro.nunique()).fillna(0).astype(int)
        base[f"socios_{b}"] = base.locus.map(g.gen_otro.agg(unir))
    base["simbolo_EgrNAC_ortologo"] = base.ortologo_Egr.map(sim_de_egr)

    def contraste(r):
        if pd.isna(r.coortologo_EgrNAC):
            return "sin co-ortólogo en T01"
        if pd.isna(r.coortologo_v2) or r.coortologo_v2 not in bed["Egr"].index:
            return "co-ortólogo fuera de los cromosomas v2"
        if r.n_socios_Egr == 0:
            return "sin par colineal"
        if r.coortologo_v2 == r.ortologo_Egr:
            return "coincide con el ortólogo colineal"
        if r.coortologo_v2 in r.socios_Egr.split(";"):
            return "co-ortólogo entre los socios colineales"
        return "par colineal distinto"

    base["contraste_T01"] = base.apply(contraste, axis=1)
    ie = inv["Egr"]
    base["mejor_socio_Egl_del_coortologo"] = base.coortologo_v2.map(
        lambda e: np.nan if not isinstance(e, str) else (
            "sin par colineal" if e not in ie.index else sim_egl.get(ie.gen_Egl[e], ie.gen_Egl[e] + " (no NAC)")))
    orden = ["simbolo", "locus", "cromosoma", "inicio", "subfamilia", "modo_duplicacion", "regla", "coortologo_EgrNAC",
             "coortologo_v2", "contraste_T01", "mejor_socio_Egl_del_coortologo"]
    for b in OTRAS:
        orden += [f"ortologo_{b}", f"ortologo_{b}_es_NAC", f"ortologo_{b}_puntaje", f"ortologo_{b}_rescatado",
                  f"ortologo_{b}_reciproco", f"cromosoma_{b}", f"bloque_{b}", f"anclas_bloque_{b}", f"n_socios_{b}",
                  f"n_socios_{b}_NAC", f"socios_{b}"] + (["simbolo_EgrNAC_ortologo"] if b == "Egr" else [])
    orden.append("simbolo_manuscrito_2026")
    t48 = base[orden].sort_values(["cromosoma", "inicio"])
    t48.to_csv(T48, index=False)

    # Grupos de parálogos de E. globulus que comparten co-ortólogo en T01: ¿cuántos lo tienen como ortólogo colineal?
    gr = pd.DataFrame([dict(coortologo=e, miembros=len(x), con_el_como_ortologo=int((x.ortologo_Egr == e).sum()))
                       for e, x in base[base.coortologo_v2.notna()].groupby("coortologo_v2") if len(x) > 1],
                      columns=["coortologo", "miembros", "con_el_como_ortologo"])

    nac_res = {}
    for b in OTRAS:
        o, es = base[f"ortologo_{b}"], base[f"ortologo_{b}_es_NAC"]
        ib = inv[b][inv[b].otro_es_NAC]
        nac_res[b] = dict(
            NAC_hmm=len(nac[b]),
            NAC_Egl_con_ortologo=int(o.notna().sum()),
            NAC_Egl_ortologo_NAC=int((es == True).sum()),
            NAC_Egl_ortologo_no_NAC=base.loc[es == False, "simbolo"].tolist(),
            NAC_Egl_sin_par=int(o.isna().sum()),
            ortologo_reciproco=int((base[f"ortologo_{b}_reciproco"] == True).sum()),
            ortologo_rescatado=int((base[f"ortologo_{b}_rescatado"] == True).sum()),
            ortologos_distintos=int(o.nunique()),
            NAC_Egl_por_ortologo={int(k): int(v) for k, v in o.value_counts().value_counts().sort_index().items()},
            socios_por_NAC={int(k): int(v) for k, v in base[f"n_socios_{b}"].value_counts().sort_index().items()},
            NAC_otro_con_par_en_Egl=len(ib),
            NAC_otro_cuyo_mejor_socio_es_NAC=int(ib.Egl_es_NAC.sum()),
            NAC_otro_cuyo_mejor_socio_no_es_NAC=[f"{x}->{r.gen_Egl}" for x, r in ib[~ib.Egl_es_NAC].iterrows()],
            NAC_otro_sin_par_en_Egl=len(nac[b] - set(ib.index)),
            con_ortologo_por_modo={k: f"{int(v.notna().sum())}/{len(v)}" for k, v in o.groupby(base.modo_duplicacion)},
            con_ortologo_por_subfamilia={k: f"{int(v.notna().sum())}/{len(v)}" for k, v in o.groupby(base.subfamilia)})
        print(f"NAC Egl-{b}:", {k: (v if not isinstance(v, list) or len(v) <= 8 else f"{len(v)} elementos")
                                for k, v in nac_res[b].items() if not k.startswith("con_ortologo_por")})

    ct = pd.crosstab(base.regla, base.contraste_T01)
    print("contraste con T01 (regla x categoría):\n" + ct.to_string())

    # 4. Pares de anclas con un NAC en algún lado
    filas = []
    for b in OTRAS:
        d = anc[b]
        x = d[d.Egl_es_NAC | d.otro_es_NAC].copy()
        x.insert(0, "especie", b)
        x["simbolo_Egl"] = x.gen_Egl.map(sim_egl)
        x["simbolo_EgrNAC"] = x.gen_otro.map(sim_de_egr) if b == "Egr" else np.nan
        filas.append(x)
    t49 = pd.concat(filas)[["especie", "bloque", "anclas_en_bloque", "simbolo_Egl", "gen_Egl", "chr_Egl", "ini_Egl",
                            "Egl_es_NAC", "gen_otro", "chr_otro", "ini_otro", "otro_es_NAC", "simbolo_EgrNAC",
                            "es_mejor_socio_Egl", "es_mejor_socio_otro", "rescatado", "puntaje"]]
    t49.to_csv(T49, index=False)

    # 5. Contraste con el mapeo posicional anterior
    o = pd.read_csv(ANTIGUO)
    o["locus_v2"] = o.locus.str.replace("Eucgr.", "EUGRSUZ_", regex=False)
    o["chr_Egr_v2"] = o.locus_v2.map(bed["Egr"].chr)
    seq = t01.dropna(subset=["simbolo_manuscrito_2026"]).set_index("simbolo_manuscrito_2026")
    o["locus_Egl_por_secuencia"] = o.EglNAC.map(seq.locus)
    o["simbolo_189_por_secuencia"] = o.EglNAC.map(seq.simbolo)
    o["ortologo_colineal_Egl"] = o.locus_v2.map(ie.gen_Egl)
    o["simbolo_189_colineal"] = o.ortologo_colineal_Egl.map(
        lambda x: sim_egl.get(x, "no NAC") if isinstance(x, str) else np.nan)
    o["chr_Egl_colineal"] = o.ortologo_colineal_Egl.map(bed["Egl"].chr)
    o["socios_colineales_Egl"] = o.locus_v2.map(anc["Egr"].groupby("gen_otro").gen_Egl.agg(unir))
    o["cromosoma_antiguo_coincide"] = [("Egl%02d" % int(c) == x) if isinstance(x, str) and pd.notna(c) else np.nan
                                       for c, x in zip(o.chr_globulus, o.chr_Egl_colineal)]
    o["sintenia_coincide_con_secuencia"] = [(l == x) if isinstance(l, str) and isinstance(x, str) else np.nan
                                            for l, x in zip(o.locus_Egl_por_secuencia, o.ortologo_colineal_Egl)]
    o["secuencia_entre_socios_colineales"] = [(l in s.split(";")) if isinstance(l, str) and isinstance(s, str)
                                              else np.nan for l, s in zip(o.locus_Egl_por_secuencia,
                                                                          o.socios_colineales_Egl)]
    L_act, L_ant = largos(SR_ACT), largos(SR_ANT)
    nacb = bed["Egl"].loc[sorted(nac["Egl"])]
    dist = []
    for c, p in zip(o.chr_globulus, o.pos_globulus):
        n = nacb[nacb.chr == "Egl%02d" % int(c)]
        dd = np.where((n.ini <= p) & (p <= n.fin), 0, np.minimum((n.ini - p).abs(), (n.fin - p).abs()))
        dist.append(int(dd.min()) if len(dd) else np.nan)
    o["dist_al_NAC_mas_cercano_pb_ensamblado_actual"] = dist
    o.to_csv(T50, index=False)
    ev, en, cc = o.sintenia_coincide_con_secuencia, o.secuencia_entre_socios_colineales, o.cromosoma_antiguo_coincide
    colapsados = [dict(chr=int(c), pos=int(p), simbolos=";".join(x.EglNAC),
                       ortologos_colineales_distintos=int(x.ortologo_colineal_Egl.nunique()),
                       simbolos_189_colineales=";".join(x.simbolo_189_colineal.fillna("sin par")))
                  for (c, p), x in o[o.sitio_colapsado.astype(bool)].groupby(["chr_globulus", "pos_globulus"])]
    cambian = [dict(simbolo=r.EglNAC, chr_grandis=str(r.chr_grandis), chr_globulus_antiguo=int(r.chr_globulus),
                    chr_Egl_colineal=r.chr_Egl_colineal if isinstance(r.chr_Egl_colineal, str) else None,
                    simbolo_189_colineal=r.simbolo_189_colineal if isinstance(r.simbolo_189_colineal, str) else None,
                    simbolo_189_por_secuencia=r.simbolo_189_por_secuencia
                    if isinstance(r.simbolo_189_por_secuencia, str) else None)
               for r in o[~o.mismo_cromosoma.astype(bool)].itertuples()]
    d_act = o.dist_al_NAC_mas_cercano_pb_ensamblado_actual
    antiguo = dict(
        filas=len(o), en_cromosomas_Egr_v2=int(o.chr_Egr_v2.notna().sum()),
        con_ortologo_colineal=int(o.ortologo_colineal_Egl.notna().sum()),
        ortologo_colineal_es_NAC=int((o.simbolo_189_colineal.notna() & (o.simbolo_189_colineal != "no NAC")).sum()),
        con_locus_por_secuencia=int(o.locus_Egl_por_secuencia.notna().sum()),
        sintenia_coincide_con_secuencia=f"{int((ev == True).sum())} de {int(ev.notna().sum())}",
        secuencia_entre_socios_colineales=f"{int((en == True).sum())} de {int(en.notna().sum())}",
        cromosoma_antiguo_coincide=f"{int((cc == True).sum())} de {int(cc.notna().sum())}",
        discrepancias=registros(o.loc[ev == False, ["EglNAC", "simbolo_189_por_secuencia", "simbolo_189_colineal",
                                                    "secuencia_entre_socios_colineales"]]),
        cambian_de_cromosoma_en_el_mapeo_antiguo=cambian, sitios_colapsados=colapsados,
        coordenadas=dict(
            fuera_del_cromosoma_GCA_046226855=int((o.pos_globulus > o.chr_globulus.map(L_act)).sum()),
            fuera_del_cromosoma_GCA_014182545=int((o.pos_globulus > o.chr_globulus.map(L_ant)).sum()),
            dentro_de_un_NAC_de_T01=int((d_act == 0).sum()), a_menos_de_10_kb=int((d_act <= 1e4).sum()),
            mediana_dist_al_NAC_kb=round(float(np.nanmedian(d_act)) / 1e3, 1)))
    print("mapeo anterior:", {k: v for k, v in antiguo.items() if k not in (
        "discrepancias", "cambian_de_cromosoma_en_el_mapeo_antiguo", "sitios_colapsados")})

    resumen = dict(
        metodo=dict(versiones=versiones, cscore=CSCORE, tandem_Nmax=10, min_anclas=4, dist_genes=20,
                    archivo="lifted.anchors", alineador="LAST sobre CDS (lastal -u0 -i3G)",
                    ortologo="mejor socio colineal: mayor puntaje LAST (bits) en lifted.anchors; empate: ancla original",
                    reciproco="el gen de E. globulus es también el mejor socio colineal de su ortólogo"),
        NAC_por_especie={sp: len(v) for sp, v in nac.items()}, Hussey_S4_en_v2=hus_v2,
        genoma=genoma, NAC=nac_res,
        contraste_T01=dict(
            categorias=base.contraste_T01.value_counts().to_dict(),
            regla_x_categoria={r: {k: int(v) for k, v in f.items() if v} for r, f in ct.iterrows()},
            grupos_coortologo_compartido=dict(
                grupos=len(gr), genes=int(gr.miembros.sum()),
                con_un_miembro_ortologo=int((gr.con_el_como_ortologo == 1).sum()),
                con_ninguno=int((gr.con_el_como_ortologo == 0).sum()),
                con_varios=int((gr.con_el_como_ortologo > 1).sum())),
            sin_coortologo_T01=registros(base[base.contraste_T01 == "sin co-ortólogo en T01"][
                ["simbolo", "ortologo_Egr", "ortologo_Egr_es_NAC", "simbolo_EgrNAC_ortologo"]]),
            par_colineal_distinto=registros(base[base.contraste_T01 == "par colineal distinto"][
                ["simbolo", "coortologo_EgrNAC", "ortologo_Egr", "simbolo_EgrNAC_ortologo",
                 "mejor_socio_Egl_del_coortologo"]])),
        mapeo_posicional_anterior=antiguo,
        salidas=[corto(x) for x in (T48, T49, T50, RESUMEN)] + [corto(JC) + "/"],
        segundos=round(time.time() - t0))
    json.dump(resumen, open(RESUMEN, "w"), ensure_ascii=False, indent=1,
              default=lambda v: v.item() if hasattr(v, "item") else str(v))
    print(f"escritos: {corto(T48)} ({len(t48)} filas), {corto(T49)} ({len(t49)} filas), {corto(T50)} ({len(o)} filas), "
          f"{corto(RESUMEN)} | {time.time() - t0:.0f} s")

if __name__ == "__main__":
    main()
