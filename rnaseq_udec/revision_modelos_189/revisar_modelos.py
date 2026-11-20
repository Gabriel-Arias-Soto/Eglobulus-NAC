# -*- coding: utf-8 -*-
"""Revisión de los modelos génicos de los 189 NAC con dos evidencias independientes (paso 2).

  - Proteínas: miniprot 0.18 alineó al genoma las proteínas NAC de E. grandis v2 y de A. thaliana (miniprot.gff,
    con la traducción de cada alineamiento). En cada locus se toma el mejor alineamiento de E. grandis (y, si no hay,
    de A. thaliana) que solapa la CDS anotada en la misma hebra; se comparan sus intrones con los de la CDS anotada,
    se buscan desplazamientos de marco y codones de parada internos, y se mide con hmmsearch qué fracción del perfil
    NAM (PF02365) cubre la proteína alineada.
  - RNA-seq: 5 bibliotecas de hoja completas alineadas con HISAT2 al genoma entero (SRR*.regiones.bam, filtrados a
    los 189 loci ±5 kb); solo alineamientos únicos y primarios. Uniones fiables con el mismo criterio que el análisis
    de los 6 genes de dominio parcial: canónicas, en la hebra del gen, ≥3 lecturas en ≥2 bibliotecas, intrón ≤20 kb.

Reglas de clasificación (declaradas aquí, antes de ver los resultados; el orden es la prioridad):
  1. «intrón no anotado en la CDS»: una unión fiable no anotada con ambos extremos dentro de la CDS y usada en al menos
     la mitad de las lecturas de la unión anotada más apoyada del gen (si es menor, se informa como empalme alternativo
     minoritario y no cambia el veredicto).
  2. «posible pseudogén o error de secuencia»: el mejor alineamiento de proteína tiene desplazamientos de marco o paradas
     internas y ningún otro alineamiento limpio cubre ≥80 % de su proteína.
  3. «modelo incompleto según proteínas»: un alineamiento limpio cubre ≥99 % del perfil NAM y la proteína anotada no (T16).
  4. «estructura distinta según proteínas»: los intrones del mejor alineamiento limpio no coinciden con los de la CDS
     anotada en el tramo que comparten; se indica si el RNA-seq apoya los intrones de la proteína o los anotados.
  5. «modelo respaldado»: el mejor alineamiento limpio (≥80 % de su proteína) reproduce los intrones anotados.
  6. «sin evidencia suficiente»: ningún alineamiento cubre ≥50 % de su proteína y el gen no tiene uniones fiables.
Controles: EglNAC66 (modelo confirmado) y EglNAC87 (intrón de 93 pb no anotado) deben salir con esos veredictos.

Parámetros y reglas:
  - miniprot se corre con -p 0.3 --outs=0.3 -N 100: con los valores por omisión cada proteína de E. grandis solo se
    alineaba a su mejor locus, y 18 genes de los grandes bloques en tándem (IVa, VII) quedaban sin ninguna proteína.
  - El desplazamiento de marco o la parada se ubica en la proteína (cadena cs) y se compara con el dominio NAM de esa
    proteína (hmmsearch --cut_ga): la regla 2 exige que caiga dentro del dominio y que la proteína anotada no tenga el
    dominio completo. Sin esa condición, genes muy expresados y con dominio completo (EglNAC19, 39 y 68) saldrían como
    pseudogenes por una diferencia con el ortólogo fuera del dominio; con ella quedan como «diferencia con el ortólogo
    fuera del dominio».
  - En la regla 4, el RNA-seq arbitra: si las uniones fiables confirman los intrones anotados y no los de la proteína,
    el veredicto es «modelo respaldado por RNA-seq (la proteína homóloga difiere)»; si confirman los de la proteína, es
    «estructura distinta (proteína y RNA-seq)»; si no hay uniones, «estructura distinta según proteínas, sin RNA-seq».
    Los intrones que difieren en ≤6 pb en un extremo se informan como «sitio de empalme desplazado».
  - La regla 1 exige además ≥10 lecturas en la unión no anotada.
  - Las proteínas se alinean a cada locus ±5 kb por separado (regiones.fa; miniprot -p 0.01 --outs=0.01 -N 500), no al
    genoma entero: alineadas al genoma, las proteínas de E. grandis solo se informaban en sus loci más parecidos, y 8
    genes de la subfamilia VII con el dominio completo quedaban sin evidencia aunque su proteína sí alinea a su locus.
  - Ajuste posterior a ver los resultados (se declara como tal): cuando el RNA-seq apoya a la vez un intrón anotado y uno
    de la proteína que lo reemplaza, el veredicto es «isoformas alternativas» y no «estructura distinta».
  - Con proteínas de E. grandis, varios genes de VII no tienen evidencia de proteína: sus co-ortólogos alinean a su locus
    con ~35 % de identidad (EglNAC53 frente a EUGRSUZ_E00541), por debajo de lo que miniprot informa; la proteína de
    E. globulus del propio gen sí alinea a su locus (control de que la región y el método funcionan).
"""
import bisect, collections, csv, glob, json, os, re, subprocess
import pandas as pd

AQUI = os.path.dirname(os.path.abspath(__file__))
RNA = os.path.dirname(AQUI)
RAIZ = os.path.dirname(RNA)
A = os.path.join(RAIZ, "ANALISIS_Eglobulus")
OUT = os.path.join(AQUI, "resultados"); os.makedirs(OUT, exist_ok=True)
HOME = os.path.expanduser("~")
SAM = (os.environ.get("CONDA_ENVS", HOME + "/conda/envs") + "/rnaseq/bin/samtools")
HMMSEARCH = (os.environ.get("CONDA_ENVS", HOME + "/conda/envs") + "/hmmer/bin/hmmsearch")
HMM = os.path.join(A, "referencia", "PF02365.hmm")
GENOMA = os.path.join(RNA, "ref", "genoma.fa")
LIBS = ["SRR7726470", "SRR7726471", "SRR7726474", "SRR7726479", "SRR9887586"]
MIN_LECT, MIN_BIBL, MAX_INTRON = 3, 2, 20000
FRAC_DOMINANTE, MIN_LECT_DOMINANTE, DESPLAZ_MAX = 0.5, 10, 6
QCOV_LIMPIO, QCOV_MIN, COMPLETO = 0.80, 0.50, 0.99
CANON = {"GT-AG": "+", "GC-AG": "+", "AT-AC": "+", "CT-AC": "-", "CT-GC": "-", "GT-AT": "-"}
sh = lambda cmd: subprocess.run(cmd, capture_output=True, text=True, check=True).stdout

# --- genoma -------------------------------------------------------------------------------------------------
fai = {l.split("\t")[0]: list(map(int, l.split("\t")[1:5])) for l in open(GENOMA + ".fai")}


def subseq(seq, s, e):
    ln, off, lb, lB = fai[seq]
    o = lambda p: off + (p - 1) // lb * lB + (p - 1) % lb
    with open(GENOMA, "rb") as fh:
        fh.seek(o(s)); return fh.read(o(e) - o(s) + 1).decode().replace("\n", "").upper()


# --- entradas -----------------------------------------------------------------------------------------------
T = lambda c: pd.read_csv(glob.glob(os.path.join(A, "tablas", f"{c}_*.csv"))[0])
t01 = T("T01").set_index("simbolo"); t16 = T("T16").set_index("simbolo")
mod = pd.read_csv(os.path.join(AQUI, "modelos_anotados_189.tsv"), sep="\t")
reg = pd.read_csv(os.path.join(AQUI, "regiones_5kb.bed"), sep="\t", header=None, names=["seq", "a", "b", "gen"]).set_index("gen")
ss = pd.read_csv(os.path.join(RNA, "hisat2_idx", "uniones_anotadas.txt"), sep="\t", header=None, names=["seq", "d", "a", "h"])
SS = set(zip(ss.seq, ss.d + 1, ss.a + 1))           # (seq, última base exónica, primera base del exón siguiente), 1-based
GENES = sorted(t01.index, key=lambda g: (reg.loc[g, "seq"], reg.loc[g, "a"]))
assert len(GENES) == 189

INFO = {}
for g in GENES:
    m = mod[mod.gen == g]
    cds = m[m.tipo == "CDS"].sort_values("ini"); ex = m[m.tipo == "exon"].sort_values("ini")
    h = cds.hebra.iloc[0]; s = cds.seq.iloc[0]
    INFO[g] = dict(seq=s, hebra=h, cds=list(zip(cds.ini, cds.fin)), cds_ini=int(cds.ini.min()), cds_fin=int(cds.fin.max()),
                   tx_ini=int(ex.ini.min()) if len(ex) else int(cds.ini.min()), tx_fin=int(ex.fin.max()) if len(ex) else int(cds.fin.max()),
                   intr_cds={(a_[1], b_[0]) for a_, b_ in zip(list(zip(cds.ini, cds.fin))[:-1], list(zip(cds.ini, cds.fin))[1:])},
                   intr_tx={(a_[1], b_[0]) for a_, b_ in zip(list(zip(ex.ini, ex.fin))[:-1], list(zip(ex.ini, ex.fin))[1:])},
                   exones=list(zip(ex.ini, ex.fin)))

# --- miniprot -------------------------------------------------------------------------------------------------
REC, cur = [], None
for l in open(os.path.join(AQUI, "miniprot.gff")):
    if l.startswith("##PAF"):
        f = l.rstrip("\n").split("\t")
        tags = dict(x.split(":", 2)[0::2] for x in f[13:])
        cur = dict(prot=f[1], qlen=int(f[2]), qs=int(f[3]), qe=int(f[4]), hebra=f[5], seq=f[6], ts=int(f[8]) + 1, te=int(f[9]),
                   ms=int(tags["ms"]), fs=int(tags["fs"]), st=int(tags["st"]), cs=tags.get("cs", ""), cg=tags.get("cg", ""), cds=[]); REC.append(cur)
    elif l.startswith("##STA"):
        cur["sta"] = l.rstrip("\n").split("\t")[1]
    elif not l.startswith("#"):
        f = l.split("\t")
        if f[2] == "mRNA":
            cur["id"] = re.search(r"ID=([^;]+)", f[8]).group(1); cur["identidad"] = float(re.search(r"Identity=([\d.]+)", f[8]).group(1))
        elif f[2] == "CDS":
            cur["cds"].append((int(f[3]), int(f[4])))
CG = re.compile(r"(\d+)([MIDNUVFG])")
CS_SUB = re.compile(r":(\d+)|\*([acgtn]{3})([A-Z*])|\*([acgtn]{1,2})([A-Z*])|\+([A-Z*]+)|-([acgtn]+)|~([acgtn]{2})(\d+)([acgtn]{2})")
STOP = {"taa", "tag", "tga"}


def eventos(r):
    """Posiciones (en la proteína) de los desplazamientos de marco (operaciones F/G de la cadena cg) y de las paradas
    internas (codones de parada en las sustituciones de la cadena cs). Los codones partidos por un intrón (U/V en cg;
    «*g…» y «-aa» en cs) no son desplazamientos de marco."""
    ev, q = [], r["qs"]
    for n, op in CG.findall(r["cg"]):
        n = int(n)
        if op in "MI":
            q += n
        elif op in "FG":
            ev.append(("marco", q))
    q = r["qs"]
    for m in CS_SUB.finditer(r["cs"]):
        if m.group(1):
            q += int(m.group(1))
        elif m.group(2):
            if m.group(2) in STOP:
                ev.append(("parada", q))
            q += 1
        elif m.group(4):
            q += 1
        elif m.group(6):
            q += len(m.group(6))
    assert sum(t == "marco" for t, _ in ev) == r["fs"], (r["id"], ev, r["fs"])
    return sorted(ev, key=lambda e: e[1])


DOMQ = collections.defaultdict(list)
for f_ in ("Egr.dom", "Ath.dom"):
    sp = f_.split(".")[0]
    for l in open(os.path.join(AQUI, f_)):
        if not l.startswith("#"):
            x = l.split(); DOMQ[f"{sp}|{x[0]}"].append((int(x[17]) - 1, int(x[18])))   # ali_from-1, ali_to (0-based, semiabierto)
for r in REC:                                       # coordenadas locales (región ±5 kb) -> genoma
    if ":" in r["seq"]:
        r["region"] = r["seq"]; sq, rg = r["seq"].rsplit(":", 1); off = int(rg.split("-")[0]) - 1
        r["seq"], r["ts"], r["te"] = sq, r["ts"] + off, r["te"] + off
        r["cds"] = [(a_ + off, b_ + off) for a_, b_ in r["cds"]]
for r in REC:
    r["cds"].sort(); r["qcov"] = (r["qe"] - r["qs"]) / r["qlen"]
    r["eventos"] = eventos(r)
    r["eventos_dominio"] = [e for e in r["eventos"] if any(a_ <= e[1] < b_ for a_, b_ in DOMQ.get(r["prot"], []))]
    r["intr"] = {(a_[1], b_[0]) for a_, b_ in zip(r["cds"][:-1], r["cds"][1:])}
    r["especie"] = r["prot"].split("|")[0]
# dominio NAM de cada proteína alineada (traducción del genoma; las paradas internas se cambian por X)
fa = os.path.join(OUT, "_mp.faa")
with open(fa, "w") as fh:
    fh.write("".join(f">{r['id']}\n{r['sta'].replace('*', 'X')}\n" for r in REC))
tb = os.path.join(OUT, "_mp.tbl")
sh([HMMSEARCH, "--domtblout", tb, "--noali", "-E", "10", "--domE", "10", "--cpu", "8", HMM, fa])
HMMLEN = int(re.search(r"^LENG\s+(\d+)", open(HMM).read(), re.M).group(1))
cobmp = collections.defaultdict(list)
for l in open(tb):
    if not l.startswith("#"):
        x = l.split(); cobmp[x[0]].append((int(x[15]), int(x[16])))
os.remove(fa); os.remove(tb)


def cobertura(intervalos):
    """Fracción del perfil cubierta por la unión de los tramos del perfil (hmm_ini, hmm_fin)."""
    pos = set()
    for a_, b_ in intervalos:
        pos.update(range(a_, b_ + 1))
    return len(pos) / HMMLEN


for r in REC:
    r["cob_NAM"] = cobertura(cobmp.get(r["id"], []))

# --- RNA-seq: uniones y lecturas por gen ----------------------------------------------------------------------
CIG = re.compile(r"(\d+)([MIDNSHP=X])")
POR_SEQ = collections.defaultdict(list)
for g in GENES:
    POR_SEQ[reg.loc[g, "seq"]].append((int(reg.loc[g, "a"]), int(reg.loc[g, "b"]), g))
for s in POR_SEQ:
    POR_SEQ[s].sort()
INI = {s: [x[0] for x in v] for s, v in POR_SEQ.items()}
UN = collections.Counter(); UN_LIB = collections.defaultdict(set)
LECT = collections.Counter()
for L in LIBS:
    out = subprocess.Popen([SAM, "view", "-F", "0x904", "-e", "[NH]==1", os.path.join(AQUI, f"{L}.regiones.bam")], stdout=subprocess.PIPE, text=True)
    for ln in out.stdout:
        f = ln.split("\t", 6)
        s, p = f[2], int(f[3])
        bl, jn = [], []
        for n, op in CIG.findall(f[5]):
            n = int(n)
            if op in "M=X":
                bl.append((p, p + n - 1)); p += n
            elif op == "D":
                p += n
            elif op == "N":
                jn.append((p - 1, p + n)); p += n
        for d, a_ in jn:
            UN[(s, d, a_)] += 1; UN_LIB[(s, d, a_)].add(L)
        if s in INI:
            k = bisect.bisect_right(INI[s], bl[-1][1])
            for a0, b0, g in POR_SEQ[s][max(0, k - 6):k]:
                if b0 < bl[0][0]:
                    continue
                if any(x <= e and y >= s_ for x, y in bl for s_, e in INFO[g]["exones"]):
                    LECT[g] += 1
    out.wait(); assert out.returncode == 0, L
J = pd.DataFrame([(s, d, a_, n, len(UN_LIB[(s, d, a_)])) for (s, d, a_), n in UN.items()],
                 columns=["seq", "donante_fin", "aceptor_ini", "lecturas", "bibliotecas"])
J["intron_pb"] = J.aceptor_ini - J.donante_fin - 1
J["anotada"] = [(s, d, a_) in SS for s, d, a_ in zip(J.seq, J.donante_fin, J.aceptor_ini)]
J = J[(J.lecturas >= 2) | J.anotada].copy()
J["motivo"] = [subseq(s, d + 1, d + 2) + "-" + subseq(s, a_ - 2, a_ - 1) for s, d, a_ in zip(J.seq, J.donante_fin, J.aceptor_ini)]
J["hebra_motivo"] = J.motivo.map(CANON).fillna("no canónico")
J["fiable"] = (J.lecturas >= MIN_LECT) & (J.bibliotecas >= MIN_BIBL) & (J.hebra_motivo != "no canónico") & (J.intron_pb <= MAX_INTRON)

# --- veredicto por gen ----------------------------------------------------------------------------------------
filas, uniones_nuevas = [], []
for g in GENES:
    I = INFO[g]
    # proteínas
    region = f"{reg.loc[g, 'seq']}:{int(reg.loc[g, 'a']) + 1}-{int(reg.loc[g, 'b'])}"
    hits = [r for r in REC if r.get("region") == region and r["hebra"] == I["hebra"]
            and min(r["te"], I["cds_fin"]) - max(r["ts"], I["cds_ini"]) + 1 >= 0.5 * min(r["te"] - r["ts"] + 1, I["cds_fin"] - I["cds_ini"] + 1)]
    egr = [r for r in hits if r["especie"] == "Egr"] or hits
    mejor = max(egr, key=lambda r: r["ms"]) if egr else None
    limpios = [r for r in egr if r["fs"] == 0 and r["st"] == 0 and r["qcov"] >= QCOV_LIMPIO]
    mejor_l = max(limpios, key=lambda r: r["ms"]) if limpios else None
    dom_l = max((r["cob_NAM"] for r in limpios), default=0.0)
    cob_anot = float(t16.loc[g, "cobertura_hmm"])
    # intrones de la proteína frente a los anotados, en el tramo que comparten
    def comparar(r):
        a0, b0 = max(r["ts"], I["cds_ini"]), min(r["te"], I["cds_fin"])
        dentro = lambda iv: a0 <= iv[0] and iv[1] <= b0
        an = {iv for iv in I["intr_cds"] if dentro(iv)}; mp = {iv for iv in r["intr"] if dentro(iv)}
        return an, mp
    if mejor_l:
        an, mp = comparar(mejor_l); igual = an == mp
    else:
        an, mp, igual = set(), set(), None
    # RNA-seq
    jg = J[(J.seq == I["seq"]) & (J.donante_fin >= I["tx_ini"] - 5000) & (J.aceptor_ini <= I["tx_fin"] + 5000)]
    jf = jg[jg.fiable & (jg.hebra_motivo == I["hebra"])]
    k = lambda d, a_: (d, a_)
    anot_cds = jf[[k(d, a_) in I["intr_cds"] for d, a_ in zip(jf.donante_fin, jf.aceptor_ini)]]
    sub = jg[[k(d, a_) in I["intr_tx"] for d, a_ in zip(jg.donante_fin, jg.aceptor_ini)]] if len(jg) else jg
    max_anot = int(sub.lecturas.max()) if len(sub) else 0
    nuevas = jf[~jf.anotada & (jf.donante_fin >= I["cds_ini"]) & (jf.aceptor_ini <= I["cds_fin"])]
    dominantes = nuevas[(nuevas.lecturas >= FRAC_DOMINANTE * max(max_anot, 1)) & (nuevas.lecturas >= MIN_LECT_DOMINANTE)]
    for _, u in nuevas.iterrows():
        uniones_nuevas.append(dict(gen=g, seq=u.seq, donante_fin=u.donante_fin, aceptor_ini=u.aceptor_ini, intron_pb=u.intron_pb,
                                   lecturas=u.lecturas, bibliotecas=u.bibliotecas, uso_relativo=round(u.lecturas / max(max_anot, 1), 2),
                                   en_intrones_proteina=bool(mejor_l and (u.donante_fin, u.aceptor_ini) in mejor_l["intr"])))
    apoyo = lambda ivs: sum(1 for iv in ivs if ((jf.donante_fin == iv[0]) & (jf.aceptor_ini == iv[1])).any())
    # intrones que difieren solo en un desplazamiento pequeño de un extremo
    solo_an, solo_mp = an - mp, mp - an
    desplazados = sum(1 for x in solo_an if any((x[0] == y[0] and abs(x[1] - y[1]) <= DESPLAZ_MAX) or
                                                 (x[1] == y[1] and abs(x[0] - y[0]) <= DESPLAZ_MAX) for y in solo_mp))
    rna_prot, rna_anot = (apoyo(solo_mp), apoyo(solo_an)) if mejor_l else (0, 0)
    ev_dom = mejor["eventos_dominio"] if mejor else []
    # clasificación
    if len(dominantes):
        v = "intrón no anotado en la CDS"
    elif mejor and ev_dom and not limpios and cob_anot < COMPLETO:
        v = "posible pseudogén o error de secuencia"
    elif dom_l >= COMPLETO and cob_anot < COMPLETO:
        v = "modelo incompleto según proteínas"
    elif mejor and (mejor["fs"] or mejor["st"]) and not limpios:
        v = "diferencia con el ortólogo" + ("" if ev_dom else " fuera del dominio")
    elif mejor_l and not igual:
        if rna_prot and rna_anot:
            v = "isoformas alternativas (RNA-seq apoya la anotada y la de la proteína)"
        elif rna_prot:
            v = "estructura distinta (proteína y RNA-seq)"
        elif rna_anot:
            v = "modelo respaldado por RNA-seq (la proteína homóloga difiere)"
        elif desplazados == len(solo_an) == len(solo_mp):
            v = "sitio de empalme desplazado"
        else:
            v = "estructura distinta según proteínas, sin RNA-seq"
    elif mejor_l and igual:
        v = "modelo respaldado"
    elif (not mejor or mejor["qcov"] < QCOV_MIN) and not len(jf):
        v = "sin evidencia suficiente"
    else:
        v = "evidencia parcial"
    filas.append(dict(
        gen=g, locus=t01.loc[g, "locus"], cromosoma=t01.loc[g, "cromosoma"], hebra=I["hebra"], subfamilia=t01.loc[g, "subfamilia_final"],
        cobertura_NAM_anotada=round(cob_anot, 3), intrones_CDS_anotados=len(I["intr_cds"]),
        prot_mejor=mejor["prot"] if mejor else "", prot_identidad=mejor["identidad"] if mejor else None,
        prot_cobertura=round(mejor["qcov"], 3) if mejor else None, prot_desplazamientos=mejor["fs"] if mejor else None,
        prot_paradas=mejor["st"] if mejor else None, prot_limpia=mejor_l["prot"] if mejor_l else "",
        cobertura_NAM_prot_limpia=round(dom_l, 3), intrones_comparados_anot=len(an), intrones_comparados_prot=len(mp),
        intrones_coinciden=len(an & mp), intrones_iguales=igual, intrones_desplazados=desplazados,
        eventos_marco_parada=";".join(f"{t}@{q}" for t, q in mejor["eventos"]) if mejor else "",
        eventos_en_dominio=len(ev_dom),
        intrones_prot_con_RNA=apoyo(mp - an) if mejor_l else 0, intrones_anot_sin_prot_con_RNA=apoyo(an - mp) if mejor_l else 0,
        lecturas_exones=LECT[g], uniones_fiables=len(jf), intrones_CDS_confirmados_RNA=len(anot_cds),
        uniones_nuevas_CDS=len(nuevas), uniones_nuevas_dominantes=len(dominantes), veredicto=v))
V = pd.DataFrame(filas)
V.to_csv(os.path.join(OUT, "revision_modelos_189.csv"), index=False)
pd.DataFrame(uniones_nuevas).to_csv(os.path.join(OUT, "uniones_no_anotadas_CDS.csv"), index=False)
J.to_csv(os.path.join(OUT, "uniones_189.csv"), index=False)
# controles
c66, c87 = V.set_index("gen").loc["EglNAC66"], V.set_index("gen").loc["EglNAC87"]
assert c66.veredicto == "modelo respaldado", c66.to_dict()
assert c87.veredicto == "intrón no anotado en la CDS", c87.to_dict()
res = dict(libs=LIBS, proteinas=len({r["prot"] for r in REC}), alineamientos=len(REC), veredictos=V.veredicto.value_counts().to_dict(),
           parametros=dict(MIN_LECT=MIN_LECT, MIN_BIBL=MIN_BIBL, MAX_INTRON=MAX_INTRON, FRAC_DOMINANTE=FRAC_DOMINANTE,
                           QCOV_LIMPIO=QCOV_LIMPIO, QCOV_MIN=QCOV_MIN, COMPLETO=COMPLETO),
           versiones=dict(miniprot=sh(["miniprot", "--version"]).strip(),
                          samtools=sh([SAM, "--version"]).split("\n")[0].split()[-1]))
json.dump(res, open(os.path.join(OUT, "resumen_revision.json"), "w", encoding="utf-8"), ensure_ascii=False, indent=1)
print(json.dumps(res["veredictos"], ensure_ascii=False))
