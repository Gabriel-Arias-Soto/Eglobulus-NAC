#!/usr/bin/env python3
"""Genes NAC de dominio parcial: ¿modelo génico incompleto o gen truncado?

Lee los alineamientos HISAT2 contra el GENOMA COMPLETO (rnaseq_udec/alineamiento_parciales/{SRR}.regiones.bam,
filtrados a regiones_5kb.bed), usa solo alineamientos únicos (NH:i:1, primarios) y, por locus ±5 kb:
cobertura, uniones de empalme (anotadas o no), isoforma dominante y todas las isoformas compatibles con uniones
fiables, traducción y búsqueda del dominio NAM (Pfam PF02365) con hmmsearch. Escribe resultados/ y para si
alguna comprobación falla. No se teclea ninguna cifra de resultados: todo sale de los BAM, del genoma, de la
anotación y de T01/T16.
Uso: python analisis_parciales.py   (desde cualquier carpeta)
"""
import os, re, json, gzip, itertools, subprocess, sys
import numpy as np, pandas as pd

AQUI = os.path.dirname(os.path.abspath(__file__))
RAIZ = os.path.abspath(os.path.join(AQUI, "..", ".."))            # Revisión NAC
RNA = os.path.join(RAIZ, "rnaseq_udec")
OUT = os.path.join(AQUI, "resultados"); os.makedirs(OUT, exist_ok=True)
HOME = os.path.expanduser("~")
SAM = (os.environ.get("CONDA_ENVS", HOME + "/conda/envs") + "/rnaseq/bin/samtools")
HMMSEARCH = (os.environ.get("CONDA_ENVS", HOME + "/conda/envs") + "/hmmer/bin/hmmsearch")
HMM = os.path.join(RAIZ, "ANALISIS_Eglobulus", "referencia", "PF02365.hmm")
T01 = os.path.join(RAIZ, "ANALISIS_Eglobulus", "tablas", "T01_censo_maestro_189.csv")
T16 = os.path.join(RAIZ, "ANALISIS_Eglobulus", "tablas", "T16_arquitectura_dominios.csv")
GENOMA = os.path.join(RNA, "ref", "genoma.fa")
PROTEOMA = os.path.join(RNA, "ref", "GCA_046226855.1_ASM4622685v1_protein.faa.gz")
UNIONES_ANOT = os.path.join(RNA, "hisat2_idx", "uniones_anotadas.txt")

# --- diseño y reglas (declaradas aquí; ninguna se ajusta a posteriori) ---------------------------------------
LIBS = ["SRR7726479", "SRR7726470", "SRR7726471", "SRR9887586", "SRR7726474"]
PARCIALES = ["EglNAC37", "EglNAC87", "EglNAC195", "EglNAC197", "EglNAC214", "EglNAC249"]
CONTROL = "EglNAC66"
MIN_LECT_UNION, MIN_BIBL_UNION, MAX_INTRON = 3, 2, 20000      # unión fiable: canónica, ≥3 lecturas en ≥2 bibliotecas, intrón ≤20 kb
FRAC_DOMINANTE = 0.10                                          # isoforma dominante: uniones con ≥10 % de la más apoyada
UMBRAL_ISLA = 5                                                # isla de cobertura: ≥5 lecturas únicas (suma de bibliotecas)
MIN_COB_EXONES = 5.0                                           # "sin datos suficientes" si la cobertura media de los exones anotados es menor
COBERTURA_COMPLETA = 0.99                                      # criterio del censo: dominio completo si cubre ≥99 % del perfil (T16)
MIN_FRAG_AA = 30
CANON = {"GT-AG": "+", "GC-AG": "+", "AT-AC": "+", "CT-AC": "-", "CT-GC": "-", "GT-AT": "-"}

def sh(cmd): return subprocess.run(cmd, capture_output=True, text=True, check=True).stdout

# --- genoma por acceso aleatorio ------------------------------------------------------------------------------
fai = pd.read_csv(GENOMA + ".fai", sep="\t", header=None, names=["seq", "len", "off", "lb", "lB"]).set_index("seq")
def subseq(seq, s, e):
    r = fai.loc[seq]; o = lambda p: int(r.off + (p - 1) // r.lb * r.lB + (p - 1) % r.lb)
    with open(GENOMA, "rb") as fh:
        fh.seek(o(s)); return fh.read(o(e) - o(s) + 1).decode().replace("\n", "").upper()
RCT = str.maketrans("ACGTN", "TGCAN"); rc = lambda s: s.translate(RCT)[::-1]
B_, A_ = "TCAG", "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
COD = {a + b + c: A_[16 * i + 4 * j + k] for i, a in enumerate(B_) for j, b in enumerate(B_) for k, c in enumerate(B_)}
trad = lambda s: "".join(COD.get(s[i:i + 3], "X") for i in range(0, len(s) - 2, 3))
HMMLEN = int(re.search(r"^LENG\s+(\d+)", open(HMM).read(), re.M).group(1))

def hmmsearch(seqs, tag):
    fa = os.path.join(OUT, f"_{tag}.faa"); tb = os.path.join(OUT, f"_{tag}.tbl")
    with open(fa, "w") as fh: fh.write("".join(f">{k}\n{v}\n" for k, v in seqs.items()))
    sh([HMMSEARCH, "--domtblout", tb, "--noali", "-E", "10", "--domE", "10", HMM, fa])
    rows = [l.split() for l in open(tb) if not l.startswith("#")]
    d = pd.DataFrame([(r[0], int(r[2]), float(r[12]), float(r[13]), int(r[15]), int(r[16]), int(r[17]), int(r[18])) for r in rows],
                     columns=["consulta", "largo_aa", "iEvalue", "score", "hmm_ini", "hmm_fin", "ali_ini", "ali_fin"])
    d["cobertura_hmm"] = (d.hmm_fin - d.hmm_ini + 1) / HMMLEN
    os.remove(fa); os.remove(tb)
    return d

def fragmentos(s, prefijo):
    """Tramos sin codón de parada (≥MIN_FRAG_AA) en los 3 marcos sentido."""
    out = {}
    for f in range(3):
        for m in re.finditer(r"[^*]{%d,}" % MIN_FRAG_AA, trad(s[f:])): out[f"{prefijo}|f{f}|{m.start()}"] = m.group()
    return out

# --- entradas -------------------------------------------------------------------------------------------------
t01 = pd.read_csv(T01).set_index("simbolo"); t16 = pd.read_csv(T16).set_index("simbolo")
GENES = PARCIALES + [CONTROL]
reg = pd.read_csv(os.path.join(AQUI, "regiones_5kb.bed"), sep="\t", header=None, names=["seq", "a", "b", "gen"]).set_index("gen")
mod = pd.read_csv(os.path.join(AQUI, "modelos_anotados_7loci.tsv"), sep="\t")
ss = pd.read_csv(UNIONES_ANOT, sep="\t", header=None, names=["seq", "d", "a", "h"]); SS = set(zip(ss.seq, ss.d, ss.a))
assert set(reg.index) == set(GENES) and set(mod.gen) == set(GENES)

# totales por biblioteca (resúmenes de HISAT2 de las tandas; la prueba piloto no entra)
tot = {}
for L in LIBS:
    fs = sorted(f for f in os.listdir(os.path.join(AQUI, "tandas")) if f.startswith(L + "_t") and f.endswith(".resumen.txt"))
    tx = "".join(open(os.path.join(AQUI, "tandas", f)).read() for f in fs)
    tot[L] = dict(tandas=len(fs), lecturas=sum(int(x) for x in re.findall(r"Total reads: (\d+)", tx)),
                  unicas_genoma=sum(int(x) for x in re.findall(r"Aligned 1 time: (\d+)", tx)),
                  multiples_genoma=sum(int(x) for x in re.findall(r"Aligned >1 times: (\d+)", tx)))
diseno = pd.read_csv(os.path.join(RNA, "ena_diseno.tsv"), sep="\t").set_index("run_accession")
for L in LIBS: assert tot[L]["lecturas"] == int(diseno.loc[L, "read_count"]), (L, tot[L], diseno.loc[L, "read_count"])

CIG = re.compile(r"(\d+)([MIDNSHP=X])")
def bloques(pos, cigar):
    p = pos; bl = []; jn = []
    for n, op in CIG.findall(cigar):
        n = int(n)
        if op in "M=X": bl.append((p, p + n - 1)); p += n
        elif op == "D": p += n
        elif op == "N": jn.append((p - 1, p + n)); p += n          # última base exónica, primera del exón siguiente
    return bl, jn

def leer(L, gen, expr):
    r = reg.loc[gen]
    out = sh([SAM, "view", "-F", "0x904", "-e", expr, os.path.join(AQUI, f"{L}.regiones.bam"), f"{r.seq}:{r.a + 1}-{r.b}"])
    for ln in out.splitlines():
        f = ln.split("\t"); yield f[0], int(f[3]), f[5]

# --- cobertura, lecturas y uniones ----------------------------------------------------------------------------
filas_loc, filas_un, filas_ex = [], [], []; COB = {}
for gen in GENES:
    r = reg.loc[gen]; Lr = r.b - r.a; c = np.zeros((len(LIBS), Lr), dtype=np.int32)
    exa = mod[(mod.gen == gen) & (mod.tipo == "exon")].sort_values("ini")
    en_ex = lambda bl: any(s <= e2 and e >= s2 for s, e in bl for s2, e2 in zip(exa.ini, exa.fin))
    fila = dict(gen=gen)
    for k, L in enumerate(LIBS):
        n_reg = n_ex = 0
        for _, pos, cig in leer(L, gen, "[NH]==1"):
            bl, jn = bloques(pos, cig); n_reg += 1; n_ex += en_ex(bl)
            for s, e in bl:
                s0, e0 = max(s - 1 - r.a, 0), min(e - r.a, Lr)
                if e0 > s0: c[k, s0:e0] += 1
            for d, a in jn: filas_un.append((gen, L, r.seq, d, a))
        n_mul_ex = sum(en_ex(bloques(pos, cig)[0]) for _, pos, cig in leer(L, gen, "[NH]>1"))
        fila[f"unicas_region_{L}"] = n_reg; fila[f"unicas_exones_{L}"] = n_ex; fila[f"multiples_exones_{L}"] = n_mul_ex
    COB[gen] = c; csum = c.sum(0)
    for _, e in exa.iterrows():
        v = csum[e.ini - 1 - r.a:e.fin - r.a]
        filas_ex.append(dict(gen=gen, ini=e.ini, fin=e.fin, largo=e.fin - e.ini + 1, cobertura_media=round(float(v.mean()), 1), cobertura_min=int(v.min())))
    cob_ex = np.concatenate([csum[e.ini - 1 - r.a:e.fin - r.a] for _, e in exa.iterrows()])
    fila["cobertura_media_exones"] = round(float(cob_ex.mean()), 2)
    fila["unicas_exones"] = sum(fila[f"unicas_exones_{L}"] for L in LIBS)
    fila["multiples_exones"] = sum(fila[f"multiples_exones_{L}"] for L in LIBS)
    filas_loc.append(fila)
    np.save(os.path.join(OUT, f"cobertura_{gen}.npy"), c)

J = pd.DataFrame(filas_un, columns=["gen", "lib", "seq", "donante_fin", "aceptor_ini"])
JJ = J.groupby(["gen", "seq", "donante_fin", "aceptor_ini"]).agg(lecturas=("lib", "size"), bibliotecas=("lib", "nunique")).reset_index()
JJ["anotada"] = [(s, d - 1, a - 1) in SS for s, d, a in zip(JJ.seq, JJ.donante_fin, JJ.aceptor_ini)]
JJ["intron_pb"] = JJ.aceptor_ini - JJ.donante_fin - 1
JJ["motivo"] = [subseq(s, d + 1, d + 2) + "-" + subseq(s, a - 2, a - 1) for s, d, a in zip(JJ.seq, JJ.donante_fin, JJ.aceptor_ini)]
JJ["hebra_motivo"] = JJ.motivo.map(CANON).fillna("no canónico")
JJ["fiable"] = (JJ.lecturas >= MIN_LECT_UNION) & (JJ.bibliotecas >= MIN_BIBL_UNION) & (JJ.hebra_motivo != "no canónico") & (JJ.intron_pb <= MAX_INTRON)
JJ["hebra_gen"] = JJ.gen.map(lambda g: mod[mod.gen == g].hebra.iloc[0])
JJ["toca_modelo"] = [any(s <= x <= e for x in (d, a) for s, e in zip(mod[(mod.gen == g) & (mod.tipo == "exon")].ini, mod[(mod.gen == g) & (mod.tipo == "exon")].fin))
                     for g, d, a in zip(JJ.gen, JJ.donante_fin, JJ.aceptor_ini)]

def islas(gen, umbral):
    r = reg.loc[gen]; on = COB[gen].sum(0) >= umbral; out = []; i = 0
    while i < len(on):
        if on[i]:
            j = i
            while j < len(on) and on[j]: j += 1
            out.append((int(r.a + i + 1), int(r.a + j))); i = j
        else: i += 1
    return out
def isla_de(gen, p, lado, umbral=UMBRAL_ISLA):
    for s, e in islas(gen, umbral):
        if s <= p <= e: return s if lado == "ini" else e
    return None

def transcrito(gen, cj):
    """Exones a partir de una cadena de uniones: el primero y el último se extienden hasta el borde de su isla."""
    r = reg.loc[gen]; h = mod[mod.gen == gen].hebra.iloc[0]
    s0 = isla_de(gen, int(cj.donante_fin.iloc[0]), "ini", 3) or int(cj.donante_fin.iloc[0]) - 150
    e1 = isla_de(gen, int(cj.aceptor_ini.iloc[-1]), "fin", 3) or int(cj.aceptor_ini.iloc[-1]) + 150
    exs = [(s0, int(cj.donante_fin.iloc[0]))] + [(int(cj.aceptor_ini.iloc[k]), int(cj.donante_fin.iloc[k + 1])) for k in range(len(cj) - 1)] + [(int(cj.aceptor_ini.iloc[-1]), e1)]
    s = "".join(subseq(r.seq, a, b) for a, b in exs)
    return exs, (rc(s) if h == "-" else s)

def compatible(cj): return all(cj.donante_fin.iloc[k + 1] >= cj.aceptor_ini.iloc[k] + 20 for k in range(len(cj) - 1))

# --- proteínas: anotada (con comprobación contra el proteoma), isoforma dominante y todas las isoformas -------
prot_ref = {}; want = {t01.loc[g, "proteina"]: g for g in GENES}; cur = None
with gzip.open(PROTEOMA, "rt") as fh:
    for ln in fh:
        if ln.startswith(">"): cur = want.get(ln[1:].split()[0]); prot_ref.setdefault(cur, "") if cur else None
        elif cur: prot_ref[cur] += ln.strip()

consultas, iso_info, dom_info, CDS_INFO = {}, {}, {}, {}
for gen in GENES:
    r = reg.loc[gen]; h = mod[mod.gen == gen].hebra.iloc[0]
    cd = mod[(mod.gen == gen) & (mod.tipo == "CDS")].sort_values("ini")
    s = "".join(subseq(r.seq, a, b) for a, b in zip(cd.ini, cd.fin)); s = rc(s) if h == "-" else s
    pa = trad(s).rstrip("*")
    # NCBI traduce también un codón final incompleto (CDS parcial en 3'): se acepta un residuo de más solo en ese caso
    parcial3 = len(s) % 3 != 0
    assert pa == prot_ref[gen] or (parcial3 and prot_ref[gen][:len(pa)] == pa and len(prot_ref[gen]) - len(pa) == 1), \
        f"{gen}: la traducción del CDS anotado no coincide con {t01.loc[gen, 'proteina']}"
    fin3 = (cd.fin.max() + 1, cd.fin.max() + 3) if h == "+" else (cd.ini.min() - 3, cd.ini.min() - 1)
    cod_sig = subseq(r.seq, *fin3); cod_sig = rc(cod_sig) if h == "-" else cod_sig
    termina_en_stop = (trad(s[-3:]) == "*") or (not parcial3 and COD.get(cod_sig) == "*")
    CDS_INFO[gen] = dict(largo_cds=len(s), cds_multiplo_3=not parcial3, cds_termina_en_stop=bool(termina_en_stop))
    consultas[f"{gen}|anotada"] = pa
    j = JJ[(JJ.gen == gen) & JJ.fiable & (JJ.hebra_motivo == h)].sort_values("aceptor_ini").reset_index(drop=True)
    info = dict(n_uniones_fiables=len(j), dominante=None, n_isoformas=0)
    if len(j):
        jd = j[j.lecturas >= FRAC_DOMINANTE * j.lecturas.max()].reset_index(drop=True)
        n = len(jd); best = list(jd.lecturas); prev = [-1] * n            # cadena compatible de máximo apoyo
        for i in range(n):
            for k in range(i):
                if jd.donante_fin[i] >= jd.aceptor_ini[k] + 20 and best[k] + jd.lecturas[i] > best[i]: best[i] = best[k] + jd.lecturas[i]; prev[i] = k
        i = int(np.argmax(best)); ch = []
        while i != -1: ch.append(i); i = prev[i]
        cj = jd.loc[ch[::-1]]; exs, st = transcrito(gen, cj)
        info["dominante"] = dict(exones=exs, uniones=cj[["donante_fin", "aceptor_ini", "lecturas", "anotada"]].values.tolist())
        consultas.update(fragmentos(st, f"{gen}|dominante"))
        assert len(j) <= 16, f"{gen}: demasiadas uniones fiables para enumerar isoformas"
        for nn in range(1, len(j) + 1):
            for comb in itertools.combinations(range(len(j)), nn):
                cjj = j.loc[list(comb)]
                if not compatible(cjj): continue
                info["n_isoformas"] += 1
                _, st = transcrito(gen, cjj)
                consultas.update(fragmentos(st, f"{gen}|iso{'-'.join(map(str, comb))}"))
    iso_info[gen] = info
    # uniones observadas de cualquier apoyo (canónicas, misma hebra, ≤20 kb): ¿unen tramos del dominio no anotados?
    jo = JJ[(JJ.gen == gen) & (JJ.hebra_motivo == h) & (JJ.intron_pb <= MAX_INTRON)]
    for _, u in jo.iterrows():
        st = subseq(r.seq, int(u.donante_fin) - 600, int(u.donante_fin)) + subseq(r.seq, int(u.aceptor_ini), int(u.aceptor_ini) + 600)
        consultas.update(fragmentos(rc(st) if h == "-" else st, f"{gen}|union{int(u.donante_fin)}-{int(u.aceptor_ini)}"))
    # genoma ±5 kb en 6 fases (independiente del ARN)
    g = subseq(r.seq, r.a + 1, r.b)
    for hh, sq in (("+", g), ("-", rc(g))):
        for f in range(3):
            for m in re.finditer(r"[^*]{15,}", trad(sq[f:])): consultas[f"{gen}|genoma{hh}|f{f}|{m.start()}"] = m.group()

D = hmmsearch(consultas, "todo")
D["gen"] = D.consulta.str.split("|").str[0]; D["tipo"] = D.consulta.str.split("|").str[1].str.replace(r"\d.*$", "", regex=True)
D["con_M_antes"] = [("M" in consultas[q][:a - 1]) or consultas[q].startswith("M") for q, a in zip(D.consulta, D.ali_ini)]

# coordenadas genómicas de los tramos del dominio en 6 fases
def a_genoma(q, ai, af):
    gen, _, hf, st = q.split("|"); hh, f = hf[-1], int(st)  # placeholder (sobrescrito abajo)
gen6 = []
for _, x in D[D.tipo.str.startswith("genoma")].iterrows():
    gen, gh, fr, st = x.consulta.split("|"); hh = gh[-1]; f = int(fr[1:]); st = int(st); r = reg.loc[gen]; Lg = r.b - r.a
    nt0 = f + 3 * (st + x.ali_ini - 1); nt1 = f + 3 * (st + x.ali_fin) - 1
    g0, g1 = (r.a + 1 + nt0, r.a + 1 + nt1) if hh == "+" else (r.a + Lg - nt1, r.a + Lg - nt0)
    gen6.append(dict(gen=gen, hebra=hh, fase=f, hmm_ini=x.hmm_ini, hmm_fin=x.hmm_fin, score=x.score, iEvalue=x.iEvalue, g_ini=g0, g_fin=g1))
G6 = pd.DataFrame(gen6).query("score >= 5").sort_values(["gen", "g_ini"])

# --- veredictos -----------------------------------------------------------------------------------------------
LOC = pd.DataFrame(filas_loc).set_index("gen")
ver = []
for gen in GENES:
    d = D[D.gen == gen]; h = mod[mod.gen == gen].hebra.iloc[0]
    best = lambda tipo: d[d.tipo == tipo].cobertura_hmm.max() if (d.tipo == tipo).any() else 0.0
    cob_anot = best("anotada"); cob_dom = best("dominante")
    iso = d[d.tipo == "iso"]; cob_iso = iso.cobertura_hmm.max() if len(iso) else 0.0
    iso_ok = iso[(iso.cobertura_hmm >= COBERTURA_COMPLETA) & iso.con_M_antes]
    un = d[d.tipo == "union"]; un_ok = un[un.cobertura_hmm >= COBERTURA_COMPLETA]
    un_ok_info = []
    for q in un_ok.consulta.unique():
        dd, aa = map(int, q.split("|")[1].replace("union", "").split("-"))
        row = JJ[(JJ.gen == gen) & (JJ.donante_fin == dd) & (JJ.aceptor_ini == aa)].iloc[0]
        un_ok_info.append(dict(union=f"{dd}-{aa}", lecturas=int(row.lecturas), bibliotecas=int(row.bibliotecas), fiable=bool(row.fiable), toca_modelo=bool(row.toca_modelo)))
    cob_ex = LOC.loc[gen, "cobertura_media_exones"]
    if gen == CONTROL:
        v = "control: modelo confirmado" if (cob_anot >= COBERTURA_COMPLETA and cob_dom >= COBERTURA_COMPLETA) else "control: FALLA"
    elif cob_ex < MIN_COB_EXONES: v = "sin datos suficientes"
    elif len(iso_ok): v = "modelo incompleto"
    else: v = "compatible con truncado"
    ver.append(dict(gen=gen, locus=t01.loc[gen, "locus"], hebra=h, cobertura_T16=round(float(t16.loc[gen, "cobertura_hmm"]), 3),
                    cds_anotado_completo_3p=CDS_INFO[gen]['cds_multiplo_3'] and CDS_INFO[gen]['cds_termina_en_stop'],
                    cobertura_media_exones=cob_ex, unicas_exones=int(LOC.loc[gen, "unicas_exones"]), multiples_exones=int(LOC.loc[gen, "multiples_exones"]),
                    uniones_fiables=iso_info[gen]["n_uniones_fiables"], uniones_fiables_no_anotadas=int((JJ[(JJ.gen == gen) & JJ.fiable & ~JJ.anotada & JJ.toca_modelo]).shape[0]),
                    isoformas_probadas=iso_info[gen]["n_isoformas"], cob_dominio_anotada=round(float(cob_anot), 3),
                    cob_dominio_isoforma_dominante=round(float(cob_dom), 3), cob_dominio_mejor_isoforma=round(float(cob_iso), 3),
                    dominio_completo_via_union_no_fiable=json.dumps(un_ok_info, ensure_ascii=False) if un_ok_info and not len(iso_ok) else "",
                    tramos_genomicos=";".join(f"{a}-{b}" for a, b in zip(G6[G6.gen == gen].hmm_ini, G6[G6.gen == gen].hmm_fin)),
                    veredicto=v))
V = pd.DataFrame(ver)
assert V.loc[V.gen == CONTROL, "veredicto"].item() == "control: modelo confirmado", "el control positivo falla"

# --- salidas --------------------------------------------------------------------------------------------------
V.to_csv(os.path.join(OUT, "veredictos_parciales.csv"), index=False)
LOC.reset_index().to_csv(os.path.join(OUT, "lecturas_por_locus.csv"), index=False)
pd.DataFrame(filas_ex).to_csv(os.path.join(OUT, "cobertura_exones_anotados.csv"), index=False)
JJ.sort_values(["gen", "donante_fin"]).to_csv(os.path.join(OUT, "uniones_empalme.csv"), index=False)
D.drop(columns=[]).sort_values(["gen", "tipo", "cobertura_hmm"], ascending=[True, True, False]).to_csv(os.path.join(OUT, "dominio_hmmsearch.csv"), index=False)
G6.to_csv(os.path.join(OUT, "dominio_tramos_genomicos.csv"), index=False)
json.dump(dict(bibliotecas=tot, parametros=dict(MIN_LECT_UNION=MIN_LECT_UNION, MIN_BIBL_UNION=MIN_BIBL_UNION, MAX_INTRON=MAX_INTRON,
               FRAC_DOMINANTE=FRAC_DOMINANTE, UMBRAL_ISLA=UMBRAL_ISLA, MIN_COB_EXONES=MIN_COB_EXONES, COBERTURA_COMPLETA=COBERTURA_COMPLETA,
               filtro_unicos="NH:i:1, -F 0x904", HMM_LENG=HMMLEN), isoformas={g: iso_info[g] for g in GENES},
               versiones=dict(samtools=sh([SAM, "--version"]).splitlines()[0], hmmsearch=sh([HMMSEARCH, "-h"]).splitlines()[1])),
          open(os.path.join(OUT, "resumen_parciales.json"), "w"), ensure_ascii=False, indent=1, default=int)
print(V[["gen", "cobertura_media_exones", "unicas_exones", "uniones_fiables_no_anotadas", "cob_dominio_anotada",
         "cob_dominio_isoforma_dominante", "cob_dominio_mejor_isoforma", "veredicto"]].to_string(index=False))
