#!/usr/bin/env python
"""Revisión del árbol de la familia NAC: EglNAC37 y EglNAC87 entran con su modelo corregido.

Pasos (los del Methods S1, «Identification of NAC genes» y «Review of the gene models»):
  1. Proteínas corregidas.
     EglNAC87: CDS anotada sin el intrón no anotado que respaldan las lecturas de RNA-seq
               (rnaseq_udec/revision_modelos_189/resultados/uniones_no_anotadas_CDS.csv).
     EglNAC37: traducción del genoma a lo largo del alineamiento limpio (sin desplazamientos ni paradas) de la proteína
               de E. grandis que eligió la revisión de modelos (columna prot_limpia; línea ##STA de miniprot.gff).
  2. Dominio NAM de cada proteína corregida con hmmsearch --cut_ga (envolvente del dominio, como en el resto).
  3. Los dos dominios se añaden al alineamiento recortado de los otros 339 con MAFFT --add --keeplength.
     El árbol se infiere después con IQ-TREE (ver construir_arbol_185.sh).
  4. Lecturas de hoja sobre los desplazamientos de marco de los posibles pseudogenes con lecturas (EglNAC193, EglNAC196).
  5. Escaneo de los siete motivos MEME de los subdominios A-E sobre las proteínas corregidas, con matrices de pesos
     construidas con los sitios MEME de T05 (umbral = percentil 5 de las puntuaciones de esos sitios).
Salidas en arbol_185/: corregidas.faa, corregidas.dom, nuevos_dom.fa, NAC_dom.trim_341.fa,
  pseudogenes_lecturas/resumen_pseudogenes_lecturas.json, motivos_modelos_corregidos.json
Uso: NAC_ROOT=<carpeta del proyecto> python revision_arbol_185.py <alineamiento recortado de 339 secuencias>
Requiere en el PATH: hmmsearch (HMMER 3.4), mafft (7.526), miniprot (0.18, solo para el paso 4), samtools (1.24).
"""
import collections, itertools, json, os, re, shlex, subprocess, sys
import numpy as np, pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

RAIZ = os.environ.get("NAC_ROOT", ".")
A = os.path.join(RAIZ, "ANALISIS_Eglobulus"); TB = os.path.join(A, "tablas")
RM = os.path.join(RAIZ, "rnaseq_udec", "revision_modelos_189")
GEN = os.path.join(RAIZ, "rnaseq_udec", "ref", "genoma.fa")
HMM = os.path.join(A, "referencia", "PF02365.hmm")
OUT = os.environ.get("NAC_OUT", os.path.join(RAIZ, "arbol_185")); os.makedirs(OUT, exist_ok=True)
sh = lambda cmd, **k: subprocess.run(cmd, capture_output=True, text=True, check=True, **k).stdout

fai = {l.split("\t")[0]: list(map(int, l.split("\t")[1:5])) for l in open(GEN + ".fai")}


def subseq(seq, s, e):
    ln, off, lb, lB = fai[seq]
    o = lambda p: off + (p - 1) // lb * lB + (p - 1) % lb
    with open(GEN, "rb") as fh:
        fh.seek(o(s)); return fh.read(o(e) - o(s) + 1).decode().replace("\n", "").upper()


def traducir(exones, seq, hebra):
    s = "".join(subseq(seq, a, b) for a, b in sorted(exones))
    s = str(Seq(s).reverse_complement()) if hebra == "-" else s
    return str(Seq(s[: len(s) // 3 * 3]).translate()).rstrip("*")


rev = pd.read_csv(os.path.join(RM, "resultados", "revision_modelos_189.csv"))
mod = pd.read_csv(os.path.join(RM, "modelos_anotados_189.tsv"), sep="\t")
reg = pd.read_csv(os.path.join(RM, "regiones_5kb.bed"), sep="\t", header=None, names=["seq", "a", "b", "gen"]).set_index("gen")
t01 = pd.read_csv(os.path.join(TB, "T01_censo_maestro_189.csv")).set_index("simbolo")

# --- miniprot ------------------------------------------------------------------------------------------------
REC, cur = [], None
for l in open(os.path.join(RM, "miniprot.gff")):
    if l.startswith("##PAF"):
        f = l.rstrip("\n").split("\t"); tg = dict(x.split(":", 2)[0::2] for x in f[13:])
        cur = dict(prot=f[1], qlen=int(f[2]), qs=int(f[3]), qe=int(f[4]), seq=f[6], ts=int(f[8]) + 1, te=int(f[9]),
                   fs=int(tg["fs"]), st=int(tg["st"]), cg=tg.get("cg", ""), cds=[]); REC.append(cur)
    elif l.startswith("##STA"):
        cur["sta"] = l.rstrip("\n").split("\t")[1]
    elif not l.startswith("#") and l.split("\t")[2] == "CDS":
        f = l.split("\t"); cur["cds"].append((int(f[3]), int(f[4])))


def registros(g):
    s, a = reg.loc[g, "seq"], int(reg.loc[g, "a"])
    return [r for r in REC if ":" in r["seq"] and r["seq"].rsplit(":", 1)[0] == s
            and abs(int(r["seq"].rsplit(":", 1)[1].split("-")[0]) - 1 - a) <= 1]


# --- 1. proteínas corregidas ------------------------------------------------------------------------------------
cds = lambda g: mod[(mod.gen == g) & (mod.tipo == "CDS")].sort_values("ini")
un = pd.read_csv(os.path.join(RM, "resultados", "uniones_no_anotadas_CDS.csv"))
u87 = un[un.gen == "EglNAC87"].sort_values("lecturas").iloc[-1]
c87 = cds("EglNAC87"); ex = []
for a_, b_ in zip(c87.ini, c87.fin):
    if a_ <= u87.donante_fin < b_ and a_ < u87.aceptor_ini <= b_:
        ex += [(a_, int(u87.donante_fin)), (int(u87.aceptor_ini), b_)]
    else:
        ex.append((a_, b_))
p87 = traducir(ex, c87.seq.iloc[0], c87.hebra.iloc[0])
assert "*" not in p87 and p87.startswith("M")
pl = rev.loc[rev.gen == "EglNAC37", "prot_limpia"].iloc[0]
r37 = sorted([r for r in registros("EglNAC37") if r["prot"] == pl and r["fs"] == 0 and r["st"] == 0],
             key=lambda r: r["qe"] - r["qs"])[-1]
p37 = r37["sta"].rstrip("*")
assert "*" not in p37 and p37.startswith("M")
CORR = {"EglNAC37": p37, "EglNAC87": p87}
with open(os.path.join(OUT, "corregidas.faa"), "w") as fh:
    fh.write("".join(f">{g}\n{p}\n" for g, p in CORR.items()))

# --- 2. dominio NAM (envolvente de hmmsearch --cut_ga) -------------------------------------------------------
dom = os.path.join(OUT, "corregidas.dom")
sh(["hmmsearch", "--cut_ga", "--domtblout", dom, "--noali", HMM, os.path.join(OUT, "corregidas.faa")])
D = {}
for l in open(dom):
    if not l.startswith("#"):
        x = l.split(); D.setdefault(x[0], (int(x[19]), int(x[20])))       # primer dominio: env_from, env_to
assert set(D) == set(CORR), D
with open(os.path.join(OUT, "nuevos_dom.fa"), "w") as fh:
    for g, (a_, b_) in D.items():
        fh.write(f">GL_{t01.loc[g, 'locus']}\n{CORR[g][a_ - 1:b_]}\n")

# --- 3. MAFFT --add --keeplength sobre el alineamiento recortado de 339 secuencias -----------------------------
if len(sys.argv) > 1:
    aln = sh(["mafft", "--quiet", "--add", os.path.join(OUT, "nuevos_dom.fa"), "--keeplength", "--anysymbol", sys.argv[1]])
    open(os.path.join(OUT, "NAC_dom.trim_341.fa"), "w").write(aln)
    largos = {len(str(r.seq)) for r in SeqIO.parse(os.path.join(OUT, "NAC_dom.trim_341.fa"), "fasta")}
    assert len(largos) == 1, largos

# --- 4. lecturas sobre los desplazamientos de marco de los posibles pseudogenes con lecturas ---------------------
BAMS = sorted(os.path.join(RM, f) for f in os.listdir(RM) if f.endswith(".regiones.bam"))
QB = " ".join(shlex.quote(b) for b in BAMS)
flujo = lambda r_, cola: (f"samtools merge -f -R {r_} -o - {QB} | samtools view -u -e '[NH]==1' -F 0x904 - | {cola}")


def bases(s):
    i = n = ind = 0
    while i < len(s):
        ch = s[i]
        if ch == "^": i += 2; continue
        if ch == "$": i += 1; continue
        if ch in "+-":
            j = i + 1; k = ""
            while s[j].isdigit(): k += s[j]; j += 1
            ind += 1; i = j + int(k); continue
        n += ch in ".,ACGTNacgtn*"; i += 1
    return n, ind


PS = {}
d_ps = os.path.join(OUT, "pseudogenes_lecturas"); os.makedirs(d_ps, exist_ok=True)
uj = pd.read_csv(os.path.join(RM, "resultados", "uniones_189.csv"))
for g in ("EglNAC193", "EglNAC196"):
    pm = rev.loc[rev.gen == g, "prot_mejor"].iloc[0]
    r = max([r for r in registros(g) if r["prot"] == pm], key=lambda r: (r["qe"] - r["qs"], -r["fs"] - r["st"]))
    s, off = r["seq"].rsplit(":", 1)[0], int(r["seq"].rsplit(":", 1)[1].split("-")[0]) - 1
    t, sitios = r["ts"] + off - 1, []
    for n_, op in re.findall(r"(\d+)([MIDNUVFG])", r["cg"]):
        n_ = int(n_)
        if op in "MD": t += 3 * n_
        elif op in "FG": sitios.append(dict(pos=t, pb=n_)); t += n_
        elif op in "NUV": t += n_
    assert len(sitios) == r["fs"]
    for x in sitios:
        r_ = f"{s}:{x['pos'] - 2}-{x['pos'] + x['pb'] + 2}"
        out = subprocess.run(["bash", "-o", "pipefail", "-c", flujo(r_, f"samtools mpileup -f {shlex.quote(GEN)} -B -Q 0 -q 0 -d 1000000 -")],
                             capture_output=True, text=True, check=True).stdout
        pb = [bases(f[4]) for f in (l.split("\t") for l in out.splitlines()) if x["pos"] - 2 <= int(f[1]) <= x["pos"] + x["pb"] + 2]
        c_ = cds(g)
        x.update(lecturas_exonicas=min(v[0] for v in pb) if pb else 0, lecturas_con_indel=sum(v[1] for v in pb),
                 en_CDS_anotada=bool(any(a_ <= x["pos"] <= b_ for a_, b_ in zip(c_.ini, c_.fin))))
    PS[g] = dict(locus=t01.loc[g, "locus"], desplazamientos=r["fs"], paradas=r["st"], sitios=sitios,
                 nota="posiciones 1-based en el genoma; lecturas = alineamientos únicos (NH=1) de las 5 bibliotecas de hoja "
                      "que cubren el sitio como exón")
c193 = cds("EglNAC193"); s193 = PS["EglNAC193"]["sitios"][0]
pm = rev.loc[rev.gen == "EglNAC193", "prot_mejor"].iloc[0]
r193 = max([r for r in registros("EglNAC193") if r["prot"] == pm], key=lambda r: (r["qe"] - r["qs"], -r["fs"] - r["st"]))
off = int(r193["seq"].rsplit(":", 1)[1].split("-")[0]) - 1
ex_ps = [(a_ + off, b_ + off) for a_, b_ in r193["cds"] if a_ + off <= s193["pos"] <= b_ + off]
assert len(ex_ps) == 1                      # exón de la proteína homóloga que contiene el desplazamiento
j = uj[(uj.seq == reg.loc["EglNAC193", "seq"]) & (uj.donante_fin == ex_ps[0][1]) & ~uj.anotada & uj.fiable]
assert len(j) == 1, j
j = j.iloc[0]
PS["EglNAC193"]["exon_no_anotado"] = dict(aguas_arriba_de_la_CDS=bool(s193["pos"] < c193.ini.min()),
    union=dict(donante_fin=int(j.donante_fin), aceptor_ini=int(j.aceptor_ini), lecturas=int(j.lecturas),
               bibliotecas=int(j.bibliotecas), anotada=bool(j.anotada), fiable=bool(j.fiable)))
json.dump(PS, open(os.path.join(d_ps, "resumen_pseudogenes_lecturas.json"), "w", encoding="utf-8"), ensure_ascii=False, indent=1)

# --- 5. motivos A-E sobre las proteínas corregidas -------------------------------------------------------------
anot = {}
for g, m in mod[mod.tipo == "CDS"].groupby("gen"):
    m = m.sort_values("ini"); anot[g] = traducir(list(zip(m.ini, m.fin)), m.seq.iloc[0], m.hebra.iloc[0])
t04 = pd.read_csv(os.path.join(TB, "T04_motivos_MEME.csv")); t05 = pd.read_csv(os.path.join(TB, "T05_motivos_sitios.csv"))
AA = "ACDEFGHIKLMNPQRSTVWY"; bg = collections.Counter("".join(anot.values()))
bgv = np.array([bg[a_] for a_ in AA], float); bgv /= bgv.sum()
SUB = t04.set_index("motivo").subdominio; DOM7 = t04[t04.subdominio.isin(list("ABCDE"))].motivo.tolist()
llamados = {(r.simbolo_EglNAC, r.motivo) for r in t05.itertuples()}
puntua = lambda L, x: sum(L[i, AA.index(c)] for i, c in enumerate(x) if c in AA)
mejor = lambda L, p: max(puntua(L, p[i:i + L.shape[0]]) for i in range(len(p) - L.shape[0] + 1))
MC = {"metodo": "matriz de pesos por posición construida con los sitios MEME de T05 (seudoconteo = frecuencia de fondo); "
                "umbral = percentil 5 de las puntuaciones de esos sitios", "motivos": {}, "genes": {}}
for mtf in DOM7:
    w = int(t04.set_index("motivo").loc[mtf, "ancho"])
    sit = [anot[r.simbolo_EglNAC][r.inicio - 1:r.inicio - 1 + w] for r in t05[t05.motivo == mtf].itertuples()]
    sit = [x for x in sit if len(x) == w]
    M = np.zeros((w, 20))
    for x in sit:
        for i, c in enumerate(x):
            if c in AA: M[i, AA.index(c)] += 1
    L = np.log2(((M + bgv) / (M.sum(1, keepdims=True) + 1)) / bgv)
    umbral = float(np.percentile([puntua(L, x) for x in sit], 5))
    b_ = {g: mejor(L, p) for g, p in anot.items()}
    MC["motivos"][mtf] = dict(subdominio=SUB[mtf], umbral=round(umbral, 1),
                              coinciden_MEME=sum(b_[g] >= umbral and (g, mtf) in llamados for g in anot),
                              solo_matriz=sum(b_[g] >= umbral and (g, mtf) not in llamados for g in anot),
                              solo_MEME=sum(b_[g] < umbral and (g, mtf) in llamados for g in anot))
    for g, p in CORR.items():
        d_ = MC["genes"].setdefault(g, {"MEME_anotada": [], "matriz_corregida": []})
        if (g, mtf) in llamados: d_["MEME_anotada"].append(mtf)
        if mejor(L, p) >= umbral: d_["matriz_corregida"].append(mtf)
for g, d_ in MC["genes"].items():
    d_["union"] = sorted(set(d_["MEME_anotada"]) | set(d_["matriz_corregida"]))
    d_["n_MEME_anotada"], d_["n_union"] = len(d_["MEME_anotada"]), len(d_["union"])
    d_["nuevos"] = sorted(set(d_["matriz_corregida"]) - set(d_["MEME_anotada"]))
    d_["subdominios_nuevos"] = sorted({SUB[m_] for m_ in d_["nuevos"]})
json.dump(MC, open(os.path.join(OUT, "motivos_modelos_corregidos.json"), "w", encoding="utf-8"), ensure_ascii=False, indent=1,
          default=lambda o: o.item() if hasattr(o, "item") else str(o))
print({g: len(p) for g, p in CORR.items()}, D, {g: [(x["pos"], x["lecturas_exonicas"]) for x in v["sitios"]] for g, v in PS.items()})
