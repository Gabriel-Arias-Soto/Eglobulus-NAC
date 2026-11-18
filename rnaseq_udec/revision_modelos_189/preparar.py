# -*- coding: utf-8 -*-
"""Revisión de los modelos génicos de los 189 NAC (paso 1: entradas).

Escribe, en esta carpeta:
  - modelos_anotados_189.tsv: exones y CDS del transcrito representativo de cada NAC (T01.transcrito), de anot.gtf;
  - regiones_5kb.bed: cada locus ±5 kb (filtra los alineamientos de RNA-seq);
  - proteinas_evidencia.faa: proteínas NAC de E. grandis v2 (co-ortólogos y ortólogos colineales NAC de T48) y de
    A. thaliana (NAC de UniProt, Araport11), más todas las que tienen el dominio NAM en hmmsearch --cut_ga, para miniprot. Las de E. globulus no se usan (serían circulares).
Para si falta un transcrito, si un locus no tiene CDS o si alguna proteína pedida no está en los FASTA.
"""
import csv, glob, os, re, collections

AQUI = os.path.dirname(os.path.abspath(__file__))
R = os.path.dirname(os.path.dirname(AQUI))
A = os.path.join(R, "ANALISIS_Eglobulus")
T = lambda c: list(csv.DictReader(open(glob.glob(os.path.join(A, "tablas", f"{c}_*.csv"))[0], encoding="utf-8")))
t01 = T("T01"); t48 = T("T48")
assert len(t01) == 189
TX = {x["transcrito"]: x["simbolo"] for x in t01}
assert len(TX) == 189

# --- modelos anotados -----------------------------------------------------------------------------------
filas = []
rx = re.compile(r'transcript_id "([^"]+)"')
with open(os.path.join(R, "rnaseq_udec", "ref", "anot.gtf"), encoding="utf-8") as fh:
    for l in fh:
        if l.startswith("#"):
            continue
        c = l.rstrip("\n").split("\t")
        if c[2] not in ("exon", "CDS"):
            continue
        m = rx.search(c[8])
        if m and m.group(1) in TX:
            filas.append((TX[m.group(1)], c[0], c[2], int(c[3]), int(c[4]), c[6]))
por = collections.defaultdict(collections.Counter)
for g, s, tp, a, b, h in filas:
    por[g][tp] += 1
assert len(por) == 189 and all(v["CDS"] > 0 for v in por.values()), [g for g in TX.values() if por[g]["CDS"] == 0]
with open(os.path.join(AQUI, "modelos_anotados_189.tsv"), "w") as fo:
    fo.write("gen\tseq\ttipo\tini\tfin\thebra\n")
    for f in sorted(filas, key=lambda x: (x[0], x[2], x[3])):
        fo.write("\t".join(map(str, f)) + "\n")

# --- regiones ±5 kb ---------------------------------------------------------------------------------------
fai = {l.split("\t")[0]: int(l.split("\t")[1]) for l in open(os.path.join(R, "rnaseq_udec", "ref", "genoma.fa.fai"))}
reg = {}
for g, s, tp, a, b, h in filas:
    x = reg.setdefault(g, [s, a, b])
    x[1], x[2] = min(x[1], a), max(x[2], b)
with open(os.path.join(AQUI, "regiones_5kb.bed"), "w") as fo:
    for g, (s, a, b) in sorted(reg.items(), key=lambda kv: (kv[1][0], kv[1][1])):
        fo.write(f"{s}\t{max(0, a - 1 - 5000)}\t{min(fai[s], b + 5000)}\t{g}\n")


# --- proteínas de evidencia ------------------------------------------------------------------------------
def fasta(f):
    d, k = {}, None
    for l in open(f, encoding="utf-8"):
        if l.startswith(">"):
            k = l[1:].split()[0]; d[k] = []
        else:
            d[k].append(l.strip())
    return {k: "".join(v) for k, v in d.items()}


egr = fasta(os.path.join(A, "datos", "entrega3", "isoformas", "Egr.pep"))
ath = fasta(os.path.join(A, "datos", "entrega3", "isoformas", "Ath.pep"))
ids_egr = {x["coortologo_v2"] for x in t48 if x["coortologo_v2"].startswith("EUGRSUZ_")}
ids_egr |= {x["ortologo_Egr"] for x in t48 if x["ortologo_Egr"].startswith("EUGRSUZ_") and x["ortologo_Egr_es_NAC"] == "True"}
up = list(csv.DictReader(open(os.path.join(A, "datos", "entrega3", "go", "uniprot_NAC_Ath_UP000006548.tsv"), encoding="utf-8"), delimiter="\t"))
ids_ath = {a.strip() for x in up for a in x["Araport"].split(";") if a.strip()}
# además, todas las proteínas con el dominio NAM sobre el umbral de Pfam (hmmsearch --cut_ga; Egr.dom y Ath.dom, en esta carpeta)
dom = lambda f: {l.split()[0] for l in open(os.path.join(AQUI, f)) if not l.startswith("#")}
n_t48, n_up = len(ids_egr), len(ids_ath)
ids_egr |= dom("Egr.dom"); ids_ath |= dom("Ath.dom")
falta = sorted((ids_egr - set(egr)) | (ids_ath - set(ath)))
assert not falta, falta[:10]
with open(os.path.join(AQUI, "proteinas_evidencia.faa"), "w") as fo:
    for k in sorted(ids_egr):
        fo.write(f">Egr|{k}\n{egr[k]}\n")
    for k in sorted(ids_ath):
        fo.write(f">Ath|{k}\n{ath[k]}\n")
print("modelos:", len(por), "| filas:", len(filas), "| regiones:", len(reg), "| proteínas Egr:", len(ids_egr), f"(T48: {n_t48})", "Ath:", len(ids_ath), f"(UniProt: {n_up})")
