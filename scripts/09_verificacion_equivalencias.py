#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
09_verificacion_equivalencias.py — Revisión NAC (E. globulus)

Verifica POR SECUENCIA, no por nombre, la cadena de equivalencias que sostiene
el capítulo de expresión y el paso estructural:

    partidor qPCR -> transcrito (genoma E. globulus 2024, GCA_046226855.1)
                  -> locus ACJRO7_* -> símbolo del censo (T02)
                  -> co-ortólogo EgrNAC (Hussey 2015 Tabla S4)
                  -> locus Eucgr.* de E. grandis v1 -> fila del DE

Motivo: el manuscrito está anclado a E. grandis v1 (los DE son de globulus pero
referenciados a Eucgr.*), el censo al genoma conespecífico de 2024, y los
partidores se diseñaron sobre la asamblea propia del manuscrito. Tres sistemas
de coordenadas simultáneos y cuatro nomenclaturas para los mismos genes.

Salidas (en ANALISIS_Eglobulus/tablas/):
  T24_verificacion_equivalencias.csv  una fila por gen medido por qPCR
  T25_trazabilidad_T10.csv            una fila por valor de RNA-seq de T10

Uso:  python3 scripts/09_verificacion_equivalencias.py
"""
import csv, gzip, os, re, sys, collections

R  = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TB = os.path.join(R, "ANALISIS_Eglobulus", "tablas")
PA = os.path.join(R, "1ra entrega data paula")

# ---------------------------------------------------------------- utilidades
def rc(s):
    return s[::-1].translate(str.maketrans("ACGTacgt", "TGCAtgca"))

def leer_fasta(ruta, op=open):
    d, k = {}, None
    with op(ruta, "rt", encoding="utf-8", errors="replace") as fh:
        for ln in fh:
            if ln.startswith(">"):
                k = ln[1:].strip(); d[k] = []
            else:
                d[k].append(ln.strip())
    return {k: "".join(v).upper() for k, v in d.items()}

def desap(a, b):
    return sum(1 for x, y in zip(a, b) if x != y)

def mejor_ventana(p, s):
    """Mínimo de desapareamientos del partidor p en la secuencia s."""
    if len(s) < len(p):
        return len(p), -1
    return min((desap(p, s[i:i+len(p)]), i) for i in range(len(s) - len(p) + 1))

def identidad(a, b, abre=-8, coin=5, dif=-4):
    """Needleman-Wunsch con penalización lineal; devuelve (idénticas, columnas, %)."""
    n, m = len(a), len(b)
    prev = [abre * j for j in range(m + 1)]
    ptr = []
    for i in range(1, n + 1):
        cur = [abre * i] + [0] * m
        fila = []
        for j in range(1, m + 1):
            d = prev[j-1] + (coin if a[i-1] == b[j-1] else dif)
            cur[j] = max(d, prev[j] + abre, cur[j-1] + abre)
            fila.append(0 if cur[j] == d else (1 if cur[j] == prev[j] + abre else 2))
        ptr.append(fila); prev = cur
    i, j, ig, al = n, m, 0, 0
    while i > 0 and j > 0:
        t = ptr[i-1][j-1]
        if t == 0:
            ig += (a[i-1] == b[j-1]); al += 1; i -= 1; j -= 1
        elif t == 1:
            i -= 1
        else:
            j -= 1
    return ig, al, (100.0 * ig / al if al else 0.0)

# ---------------------------------------------------------------- entradas
t02 = list(csv.DictReader(open(os.path.join(TB, "T02_equivalencia_nomenclaturas_189.csv"), encoding="utf-8")))
t01 = {r["locus"]: r for r in csv.DictReader(open(os.path.join(TB, "T01_censo_maestro_189.csv"), encoding="utf-8"))}
t10 = list(csv.DictReader(open(os.path.join(TB, "T10_correlacion_qpcr_rnaseq.csv"), encoding="utf-8")))
s4  = {r["simbolo_EgrNAC"]: r for r in csv.DictReader(open(os.path.join(R, "Hussey_TablaS4_189simbolos.csv"), encoding="utf-8"))}
aud = {r["locus"]: r for r in csv.DictReader(open(os.path.join(R, "_respaldo_app", "auditoria_loci_simbolos.csv"), encoding="utf-8"))}
prim = [l.split("\t") for l in open(os.path.join(R, "rnaseq_udec", "primers.tsv"), encoding="utf-8").read().splitlines() if l.strip()]

trans = leer_fasta(os.path.join(R, "rnaseq_udec", "ref", "transcritos.fa"))
prot  = leer_fasta(os.path.join(R, "rnaseq_udec", "ref", "GCA_046226855.1_ASM4622685v1_protein.faa.gz"), gzip.open)
manu  = leer_fasta(os.path.join(R, "_respaldo_app", "EglNAC_proteinas.fasta"))

# DE del proyecto (reanálisis propio, 24 muestras, loci Eucgr de E. grandis v1)
de = collections.defaultdict(dict)
for r in csv.DictReader(open(os.path.join(R, "DE_NAC_modelo_conjunto.csv"), encoding="utf-8")):
    de[r["locus"]][r["contraste"]] = (float(r["logFC"]), float(r["FDR"]))

# edgeR original de la colaboradora: de aquí salieron las cifras publicadas
PAULA = {"C_DS": "verdadero.txt",
         "NA_CAAF": "Tdos_alltest4output_NA_CAAF_verdadero.txt",
         "NA_CABF": "todos_alltest4output_NA_CABF_tagwise.txt",
         "CABF_CAAF": "todos_alltest4output_CABF_CAAF_verdadero.txt"}
orig = {}
for cc, f in PAULA.items():
    ruta = os.path.join(PA, f)
    if not os.path.exists(ruta):
        continue
    for ln in open(ruta, encoding="utf-8"):
        p = ln.rstrip("\n").split("\t")
        if p[0].startswith('"Eucgr.'):
            orig[(p[0].strip('"'), cc)] = (float(p[1]), float(p[3]))

# ---------------------------------------------------------------- 1) PCR in silico
def pcr(f, r, tol=0):
    """Amplicones (transcrito, largo, desap_F, desap_R) en todo el transcriptoma."""
    F, Rc, out = f.upper(), rc(r.upper()), []
    for t, s in trans.items():
        i = s.find(F) if tol == 0 else -1
        if tol == 0:
            if i < 0:
                continue
            j = s.find(Rc, i)
            if j < 0:
                continue
            out.append((t.split("|")[-1], j + len(Rc) - i, 0, 0))
        else:
            eF, iF = mejor_ventana(F, s)
            eR, iR = mejor_ventana(Rc, s)
            if eF <= tol and eR <= tol and 50 < iR + len(Rc) - iF < 400:
                out.append((t.split("|")[-1], iR + len(Rc) - iF, eF, eR))
    return out

# ---------------------------------------------------------------- 2) verificación por gen
qpcr = [r for r in t02 if r["medido_por_qPCR"] == "True"]
filas = []
for r in qpcr:
    sim, loc = r["simbolo"], r["locus"]
    tr_decl = r["transcrito"].split("|")[-1]
    par = next((p for p in prim if p[0] in (r["etiqueta_FileS1"], r["simbolo_manuscrito_2026"])), None)
    # amplicón declarado
    amp_ex = pcr(par[1], par[2]) if par else []
    tr_obs = amp_ex[0][0] if len(amp_ex) == 1 else ""
    largo_obs = amp_ex[0][1] if len(amp_ex) == 1 else ""
    # si no pega exacto, medir desapareamientos contra el transcrito declarado
    clave = next((k for k in trans if tr_decl in k), None)   # la cabecera lleva descripción
    dF = dR = ""
    if par and clave:
        dF = mejor_ventana(par[1].upper(), trans[clave])[0]
        dR = mejor_ventana(rc(par[2].upper()), trans[clave])[0]
        if not tr_obs and dF <= 2 and dR <= 2:
            iF = mejor_ventana(par[1].upper(), trans[clave])[1]
            iR = mejor_ventana(rc(par[2].upper()), trans[clave])[1]
            tr_obs, largo_obs = tr_decl + " (con desapareamiento)", iR + len(par[2]) - iF
    # proteína: locus del censo vs proteína del manuscrito
    pk = next((k for k in prot if loc.split("_")[-1] in k), None)
    mk = next((k for k in manu if re.search(r"%s\b" % re.escape(r["simbolo_manuscrito_2026"]), k)), None)
    pid = identidad(prot[pk], manu[mk])[2] if pk and mk else float("nan")
    # eslabón a E. grandis v1
    eg = r["coortologo_EgrNAC"]
    euc = s4.get(eg, {}).get("locus_v1_0", "")
    veredicto = []
    if tr_obs.startswith(tr_decl):
        veredicto.append("partidor confirma el transcrito declarado")
    elif tr_obs:
        veredicto.append("el partidor amplifica OTRO transcrito: " + tr_obs)
    else:
        veredicto.append("partidor sin amplicón identificable")
    if pid == pid and pid >= 95:
        veredicto.append("proteína del censo = proteína del manuscrito (%.1f%%)" % pid)
    elif pid == pid:
        veredicto.append("proteína discordante (%.1f%%)" % pid)
    if r["etiqueta_FileS1"] and r["etiqueta_FileS1"] != r["simbolo_manuscrito_2026"]:
        veredicto.append("ATENCIÓN: File S1 lo rotula %s pero el locus es %s"
                         % (r["etiqueta_FileS1"], r["simbolo_manuscrito_2026"]))
    filas.append({
        "simbolo_censo": sim, "locus_2024": loc, "transcrito_declarado": tr_decl,
        "partidor": par[0] if par else "", "amplicon_transcrito": tr_obs,
        "amplicon_pb": largo_obs, "desap_F_en_declarado": dF, "desap_R_en_declarado": dR,
        "n_amplicones_exactos": len(amp_ex),
        "simbolo_manuscrito": r["simbolo_manuscrito_2026"], "etiqueta_FileS1": r["etiqueta_FileS1"],
        "detector_placa": r["detector_placa"], "coortologo_EgrNAC": eg, "locus_Egrandis_v1": euc,
        "largo_aa_2024": len(prot[pk]) if pk else "", "largo_aa_manuscrito": len(manu[mk]) if mk else "",
        "largo_aa_S4_EgrNAC": s4.get(eg, {}).get("largo_aa", ""),
        "identidad_censo_vs_manuscrito_pct": round(pid, 1) if pid == pid else "",
        "veredicto": "; ".join(veredicto)})

with open(os.path.join(TB, "T24_verificacion_equivalencias.csv"), "w", newline="", encoding="utf-8") as fh:
    w = csv.DictWriter(fh, fieldnames=list(filas[0])); w.writeheader(); w.writerows(filas)

# ---------------------------------------------------------------- 3) trazabilidad de T10
CM = {"C_DS": "Leaf_NA_DS", "NA_CABF": "Leaf_NA_CABF", "NA_CAAF": "Leaf_NA_CAAF"}
mapa = {r["simbolo"]: r for r in t02}
tr10 = []
for r in t10:
    sim, cc = r["simbolo_unico_EglNAC"], r["contraste"]
    val = float(r["rnaseq_log2FC"])
    eg = mapa.get(sim, {}).get("coortologo_EgrNAC", "")
    euc = s4.get(eg, {}).get("locus_v1_0", "")
    vo = orig.get((euc, cc), (None, None))          # edgeR original (publicado)
    vp = de.get(euc, {}).get(CM.get(cc, ""), (None, None))  # reanálisis del proyecto
    # ¿de qué gen es realmente el valor publicado?
    dueno = [L for (L, c), (fc, _) in orig.items() if c == cc and abs(fc - val) < 1e-4]
    dueno_sim = sorted({aud.get(L, {}).get("EglNAC", L) for L in dueno})
    if vo[0] is not None and abs(vo[0] - val) < 1e-3:
        ver = "correcto: coincide con el edgeR original del locus"
    elif vo[0] is None:
        ver = "el locus no está en ese contraste (filtrado); el valor de T10 no puede venir de aquí"
    else:
        ver = "INCORRECTO: el locus vale %.4f en ese contraste" % vo[0]
        if dueno_sim:
            ver += "; el valor publicado pertenece a %s" % ", ".join(dueno_sim[:3])
    tr10.append({"simbolo_censo": sim, "contraste": cc, "valor_en_T10": val,
                 "coortologo_EgrNAC": eg, "locus_Egrandis_v1": euc,
                 "edgeR_original_logFC": "" if vo[0] is None else round(vo[0], 4),
                 "edgeR_original_FDR": "" if vo[1] is None else "%.3g" % vo[1],
                 "reanalisis_proyecto_logFC": "" if vp[0] is None else round(vp[0], 4),
                 "reanalisis_proyecto_FDR": "" if vp[1] is None else "%.3g" % vp[1],
                 "veredicto": ver})

with open(os.path.join(TB, "T25_trazabilidad_T10.csv"), "w", newline="", encoding="utf-8") as fh:
    w = csv.DictWriter(fh, fieldnames=list(tr10[0])); w.writeheader(); w.writerows(tr10)

# ---------------------------------------------------------------- 4) auditoría de la regla de T02
num = lambda s: int(re.sub(r"\D", "", s) or 0)
reglas = collections.Counter(r["regla"].split(":")[0] for r in t02)
hereda = [r for r in t02 if r["regla"].startswith("hereda")]
mal_h = [r["simbolo"] for r in hereda if num(r["simbolo"]) != num(r["coortologo_EgrNAC"])]
comp = collections.Counter(r["coortologo_EgrNAC"] for r in t02 if r["coortologo_EgrNAC"].strip())
exp = [r for r in t02 if r["regla"].startswith("parálogo")]
mal_e = [r["simbolo"] for r in exp if comp[r["coortologo_EgrNAC"]] < 2]

print("T02: %d filas | reglas: %s" % (len(t02), dict(reglas)))
print("  heredan número y lo cumplen: %d/%d %s" % (len(hereda) - len(mal_h), len(hereda),
                                                   "" if not mal_h else "fallan: " + ", ".join(mal_h)))
print("  parálogos de expansión con co-ortólogo compartido: %d/%d %s" % (len(exp) - len(mal_e), len(exp),
                                                   "" if not mal_e else "fallan: " + ", ".join(mal_e)))
print("T24: %d genes verificados por secuencia -> %s" % (len(filas), "T24_verificacion_equivalencias.csv"))
for f in filas:
    print("   %-9s %-14s %s" % (f["simbolo_censo"], f["locus_2024"], f["veredicto"]))
ok = sum(1 for t in tr10 if t["veredicto"].startswith("correcto"))
print("T25: %d/%d valores de RNA-seq de T10 trazables al edgeR original -> T25_trazabilidad_T10.csv" % (ok, len(tr10)))
for t in tr10:
    if not t["veredicto"].startswith("correcto"):
        print("   %-9s %-9s T10=%-9s %s" % (t["simbolo_censo"], t["contraste"], t["valor_en_T10"], t["veredicto"]))
