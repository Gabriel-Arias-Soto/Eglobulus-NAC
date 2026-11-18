#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""14 - ¿de quien es ortologo el gen focal: de ANAC019 o de ANAC072/RD26?

    Uso:  cd <raiz del proyecto> && python scripts/14_ortologia_SNAC.py
    Entorno conda `filo` (mafft, trimal, iqtree, biopython).

El manuscrito afirma la ortologia con ANAC072/RD26 apoyandose en un porcentaje de
identidad. MAPA_revisores.md pide rehacerlo con arbol de genes y soporte de clado,
y tiene razon por un motivo que conviene dejar escrito: la identidad no distingue
ortologia de paralogia. ANAC019, ANAC055 y ANAC072 son parientes muy proximos
entre si, de modo que los tres dan identidades parecidas contra cualquier NAC de
estres de otra especie, y elegir al mayor es elegir ruido.

Construye un arbol enfocado: los 10 NAC de subfamilia Vb de E. globulus (T02),
el trio de estres de Arabidopsis, sus dos hermanos ATAF, SNAC1 de arroz y dos NAC
de otra subfamilia como grupo externo. MAFFT L-INS-i -> trimAl -> IQ-TREE con
ModelFinder, 1000 ultrafast bootstrap y 1000 replicas SH-aLRT.

Produce ademas dos cosas que el arbol por si solo no contesta:

  T30  identidades dominio/completa contra cada ancla, declarando isoforma,
       matriz, penalizaciones y region. Cierra el punto C2 de
       DATOS_a_solicitar_Paula.md (el manuscrito dice 62,3 % y no reproduce).
  T31  que residuo lleva cada ortologo en las ocho posiciones que leen bases
       nitrogenadas en el co-cristal 3SWP (T18). Si las ocho son identicas en los
       tres parologos de Arabidopsis, la cara de union al ADN no puede decidir
       entre ellos y hay que decirlo en vez de insinuar lo contrario.
"""
import os, re, csv, glob, json, subprocess, urllib.request, urllib.parse

# --- rutas y contrato de estilo -------------------------------------------
RAIZ = os.getcwd()
assert os.path.isfile(os.path.join(RAIZ, "scripts", "00_estilo.R")), \
    "ejecutar desde la raiz del proyecto: python scripts/14_ortologia_SNAC.py"
ANALISIS = os.path.join(RAIZ, "ANALISIS_Eglobulus")
TABLAS   = os.path.join(ANALISIS, "tablas")
DATOS    = os.path.join(ANALISIS, "datos", "snac")
CIFRAS   = os.path.join(TABLAS, "_cifras")
for d in (DATOS, CIFRAS):
    os.makedirs(d, exist_ok=True)

_estilo = open(os.path.join(RAIZ, "scripts", "00_estilo.R"), encoding="utf-8").read()
def cte(nombre):
    m = re.search(r'^%s\s*<-\s*"?([^"#\n]+?)"?\s*(?:#.*)?$' % nombre, _estilo, re.M)
    assert m, "no encuentro %s en 00_estilo.R" % nombre
    return m.group(1).strip()
GEN    = cte("GEN_FOCAL")
LOCUS  = cte("GEN_FOCAL_LOCUS")
SUBFAM = cte("GEN_FOCAL_SUBFAM")

def leer(tabla):
    with open(os.path.join(TABLAS, tabla), encoding="utf-8") as fh:
        return list(csv.DictReader(fh))

# --- 1. conjunto de secuencias --------------------------------------------
# Anclas de Arabidopsis. El trio de estres, sus dos hermanos ATAF (el grupo
# hermano inmediato, sin el cual el trio parece mas aislado de lo que esta) y dos
# NAC de otra subfamilia que hacen de grupo externo para enraizar.
ANCLAS = [("AT1G52890", "ANAC019",  "trio de estres"),
          ("AT3G15500", "ANAC055",  "trio de estres"),
          ("AT4G27410", "RD26",     "trio de estres"),
          ("AT1G01720", "ATAF1",    "hermano ATAF"),
          ("AT5G08790", "ATAF2",    "hermano ATAF"),
          ("AT5G39610", "ORE1",     "grupo externo"),
          ("AT1G69490", "NAP",      "grupo externo")]
ARROZ = ("Q75HE5", "SNAC1", "SNAC-A de arroz")

def bajar_uniprot(destino):
    """Isoforma canonica de SwissProt. Se guarda en disco: la cifra de identidad
    depende de la isoforma y el proyecto ya tiene una discrepancia por no saber
    cual uso el manuscrito."""
    filas = []
    with open(destino, "w", encoding="utf-8") as fh:
        for locus, nombre, _ in ANCLAS:
            q = "(%s)AND(organism_id:3702)AND(reviewed:true)" % locus
            url = "https://rest.uniprot.org/uniprotkb/search?format=fasta&query=" + \
                  urllib.parse.quote(q, safe="()&:")
            ent = [e for e in urllib.request.urlopen(url).read().decode().split(">") if e.strip()]
            assert len(ent) == 1, "%s devuelve %d entradas revisadas" % (locus, len(ent))
            cab, seq = ent[0].split("\n", 1); seq = "".join(seq.split())
            fh.write(">%s\n%s\n" % (nombre, seq))
            filas.append((nombre, cab.split("|")[1], locus, len(seq)))
        txt = urllib.request.urlopen(
            "https://rest.uniprot.org/uniprotkb/%s.fasta" % ARROZ[0]).read().decode()
        seq = "".join(txt.split("\n")[1:])
        fh.write(">%s\n%s\n" % (ARROZ[1], seq))
        filas.append((ARROZ[1], ARROZ[0], "Os03g0815100", len(seq)))
    return filas

def bajar_ncbi(destino, acc):
    url = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein"
           "&rettype=fasta&retmode=text&id=" + ",".join(acc.values()))
    fa = urllib.request.urlopen(url).read().decode()
    inv = {v: k for k, v in acc.items()}
    filas = []
    with open(destino, "w", encoding="utf-8") as fh:
        for e in (x for x in fa.split(">") if x.strip()):
            cab, seq = e.split("\n", 1); seq = "".join(seq.split())
            sim = inv[cab.split()[0]]
            fh.write(">%s\n%s\n" % (sim, seq))
            filas.append((sim, cab.split()[0], len(seq)))
    assert len(filas) == len(acc), "NCBI devolvio %d de %d" % (len(filas), len(acc))
    return filas

t02 = leer("T02_equivalencia_nomenclaturas_189.csv")
vb  = [r for r in t02 if r["subfamilia_final"] == SUBFAM]
assert any(r["simbolo"] == GEN for r in vb), "%s no aparece en la subfamilia %s" % (GEN, SUBFAM)
l2p = dict(l.split("\t") for l in
           open(os.path.join(ANALISIS, "datos", "locus2prot.tsv"), encoding="utf-8").read().splitlines()
           if "\t" in l)
acc_vb = {r["simbolo"]: l2p[r["locus"]] for r in vb}

f_anc, f_euc = os.path.join(DATOS, "anclas.faa"), os.path.join(DATOS, "vb.faa")
meta_anc = bajar_uniprot(f_anc) if not os.path.isfile(f_anc) else None
meta_euc = bajar_ncbi(f_euc, acc_vb) if not os.path.isfile(f_euc) else None

def leer_faa(path):
    d, k = {}, None
    for l in open(path, encoding="utf-8"):
        if l.startswith(">"): k = l[1:].strip().split()[0]; d[k] = ""
        else: d[k] += l.strip()
    return d
seqs = leer_faa(f_anc); seqs.update(leer_faa(f_euc))
assert GEN in seqs and "RD26" in seqs and "ANAC019" in seqs

# El gen focal en disco debe ser el mismo que el modelado en la seccion estructural.
disco = "".join(l.strip() for l in open(os.path.join(ANALISIS, "datos", "%s.faa" % GEN))
                if not l.startswith(">"))
assert seqs[GEN] == disco, "la secuencia descargada de %s no coincide con datos/%s.faa" % (GEN, GEN)

# Metadatos de cada punta, para que la figura no teclee grupos ni especies.
ETQ = os.path.join(DATOS, "etiquetas.csv")
with open(ETQ, "w", newline="", encoding="utf-8") as fh:
    w = csv.writer(fh); w.writerow(["id", "grupo", "especie", "referencia", "focal"])
    for locus, nombre, grupo in ANCLAS:
        w.writerow([nombre, grupo, "A. thaliana", locus, False])
    w.writerow([ARROZ[1], ARROZ[2], "O. sativa", ARROZ[0], False])
    for r in vb:
        w.writerow([r["simbolo"], "subfamilia %s de E. globulus" % SUBFAM, "E. globulus",
                    acc_vb[r["simbolo"]], r["simbolo"] == GEN])

CONJ = os.path.join(DATOS, "conjunto.faa")
with open(CONJ, "w", encoding="utf-8") as fh:
    for k, v in seqs.items():
        fh.write(">%s\n%s\n" % (k, v))
print("  %d secuencias -> %s" % (len(seqs), os.path.relpath(CONJ, RAIZ)))

# --- 2. alineamiento -------------------------------------------------------
MSA = os.path.join(DATOS, "conjunto_mafft.faa")
with open(MSA, "w") as fh:
    subprocess.run(["mafft", "--localpair", "--maxiterate", "1000", "--quiet", CONJ],
                   stdout=fh, check=True)
ali = leer_faa(MSA)
ancho = len(next(iter(ali.values())))
assert all(len(v) == ancho for v in ali.values())
print("  alineamiento L-INS-i: %d columnas" % ancho)

# Columnas del alineamiento que ocupa el dominio NAM del gen focal, segun los
# limites del HMM en T16 (no los del modelo 3D, que son algo mas anchos).
t16 = {r["simbolo"]: r for r in leer("T16_arquitectura_dominios.csv")}
dom_ini, dom_fin = int(t16[GEN]["ali_ini"]), int(t16[GEN]["ali_fin"])
col_de_pos, pos = {}, 0
for c, ch in enumerate(ali[GEN]):
    if ch != "-":
        pos += 1; col_de_pos[pos] = c
cols_dom = set(range(col_de_pos[dom_ini], col_de_pos[dom_fin] + 1))
print("  dominio NAM de %s: %d-%d aa (T16) = %d columnas" % (GEN, dom_ini, dom_fin, len(cols_dom)))

# --- 3. identidades, con el metodo declarado -------------------------------
from Bio import Align
from Bio.Align import substitution_matrices
aln = Align.PairwiseAligner(mode="global", substitution_matrix=substitution_matrices.load("BLOSUM62"),
                            open_gap_score=-10.0, extend_gap_score=-0.5)
SEMILLA = "20260919"   # IQ-TREE remuestrea: sin semilla fija el soporte y el AU se mueven
METODO = "Needleman-Wunsch global, BLOSUM62, apertura 10, extension 0,5; identidad = identicos / columnas sin hueco"

def identidad(a, b, region=None):
    """region: conjunto de posiciones (1-based) de `a` a las que restringir el
    recuento. Se usa la misma alineacion global para dominio y proteina completa,
    asi las dos cifras son comparables entre si."""
    par = aln.align(a, b)[0]
    A, B = par[0], par[1]
    ia = 0; ident = tot = 0
    for x, y in zip(A, B):
        if x != "-": ia += 1
        if x == "-" or y == "-": continue
        if region is not None and ia not in region: continue
        tot += 1; ident += (x == y)
    return ident, tot

dom_focal = set(range(dom_ini, dom_fin + 1))
filas30 = []
for nombre in [a[1] for a in ANCLAS] + [ARROZ[1]]:
    i_f, n_f = identidad(seqs[GEN], seqs[nombre])
    i_d, n_d = identidad(seqs[GEN], seqs[nombre], dom_focal)
    grupo = dict((a[1], a[2]) for a in ANCLAS); grupo[ARROZ[1]] = ARROZ[2]
    filas30.append(dict(ancla=nombre, grupo=grupo[nombre], largo_aa=len(seqs[nombre]),
                        ident_dominio_pct=round(100*i_d/n_d, 1), n_dominio=n_d,
                        ident_completa_pct=round(100*i_f/n_f, 1), n_completa=n_f))
filas30.sort(key=lambda r: -r["ident_dominio_pct"])
with open(os.path.join(TABLAS, "T30_identidad_%s_vs_anclas.csv" % GEN), "w", newline="", encoding="utf-8") as fh:
    w = csv.DictWriter(fh, fieldnames=list(filas30[0]) + ["metodo"])
    w.writeheader()
    for r in filas30: w.writerow(dict(r, metodo=METODO))
print("  T30: %d anclas" % len(filas30))

# --- 4. los ocho lectores de base, en cada ortologo ------------------------
# T18 numera sobre ANAC019 (la cadena del co-cristal). Se pasa esa posicion a
# columna del alineamiento y se lee que tiene cada secuencia en esa columna.
t18 = leer("T18_lectura_bases_3SWP.csv")
col_anac019, p = {}, 0
for c, ch in enumerate(ali["ANAC019"]):
    if ch != "-":
        p += 1; col_anac019[p] = c
orden = ["ANAC019", "ANAC055", "RD26", "ATAF1", "ATAF2", "SNAC1", GEN] + \
        [r["simbolo"] for r in vb if r["simbolo"] != GEN]
filas31 = []
for r in t18:
    pos19 = int(r["pos_ANAC019"]); c = col_anac019[pos19]
    fila = dict(pos_ANAC019=pos19, aa_ANAC019=r["aa_ANAC019"],
                n_contactos_base=int(r["n_contactos_base"]))
    for s in orden: fila[s] = ali[s][c]
    fila["identico_en_el_trio"] = len({ali[s][c] for s in ("ANAC019", "ANAC055", "RD26")}) == 1
    assert fila["aa_ANAC019"] == ali["ANAC019"][c], \
        "la columna %d no lleva el aa que dice T18 para la posicion %d" % (c, pos19)
    filas31.append(fila)
with open(os.path.join(TABLAS, "T31_lectores_base_en_ortologos.csv"), "w", newline="", encoding="utf-8") as fh:
    w = csv.DictWriter(fh, fieldnames=list(filas31[0])); w.writeheader(); w.writerows(filas31)
n_trio = sum(f["identico_en_el_trio"] for f in filas31)
n_focal = sum(f[GEN] == f["aa_ANAC019"] for f in filas31)
print("  T31: %d/%d posiciones identicas en el trio; %d/%d conservadas en %s"
      % (n_trio, len(filas31), n_focal, len(filas31), GEN))

# --- 5. arboles ------------------------------------------------------------
# Dos alineamientos, a proposito. El dominio NAM es la practica establecida en la
# familia NAC porque la region C-terminal no es alineable entre subfamilias; la
# version completa se conserva como analisis de sensibilidad, para poder decir si
# la conclusion depende de incluir una cola que el propio gen focal tiene larga.
from Bio import Phylo
iq = "iqtree3" if subprocess.run(["which", "iqtree3"], capture_output=True).returncode == 0 else "iqtree"

ALI_DOM = os.path.join(DATOS, "conjunto_dominio.faa")
cols = sorted(cols_dom)
with open(ALI_DOM, "w", encoding="utf-8") as fh:
    for k, v in ali.items():
        fh.write(">%s\n%s\n" % (k, "".join(v[c] for c in cols)))

TRIM = os.path.join(DATOS, "conjunto_trim.faa")
subprocess.run(["trimal", "-in", MSA, "-out", TRIM, "-gappyout"], check=True)
n_trim = len(next(iter(leer_faa(TRIM).values())))
print("  trimAl -gappyout sobre la completa: %d -> %d columnas" % (ancho, n_trim))

def inferir(etiqueta, alineamiento):
    pref = os.path.join(DATOS, "arbol_%s" % etiqueta)
    # KMP_AFFINITY=disabled: el runtime OpenMP de este contenedor aborta en
    # __kmp_middle_initialize al leer la mascara de afinidad. No cambia el resultado.
    subprocess.run([iq, "-s", alineamiento, "-m", "MFP", "-B", "1000", "--alrt", "1000",
                    "-o", "ORE1,NAP", "-T", "2", "--seed", SEMILLA, "--prefix", pref, "-redo", "--quiet"],
                   check=True, env=dict(os.environ, KMP_AFFINITY="disabled"))
    arbol = Phylo.read(pref + ".treefile", "newick")
    arbol.root_with_outgroup(["ORE1", "NAP"])
    d = dict(etiqueta=etiqueta, treefile=pref + ".treefile",
             modelo=re.search(r"Best-fit model.*?: (\S+)", open(pref + ".log").read()).group(1),
             n_cols=len(next(iter(leer_faa(alineamiento).values()))))

    def sop(c):   # IQ-TREE guarda "SH-aLRT/UFboot" como confianza del nodo interno
        return c.confidence if c.confidence is not None else c.name
    def minimo(hojas):
        return arbol.common_ancestor([arbol.find_any(name=h) for h in hojas])
    def hojas_de(c):
        return sorted(t.name for t in c.get_terminals())

    trio = minimo(["ANAC019", "ANAC055", "RD26"])
    d["trio_hojas"] = hojas_de(trio)
    d["trio_monofiletico"] = d["trio_hojas"] == ["ANAC019", "ANAC055", "RD26"]
    d["trio_soporte"] = sop(trio)
    # clado minimo que contiene al gen focal y a los tres de Arabidopsis
    con_focal = minimo(["ANAC019", "ANAC055", "RD26", GEN])
    d["clado_estres_hojas"] = hojas_de(con_focal)
    d["clado_estres_soporte"] = sop(con_focal)
    d["eucaliptos_en_clado_estres"] = [h for h in d["clado_estres_hojas"] if h.startswith("EglNAC")]
    # hermano del gen focal = el clado mas pequeno que lo contiene con alguien mas
    padre = min((c for c in arbol.get_nonterminals()
                 if GEN in hojas_de(c) and len(c.get_terminals()) > 1),
                key=lambda c: len(c.get_terminals()))
    d["hermano"] = [h for h in hojas_de(padre) if h != GEN]
    d["hermano_soporte"] = sop(padre)
    print("  [%s] %s, %d columnas" % (etiqueta, d["modelo"], d["n_cols"]))
    print("      trio monofiletico: %-5s (%s) soporte %s"
          % (d["trio_monofiletico"], ",".join(d["trio_hojas"]), d["trio_soporte"]))
    print("      clado de estres: %s  soporte %s"
          % (",".join(d["clado_estres_hojas"]), d["clado_estres_soporte"]))
    print("      hermano de %s: %s  soporte %s" % (GEN, ",".join(d["hermano"]), d["hermano_soporte"]))
    return d

arb = {a["etiqueta"]: a for a in (inferir("dominio", ALI_DOM), inferir("completa", TRIM))}
pri = arb["dominio"]          # el arbol que se reporta

# --- 5b. ¿rechazan los datos la topologia esperada? ------------------------
# La maxima verosimilitud mete al gen focal DENTRO del trio, lo que contradice la
# filogenia de especies: si las duplicaciones que dieron 019/055/072 ocurrieron en
# Brassicaceae, ningun gen de Eucalyptus puede quedar entre medio. Con soporte por
# debajo del umbral eso es ruido, pero conviene demostrarlo en vez de afirmarlo:
# se compara la topologia ML contra otra con el trio forzado a ser monofiletico
# (la que predice la filogenia de especies) mediante la prueba AU.
CONS = os.path.join(DATOS, "constrainte.nwk")
resto = [h for h in sorted(seqs) if h not in ("ANAC019", "ANAC055", "RD26")]
open(CONS, "w").write("((ANAC019,ANAC055,RD26),%s);\n" % ",".join(resto))
PREF_C = os.path.join(DATOS, "arbol_dominio_constr")
subprocess.run([iq, "-s", ALI_DOM, "-m", pri["modelo"], "-g", CONS, "-T", "2", "--seed", SEMILLA,
                "--prefix", PREF_C, "-redo", "--quiet"],
               check=True, env=dict(os.environ, KMP_AFFINITY="disabled"))
CAND = os.path.join(DATOS, "candidatas.nwk")
with open(CAND, "w") as fh:
    fh.write(open(pri["treefile"]).read().strip() + "\n")
    fh.write(open(PREF_C + ".treefile").read().strip() + "\n")
PREF_AU = os.path.join(DATOS, "au_test")
subprocess.run([iq, "-s", ALI_DOM, "-m", pri["modelo"], "-z", CAND, "-n", "0",
                "-zb", "100000", "-au", "-T", "2", "--seed", SEMILLA, "--prefix", PREF_AU, "-redo", "--quiet"],
               check=True, env=dict(os.environ, KMP_AFFINITY="disabled"))
# Tabla USER TREES: primer campo = indice del arbol, segundo = logL, ultimo
# numerico = p-AU. Se parsea por campos y no con un patron de decimales porque
# IQ-TREE imprime "1" y no "1.000" cuando la probabilidad es la unidad.
p_au = {}
for linea in open(PREF_AU + ".iqtree"):
    campos = [c for c in linea.split() if c not in ("+", "-")]
    if len(campos) >= 4 and campos[0] in ("1", "2") and campos[1].startswith("-"):
        p_au[int(campos[0])] = (float(campos[1]), float(campos[-1]))
assert len(p_au) == 2, "no pude leer la tabla AU en %s.iqtree" % PREF_AU
print("  AU: ML logL=%.2f p=%.3f | trio forzado logL=%.2f p=%.3f"
      % (p_au[1][0], p_au[1][1], p_au[2][0], p_au[2][1]))
pri["au_p_ML"], pri["au_p_constr"] = p_au[1][1], p_au[2][1]
with open(os.path.join(TABLAS, "T32_arboles_SNAC.csv"), "w", newline="", encoding="utf-8") as fh:
    campos = ["etiqueta", "modelo", "n_cols", "trio_monofiletico", "trio_soporte",
              "clado_estres_soporte", "hermano_soporte"]
    w = csv.writer(fh)
    w.writerow(campos + ["trio_hojas", "clado_estres_hojas", "hermano", "eucaliptos_en_clado_estres"])
    for a in arb.values():
        w.writerow([a[c] for c in campos] + [";".join(a[k]) for k in
                   ("trio_hojas", "clado_estres_hojas", "hermano", "eucaliptos_en_clado_estres")])
print("  T32: 2 arboles")

# --- 6. cifras para las leyendas ------------------------------------------
d30 = {r["ancla"]: r for r in filas30}
cif = [
    ("ort_n_secuencias", len(seqs), "conjunto.faa", "10 %s + 7 Arabidopsis + SNAC1" % SUBFAM),
    ("ort_modelo", pri["modelo"], "IQ-TREE ModelFinder", "arbol del dominio NAM"),
    ("ort_cols_dominio", pri["n_cols"], "T16 %d-%d aa" % (dom_ini, dom_fin), "columnas del dominio"),
    ("ort_cols_completa", n_trim, "trimAl -gappyout", "de %d columnas del L-INS-i" % ancho),
    ("ort_clado_estres", ", ".join(pri["clado_estres_hojas"]), pri["treefile"],
     "clado minimo con %s y el trio" % GEN),
    ("ort_clado_estres_soporte", pri["clado_estres_soporte"], pri["treefile"], "SH-aLRT/UFboot"),
    ("ort_eucaliptos_en_clado", ", ".join(pri["eucaliptos_en_clado_estres"]), pri["treefile"],
     "unicos Eucalyptus dentro del clado de estres"),
    ("ort_trio_monofiletico", pri["trio_monofiletico"], pri["treefile"],
     "el trio de Arabidopsis excluye a %s" % GEN),
    ("ort_hermano_focal", ", ".join(pri["hermano"]), pri["treefile"], "hermano de %s" % GEN),
    ("ort_hermano_soporte", pri["hermano_soporte"], pri["treefile"],
     "SH-aLRT/UFboot; <80/<95 = no resuelto"),
    ("ort_hermano_completa", ", ".join(arb["completa"]["hermano"]), arb["completa"]["treefile"],
     "sensibilidad: hermano en el arbol de proteina completa, soporte %s" % arb["completa"]["hermano_soporte"]),
    ("ort_au_p_ML", "%.3f" % pri["au_p_ML"], PREF_AU + ".iqtree", "prueba AU, topologia de maxima verosimilitud"),
    ("ort_au_p_trio", "%.3f" % pri["au_p_constr"], PREF_AU + ".iqtree",
     "prueba AU (100.000 replicas RELL); >0,05 = no se rechaza la coortologia 1:3. Es una estimacion Monte Carlo: cuatro corridas con -zb 10000 dieron 0,0570, 0,0590, 0,0620 y 0,0657, o sea justo en el limite. Se reporta como ~0,06, no como si fuera exacta."),
    ("ort_ident_dom_ANAC019", "%.1f %%" % d30["ANAC019"]["ident_dominio_pct"], "T30", METODO),
    ("ort_ident_dom_RD26", "%.1f %%" % d30["RD26"]["ident_dominio_pct"], "T30", METODO),
    ("ort_ident_comp_ANAC019", "%.1f %%" % d30["ANAC019"]["ident_completa_pct"], "T30", METODO),
    ("ort_ident_comp_RD26", "%.1f %%" % d30["RD26"]["ident_completa_pct"], "T30", METODO),
    ("ort_lectores_identicos_trio", "%d/%d" % (n_trio, len(filas31)), "T31",
     "posiciones que leen base identicas en ANAC019, ANAC055 y RD26"),
]
with open(os.path.join(CIFRAS, "ortologia.csv"), "w", newline="", encoding="utf-8") as fh:
    w = csv.writer(fh); w.writerow(["clave", "valor", "fuente", "nota"]); w.writerows(cif)
print("  %d cifras -> _cifras/ortologia.csv" % len(cif))
