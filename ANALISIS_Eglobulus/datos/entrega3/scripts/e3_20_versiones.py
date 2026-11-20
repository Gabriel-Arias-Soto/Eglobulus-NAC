#!/usr/bin/env python
"""Entrega 3, paso «Informe de validación y métodos»: tabla de versiones T61.

Una fila por programa, paquete o base de datos que usan los guiones e3_01-e3_18, con la versión
tal como quedó registrada al correr cada paso. No se teclea ninguna versión: cada una sale de

  - datos/entrega3/entorno_nac3_paquetes.txt      lista exportada del entorno nac3 (micromamba);
  - los resumen_*.json de datos/entrega3            versiones registradas en tiempo de ejecución (HMMER y
                                                     el perfil Pfam en JCVI, FIMO, UniProt, GO, JASPAR);
  - los logs con sessionInfo()                      R y paquetes del entorno r-nac (figuras);
  - datos/entrega3/validacion/resumen_controles.json  ensamblados y anotaciones de NCBI (e3_19).

Salida: tablas/T61_versiones_entrega3.csv
  componente, tipo (programa, paquete de Python, paquete de R, base de datos), version, entorno,
  pasos (guiones e3_ que lo usan), fuente (de dónde sale la versión), cita (clave BibTeX o vacío)
Uso (python de nac3): python e3_20_versiones.py
"""
import ast
import datetime
import glob
import json
import os
import re
import sys

import pandas as pd


def _raiz():
    d = os.path.dirname(os.path.abspath(__file__))
    while True:
        if os.path.basename(d) == "ANALISIS_Eglobulus":
            return d
        if os.path.isdir(os.path.join(d, "ANALISIS_Eglobulus")):
            return os.path.join(d, "ANALISIS_Eglobulus")
        if os.path.dirname(d) == d:
            sys.exit("No se encuentra la carpeta ANALISIS_Eglobulus por encima del guion")
        d = os.path.dirname(d)


A = _raiz()
E3 = os.path.join(A, "datos", "entrega3")
LOG = os.path.join(A, "logs")


def js(f):
    return json.load(open(os.path.join(E3, f), encoding="utf-8"))


# --- entorno nac3 --------------------------------------------------------------------------------
fpaq = os.path.join(E3, "entorno_nac3_paquetes.txt")
PAQ = {}
for l in open(fpaq, encoding="utf-8"):
    c = l.split()
    if len(c) >= 4 and c[0] not in ("List", "Name") and not c[0].startswith("─"):
        PAQ[c[0]] = c[1]
FPAQ = "datos/entrega3/entorno_nac3_paquetes.txt"
FILAS = []


def fila(componente, tipo, version, entorno, pasos, fuente, cita=""):
    assert version not in (None, ""), componente
    FILAS.append(dict(componente=componente, tipo=tipo, version=str(version), entorno=entorno,
                      pasos=pasos, fuente=fuente, cita=cita))


# Pasos que usan cada componente, leídos de los propios guiones: import/library para los paquetes y
# el nombre del ejecutable para los programas externos.
SCR = os.path.join(E3, "scripts")
def _codigo(f):
    """Código sin docstring de módulo ni comentarios, para no contar menciones en el texto."""
    t = open(f, encoding="utf-8").read()
    if f.endswith(".py"):
        doc = ast.get_docstring(ast.parse(t), clean=False)
        if doc:
            t = t.replace(doc, "", 1)
        return "\n".join(l for l in t.splitlines() if not l.lstrip().startswith("#"))
    return "\n".join(re.sub(r"(^|\s)#.*$", "", l) for l in t.splitlines())


# e3_19 y e3_20 (controles y versiones) no analizan nada: no cuentan como uso
FUENTES = {os.path.basename(f)[:5]: _codigo(f)
           for f in sorted(glob.glob(os.path.join(SCR, "e3_*.py")) + glob.glob(os.path.join(SCR, "e3_*.R")))
           if os.path.basename(f)[:5] not in ("e3_19", "e3_20")}


def usan(patron):
    return [k for k, t in FUENTES.items() if re.search(patron, t, re.M)]


def rango(ps):
    return ", ".join(ps) if ps else "(ningún guion e3_)"


def nac3(paquete, componente, tipo, patron, cita="", nota=""):
    ps = usan(patron)
    fila(componente, tipo, PAQ[paquete], "nac3", rango(ps) + nota, f"{FPAQ} ({paquete})", cita)


PY_IMP = lambda m: rf"^\s*(import|from)\s+{m}\b"
R_LIB = lambda m: rf"library\({re.escape(m)}\)"
nac3("python", "Python", "programa", r"^#!/usr/bin/env python|^\s*import ")
nac3("diamond", "DIAMOND", "programa", r"/diamond[\"']", "Buchfink2021")
nac3("mcscanx", "MCScanX", "programa", r"/(MCScanX|duplicate_gene_classifier)[\"']", "Wang2012")
nac3("mafft", "MAFFT", "programa", r"/mafft[\"']", "Katoh2013")
nac3("paml", "PAML (yn00, codeml)", "programa", r"/(yn00|codeml)[\"']", "Yang2007; Yang2000")
nac3("last", "LAST", "programa", r"jcvi\.compara\.catalog", "Kielbasa2011")
nac3("jcvi", "JCVI (MCscan)", "paquete de Python", r"jcvi\.(compara|graphics)", "Tang2024")
nac3("biopython", "Biopython", "paquete de Python", PY_IMP("Bio"), "Cock2009")
nac3("numpy", "NumPy", "paquete de Python", PY_IMP("numpy"))
nac3("pandas", "pandas", "paquete de Python", PY_IMP("pandas"))
nac3("scipy", "SciPy", "paquete de Python", PY_IMP("scipy"), "Virtanen2020")
nac3("statsmodels", "statsmodels", "paquete de Python", PY_IMP("statsmodels"), "Seabold2010")
nac3("matplotlib-base", "Matplotlib", "paquete de Python", PY_IMP("matplotlib"))
nac3("pycirclize", "pyCirclize", "paquete de Python", PY_IMP("pycirclize"))
nac3("pygenomeviz", "pyGenomeViz", "paquete de Python", PY_IMP("pygenomeviz"))
nac3("ete4", "ete4", "paquete de Python", r"jcvi\.graphics\.karyotype", nota=" (dependencia de jcvi.graphics.karyotype)")
nac3("r-base", "R (nac3)", "programa", R_LIB("clusterProfiler"), "RCore2026")
nac3("bioconductor-clusterprofiler", "clusterProfiler", "paquete de R", R_LIB("clusterProfiler"), "Yu2012; Wu2021")
nac3("bioconductor-enrichplot", "enrichplot", "paquete de R", R_LIB("enrichplot") + r"|enrichplot::")
nac3("bioconductor-gosemsim", "GOSemSim", "paquete de R", R_LIB("GOSemSim") + r"|GOSemSim::", "Yu2010; Wang2007")
nac3("bioconductor-org.at.tair.db", "org.At.tair.db", "paquete de R", R_LIB("org.At.tair.db"))
nac3("bioconductor-go.db", "GO.db", "paquete de R", R_LIB("GO.db"))
nac3("bioconductor-annotationdbi", "AnnotationDbi", "paquete de R", R_LIB("AnnotationDbi"))

# --- versiones registradas al correr ----------------------------------------------------------------
INT = "colinealidad/resumen_colinealidad_interespecies.json"
vi = js(INT)["metodo"]["versiones"]
fila("HMMER (hmmsearch)", "programa", re.search(r"HMMER ([\d.]+)", vi["hmmer"]).group(1), "hmmer", rango(usan(r"/hmmsearch[\"']")),
     f"datos/entrega3/{INT} (metodo.versiones.hmmer)", "Eddy2011")
assert vi["jcvi"].split()[-1] == PAQ["jcvi"] and vi["last"].split()[-1] == PAQ["last"], vi
prep = js("cis_jaspar/preparacion_escaneo.json")
fila("MEME Suite (FIMO, fasta-get-markov)", "programa", prep["fimo"], "rnaseq", rango(usan(r"MEMEBIN, [\"'](fimo|fasta-get-markov)[\"']|/(fimo|fasta-get-markov)[\"']")),
     "datos/entrega3/cis_jaspar/preparacion_escaneo.json (fimo)", "Grant2011; Bailey2015")

# --- R del entorno r-nac (figuras), desde sessionInfo() -------------------------------------------
lg = open(os.path.join(LOG, "entrega3_09_go_figuras.log"), encoding="utf-8").read()
rv = re.search(r"^R version ([\d.]+)", lg, re.M).group(1)
fila("R (r-nac)", "programa", rv, "r-nac", rango(usan(R_LIB("patchwork"))), "logs/entrega3_09_go_figuras.log (sessionInfo)", "RCore2026")
for p, cita in (("ggplot2", "Wickham2016"), ("ggtree", ""), ("patchwork", ""), ("ape", "")):
    m = re.search(rf"\b{p}_([\d.\-]+)", lg)
    fila(p, "paquete de R", m.group(1), "r-nac", rango([k for k in usan(R_LIB(p) + r"|\bggplot\(|tidyverse" if p == "ggplot2" else R_LIB(p)) if k in ("e3_08", "e3_18")]),
         "logs/entrega3_09_go_figuras.log (sessionInfo)", cita)
lg5 = open(os.path.join(LOG, "entrega3_06_go.log"), encoding="utf-8").read()
for p in ("clusterProfiler", "enrichplot", "GOSemSim", "org.At.tair.db", "GO.db"):
    v = re.search(rf"\b{re.escape(p)}_([\d.\-]+)", lg5).group(1)
    ref = next(f for f in FILAS if f["componente"] == p)
    assert v == ref["version"].split("-")[0], (p, v, ref["version"])

# --- bases de datos y servicios ---------------------------------------------------------------------
rc = json.load(open(os.path.join(E3, "validacion", "resumen_controles.json"), encoding="utf-8"))
CITA_GEN = {"Egl": "Li_2025a", "Egr": "Myburg2014", "Ptr": "Tuskan2006", "Ath": "Lamesch2012; Cheng2017"}
NOMBRE = {"Egl": "E. globulus", "Egr": "E. grandis", "Ptr": "P. trichocarpa", "Ath": "A. thaliana"}
for sp in ("Egl", "Egr", "Ptr", "Ath"):
    r = rc["referencias"][sp]
    fila(f"Genoma de {NOMBRE[sp]} ({r['accesion']})", "base de datos",
         f"{r['ensamblado']}; anotación «{r['anotacion']}» ({r['fecha_anotacion']})", "NCBI Datasets",
         rango(usan(re.escape(r["accesion"].split(".")[0]) + "|" + {"Egl": "Eglobulus_ASM", "Egr": "Egrandis_v2", "Ptr": "Ptrichocarpa_v4", "Ath": "Athaliana_TAIR10"}[sp])),
         "datos/entrega3/validacion/resumen_controles.json (referencias)",
         CITA_GEN[sp] + ("; OLeary2024" if sp != "Egl" else ""))
uni = js("go/uniprot_NAC_Ath_UP000006548.version.json")
fila("UniProtKB (proteoma UP000006548)", "base de datos", f"release {uni['cabeceras']['X-UniProt-Release']} "
     f"(descarga {uni['fecha_descarga']})", "API REST de UniProt", rango(usan(r"UP000006548")),
     "datos/entrega3/go/uniprot_NAC_Ath_UP000006548.version.json", "UniProt2025")
go = js("go/resumen_go.json")["fuentes"]
fila("Gene Ontology (GO.db) y anotación GO de TAIR (org.At.tair.db)", "base de datos",
     f"GO {go['GO_fecha']}; TAIR {go['TAIR_GO_fecha'].split(';')[0]}", "nac3", rango(usan(R_LIB("org.At.tair.db"))),
     "datos/entrega3/go/resumen_go.json (fuentes)", "Lamesch2012")
fila("Pfam NAM (PF02365)", "base de datos", vi["perfil"].split()[1] + " (umbral GA " + vi["perfil"].split()[3] + " bits)",
     "referencia/PF02365.hmm", rango(usan(r"PF02365\.hmm")), f"datos/entrega3/{INT} (metodo.versiones.perfil)", "Mistry2021")
fila("JASPAR CORE plants non-redundant", "base de datos",
     f"{prep['jaspar'].split(' CORE')[0]} (descarga {prep['fecha_descarga_jaspar'][:10]})", "descarga local", rango(usan(r"JASPAR2026_CORE")),
     "datos/entrega3/cis_jaspar/preparacion_escaneo.json (jaspar)", "JASPAR2026")
res = sorted(glob.glob(os.path.join(E3, "cis_plantcare", "resultados", "PlantCARE_*_plantCARE.tar.gz")))
fechas = sorted(datetime.datetime.fromtimestamp(os.path.getmtime(f)).strftime("%Y-%m-%d") for f in res)
lec = js("cis_plantcare/resumen_lectura_PlantCARE.json")
fila("PlantCARE", "servicio web", f"resultados de {len(res)} lotes guardados el {fechas[0]}"
     + (f" a {fechas[-1]}" if fechas[-1] != fechas[0] else "") + f"; leídos el {lec['fecha'][:10]}",
     "bioinformatics.psb.ugent.be", rango(usan(r"plantcare/html|plantCARE\.tar\.gz")), "datos/entrega3/cis_plantcare/resultados (fechas de los .tar.gz)",
     "Lescot2002")

t61 = pd.DataFrame(FILAS)
assert t61.componente.is_unique
t61.to_csv(os.path.join(A, "tablas", "T61_versiones_entrega3.csv"), index=False, encoding="utf-8")
print(f"T61: {len(t61)} componentes; {t61.tipo.value_counts().to_dict()}")
