#!/usr/bin/env python
"""Entrega 3, paso «Informe de validación y métodos»: tabla de controles T60.

Reúne en una tabla cada control del plan de la entrega 3 con su criterio, la cifra observada y el
resultado. Las cifras se leen de lo que dejaron los guiones e3_01-e3_18 (resumen_*.json de
datos/entrega3, tablas T40-T59, figuras F13-F19, logs) y de comprobaciones propias de este guion:

  - md5 de los ficheros de los paquetes de NCBI Datasets (E. grandis, P. trichocarpa, A. thaliana)
    frente al md5sum.txt de cada paquete: la copia extraída que leen los guiones, si existe, y si no
    el miembro del .zip;
  - md5 de los ficheros de GCA_046226855.1 (E. globulus) que leen los guiones frente al
    md5checksums.txt de NCBI, y genes codificantes frente al informe de NCBI Datasets (ambos en
    datos/entrega3/validacion, descargados para este paso);
  - genes y proteínas por especie;
  - figuras: PNG a 600 ppp, 174 mm de ancho y <= 234 mm de alto, página del PDF del mismo tamaño y
    solo fuentes Liberation Sans incrustadas;
  - filas de las tablas frente a las cifras de los resúmenes.

Resultado: PASA (se cumple el criterio), FALLA (no se cumple) o INFORMATIVO (cifra descriptiva sin
criterio de aceptación fijado de antemano; se informa porque hace falta para leer otros controles).

Calibraciones (falsos «NAC» sacados del fondo): la fracción de p < 0,05 no debe superar 0,05 (prueba
binomial unilateral de exceso con todas las pruebas de los conjuntos de calibración; las pruebas de un
mismo conjunto están correlacionadas, así que el criterio es más exigente que el nominal) y la fracción
de conjuntos con algún q < 0,05 tampoco (con la hipótesis nula global, BH acota esa fracción en 0,05).

Barreras: si falta una clave de un resumen, el guion para (KeyError); si algún booleano de los
resumen_*.json no lo cita ningún control, el guion para y lo lista.

Salidas
  tablas/T60_controles_entrega3.csv                  id, paso, control, criterio, observado, resultado, fuente
  datos/entrega3/validacion/resumen_controles.json   recuentos, cifras propias y versiones
  datos/entrega3/validacion/md5_referencias.json     md5 por fichero (también es la caché)
Uso (python de nac3, desde cualquier carpeta):
  python e3_19_controles.py > ANALISIS_Eglobulus/logs/entrega3_21_controles.log 2>&1
"""
import datetime
import glob
import gzip
import hashlib
import json
import os
import platform
import re
import sys
import zipfile
import zlib

import pandas as pd
import PIL
from PIL import Image
import scipy
from scipy.stats import binom


# --- rutas -------------------------------------------------------------------------------------
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
BASE = os.path.dirname(A)                      # Revisión NAC
E3 = os.path.join(A, "datos", "entrega3")
TAB = os.path.join(A, "tablas")
FIG = os.path.join(A, "figuras")
LOG = os.path.join(A, "logs")
REF = os.path.join(A, "referencia")
VAL = os.path.join(E3, "validacion")
RQ = os.path.join(BASE, "rnaseq_udec", "ref")
os.makedirs(VAL, exist_ok=True)

PASOS = {
    "01": "Crear entorno nac3 en el disco Z",
    "02": "Descargar referencias de NCBI a Z",
    "03": "Elegir una isoforma por gen",
    "04": "Extraer promotores de 2 kb",
    "05": "Enviar 4 lotes a PlantCARE",
    "06": "Asignar homólogos de Arabidopsis",
    "07": "Enriquecimiento GO con clusterProfiler",
    "08": "Figuras y tablas de GO",
    "09": "Colinealidad dentro de E. globulus",
    "10": "Ka/Ks y datación de duplicados",
    "11": "Colinealidad con E. grandis y P. trichocarpa",
    "12": "Figuras de colinealidad",
    "13": "Leer y verificar la salida de PlantCARE",
    "14": "Contraste con fondo emparejado por GC",
    "14J": "Análisis cis local con JASPAR 2026 y FIMO (añadido al plan)",
    "15": "Figuras y tablas de elementos cis",
}

# --- resúmenes JSON: accesor que registra lo que se usa ------------------------------------------
RES = sorted(os.path.relpath(p, E3) for p in glob.glob(os.path.join(E3, "**", "*.json"), recursive=True)
             if not os.path.relpath(p, E3).startswith(("validacion", "scripts")))
JS = {f: json.load(open(os.path.join(E3, f), encoding="utf-8")) for f in RES}
USADAS = set()


def _hojas(x, ruta=()):
    if isinstance(x, dict):
        for k, v in x.items():
            yield from _hojas(v, ruta + (k,))
    elif isinstance(x, list):
        for i, v in enumerate(x):
            yield from _hojas(v, ruta + (i,))
    else:
        yield ruta, x


def J(f, *ruta):
    """Valor de un resumen; KeyError si falta. Marca como usadas todas las hojas leídas."""
    x = JS[f]
    for p in ruta:
        x = x[p]
    for r, _ in _hojas(x, tuple(ruta)):
        USADAS.add((f, r))
    return x


def fu(f, *ruta):
    return f"{f}:" + ".".join(str(p) for p in ruta)


ISO = "isoformas/resumen_isoformas.json"
PROM = "promotores/resumen_promotores.json"
HOM = "go/resumen_homologos.json"
UNI = "go/uniprot_NAC_Ath_UP000006548.version.json"
GO = "go/resumen_go.json"
REC = "go/control_recalculo_go.json"
GOF = "go/resumen_go_figuras.json"
COL = "colinealidad/resumen_colinealidad_Egl.json"
KK = "colinealidad/resumen_kaks_Egl.json"
INT = "colinealidad/resumen_colinealidad_interespecies.json"
LEC = "cis_plantcare/resumen_lectura_PlantCARE.json"
FGC = "cis_plantcare/fondo/resumen_fondo_GC.json"
FGC10 = "cis_plantcare/fondo/resumen_fondo_GC_10000.json"
CPC = "cis_plantcare/resumen_contraste_PlantCARE.json"
JPREP = "cis_jaspar/preparacion_escaneo.json"
JESC = "cis_jaspar/resumen_escaneo_JASPAR2026.json"
JCON = "cis_jaspar/resumen_contraste_JASPAR2026.json"

# --- formato (texto en español, como mil() de metodos_bioinformatica.Rmd: coma decimal, punto de miles) ---


def n(x):
    x = int(x)
    s = f"{abs(x):,d}".replace(",", ".")
    return ("-" if x < 0 else "") + s


def d(x, k=2):
    return f"{x:.{k}f}".replace(".", ",")


def pc(a, b, k=1):
    return d(100 * a / b, k) + " %"


def ec(p, k=1):
    m, e = f"{p:.{k}e}".split("e")
    return f"{m.replace('.', ',')}e{int(e)}"


def es_v(s):
    return s.astype(str).str.lower().isin(["true", "1"])


def tabla(nombre):
    f = glob.glob(os.path.join(TAB, nombre + "_*.csv"))
    assert len(f) == 1, (nombre, f)
    return pd.read_csv(f[0]), "tablas/" + os.path.basename(f[0])


# --- filas de T60 --------------------------------------------------------------------------------
FILAS = []


def ctl(paso, control, criterio, observado, resultado, fuentes):
    assert paso in PASOS, paso
    assert resultado in ("PASA", "FALLA", "INFORMATIVO"), resultado
    assert resultado == "INFORMATIVO" or criterio, control
    k = sum(1 for f in FILAS if f["paso_id"] == paso) + 1
    FILAS.append(dict(id=f"C{paso}.{k}", paso_id=paso, paso=PASOS[paso], control=control,
                      criterio=criterio, observado=observado, resultado=resultado,
                      fuente="; ".join(fuentes)))


def ok(c):
    return "PASA" if c else "FALLA"


PROPIO = {}          # cifras calculadas aquí, para resumen_controles.json

# --- md5 con caché ------------------------------------------------------------------------------
CACHE_F = os.path.join(VAL, "md5_referencias.json")
CACHE = json.load(open(CACHE_F, encoding="utf-8")) if os.path.exists(CACHE_F) else {}


def _guardar_cache():
    tmp = CACHE_F + ".tmp"
    json.dump(CACHE, open(tmp, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    os.replace(tmp, CACHE_F)


def _md5_flujo(fh, bloque=1 << 22):
    h = hashlib.md5()
    for b in iter(lambda: fh.read(bloque), b""):
        h.update(b)
    return h.hexdigest()


def _md5(clave, st, abrir):
    c = CACHE.get(clave)
    if c and c["bytes"] == st.st_size and c["mtime"] == int(st.st_mtime):
        return c["md5"]
    with abrir() as fh:
        m = _md5_flujo(fh)
    CACHE[clave] = dict(bytes=st.st_size, mtime=int(st.st_mtime), md5=m)
    _guardar_cache()
    print(f"  md5 {clave}: {m}", flush=True)
    return m


def md5_fichero(p):
    return _md5(os.path.relpath(p, BASE), os.stat(p), lambda: open(p, "rb"))


def md5_miembro(zp, z, nombre):
    return _md5(os.path.relpath(zp, BASE) + "::" + nombre, os.stat(zp), lambda: z.open(nombre))


def n_fasta(p):
    with (gzip.open(p, "rb") if p.endswith(".gz") else open(p, "rb")) as fh:
        return sum(1 for l in fh if l.startswith(b">"))


# =================================================================================================
print(f"e3_19_controles.py  {datetime.datetime.now():%Y-%m-%d %H:%M:%S}  raíz: {A}", flush=True)
print(f"{len(RES)} resúmenes JSON", flush=True)

# --- 01 entorno nac3 --------------------------------------------------------------------------------
fpaq = os.path.join(E3, "entorno_nac3_paquetes.txt")
lin = open(fpaq, encoding="utf-8").read().splitlines()
prefijo = re.search(r'environment: "([^"]+)"', lin[0]).group(1)
PAQ = {}
for l in lin:
    c = l.split()
    if len(c) >= 4 and c[0] not in ("List", "Name") and not c[0].startswith("─"):
        PAQ[c[0]] = c[1]
cache_mamba = os.path.join(os.path.dirname(os.path.dirname(prefijo)), "pkgs")
ctl("01", "Prefijo del entorno y caché de paquetes en el disco Z",
    "prefijo del entorno y caché de paquetes bajo /mnt/z (disco Z)",
    f"prefijo {prefijo}; caché {cache_mamba} ({'existe' if os.path.isdir(cache_mamba) else 'no existe'})",
    ok(prefijo.startswith("/mnt/z/") and os.path.isdir(cache_mamba)),
    [os.path.relpath(fpaq, A)])
PLAN = {"diamond": "2.2.8", "jcvi": "1.6.7", "pycirclize": "1.10", "pygenomeviz": "1.7", "biopython": "1.88",
        "bioconductor-clusterprofiler": "4.18", "bioconductor-enrichplot": "1.30",
        "bioconductor-org.at.tair.db": "3.22", "bioconductor-complexheatmap": "2.26",
        "mcscanx": "", "last": "", "mafft": ""}
vers_ok = {k: (k in PAQ and (PAQ[k] == v or PAQ[k].startswith(v + ".") or v == "")) for k, v in PLAN.items()}
ctl("01", "Versiones de los programas del plan en nac3",
    "DIAMOND 2.2.8, JCVI 1.6.7, pyCirclize 1.10, pyGenomeViz 1.7, Biopython 1.88, clusterProfiler 4.18, "
    "enrichplot 1.30, org.At.tair.db 3.22, ComplexHeatmap 2.26; MCScanX, LAST y MAFFT instalados",
    "; ".join(f"{k.replace('bioconductor-', '')} {PAQ.get(k, 'ausente')}" for k in PLAN),
    ok(all(vers_ok.values())), [os.path.relpath(fpaq, A)])
fyml = os.path.join(E3, "entorno_nac3.yml")
ctl("01", "Lista de versiones exportada",
    "entorno_nac3.yml y lista de paquetes en datos/entrega3",
    f"{len(PAQ)} paquetes en la lista; yml {'presente' if os.path.exists(fyml) else 'ausente'}",
    ok(os.path.exists(fyml) and len(PAQ) > 0), [os.path.relpath(fpaq, A), os.path.relpath(fyml, A)])
PROPIO["nac3"] = dict(prefijo=prefijo, paquetes=len(PAQ),
                      versiones={k: PAQ.get(k) for k in sorted(set(PLAN) | {
                          "python", "r-base", "paml", "numpy", "pandas", "scipy", "matplotlib-base",
                          "statsmodels", "bioconductor-go.db", "bioconductor-gosemsim", "bioconductor-dose",
                          "bioconductor-annotationdbi", "ete4", "r-ggplot2", "bedtools", "r-jsonlite"})})

# --- 02 referencias de NCBI -------------------------------------------------------------------------
NCBI = {"Egr": ("Egrandis_v2", ["GCA_000612305.2_ncbi_dataset.zip"]),
        "Ptr": ("Ptrichocarpa_v4.1", ["GCF_000002775.5_ncbi_dataset.zip"]),
        "Ath": ("Athaliana_TAIR10.1", ["GCF_000001735.4_ncbi_dataset.zip", "GCF_000001735.4_cds.zip"])}
NOMBRE = {"Egl": "E. globulus", "Egr": "E. grandis", "Ptr": "P. trichocarpa", "Ath": "A. thaliana"}
PROPIO["referencias"] = {}
for sp, (carpeta, zips) in NCBI.items():
    dsp = os.path.join(REF, carpeta)
    filas_md5, informe = [], None
    for zn in zips:
        zp = os.path.join(dsp, zn)
        z = zipfile.ZipFile(zp)
        for l in z.read("ncbi_dataset/data/md5sum.txt" if "ncbi_dataset/data/md5sum.txt" in z.namelist()
                        else "md5sum.txt").decode().splitlines():
            if not l.strip():
                continue
            esperado, miembro = l.split(None, 1)
            miembro = miembro.strip()
            ext = os.path.join(dsp, os.path.basename(miembro))
            if os.path.exists(ext):
                obs, origen = md5_fichero(ext), "extraído"
            else:
                obs, origen = md5_miembro(zp, z, miembro), "zip"
            filas_md5.append(dict(zip=zn, fichero=os.path.basename(miembro), origen=origen,
                                  esperado=esperado, observado=obs, igual=obs == esperado))
            if miembro.endswith("assembly_data_report.jsonl") and informe is None:
                informe = json.loads(z.read(miembro).decode().splitlines()[0])
    n_ig = sum(r["igual"] for r in filas_md5)
    n_ext = sum(r["origen"] == "extraído" for r in filas_md5)
    ctl("02", f"md5 del paquete NCBI Datasets de {NOMBRE[sp]} ({informe['accession']})",
        "md5 de cada fichero = md5sum.txt del paquete (copia extraída si existe; si no, miembro del zip)",
        f"{n_ig}/{len(filas_md5)} iguales ({n_ext} extraídos, {len(filas_md5) - n_ext} dentro del zip)",
        ok(n_ig == len(filas_md5)), [f"referencia/{carpeta}/{z_}:md5sum.txt" for z_ in zips])
    gc = informe["annotationInfo"]["stats"]["geneCounts"]
    prot = n_fasta(os.path.join(dsp, "protein.faa"))
    gp = J(ISO, sp, "genes_con_proteina")
    ctl("02", f"Genes y proteínas de {NOMBRE[sp]}",
        "genes codificantes del informe de NCBI = genes con proteína leídos del GFF y protein.faa",
        f"NCBI {n(gc['proteinCoding'])} codificantes ({n(gc['total'])} genes en total); con proteína "
        f"{n(gp)}; {n(prot)} proteínas en protein.faa",
        ok(gc["proteinCoding"] == gp), [fu(ISO, sp, "genes_con_proteina"),
                                        f"referencia/{carpeta}/{zips[0]}:assembly_data_report.jsonl"])
    ai = informe["annotationInfo"]
    PROPIO["referencias"][sp] = dict(
        accesion=informe["accession"], ensamblado=informe["assemblyInfo"]["assemblyName"],
        anotacion=ai.get("name"), fecha_anotacion=ai.get("releaseDate"),
        genes_codificantes=gc["proteinCoding"], genes_total=gc["total"], proteinas=prot,
        descarga=datetime.datetime.fromtimestamp(os.path.getmtime(os.path.join(dsp, zips[0]))).strftime("%Y-%m-%d"),
        md5=filas_md5)
# E. globulus (GCA_046226855.1): ficheros que leen e3_01, e3_02 y e3_04
md5ncbi = {}
for l in open(os.path.join(VAL, "md5checksums_GCA_046226855.1.txt"), encoding="utf-8"):
    if l.strip():
        h, f = l.split(None, 1)
        md5ncbi[f.strip().lstrip("./")] = h
egl_f = [f"GCA_046226855.1_ASM4622685v1_{s}" for s in
         ("genomic.fna.gz", "protein.faa.gz", "cds_from_genomic.fna.gz", "genomic.gtf.gz")]
egl_ok = {f: md5_fichero(os.path.join(RQ, f)) == md5ncbi[f] for f in egl_f}
fecha_md5 = open(os.path.join(VAL, "md5checksums_GCA_046226855.1.fecha.txt")).read().strip()[:10]
ctl("02", "md5 de los ficheros de E. globulus GCA_046226855.1 que leen los guiones",
    "md5 = md5checksums.txt publicado por NCBI (genoma, proteínas, CDS y GTF)",
    f"{sum(egl_ok.values())}/{len(egl_ok)} iguales (md5checksums.txt descargado el {fecha_md5})",
    ok(all(egl_ok.values())), ["datos/entrega3/validacion/md5checksums_GCA_046226855.1.txt"])
rep = json.load(open(os.path.join(VAL, "dataset_report_GCA_046226855.1.json"), encoding="utf-8"))["reports"][0]
gce = rep["annotation_info"]["stats"]["gene_counts"]
prot_e = n_fasta(os.path.join(RQ, egl_f[1]))
gp = J(ISO, "Egl", "genes_con_proteina")
ctl("02", "Genes y proteínas de E. globulus",
    "genes codificantes del informe de NCBI Datasets = genes con proteína leídos del GTF y protein.faa",
    f"NCBI {n(gce['protein_coding'])} codificantes; con proteína {n(gp)}; {n(prot_e)} proteínas en protein.faa",
    ok(gce["protein_coding"] == gp), [fu(ISO, "Egl", "genes_con_proteina"),
                                      "datos/entrega3/validacion/dataset_report_GCA_046226855.1.json"])
PROPIO["referencias"]["Egl"] = dict(
    accesion=rep["accession"], ensamblado=rep["assembly_info"]["assembly_name"],
    anotacion=rep["annotation_info"].get("name"), fecha_anotacion=rep["annotation_info"].get("release_date"),
    genes_codificantes=gce["protein_coding"], proteinas=prot_e, md5_iguales=egl_ok, md5_fecha=fecha_md5)

# --- 03 isoforma representativa -------------------------------------------------------------------
t01 = pd.read_csv(os.path.join(TAB, "T01_censo_maestro_189.csv"), dtype=str, keep_default_na=False)
k189 = ["encontrados", "proteina_igual_T01", "largo_aa_igual_T01", "traduccion_ok", "empiezan_ATG", "en_cromosomas"]
v189 = {k: J(ISO, "Egl", "NAC_189", k) for k in k189}
ctl("03", "Los 189 NAC: proteína fijada en T01, CDS que la traduce, ATG inicial y cromosoma",
    "189/189 en cada comprobación (identificador y secuencia = T01)",
    "; ".join(f"{k.replace('_', ' ')} {v}" for k, v in v189.items()),
    ok(all(v == 189 for v in v189.values()) and len(t01) == 189),
    [fu(ISO, "Egl", "NAC_189", k) for k in k189] + ["tablas/T01_censo_maestro_189.csv"])
rep_egl = pd.read_csv(os.path.join(E3, "isoformas", "Egl_representativas.tsv"), sep="\t")
parc = rep_egl[rep_egl.gen.isin(t01.locus) & es_v(rep_egl.parcial)].merge(t01[["locus", "simbolo"]],
                                                                          left_on="gen", right_on="locus")
assert len(parc) == J(ISO, "Egl", "NAC_189", "parciales"), parc
ctl("03", "NAC con CDS parcial", "",
    "; ".join(f"{r.simbolo} ({r.gen}, parcial en {r.parcial_tipo}, {r.largo_aa} aa)" for r in parc.itertuples()),
    "INFORMATIVO", [fu(ISO, "Egl", "NAC_189", "parciales"), "datos/entrega3/isoformas/Egl_representativas.tsv"])
for sp in ("Egl", "Egr", "Ptr", "Ath"):
    g = {k: J(ISO, sp, k) for k in ("genes_con_proteina", "en_cromosomas", "cromosomas", "traduccion_ok",
                                    "traduccion_falla", "inicio_alternativo", "con_paradas_internas",
                                    "con_selenocisteina", "completas_5p", "completas_5p_con_ATG", "parciales",
                                    "genes_con_varias_isoformas")}
    rsp = pd.read_csv(os.path.join(E3, "isoformas", f"{sp}_representativas.tsv"), sep="\t")
    tok = es_v(rsp.traduccion_ok)
    par = rsp.paradas_internas.fillna(0).astype(float) > 0
    alt = es_v(rsp.inicio_alternativo)
    assert len(rsp) == g["genes_con_proteina"] and (~tok).sum() == g["traduccion_falla"] and par.sum() == g["con_paradas_internas"]
    ctl("03", f"Genoma de {NOMBRE[sp]}: traducción de la CDS representativa", "",
        f"{n(g['traduccion_ok'])}/{n(g['genes_con_proteina'])} traducen a la proteína "
        f"({pc(g['traduccion_ok'], g['genes_con_proteina'], 2)}); no traducen {g['traduccion_falla']} "
        f"({(~tok & par).sum()} con paradas internas, {(~tok & ~par).sum()} por CDS parcial o marco); "
        f"{(tok & alt).sum()} traducen con un codón de inicio alternativo; completas en 5' con ATG "
        f"{n(g['completas_5p_con_ATG'])}; parciales {g['parciales']}; {n(g['en_cromosomas'])} en {g['cromosomas']} "
        f"cromosomas; {n(g['genes_con_varias_isoformas'])} con varias isoformas",
        "INFORMATIVO", [fu(ISO, sp, k) for k in g] + [f"datos/entrega3/isoformas/{sp}_representativas.tsv"])

# --- 04 promotores -------------------------------------------------------------------------------
t40, f40 = tabla("T40")
nac = {k: J(PROM, "NAC", k) for k in ("n", "ATG_verificado", "largo_2000", "con_N", "max_pct_N", "atg_dentro_del_gen")}
t40_ok = (len(t40) == 189 and es_v(t40.atg_verificado).all() and (t40.largo == 2000).all()
          and (t40.pct_N == 0).all() and es_v(t40.atg_dentro_del_gen).all())
ctl("04", "189 promotores NAC: ATG tras la ventana, 2000 pb, sin N, ATG dentro del gen",
    "189/189 en cada comprobación, en el resumen y en T40",
    f"n {nac['n']}; ATG {nac['ATG_verificado']}; 2000 pb {nac['largo_2000']}; con N {nac['con_N']} "
    f"(máximo {d(nac['max_pct_N'], 1)} %); ATG dentro del gen {nac['atg_dentro_del_gen']}; T40 {len(t40)} filas",
    ok(all(nac[k] == 189 for k in ("n", "ATG_verificado", "largo_2000", "atg_dentro_del_gen"))
       and nac["con_N"] == 0 and t40_ok),
    [fu(PROM, "NAC", k) for k in nac] + [f40])
gen = {k: J(PROM, k) for k in ("genes_en_cromosomas", "promotores", "truncados", "con_N", "ATG_verificado")}
ctl("04", "Genoma: promotores sin truncar y sin N",
    "ningún promotor truncado en el extremo del cromosoma ni con N",
    f"{n(gen['promotores'])} promotores de {n(gen['genes_en_cromosomas'])} genes en cromosomas "
    f"(genes con CDS completa en 5' en los {J(ISO, 'Egl', 'cromosomas')} cromosomas); truncados {gen['truncados']}; con N {gen['con_N']}",
    ok(gen["truncados"] == 0 and gen["con_N"] == 0),
    [fu(PROM, k) for k in ("genes_en_cromosomas", "promotores", "truncados", "con_N")])
prt = pd.read_csv(os.path.join(E3, "promotores", "promotores_2kb_todos.tsv"), sep="\t")
cod = prt.loc[~es_v(prt.atg_verificado), "codon_siguiente"].value_counts()
assert len(prt) == gen["promotores"] and es_v(prt.atg_verificado).sum() == gen["ATG_verificado"]
ctl("04", "Genoma: codón inmediatamente después de la ventana", "",
    f"ATG en {n(gen['ATG_verificado'])} de {n(gen['promotores'])}; otros: "
    + ", ".join(f"{c} {k}" for c, k in cod.items()) + " (ninguno en los NAC)",
    "INFORMATIVO", [fu(PROM, "ATG_verificado"), "datos/entrega3/promotores/promotores_2kb_todos.tsv"])
sol = {k: J(PROM, "NAC", k) for k in ("solapan_otro_gen", "solape_mediana_bp")}
fon = {k: J(PROM, "fondo", k) for k in ("n", "solapan_otro_gen", "GC_mediana")}
s40 = t40.solape_otros_genes_bp
assert (s40 > 0).sum() == sol["solapan_otro_gen"] and s40.median() == sol["solape_mediana_bp"]
ctl("04", "Ventanas que caen dentro de otros genes", "",
    f"NAC {sol['solapan_otro_gen']}/189 (solape en esos {sol['solapan_otro_gen']}: mediana {d(s40[s40 > 0].median(), 1)} pb, "
    f"{s40[s40 > 0].min()}-{s40[s40 > 0].max()} pb; en los 189: mediana {d(sol['solape_mediana_bp'], 0)} pb); "
    f"fondo {n(fon['solapan_otro_gen'])}/{n(fon['n'])} ({pc(fon['solapan_otro_gen'], fon['n'])})",
    "INFORMATIVO", [fu(PROM, "NAC", k) for k in sol] + [fu(PROM, "fondo", "solapan_otro_gen"), fu(PROM, "fondo", "n"), f40])
gcn, gcr = J(PROM, "NAC", "GC_mediana"), J(PROM, "NAC", "GC_rango")
ctl("04", "Contenido en GC de los promotores", "",
    f"mediana NAC {d(gcn, 3)} (rango {d(gcr[0], 3)}-{d(gcr[1], 3)}); mediana del fondo {d(fon['GC_mediana'], 4)}: "
    "motiva el fondo emparejado por GC",
    "INFORMATIVO", [fu(PROM, "NAC", "GC_mediana"), fu(PROM, "NAC", "GC_rango"), fu(PROM, "fondo", "GC_mediana")])

# --- 05 envío a PlantCARE ------------------------------------------------------------------------------
ENV = os.path.join(E3, "cis_plantcare", "envios")
reg = pd.read_csv(os.path.join(ENV, "registro_envios.tsv"), sep="\t")
trab = pd.read_csv(os.path.join(ENV, "trabajos_plantcare.tsv"), sep="\t")
lot = J(PROM, "lotes")
md5_lotes = [md5_fichero(os.path.join(ENV, f)) for f in reg.fichero]
ctl("05", "Cuatro lotes de menos de 100 000 bytes con los 189 promotores; md5 registrado",
    "4 lotes < 100 000 bytes; 189 secuencias en total; md5 actual del lote = md5 registrado al enviar",
    f"bytes {', '.join(n(b) for b in reg.bytes)}; secuencias {'+'.join(str(s) for s in reg.secuencias)} = "
    f"{reg.secuencias.sum()}; md5 iguales {sum(a == b for a, b in zip(md5_lotes, reg.md5))}/4",
    ok(len(reg) == 4 and (reg.bytes < 100000).all() and reg.secuencias.sum() == 189
       and all(a == b for a, b in zip(md5_lotes, reg.md5))
       and [l_["bytes"] for l_ in lot] == reg.bytes.tolist() and [l_["n"] for l_ in lot] == reg.secuencias.tolist()),
    ["datos/entrega3/cis_plantcare/envios/registro_envios.tsv", fu(PROM, "lotes")])
horas = pd.to_datetime(reg.hora_envio)
gaps = horas.diff().dt.total_seconds().dropna().astype(int).tolist()
ctl("05", "Espera entre envíos", ">= 60 s entre envíos consecutivos",
    f"{', '.join(str(g) for g in gaps)} s (envíos del {horas.iloc[0]:%d-%m-%Y}, {horas.iloc[0]:%H:%M:%S}-{horas.iloc[-1]:%H:%M:%S})",
    ok(all(g >= 60 for g in gaps)), ["datos/entrega3/cis_plantcare/envios/registro_envios.tsv"])
ctl("05", "Respuesta del servidor de PlantCARE",
    "HTTP 200 y secuencias aceptadas = enviadas en los 4 trabajos",
    f"HTTP {', '.join(str(s) for s in reg.estado_http)}; aceptadas {'+'.join(str(s) for s in trab.secuencias_aceptadas)}; "
    f"trabajos {', '.join(trab.trabajo)}",
    ok((reg.estado_http == 200).all() and (trab.secuencias_aceptadas == trab.secuencias_enviadas).all()
       and trab.md5.tolist() == reg.md5.tolist()),
    ["datos/entrega3/cis_plantcare/envios/trabajos_plantcare.tsv"])
tr_res = [J(LEC, "lotes", str(i), "trabajo_resultado") for i in range(1, 5)]
tr_reg = [J(LEC, "lotes", str(i), "trabajo_registro") for i in range(1, 5)]
ctl("05", "Resultados recibidos", "",
    f"los trabajos del envío automático ({', '.join(tr_reg)}) no llegaron por correo (dominio del correo "
    f"registrado mal escrito); los mismos 4 lotes se subieron desde la web ({', '.join(tr_res)}); la identidad se "
    "comprueba por contenido en el paso «Leer y verificar»",
    "INFORMATIVO", [fu(LEC, "lotes", str(i), k) for i in range(1, 5) for k in ("trabajo_resultado", "trabajo_registro")]
    + [fu(LEC, "entrada", "nota")])
J(LEC, "entrada", "nota")

# --- 06 homólogos de Arabidopsis ---------------------------------------------------------------------
h = {k: J(HOM, k) for k in ("genes_Egl", "genes_Ath_nucleares", "genes_Ath_organulos_excluidos",
                            "consultas_con_acierto", "genes_Egl_con_homologo", "AGI_distintos_universo",
                            "NAC_Ath_UniProt", "NAC_Ath_con_numero_ANAC", "NAC189_con_homologo",
                            "NAC189_mejor_es_NAC", "NAC189_AGI_distintos")}
fil = J(HOM, "filtros")
t41, f41 = tabla("T41")
ctl("06", "Proporción de los 189 cuyo mejor homólogo es un NAC de Arabidopsis",
    "todo NAC con homólogo que pasa los filtros tiene por mejor acierto un NAC de Arabidopsis (lista de UniProt)",
    f"{h['NAC189_mejor_es_NAC']}/{h['NAC189_con_homologo']} ({pc(h['NAC189_mejor_es_NAC'], h['NAC189_con_homologo'])}); "
    f"con homólogo {h['NAC189_con_homologo']}/189 ({pc(h['NAC189_con_homologo'], 189)}); {h['NAC189_AGI_distintos']} AGI "
    f"distintos; T41: {es_v(t41.pasa_filtros).sum()} pasan filtros, {es_v(t41.es_NAC_Ath).sum()} NAC",
    ok(h["NAC189_mejor_es_NAC"] == h["NAC189_con_homologo"] and len(t41) == 189
       and es_v(t41.pasa_filtros).sum() == h["NAC189_con_homologo"]
       and (es_v(t41.pasa_filtros) & es_v(t41.es_NAC_Ath)).sum() == h["NAC189_mejor_es_NAC"]),
    [fu(HOM, "NAC189_mejor_es_NAC"), fu(HOM, "NAC189_con_homologo"), fu(HOM, "NAC189_AGI_distintos"), f41])
ctl("06", "Filtros y universo de homólogos", "",
    f"e <= {ec(fil['evalue'], 0)}, identidad >= {d(fil['identidad_min'], 0)} %, cobertura de la consulta >= "
    f"{d(fil['cobertura_consulta_min'], 0)} %; {n(h['genes_Egl_con_homologo'])}/{n(h['genes_Egl'])} genes de "
    f"E. globulus con homólogo; {n(h['AGI_distintos_universo'])} AGI distintos; {n(h['genes_Ath_nucleares'])} genes "
    f"nucleares de Arabidopsis ({h['genes_Ath_organulos_excluidos']} de orgánulos excluidos)",
    "INFORMATIVO", [fu(HOM, "filtros")] + [fu(HOM, k) for k in ("genes_Egl_con_homologo", "genes_Egl",
                                                                  "AGI_distintos_universo", "genes_Ath_nucleares",
                                                                  "genes_Ath_organulos_excluidos", "consultas_con_acierto")])
u = {k: J(UNI, k) for k in ("entradas", "agi_unicos", "fecha_descarga")}
ctl("06", "Lista de NAC de Arabidopsis (UniProt)", "",
    f"UniProt {J(UNI, 'cabeceras', 'X-UniProt-Release')}: {u['entradas']} entradas, {u['agi_unicos']} AGI "
    f"(PF02365, IPR003441 o PS51005 en UP000006548); {h['NAC_Ath_con_numero_ANAC']} con número ANAC; descarga {u['fecha_descarga']}",
    "INFORMATIVO", [fu(UNI, "cabeceras", "X-UniProt-Release"), fu(UNI, "entradas"), fu(UNI, "agi_unicos"),
                    fu(HOM, "NAC_Ath_con_numero_ANAC")])
J(UNI, "consulta")

# --- 07 enriquecimiento GO ----------------------------------------------------------------------------
TOL = 1e-6   # criterio de e3_06_go_recalculo.py
for ont in ("BP", "MF", "CC"):
    for un, et in (("AGI", "unidad = AGI"), ("Egl", "unidad = gen de E. globulus")):
        r = {k: J(REC, ont, un, k) for k in ("terminos_py", "terminos_R", "mismos_terminos", "k_igual",
                                             "dif_rel_max_p", "dif_rel_max_padj", "pasa")}
        ctl("07", f"Recálculo hipergeométrico independiente en Python = clusterProfiler ({ont}, {et})",
            "mismos términos, mismo k y diferencia relativa < 1e-6 en p y en p.adjust",
            f"términos {r['terminos_py']} (Python) y {r['terminos_R']} (R); k igual: {'sí' if r['k_igual'] else 'no'}; "
            f"dif. rel. máx. p {ec(r['dif_rel_max_p'])}, p.adjust {ec(r['dif_rel_max_padj'])}",
            ok(r["mismos_terminos"] and r["k_igual"] and r["pasa"] and r["dif_rel_max_p"] < TOL
               and r["dif_rel_max_padj"] < TOL),
            [fu(REC, ont, un, k) for k in r] + ([fu(REC, "pasa_todo")] if (ont, un) == ("BP", "AGI") else []))
J(REC, "pasa_todo")
t43, f43 = tabla("T43")
gsig = {o: J(GO, o, "significativos") for o in ("BP", "MF", "CC")}
gegl = {o: J(GO, o, "significativos_contraste_Egl") for o in ("BP", "MF", "CC")}
gpro = {o: J(GO, o, "terminos_probados") for o in ("BP", "MF", "CC")}
t43c = {o: (es_v(g.significativo).sum(), es_v(g.significativo_Egl).sum(), (~es_v(g.supera_maxGSSize)).sum())
        for o, g in t43.groupby("ontologia")}
ctl("07", "Términos probados y significativos: resumen = T43",
    "en cada ontología, significativos (unidad AGI y unidad gen) y términos probados de T43 = resumen",
    "; ".join(f"{o}: probados {gpro[o]}, significativos {gsig[o]} (AGI) / {gegl[o]} (genes)" for o in gsig)
    + f"; BP tras simplify {J(GO, 'BP', 'tras_simplify')}",
    ok(all(t43c[o] == (gsig[o], gegl[o], gpro[o]) for o in gsig)
       and es_v(t43.representativo).sum() == J(GO, "BP", "tras_simplify")),
    [fu(GO, o, k) for o in gsig for k in ("significativos", "significativos_contraste_Egl", "terminos_probados")]
    + [fu(GO, "BP", "tras_simplify"), f43])
pos = []
for o in ("BP", "MF", "CC"):
    p0 = J(GOF, f"sensibilidad_{o}", "primeros_mayores_500", 0)
    pos.append((o, p0["ID"], p0["termino"], p0["GeneRatio"], p0["p.adjust"]))
ctl("07", "Control positivo de la transferencia de anotación (sin límite superior de tamaño)",
    "con maxGSSize = Inf, los términos propios de un factor de transcripción salen significativos",
    "; ".join(f"{o} {i} «{t}» {g}, p.adjust {ec(p)}" for o, i, t, g, p in pos),
    ok(all(p < 0.05 for *_, p in pos)
       and {i for _, i, *_ in pos} == {"GO:0006355", "GO:0003700", "GO:0005634"}),
    [fu(GOF, f"sensibilidad_{o}", "primeros_mayores_500", 0) for o in ("BP", "MF", "CC")])
daga = int((es_v(t43.representativo) & ~es_v(t43.significativo_Egl)).sum())
ctl("07", "Contraste con los genes de E. globulus como unidades (BP)", "",
    f"significativos {gsig['BP']} (AGI) frente a {gegl['BP']} (genes); compartidos {J(GO, 'BP', 'compartidos')}; "
    f"{gsig['BP'] - J(GO, 'BP', 'compartidos')} de los {gsig['BP']} no lo son con genes como unidades; "
    f"{daga} de los {J(GO, 'BP', 'tras_simplify')} representativos de F13 {'lleva' if daga == 1 else 'llevan'} daga por eso",
    "INFORMATIVO", [fu(GO, "BP", "compartidos"), f43])
PROPIO["go"] = dict(representativos_con_daga=daga)
ctl("07", "Universo anotado", "",
    "; ".join(f"{o}: N = {n(J(GO, o, 'N_universo_anotado'))}" for o in ("BP", "MF", "CC"))
    + f"; universo AGI {n(J(GO, 'parametros', 'universo_AGI'))}, universo de genes {n(J(GO, 'parametros', 'universo_Egl'))}",
    "INFORMATIVO", [fu(GO, o, "N_universo_anotado") for o in ("BP", "MF", "CC")]
    + [fu(GO, "parametros", "universo_AGI"), fu(GO, "parametros", "universo_Egl")])


# --- figuras (pasos 08, 12, 15) ----------------------------------------------------------------------
def pdf_info(p):
    b = open(p, "rb").read()
    textos = [b]
    for m in re.finditer(rb"stream\r?\n", b):
        fin = b.find(b"endstream", m.end())
        if fin > 0:
            try:
                textos.append(zlib.decompressobj().decompress(b[m.end():fin]))
            except zlib.error:
                pass
    fuentes, cajas = set(), set()
    for t in textos:
        fuentes |= {f.decode("latin-1").split("+")[-1] for f in re.findall(rb"/BaseFont\s*/([^\s/\[\]<>()]+)", t)}
        cajas |= {tuple(round(float(v), 2) for v in c) for c in
                  re.findall(rb"/MediaBox\s*\[\s*([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s*\]", t)}
    return fuentes, cajas


PROPIO["figuras"] = {}


def control_figura(paso, base):
    png, pdf = os.path.join(FIG, base + ".png"), os.path.join(FIG, base + ".pdf")
    im = Image.open(png)
    dpi = im.info.get("dpi", (0, 0))
    mm = [s / r * 25.4 for s, r in zip(im.size, dpi)]
    fuentes, cajas = pdf_info(pdf)
    caja = list(cajas)[0] if len(cajas) == 1 else None
    mm_pdf = [(caja[2] - caja[0]) / 72 * 25.4, (caja[3] - caja[1]) / 72 * 25.4] if caja else [float("nan")] * 2
    c = (abs(dpi[0] - 600) < 0.5 and abs(mm[0] - 174) <= 0.3 and mm[1] <= 234 and caja is not None
         and all(abs(a - b) <= 0.5 for a, b in zip(mm, mm_pdf))
         and fuentes and all(f.startswith("LiberationSans") for f in fuentes))
    ctl(paso, f"Figura {base}: tamaño, resolución y fuentes",
        "PNG a 600 ppp, 174 mm de ancho y <= 234 mm de alto; página del PDF del mismo tamaño (+/- 0,5 mm); "
        "solo Liberation Sans incrustada",
        f"PNG {im.size[0]} x {im.size[1]} px a {d(dpi[0], 0)} ppp = {d(mm[0], 1)} x {d(mm[1], 1)} mm; PDF "
        f"{d(mm_pdf[0], 1)} x {d(mm_pdf[1], 1)} mm; fuentes: {', '.join(sorted(fuentes))}",
        ok(c), [f"figuras/{base}.png", f"figuras/{base}.pdf"])
    PROPIO["figuras"][base] = dict(px=list(im.size), mm=[round(x, 1) for x in mm], fuentes=sorted(fuentes))


def control_log(paso, nombre, marca, que):
    t = open(os.path.join(LOG, nombre), encoding="utf-8", errors="replace").read()
    err = [l for l in t.splitlines() if re.match(r"\s*(Error|Traceback|Execution halted)", l)]
    hallado = marca in t[-6000:]
    ctl(paso, f"Comprobaciones internas de {que} (el guion para si alguna falla)",
        f"log completo (aparece «{marca}» al final) y sin errores",
        f"{'completo' if hallado else 'incompleto'}; líneas de error: {len(err)}",
        ok(hallado and not err), [f"logs/{nombre}"])


for b_ in ("F13_GO_dotplot", "F13S_GO_arbol", "F14_GO_subfamilias"):
    control_figura("08", b_)
control_log("08", "entrega3_09_go_figuras.log", "cifras -> _cifras/e3_go.csv", "e3_08_go_figuras.R (recuentos, representativos, "
            "hojas del árbol, listas de compareCluster, EglNAC por término)")
t42, f42 = tabla("T42")
t44, f44 = tabla("T44")
cc_prob = sum(J(GOF, "compareCluster_subfamilia", "probados").values()) + sum(J(GOF, "compareCluster_grupo", "probados").values())
ctl("08", "Tablas de GO: filas",
    "T42 = 189 genes; T44 = términos probados por subfamilia y por grupo en compareCluster",
    f"T42 {len(t42)}; T44 {n(len(t44))} = {n(cc_prob)} probados", ok(len(t42) == 189 and len(t44) == cc_prob),
    [f42, f44, fu(GOF, "compareCluster_subfamilia", "probados"), fu(GOF, "compareCluster_grupo", "probados")])

# --- 09 colinealidad intraespecífica -----------------------------------------------------------------
mc = {k: J(COL, "genoma", "MCScanX", k) for k in ("bloques", "pares_ancla", "genes_ancla", "clase4_igual_a_anclas")}
cls_g = J(COL, "genoma", "clases")
ctl("09", "Clase WGD/segmental de duplicate_gene_classifier = genes ancla de los bloques de MCScanX",
    "mismos genes", f"{n(cls_g['WGD/segmental'])} genes en la clase 4 y {n(mc['genes_ancla'])} genes ancla "
                    f"({mc['bloques']} bloques, {n(mc['pares_ancla'])} pares ancla)",
    ok(mc["clase4_igual_a_anclas"] and cls_g["WGD/segmental"] == mc["genes_ancla"]),
    [fu(COL, "genoma", "MCScanX", k) for k in mc] + [fu(COL, "genoma", "clases", "WGD/segmental")])
cls_n = J(COL, "NAC189", "clases")
ctl("09", "Las clases suman todos los genes",
    "genoma: suma de las 5 clases = genes en cromosomas; NAC: suma = 189",
    f"genoma: {', '.join(f'{k} {n(v)}' for k, v in cls_g.items())} = {n(sum(cls_g.values()))} de {n(J(COL, 'genoma', 'genes'))}; "
    f"NAC: {', '.join(f'{k} {v}' for k, v in cls_n.items())} = {sum(cls_n.values())}",
    ok(sum(cls_g.values()) == J(COL, "genoma", "genes") == J(ISO, "Egl", "en_cromosomas") and sum(cls_n.values()) == 189),
    [fu(COL, "genoma", "clases"), fu(COL, "NAC189", "clases"), fu(COL, "genoma", "genes")])
ct = lambda *r: J(COL, "NAC189", "contraste_tandem", *r)
ctl("09", "Pares en tándem de MCScanX = pares NAC adyacentes con acierto DIAMOND",
    "mismos pares", f"{ct('MCScanX', 'pares_tandem')} pares; coinciden: "
                    f"{'sí' if ct('MCScanX', 'pares_tandem_igual_a_adyacentes_con_acierto') else 'no'}",
    ok(ct("MCScanX", "pares_tandem_igual_a_adyacentes_con_acierto")),
    [fu(COL, "NAC189", "contraste_tandem", "MCScanX", k) for k in ("pares_tandem", "pares_tandem_igual_a_adyacentes_con_acierto")])
ctl("09", "Regla de 200 kb (artículo 3) frente a T01.bloque_tandem",
    "todos los genes de los bloques de T01 (<= 5 genes intercalados y <= 200 kb) están en grupos de 200 kb",
    f"grupos de 200 kb: {ct('art3', 'grupos')} con {ct('art3', 'genes')} genes ({ct('art3', 'pares')} pares); T01: "
    f"{ct('T01', 'bloques')} bloques con {ct('T01', 'genes')} genes, {ct('T01', 'genes_en_art3')} en grupos de 200 kb, "
    f"{ct('T01', 'T01_fuera_de_art3')} fuera; {ct('T01', 'art3_fuera_de_T01')} genes de los grupos fuera de T01; "
    f"{ct('T01', 'bloques_identicos_a_un_grupo_art3')} bloques idénticos a un grupo",
    ok(ct("T01", "T01_fuera_de_art3") == 0 and ct("T01", "genes_en_art3") == ct("T01", "genes")),
    [fu(COL, "NAC189", "contraste_tandem", "art3", k) for k in ("grupos", "genes", "pares")]
    + [fu(COL, "NAC189", "contraste_tandem", "T01", k) for k in ("bloques", "genes", "genes_en_art3", "T01_fuera_de_art3",
                                                                  "art3_fuera_de_T01", "bloques_identicos_a_un_grupo_art3")])
ctl("09", "Bloques en tándem de T15 dentro de los grupos de 200 kb y de T01",
    "los genes de T15 están en ambos y sus símbolos = los listados en T15",
    f"T15: {ct('T15', 'bloques')} bloques, {ct('T15', 'genes')} genes; en grupos de 200 kb {ct('T15', 'genes_en_art3')}; "
    f"en T01 {ct('T15', 'genes_en_T01')}; símbolos iguales: {'sí' if ct('T15', 'coincide_con_simbolos_listados') else 'no'}",
    ok(ct("T15", "coincide_con_simbolos_listados") and ct("T15", "genes_en_art3") == ct("T15", "genes") == ct("T15", "genes_en_T01")),
    [fu(COL, "NAC189", "contraste_tandem", "T15", k) for k in ("bloques", "genes", "genes_en_art3", "genes_en_T01",
                                                                "coincide_con_simbolos_listados")])
t45, f45 = tabla("T45")
t46, f46 = tabla("T46")
ctl("09", "Tablas de duplicación: filas", "T45 = 189 genes; T46 = pares NAC del resumen",
    f"T45 {len(t45)}; T46 {len(t46)} = {J(COL, 'NAC189', 'pares', 'T46')}",
    ok(len(t45) == 189 and len(t46) == J(COL, "NAC189", "pares", "T46")), [f45, f46, fu(COL, "NAC189", "pares", "T46")])

# --- 10 Ka/Ks y datación --------------------------------------------------------------------------------
cc_ = J(KK, "conciliacion_CDS_proteina")
pn = {k: J(KK, "pares_NAC", k) for k in ("analizados", "Ks_saturado")}
ok_pares = J(KK, "pares_NAC", "estados", "ok")
an = {k: J(KK, "anclas_genoma", k) for k in ("pares", "estimados")}
ctl("10", "Conciliación CDS-proteína antes de alinear",
    "todos los pares NAC con CDS conciliada (traducción = proteína) y estimación válida",
    f"genes: {n(cc_['exactos'])} exactos de {n(cc_['genes'])}, {cc_['con_diferencias']} con diferencias, "
    f"{cc_['sin_conciliar']} sin conciliar (fuera: {', '.join(cc_['genes_sin_conciliar'])}); pares NAC "
    f"{ok_pares}/{pn['analizados']}; anclas estimadas {n(an['estimados'])}/{n(an['pares'])}",
    ok(ok_pares == pn["analizados"] == J(COL, "NAC189", "pares", "T46") and cc_["sin_cds"] == 0),
    [fu(KK, "conciliacion_CDS_proteina"), fu(KK, "pares_NAC", "estados", "ok"), fu(KK, "pares_NAC", "analizados"),
     fu(KK, "anclas_genoma", "pares"), fu(KK, "anclas_genoma", "estimados")])
rc = J(KK, "control_recalculo_5_pares")
ctl("10", "Recálculo independiente de 5 pares (Ka y Ks con yn00)", "diferencia absoluta máxima = 0",
    f"pares {', '.join(rc['claves'])}; diferencia máxima {d(rc['dif_abs_max_Ka_Ks_YN00'], 1)}",
    ok(rc["dif_abs_max_Ka_Ks_YN00"] == 0), [fu(KK, "control_recalculo_5_pares")])
t47, f47 = tabla("T47")
LAMBDA = J(KK, "parametros", "lambda_sitio_anio")
KSAT = J(KK, "parametros", "Ks_saturado")
sat = es_v(t47.Ks_saturado)
dat = ~sat & t47.T_Mya_YN00.notna()
dmax = (t47.loc[dat, "T_Mya_YN00"] - t47.loc[dat, "Ks_YN00"] / (2 * LAMBDA) / 1e6).abs().max()
ctl("10", "Ks > 2 marcados y sin datar; T = Ks / 2λ en el resto (T47)",
    f"saturado si y solo si Ks > {d(KSAT, 0)}; sin T en los saturados; |T - Ks/2λ| < 0,01 Ma en los demás",
    f"{sat.sum()} saturados de {len(t47)} (resumen: {pn['Ks_saturado']}); con Ks > 2: {(t47.Ks_YN00 > KSAT).sum()}; "
    f"saturados con T: {t47.loc[sat, 'T_Mya_YN00'].notna().sum()}; no saturados sin T: {(~sat & t47.T_Mya_YN00.isna()).sum()}; "
    f"|T - Ks/2λ| máx. {d(dmax, 4)} Ma (λ = {ec(LAMBDA)})",
    ok(len(t47) == pn["analizados"] and sat.sum() == pn["Ks_saturado"] == (t47.Ks_YN00 > KSAT).sum()
       and t47.loc[sat, "T_Mya_YN00"].isna().all() and t47.loc[~sat, "T_Mya_YN00"].notna().all() and dmax < 0.01),
    [f47, fu(KK, "pares_NAC", "Ks_saturado"), fu(KK, "parametros", "lambda_sitio_anio"), fu(KK, "parametros", "Ks_saturado")])
pk = J(KK, "anclas_genoma", "picos_KDE_anclas_YN00_PAML")
ks_pico = pk["picos"][0]["Ks"]
t_pico = J(KK, "anclas_genoma", "T_Mya_del_pico_principal")
ks_ng = J(KK, "anclas_genoma", "picos_KDE_anclas_NG86")[0]["Ks"]
ks_bl = J(KK, "anclas_genoma", "picos_KDE_medianas_por_bloque", "picos")[0]["Ks"]
J(KK, "anclas_genoma", "picos_KDE_medianas_por_bloque", "bloques")
MYB = dict(moda=109.93, ic=(105.96, 113.91))   # Myburg et al. 2014, Nature 510:356, Extended Data Fig. 4
ctl("10", "Pico de la duplicación antigua de Eucalyptus en los Ks de las anclas del genoma",
    "pico principal de la densidad de Ks (yn00) con T = Ks/2λ dentro del IC 95 % de la WGD de Eucalyptus "
    f"(Myburg et al. 2014: {d(MYB['moda'])} Ma, {d(MYB['ic'][0])}-{d(MYB['ic'][1])} Ma)",
    f"{len(pk['picos'])} pico en {n(pk['n'])} anclas no saturadas: Ks {d(ks_pico, 3)} (mitad de altura "
    f"{d(pk['picos'][0]['mitad_altura'][0], 3)}-{d(pk['picos'][0]['mitad_altura'][1], 3)}) -> T {d(t_pico, 1)} Ma; "
    f"NG86 {d(ks_ng, 3)} ({d(ks_ng / (2 * LAMBDA) / 1e6, 1)} Ma); medianas por bloque {d(ks_bl, 3)} "
    f"({d(ks_bl / (2 * LAMBDA) / 1e6, 1)} Ma); Myburg dató por filogenómica, no por Ks/2λ",
    ok(len(pk["picos"]) >= 1 and abs(t_pico - ks_pico / (2 * LAMBDA) / 1e6) < 0.1 and MYB["ic"][0] <= t_pico <= MYB["ic"][1]),
    [fu(KK, "anclas_genoma", "picos_KDE_anclas_YN00_PAML"), fu(KK, "anclas_genoma", "T_Mya_del_pico_principal"),
     fu(KK, "anclas_genoma", "picos_KDE_anclas_NG86"), fu(KK, "anclas_genoma", "picos_KDE_medianas_por_bloque")])
PROPIO["datacion"] = dict(
    Ks_pico=ks_pico, T_pico_Ma=t_pico, lambda_sitio_anio=LAMBDA,
    Myburg2014=dict(moda_Ma=MYB["moda"], IC95_Ma=list(MYB["ic"]), metodo="datación filogenómica de anclas",
                    fuente="Myburg et al. 2014, Nature 510:356-362, p. 2 y leyenda de Extended Data Fig. 4"),
    dentro_IC=MYB["ic"][0] <= t_pico <= MYB["ic"][1])
cn = J(KK, "concordancia")
ctl("10", "Concordancia entre estimadores de Ks", "",
    "; ".join(f"{g}: yn00 de PAML frente a YN00 de Biopython rho {d(v['rho_Ks_YN00_PAML_vs_Biopython'], 3)} "
              f"(Biopython no estimable en {n(v['YN00_Biopython_no_estimable'])} de {n(v['de'])}), frente a NG86 "
              f"rho {d(v['rho_Ks_YN00_PAML_vs_NG86'], 3)} (mediana de la diferencia {d(v['mediana_dif_rel_Ks_PAML_vs_NG86_pct'], 1)} %)"
              for g, v in cn.items()),
    "INFORMATIVO", [fu(KK, "concordancia")])

# --- 11 colinealidad interespecífica -------------------------------------------------------------------
hmm = {l.split()[0] for l in open(os.path.join(E3, "colinealidad", "jcvi", "hmm", "Egl_PF02365_cutga.domtbl"))
       if not l.startswith("#")}
ctl("11", "hmmsearch --cut_ga con PF02365 en E. globulus = los 189 loci de T01",
    "mismo conjunto de loci", f"{len(hmm)} loci con dominio; diferencia simétrica con T01: {len(hmm ^ set(t01.locus))}; "
                              f"resumen: {J(INT, 'NAC_por_especie', 'Egl')}",
    ok(hmm == set(t01.locus) and J(INT, "NAC_por_especie", "Egl") == 189),
    ["datos/entrega3/colinealidad/jcvi/hmm/Egl_PF02365_cutga.domtbl", fu(INT, "NAC_por_especie", "Egl")])
for sp in ("Egr", "Ptr"):
    ej = JS[f"colinealidad/jcvi/Egl.{sp}.ejecucion.json"]
    for k in ej:
        J(f"colinealidad/jcvi/Egl.{sp}.ejecucion.json", k)
    x = {k: J(INT, "NAC", sp, k) for k in ("NAC_hmm", "NAC_otro_con_par_en_Egl", "NAC_otro_cuyo_mejor_socio_es_NAC",
                                           "NAC_otro_cuyo_mejor_socio_no_es_NAC", "NAC_Egl_con_ortologo",
                                           "NAC_Egl_ortologo_NAC", "NAC_Egl_ortologo_no_NAC", "ortologo_reciproco")}
    ctl("11", f"JCVI E. globulus - {NOMBRE[sp]}: ejecución y coherencia de los NAC",
        "código de salida 0; todo NAC de la otra especie con socio colineal en E. globulus tiene por mejor socio un NAC",
        f"código {ej['codigo']}; NAC (hmmsearch) {x['NAC_hmm']}; con socio en E. globulus {x['NAC_otro_con_par_en_Egl']}, "
        f"mejor socio NAC {x['NAC_otro_cuyo_mejor_socio_es_NAC']}; NAC de E. globulus con ortólogo {x['NAC_Egl_con_ortologo']}/189 "
        f"(NAC {x['NAC_Egl_ortologo_NAC']}; sin PF02365 sobre el umbral GA: {len(x['NAC_Egl_ortologo_no_NAC'])}; "
        f"recíprocos {x['ortologo_reciproco']})",
        ok(ej["codigo"] == 0 and not x["NAC_otro_cuyo_mejor_socio_no_es_NAC"]
           and x["NAC_otro_cuyo_mejor_socio_es_NAC"] == x["NAC_otro_con_par_en_Egl"]),
        [fu(f"colinealidad/jcvi/Egl.{sp}.ejecucion.json", "codigo")] + [fu(INT, "NAC", sp, k) for k in x])
cat = J(INT, "contraste_T01", "categorias")
sinco = J(INT, "contraste_T01", "sin_coortologo_T01")
ctl("11", "Contraste (a): ortólogo colineal en E. grandis frente al co-ortólogo EgrNAC de T01", "",
    "; ".join(f"{k} {v}" for k, v in cat.items()) + f" (total {sum(cat.values())}); sin co-ortólogo en T01: "
    + ", ".join(f"{x['simbolo']} -> {x['simbolo_EgrNAC_ortologo'] or x['ortologo_Egr'] or 'sin ortólogo colineal'}"
                for x in sinco),
    "INFORMATIVO", [fu(INT, "contraste_T01", "categorias"), fu(INT, "contraste_T01", "sin_coortologo_T01")])
mp = {k: J(INT, "mapeo_posicional_anterior", k) for k in ("filas", "en_cromosomas_Egr_v2", "con_ortologo_colineal",
                                                           "cromosoma_antiguo_coincide", "sintenia_coincide_con_secuencia",
                                                           "secuencia_entre_socios_colineales")}
t50, f50 = tabla("T50")
disc = J(INT, "mapeo_posicional_anterior", "discrepancias")
ctl("11", "Contraste (b): mapeo posicional anterior", "",
    f"{mp['filas']} símbolos (T50 {len(t50)}); en cromosomas de la v2 {mp['en_cromosomas_Egr_v2']}; con ortólogo colineal "
    f"{mp['con_ortologo_colineal']}; cromosoma anterior coincide {mp['cromosoma_antiguo_coincide']}; sintenia = secuencia "
    f"{mp['sintenia_coincide_con_secuencia']}; en {len(disc)} discrepancias el locus por secuencia está entre los socios "
    f"colineales en {sum(x['secuencia_entre_socios_colineales'] is True for x in disc)} "
    f"(total: {mp['secuencia_entre_socios_colineales']})",
    "INFORMATIVO", [fu(INT, "mapeo_posicional_anterior", k) for k in mp]
    + [fu(INT, "mapeo_posicional_anterior", "discrepancias"), f50])
cob = {sp: (J(INT, "genoma", sp, "pct_genes_Egl"), J(INT, "genoma", sp, "pct_genes_otro"),
            J(INT, "genoma", sp, "bloques"), J(INT, "genoma", sp, "pares_ancla")) for sp in ("Egr", "Ptr")}
ctl("11", "Cobertura colineal del genoma", "",
    "; ".join(f"{NOMBRE[sp]}: {c[2]} bloques, {n(c[3])} pares ancla; genes de E. globulus en bloques {d(c[0], 1)} %, "
              f"de la otra especie {d(c[1], 1)} %" for sp, c in cob.items()),
    "INFORMATIVO", [fu(INT, "genoma", sp, k) for sp in cob for k in ("pct_genes_Egl", "pct_genes_otro", "bloques", "pares_ancla")])
t48, f48 = tabla("T48")
t49, f49 = tabla("T49")
ctl("11", "Tablas de ortólogos: filas", "T48 = 189 genes; T50 = símbolos del mapeo anterior",
    f"T48 {len(t48)}; T49 {len(t49)} pares; T50 {len(t50)}", ok(len(t48) == 189 and len(t50) == mp["filas"]),
    [f48, f49, f50])

# --- 12 figuras de colinealidad --------------------------------------------------------------------------
for b_ in ("F15_circos_colinealidad_Egl", "F16_cariotipo_Egr_Egl_Ptr", "F17_KaKs_datacion"):
    control_figura("12", b_)
control_log("12", "entrega3_16_figuras_colinealidad.log", "exit=0", "e3_12_figuras_colinealidad.py (tamaño del PNG, leyendas)")

# --- 13 lectura de PlantCARE -----------------------------------------------------------------------------
lotes = {i: {k: J(LEC, "lotes", str(i), k) for k in ("md5_igual_al_registro", "secuencias", "html", "filas_tab", "md5_lote")}
         for i in range(1, 5)}
ctl("13", "Los resultados corresponden a los lotes registrados",
    "md5 del lote = md5 registrado al enviar, en los 4 lotes; un HTML por gen",
    "; ".join(f"lote {i}: md5 {'igual' if v['md5_igual_al_registro'] else 'distinto'}, {v['secuencias']} secuencias, "
              f"{v['html']} HTML, {n(v['filas_tab'])} filas .tab" for i, v in lotes.items()),
    ok(all(v["md5_igual_al_registro"] and v["secuencias"] == v["html"] for v in lotes.values())
       and [v["md5_lote"] for v in lotes.values()] == reg.md5.tolist()),
    [fu(LEC, "lotes", str(i), k) for i in lotes for k in lotes[i]])
en = {k: J(LEC, "entrada", k) for k in ("genes", "secuencias_lotes_iguales_a_NAC189_fa", "largo", "solo_ACGT")}
ctl("13", "Secuencias de los lotes = FASTA de los 189 promotores", "189 genes; secuencias idénticas; 2000 pb; solo ACGT",
    f"genes {en['genes']}; idénticas: {'sí' if en['secuencias_lotes_iguales_a_NAC189_fa'] else 'no'}; largo {en['largo']}; "
    f"solo ACGT: {'sí' if en['solo_ACGT'] else 'no'}",
    ok(en["genes"] == 189 and en["secuencias_lotes_iguales_a_NAC189_fa"] and en["largo"] == 2000 and en["solo_ACGT"]),
    [fu(LEC, "entrada", k) for k in en])
si = {k: J(LEC, "sitios", k) for k in ("filas_con_nombre", "dentro_de_la_secuencia", "puntaje_calculado_igual_al_declarado",
                                       "puntaje_>=0.9_del_maximo", "coincidencia_exacta", "fraccion_exacta",
                                       "con_un_fallo_o_degenerada_incompleta")}
ctl("13", "Cada sitio en la posición y la hebra declaradas, con el puntaje declarado",
    "100 % de las filas con nombre: sitio dentro de la secuencia enviada, puntaje calculado = declarado y >= 0,9 del máximo",
    f"de {n(si['filas_con_nombre'])} filas con nombre: dentro de la secuencia {n(si['dentro_de_la_secuencia'])}, puntaje igual "
    f"{n(si['puntaje_calculado_igual_al_declarado'])}, >= 0,9 del máximo {n(si['puntaje_>=0.9_del_maximo'])}; exactos {n(si['coincidencia_exacta'])} ({d(100 * si['fraccion_exacta'], 2)} %), "
    f"{n(si['con_un_fallo_o_degenerada_incompleta'])} con un fallo o con una degenerada incompleta",
    ok(si["filas_con_nombre"] == si["dentro_de_la_secuencia"] == si["puntaje_calculado_igual_al_declarado"]
       == si["puntaje_>=0.9_del_maximo"]),
    [fu(LEC, "sitios", k) for k in si])
rg = {k: J(LEC, "regla", k) for k in ("registros", "elementos", "posiciones_PlantCARE", "posiciones_escaneo", "comunes",
                                      "solo_PlantCARE", "solo_escaneo", "mayor_fraccion_no_informada",
                                      "menor_fraccion_informada")}
J(LEC, "regla", "descripcion")
ctl("13", "La regla de búsqueda reproduce PlantCARE en los 189",
    "mismo conjunto (gen, elemento, posición) en el escaneo y en PlantCARE; umbral acotado: mayor fracción no informada "
    "< 0,9 <= menor fracción informada",
    f"{rg['registros']} registros, {rg['elementos']} elementos; posiciones {n(rg['posiciones_PlantCARE'])} (PlantCARE) y "
    f"{n(rg['posiciones_escaneo'])} (escaneo), comunes {n(rg['comunes'])}, solo en una {rg['solo_PlantCARE']} y "
    f"{rg['solo_escaneo']}; fracciones {d(rg['mayor_fraccion_no_informada'], 4)} < 0,9 <= {d(rg['menor_fraccion_informada'], 4)}",
    ok(rg["solo_PlantCARE"] == 0 and rg["solo_escaneo"] == 0 and rg["comunes"] == rg["posiciones_PlantCARE"]
       and rg["mayor_fraccion_no_informada"] < 0.9 <= rg["menor_fraccion_informada"]),
    [fu(LEC, "regla", k) for k in rg] + [fu(LEC, "regla", "descripcion")])
sn = {k: J(LEC, "filas_sin_nombre", k) for k in ("filas", "genes", "puntaje", "tetrameros_distintos", "decision")}
J(LEC, "filas_sin_nombre", "tetrameros")
ctl("13", "Filas sin nombre de PlantCARE", "",
    f"{n(sn['filas'])} filas en {sn['genes']} genes, {sn['tetrameros_distintos']} tetrámeros distintos, puntaje "
    f"{', '.join(d(p, 0) for p in sn['puntaje'])}; decisión: {sn['decision']}",
    "INFORMATIVO", [fu(LEC, "filas_sin_nombre", k) for k in sn] + [fu(LEC, "filas_sin_nombre", "tetrameros")])
t54, f54 = tabla("T54")
CLASES5 = ["Luz", "Desarrollo", "Fitohormonas", "Estrés", "Otros"]
cl54 = t54.clase.value_counts()
ctl("13", "Cada unidad tiene clase (5 clases del artículo 1, núcleo del promotor o sin función anotada)",
    "ninguna unidad sin clase; T54 = unidades del contraste; núcleo del promotor = TATA-box y CAAT-box",
    "; ".join(f"{c} {cl54.get(c, 0)}" for c in CLASES5 + ["Núcleo del promotor", "Sin función anotada"])
    + f"; total {len(t54)} (contraste: {J(CPC, 'n_unidades')}); núcleo: {', '.join(sorted(t54.unidad[t54.clase == 'Núcleo del promotor']))}",
    ok(t54.clase.notna().all() and set(t54.clase) <= set(CLASES5) | {"Núcleo del promotor", "Sin función anotada"}
       and len(t54) == J(CPC, "n_unidades")
       and set(t54.unidad[t54.clase == "Núcleo del promotor"]) == {"TATA-box", "CAAT-box"}),
    [f54, fu(CPC, "n_unidades")])
t55, f55 = tabla("T55")
t56, f56 = tabla("T56")
fus = {k: J(LEC, "clases", k) for k in ("filas_fundidas_misma_posicion", "sitios_tras_fusion")}
ctl("13", "Tablas de sitios y recuentos: filas",
    "T55 = filas con nombre menos las fundidas (misma unidad en la misma posición) = sitios del escaneo en los NAC; T56 = 189 genes",
    f"T55 {n(len(t55))} = {n(si['filas_con_nombre'])} - {fus['filas_fundidas_misma_posicion']} fundidas "
    f"(resumen: {n(fus['sitios_tras_fusion'])}); T56 {len(t56)}",
    ok(len(t55) == fus["sitios_tras_fusion"] == si["filas_con_nombre"] - fus["filas_fundidas_misma_posicion"]
       == J(CPC, "controles", "sitios_T55") and len(t56) == 189),
    [f55, f56, fu(LEC, "clases", "filas_fundidas_misma_posicion"), fu(LEC, "clases", "sitios_tras_fusion"),
     fu(CPC, "controles", "sitios_T55")])

# --- 14 contraste con el fondo emparejado por GC ---------------------------------------------------------
for f_, et in ((FGC, "1000"), (FGC10, "10 000")):
    c_ = {k: J(f_, "controles", k) for k in ("GC_tabla_igual_FASTA", "NAC_identicos_al_genomico", "intervalo_correcto",
                                              "sin_repetidos_en_conjunto")}
    p_ = {k: J(f_, k) for k in ("n_conjuntos", "semilla", "ancho_GC", "min_razon_fondo_NAC", "intervalos_con_NAC",
                                "fusiones", "razon_minima", "promotores_fondo_usados", "usos_max", "n_fondo")}
    ctl("14", f"Fondo emparejado por GC ({et} conjuntos de 189): construcción",
        "GC de la tabla = GC del FASTA; 189 NAC idénticos al FASTA genómico; cada promotor en el intervalo de su NAC; "
        "sin repetidos dentro de un conjunto",
        f"GC igual: {'sí' if c_['GC_tabla_igual_FASTA'] else 'no'}; NAC idénticos {c_['NAC_identicos_al_genomico']}; "
        f"intervalo correcto: {'sí' if c_['intervalo_correcto'] else 'no'}; sin repetidos: {'sí' if c_['sin_repetidos_en_conjunto'] else 'no'}; "
        f"semilla {p_['semilla']}; {p_['intervalos_con_NAC']} intervalos de {d(p_['ancho_GC'], 2)} ({p_['fusiones']} fusiones, razón "
        f"mínima {d(p_['razon_minima'], 1)}); {n(p_['promotores_fondo_usados'])} de {n(p_['n_fondo'])} promotores usados (máx. {p_['usos_max']} usos)",
        ok(c_["GC_tabla_igual_FASTA"] and c_["NAC_identicos_al_genomico"] == 189 and c_["intervalo_correcto"]
           and c_["sin_repetidos_en_conjunto"]),
        [fu(f_, "controles", k) for k in c_] + [fu(f_, k) for k in p_])
ctl("14", "La colección de 10 000 conjuntos amplía la original de 1000", "los primeros 1000 conjuntos son los originales",
    "iguales" if J(FGC10, "primeros_conjuntos_iguales_a_la_coleccion_original") else "distintos",
    ok(J(FGC10, "primeros_conjuntos_iguales_a_la_coleccion_original")),
    [fu(FGC10, "primeros_conjuntos_iguales_a_la_coleccion_original")])
ks = {k: J(FGC10, k) for k in ("KS_NAC_vs_fondo_elegido", "KS_NAC_vs_fondo_completo", "GC_medio", "dif_GC_medio_por_conjunto")}
J(FGC, "KS_NAC_vs_fondo_elegido"); J(FGC, "KS_NAC_vs_fondo_completo"); J(FGC, "GC_medio"); J(FGC, "dif_GC_medio_por_conjunto")
ctl("14", "Emparejamiento por GC (10 000 conjuntos)",
    "distribución de GC del fondo elegido indistinguible de la de los NAC (Kolmogorov-Smirnov p >= 0,05)",
    f"KS NAC frente al fondo elegido D {d(ks['KS_NAC_vs_fondo_elegido']['D'], 4)}, p {d(ks['KS_NAC_vs_fondo_elegido']['p'], 4)}; "
    f"frente al fondo completo D {d(ks['KS_NAC_vs_fondo_completo']['D'], 4)}, p {ec(ks['KS_NAC_vs_fondo_completo']['p'])}; "
    f"GC medio NAC {d(ks['GC_medio']['NAC'], 4)}, conjuntos {d(ks['GC_medio']['conjuntos'], 4)}, fondo completo "
    f"{d(ks['GC_medio']['fondo_completo'], 4)}",
    ok(ks["KS_NAC_vs_fondo_elegido"]["p"] >= 0.05),
    [fu(FGC10, k) for k in ks] + [fu(FGC, "KS_NAC_vs_fondo_elegido")])
c = {k: J(CPC, "controles", k) for k in ("secuencias_NAC_iguales_en_los_dos_FASTA", "sitios_NAC_escaneo", "sitios_T55",
                                         "sitios_iguales_a_T55", "recuentos_NAC_iguales_a_T56_unidades",
                                         "recuentos_NAC_iguales_a_T56_clases", "T54_sitios_y_genes_iguales",
                                         "todos_los_registros_con_sitio_en_los_NAC")}
ctl("14", "El escaneo del genoma reproduce PlantCARE en los 189 NAC",
    "mismos sitios que T55 (gen, unidad, posición), mismos recuentos que T56 por unidad y por clase, mismos sitios y genes "
    "por unidad que T54; todos los registros con sitio en los NAC",
    f"secuencias NAC iguales en los dos FASTA {c['secuencias_NAC_iguales_en_los_dos_FASTA']}; sitios {n(c['sitios_NAC_escaneo'])} "
    f"(escaneo) y {n(c['sitios_T55'])} (T55); " + "; ".join(
        f"{k.replace('_', ' ')}: {'sí' if c[k] else 'no'}" for k in ("sitios_iguales_a_T55", "recuentos_NAC_iguales_a_T56_unidades",
                                                                     "recuentos_NAC_iguales_a_T56_clases", "T54_sitios_y_genes_iguales",
                                                                     "todos_los_registros_con_sitio_en_los_NAC")),
    ok(c["secuencias_NAC_iguales_en_los_dos_FASTA"] == 189 and c["sitios_NAC_escaneo"] == c["sitios_T55"]
       and all(c[k] for k in ("sitios_iguales_a_T55", "recuentos_NAC_iguales_a_T56_unidades",
                              "recuentos_NAC_iguales_a_T56_clases", "T54_sitios_y_genes_iguales",
                              "todos_los_registros_con_sitio_en_los_NAC"))),
    [fu(CPC, "controles", k) for k in c])
e = {k: J(CPC, "controles", k) for k in ("exactos_T55", "exactos_escaneo", "exactos_T55_incluidos_en_el_escaneo",
                                         "exactos_escaneo_no_marcados_en_T55",
                                         "exactos_no_marcados_explicados_por_otro_registro_de_la_unidad",
                                         "exactos_matriz_igual_al_escaneo", "T54_sitios_exactos_iguales_salvo_esos_sitios",
                                         "exactos_no_marcados_por_unidad")}
J(CPC, "controles", "exactos_no_marcados_en_T55")
ctl("14", "Sitios exactos del escaneo frente a T55",
    "los exactos de T55 están en el escaneo; los que añade el escaneo son sitios donde PlantCARE informó un registro no "
    "exacto y otro registro de la unidad es exacto; matriz de exactos = escaneo; T54 igual salvo esos sitios",
    f"exactos {n(e['exactos_T55'])} (T55) y {n(e['exactos_escaneo'])} (escaneo); añadidos {e['exactos_escaneo_no_marcados_en_T55']} "
    f"(" + ", ".join(f"{k} {v}" for k, v in e["exactos_no_marcados_por_unidad"].items()) + "), explicados: "
    f"{'sí' if e['exactos_no_marcados_explicados_por_otro_registro_de_la_unidad'] else 'no'}",
    ok(e["exactos_T55_incluidos_en_el_escaneo"] and e["exactos_no_marcados_explicados_por_otro_registro_de_la_unidad"]
       and e["exactos_matriz_igual_al_escaneo"] and e["T54_sitios_exactos_iguales_salvo_esos_sitios"]
       and e["exactos_escaneo"] - e["exactos_T55"] == e["exactos_escaneo_no_marcados_en_T55"]),
    [fu(CPC, "controles", k) for k in e] + [fu(CPC, "controles", "exactos_no_marcados_en_T55")])
esc = {k: J(CPC, "escaneo", k) for k in ("reutilizado", "promotores", "trozos", "sitios_genoma", "exactos_genoma")}
for k in JS[CPC]["escaneo"]:
    J(CPC, "escaneo", k)
ctl("14", "Escaneo de PlantCARE en el genoma", "",
    f"{n(esc['promotores'])} promotores en {esc['trozos']} trozos; {n(esc['sitios_genoma'])} sitios ({n(esc['exactos_genoma'])} "
    f"exactos); {'reutilizado de la corrida anterior (mismos registros; lo validan los dos controles previos)' if esc['reutilizado'] else 'calculado en esta corrida'}",
    "INFORMATIVO", [fu(CPC, "escaneo", k) for k in esc])


def calib(paso, que, frac_media, frac_rango, n_pruebas, lista, fuentes, extra=""):
    N = len(lista) * n_pruebas
    S = int(round(frac_media * N))
    p_exc = binom.sf(S - 1, N, 0.05) if S > 0 else 1.0
    k = sum(1 for x in lista if x > 0)
    p_q = binom.sf(k - 1, len(lista), 0.05) if k > 0 else 1.0
    ctl(paso, f"Calibración: {que}",
        "fracción de p < 0,05 no mayor que 0,05 (binomial unilateral de exceso, p >= 0,05) y conjuntos falsos con algún "
        "q < 0,05 no más del 5 % (binomial unilateral, p >= 0,05)",
        f"fracción de p < 0,05 {d(frac_media, 4)}" + (f" (por conjunto {d(frac_rango[0], 3)}-{d(frac_rango[1], 3)})" if frac_rango else "")
        + f" en ~{n(N)} pruebas (exceso p = {d(p_exc, 3)}); conjuntos con algún q < 0,05: {k}/{len(lista)} (p = {d(p_q, 3)})" + extra,
        ok(p_exc >= 0.05 and p_q >= 0.05), fuentes)
    return dict(frac=frac_media, pruebas=N, p_exceso=p_exc, conjuntos_q=k, de=len(lista), p_q=p_q)


PROPIO["calibracion"] = {}
cs = J(CPC, "calibracion", "u_sitios_sin_seleccion")
PROPIO["calibracion"]["PlantCARE_unidades"] = calib(
    "14", f"PlantCARE, unidades sin selección ({J(CPC, 'n_calibracion')} falsos NAC)", cs["fraccion_p<0.05_media"],
    cs["fraccion_p<0.05_rango"], J(CPC, "n_unidades"), cs["n_q<0.05_por_conjunto"],
    [fu(CPC, "calibracion", "u_sitios_sin_seleccion"), fu(CPC, "n_calibracion")])
cs = J(CPC, "calibracion", "c_sitios_sin_seleccion")
PROPIO["calibracion"]["PlantCARE_clases"] = calib(
    "14", "PlantCARE, clases sin selección", cs["fraccion_p<0.05_media"], cs["fraccion_p<0.05_rango"],
    J(CPC, "resolucion", "c_sitios", "evaluables"), cs["n_q<0.05_por_conjunto"],
    [fu(CPC, "calibracion", "c_sitios_sin_seleccion"), fu(CPC, "resolucion", "c_sitios", "evaluables")])
cs = J(CPC, "calibracion", "u_sitios_con_seleccion")
estr = cs["por_estrato_de_sitios"]
PROPIO["calibracion"]["PlantCARE_unidades_con_seleccion"] = calib(
    "14", "PlantCARE, unidades con la misma selección por presencia que el análisis real", cs["fraccion_p<0.05"], None,
    cs["unidades_probadas_media"], cs["n_q<0.05_por_conjunto"], [fu(CPC, "calibracion", "u_sitios_con_seleccion")],
    extra=f"; condicionada a presencia {d(cs['fraccion_p_condicional<0.05'], 4)}, conjuntos con algún q < 0,05 "
          f"{sum(1 for x in cs['n_q_condicional<0.05_por_conjunto'] if x > 0)}/{len(cs['n_q_condicional<0.05_por_conjunto'])}; por estrato "
          "de sitios en el falso NAC: " + ", ".join(f"{k} {d(v['fraccion_p<0.05'], 3)}" for k, v in estr.items()))
rv = {k: J(CPC, "resolucion", k) for k in JS[CPC]["resolucion"]}
c1 = {k: J(CPC, "control_N1000", k) for k in JS[CPC]["control_N1000"]}
ctl("14", "Resolución: 1000 frente a 10 000 conjuntos de fondo", "",
    f"unidades con q < 0,05: {c1['u_sitios']['q<0.05_N1000']} con 1000 y {c1['u_sitios']['q<0.05_N10000']} con 10 000 "
    f"(Spearman de p {d(c1['u_sitios']['spearman_p'], 4)}); clases: {c1['c_sitios']['q<0.05_N1000']} y {c1['c_sitios']['q<0.05_N10000']}; "
    f"con 10 000, q mínima de unidades {d(rv['u_sitios']['q_minima'], 4)} y de clases {d(rv['c_sitios']['q_minima'], 4)}",
    "INFORMATIVO", [fu(CPC, "control_N1000"), fu(CPC, "resolucion")])

# --- 14J JASPAR 2026 + FIMO ---------------------------------------------------------------------------
md5j = J(JPREP, "md5")
for k in JS[JESC]["md5"]:
    J(JESC, "md5", k)
rutas_md5 = {"promotores_2kb_todos.fa": os.path.join(E3, "promotores", "promotores_2kb_todos.fa"),
             "NAC189_promotores_2kb.fa": os.path.join(E3, "promotores", "NAC189_promotores_2kb.fa"),
             "fondo_markov0.txt": os.path.join(E3, "cis_jaspar", "fondo_markov0.txt"),
             "JASPAR2026_CORE_plants_non-redundant_pfms_meme.txt":
                 os.path.join(E3, "cis_jaspar", "matrices", "JASPAR2026_CORE_plants_non-redundant_pfms_meme.txt"),
             "JASPAR2026_CORE_plants_non-redundant_pfms_transfac.txt":
                 os.path.join(E3, "cis_jaspar", "matrices", "JASPAR2026_CORE_plants_non-redundant_pfms_transfac.txt")}
iguales = {f: md5_fichero(p) == md5j[f] == JS[JESC]["md5"][f] for f, p in rutas_md5.items()}
ctl("14J", "Entradas del escaneo sin cambios desde que se escanearon",
    "md5 actual de los promotores, matrices y fondo de Markov = md5 registrado al preparar el escaneo",
    f"{sum(iguales.values())}/{len(iguales)} iguales ({', '.join(iguales)})", ok(all(iguales.values())),
    [fu(JPREP, "md5"), fu(JESC, "md5")])
jx = {k: J(JESC, k) for k in ("control_NAC_identico", "control_NAC_nota", "NAC_lineas_con_hebra_distinta_entre_corridas",
                              "NAC_de_ellas_palindromicas", "NAC_sitios_palindromicos", "sitios_NAC", "n_trozos",
                              "promotores_sin_sitios")}
ctl("14J", "Sitios NAC del escaneo por trozos = escaneo aparte del FASTA de los 189",
    "idénticos línea a línea (salvo la hebra elegida al azar por FIMO --max-strand en posiciones con empate)",
    f"idénticos: {'sí' if jx['control_NAC_identico'] else 'no'}; {n(jx['sitios_NAC'])} sitios NAC; líneas con hebra distinta "
    f"{jx['NAC_lineas_con_hebra_distinta_entre_corridas']} ({jx['NAC_de_ellas_palindromicas']} palindrómicas); "
    f"{jx['n_trozos']} trozos; promotores sin sitios {jx['promotores_sin_sitios']}",
    ok(jx["control_NAC_identico"] and jx["promotores_sin_sitios"] == 0), [fu(JESC, k) for k in jx])
cj = J(JCON, "calibracion", "mot_sitios")
PROPIO["calibracion"]["JASPAR_motivos"] = calib(
    "14J", f"JASPAR, motivos ({J(JCON, 'n_calibracion')} falsos NAC)", cj["fraccion_p<0.05_media"], cj["fraccion_p<0.05_rango"],
    J(JCON, "resolucion", "mot_sitios", "evaluables"), cj["n_q<0.05_por_conjunto"],
    [fu(JCON, "calibracion", "mot_sitios"), fu(JCON, "resolucion", "mot_sitios", "evaluables"), fu(JCON, "n_calibracion")])
cj = J(JCON, "calibracion", "fam_regiones")
PROPIO["calibracion"]["JASPAR_familias"] = calib(
    "14J", "JASPAR, familias", cj["fraccion_p<0.05_media"], cj["fraccion_p<0.05_rango"],
    J(JCON, "resolucion", "fam_regiones", "evaluables"), cj["n_q<0.05_por_conjunto"],
    [fu(JCON, "calibracion", "fam_regiones"), fu(JCON, "resolucion", "fam_regiones", "evaluables")])
t51, f51 = tabla("T51")
t52, f52 = tabla("T52")
t53, f53 = tabla("T53")
ctl("14J", "Tablas de JASPAR: filas", "T51 = motivos; T52 = familias; T53 = 189 genes x familias",
    f"T51 {len(t51)} = {J(JCON, 'n_motivos')}; T52 {len(t52)} = {J(JCON, 'n_familias')}; T53 {n(len(t53))} = 189 x {len(t52)}",
    ok(len(t51) == J(JCON, "n_motivos") == J(JESC, "n_motivos") and len(t52) == J(JCON, "n_familias")
       and len(t53) == 189 * len(t52)),
    [f51, f52, f53, fu(JCON, "n_motivos"), fu(JCON, "n_familias"), fu(JESC, "n_motivos")])

# --- 15 figuras y tablas de elementos cis --------------------------------------------------------------
for b_ in ("F18_cis_PlantCARE_por_gen", "F19_cis_enriquecimiento_fondo_GC"):
    control_figura("15", b_)
control_log("15", "entrega3_20_figuras_cis.log", "cifras -> _cifras/e3_cis.csv",
            "e3_18_figuras_cis.R (genes y subfamilias de T56, orden del árbol de F06, medias del fondo)")
t57, f57 = tabla("T57")
t58, f58 = tabla("T58")
t59, f59 = tabla("T59")
ctl("15", "Tablas de enriquecimiento de PlantCARE: filas", "T57 = unidades; T58 = 7 grupos; T59 = 189 genes x 7 grupos",
    f"T57 {len(t57)}; T58 {len(t58)}; T59 {n(len(t59))}",
    ok(len(t57) == J(CPC, "n_unidades") and len(t58) == 7 and len(t59) == 189 * 7), [f57, f58, f59])

# =================================================================================================
# barrera: todo booleano de los resúmenes tiene que estar citado
booleanos = {(f, r): v for f in RES for r, v in _hojas(JS[f]) if isinstance(v, bool)}
faltan = sorted(set(booleanos) - USADAS, key=str)
if faltan:
    sys.exit("Booleanos sin control en T60:\n" + "\n".join(f"  {f}: {'.'.join(map(str, r))}" for f, r in faltan))

t60 = pd.DataFrame(FILAS)[["id", "paso", "control", "criterio", "observado", "resultado", "fuente"]]
t60.to_csv(os.path.join(TAB, "T60_controles_entrega3.csv"), index=False, encoding="utf-8")
cuenta = t60.resultado.value_counts().to_dict()
res = dict(fecha=f"{datetime.datetime.now():%Y-%m-%d %H:%M:%S}", guion="e3_19_controles.py",
           controles=len(t60), por_resultado=cuenta, booleanos_en_resumenes=len(booleanos),
           booleanos_falsos=sorted(fu(f, *r) for (f, r), v in booleanos.items() if v is False),
           resumenes=RES, versiones=dict(python=platform.python_version(), pandas=pd.__version__,
                                         scipy=scipy.__version__, pillow=PIL.__version__),
           **PROPIO)
json.dump(res, open(os.path.join(VAL, "resumen_controles.json"), "w", encoding="utf-8"), ensure_ascii=False, indent=1,
          default=lambda o: o.item() if hasattr(o, "item") else str(o))
print(f"T60: {len(t60)} controles; {cuenta}; {len(booleanos)} booleanos de {len(RES)} resúmenes, todos citados", flush=True)
for r in t60[t60.resultado == "FALLA"].itertuples():
    print(f"  FALLA {r.id} {r.control}: {r.observado}", flush=True)
