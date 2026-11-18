"""Entrega 3, paso 3.1: homólogos de Arabidopsis para el análisis GO.

Entradas
  datos/entrega3/isoformas/{Egl,Ath}_representativas.tsv   una proteína por gen (paso 1.2)
  rnaseq_udec/ref/GCA_046226855.1_ASM4622685v1_protein.faa.gz   proteínas de E. globulus
  referencia/Athaliana_TAIR10.1/protein.faa                      proteínas de TAIR10.1 (GCF_000001735.4)
  datos/entrega3/go/uniprot_NAC_Ath_UP000006548.tsv             NAC de Arabidopsis en UniProt 2026_03
      (proteoma de referencia UP000006548 con PF02365, IPR003441 o PS51005)
  tablas/T01_censo_maestro_189.csv

Proceso
  1. Proteomas representativos: 36 387 genes de E. globulus y 27 444 genes nucleares de
     Arabidopsis. Se excluyen los 118 genes de cloroplasto (85, ArthCp...) y mitocondria
     (33, DA397_mgp...), que no tienen código AGI en RefSeq ni entran en org.At.tair.db.
  2. DIAMOND blastp --ultra-sensitive, e <= 1e-5, hasta 5 dianas por consulta.
  3. Mejor acierto por gen de E. globulus (mayor bitscore) entre los que cumplen
     e <= 1e-5, identidad >= 30 % y cobertura de la consulta >= 50 %.

Salidas
  datos/entrega3/go/Egl_vs_Ath.tsv                 salida de DIAMOND (formato 6 ampliado)
  datos/entrega3/go/Egl_mejor_homologo_Ath.tsv     mejor homólogo de cada gen de E. globulus
  tablas/T41_homologos_Ath_NAC189.csv              EglNAC -> AGI -> ANAC
  datos/entrega3/go/resumen_homologos.json
"""
import gzip, json, os, re, subprocess, sys
import pandas as pd

RAIZ = (os.environ.get("NAC_ROOT", ".") + "")
A = f"{RAIZ}/ANALISIS_Eglobulus"
E3 = f"{A}/datos/entrega3"
GO = f"{E3}/go"
ISO = f"{E3}/isoformas"
BIN = (os.environ.get("CONDA_ENVS", os.path.expanduser("~/conda/envs")) + "/nac3/bin")
PROT_EGL = f"{RAIZ}/rnaseq_udec/ref/GCA_046226855.1_ASM4622685v1_protein.faa.gz"
PROT_ATH = f"{A}/referencia/Athaliana_TAIR10.1/protein.faa"
UNIPROT_NAC = f"{GO}/uniprot_NAC_Ath_UP000006548.tsv"
T01 = f"{A}/tablas/T01_censo_maestro_189.csv"
T41 = f"{A}/tablas/T41_homologos_Ath_NAC189.csv"

EVALUE, IDENT, COB = 1e-5, 30.0, 50.0
COLS = "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovhsp scovhsp".split()
AGI = re.compile(r"^AT[1-5]G\d{5}$")


def fasta(p):
    d, h = {}, None
    with (gzip.open if p.endswith(".gz") else open)(p, "rt") as fh:
        for l in fh:
            if l.startswith(">"):
                h = l[1:].split()[0]
                d[h] = []
            else:
                d[h].append(l.strip())
    return {k: "".join(v) for k, v in d.items()}


def proteomas():
    egl = pd.read_csv(f"{ISO}/Egl_representativas.tsv", sep="\t")
    ath = pd.read_csv(f"{ISO}/Ath_representativas.tsv", sep="\t")
    org = ath[~ath.gen.str.match(AGI)]
    # cloroplasto NC_000932.1 (ArthCp...) y mitocondria NC_037304.1 (DA397_mgp...)
    assert org.seqid.isin(["NC_000932.1", "NC_037304.1"]).all(), org.gen.head().tolist()
    ath = ath[ath.gen.str.match(AGI)]
    for nombre, rep, prot in (("Egl_todas.pep", egl, fasta(PROT_EGL)), ("Ath_nucleares.pep", ath, fasta(PROT_ATH))):
        with open(f"{GO}/{nombre}", "w") as fh:
            for g, p in zip(rep.gen, rep.proteina):
                fh.write(f">{g}\n{prot[p].rstrip('*')}\n")
    return egl, ath, len(org)


def diamond(forzar=False):
    salida = f"{GO}/Egl_vs_Ath.tsv"
    if os.path.exists(salida) and not forzar:
        return salida
    os.makedirs(f"{GO}/tmp", exist_ok=True)
    subprocess.run([f"{BIN}/diamond", "makedb", "--in", f"{GO}/Ath_nucleares.pep", "-d", f"{GO}/Ath_nucleares",
                    "--threads", "12"], check=True)
    subprocess.run([f"{BIN}/diamond", "blastp", "-q", f"{GO}/Egl_todas.pep", "-d", f"{GO}/Ath_nucleares", "-o", salida,
                    "--ultra-sensitive", "--evalue", str(EVALUE), "--max-target-seqs", "5", "--threads", "12",
                    "--tmpdir", f"{GO}/tmp", "--outfmt", "6", *COLS], check=True)
    return salida


def nac_arabidopsis():
    up = pd.read_csv(UNIPROT_NAC, sep="\t", dtype=str, keep_default_na=False)
    texto = up["Araport"] + " " + up["Gene Names (ordered locus)"] + " " + up["Gene Names"] + " " + up["Protein names"]
    up["agi"] = texto.map(lambda s: sorted({a.upper() for a in re.findall(r"\b(At[1-5CM]g\d{5})\b", s, flags=re.I)}))
    assert (up.agi.str.len() == 1).all()
    up["agi"] = up.agi.str[0]
    up["revisada"] = up.Reviewed.eq("reviewed")
    up = up.sort_values(["agi", "revisada"], ascending=[True, False]).drop_duplicates("agi")

    def anac(nombres):
        # Numeración de Ooka et al. (2003). UniProt la escribe como ANAC059, NAC059 o NAC59;
        # NAC1, NAC2, NAC3 y NAC6 son nombres históricos, no números ANAC.
        m = re.search(r"\bA?NAC(\d{3})\b", nombres) or re.search(r"\bA?NAC([1-9]\d)\b", nombres)
        return f"ANAC{int(m.group(1)):03d}" if m else ""

    def principal(nombres):
        # primer nombre que no sea código AGI, nombre de clon (F14D16.19, T13D8_13) ni A?NAC + número
        for t in nombres.split():
            if not re.match(r"^(At[1-5]g\d{5}|[A-Z0-9]{2,8}[._]\d+|A?NAC\d+)$", t, flags=re.I):
                return t
        return ""

    up["ANAC"] = up["Gene Names"].map(anac)
    up["nombre_Ath"] = up["Gene Names"].map(principal)
    return up[["agi", "ANAC", "nombre_Ath", "Gene Names", "Entry", "revisada"]].rename(
        columns={"Gene Names": "nombres_UniProt", "Entry": "UniProt"})


def mejor_acierto(hits):
    ok = hits[(hits.evalue <= EVALUE) & (hits.pident >= IDENT) & (hits.qcovhsp >= COB)]
    orden = ["qseqid", "bitscore", "evalue", "pident"]
    return ok.sort_values(orden, ascending=[True, False, True, False]).drop_duplicates("qseqid")


def main():
    forzar = "--forzar" in sys.argv
    egl, ath, n_org = proteomas()
    hits = pd.read_csv(diamond(forzar), sep="\t", names=COLS)
    nac = nac_arabidopsis()
    nac_set = set(nac.agi)
    mejor = mejor_acierto(hits)
    mejor = mejor.assign(es_NAC_Ath=mejor.sseqid.isin(nac_set))
    mejor[["qseqid", "sseqid", "pident", "length", "evalue", "bitscore", "qlen", "slen", "qcovhsp", "scovhsp",
           "es_NAC_Ath"]].rename(columns={"qseqid": "gen_Egl", "sseqid": "AGI"}).to_csv(
        f"{GO}/Egl_mejor_homologo_Ath.tsv", sep="\t", index=False)

    t01 = pd.read_csv(T01, dtype=str, keep_default_na=False)
    assert t01.locus.isin(egl.gen).all()
    rep = egl.set_index("gen")
    b = mejor.set_index("qseqid")
    sin_filtro = hits.sort_values(["qseqid", "bitscore"], ascending=[True, False]).drop_duplicates("qseqid").set_index("qseqid")
    hits_nac = hits[hits.sseqid.isin(nac_set)].sort_values(["qseqid", "bitscore"], ascending=[True, False]).drop_duplicates("qseqid").set_index("qseqid")
    info = nac.set_index("agi")
    filas = []
    for _, r in t01.iterrows():
        g = r.locus
        f = dict(simbolo=r.simbolo, locus=g, subfamilia=r.get("subfamilia_final", ""), proteina_Egl=rep.at[g, "proteina"],
                 largo_aa=int(rep.at[g, "largo_aa"]))
        if g in b.index:
            h = b.loc[g]
            f.update(AGI=h.sseqid, ANAC=info.ANAC.get(h.sseqid, ""), nombre_Ath=info.nombre_Ath.get(h.sseqid, ""),
                     es_NAC_Ath=bool(h.sseqid in nac_set), identidad=h.pident, cobertura_consulta=h.qcovhsp,
                     cobertura_diana=h.scovhsp, evalue=h.evalue, bitscore=h.bitscore, pasa_filtros=True)
        else:
            f.update(AGI="", ANAC="", nombre_Ath="", es_NAC_Ath=False, identidad=None, cobertura_consulta=None,
                     cobertura_diana=None, evalue=None, bitscore=None, pasa_filtros=False)
        s = sin_filtro.loc[g] if g in sin_filtro.index else None
        f.update(mejor_AGI_sin_filtros="" if s is None else s.sseqid,
                 identidad_sin_filtros=None if s is None else s.pident,
                 cobertura_sin_filtros=None if s is None else s.qcovhsp)
        n = hits_nac.loc[g] if g in hits_nac.index else None
        f.update(mejor_NAC_Ath_entre_dianas="" if n is None else n.sseqid,
                 ANAC_entre_dianas="" if n is None else info.ANAC.get(n.sseqid, ""),
                 identidad_NAC_entre_dianas=None if n is None else n.pident)
        filas.append(f)
    t41 = pd.DataFrame(filas)
    t41.to_csv(T41, index=False)

    res = dict(genes_Egl=len(egl), genes_Ath_nucleares=len(ath), genes_Ath_organulos_excluidos=n_org,
               consultas_con_acierto=int(hits.qseqid.nunique()), genes_Egl_con_homologo=int(len(mejor)),
               AGI_distintos_universo=int(mejor.sseqid.nunique()), NAC_Ath_UniProt=len(nac_set),
               NAC_Ath_con_numero_ANAC=int(nac.ANAC.ne("").sum()),
               NAC189_con_homologo=int(t41.pasa_filtros.sum()), NAC189_mejor_es_NAC=int(t41.es_NAC_Ath.sum()),
               NAC189_AGI_distintos=int(t41.loc[t41.pasa_filtros, "AGI"].nunique()),
               filtros=dict(evalue=EVALUE, identidad_min=IDENT, cobertura_consulta_min=COB))
    json.dump(res, open(f"{GO}/resumen_homologos.json", "w"), ensure_ascii=False, indent=1)
    print(json.dumps(res, ensure_ascii=False))


if __name__ == "__main__":
    main()
