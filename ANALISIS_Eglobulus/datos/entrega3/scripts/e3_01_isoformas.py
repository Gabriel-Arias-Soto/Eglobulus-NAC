"""Entrega 3, paso 1.3: una isoforma representativa por gen en las cuatro especies.

Regla: en los 189 NAC de E. globulus se usa la proteína fijada en T01; en el resto de
genes, la proteína más larga (empate: CDS completa antes que parcial; después, el
identificador menor). Cada representativa se verifica traduciendo su CDS.

Salidas en datos/entrega3/isoformas/:
  <sp>_representativas.tsv  tabla por gen (coordenadas, proteína, controles)
  <sp>.pep / <sp>.cds       FASTA con el identificador del gen
  <sp>.bed                  genes en cromosomas nucleares (formato JCVI)
  resumen_isoformas.json    controles por especie
"""
import gzip, os, re, json
import pandas as pd

BASE = (os.environ.get("NAC_ROOT", ".") + "")
A = os.path.join(BASE, "ANALISIS_Eglobulus")
REF = os.path.join(A, "referencia")
OUT = os.path.join(A, "datos", "entrega3", "isoformas")
RQ = os.path.join(BASE, "rnaseq_udec", "ref")

ESPECIES = {
    "Egl": dict(nombre="Eucalyptus globulus", acc="GCA_046226855.1",
                gtf=os.path.join(RQ, "GCA_046226855.1_ASM4622685v1_genomic.gtf.gz"),
                prot=os.path.join(RQ, "GCA_046226855.1_ASM4622685v1_protein.faa.gz"),
                cds=os.path.join(RQ, "GCA_046226855.1_ASM4622685v1_cds_from_genomic.fna.gz"),
                genes=os.path.join(A, "datos", "genes_todos.tsv"),
                sr=os.path.join(REF, "Eglobulus_ASM4622685v1", "sequence_report.json")),
    "Egr": dict(nombre="Eucalyptus grandis", acc="GCA_000612305.2", carpeta="Egrandis_v2"),
    "Ptr": dict(nombre="Populus trichocarpa", acc="GCF_000002775.5", carpeta="Ptrichocarpa_v4.1"),
    "Ath": dict(nombre="Arabidopsis thaliana", acc="GCF_000001735.4", carpeta="Athaliana_TAIR10.1"),
}

BASES = "TCAG"
AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON = {a + b + c: AA[16 * i + 4 * j + k] for i, a in enumerate(BASES) for j, b in enumerate(BASES) for k, c in enumerate(BASES)}
INICIOS_ALT = {"TTG", "CTG", "GTG", "ATT", "ATC", "ATA", "ACG"}


def abrir(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)


def leer_fasta(p):
    seqs, h, buf = {}, None, []
    with abrir(p) as fh:
        for l in fh:
            if l.startswith(">"):
                if h is not None:
                    seqs[h] = "".join(buf)
                h, buf = l[1:].rstrip("\n"), []
            else:
                buf.append(l.strip())
    if h is not None:
        seqs[h] = "".join(buf)
    return seqs


def traducir(nt, inicio=1):
    s = nt[inicio - 1:].upper()
    return "".join(CODON.get(s[i:i + 3], "X") for i in range(0, len(s) - 2, 3))


def verificar(cds, prot, marco=1):
    """(traducción coincide, inicio alternativo, paradas internas, selenocisteínas)"""
    tr = traducir(cds, marco)
    if tr.endswith("*"):
        tr = tr[:-1]
    sec = sum(1 for x, y in zip(tr, prot) if x == "*" and y == "U")
    internas = tr.count("*") - sec
    if tr == prot:
        return True, False, internas, sec
    if len(tr) != len(prot):
        return False, False, internas, sec
    dif = [i for i, (x, y) in enumerate(zip(tr, prot)) if x != y]
    alt = dif[:1] == [0] and prot[0] == "M" and cds[marco - 1:marco + 2].upper() in INICIOS_ALT
    resto = dif[1:] if alt else dif
    ok = all(tr[i] == "X" or prot[i] == "X" or (tr[i] == "*" and prot[i] == "U") for i in resto)
    return ok, alt, internas, sec


def attrs_cabecera(h):
    return dict(re.findall(r"\[(\w+)=([^\]]*)\]", h))


def cds_por_proteina(p):
    """FASTA de CDS de NCBI -> {protein_id: (secuencia, marco, parcial)}"""
    out = {}
    for h, s in leer_fasta(p).items():
        a = attrs_cabecera(h)
        pid = a.get("protein_id")
        if pid and a.get("pseudo") != "true" and pid not in out:
            out[pid] = (s, int(a.get("frame", 1)), a.get("partial", ""))
    return out


def proteinas(p):
    return {h.split()[0]: s for h, s in leer_fasta(p).items()}


def attrs_gff3(s):
    return dict(kv.split("=", 1) for kv in s.strip().split(";") if "=" in kv)


def parse_gff3(p):
    """GFF3 de NCBI -> (genes {ID: info}, {protein_id: ID del gen})"""
    genes, padre, cds_padre = {}, {}, {}
    with open(p) as fh:
        for l in fh:
            if l.startswith("#"):
                continue
            c = l.rstrip("\n").split("\t")
            if len(c) < 9 or c[2] not in ("gene", "pseudogene", "mRNA", "transcript", "CDS"):
                continue
            a = attrs_gff3(c[8])
            if c[2] in ("gene", "pseudogene"):
                genes[a["ID"]] = dict(seqid=c[0], ini=int(c[3]), fin=int(c[4]), hebra=c[6],
                                      gen=a.get("locus_tag") or a.get("Name") or a["ID"],
                                      biotipo=a.get("gene_biotype", c[2]), nombre=a.get("Name", ""))
            elif c[2] in ("mRNA", "transcript"):
                padre[a["ID"]] = a.get("Parent", "").split(",")[0]
            else:
                pid = a.get("protein_id")
                if pid and pid not in cds_padre:
                    cds_padre[pid] = a.get("Parent", "").split(",")[0]
    p2g = {}
    for pid, par in cds_padre.items():
        g = padre.get(par, par)
        if g in genes:
            p2g[pid] = g
    return genes, p2g


def mapa_cromosomas_jsonl(p, sp):
    m = {}
    for l in open(p):
        s = json.loads(l)
        if s.get("role") == "assembled-molecule" and s.get("chrName", "").isdigit():
            for k in ("genbankAccession", "refseqAccession"):
                if s.get(k):
                    m[s[k]] = (f"{sp}{int(s['chrName']):02d}", int(s["length"]))
    return m


def elegir(df, fijadas=None):
    df = df.sort_values(["gen", "largo_aa", "parcial", "proteina"], ascending=[True, False, True, True])
    rep = df.drop_duplicates("gen").set_index("gen")
    n_iso = df.groupby("gen").size()
    if fijadas:
        for g, pid in fijadas.items():
            fila = df[(df.gen == g) & (df.proteina == pid)]
            assert len(fila) == 1, (g, pid, len(fila))
            rep.loc[g] = fila.iloc[0].drop("gen")
    rep["n_isoformas"] = n_iso.reindex(rep.index).values
    return rep.reset_index()


def tabla_ncbi(sp, cfg):
    d = os.path.join(REF, cfg["carpeta"])
    genes, p2g = parse_gff3(os.path.join(d, "genomic.gff"))
    prot = proteinas(os.path.join(d, "protein.faa"))
    cds = cds_por_proteina(os.path.join(d, "cds_from_genomic.fna"))
    crom = mapa_cromosomas_jsonl(os.path.join(d, "sequence_report.jsonl"), sp)
    filas = []
    for pid, gid in p2g.items():
        if pid not in prot or pid not in cds:
            continue
        g = genes[gid]
        s, marco, parcial = cds[pid]
        filas.append(dict(gen=g["gen"], proteina=pid, largo_aa=len(prot[pid]), parcial=int(bool(parcial)),
                          parcial_tipo=parcial, marco=marco, seqid=g["seqid"], ini=g["ini"], fin=g["fin"], hebra=g["hebra"]))
    df = pd.DataFrame(filas)
    rep = elegir(df)
    rep["cromosoma"] = rep.seqid.map(lambda x: crom.get(x, (None, 0))[0])
    return rep, prot, cds, crom


def tabla_egl(cfg, t01):
    prot = proteinas(cfg["prot"])
    cds = cds_por_proteina(cfg["cds"])
    pg, segs = {}, {}
    with abrir(cfg["gtf"]) as fh:
        for l in fh:
            c = l.rstrip("\n").split("\t")
            if len(c) < 9 or c[2] != "CDS":
                continue
            a = dict(re.findall(r'(\S+) "([^"]*)"', c[8]))
            pid = a.get("protein_id")
            if pid:
                pg.setdefault(pid, (a["gene_id"], a.get("transcript_id", "")))
                segs.setdefault(pid, []).append((int(c[3]), int(c[4])))
    gt = pd.read_csv(cfg["genes"], sep="\t", header=None, names=["seqid", "ini", "fin", "hebra", "gen", "biotipo"])
    gt = gt.drop_duplicates("gen").set_index("gen")
    sr = json.load(open(cfg["sr"]))["reports"]
    crom = {s["genbank_accession"]: (f"Egl{int(s['chr_name']):02d}", int(s["length"])) for s in sr
            if s.get("role") == "assembled-molecule" and s.get("chr_name", "").isdigit()}
    filas = []
    for pid, (g, tx) in pg.items():
        if pid not in prot or pid not in cds or g not in gt.index:
            continue
        s, marco, parcial = cds[pid]
        gi = gt.loc[g]
        sg = segs[pid]
        ini_cds = min(a for a, b in sg) if gi.hebra == "+" else max(b for a, b in sg)
        filas.append(dict(gen=g, proteina=pid, transcrito=tx, largo_aa=len(prot[pid]), parcial=int(bool(parcial)),
                          parcial_tipo=parcial, marco=marco, seqid=gi.seqid, ini=int(gi.ini), fin=int(gi.fin),
                          hebra=gi.hebra, cds_inicio=ini_cds))
    df = pd.DataFrame(filas)
    fijadas = dict(zip(t01.locus, t01.proteina))
    rep = elegir(df, fijadas)
    rep["cromosoma"] = rep.seqid.map(lambda x: crom.get(x, (None, 0))[0])
    return rep, prot, cds, crom


def controles(rep, prot, cds):
    ok, alt, internas, sec, atg = [], [], [], [], []
    for _, r in rep.iterrows():
        s, marco, parcial = cds[r.proteina]
        o, a, i, u = verificar(s, prot[r.proteina], marco)
        ok.append(o); alt.append(a); internas.append(i); sec.append(u)
        atg.append(s[:3].upper() == "ATG" if "5'" not in parcial else None)
    rep["traduccion_ok"], rep["inicio_alternativo"], rep["paradas_internas"], rep["selenocisteinas"], rep["empieza_ATG"] = ok, alt, internas, sec, atg
    return rep


def escribir(sp, rep, prot, cds):
    os.makedirs(OUT, exist_ok=True)
    rep.to_csv(os.path.join(OUT, f"{sp}_representativas.tsv"), sep="\t", index=False)
    nuc = rep[rep.cromosoma.notna()].sort_values(["cromosoma", "ini"])
    with open(os.path.join(OUT, f"{sp}.pep"), "w") as fp, open(os.path.join(OUT, f"{sp}.cds"), "w") as fc:
        for _, r in nuc.iterrows():
            fp.write(f">{r.gen}\n{prot[r.proteina].rstrip('*')}\n")
            fc.write(f">{r.gen}\n{cds[r.proteina][0][cds[r.proteina][1] - 1:]}\n")
    nuc[["cromosoma", "ini", "fin", "gen", "hebra"]].assign(ini=lambda x: x.ini - 1, puntaje=0)[
        ["cromosoma", "ini", "fin", "gen", "puntaje", "hebra"]].to_csv(os.path.join(OUT, f"{sp}.bed"), sep="\t", header=False, index=False)
    return nuc


def resumen(sp, rep, nuc, crom):
    comp5 = rep.empieza_ATG.notna()
    return dict(genes_con_proteina=int(len(rep)), en_cromosomas=int(len(nuc)), cromosomas=len({v[0] for v in crom.values()}),
                traduccion_ok=int(rep.traduccion_ok.sum()), traduccion_falla=int((~rep.traduccion_ok).sum()),
                inicio_alternativo=int(rep.inicio_alternativo.sum()), con_paradas_internas=int((rep.paradas_internas > 0).sum()),
                con_selenocisteina=int((rep.selenocisteinas > 0).sum()), completas_5p=int(comp5.sum()),
                completas_5p_con_ATG=int(rep.loc[comp5, "empieza_ATG"].astype(bool).sum()),
                parciales=int(rep.parcial.sum()), genes_con_varias_isoformas=int((rep.n_isoformas > 1).sum()))


if __name__ == "__main__":
    t01 = pd.read_csv(os.path.join(A, "tablas", "T01_censo_maestro_189.csv"), dtype=str, keep_default_na=False)
    res = {}
    for sp, cfg in ESPECIES.items():
        rep, prot, cds, crom = tabla_egl(cfg, t01) if sp == "Egl" else tabla_ncbi(sp, cfg)
        rep = controles(rep, prot, cds)
        nuc = escribir(sp, rep, prot, cds)
        res[sp] = resumen(sp, rep, nuc, crom)
        if sp == "Egl":
            nac = rep.set_index("gen").loc[t01.locus]
            res[sp]["NAC_189"] = dict(
                encontrados=int(len(nac)), proteina_igual_T01=int((nac.proteina.values == t01.proteina.values).sum()),
                largo_aa_igual_T01=int((nac.largo_aa.astype(str).values == t01.largo_aa.values).sum()),
                traduccion_ok=int(nac.traduccion_ok.sum()), empiezan_ATG=int(nac.empieza_ATG.eq(True).sum()),
                parciales=int(nac.parcial.sum()), en_cromosomas=int(nac.cromosoma.notna().sum()))
        print(sp, json.dumps(res[sp], ensure_ascii=False))
    json.dump(res, open(os.path.join(OUT, "resumen_isoformas.json"), "w"), ensure_ascii=False, indent=1)
