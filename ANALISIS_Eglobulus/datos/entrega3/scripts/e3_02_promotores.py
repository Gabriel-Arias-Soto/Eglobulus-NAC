"""Entrega 3, paso 2.1: promotores de 2 kb aguas arriba del ATG (E. globulus, GCA_046226855.1).

Definición (artículos 1 y 3): 2000 pb inmediatamente aguas arriba del primer nucleótido del
codón de inicio de la CDS representativa, según la hebra. Se extraen para todos los genes con
CDS completa en 5' y ATG en los 11 cromosomas: los 189 NAC son el conjunto de interés y el
resto, el fondo genómico.

Controles por secuencia: ATG inmediatamente después de la ventana, longitud (truncado en el
extremo del cromosoma), % de N, GC y pb de la ventana que caen dentro de otros genes.

Salidas:
  datos/entrega3/promotores/promotores_2kb_todos.fa      todos los genes (cabecera = locus)
  datos/entrega3/promotores/promotores_2kb_todos.tsv     controles por gen
  datos/entrega3/promotores/NAC189_promotores_2kb.fa     189 NAC (cabecera = símbolo EglNAC)
  datos/entrega3/cis_plantcare/envios/lote_{1..4}.fa         lotes de < 100 000 bytes para PlantCARE
  tablas/T40_promotores_NAC_2kb.csv                      controles de los 189
"""
import bisect, gzip, json, os
import pandas as pd

BASE = (os.environ.get("NAC_ROOT", ".") + "")
A = os.path.join(BASE, "ANALISIS_Eglobulus")
E3 = os.path.join(A, "datos", "entrega3")
GENOMA = os.path.join(BASE, "rnaseq_udec", "ref", "GCA_046226855.1_ASM4622685v1_genomic.fna.gz")
LARGO = 2000
TAM_MAX_LOTE = 100_000
N_LOTES = 4

COMP = str.maketrans("ACGTNacgtnRYKMSWBDHVrykmswbdhv", "TGCANtgcanYRMKSWVHDBymrkswvhdb")


def revcomp(s):
    return s.translate(COMP)[::-1]


def leer_cromosomas(p, seqids):
    seqs, h, buf = {}, None, []
    with gzip.open(p, "rt") as fh:
        for l in fh:
            if l.startswith(">"):
                if h in seqids:
                    seqs[h] = "".join(buf).upper()
                h, buf = l[1:].split()[0], []
            elif h in seqids:
                buf.append(l.strip())
    if h in seqids:
        seqs[h] = "".join(buf).upper()
    return seqs


def solape_bp(ini, fin, intervalos, inicios, excluir):
    """pb de [ini, fin] (1-based, inclusivo) cubiertos por genes distintos de `excluir`."""
    cub = []
    j = bisect.bisect_right(inicios, fin)
    for k in range(max(0, j - 400), j):
        a, b, g = intervalos[k]
        if g != excluir and b >= ini and a <= fin:
            cub.append((max(a, ini), min(b, fin)))
    cub.sort()
    total, cur_a, cur_b = 0, None, None
    for a, b in cub:
        if cur_b is None or a > cur_b + 1:
            if cur_b is not None:
                total += cur_b - cur_a + 1
            cur_a, cur_b = a, b
        else:
            cur_b = max(cur_b, b)
    if cur_b is not None:
        total += cur_b - cur_a + 1
    return total


def promotor(cromo, hebra, atg):
    """atg: coordenada 1-based del primer nucleótido del codón de inicio en la hebra del gen."""
    L = len(cromo)
    if hebra == "+":
        ini, fin = max(1, atg - LARGO), atg - 1
        s = cromo[ini - 1:fin]
        codon = cromo[atg - 1:atg + 2]
    else:
        ini, fin = atg + 1, min(L, atg + LARGO)
        s = revcomp(cromo[ini - 1:fin])
        codon = revcomp(cromo[atg - 3:atg])
    return s, ini, fin, codon


def main():
    rep = pd.read_csv(os.path.join(E3, "isoformas", "Egl_representativas.tsv"), sep="\t")
    rep = rep[rep.cromosoma.notna()].copy()
    t01 = pd.read_csv(os.path.join(A, "tablas", "T01_censo_maestro_189.csv"), dtype=str, keep_default_na=False)
    todos = pd.read_csv(os.path.join(A, "datos", "genes_todos.tsv"), sep="\t", header=None,
                        names=["seqid", "ini", "fin", "hebra", "gen", "biotipo"])
    seqids = set(rep.seqid)
    cromos = leer_cromosomas(GENOMA, seqids)
    assert set(cromos) == seqids, (len(cromos), len(seqids))
    interv = {}
    for sid, g in todos[todos.seqid.isin(seqids)].groupby("seqid"):
        g = g.sort_values("ini")
        interv[sid] = (list(zip(g.ini, g.fin, g.gen)), list(g.ini))

    filas, seqs = [], {}
    for r in rep.itertuples():
        comp5 = isinstance(r.parcial_tipo, float) or "5'" not in str(r.parcial_tipo)
        if not comp5:
            continue
        s, ini, fin, codon = promotor(cromos[r.seqid], r.hebra, int(r.cds_inicio))
        iv, ini_s = interv[r.seqid]
        n_n = s.count("N")
        gc = (s.count("G") + s.count("C")) / max(1, len(s) - n_n)
        filas.append(dict(locus=r.gen, cromosoma=r.cromosoma, seqid=r.seqid, hebra=r.hebra, atg=int(r.cds_inicio),
                          ventana_ini=ini, ventana_fin=fin, largo=len(s), codon_siguiente=codon,
                          atg_verificado=codon == "ATG", pct_N=round(100 * n_n / max(1, len(s)), 2), GC=round(gc, 4),
                          solape_otros_genes_bp=solape_bp(ini, fin, iv, ini_s, r.gen)))
        seqs[r.gen] = s
    tab = pd.DataFrame(filas)
    nac = set(t01.locus)
    tab["es_NAC"] = tab.locus.isin(nac)
    os.makedirs(os.path.join(E3, "promotores"), exist_ok=True)
    tab.to_csv(os.path.join(E3, "promotores", "promotores_2kb_todos.tsv"), sep="\t", index=False)
    with open(os.path.join(E3, "promotores", "promotores_2kb_todos.fa"), "w") as fh:
        for g in tab.locus:
            fh.write(f">{g}\n{seqs[g]}\n")

    # 189 NAC en el orden del censo
    sim = dict(zip(t01.locus, t01.simbolo))
    t40 = t01[["simbolo", "locus", "subfamilia_final", "inicio", "fin"]].rename(
        columns={"inicio": "gen_inicio", "fin": "gen_fin"}).merge(tab, on="locus", how="left", validate="1:1")
    assert t40.largo.notna().all(), "NAC sin promotor"
    t40["atg_dentro_del_gen"] = [int(a) >= int(i) and int(a) <= int(f) for a, i, f in zip(t40.atg, t40.gen_inicio, t40.gen_fin)]
    t40.drop(columns=["es_NAC"]).to_csv(os.path.join(A, "tablas", "T40_promotores_NAC_2kb.csv"), index=False)
    with open(os.path.join(E3, "promotores", "NAC189_promotores_2kb.fa"), "w") as fh:
        for g in t40.locus:
            fh.write(f">{sim[g]}\n{seqs[g]}\n")

    # lotes para PlantCARE: secuencia en una sola línea, < 100 000 bytes por fichero
    orden = sorted(t40.locus, key=lambda g: int(sim[g].replace("EglNAC", "")))
    tam = -(-len(orden) // N_LOTES)
    os.makedirs(os.path.join(E3, "cis_plantcare", "envios"), exist_ok=True)
    lotes = []
    for i in range(N_LOTES):
        parte = orden[i * tam:(i + 1) * tam]
        txt = "".join(f">{sim[g]}\n{seqs[g]}\n" for g in parte)
        b = txt.encode("ascii")
        assert len(b) < TAM_MAX_LOTE, (i, len(b))
        p = os.path.join(E3, "cis_plantcare", "envios", f"lote_{i + 1}.fa")
        open(p, "wb").write(b)
        lotes.append(dict(lote=i + 1, n=len(parte), bytes=len(b), primero=sim[parte[0]], ultimo=sim[parte[-1]]))

    fondo = tab[~tab.es_NAC]
    res = dict(genes_en_cromosomas=int(len(rep)), promotores=int(len(tab)), ATG_verificado=int(tab.atg_verificado.sum()),
               truncados=int((tab.largo < LARGO).sum()), con_N=int((tab.pct_N > 0).sum()),
               NAC=dict(n=int(len(t40)), ATG_verificado=int(t40.atg_verificado.sum()), largo_2000=int((t40.largo == LARGO).sum()),
                        con_N=int((t40.pct_N > 0).sum()), max_pct_N=float(t40.pct_N.max()),
                        atg_dentro_del_gen=int(t40.atg_dentro_del_gen.sum()),
                        solapan_otro_gen=int((t40.solape_otros_genes_bp > 0).sum()),
                        solape_mediana_bp=float(t40.solape_otros_genes_bp.median()),
                        GC_mediana=float(t40.GC.median()), GC_rango=[float(t40.GC.min()), float(t40.GC.max())]),
               fondo=dict(n=int(len(fondo)), GC_mediana=float(fondo.GC.median()),
                          solapan_otro_gen=int((fondo.solape_otros_genes_bp > 0).sum())),
               lotes=lotes)
    json.dump(res, open(os.path.join(E3, "promotores", "resumen_promotores.json"), "w"), ensure_ascii=False, indent=1)
    print(json.dumps(res, ensure_ascii=False))


if __name__ == "__main__":
    main()
