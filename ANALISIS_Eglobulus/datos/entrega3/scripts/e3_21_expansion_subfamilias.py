#!/usr/bin/env python
"""Entrega 3, complemento: ¿qué subfamilias NAC se expandieron y cómo, y se corresponden con pared
secundaria / xilema o con respuesta a sequía?

No recalcula nada: cruza por subfamilia (T01.subfamilia_final) lo que ya está en
  T01 censo, T41 homólogos de A. thaliana, T45 modo de duplicación de MCScanX, T47 Ka/Ks de pares,
  T44 compareCluster GO por subfamilia, T36 expresión diferencial en hoja y
  Hussey_FigS1_subfamilias_189.csv (subfamilias de los 189 EgrNAC de E. grandis, Hussey et al. 2015).

Definiciones
  pares locales = tipo_par tandem o proximal de MCScanX (T47) dentro de la misma subfamilia, no saturados (Ks <= 2)
  pares segmentales = tipo_par WGD/segmental (cualquier subfamilia de la pareja)
  Ks < 1: más joven que el pico de la duplicación genómica completa (Ks 1,475; resumen_kaks_Egl.json)
  GO pared/xilema: términos BP significativos de la subfamilia (T44, Count >= 2) con cell wall|xylem|tracheary|lignin
  GO estrés: ... con stress|water deprivation|abscisic|osmotic|hypoxia|wounding|abiotic
  al alza en sequía: Leaf_NA_DS, FDR < 0,05 y log2FC >= 1 (T36); solo cuentan los loci cuantificados en T36
  expandidas: IVa, IVc y VII (las de más genes; Hussey et al. 2015 señalan IVa, IVc, Vb y VII en E. grandis)

Salidas
  tablas/T62_expansion_subfamilias_NAC.csv
  datos/entrega3/expansion/resumen_expansion.json  (pruebas exactas de Fisher y de Mann-Whitney)
"""
import json, os
import numpy as np, pandas as pd
from scipy.stats import fisher_exact, mannwhitneyu

A = os.path.dirname(os.path.abspath(__file__))
while os.path.basename(A) != "ANALISIS_Eglobulus":
    A = os.path.dirname(A)
T = os.path.join(A, "tablas")
RZ = os.path.dirname(A)
NE = "no evaluada (dominio parcial)"
EXP = ["IVa", "IVc", "VII"]
SCW = r"cell wall|xylem|tracheary|lignin"
STR = r"stress|water deprivation|abscisic|osmotic|hypoxia|wounding|abiotic"
verdad = lambda s: s.astype(str).str.lower().isin(["true", "1"])

t01 = pd.read_csv(os.path.join(T, "T01_censo_maestro_189.csv"))
t41 = pd.read_csv(os.path.join(T, "T41_homologos_Ath_NAC189.csv"))
t45 = pd.read_csv(os.path.join(T, "T45_modos_duplicacion_NAC189.csv"))
t47 = pd.read_csv(os.path.join(T, "T47_KaKs_pares_NAC.csv"))
t36 = pd.read_csv(os.path.join(T, "T36_DE_NAC_globulus.csv"))
t44 = pd.read_csv(os.path.join(T, "T44_compareCluster_GO.csv"))
hz = pd.read_csv(os.path.join(RZ, "Hussey_FigS1_subfamilias_189.csv"))
ks_pico = json.load(open(os.path.join(A, "datos", "entrega3", "validacion", "resumen_controles.json")))["datacion"]["Ks_pico"]
assert len(t01) == 189 and len(hz) == 189

g = t01[["simbolo", "locus", "subfamilia_final"]].rename(columns={"subfamilia_final": "sf"})
g = g.merge(t45[["locus", "modo_duplicacion"]], on="locus", validate="1:1")
g = g.merge(t41[["locus", "AGI", "nombre_Ath", "pasa_filtros"]], on="locus", validate="1:1")
g["homologo"] = verdad(g.pasa_filtros)
de = t36[(t36.FDR < 0.05) & (t36.logFC >= 1)]
g["cuantificado"] = g.locus.isin(set(t36.locus))
g["al_alza_sequia"] = g.locus.isin(set(de[de.contraste == "Leaf_NA_DS"].locus))
g["al_alza_frio"] = g.locus.isin(set(de[de.contraste.isin(["Leaf_NA_CABF", "Leaf_NA_CAAF"])].locus))
assert set(t36.contraste) == {"Leaf_NA_DS", "Leaf_NA_CABF", "Leaf_NA_CAAF", "Leaf_CABF_CAAF"}

sat = verdad(t47.Ks_saturado)
ok = t47[~sat & t47.Ks_YN00.notna()]
loc = ok[ok.tipo_par.isin(["tandem", "proximal"]) & verdad(ok.misma_subfamilia)]
seg = t47[t47.tipo_par == "WGD/segmental"]
seg_ok = seg[~verdad(seg.Ks_saturado) & seg.Ks_YN00.notna()]
go = t44[(t44.nivel == "subfamilia") & verdad(t44.significativo) & (t44.Count >= 2)]

filas = []
for sf, x in g.groupby("sf"):
    a = loc[loc.subfamilia1 == sf]
    b = seg_ok[(seg_ok.subfamilia1 == sf) | (seg_ok.subfamilia2 == sf)]
    gs = go[go.grupo == sf]
    h = x[x.homologo]
    top = h.nombre_Ath.fillna(h.AGI).value_counts()
    filas.append(dict(
        subfamilia=sf, n_Egl=len(x), n_Egr_Hussey=int((hz.subfamilia_Hussey == sf).sum()),
        Egl_con_homologo_Ath=len(h), AGI_distintos=h.AGI.nunique(),
        Egl_por_AGI=round(len(h) / h.AGI.nunique(), 1) if len(h) else np.nan,
        homologo_principal=f"{top.index[0]} ({top.iloc[0]})" if len(top) else "",
        tandem_o_proximal=int(x.modo_duplicacion.isin(["tandem", "proximal"]).sum()),
        segmental=int((x.modo_duplicacion == "WGD/segmental").sum()),
        dispersos=int(x.modo_duplicacion.isin(["dispersed", "singleton"]).sum()),
        pares_locales=len(a), Ks_local_mediana=round(a.Ks_YN00.median(), 3) if len(a) else np.nan,
        pares_locales_Ks_menor_1=int((a.Ks_YN00 < 1).sum()),
        pares_segmentales_datables=len(b), Ks_segmental_mediana=round(b.Ks_YN00.median(), 3) if len(b) else np.nan,
        GO_pared_xilema=int(gs.termino.str.contains(SCW, case=False).sum()),
        GO_estres=int(gs.termino.str.contains(STR, case=False).sum()),
        cuantificados=int(x.cuantificado.sum()), al_alza_sequia=int(x.al_alza_sequia.sum()),
        al_alza_frio=int(x.al_alza_frio.sum())))
tab = pd.DataFrame(filas).sort_values(["n_Egl", "subfamilia"], ascending=[False, True])
tab.to_csv(os.path.join(T, "T62_expansion_subfamilias_NAC.csv"), index=False)


def fisher(m1, m2):
    t = [[int((m1 & m2).sum()), int((m1 & ~m2).sum())], [int((~m1 & m2).sum()), int((~m1 & ~m2).sum())]]
    r = fisher_exact(t)
    return dict(tabla=t, OR=round(float(r.statistic), 2), p=float(r.pvalue))


gg = g[g.sf != NE]
e = gg.sf.isin(EXP)
q = gg[gg.cuantificado]
res = dict(
    ks_pico_WGD=ks_pico, expandidas=EXP,
    n_expandidas=int(e.sum()), n_resto=int((~e).sum()),
    tandem_o_proximal=fisher(e, gg.modo_duplicacion.isin(["tandem", "proximal"])),
    segmental=fisher(e, gg.modo_duplicacion == "WGD/segmental"),
    Ks_locales_expandidas_vs_segmentales=dict(
        mediana_locales=round(float(loc[loc.subfamilia1.isin(EXP)].Ks_YN00.median()), 3),
        n_locales=int(loc.subfamilia1.isin(EXP).sum()),
        mediana_segmentales=round(float(seg_ok.Ks_YN00.median()), 3), n_segmentales=len(seg_ok),
        segmentales_total=len(seg), segmentales_saturados=int(verdad(seg.Ks_saturado).sum()),
        p_MannWhitney=float(mannwhitneyu(loc[loc.subfamilia1.isin(EXP)].Ks_YN00, seg_ok.Ks_YN00).pvalue)),
    cuantificados_expandidas=int(q.sf.isin(EXP).sum()),
    sequia_expandidas=fisher(q.sf.isin(EXP), q.al_alza_sequia),
    sequia_Vb=fisher(q.sf == "Vb", q.al_alza_sequia),
    al_alza_sequia_total=int(q.al_alza_sequia.sum()), cuantificados_total=len(q))
os.makedirs(os.path.join(A, "datos", "entrega3", "expansion"), exist_ok=True)
json.dump(res, open(os.path.join(A, "datos", "entrega3", "expansion", "resumen_expansion.json"), "w"), ensure_ascii=False, indent=1)
print(tab.to_string(index=False))
print(json.dumps(res, ensure_ascii=False))
