# -*- coding: utf-8 -*-
"""Tabla suplementaria (en inglés) de la revisión de los modelos génicos de los 189 NAC y cifras para Métodos S1.

  - ANALISIS_Eglobulus/tablas/T65_revision_modelos_NAC189.csv: revision_modelos_189.csv con encabezados y veredictos
    en inglés.
  - resultados/resumen_S1.json: las cifras que cita Métodos S1, leídas de los resultados y de las entradas.
No se teclea ninguna cifra; para si aparece un veredicto sin traducción.
"""
import json, os, re, subprocess
import pandas as pd

AQUI = os.path.dirname(os.path.abspath(__file__))
RAIZ = os.path.dirname(os.path.dirname(AQUI))
TAB = os.path.join(RAIZ, "ANALISIS_Eglobulus", "tablas")
V = pd.read_csv(os.path.join(AQUI, "resultados", "revision_modelos_189.csv"))
RES = json.load(open(os.path.join(AQUI, "resultados", "resumen_revision.json"), encoding="utf-8"))
assert len(V) == 189

VER_EN = {
    "modelo respaldado": "model supported",
    "modelo respaldado por RNA-seq (la proteína homóloga difiere)": "model supported by RNA-seq (homologous protein differs)",
    "isoformas alternativas (RNA-seq apoya la anotada y la de la proteína)": "alternative isoforms (RNA-seq supports both)",
    "intrón no anotado en la CDS": "unannotated intron in the CDS",
    "posible pseudogén o error de secuencia": "possible pseudogene or sequence error",
    "modelo incompleto según proteínas": "incomplete model (protein evidence)",
    "estructura distinta (proteína y RNA-seq)": "different structure (protein and RNA-seq)",
    "estructura distinta según proteínas, sin RNA-seq": "different structure (protein only; no RNA-seq)",
    "sitio de empalme desplazado": "shifted splice site",
    "diferencia con el ortólogo fuera del dominio": "difference from ortholog outside the domain",
    "diferencia con el ortólogo": "difference from ortholog",
    "evidencia parcial": "partial evidence",
    "sin evidencia suficiente": "insufficient evidence",
}
assert set(V.veredicto) <= set(VER_EN), set(V.veredicto) - set(VER_EN)
SUB_EN = {"no evaluada (dominio parcial)": "not assessed (partial domain)"}
t = pd.DataFrame({
    "symbol": V.gen, "locus": V.locus, "chromosome": V.cromosoma, "strand": V.hebra,
    "subfamily": V.subfamilia.replace(SUB_EN),
    "NAM_profile_coverage_annotated": V.cobertura_NAM_anotada,
    "annotated_CDS_introns": V.intrones_CDS_anotados,
    "best_protein": V.prot_mejor.fillna("").str.replace("Egr|", "E. grandis ", regex=False).str.replace("Ath|", "A. thaliana ", regex=False),
    "best_protein_identity": V.prot_identidad, "best_protein_coverage": V.prot_cobertura,
    "frameshifts_and_stops": V.eventos_marco_parada.fillna("").str.replace("marco@", "frameshift@", regex=False).str.replace("parada@", "stop@", regex=False),
    "events_in_NAM_domain": V.eventos_en_dominio,
    "NAM_profile_coverage_clean_protein": V.cobertura_NAM_prot_limpia,
    "introns_compared_annotated": V.intrones_comparados_anot, "introns_compared_protein": V.intrones_comparados_prot,
    "introns_shared": V.intrones_coinciden,
    "protein_only_introns_with_RNAseq": V.intrones_prot_con_RNA, "annotated_only_introns_with_RNAseq": V.intrones_anot_sin_prot_con_RNA,
    "unique_reads_in_exons": V.lecturas_exones, "reliable_junctions": V.uniones_fiables,
    "annotated_CDS_introns_confirmed_RNAseq": V.intrones_CDS_confirmados_RNA,
    "unannotated_junctions_in_CDS": V.uniones_nuevas_CDS, "dominant_unannotated_junctions": V.uniones_nuevas_dominantes,
    "verdict": V.veredicto.map(VER_EN)})
t.to_csv(os.path.join(TAB, "T65_revision_modelos_NAC189.csv"), index=False)

# --- cifras para Métodos S1 ---------------------------------------------------------------------------------
cnt = V.veredicto.value_counts()
c = lambda v: int(cnt.get(v, 0))
prot = [l for l in open(os.path.join(AQUI, "proteinas_evidencia.faa")) if l.startswith(">")]
genes = lambda v: sorted(V[V.veredicto == v].gen, key=lambda g: int(re.sub(r"\D", "", g)))
p = RES["parametros"]
sin = V[V.veredicto == "sin evidencia suficiente"]
S = dict(
    miniprot=RES["versiones"]["miniprot"].split("-")[0], n_egr=sum(l.startswith(">Egr|") for l in prot), n_ath=sum(l.startswith(">Ath|") for l in prot),
    n_bibliotecas=len(RES["libs"]), qcov_limpio_pct=round(p["QCOV_LIMPIO"] * 100), frac_dominante_pct=round(p["FRAC_DOMINANTE"] * 100),
    respaldado=c("modelo respaldado"), respaldado_rna=c("modelo respaldado por RNA-seq (la proteína homóloga difiere)"),
    isoformas=c("isoformas alternativas (RNA-seq apoya la anotada y la de la proteína)"),
    intron_nuevo=c("intrón no anotado en la CDS"), intron_nuevo_genes=genes("intrón no anotado en la CDS"),
    pseudo=c("posible pseudogén o error de secuencia"), pseudo_genes=genes("posible pseudogén o error de secuencia"),
    incompleto=c("modelo incompleto según proteínas"), incompleto_genes=genes("modelo incompleto según proteínas"),
    estructura_rna=c("estructura distinta (proteína y RNA-seq)"), estructura_rna_genes=genes("estructura distinta (proteína y RNA-seq)"),
    estructura_sin=c("estructura distinta según proteínas, sin RNA-seq"),
    otros=c("sitio de empalme desplazado") + c("diferencia con el ortólogo fuera del dominio") + c("diferencia con el ortólogo") + c("evidencia parcial"),
    sin=c("sin evidencia suficiente"), sin_VII=int((sin.subfamilia == "VII").sum()),
    sin_lecturas=int((sin.lecturas_exones < 50).sum()),
)
assert sum(v for k, v in S.items() if k in ("respaldado", "respaldado_rna", "isoformas", "intron_nuevo", "pseudo", "incompleto",
                                              "estructura_rna", "estructura_sin", "otros", "sin")) == 189
# el umbral MIN_LECT_DOMINANTE y DESPLAZ_MAX se leen del guion de la revisión (no están en el JSON)
src = open(os.path.join(AQUI, "revisar_modelos.py"), encoding="utf-8").read()
m = re.search(r"FRAC_DOMINANTE, MIN_LECT_DOMINANTE, DESPLAZ_MAX = ([\d.]+), (\d+), (\d+)", src)
S["min_lect_dominante"], S["desplaz_max"] = int(m.group(2)), int(m.group(3))
json.dump(S, open(os.path.join(AQUI, "resultados", "resumen_S1.json"), "w", encoding="utf-8"), ensure_ascii=False, indent=1)
print(json.dumps(S, ensure_ascii=False))
