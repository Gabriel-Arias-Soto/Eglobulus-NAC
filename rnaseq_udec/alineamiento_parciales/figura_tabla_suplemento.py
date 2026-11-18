# -*- coding: utf-8 -*-
"""Figura y tabla suplementarias (en inglés) del chequeo de RNA-seq de los 6 NAC de dominio parcial.

  - F23_ARNseq_dominio_parcial.png (ANALISIS_Eglobulus/figuras): EglNAC87 (paneles a-d) y EglNAC249 (e-h), con las
    mismas funciones de dibujo de las figuras por gen (figuras_parciales_base.py) y los textos traducidos.
  - T64_ARNseq_dominio_parcial.csv (ANALISIS_Eglobulus/tablas): veredictos_parciales.csv con encabezados y veredictos en
    inglés y la posición de cada gen (T01).
  - resumen_suplemento.json (resultados/): las cifras que citan Métodos S1 y el manuscrito, leídas de los resultados.
No se teclea ninguna cifra; para si falta un texto por traducir.
"""
import json, os, re, subprocess
import pandas as pd

AQUI = os.path.dirname(os.path.abspath(__file__))
NS = {"__name__": "base"}
exec(open(os.path.join(AQUI, "figuras_parciales_base.py"), encoding="utf-8").read(), NS)
plt, mpl, V, J, RES, t01, R = NS["plt"], NS["mpl"], NS["V"], NS["J"], NS["RES"], NS["t01"], NS["R"]
B = NS["B"]
FIGD, TABD = f"{B}/ANALISIS_Eglobulus/figuras", f"{B}/ANALISIS_Eglobulus/tablas"

VER_EN = {"modelo incompleto": "incomplete gene model", "compatible con truncado": "consistent with a truncated domain",
          "sin datos suficientes": "insufficient data", "control: modelo confirmado": "control: model confirmed"}
assert set(V.veredicto) <= set(VER_EN), set(V.veredicto)


def nueva_union(gen):
    jj = J[(J.gen == gen) & J.fiable & ~J.anotada & J.toca_modelo & (J.hebra_motivo == J.hebra_gen)]
    return jj.sort_values("lecturas", ascending=False).iloc[0]


def titulo_en(gen):
    v = V.loc[gen]
    if v.veredicto == "modelo incompleto":
        u = nueva_union(gen)
        return (f"{gen}: an unannotated {u.intron_pb}-bp intron ({u.lecturas} reads) completes the NAM domain "
                f"({v.cob_dominio_mejor_isoforma:.0%} of the profile)")
    if v.veredicto == "compatible con truncado":
        return f"{gen}: no read-supported isoform completes the NAM domain (at most {v.cob_dominio_mejor_isoforma:.0%} of the profile)"
    raise ValueError(gen)


NS["titulo"] = titulo_en
EXACTO = {"Modelo anotado": "Annotated model", "Isoforma dominante (ARN)": "Dominant isoform (RNA-seq)",
          "Dominio NAM en el genoma": "NAM domain in the genome", "sin isoforma: no hay uniones fiables": "no isoform: no reliable junctions",
          "ninguna unión fiable en la región": "no reliable junction in the region", "unión anotada": "annotated junction",
          "unión no anotada": "unannotated junction", "unión débil que completa un dominio": "weak junction completing a domain",
          "Cobertura\n(lecturas únicas)": "Coverage\n(unique reads)", "Proteína anotada": "Annotated protein",
          "Mejor isoforma (ARN)": "Best isoform (RNA-seq)", "Posición en el perfil NAM (PF02365)": "Position in the NAM profile (PF02365)"}


def tr(s):
    if s in EXACTO:
        return EXACTO[s]
    s2 = re.sub(r"\(marco (\d)\)", r"(frame \1)", s)
    s2 = re.sub(r"^(chr\d+) \(Mb\), gen en la hebra (.)$", r"\1 (Mb), gene on the \2 strand", s2)
    assert not re.search(r"[áéíóúñ]|\b(gen|hebra|lecturas|unión|ninguna|marco|isoforma)\b", s2), s2
    return s2


def en_ingles(fig, letras):
    for ax in fig.axes:
        ax.set_yticklabels([tr(t.get_text()) for t in ax.get_yticklabels()]) if ax.get_yticks().size and any(
            t.get_text() for t in ax.get_yticklabels()) else None
        ax.set_xlabel(tr(ax.get_xlabel())); ax.set_ylabel(tr(ax.get_ylabel()))
        lg = ax.get_legend()
        if lg:
            for t in lg.get_texts():
                t.set_text(tr(t.get_text()))
        for t in ax.texts:
            s = t.get_text()
            if s in ("a", "b", "c", "d"):
                t.set_text(letras["abcd".index(s)])
            else:
                t.set_text(tr(s))
    return fig


GENES = [("EglNAC87", "abcd"), ("EglNAC249", "efgh")]
os.makedirs("/tmp/f23", exist_ok=True)
pngs = []
for gen, letras in GENES:
    fig = en_ingles(NS["dibujar4"](gen), letras)
    f = f"/tmp/f23/{gen}.png"
    fig.savefig(f, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    pngs.append(f)
from PIL import Image
ims = [Image.open(f).convert("RGB") for f in pngs]
gap = 90
W = max(i.width for i in ims)
out = Image.new("RGB", (W, sum(i.height for i in ims) + gap), "white")
y = 0
for i in ims:
    out.paste(i, (0, y)); y += i.height + gap
F23 = f"{FIGD}/F23_ARNseq_dominio_parcial.png"
out.save(F23, dpi=(300, 300))

# --- tabla T64 -----------------------------------------------------------------------------------------------
v = V.reset_index()
t = pd.DataFrame({
    "gene": v.gen, "locus": v.locus, "chromosome": v.gen.map(t01.cromosoma), "start": v.gen.map(t01.inicio),
    "end": v.gen.map(t01.fin), "strand": v.hebra,
    "profile_coverage_annotated_protein": v.cob_dominio_anotada,
    "unique_reads_annotated_exons": v.unicas_exones, "multimapping_reads_annotated_exons": v.multiples_exones,
    "mean_coverage_annotated_exons": v.cobertura_media_exones,
    "reliable_junctions": v.uniones_fiables, "reliable_unannotated_junctions": v.uniones_fiables_no_anotadas,
    "isoforms_tested": v.isoformas_probadas,
    "profile_coverage_dominant_isoform": v.cob_dominio_isoforma_dominante,
    "profile_coverage_best_isoform": v.cob_dominio_mejor_isoforma,
    "verdict": v.veredicto.map(VER_EN)})
t.to_csv(f"{TABD}/T64_ARNseq_dominio_parcial.csv", index=False)

# --- cifras para el texto ------------------------------------------------------------------------------------
u87 = nueva_union("EglNAC87")
bib = RES["bibliotecas"]
hv = subprocess.run([(os.environ.get("CONDA_ENVS", os.path.expanduser("~/conda/envs")) + "/rnaseq/bin/hisat2"), "--version"], capture_output=True,
                    text=True, env={**os.environ, "PATH": (os.environ.get("CONDA_ENVS", os.path.expanduser("~/conda/envs")) + "/rnaseq/bin") + ":" + os.environ["PATH"]}).stdout
hisat_v = re.search(r"version (\d+\.\d+\.\d+)", hv).group(1)
sin = [g for g in V.index if V.loc[g, "veredicto"] == "sin datos suficientes"]
RS = dict(
    n_bibliotecas=len(bib), bibliotecas=sorted(bib), lecturas_M=round(sum(b["lecturas"] for b in bib.values()) / 1e6),
    hisat2=hisat_v, samtools=RES["versiones"]["samtools"].split()[-1], hmmer=re.search(r"HMMER (\S+)", RES["versiones"]["hmmsearch"]).group(1),
    min_lecturas_union=RES["parametros"]["MIN_LECT_UNION"], min_bibliotecas_union=RES["parametros"]["MIN_BIBL_UNION"],
    max_intron_kb=RES["parametros"]["MAX_INTRON"] // 1000, cobertura_completa_pct=round(RES["parametros"]["COBERTURA_COMPLETA"] * 100),
    EglNAC87=dict(intron_pb=int(u87.intron_pb), lecturas=int(u87.lecturas), bibliotecas=int(u87.bibliotecas),
                  cob_anotada_pct=round(V.loc["EglNAC87", "cob_dominio_anotada"] * 100),
                  cob_isoforma_pct=round(V.loc["EglNAC87", "cob_dominio_mejor_isoforma"] * 100)),
    EglNAC249=dict(isoformas=int(V.loc["EglNAC249", "isoformas_probadas"]), cob_max_pct=round(V.loc["EglNAC249", "cob_dominio_mejor_isoforma"] * 100)),
    EglNAC66=dict(uniones=[int(x) for x in J[(J.gen == "EglNAC66") & J.anotada & J.fiable].sort_values("donante_fin").lecturas],
                  cob_pct=round(V.loc["EglNAC66", "cob_dominio_mejor_isoforma"] * 100)),
    sin_datos=sin, veredictos={g: VER_EN[V.loc[g, "veredicto"]] for g in V.index})
assert RS["EglNAC87"]["cob_isoforma_pct"] >= RS["cobertura_completa_pct"] and RS["EglNAC249"]["cob_max_pct"] < RS["cobertura_completa_pct"]
assert len(sin) == 4 and V.loc["EglNAC66", "veredicto"].startswith("control")
json.dump(RS, open(f"{R}/resumen_suplemento.json", "w", encoding="utf-8"), ensure_ascii=False, indent=1)
print(out.size, len(t), json.dumps(RS, ensure_ascii=False))
