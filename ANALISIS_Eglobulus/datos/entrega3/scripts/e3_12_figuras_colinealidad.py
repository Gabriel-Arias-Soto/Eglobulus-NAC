#!/usr/bin/env python
# e3_12_figuras_colinealidad.py -- figuras de colinealidad y Ka/Ks de la tercera entrega.
#
#   F15_circos_colinealidad_Egl   grafico circular de E. globulus: 11 cromosomas, densidad
#                                 genica, los 189 NAC por modo de duplicacion (MCScanX), los
#                                 grupos de NAC a <= 200 kb, los bloques colineales del genoma
#                                 en gris y los pares segmentales NAC en rojo.
#   F16_cariotipo_Egr_Egl_Ptr     cariotipo JCVI E. grandis - E. globulus - P. trichocarpa con
#                                 los bloques colineales de fondo y cada NAC unido a su mejor
#                                 socio colineal (T48).
#   F17_KaKs_datacion             Ks de las anclas del genoma (pico y datacion) y Ks y Ka/Ks
#                                 de los 615 pares NAC por clase.
#
# Entradas: tablas T45, T47, T48; datos/entrega3/colinealidad/{Egl_pares_colineales.tsv,
#   Ks_anclas_genoma.tsv, resumen_kaks_Egl.json, jcvi/Egl.{Egr,Ptr}.lifted.anchors};
#   datos/entrega3/isoformas/{Egl,Egr,Ptr}.bed; referencia/.../sequence_report.json(l).
# Salidas: figuras/F15-F17 (.pdf y .png), tablas/_cifras/e3_colinealidad.csv y, para JCVI,
#   datos/entrega3/colinealidad/jcvi/figura/ (bloques filtrados, simple, layout, seqids).
# Ejecutar (con el python del entorno nac3, desde cualquier carpeta):
#   nac3/bin/python datos/entrega3/scripts/e3_12_figuras_colinealidad.py
#
# Reglas de 00_estilo.R que se repiten aqui porque matplotlib no pasa por theme_nac():
# sin titulos dentro de la imagen, 174 mm de ancho, 600 ppp, Arial o Liberation Sans de
# 8 pt como minimo, y cifras calculadas (nunca escritas a mano) volcadas a _cifras.

import glob
import json
import logging
import os
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from PIL import Image
from scipy.stats import gaussian_kde
# pyCirclize y JCVI se importan aqui, antes de fijar rcParams: pyCirclize pone al
# importarse savefig.bbox = "tight" y savefig.pad_inches = 0.5, que anaden 1 in a cada
# figura y rompen el ancho de 174 mm; JCVI deja el registro en DEBUG.
from pycirclize import Circos
from jcvi.graphics.karyotype import Karyotype
for _n in [n for n in logging.root.manager.loggerDict if n.split(".")[0] in ("jcvi", "fontTools")] + ["jcvi", "fontTools"]:
    logging.getLogger(_n).setLevel(logging.WARNING)
logging.getLogger().setLevel(logging.WARNING)

# --- rutas -------------------------------------------------------------------
AQUI = Path(__file__).resolve()
ANALISIS = AQUI.parents[3]
RAIZ = AQUI.parents[4]
TABLAS = ANALISIS / "tablas"
FIGURAS = ANALISIS / "figuras"
CIFRAS = TABLAS / "_cifras"
COL = ANALISIS / "datos" / "entrega3" / "colinealidad"
ISO = ANALISIS / "datos" / "entrega3" / "isoformas"
REF = ANALISIS / "referencia"
JF = COL / "jcvi" / "figura"
JF.mkdir(parents=True, exist_ok=True)
assert (RAIZ / "scripts" / "00_estilo.R").is_file(), RAIZ

# --- contrato de estilo (00_estilo.R) ----------------------------------------
_estilo = (RAIZ / "scripts" / "00_estilo.R").read_text(encoding="utf-8")


def cte(nombre):
    # valor entre comillas ("#B2182B", "grey88") o expresion sin comillas (174 * MM, 600)
    m = re.search(r'^%s\s*<-\s*(?:"([^"]*)"|([^#\n]+?))\s*(?:#.*)?$' % nombre, _estilo, re.M)
    assert m, "no encuentro %s en 00_estilo.R" % nombre
    v = (m.group(1) if m.group(1) is not None else m.group(2)).strip()
    mm = re.fullmatch(r"([\d.]+)\s*\*\s*MM", v)
    return str(float(mm.group(1)) / 25.4) if mm else v


def es_verdad(s):
    # columnas logicas que pueden llegar como bool o como texto "True"/"TRUE"
    return s.astype(str).str.strip().str.lower().isin(["true", "1"])


def color_r(v):
    # 00_estilo.R usa nombres de R ("grey88"); R los resuelve como round(n * 255 / 100).
    m = re.fullmatch(r"gr[ae]y(\d+)", v)
    if m:
        n = int(round(int(m.group(1)) * 255 / 100))
        return "#%02X%02X%02X" % (n, n, n)
    assert v.startswith("#"), v
    return v


ANCHO = float(cte("ANCHO_COLUMNA"))
ALTO_MAX = float(cte("ALTO_MAXIMO"))
DPI = int(float(cte("DPI")))
COL_DESTACADO = color_r(cte("COL_DESTACADO"))
COL_NEUTRO = color_r(cte("COL_NEUTRO"))
COL_IDEOGRAMA = color_r(cte("COL_IDEOGRAMA"))
COL_NAC = color_r(cte("COL_NAC"))
COL_BLOQUE = color_r(cte("COL_BLOQUE"))
MM = 1 / 25.4

# Colores propios de estas figuras. Rojo = par NAC dentro de un bloque colineal (segmental
# en F15, ortologo colineal en F16), como en Zhang et al. 2026 (Fig. 4); los modos locales
# van en azules (el oscuro es COL_BLOQUE, el de los bloques en tandem de F01 y F09).
COL_PROXIMAL = "#67A9CF"
COL_FONDO_BLOQUE = "#D4D4D4"     # ribbons del genoma en F15
COL_GAINSBORO = "#DCDCDC"        # fondo que JCVI pinta en F16 (fijo en ShadeManager)
COL_TEXTO_SEC = "#525252"

# --- tipografia ----------------------------------------------------------------
# Trees pide Helvetica o Arial. Liberation Sans (metricamente identica a Arial) viaja en
# el paquete R fontLiberation; se registra desde los entornos que la tengan.
for _f in glob.glob(os.path.join(os.environ.get("CONDA_ENVS", os.path.expanduser("~/conda/envs")),
                                 "*/lib/R/library/fontLiberation/fonts/liberation-fonts/LiberationSans-*.ttf")):
    try:
        fm.fontManager.addfont(_f)
    except Exception:
        pass
_disponibles = {f.name for f in fm.fontManager.ttflist}
_preferidas = [f for f in ("Arial", "Helvetica", "Liberation Sans") if f in _disponibles]
assert _preferidas, "no hay Arial, Helvetica ni Liberation Sans"
plt.rcParams.update({
    "font.family": "sans-serif",
    # solo las instaladas, para que matplotlib no avise por cada familia que falta
    "font.sans-serif": _preferidas + ["DejaVu Sans"],
    "savefig.bbox": "standard", "savefig.pad_inches": 0.0,
    "font.size": 8, "axes.labelsize": 8, "xtick.labelsize": 8, "ytick.labelsize": 8,
    "legend.fontsize": 8, "axes.linewidth": 0.5, "xtick.major.width": 0.5,
    "ytick.major.width": 0.5, "xtick.major.size": 2.5, "ytick.major.size": 2.5,
    "xtick.major.pad": 2, "ytick.major.pad": 2, "axes.labelpad": 3,
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
    "axes.spines.top": False, "axes.spines.right": False,
})
FUENTE = fm.FontProperties(fname=fm.findfont(fm.FontProperties(family=["sans-serif"]))).get_name()
assert FUENTE in ("Arial", "Helvetica", "Liberation Sans"), FUENTE
PT_MIN = 8
PT_LETRA = 11      # plot.tag de theme_nac(): negrita, base + 3

CIF = []


def cifra(clave, valor, fuente, nota="NA"):
    if isinstance(valor, (np.integer,)):
        valor = int(valor)
    if isinstance(valor, (float, np.floating)):
        valor = float(round(valor, 4))
    CIF.append(dict(clave=clave, valor=valor, fuente=fuente, nota=nota))


def guardar(fig, nombre):
    rutas = []
    for ext in ("pdf", "png"):
        p = FIGURAS / f"{nombre}.{ext}"
        fig.savefig(p, dpi=DPI)
        rutas.append(p)
    w, h = fig.get_size_inches()
    # El PNG debe medir exactamente lo pedido (sin recorte ni margen anadido) y caber en la caja.
    px = Image.open(rutas[1]).size
    assert abs(px[0] - w * DPI) <= 1 and abs(px[1] - h * DPI) <= 1, f"{nombre}: PNG {px} != {w * DPI:.0f} x {h * DPI:.0f}"
    assert w <= ANCHO + 1e-6 and h <= ALTO_MAX + 1e-6, f"{nombre}: {w / MM:.1f} x {h / MM:.1f} mm fuera de la caja"
    crudo = rutas[0].read_bytes()
    fuentes = sorted({f.decode().split("+")[-1] for f in re.findall(rb"/BaseFont\s*/([^\s/\]>]+)", crudo)})
    malas = [f for f in fuentes if not re.match(r"(Arial|Helvetica|LiberationSans)", f)]
    assert not malas, f"{nombre}: fuentes inesperadas en el PDF {malas}"
    print(f"  {nombre}: {w / MM:.0f} x {h / MM:.0f} mm, {DPI} ppp, fuentes {fuentes}")
    plt.close(fig)


def texto_minimo(fig):
    tam = [t.get_fontsize() for t in fig.findobj(matplotlib.text.Text) if t.get_text().strip()]
    assert min(tam) >= PT_MIN, f"texto de {min(tam)} pt"


# --- datos ---------------------------------------------------------------------
print("versiones:", {m: __import__(m).__version__ for m in ("numpy", "pandas", "matplotlib", "scipy",
                                                            "pycirclize", "jcvi")})
t45 = pd.read_csv(TABLAS / "T45_modos_duplicacion_NAC189.csv")
t47 = pd.read_csv(TABLAS / "T47_KaKs_pares_NAC.csv")
t48 = pd.read_csv(TABLAS / "T48_ortologos_colineales_NAC189.csv")
pc = pd.read_csv(COL / "Egl_pares_colineales.tsv", sep="\t")
an = pd.read_csv(COL / "Ks_anclas_genoma.tsv", sep="\t")
rk = json.load(open(COL / "resumen_kaks_Egl.json", encoding="utf-8"))
LAMBDA = rk["parametros"]["lambda_sitio_anio"]
KS_SAT = rk["parametros"]["Ks_saturado"]
NOMBRES_BED = ["crom", "ini", "fin", "gen", "puntaje", "hebra"]
bed = {sp: pd.read_csv(ISO / f"{sp}.bed", sep="\t", header=None, names=NOMBRES_BED)
       for sp in ("Egl", "Egr", "Ptr")}
for sp, b in bed.items():
    assert b.gen.is_unique, sp

# Longitudes de los cromosomas de E. globulus (misma fuente que e3_01_isoformas.py).
sr = json.load(open(REF / "Eglobulus_ASM4622685v1" / "sequence_report.json", encoding="utf-8"))["reports"]
LONG = {f"Egl{int(s['chr_name']):02d}": int(s["length"]) for s in sr
        if s.get("role") == "assembled-molecule" and s.get("chr_name", "").isdigit()}
LONG = dict(sorted(LONG.items()))
assert len(LONG) == 11, LONG
egl = bed["Egl"][bed["Egl"].crom.isin(LONG)].copy()
egl["mid"] = (egl.ini + egl.fin) / 2
assert (egl.fin <= egl.crom.map(LONG)).all()
POS = egl.set_index("gen")[["crom", "ini", "fin", "mid"]]
cifra("col_cromosomas_Egl", len(LONG), "sequence_report.json")
cifra("col_genes_Egl_en_cromosomas", len(egl), "isoformas/Egl.bed", "un modelo representativo por gen")
cifra("col_genes_Egl_fuera_de_cromosomas", int((~bed["Egl"].crom.isin(LONG)).sum()), "isoformas/Egl.bed")

# NAC: modo MCScanX por gen (T45) y posicion.
assert len(t45) == 189 and t45.locus.isin(POS.index).all()
nac = (t45.join(POS.add_prefix("bed_"), on="locus")
          .rename(columns={"bed_crom": "crom", "bed_mid": "mid", "bed_ini": "ini_bed", "bed_fin": "fin_bed"}))
assert nac.groupby("cromosoma").crom.nunique().eq(1).all() and nac.groupby("crom").cromosoma.nunique().eq(1).all()
MODOS = {"WGD/segmental": COL_DESTACADO, "tandem": COL_BLOQUE, "proximal": COL_PROXIMAL,
         "dispersed": COL_NAC, "singleton": COL_NAC}
assert set(nac.modo_duplicacion) <= set(MODOS), set(nac.modo_duplicacion)
for m, n in nac.modo_duplicacion.value_counts().items():
    cifra(f"col_modo_{m.replace('/', '_')}", n, "T45.modo_duplicacion", "genes NAC (MCScanX)")

# Pares segmentales NAC: T47 frente a las anclas NAC-NAC de los bloques del genoma.
seg = t47[t47.tipo_par == "WGD/segmental"]
par = lambda a, b: tuple(sorted((a, b)))
seg_pares = {par(a, b) for a, b in zip(seg.locus1, seg.locus2)}
nn = pc[es_verdad(pc.NAC1) & es_verdad(pc.NAC2)]
assert seg_pares == {par(a, b) for a, b in zip(nn.gen1, nn.gen2)}, "T47 y las anclas NAC-NAC no coinciden"
gen_seg = sorted(set(seg.locus1) | set(seg.locus2))
assert set(gen_seg) == set(nac.loc[nac.modo_duplicacion == "WGD/segmental", "locus"])
cifra("col_pares_segmentales_NAC", len(seg_pares), "T47.tipo_par")
cifra("col_genes_en_pares_segmentales_NAC", len(gen_seg), "T47")
cifra("col_pares_segmentales_misma_subfamilia", int(es_verdad(seg.misma_subfamilia).sum()), "T47.misma_subfamilia")
cifra("col_pares_segmentales_intercromosomicos", int((seg.cromosoma1 != seg.cromosoma2).sum()), "T47")
cifra("col_bloques_con_par_segmental_NAC", nn.bloque.nunique(), "Egl_pares_colineales.tsv")

# Bloques colineales del genoma (autocomparacion de E. globulus). Los nombres de
# cromosoma de MCScanX se sustituyen por los del BED, comprobando que la equivalencia
# es uno a uno.
pcx = (pc.rename(columns={"crom1": "crom1_mcscanx", "crom2": "crom2_mcscanx"})
         .join(POS.add_suffix("1"), on="gen1").join(POS.add_suffix("2"), on="gen2"))
assert pcx[["crom1", "crom2"]].notna().all().all()
for k in ("1", "2"):
    eq = pcx.groupby("crom%s_mcscanx" % k)["crom" + k].nunique()
    assert (eq == 1).all() and pcx.groupby("crom" + k)["crom%s_mcscanx" % k].nunique().eq(1).all()
bl = pcx.groupby("bloque").agg(c1=("crom1", "first"), s1=("ini1", "min"), e1=("fin1", "max"),
                               c2=("crom2", "first"), s2=("ini2", "min"), e2=("fin2", "max"),
                               n=("gen1", "size"))
assert (pcx.groupby("bloque").crom1.nunique() == 1).all() and (pcx.groupby("bloque").crom2.nunique() == 1).all()
cifra("col_bloques_intragenomicos", len(bl), "Egl_pares_colineales.tsv")
cifra("col_anclas_intragenomicas", len(pc), "Egl_pares_colineales.tsv")
cifra("col_bloques_intracromosomicos", int((bl.c1 == bl.c2).sum()), "Egl_pares_colineales.tsv")

# Grupos de NAC a <= 200 kb (criterio del articulo de referencia; T45.grupo_art3).
gr = nac[nac.n_genes_grupo_art3 > 1].groupby("grupo_art3").agg(crom=("crom", "first"), ini=("ini_bed", "min"),
                                                                fin=("fin_bed", "max"), n=("locus", "size"))
assert (nac[nac.n_genes_grupo_art3 > 1].groupby("grupo_art3").crom.nunique() == 1).all()
cifra("col_grupos_200kb", len(gr), "T45.grupo_art3", "grupos con 2 o mas NAC")
cifra("col_genes_en_grupos_200kb", int(gr.n.sum()), "T45.grupo_art3")

# ==================================================================================
# F15 -- grafico circular intraespecifico (pyCirclize)
# ==================================================================================
# Los grupos de <= 200 kb (criterio del art. 3) van en negro, fuera de la paleta de modos:
# COL_BLOQUE es el tandem (marca de F01, por adyacencia en orden genico) y aqui lo usa el
# modo tandem de MCScanX; con el mismo azul el anillo y las marcas se confundian.
COL_GRUPO = "#000000"
TICK_MB = 10               # Mb entre marcas del ideograma (se declara en la leyenda)
R_LINK = 52.0
R_DENS = (54.0, 62.0)
R_GRUPO = (63.0, 64.4)
R_NAC = (65.4, 70.6)
R_IDEO = (72.0, 76.6)
R_LAB = 79.4
ANCHO_GRUPO_MIN = 0.45     # Mb: los grupos de 200 kb se ensanchan hasta aqui para verse
VENTANA = 1.0              # Mb, densidad genica

circos = Circos({c: L / 1e6 for c, L in LONG.items()}, space=2.2)
vent = []
for c, L in LONG.items():
    Lm = L / 1e6
    bordes = np.arange(0, Lm, VENTANA)
    if Lm - bordes[-1] < 0.5 * VENTANA and len(bordes) > 1:
        bordes = bordes[:-1]              # ultima ventana corta: se une a la anterior
    bordes = np.append(bordes, Lm)
    n = np.histogram(egl.loc[egl.crom == c, "mid"] / 1e6, bins=bordes)[0]
    vent.append(pd.DataFrame(dict(crom=c, a=bordes[:-1], b=bordes[1:], n=n, d=n / np.diff(bordes))))
vent = pd.concat(vent)
DMAX = float(vent.d.max())
cifra("col_densidad_ventana_Mb", VENTANA, "e3_12")
cifra("col_densidad_max_genes_por_Mb", round(DMAX, 1), "isoformas/Egl.bed")
cifra("col_densidad_mediana_genes_por_Mb", round(float(vent.d.median()), 1), "isoformas/Egl.bed")

ORDEN_MODO = ["dispersed", "singleton", "proximal", "tandem", "WGD/segmental"]
for s in circos.sectors:
    t = s.add_track(R_IDEO)
    t.axis(fc=COL_IDEOGRAMA, ec=COL_NAC, lw=0.4)
    t.text(f"chr{int(s.name[3:])}", size=8, color="black")
    t.xticks_by_interval(TICK_MB, tick_length=0.9, outer=False, show_label=False,
                         line_kws=dict(lw=0.4, color=COL_NAC))
    v = vent[vent.crom == s.name]
    x = np.r_[0, (v.a + v.b).to_numpy() / 2, s.size]
    y = np.r_[v.d.iloc[0], v.d.to_numpy(), v.d.iloc[-1]]
    td = s.add_track(R_DENS)
    td.fill_between(x, y, vmin=0, vmax=DMAX, fc=COL_NEUTRO, ec="none", lw=0)
    tg = s.add_track(R_GRUPO)
    for _, g in gr[gr.crom == s.name].iterrows():
        a, b = g.ini / 1e6, g.fin / 1e6
        if b - a < ANCHO_GRUPO_MIN:
            m = (a + b) / 2
            a, b = m - ANCHO_GRUPO_MIN / 2, m + ANCHO_GRUPO_MIN / 2
        tg.rect(max(a, 0), min(b, s.size), fc=COL_GRUPO, ec="none", lw=0)
    tn = s.add_track(R_NAC)
    sub = nac[nac.crom == s.name].copy()
    sub["o"] = sub.modo_duplicacion.map(ORDEN_MODO.index)
    for _, g in sub.sort_values("o").iterrows():
        p = g.mid / 1e6
        tn.line([p, p], [0, 1], vmin=0, vmax=1, color=MODOS[g.modo_duplicacion], lw=0.6,
                solid_capstyle="butt")

for _, b in bl.iterrows():
    circos.link((b.c1, b.s1 / 1e6, b.e1 / 1e6), (b.c2, b.s2 / 1e6, b.e2 / 1e6), r1=R_LINK, r2=R_LINK,
                color=COL_FONDO_BLOQUE, alpha=0.6, lw=0, allow_twist=False)
for a, b in sorted(seg_pares):
    circos.link_line((POS.at[a, "crom"], POS.at[a, "mid"] / 1e6), (POS.at[b, "crom"], POS.at[b, "mid"] / 1e6),
                     r1=R_LINK, r2=R_LINK, color=COL_DESTACADO, lw=0.9, zorder=5)

H15 = ANCHO + 15 * MM      # clave de 3 filas al pie
fig = plt.figure(figsize=(ANCHO, H15))
axp = fig.add_axes([0, (H15 - ANCHO) / H15, 1, ANCHO / H15], projection="polar")
circos.plotfig(ax=axp)


def pava(y):
    # regresion isotonica (no decreciente) por bloques adyacentes
    bloques = []
    for v in y:
        bloques.append([v, 1])
        while len(bloques) > 1 and bloques[-2][0] > bloques[-1][0]:
            m2, w2 = bloques.pop()
            m1, w1 = bloques.pop()
            bloques.append([(m1 * w1 + m2 * w2) / (w1 + w2), w1 + w2])
    return np.concatenate([[m] * w for m, w in bloques])


def repartir(ang, sep):
    # posiciones angulares lo mas cercanas posible a 'ang' (minimos cuadrados) con una
    # separacion minima 'sep' entre rotulos vecinos; el circulo se corta en el hueco mayor.
    ang = np.asarray(ang, float) % 360
    orden = np.argsort(ang)
    a = ang[orden]
    huecos = np.diff(np.r_[a, a[0] + 360])
    k = int(np.argmax(huecos))
    idx = np.r_[orden[k + 1:], orden[:k + 1]]
    rot = np.r_[a[k + 1:], a[:k + 1] + 360]
    i = np.arange(len(rot))
    p = pava(rot - sep * i) + sep * i
    assert p[-1] - p[0] <= 360 - sep, "no caben los rotulos"
    out = np.empty(len(ang))
    out[idx] = p % 360
    return out


r_mm = R_LAB / 100 * (ANCHO / MM) / 2
SEP = np.degrees(8 * 1.2 / 72 * 25.4 / r_mm)    # interlinea de 8 pt a ese radio
lab = nac[nac.locus.isin(gen_seg)].copy()
lab["th"] = [np.degrees(circos.get_sector(c).x_to_rad(m / 1e6)) for c, m in zip(lab.crom, lab.mid)]
lab["th_l"] = repartir(lab.th, SEP)
for _, g in lab.iterrows():
    t0, t1 = np.radians(g.th), np.radians(g.th_l)
    axp.plot([t0, t0, t1], [R_IDEO[1] + 0.3, R_IDEO[1] + 1.1, R_LAB - 0.7], color=COL_NAC, lw=0.4,
             solid_capstyle="butt", clip_on=False)
    phi = (90 - g.th_l) % 360
    rot, ha = (phi - 180, "right") if 90 < phi < 270 else ((phi if phi <= 90 else phi - 360), "left")
    axp.text(t1, R_LAB, g.simbolo, rotation=rot, rotation_mode="anchor", ha=ha, va="center", fontsize=8,
             clip_on=False)
cifra("col_F15_rotulos", len(lab), "T47", "genes de pares segmentales rotulados")
cifra("col_F15_desplazamiento_max_grados", round(float(np.max(np.abs(((lab.th_l - lab.th + 180) % 360) - 180))), 2),
      "e3_12", "reparto de rotulos sin solaparse")

axl = fig.add_axes([0.02, 0.0, 0.96, (H15 - ANCHO) / H15])
axl.set_axis_off()
# Clave en el orden radial (de fuera adentro). La escala de las marcas y el techo de la
# pista de densidad salen de las constantes con que se dibujan, no se escriben a mano.
marca = dict(lw=0, marker="|", ms=7, mew=1.2)
h = [Patch(fc=COL_IDEOGRAMA, ec=COL_NAC, lw=0.4, label=f"Chromosome, ticks every {TICK_MB} Mb"),
     Line2D([], [], color=COL_DESTACADO, label="NAC: segmental (WGD)", **marca),
     Line2D([], [], color=COL_BLOQUE, label="NAC: tandem", **marca),
     Line2D([], [], color=COL_PROXIMAL, label="NAC: proximal", **marca),
     Line2D([], [], color=COL_NAC, label="NAC: dispersed or singleton", **marca),
     # en el anillo son trazos radiales cortos (1,4 de 5,2 unidades de las marcas NAC)
     Line2D([], [], color=COL_GRUPO, lw=0, marker="|", ms=3.5, mew=1.6, label="NAC genes within 200 kb"),
     Patch(fc=COL_NEUTRO, ec="none", label=f"Gene density, 0–{DMAX:.0f} genes per Mb"),
     Patch(fc=COL_FONDO_BLOQUE, ec="none", label="Collinear block (genome)"),
     Line2D([], [], color=COL_DESTACADO, lw=1.0, label="Segmental NAC pair")]
axl.legend(handles=h, ncol=3, loc="center", frameon=False, handlelength=1.6, columnspacing=1.6,
           handletextpad=0.5, borderaxespad=0, labelspacing=0.45)
cifra("col_F15_marcas_cada_Mb", TICK_MB, "e3_12", "marcas del ideograma")
texto_minimo(fig)
guardar(fig, "F15_circos_colinealidad_Egl")

# ==================================================================================
# F16 -- cariotipo E. grandis - E. globulus - P. trichocarpa (JCVI)
# ==================================================================================
MINSPAN = 30
ESP = {"Egr": "E. grandis", "Egl": "E. globulus", "Ptr": "P. trichocarpa"}
CROMS = {}
for sp in ("Egr", "Egl", "Ptr"):
    cs = sorted({c for c in bed[sp].crom if re.fullmatch(sp + r"\d{2}", c)})
    CROMS[sp] = cs
assert len(CROMS["Egr"]) == 11 and len(CROMS["Egl"]) == 11 and len(CROMS["Ptr"]) == 19, {k: len(v) for k, v in CROMS.items()}

rango = {sp: {g: (c, i) for c, gg in bed[sp].groupby("crom", sort=False) for i, g in enumerate(gg.sort_values("ini").gen)}
         for sp in bed}
enl_res = {}
for sp in ("Egr", "Ptr"):
    nuevo, simple = JF / f"Egl.{sp}.min{MINSPAN}.anchors", JF / f"Egl.{sp}.min{MINSPAN}.simple"
    r = subprocess.run([sys.executable, "-m", "jcvi.compara.synteny", "screen",
                        str(COL / "jcvi" / f"Egl.{sp}.lifted.anchors"), str(nuevo),
                        f"--qbed={ISO / 'Egl.bed'}", f"--sbed={ISO / (sp + '.bed')}",
                        f"--minspan={MINSPAN}", "--simple"], capture_output=True, text=True)
    assert r.returncode == 0 and simple.is_file(), r.stderr[-800:]
    fondo = [l for l in simple.read_text().splitlines() if l.strip()]
    bloques = [l.split("\t") for l in fondo]
    # Pares NAC: cada NAC de E. globulus con su mejor socio colineal (T48).
    o = t48.dropna(subset=[f"ortologo_{sp}"]).copy()
    o["es_nac"] = o[f"ortologo_{sp}_es_NAC"].astype(str).str.lower().eq("true")
    assert o.locus.isin(rango["Egl"]).all() and o[f"ortologo_{sp}"].isin(rango[sp]).all()

    def dentro(g1, g2):
        c1, i1 = rango["Egl"][g1]
        c2, i2 = rango[sp][g2]
        for a, b, c, d, *_ in bloques:
            ca, ia = rango["Egl"][a]
            cb, ib = rango["Egl"][b]
            cc, ic = rango[sp][c]
            cd, id_ = rango[sp][d]
            if ca == c1 and cc == c2 and min(ia, ib) <= i1 <= max(ia, ib) and min(ic, id_) <= i2 <= max(ic, id_):
                return True
        return False

    o["en_bloque"] = [dentro(a, b) for a, b in zip(o.locus, o[f"ortologo_{sp}"])]
    lineas = fondo[:]
    for es in (False, True):            # rojo encima del gris
        for _, g in o[o.es_nac == es].iterrows():
            col = COL_DESTACADO if es else COL_NAC
            q = g[f"ortologo_{sp}"]
            lineas.append(f"{col}*{g.locus}\t{g.locus}\t{q}\t{q}\t1\t+")
    (JF / f"Egl.{sp}.figura.simple").write_text("\n".join(lineas) + "\n")
    enl_res[sp] = o
    cifra(f"col_{sp}_bloques_min{MINSPAN}", len(fondo), f"jcvi screen --minspan={MINSPAN}",
          "bloques de fondo en F16")
    cifra(f"col_{sp}_enlaces_NAC", len(o), f"T48.ortologo_{sp}", "NAC de E. globulus con mejor socio colineal")
    cifra(f"col_{sp}_enlaces_socio_NAC", int(o.es_nac.sum()), f"T48.ortologo_{sp}_es_NAC")
    cifra(f"col_{sp}_enlaces_en_bloques_mostrados", int(o.en_bloque.sum()), "e3_12",
          f"enlaces dentro de un bloque de >= {MINSPAN} genes")

(JF / "seqids.txt").write_text("\n".join(",".join(CROMS[sp]) for sp in ("Egr", "Egl", "Ptr")) + "\n")
# Disposicion en mm. Con fracciones fijas de un eje de 85 mm quedaban 12,6 mm vacios sobre
# E. grandis; ahora la altura sale de la barra, la separacion entre pistas (sitio para las
# cintas), la clave al pie y un margen superior pequeno.
BARRA16, SEP16, MARGEN16, LEY16, HUECO16 = 3.4 * MM, 28 * MM, 2.5 * MM, 7 * MM, 3 * MM
H16 = LEY16 + HUECO16 + BARRA16 + 2 * SEP16 + MARGEN16
y_ptr = LEY16 + HUECO16 + BARRA16 / 2
Y = {"Ptr": y_ptr / H16, "Egl": (y_ptr + SEP16) / H16, "Egr": (y_ptr + 2 * SEP16) / H16}
ALTO_BARRA = BARRA16 / H16
X0, X1 = 0.135, 0.99
with open(JF / "layout.csv", "w") as f:
    f.write("# y, xstart, xend, rotation, color, label, va, bed\n")
    for sp in ("Egr", "Egl", "Ptr"):
        f.write(f"{Y[sp]}, {X0}, {X1}, 0, {COL_IDEOGRAMA}, , center, {ISO / (sp + '.bed')}\n")
    f.write("# edges\n")
    f.write(f"e, 1, 0, {JF / 'Egl.Egr.figura.simple'}\n")
    f.write(f"e, 1, 2, {JF / 'Egl.Ptr.figura.simple'}\n")

fig = plt.figure(figsize=(ANCHO, H16))
root = fig.add_axes([0, 0, 1, 1])
root.set_xlim(0, 1)
root.set_ylim(0, 1)
root.set_axis_off()
# generank=False: cromosomas a escala fisica (pb) y no por rango genico (valor por defecto
# de JCVI); sin rotulos ni circulos de JCVI, que se dibujan aqui a 8 pt. chrstyle="rect":
# con "auto" JCVI mezcla rectangulos (cromosomas cortos) y extremos redondeados cuyo radio va
# en unidades de datos de un eje no cuadrado, y salian deformados.
K = Karyotype(root, str(JF / "seqids.txt"), str(JF / "layout.csv"), height=ALTO_BARRA, lw=0.4,
              generank=False, heightpad=ALTO_BARRA / 2, plot_label=False, plot_circles=False,
              chrstyle="rect")
assert not root.texts, "JCVI ha dibujado texto propio"
for sp, tr in zip(("Egr", "Egl", "Ptr"), K.tracks):
    assert list(tr.seqids) == CROMS[sp]
    for sid in tr.seqids:
        x0 = tr.offsets[sid]
        root.text(x0 + tr.ratio * tr.sizes[sid] / 2, tr.y, str(int(sid[3:])), ha="center", va="center",
                  fontsize=8, color="black", zorder=10)
    root.text(X0 - 0.012, tr.y, ESP[sp], ha="right", va="center", fontsize=8, style="italic")
axl = fig.add_axes([0.02, 0.0, 0.96, LEY16 / H16])
axl.set_axis_off()
h = [Patch(fc=COL_GAINSBORO, ec="none", label=f"Collinear block (≥ {MINSPAN} genes)"),
     Line2D([], [], color=COL_DESTACADO, lw=1.0, label="NAC gene and its collinear NAC partner"),
     Line2D([], [], color=COL_NAC, lw=1.0, label="Collinear partner not detected as NAC")]
axl.legend(handles=h, ncol=3, loc="center", frameon=False, handlelength=1.8, columnspacing=1.6,
           handletextpad=0.5, borderaxespad=0)
texto_minimo(fig)
guardar(fig, "F16_cariotipo_Egr_Egl_Ptr")

# ==================================================================================
# F17 -- Ks de las anclas del genoma y Ks y Ka/Ks de los pares NAC
# ==================================================================================
# Pico de la densidad: misma funcion que e3_10_kaks_Egl.py (picos), comprobada contra
# su resumen para que la figura no muestre un pico distinto del que se publica.
ks_an = an.loc[an.estado == "ok", "Ks_YN00"].to_numpy(float)
ks_an = ks_an[np.isfinite(ks_an) & (ks_an > 0)]
x3 = ks_an[ks_an <= 3.0]
kde = gaussian_kde(x3)
rejilla = np.linspace(0, 3.0, 601)
dens = kde(rejilla)
ref = rk["anclas_genoma"]["picos_KDE_anclas_YN00_PAML"]
ipico = int(np.argmax(dens))
PICO = float(rejilla[ipico])
assert len(x3) == ref["n"] and round(float(kde.factor), 4) == ref["factor_banda"], (len(x3), kde.factor)
assert abs(PICO - ref["picos"][0]["Ks"]) < 1e-9 and len(ref["picos"]) == 1
T_PICO = PICO / (2 * LAMBDA) / 1e6
assert round(T_PICO, 1) == rk["anclas_genoma"]["T_Mya_del_pico_principal"]
XMAX = 6.0
cifra("kaks_anclas_estimadas", int(len(ks_an)), "Ks_anclas_genoma.tsv", "Ks > 0")
cifra("kaks_anclas_KDE_n", int(len(x3)), "Ks_anclas_genoma.tsv", "0 < Ks <= 3")
cifra("kaks_KDE_factor_banda", float(kde.factor), "scipy gaussian_kde (Scott)")
cifra("kaks_pico_Ks", PICO, "resumen_kaks_Egl.json")
cifra("kaks_pico_T_Mya", round(T_PICO, 1), "resumen_kaks_Egl.json", f"T = Ks / (2 x {LAMBDA:g})")
cifra("kaks_pico_mitad_altura", "–".join(f"{v:.3f}" for v in ref["picos"][0]["mitad_altura"]), "resumen_kaks_Egl.json")
cifra("kaks_anclas_Ks_mayor_xmax", int((ks_an > XMAX).sum()), "Ks_anclas_genoma.tsv", f"fuera de F17 (Ks > {XMAX:g})")
cifra("kaks_segmentales_en_pico", rk["segmentales_NAC"]["en_pico_principal_del_genoma"], "resumen_kaks_Egl.json",
      "dentro de la mitad de altura del pico")

CLASES = [("WGD/segmental", "Segmental", COL_DESTACADO),
          ("tandem", "Tandem", COL_BLOQUE),
          ("proximal", "Proximal", COL_PROXIMAL),
          ("a_200kb_sin_modo", "Within 200 kb,\nDIAMOND hit", COL_NAC),
          ("a_200kb_sin_acierto", "Within 200 kb,\nno DIAMOND hit", COL_NEUTRO)]
assert set(t47.tipo_par) == {c for c, *_ in CLASES} and len(t47) == rk["pares_NAC"]["analizados"]
assert t47.Ks_YN00.max() <= XMAX
YC = {c: len(CLASES) - 1 - i for i, (c, *_) in enumerate(CLASES)}
t47["fiable"] = (t47.Ks_YN00 >= 0.01) & (t47.Ks_YN00 <= KS_SAT) & np.isfinite(t47.KaKs_YN00)
rng = np.random.default_rng(12)

# Rotulos de fila de 3 lineas (7,9 mm a 8 pt): el panel inferior necesita ~9 mm por clase
# para que no se toquen. Los margenes van en mm para que no crezcan con la altura.
H17 = 104 * MM
fig = plt.figure(figsize=(ANCHO, H17))
gs = fig.add_gridspec(2, 2, width_ratios=[1.72, 1], height_ratios=[1, 1.55], left=0.165, right=0.985,
                      bottom=10.3 * MM / H17, top=1 - 12.3 * MM / H17, wspace=0.07, hspace=0.07)
axd = fig.add_subplot(gs[0, 0])
axs = fig.add_subplot(gs[1, 0], sharex=axd)
axk = fig.add_subplot(gs[1, 1], sharey=axs)
axl = fig.add_subplot(gs[0, 1])
axl.set_axis_off()

for ax in (axd, axs):
    ax.axvspan(KS_SAT, XMAX, fc="#F2F2F2", ec="none", zorder=0)
    ax.axvline(PICO, color=COL_TEXTO_SEC, lw=0.6, ls=(0, (3, 2)), zorder=1)
anch = 0.1
pesos = np.full((ks_an <= XMAX).sum(), 1 / (len(ks_an) * anch))
alturas, bordes_h, _ = axd.hist(ks_an[ks_an <= XMAX], bins=np.arange(0, XMAX + anch / 2, anch),
                                weights=pesos, fc="#CFCFCF", ec="white", lw=0.25, zorder=2)
curva = dens * len(x3) / len(ks_an)
axd.plot(rejilla, curva, color="black", lw=0.8, zorder=3)
axd.set_ylabel("Density")
# Techo con holgura para que los rotulos queden encima de la curva (se comprueba al final).
YTOP = 1.45 * max(float(alturas.max()), float(curva.max()))
axd.set_ylim(0, YTOP)
axd.tick_params(labelbottom=False)
rot_d = [axd.text(PICO + 0.07, YTOP * 0.98, f"Ks = {PICO:.2f}\n{T_PICO:.1f} Mya", ha="left", va="top",
                  fontsize=8, color="black"),
         axd.text(XMAX - 0.06, YTOP * 0.98, f"Ks > {KS_SAT:g}: not dated", ha="right", va="top",
                  fontsize=8, color=COL_TEXTO_SEC)]
sec = axd.secondary_xaxis("top", functions=(lambda k: k / (2 * LAMBDA) / 1e6, lambda t: t * 2 * LAMBDA * 1e6))
sec.set_xticks([0, 50, 100, 150])
sec.set_xlabel("Divergence time (Mya)")
sec.spines["top"].set_linewidth(0.5)

for c, etiqueta, col in CLASES:
    d = t47[t47.tipo_par == c]
    y = YC[c] + rng.uniform(-0.3, 0.3, len(d))
    # Mismo codigo que en b (la leyenda es comun): relleno = 0.01 <= Ks <= KS_SAT, hueco = fuera.
    m = d.fiable.to_numpy()
    axs.scatter(d.Ks_YN00[m], y[m], s=7, c=col, lw=0, alpha=0.8, zorder=3)
    axs.scatter(d.Ks_YN00[~m], y[~m], s=7, facecolors="none", edgecolors=col, lw=0.5, alpha=0.8, zorder=3)
    cifra(f"kaks_{c.replace('/', '_')}_n", len(d), "T47.tipo_par")
    cifra(f"kaks_{c.replace('/', '_')}_Ks_mediana", float(d.Ks_YN00.median()), "T47.Ks_YN00")
    cifra(f"kaks_{c.replace('/', '_')}_Ks_saturado", int((d.Ks_YN00 > KS_SAT).sum()), "T47.Ks_YN00", f"Ks > {KS_SAT:g}")
    cifra(f"kaks_{c.replace('/', '_')}_purificadora", int((d.seleccion_LRT == "purificadora").sum()), "T47.seleccion_LRT")
axs.set_yticks([YC[c] for c, *_ in CLASES])
axs.set_yticklabels([f"{e}\n(n = {(t47.tipo_par == c).sum()})" for c, e, _ in CLASES], linespacing=1.0)
axs.set_ylim(-0.6, len(CLASES) - 0.4)
axs.set_xlim(0, XMAX)
axs.set_xlabel("Ks (YN00)")
# pad de 7 pt: los puntos con Ks ~ 0 sobresalen a la izquierda del origen y tocaban el rotulo.
axs.tick_params(axis="y", length=0, pad=7)
axs.spines["left"].set_visible(False)

for c, etiqueta, col in CLASES:
    d = t47[(t47.tipo_par == c) & np.isfinite(t47.KaKs_YN00)]
    f_ = d[d.fiable]
    y0 = YC[c]
    if len(f_) >= 5:
        # manage_ticks=False: si no, cada caja anade su posicion como marca y rotulo al eje y
        # compartido con axs, y aparecia un numero suelto sobre el nombre de la clase.
        axk.boxplot(f_.KaKs_YN00, positions=[y0], orientation="horizontal", manage_ticks=False,
                    widths=0.55, showfliers=False, whis=1.5,
                    patch_artist=True, boxprops=dict(fc="none", ec="black", lw=0.6),
                    medianprops=dict(color="black", lw=0.9), whiskerprops=dict(lw=0.6), capprops=dict(lw=0.6),
                    zorder=4)
    yj = y0 + rng.uniform(-0.3, 0.3, len(d))
    m = d.fiable.to_numpy()
    axk.scatter(d.KaKs_YN00[m], yj[m], s=7, c=col, lw=0, alpha=0.8, zorder=3)
    axk.scatter(d.KaKs_YN00[~m], yj[~m], s=7, facecolors="none", edgecolors=col, lw=0.5, alpha=0.8, zorder=3)
    cifra(f"kaks_{c.replace('/', '_')}_KaKs_n_fiable", len(f_), "T47", f"0.01 <= Ks <= {KS_SAT:g}")
    cifra(f"kaks_{c.replace('/', '_')}_KaKs_mediana_fiable", float(f_.KaKs_YN00.median()) if len(f_) else "NA", "T47.KaKs_YN00")
axk.axvline(1, color=COL_TEXTO_SEC, lw=0.6, ls=(0, (3, 2)), zorder=1)
axk.set_xlim(0, max(3.0, float(np.nanmax(t47.KaKs_YN00)) + 0.1))
axk.set_xlabel("Ka/Ks (YN00)")
axk.tick_params(axis="y", left=False, labelleft=False)
axk.spines["left"].set_visible(False)

h = [Patch(fc="#CFCFCF", ec="none", label="Genome collinear anchors"),
     Line2D([], [], color="black", lw=0.8, label="Kernel density (Ks ≤ 3)"),
     Line2D([], [], color=COL_TEXTO_SEC, lw=0.6, ls=(0, (3, 2)), label="Density peak; Ka/Ks = 1"),
     Line2D([], [], color=COL_NAC, lw=0, marker="o", ms=3.2, label=f"NAC pair, 0.01 ≤ Ks ≤ {KS_SAT:g}"),
     Line2D([], [], color=COL_NAC, lw=0, marker="o", ms=3.2, mfc="none", mew=0.6, label="NAC pair, Ks outside"),
     Patch(fc="none", ec="black", lw=0.6, label="Box plot: filled pairs only")]
axl.legend(handles=h, loc="center left", frameon=False, handlelength=1.8, handletextpad=0.6,
           borderaxespad=0, bbox_to_anchor=(0.04, 0.5), labelspacing=0.45)
for ax, letra in ((axd, "a"), (axk, "b")):
    bb = ax.get_position()
    fig.text(0.012 if letra == "a" else bb.x0 - 0.01, bb.y1 + (0.075 if letra == "a" else 0.012), letra,
             fontsize=PT_LETRA, fontweight="bold", ha="left", va="bottom")
# Comprobaciones de maquetacion: rotulos del pico por encima de la densidad, una etiqueta
# por clase y rotulos de fila sin solaparse.
fig.canvas.draw()
_r = fig.canvas.get_renderer()
for t in rot_d:
    bb = t.get_window_extent(_r).transformed(axd.transData.inverted())
    en = (rejilla >= bb.x0) & (rejilla <= bb.x1)
    barras = alturas[(bordes_h[1:] > bb.x0) & (bordes_h[:-1] < bb.x1)]
    techo = max([0.0] + ([float(curva[en].max())] if en.any() else []) + ([float(barras.max())] if len(barras) else []))
    assert techo < bb.y0, f"F17: {t.get_text()!r} pisa la densidad ({techo:.3f} >= {bb.y0:.3f})"
et = [t for t in axs.get_yticklabels() if t.get_text()]
assert len(axs.get_yticks()) == len(CLASES) and len(et) == len(CLASES), [t.get_text() for t in et]
cajas = sorted((t.get_window_extent(_r) for t in et), key=lambda b: b.y0)
holgura = min(b.y0 - a.y1 for a, b in zip(cajas, cajas[1:])) / fig.dpi * 25.4
assert holgura > 0.5, f"F17: rotulos de fila a {holgura:.2f} mm"
# Rotulos frente a los puntos: el punto mas a la izquierda (s = 7 pt^2) menos su radio.
x_izq = axs.transData.transform((float(t47.Ks_YN00.min()), 0))[0] - np.sqrt(7) / 2 * fig.dpi / 72
hueco_x = (x_izq - max(b.x1 for b in cajas)) / fig.dpi * 25.4
assert hueco_x > 0.8, f"F17: rotulos de fila a {hueco_x:.2f} mm de los puntos"
bl = axl.get_legend().get_window_extent(_r)
for t in fig.texts:
    if t.get_text() in ("a", "b"):
        assert not bl.overlaps(t.get_window_extent(_r)), f"F17: la leyenda pisa la letra {t.get_text()}"
m_ab = t47.fiable.to_numpy()
assert m_ab.sum() == 418 and (~m_ab).sum() == 197, (m_ab.sum(), (~m_ab).sum())
print(f"  F17: hueco minimo entre rotulos de fila {holgura:.2f} mm; rotulo-punto {hueco_x:.2f} mm; "
      f"pares rellenos {m_ab.sum()}, huecos {(~m_ab).sum()}")
texto_minimo(fig)
guardar(fig, "F17_KaKs_datacion")

# --- cifras ---------------------------------------------------------------------------
cf = pd.DataFrame(CIF)
assert cf.clave.is_unique
cf.to_csv(CIFRAS / "e3_colinealidad.csv", index=False)
print(f"  {len(cf)} cifras -> _cifras/e3_colinealidad.csv")
print(cf.to_string(index=False, max_colwidth=60))
