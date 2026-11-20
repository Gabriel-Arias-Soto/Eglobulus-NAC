#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""F10 - figura estructural 3D del gen focal: tres paneles trazados con PyMOL.

    Uso:  cd <raiz del proyecto> && pymol -cq scripts/13_estructura3D.py

Es la version mantenible de estructura_EglNAC34/componer_figura.py + paneles.pml,
que fue lo que produjo figura_modelo_EglNAC34.png. Tres cambios respecto de aquella:

  1. Nomenclatura. Los rotulos salen de GEN_FOCAL (00_estilo.R, regla 4), no del
     simbolo del manuscrito. Ver NOMENCLATURA_decision.md.
  2. Las selecciones de residuos ya no van tecleadas dentro del .pml: se derivan
     de T19 (contactos atomo a atomo en 3SWP) y de qc_perfil_dos_referencias.csv
     (tramos divergentes). Si cambia la tabla, cambia la figura.
  3. El panel c separa dos cosas que antes iban del mismo color: los residuos que
     leen bases nitrogenadas -los que dan especificidad de secuencia- y los que
     solo agarran el esqueleto de fosfatos. Esa distincion es el resultado; el
     resto del panel es contexto.

El ADN no es un modelo: es el duplex del co-cristal ANAC019-ADN (3SWP), colocado
por superposicion del dominio sobre la cadena homologa. No hay prediccion aqui.
"""
import os, re, csv
import numpy as np, pandas as pd, gemmi
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from pymol import cmd

# --- rutas y contrato de estilo -------------------------------------------
RAIZ = os.getcwd()
assert os.path.isfile(os.path.join(RAIZ, "scripts", "00_estilo.R")), \
    "ejecutar desde la raiz del proyecto: pymol -cq scripts/13_estructura3D.py"
ANALISIS = os.path.join(RAIZ, "ANALISIS_Eglobulus")
TABLAS   = os.path.join(ANALISIS, "tablas")
FIGURAS  = os.path.join(ANALISIS, "figuras")
EST      = os.path.join(RAIZ, "estructura_EglNAC34")   # modelo, cristal recortado y vistas de camara
REN      = os.path.join(FIGURAS, "_render3D")
os.makedirs(REN, exist_ok=True)

_estilo = open(os.path.join(RAIZ, "scripts", "00_estilo.R"), encoding="utf-8").read()
def cte(nombre):
    m = re.search(r'^%s\s*<-\s*"?([^"#\n]+?)"?\s*(?:#.*)?$' % nombre, _estilo, re.M)
    assert m, "no encuentro %s en 00_estilo.R" % nombre
    v = m.group(1).strip()
    # Las medidas de la revista se declaran en mm ("174 * MM") para que la
    # cifra de referencia sea la del reglamento y no su conversion. Aqui se
    # resuelve esa forma, porque float() sola no la entiende.
    mm = re.fullmatch(r"([\d.]+)\s*\*\s*MM", v)
    return str(float(mm.group(1)) / 25.4) if mm else v
ANCHO = float(cte("ANCHO_COLUMNA"))
DPI   = int(float(cte("DPI")))
GEN   = cte("GEN_FOCAL")

# Springer pide Helvetica o Arial y rotulacion de 8 a 12 pt. matplotlib no pasa
# por theme_nac(), asi que la misma regla se repite aqui. Liberation Sans es
# el sustituto metricamente identico a Arial cuando esta no existe.
import glob as _glob
import matplotlib.font_manager as _fm
# Liberation Sans viaja dentro del paquete R fontLiberation, que no es un
# directorio de fuentes del sistema: matplotlib no la ve si no se le da.
for _f in _glob.glob(os.path.expanduser(
        "~/conda/envs/*/lib/R/library/fontLiberation/"
        "fonts/liberation-fonts/LiberationSans-*.ttf")):
    try:
        _fm.fontManager.addfont(_f)
    except Exception:
        pass
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "Liberation Sans",
                                   "DejaVu Sans"]

# Paleta del render. No vive en 00_estilo.R porque PyMOL no lee R; se mantiene
# aparte de PAL_PLDDT a proposito: aqui el rojo significa confianza baja (panel a)
# y por eso los residuos de contacto van en azul y nunca en rojo, para que el
# lector no lea "contacto" donde dice "mal modelado".
C_LO, C_HI = "#C1440E", "#1B6CA8"    # extremos de la escala de pLDDT
C_BASE     = "#08519C"               # residuos que leen bases
C_ESQ      = "#9ECAE1"               # residuos que solo tocan el esqueleto
C_GRIS     = "#B8B8B8"
C_ADN      = "#B08C69"
def rgb(h): return [int(h[i:i+2], 16) / 255 for i in (1, 3, 5)]

# --- de donde sale cada seleccion -----------------------------------------
t19 = pd.read_csv(os.path.join(TABLAS, "T19_contactos_proteina_ADN_3SWP.csv"))
t18 = pd.read_csv(os.path.join(TABLAS, "T18_lectura_bases_3SWP.csv"))
dos = pd.read_csv(os.path.join(RAIZ, "qc_perfil_dos_referencias.csv"))
qc  = pd.read_csv(os.path.join(RAIZ, "qc_rmsd_modelo_vs_plantillas.csv"))

res_base = sorted(t19.loc[t19.parte == "base", "pri"].astype(int).unique())
res_todos = sorted(t19.pri.astype(int).unique())
res_esq  = [r for r in res_todos if r not in res_base]
# T18 se calculo por otra via (posiciones de ANAC019 mapeadas al gen focal).
# Que las dos rutas den el mismo conjunto es la comprobacion de que el mapeo
# de numeracion es correcto; si algun dia difieren, la figura no debe salir.
assert res_base == sorted(int(re.sub(r"\D", "", s)) for s in t18.residuo_EglNAC), \
    "T18 y T19 no coinciden en que residuos leen bases"
div = dos.loc[dos.tramo_divergente, "residuo"].astype(int).tolist()
sel = lambda v: "+".join(str(x) for x in v)

# --- panel a, b, c: trazado ------------------------------------------------
cmd.bg_color("white")
for k, v in dict(ray_opaque_background=1, cartoon_transparency=0,
                 antialias=2, ray_trace_mode=0, specular=0.15, ambient=0.25,
                 direct=0.45, reflect=0.35, cartoon_fancy_helices=1,
                 cartoon_smooth_loops=0, stick_radius=0.16).items():
    cmd.set(k, v)
for nom, hexa in [("plo", C_LO), ("phi", C_HI), ("pgris", C_GRIS), ("padn", C_ADN),
                  ("pbase", C_BASE), ("pesq", C_ESQ)]:
    cmd.set_color(nom, rgb(hexa))
# Las dos matrices de camara se heredan de la figura anterior para que esta sea
# comparable con ella: misma orientacion, mismo encuadre.
VISTA_PROT = open(os.path.join(EST, "vista_prot.pml")).read()
VISTA_ADN  = open(os.path.join(EST, "vista_adn.pml")).read()
MODELO = os.path.join(EST, "modelo.pdb")

def panel_a():
    cmd.delete("all"); cmd.load(MODELO, "mod"); cmd.hide("everything")
    cmd.show("cartoon", "mod")
    cmd.spectrum("b", "plo_white_phi", "mod", 50, 100)   # B-factor = pLDDT
    cmd.do(VISTA_PROT); cmd.zoom("mod", 2); cmd.deselect()
    cmd.png(os.path.join(REN, "panel_a.png"), width=2250, height=2250, dpi=300, ray=1)

def panel_b():
    cmd.delete("all"); cmd.load(MODELO, "mod")
    cmd.load(os.path.join(EST, "1UT7_B_recortado.pdb"), "cris")
    cmd.hide("everything"); cmd.show("cartoon", "mod or cris")
    cmd.color("pgris", "cris"); cmd.color("phi", "mod")
    cmd.color("plo", "mod and resi %s" % sel(div))
    cmd.do(VISTA_PROT); cmd.zoom("mod", 2); cmd.deselect()
    cmd.png(os.path.join(REN, "panel_b.png"), width=2250, height=2250, dpi=300, ray=1)

def panel_c():
    cmd.delete("all"); cmd.load(MODELO, "mod")
    cmd.load(os.path.join(EST, "ADN_recortado.pdb"), "adn")
    cmd.hide("everything"); cmd.show("cartoon", "mod"); cmd.color("pgris", "mod")
    for nombre, lista, color in (("cb", res_base, "pbase"), ("ce", res_esq, "pesq")):
        cmd.select(nombre, "mod and resi %s" % sel(lista))
        cmd.show("sticks", "%s and not (name C+N+O)" % nombre)
        cmd.color(color, nombre)
    cmd.show("cartoon", "adn"); cmd.show("sticks", "adn")
    cmd.set("cartoon_nucleic_acid_mode", 4)
    cmd.set("cartoon_ring_mode", 3); cmd.set("cartoon_ring_finder", 1)
    cmd.color("padn", "adn"); cmd.set("stick_radius", 0.12, "adn")
    cmd.do(VISTA_ADN); cmd.zoom("mod or adn", 2); cmd.deselect()
    cmd.png(os.path.join(REN, "panel_c.png"), width=2250, height=2550, dpi=300, ray=1)

panel_a(); panel_b(); panel_c()

# --- composicion -----------------------------------------------------------
def recortar(path, margen=8):
    """Quita el fondo blanco sobrante: PyMOL centra el objeto en un lienzo
    cuadrado y deja franjas que descuadran el mosaico."""
    im = plt.imread(path)
    if im.ndim == 2:
        im = np.dstack([im] * 3)
    im = im[..., :3]
    fondo = (im > 0.985).all(-1)
    fil = np.where(~fondo.all(1))[0]; col = np.where(~fondo.all(0))[0]
    y0, y1 = max(fil[0]-margen, 0), min(fil[-1]+margen+1, im.shape[0])
    x0, x1 = max(col[0]-margen, 0), min(col[-1]+margen+1, im.shape[1])
    return im[y0:y1, x0:x1]

def con_banda(im, frac):
    """Anade papel en blanco debajo del render. Las claves y las cifras van ahi
    y no encima de la molecula: superpuestas tapan cartoon y obligan a leer
    texto sobre un fondo que cambia de color de una letra a otra."""
    h = int(round(im.shape[0] * frac))
    return np.vstack([im, np.ones((h, im.shape[1], 3), im.dtype)]), h

# Alto de la franja como fraccion del render. El panel b pide mas porque apila
# la clave sobre las dos lineas de RMSD; el c pide menos porque ya es alto.
pa, ha = con_banda(recortar(os.path.join(REN, "panel_a.png")), 0.26)
pb, hb = con_banda(recortar(os.path.join(REN, "panel_b.png")), 0.38)
pc, hc = con_banda(recortar(os.path.join(REN, "panel_c.png")), 0.12)

# a y b comparten el ancho de columna del mosaico. Si sus proporciones no
# coinciden, matplotlib centra el mas ancho dentro de su caja y deja margen
# arriba y abajo: la franja dejaria de empezar donde dicen estas cuentas y los
# rotulos volverian a caer sobre el render. Se igualan con mas papel.
def igualar(im, h_banda, objetivo):
    extra = int(round(objetivo * im.shape[1])) - im.shape[0]
    if extra <= 0:
        return im, h_banda
    return np.vstack([im, np.ones((extra, im.shape[1], 3), im.dtype)]), h_banda + extra
_obj = max(pa.shape[0]/pa.shape[1], pb.shape[0]/pb.shape[1])
pa, ha = igualar(pa, ha, _obj)
pb, hb = igualar(pb, hb, _obj)

asp = [p.shape[0] / p.shape[1] for p in (pa, pb, pc)]
frac_banda = {"a": ha/pa.shape[0], "b": hb/pb.shape[0], "c": hc/pc.shape[0]}

# Cifras: calculadas, nunca escritas a mano (regla 5 de 00_estilo.R).
q1 = qc.loc[qc.referencia == "1UT7-B"].iloc[0]
st = gemmi.read_structure(MODELO)
pl = np.array([[a.b_iso for a in r if a.name == "CA"][0] for r in st[0][0]])
pl_med, pl_alta = pl.mean(), 100 * (pl >= 70).mean()
nucleo = dos.loc[~dos.tramo_divergente, "desviacion_1UT7_A"]
n_nt = sum(len(c) for c in gemmi.read_structure(os.path.join(EST, "ADN_recortado.pdb"))[0])

W = ANCHO
H = W / (1/(2*asp[0]) + 1/asp[2] + 0.02)
fig = plt.figure(figsize=(W, H), dpi=DPI)
gs = fig.add_gridspec(2, 2, width_ratios=[(H/2)/asp[0], H/asp[2]],
                      wspace=0.02, hspace=0.03, left=0.005, right=0.995,
                      top=0.995, bottom=0.005)
ax_a, ax_b = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[1, 0])
ax_c = fig.add_subplot(gs[:, 1])
for ax, im in ((ax_a, pa), (ax_b, pb), (ax_c, pc)):
    ax.imshow(im); ax.set_axis_off()
for ax, t in ((ax_a, "a"), (ax_b, "b"), (ax_c, "c")):
    ax.text(0.005, 0.995, t, transform=ax.transAxes, va="top", ha="left",
            fontsize=9, fontweight="bold", color="black")
fig.canvas.draw()   # sin renderer no se puede medir el alto real de cada eje

# La maquetacion de los rotulos esta calibrada en pulgadas para un cuerpo de
# 5,6 pt, y Springer exige un minimo de 8 pt. Subir solo la fuente haria que la
# clave del panel b pisara el bloque de RMSD y que el titulo de la barra de color
# chocara con la linea de pLDDT, asi que los desplazamientos se escalan por el
# mismo factor que el cuerpo: la geometria relativa se conserva exacta.
FS_ORIG = 5.6
FS      = 8.0
K       = FS / FS_ORIG
PASO    = 0.088 * K
def alto_in(ax):  return ax.get_window_extent().height / fig.dpi
def y_frac(ax, pulgadas): return pulgadas / alto_in(ax)
def cabe(ax, panel, pulgadas):
    """La franja es finita. Si un rotulo se sale por arriba vuelve a montarse
    sobre la molecula y nadie se entera hasta mirar la figura: mejor reventar."""
    disp = alto_in(ax) * frac_banda[panel]
    assert pulgadas <= disp, "la franja de %s necesita %.2f in y solo tiene %.2f" % (
        panel, pulgadas, disp)

def clave(ax, panel, entradas, y0=0.05 * K, x=0.04):
    cabe(ax, panel, y0 + PASO*(len(entradas)-1) + 0.05 * K)
    for i, (col, txt) in enumerate(entradas):
        y = y_frac(ax, y0 + i*PASO)
        ax.plot([x, x+0.045], [y]*2, color=col, lw=2.4, solid_capstyle="round",
                transform=ax.transAxes, clip_on=False)
        ax.text(x+0.062, y, txt, transform=ax.transAxes, va="center", ha="left",
                fontsize=FS, color="#222222")

# panel a: barra de confianza abajo, cifras encima
cabe(ax_a, "a", 0.30 * K)
cmap = LinearSegmentedColormap.from_list("plddt", [C_LO, "#F2F2F2", C_HI])
cax = ax_a.inset_axes([0.06, y_frac(ax_a, 0.080 * K), 0.34, y_frac(ax_a, 0.030 * K)])
mpl.colorbar.ColorbarBase(cax, cmap=cmap, norm=mpl.colors.Normalize(50, 100),
                          orientation="horizontal", ticks=[50, 75, 100])
# el suelo de 8 pt tambien rige para los numeros de la barra
cax.tick_params(labelsize=FS, length=2, pad=1.5)
cax.set_title("pLDDT", fontsize=FS, pad=2.0)
ax_a.text(0.06, y_frac(ax_a, 0.255 * K), "mean pLDDT %.0f; %.0f %% of residues \u2265 70"
          % (pl_med, pl_alta), transform=ax_a.transAxes, va="center", ha="left",
          fontsize=FS, color="#222222")

# panel b: cifras en las dos lineas de abajo, clave encima. Apiladas y no lado a
# lado porque juntas no caben en el ancho de columna.
# A 8 pt las cinco lineas del panel b (2 de cifras + 3 de clave) no caben en
# el hueco que deja el render, y frac_banda no es ajustable: se MIDE de la
# imagen, de modo que forzarlo montaria el texto sobre la molecula. Se gana el
# renglon fundiendo las dos lineas de cifras en una; no se pierde ningun dato.
cabe(ax_b, "b", 0.49 * K)
ax_b.text(0.04, y_frac(ax_b, 0.045 * K), "RMSD %.2f \u00c5 (%d C\u03b1); \u03b2 core median %.2f \u00c5"
          % (q1.rmsd_A, int(q1.n_pares), nucleo.median()),
          transform=ax_b.transAxes, va="bottom", ha="left", fontsize=FS, color="#222222")
clave(ax_b, "b", [(C_GRIS, "free ANAC019 (1UT7-B)"), (C_LO, "divergent segments"),
                  (C_HI, "%s model" % GEN)], y0=0.285 * K - PASO)

clave(ax_c, "c", [(C_ADN, "3SWP duplex (%d nt)" % n_nt),
                  (C_ESQ, "%d contacting the backbone" % len(res_esq)),
                  (C_BASE, "%d residues contacting bases" % len(res_base)),
                  (C_GRIS, "NAC domain of %s" % GEN)])

sal = os.path.join(FIGURAS, "F10_estructura3D_%s.png" % GEN)
fig.savefig(sal, dpi=DPI, facecolor="white")
# matplotlib escribe PNG con canal alfa; Springer pide RGB de 8 bits
# por canal y el alfa da problemas en preimpresion. Se aplana
# sobre blanco conservando la resolucion declarada.
from PIL import Image
_im = Image.open(sal)
if _im.mode != "RGB":
    _fondo = Image.new("RGB", _im.size, "white")
    _fondo.paste(_im, mask=_im.split()[-1] if "A" in _im.mode else None)
    _fondo.save(sal, dpi=(DPI, DPI))
_im.close()
print("  guardada %-34s %.2f x %.2f in" % (os.path.basename(sal), W, H))

# --- cifras para las leyendas ---------------------------------------------
cif = os.path.join(TABLAS, "_cifras")
os.makedirs(cif, exist_ok=True)
cols = list(pd.read_csv(os.path.join(cif, "f08.csv"), nrows=1).columns) \
       if os.path.isfile(os.path.join(cif, "f08.csv")) else ["clave", "valor", "fuente"]
filas = [
    ("f10_plddt_medio",     "%.0f" % pl_med,   "B-factor de modelo.pdb"),
    ("f10_plddt_alta_pct",  "%.0f %%" % pl_alta, "residuos con pLDDT >= 70"),
    ("f10_rmsd_1ut7",       "%.2f \u00c5 (%d C\u03b1)" % (q1.rmsd_A, int(q1.n_pares)),
                            "qc_rmsd_modelo_vs_plantillas.csv"),
    ("f10_nucleo_mediana",  "%.2f \u00c5" % nucleo.median(), "qc_perfil_dos_referencias.csv"),
    ("f10_contactos_n",     len(res_todos),    "T19: residuos con algun contacto"),
    ("f10_lectores_base_n", len(res_base),     "T19 parte=base; verificado contra T18"),
    ("f10_lectores_base",   ", ".join(str(r) for r in res_base), "T19 parte=base"),
    ("f10_esqueleto_n",     len(res_esq),      "T19: contacto solo con fosfato o azucar"),
    ("f10_adn_nt",          n_nt,              "ADN_recortado.pdb, duplex de 3SWP"),
]
with open(os.path.join(cif, "f10.csv"), "w", newline="", encoding="utf-8") as fh:
    w = csv.writer(fh); w.writerow(cols)
    for k, v, f in filas:
        d = {"clave": k, "valor": v, "fuente": f}
        w.writerow([d.get(c, "") for c in cols])
print("  %d cifras -> _cifras/f10.csv" % len(filas))
