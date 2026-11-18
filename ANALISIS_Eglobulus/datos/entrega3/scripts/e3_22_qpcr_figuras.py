#!/usr/bin/env python
"""Figuras de qRT-PCR rehechas desde los datos por réplica (cierre del manuscrito, 26-10-2026).

F21_qpcr_NAC: los cuatro NAC (EglNAC33, EglNAC66, EglNAC85, EglNAC115) en sequía (a) y en frío (b).
F22_qpcr_blancos: los cuatro genes blanco (ZHD11, LTI65, NCED3, CLPD) en sequía (a) y en frío (b).
Sustituyen a las Fig. 3 y 4 (NAC) y 7 y 8 (blancos) de la versión sometida, de las que no hay imagen en el Word.

Datos
  NAC: tablas/T08 (RQ = 2^-(ΔCq - media ΔCq del control), recalculado desde los Cq crudos).
  Blancos: «2da entrega de datos/graficos.xlsx», hojas ZHD11, LTI65, NCED3, CLPD y sus «-c»: valores de expresión relativa
  por réplica que usaron los autores para las figuras originales (no hay Cq crudos de los blancos). Cada bloque de genotipo
  tiene 9 filas: control, tratamiento 1 y tratamiento 2, tres réplicas cada uno.
Clave de genotipos de las placas de NAC (deducida y comprobada aquí, no supuesta)
  Las p de Kruskal-Wallis y Dunn de cada hoja de placa (T09) se cruzan con las tablas publicadas de los autores
  («Analisis K_W_FRIO (1).xlsx», «Analisis K_W_sequia.xlsx»), y la clave resultante se verifica además comparando los RQ de
  T08 con los valores de graficos.xlsx de los NAC, réplica a réplica. Si alguna de las dos comprobaciones no cuadra, el guion
  para. Las hojas sin genotipo publicado (G2 en frío; G2 y G3 en sequía) no se dibujan.
Estadística (la del texto de Métodos)
  Por gen y genotipo: Kruskal-Wallis entre los tres tratamientos y Dunn sin corrección (con corrección por empates);
  letras de grupo con p <= 0,05 en Dunn, solo cuando el Kruskal-Wallis del genotipo da p <= 0,05 (decisión de
  Gabriel, 02-11-2026; si no, todos los tratamientos llevan «a»). Las p de los NAC tienen que reproducir T09.
Salidas
  figuras/F21_qpcr_NAC.{png,pdf}, figuras/F22_qpcr_blancos.{png,pdf}
  tablas/T63_qpcr_estadistica_figuras.csv (gen, genotipo, KW, Dunn, letras, media, EE)
  datos/entrega3/qpcr/resumen_qpcr_figuras.json (clave de genotipos y comprobaciones)
"""
import itertools, json, os, re
import numpy as np, pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager

A = os.path.dirname(os.path.abspath(__file__))
while os.path.basename(A) != "ANALISIS_Eglobulus":
    A = os.path.dirname(A)
RZ = os.path.dirname(A)
T, FIG = os.path.join(A, "tablas"), os.path.join(A, "figuras")
D2 = os.path.join(RZ, "2da entrega de datos")
OUT = os.path.join(A, "datos", "entrega3", "qpcr")
os.makedirs(OUT, exist_ok=True)

TRAT = {"sequia": ["CON", "IDS", "SDS"], "frio": ["NA", "CABF", "CAAF"]}
GENOS = {"sequia": ["GDT1", "GDT2", "GDS1", "GDS2", "GDC"], "frio": ["GCT1", "GCT2", "GCS1", "GCS2", "GDC"]}
PLACA_A_TRAT = {"sequia": {"NA": "CON", "AAH": "IDS", "ADH": "SDS"}, "frio": {"NA": "NA", "AAH": "CABF", "ADH": "CAAF"}}
# símbolo de las placas (T08.gen, el del File S1) -> símbolo del censo (T02)
t02 = pd.read_csv(os.path.join(T, "T02_equivalencia_nomenclaturas_189.csv"))
q4 = t02[t02.medido_por_qPCR.astype(str).str.lower().eq("true")]
VIEJO = dict(zip(q4.etiqueta_FileS1, q4.simbolo))
assert len(VIEJO) == 4
ORDEN_NAC = ["EglNAC33", "EglNAC66", "EglNAC85", "EglNAC115"]
assert set(VIEJO.values()) == set(ORDEN_NAC)
HOJA_GRAF = {"EglNAC17": "EglNAC18", "EglNAC34": "EglNAC34", "EglNAC47": "EglNAC47", "EglNAC62": "EglNAC62"}  # graficos.xlsx
BLANCOS = ["ZHD11", "LTI65", "NCED3", "CLPD"]


# --- Dunn sin corrección -------------------------------------------------------------------------------
def kw_dunn(grupos):
    """grupos: lista de arrays (mismo orden que TRAT). Devuelve p de KW y dict de p de Dunn (bilateral, sin corrección)."""
    x = np.concatenate(grupos); n = len(x)
    r = stats.rankdata(x)
    k = [len(g) for g in grupos]
    rm, i0 = [], 0
    for kk in k:
        rm.append(r[i0:i0 + kk].mean()); i0 += kk
    _, cnt = np.unique(x, return_counts=True)
    corr = (cnt ** 3 - cnt).sum() / (12 * (n - 1))
    p_kw = stats.kruskal(*grupos).pvalue
    dunn = {}
    for a, b in itertools.combinations(range(len(grupos)), 2):
        se = np.sqrt((n * (n + 1) / 12 - corr) * (1 / k[a] + 1 / k[b]))
        z = (rm[a] - rm[b]) / se
        dunn[(a, b)] = 2 * stats.norm.sf(abs(z))
    return p_kw, dunn


def letras(medias, dunn, alfa=0.05):
    """Letras de grupo (inserción y absorción) a partir de las p de Dunn; la «a» va al grupo de mayor media."""
    m = len(medias)
    dif = {frozenset(par) for par, p in dunn.items() if p <= alfa}
    cols = [set(range(m))]
    for par in dif:
        i, j = tuple(par)
        nuevas = []
        for c in cols:
            if i in c and j in c:
                nuevas += [c - {i}, c - {j}]
            else:
                nuevas.append(c)
        cols = [c for c in nuevas if not any(c < d for d in nuevas)]
        cols = [c for k_, c in enumerate(cols) if c not in cols[:k_]]
    orden = sorted(cols, key=lambda c: -max(medias[i] for i in c))
    out = ["" for _ in range(m)]
    for li, c in enumerate(orden):
        for i in c:
            out[i] += "abcdefgh"[li]
    return out


# --- 1. clave de genotipos de las placas de NAC ---------------------------------------------------------------
t08 = pd.read_csv(os.path.join(T, "T08_qpcr_expresion_relativa.csv"))
t08["G"] = t08.hoja.str.extract(r"-(G\d+(?:_2018)?)")[0]
t08["trat"] = [PLACA_A_TRAT[e][g if isinstance(g, str) else "NA"] for e, g in zip(t08.ensayo, t08.grupo)]
# hoja NAC18-G8(S2): su detector viene rotulado NAC35-C en el fichero de Cq, pero su Cq diana es distinto del de
# NAC35-G8(S2) (misma referencia) y sus p son las de EglNAC17 / GDS en la tabla publicada: es EglNAC17.
m18 = t08.hoja == "NAC18-G8(S2)"
a_ = t08[m18].Cttgt.values; b_ = t08[t08.hoja == "NAC35-G8(S2)"].Cttgt.values
assert len(a_) == len(b_) == 9 and not np.allclose(a_, b_)
t08.loc[m18, "gen"] = "EglNAC17"
t09 = pd.read_csv(os.path.join(T, "T09_qpcr_KW_Dunn.csv"))
t09["G"] = t09.genotipo_hoja.str.extract(r"-(G\d+(?:_2018)?)")[0]
t09.loc[t09.genotipo_hoja == "NAC18-G8(S2)", "gen"] = "EglNAC17"
COLS = {"frio": ["p_KW", "CON-CABF", "CON-CAAF", "CABF-CAAF"], "sequia": ["p_KW", "CON-IDS", "CON-SDS", "IDS-SDS"]}
PUB = {"frio": ("Analisis K_W_FRIO (1).xlsx", {"EglNac18": "EglNAC17"}),
       "sequia": ("Analisis K_W_sequia.xlsx", {"EglNac18": "EglNAC17", "EglNAC62": "EglNAC62", "EglNAC34": "EglNAC34"})}
cerca = lambda a, b: all(abs(x - y) <= 6e-4 * max(1, abs(y)) + 5e-4 for x, y in zip(a, b))
CLAVE, COMPROB = {}, {}
for ens, (fx, hojas) in PUB.items():
    xl = pd.ExcelFile(os.path.join(D2, fx))
    filas = []
    for h, gen in hojas.items():
        d = xl.parse(h, header=None)
        vals = [[float(v) for v in r.values if isinstance(v, (int, float)) and not pd.isna(v)][:4] for _, r in d.iterrows()
                if any(isinstance(v, str) and re.match(r"^G[CD][TS]?\d?$|^GDS$", v.strip()) for v in r.values[:2])]
        assert len(vals) == 5, (ens, h)
        filas += [(gen, GENOS[ens][i], v) for i, v in enumerate(vals)]      # posición = orden de genotipos de Métodos
    Gs = sorted(t09[t09.ensayo == ens].G.unique())
    sols = []
    for perm in itertools.permutations(Gs, 5):
        mp = dict(zip(GENOS[ens], perm)); n_ok = 0
        for gen, gt, v in filas:
            s = t09[(t09.ensayo == ens) & (t09.gen == gen) & (t09.G == mp[gt])]
            n_ok += any(cerca(v, [r[c] for c in COLS[ens]]) for _, r in s.iterrows())
        sols.append((n_ok, mp))
    mejor = max(s[0] for s in sols)
    top = [mp for n_, mp in sols if n_ == mejor]
    assert len(top) == 1, (ens, mejor, top[:3])
    CLAVE[ens] = top[0]
    COMPROB[ens] = dict(filas_publicadas=len(filas), filas_reproducidas=mejor)
# verificación independiente: RQ de T08 contra graficos.xlsx (valores por réplica de las figuras originales)
xg = pd.ExcelFile(os.path.join(D2, "graficos.xlsx"))


def bloques(hoja, ens):
    d = xg.parse(hoja, header=None)
    out, cur = {}, None
    for _, r in d.iterrows():
        c0 = str(r[0]).strip() if not pd.isna(r[0]) else ""
        if re.match(r"^G[CD]", c0):
            cur = "GDC" if c0 == "GDS" and ens else c0; out[cur] = []
        elif cur is not None and not pd.isna(r[1]) and len(out[cur]) < 9:
            out[cur].append(float(r[1]))
    assert set(out) == set(GENOS[ens]) and all(len(v) == 9 for v in out.values()), (hoja, {k: len(v) for k, v in out.items()})
    return {g: [np.array(v[0:3]), np.array(v[3:6]), np.array(v[6:9])] for g, v in out.items()}


# graficos.xlsx usa como calibrador una réplica del control (vale 1) y T08 la media del control: los dos juegos son
# el mismo dato si, para cada gen y genotipo, log2(graficos / T08) es constante en las 9 réplicas.
coinc, CONF = {}, {e: set() for e in GENOS}
for ens in ("sequia", "frio"):
    n_ok = n_tot = 0; difieren = []
    for viejo, hg in HOJA_GRAF.items():
        gb = bloques(hg + ("-C" if ens == "frio" else ""), ens)
        for gt in GENOS[ens]:
            a = np.concatenate([np.sort(t08[(t08.ensayo == ens) & (t08.gen == viejo) & (t08.G == CLAVE[ens][gt]) & (t08.trat == tr)].RQ.values)
                                for tr in TRAT[ens]])
            b = np.concatenate([np.sort(g) for g in gb[gt]])
            n_tot += 1
            igual = len(a) == len(b) == 9 and int(np.sum(np.abs(np.log2(b / a) - np.median(np.log2(b / a))) < 0.05)) >= 8
            n_ok += igual
            if igual:
                CONF[ens].add(gt)
            else:
                difieren.append(f"{VIEJO[viejo]}/{gt}")
    coinc[ens] = dict(genes_x_genotipos=n_tot, mismo_dato=n_ok, difieren=difieren)
# La clave de genotipos se fija con las p publicadas (arriba) y queda confirmada por EglNAC33 y EglNAC115, que coinciden
# réplica a réplica en los 5 genotipos de los dos ensayos. Las celdas gen x genotipo en que los Cq crudos NO dan los valores
# de las figuras originales se declaran aquí una a una: cualquier discrepancia nueva, o una declarada que desaparezca, para
# el guion. La figura nueva usa los Cq crudos (como declara Métodos) y la leyenda informa de estas celdas.
DECLARADAS = {"sequia": ["EglNAC66/GDT1", "EglNAC66/GDC"],
              "frio": ["EglNAC66/GCT1", "EglNAC66/GCT2", "EglNAC66/GCS1", "EglNAC66/GCS2", "EglNAC66/GDC",
                       "EglNAC85/GCT1", "EglNAC85/GCT2", "EglNAC85/GCS1", "EglNAC85/GCS2"]}
for e in GENOS:
    assert sorted(coinc[e]["difieren"]) == sorted(DECLARADAS[e]), (e, coinc[e]["difieren"])
    assert CONF[e] == set(GENOS[e]), CONF
for gen_ok in ("EglNAC33", "EglNAC115"):
    for e in GENOS:
        assert not any(x.startswith(gen_ok + "/") for x in coinc[e]["difieren"]), (gen_ok, e)
COMPROB["RQ_T08_vs_graficos"] = coinc

# --- búsqueda exhaustiva del origen de las celdas discrepantes en los ficheros de Cq crudos -----------------
CRUDOS = {"frio": [f for f in os.listdir(D2) if f.startswith("FR") and f.endswith(".xlsx")][0],
          "sequia": [f for f in os.listdir(D2) if f.startswith("DS") and f.endswith(".xlsx")][0]}
MUESTRAS = [f"{t}-R{i}" for t in ("NA", "AAH", "ADH") for i in (1, 2, 3)]


def hoja_cruda(xl, h):
    d = xl.parse(h, header=None); out, det = {}, None
    for i in range(len(d)):
        r = d.iloc[i]
        if str(r[0]).strip() == "Detector" and i + 1 < len(d):
            det = str(d.iloc[i + 1][0]).strip(); out.setdefault(det, {})
        elif det and isinstance(r[1], str) and re.match(r"^(NA|AAH|ADH)-R\d$", r[1].strip()) and not pd.isna(r[5]):
            out[det][r[1].strip()] = (float(r[5]), float(r[7]) if not pd.isna(r[7]) else np.nan)
    return out


def variantes(H):
    tg = [k for k in H if k.startswith("NAC")]
    if len(tg) != 1 or any(len(H[k]) != 9 for k in H):
        return {}
    tc = np.array([H[tg[0]][s][0] for s in MUESTRAS]); tq = np.array([H[tg[0]][s][1] for s in MUESTRAS])
    refs = [k for k in H if k != tg[0]]; v = {}
    for nm, rr in [("media de " + " y ".join(refs), refs)] + [(r, [r]) for r in refs]:
        rc = np.mean([[H[r][s][0] for s in MUESTRAS] for r in rr], axis=0)
        rq = np.exp(np.mean([np.log([H[r][s][1] for s in MUESTRAS]) for r in rr], axis=0))
        v[f"2^-ΔΔCq, {nm}"] = 2.0 ** -(tc - rc); v[f"curva estándar, {nm}"] = tq / rq
    return v


BUSQ = {}
for ens in ("sequia", "frio"):
    xl = pd.ExcelFile(os.path.join(D2, CRUDOS[ens]))
    cand = {h: variantes(hoja_cruda(xl, h)) for h in xl.sheet_names}
    cand = {h: v for h, v in cand.items() if v}
    n_var = len({k for v in cand.values() for k in v})
    for celda in coinc[ens]["difieren"]:
        gen, gt = celda.split("/")
        viejo = {v: k for k, v in VIEJO.items()}[gen]
        b_ = bloques(HOJA_GRAF[viejo] + ("-C" if ens == "frio" else ""), ens)[gt]
        b_ = np.concatenate([np.sort(g) for g in b_])
        mejor = min((float(np.max(np.abs(np.log2(b_ / np.concatenate([np.sort(x[3 * i:3 * i + 3]) for i in range(3)])) -
                                       np.median(np.log2(b_ / np.concatenate([np.sort(x[3 * i:3 * i + 3]) for i in range(3)])))))), h, k)
                    for h, vs in cand.items() for k, x in vs.items())
        BUSQ[f"{ens}/{celda}"] = dict(hojas_probadas=len(cand), normalizaciones=n_var, mejor_desvio_log2=round(mejor[0], 2),
                                      mejor_hoja=mejor[1], mejor_normalizacion=mejor[2], reproduce=mejor[0] < 0.05)
COMPROB["busqueda_en_Cq_crudos"] = BUSQ

# --- 2. datos y estadística -------------------------------------------------------------------------------
REG = []


def registrar(fig, ens, gen, gt, grupos):
    p_kw, dn = kw_dunn(grupos)
    med = [float(np.mean(g)) for g in grupos]
    se = [float(np.std(g, ddof=1) / np.sqrt(len(g))) for g in grupos]
    le = letras(med, dn) if p_kw <= 0.05 else ["a"] * len(grupos)
    REG.append(dict(figura=fig, ensayo=ens, gen=gen, genotipo=gt, p_KW=p_kw,
                    **{f"Dunn_{TRAT[ens][a]}-{TRAT[ens][b]}": p for (a, b), p in dn.items()},
                    **{f"n_{t}": len(g) for t, g in zip(TRAT[ens], grupos)},
                    **{f"media_{t}": m for t, m in zip(TRAT[ens], med)},
                    **{f"EE_{t}": s for t, s in zip(TRAT[ens], se)},
                    **{f"letra_{t}": l for t, l in zip(TRAT[ens], le)}, KW_significativo=bool(p_kw <= 0.05)))
    return med, se, le


DATOS = {"F21": {}, "F22": {}}
COPIAS = []
for ens in ("sequia", "frio"):
    for viejo, nuevo in VIEJO.items():
        for gt in GENOS[ens]:
            sub = t08[(t08.ensayo == ens) & (t08.gen == viejo) & (t08.G == CLAVE[ens][gt])]
            grupos = [sub[sub.trat == tr].RQ.values for tr in TRAT[ens]]
            DATOS["F21"][(ens, nuevo, gt)] = (grupos, *registrar("F21", ens, nuevo, gt, grupos))
    for b in BLANCOS:
        gb = bloques(b + ("-c" if ens == "frio" else "") + (" " if b == "ZHD11" and ens == "sequia" else ""), ens)
        for gt in GENOS[ens]:
            copia = [f"{VIEJO[v]}" for v, hg in HOJA_GRAF.items()
                     if all(np.allclose(x, y) for x, y in zip(gb[gt], bloques(hg + ("-C" if ens == "frio" else ""), ens)[gt]))]
            if copia:                                         # bloque copiado de un NAC en graficos.xlsx: no es dato del blanco
                COPIAS.append(dict(ensayo=ens, gen=b, genotipo=gt, copia_de=copia))
                continue
            cal = np.exp(np.mean(np.log(gb[gt][0])))          # mismo calibrador que T08: media geométrica del control
            grupos = [g / cal for g in gb[gt]]
            DATOS["F22"][(ens, b, gt)] = (grupos, *registrar("F22", ens, b, gt, grupos))
reg = pd.DataFrame(REG)
# reproducción de T09 (NAC) y de las tablas publicadas de los blancos
rep = 0; tot = 0
for _, r in reg[reg.figura == "F21"].iterrows():
    viejo = {v: k for k, v in VIEJO.items()}[r.gen]
    s = t09[(t09.ensayo == r.ensayo) & (t09.gen == viejo) & (t09.G == CLAVE[r.ensayo][r.genotipo])]
    mine = [r.p_KW] + [r[f"Dunn_{a}-{b}"] for a, b in itertools.combinations(TRAT[r.ensayo], 2)]
    tot += 1; rep += any(cerca(mine, [x[c] for c in COLS[r.ensayo]]) for _, x in s.iterrows())
assert rep == tot, (rep, tot)
COMPROB["NAC_reproduce_T09"] = dict(celdas=tot, iguales=rep)
pubb = {}
for ens, fx, hojas in (("sequia", "Analisis K_W_sequia.xlsx", {"zdh11": "ZHD11", "LTI65": "LTI65", "NCED3": "NCED3"}),
                       ("frio", "Analisis K_W_FRIO (1).xlsx", {"zdh11": "ZHD11", "LTI65": "LTI65"})):
    xl = pd.ExcelFile(os.path.join(D2, fx))
    for h, gen in hojas.items():
        d = xl.parse(h, header=None)
        vals = [[float(v) for v in r.values if isinstance(v, (int, float)) and not pd.isna(v)][:4] for _, r in d.iterrows()
                if any(isinstance(v, str) and re.match(r"^G[CD][TS]?\d?$|^GDS$", v.strip()) for v in r.values[:2])]
        n_ok = 0
        for i, v in enumerate(vals[:5]):
            if not ((reg.figura == "F22") & (reg.ensayo == ens) & (reg.gen == gen) & (reg.genotipo == GENOS[ens][i])).any():
                continue
            r = reg[(reg.figura == "F22") & (reg.ensayo == ens) & (reg.gen == gen) & (reg.genotipo == GENOS[ens][i])].iloc[0]
            mine = [r.p_KW] + [r[f"Dunn_{a}-{b}"] for a, b in itertools.combinations(TRAT[ens], 2)]
            n_ok += cerca(mine, v)
        pubb[f"{ens}/{gen}"] = dict(genotipos=len(vals[:5]), reproducidos=n_ok)
COMPROB["blancos_reproducen_tabla_publicada"] = pubb
reg.to_csv(os.path.join(T, "T63_qpcr_estadistica_figuras.csv"), index=False)

# --- 3. figuras -------------------------------------------------------------------------------------------
import glob as _glob
for f_ in _glob.glob((os.environ.get("CONDA_ENVS", os.path.expanduser("~/conda/envs")) + "/r-nac/lib/R/library/fontLiberation/fonts/liberation-fonts/LiberationSans-*.ttf")):
    font_manager.fontManager.addfont(f_)
disp = {f.name for f in font_manager.fontManager.ttflist}
FAM = "Arial" if "Arial" in disp else "Liberation Sans"
assert FAM in disp, "sin Arial ni Liberation Sans"
plt.rcParams.update({"font.family": FAM, "font.size": 7, "axes.linewidth": 0.6, "xtick.major.width": 0.6,
                     "ytick.major.width": 0.6, "pdf.fonttype": 42, "axes.spines.top": False, "axes.spines.right": False})
COL = {"sequia": ["#BDBDBD", "#E6A157", "#B5461F"], "frio": ["#BDBDBD", "#8FB8DE", "#2F5D8C"]}
MM = 1 / 25.4


def figura(clave, genes, nombre, cursiva):
    fig, axs = plt.subplots(2, 4, figsize=(174 * MM, 120 * MM), sharey=False)
    for fi, ens in enumerate(("sequia", "frio")):
        for ci, gen in enumerate(genes):
            ax = axs[fi, ci]
            ymax = 0
            ax._letras = []
            for gi, gt in enumerate(GENOS[ens]):
                if (ens, gen, gt) not in DATOS[clave]:
                    ax._sin = getattr(ax, "_sin", []) + [gi]
                    continue
                grupos, med, se, le = DATOS[clave][(ens, gen, gt)]
                for ti in range(3):
                    xpos = gi + (ti - 1) * 0.27
                    ax.bar(xpos, med[ti], width=0.25, color=COL[ens][ti], edgecolor="black", linewidth=0.4, zorder=2)
                    ax.errorbar(xpos, med[ti], yerr=se[ti], color="black", lw=0.6, capsize=1.2, zorder=3)
                    ax.scatter(np.full(len(grupos[ti]), xpos), grupos[ti], s=3, color="black", zorder=4, linewidths=0)
                    top = max(med[ti] + se[ti], np.max(grupos[ti]))
                    ymax = max(ymax, top)
                    ax._letras = getattr(ax, "_letras", []) + [(xpos, top, le[ti])]
            for xpos, top, l in ax._letras:
                ax.text(xpos, top + ymax * 0.03, l, ha="center", va="bottom", fontsize=6)
            for gi in getattr(ax, "_sin", []):
                ax.text(gi, ymax * 0.04, "n.a.", ha="center", va="bottom", fontsize=6, color="0.35")
            ax.set_ylim(0, ymax * 1.15)
            ax.set_xticks(range(5)); ax.set_xticklabels(GENOS[ens], rotation=45, ha="right", fontsize=6)
            ax.set_title(gen, fontsize=7, style="italic" if cursiva else "normal", pad=3)
            if ci == 0:
                ax.set_ylabel("Relative expression (2$^{-ΔΔCq}$)", fontsize=7)
        axs[fi, 0].text(-0.42, 1.12, "ab"[fi], transform=axs[fi, 0].transAxes, fontsize=9, fontweight="bold")
        hs = [plt.Rectangle((0, 0), 1, 1, fc=c, ec="black", lw=0.4) for c in COL[ens]]
        axs[fi, 3].legend(hs, TRAT[ens], fontsize=6, frameon=False, loc="upper left", bbox_to_anchor=(1.0, 1.0))
    fig.tight_layout(w_pad=0.8, h_pad=1.2)
    fig.savefig(os.path.join(FIG, f"{nombre}.png"), dpi=600)
    fig.savefig(os.path.join(FIG, f"{nombre}.pdf"))
    plt.close(fig)


figura("F21", ORDEN_NAC, "F21_qpcr_NAC", False)
figura("F22", BLANCOS, "F22_qpcr_blancos", True)
res = dict(bloques_copiados=COPIAS, clave_genotipos={e: v for e, v in CLAVE.items()}, comprobaciones=COMPROB, simbolos=VIEJO,
           hojas_sin_genotipo_publicado={e: sorted(set(t08[t08.ensayo == e].G) - set(CLAVE[e].values())) for e in CLAVE},
           n_replicas={e: sorted({int(v) for c in reg.columns if c.startswith("n_") for v in reg[reg.ensayo == e][c].dropna()}) for e in CLAVE},
           p_KW_minima=float(reg.p_KW.min()), fuente_fuente_blancos="2da entrega de datos/graficos.xlsx")
json.dump(res, open(os.path.join(OUT, "resumen_qpcr_figuras.json"), "w"), ensure_ascii=False, indent=1)
print(json.dumps(res, ensure_ascii=False, indent=1))
