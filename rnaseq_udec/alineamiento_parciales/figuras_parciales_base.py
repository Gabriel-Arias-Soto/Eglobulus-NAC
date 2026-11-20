# Funciones de dibujo de las figuras por gen del chequeo de RNA-seq de los 6 NAC de dominio parcial.
# Solo definen; no dibujan. Las usa figura_tabla_suplemento.py.
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import glob
import os
import json
import re
import subprocess
import pandas as pd
import gzip
from matplotlib import font_manager
from matplotlib.patches import Rectangle, Arc

# skill:figure-style kernel.py (auto-injected on skill load)
META_GREY = "#888888"


def apply_figure_style(*, frame="open", font=None, sizes=(8, 7, 6), grid=False):
    import matplotlib as mpl
    if frame not in ("open", "boxed", "none"):
        raise ValueError(f"frame must be 'open'|'boxed'|'none', got {frame!r}")
    try:
        import os, sys, glob, matplotlib.font_manager as fm
        fdir = os.path.join(os.environ.get("CONDA_PREFIX") or sys.prefix, "fonts")
        if os.path.isdir(fdir):
            known = {f.fname for f in fm.fontManager.ttflist}
            for f in glob.glob(os.path.join(fdir, "*.ttf")):
                if f not in known:
                    fm.fontManager.addfont(f)
    except Exception:
        pass
    base, secondary, tick = sizes
    boxed = (frame == "boxed")
    rc = {
        "font.family": "sans-serif",
        "font.size": base,
        "axes.labelsize": base,
        "axes.titlesize": base,
        "legend.fontsize": secondary,
        "xtick.labelsize": tick,
        "ytick.labelsize": tick,
        "axes.linewidth": 0.6,
        "xtick.direction": "out", "ytick.direction": "out",
        "xtick.major.size": 3, "ytick.major.size": 3,
        "xtick.major.width": 0.6, "ytick.major.width": 0.6,
        "axes.spines.top": boxed, "axes.spines.right": boxed,
        "axes.spines.left": frame != "none", "axes.spines.bottom": frame != "none",
        "axes.grid": bool(grid),
        "legend.frameon": False,
        "figure.dpi": 200,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "axes.titleweight": "normal",
        "axes.titlelocation": "left",
        "axes.labelweight": "normal",
        "lines.linewidth": 1.2,
        "patch.linewidth": 0.6,
        "pdf.fonttype": 42, "ps.fonttype": 42,
    }
    if font:
        rc["font.sans-serif"] = [font, "DejaVu Sans"]
    mpl.rcParams.update(rc)


def set_frame(ax, style="open"):
    show = {"open": (False, False, True, True),
            "boxed": (True, True, True, True),
            "none": (False, False, False, False)}[style]
    for side, vis in zip(("top", "right", "bottom", "left"), show):
        ax.spines[side].set_visible(vis)
        if vis:
            ax.spines[side].set_linewidth(0.6)
    ax.tick_params(direction="out", length=0 if style == "none" else 3, width=0.6)


def panel_letter(ax, letter, dx=-0.18, dy=1.02, case="lower", fontsize=None):
    import matplotlib.pyplot as plt
    if fontsize is None:
        fontsize = plt.rcParams.get("font.size", 8) + 1
    s = letter.lower() if case == "lower" else letter.upper()
    ax.text(dx, dy, s, transform=ax.transAxes,
            fontweight="bold", fontsize=fontsize, va="bottom", ha="left")


B = (os.environ.get("NAC_ROOT", ".") + "")

for f in glob.glob((os.environ.get("CONDA_ENVS", os.path.expanduser("~/conda/envs")) + "/r-nac/lib/R/library/fontLiberation/fonts/liberation-fonts/LiberationSans-*.ttf")):
    font_manager.fontManager.addfont(f)
apply_figure_style(font="Liberation Sans", sizes=(8, 7, 6))

AL = f"{B}/rnaseq_udec/alineamiento_parciales"
R = f"{AL}/resultados"
FIG = f"{AL}/figuras"
os.makedirs(FIG, exist_ok=True)

reg = pd.read_csv(f"{AL}/regiones_5kb.bed", sep="\t", header=None, names=['seq', 'a', 'b', 'gen']).set_index('gen')
mod = pd.read_csv(f"{AL}/modelos_anotados_7loci.tsv", sep="\t")
J = pd.read_csv(f"{R}/uniones_empalme.csv")
D = pd.read_csv(f"{R}/dominio_hmmsearch.csv")
G6 = pd.read_csv(f"{R}/dominio_tramos_genomicos.csv")
V = pd.read_csv(f"{R}/veredictos_parciales.csv").set_index('gen')
RES = json.load(open(f"{R}/resumen_parciales.json"))
t01 = pd.read_csv(f"{B}/ANALISIS_Eglobulus/tablas/T01_censo_maestro_189.csv").set_index('simbolo')

ntot = sum(v['lecturas'] for v in RES['bibliotecas'].values())
nlib = len(RES['bibliotecas'])
C_ANOT, C_ARN, C_NUEVA, C_DOM, C_COB = "#4D4D4D", "#0072B2", "#D55E00", "#009E73", "#9ECAE1"


def titulo(gen):
    v = V.loc[gen]
    jj = J[(J.gen == gen) & J.fiable & ~J.anotada & J.toca_modelo & (J.hebra_motivo == J.hebra_gen)].sort_values('lecturas', ascending=False)
    if v.veredicto == "modelo incompleto":
        u = jj.iloc[0]
        return f"{gen}: un intrón no anotado de {u.intron_pb} pb ({u.lecturas} lecturas) completa el dominio NAM ({v.cob_dominio_mejor_isoforma:.0%} del perfil)"
    if v.veredicto == "compatible con truncado":
        return f"{gen}: ninguna isoforma apoyada por lecturas completa el dominio NAM (máximo {v.cob_dominio_mejor_isoforma:.0%} del perfil)"
    if v.veredicto == "sin datos suficientes":
        return f"{gen}: ninguna lectura única en los exones anotados ({nlib} bibliotecas de hoja)"
    ja = J[(J.gen == gen) & J.anotada & J.fiable].sort_values('donante_fin')
    return f"{gen} (control): las uniones anotadas ({' y '.join(map(str, ja.lecturas))} lecturas) reproducen el modelo y el dominio completo"


def dibujar(gen):
    r = reg.loc[gen]
    c = np.load(f"{R}/cobertura_{gen}.npy").sum(0)
    x = np.arange(r.a + 1, r.b + 1) / 1e6
    h = mod[mod.gen == gen].hebra.iloc[0]
    chrn = t01.loc[gen, 'cromosoma']
    fig = plt.figure(figsize=(6.85, 4.9))
    gs = fig.add_gridspec(4, 1, height_ratios=[1.25, 0.8, 1.0, 0.95], hspace=0.55)
    a1 = fig.add_subplot(gs[0])
    a2 = fig.add_subplot(gs[1], sharex=a1)
    a3 = fig.add_subplot(gs[2], sharex=a1)
    a4 = fig.add_subplot(gs[3])
    # a: estructura
    filas = [("Modelo anotado", 2.0), ("Isoforma dominante (ARN)", 1.0), ("Dominio NAM en el genoma", 0.0)]
    exa = mod[(mod.gen == gen) & (mod.tipo == "exon")]
    cda = mod[(mod.gen == gen) & (mod.tipo == "CDS")]
    a1.plot([exa.ini.min() / 1e6, exa.fin.max() / 1e6], [2, 2], color=C_ANOT, lw=0.8)
    for _, e in exa.iterrows():
        a1.add_patch(Rectangle((e.ini / 1e6, 1.85), (e.fin - e.ini + 1) / 1e6, 0.3, color=C_ANOT, alpha=0.45, lw=0))
    for _, e in cda.iterrows():
        a1.add_patch(Rectangle((e.ini / 1e6, 1.75), (e.fin - e.ini + 1) / 1e6, 0.5, color=C_ANOT, lw=0))
    dom = RES['isoformas'][gen].get('dominante')
    if dom:
        exr = dom['exones']
        a1.plot([exr[0][0] / 1e6, exr[-1][1] / 1e6], [1, 1], color=C_ARN, lw=0.8)
        for s, e in exr:
            a1.add_patch(Rectangle((s / 1e6, 0.75), (e - s + 1) / 1e6, 0.5, color=C_ARN, lw=0))
    else:
        a1.text(x[len(x) // 2], 1, "sin isoforma: no hay uniones fiables", ha="center", va="center", fontsize=7, color=C_ARN)
    g6 = G6[G6.gen == gen]
    for k, (_, t) in enumerate(g6.iterrows()):
        a1.add_patch(Rectangle((t.g_ini / 1e6, -0.2), (t.g_fin - t.g_ini + 1) / 1e6, 0.4, color=C_DOM, lw=0))
        a1.text((t.g_ini + t.g_fin) / 2e6, -0.32 - 0.0 * k, f"{t.hmm_ini}–{t.hmm_fin}", ha="center", va="top", fontsize=6, color=C_DOM)
    a1.set_yticks([f[1] for f in filas])
    a1.set_yticklabels([f[0] for f in filas])
    a1.set_ylim(-0.75, 2.5)
    for sp in ("left", "bottom"):
        a1.spines[sp].set_visible(False)
    a1.tick_params(axis="y", length=0)
    a1.tick_params(axis="x", labelbottom=False, length=0)
    a1.set_title(titulo(gen), loc="left", fontsize=8)
    a1.text(1.0, 1.02, f"hebra {'+' if h == '+' else '−'}", transform=a1.transAxes, ha="right", va="bottom", fontsize=7, color=C_ANOT)
    # b: uniones fiables (arcos), dentro de la región
    jj = J[(J.gen == gen) & J.fiable & (J.donante_fin > r.a) & (J.aceptor_ini < r.b)]
    ymax = 1
    for _, u in jj.iterrows():
        w = (u.aceptor_ini - u.donante_fin) / 1e6
        cx = (u.aceptor_ini + u.donante_fin) / 2e6
        hh = 0.35 + 0.65 * np.log10(u.lecturas) / np.log10(max(J[J.gen == gen].lecturas.max(), 10))
        col = C_ANOT if u.anotada else C_NUEVA
        a2.add_patch(Arc((cx, 0), w, 2 * hh, theta1=0, theta2=180, color=col, lw=0.6 + 1.6 * np.log10(u.lecturas) / 3))
        a2.text(cx, hh + 0.04, str(u.lecturas), ha="center", va="bottom", fontsize=6, color=col)
        ymax = max(ymax, hh)
    if len(jj) == 0:
        a2.text(x[len(x) // 2], 0.5, "ninguna unión fiable en la región", ha="center", va="center", fontsize=7, color=C_ANOT)
    a2.set_ylim(0, ymax + 0.35)
    a2.set_yticks([])
    a2.spines["left"].set_visible(False)
    a2.tick_params(axis="x", labelbottom=False)
    a2.plot([], [], color=C_ANOT, label="unión anotada")
    a2.plot([], [], color=C_NUEVA, label="unión no anotada")
    a2.legend(loc="upper right", fontsize=6, ncol=2, handlelength=1.2, borderaxespad=0)
    # c: cobertura
    a3.fill_between(x, c, step="mid", color=C_COB, lw=0)
    a3.plot(x, c, color=C_ARN, lw=0.5, drawstyle="steps-mid")
    a3.set_ylabel("Lecturas únicas")
    a3.set_xlabel(f"{chrn} (Mb)")
    if c.max() == 0:
        a3.set_ylim(0, 1)
        a3.text(x[len(x) // 2], 0.5, f"0 lecturas únicas en la región del modelo ({ntot / 1e6:.0f} M lecturas alineadas)", ha="center", va="center", fontsize=7)
    a3.set_xlim(x[0], x[-1])
    a3.margins(y=0.04)
    a3.xaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda v, _: f"{v:.3f}"))
    # d: cobertura del perfil NAM
    rows = [("Proteína anotada", "anotada"), ("Isoforma dominante (ARN)", "dominante"), ("Mejor isoforma (ARN)", "iso")]
    for k, (lab, tp) in enumerate(rows):
        y = len(rows) - 1 - k
        d = D[(D.gen == gen) & (D.tipo == tp)]
        a4.add_patch(Rectangle((1, y - 0.3), 121, 0.6, fill=False, ec="#BDBDBD", lw=0.6))
        if len(d):
            q = d.sort_values('score', ascending=False).consulta.iloc[0] if tp != "anotada" else d.consulta.iloc[0]
            for _, s in d[d.consulta == q].iterrows():
                a4.add_patch(Rectangle((s.hmm_ini, y - 0.3), s.hmm_fin - s.hmm_ini + 1, 0.6, color=C_ANOT if tp == "anotada" else C_ARN, lw=0))
            best = d[d.consulta == q].cobertura_hmm.max()
            a4.text(124, y, f"{best:.0%}", va="center", fontsize=6)
        else:
            a4.text(124, y, "—", va="center", fontsize=6)
    a4.axvline(0.99 * 122, color="#BDBDBD", lw=0)
    a4.set_xlim(0, 132)
    a4.set_ylim(-0.6, len(rows) - 0.4)
    a4.set_yticks(range(len(rows)))
    a4.set_yticklabels([r_[0] for r_ in rows][::-1])
    a4.set_xticks([1, 30, 60, 90, 122])
    a4.set_xlabel("Posición en el perfil NAM (PF02365)")
    a4.spines["left"].set_visible(False)
    a4.tick_params(axis="y", length=0)
    for ax, l in zip((a1, a2, a3, a4), "abcd"):
        panel_letter(ax, l)
    return fig


def grupos_dominio(gen):
    g = G6[G6.gen == gen].sort_values('g_ini')
    out = []
    for _, t in g.iterrows():
        key = (t.hebra, t.fase)
        if out and out[-1]['key'] == key and t.g_ini - out[-1]['g_fin'] <= 30:
            out[-1]['g_fin'] = max(out[-1]['g_fin'], t.g_fin)
            out[-1]['lab'].append(f"{t.hmm_ini}–{t.hmm_fin}")
        else:
            out.append(dict(key=key, g_ini=t.g_ini, g_fin=t.g_fin, lab=[f"{t.hmm_ini}–{t.hmm_fin}"]))
    marcos = {k: i + 1 for i, k in enumerate(dict.fromkeys(o['key'] for o in out))}
    for o in out:
        o['marco'] = marcos[o['key']]
    return out, len(marcos)


_dibujar_v1 = dibujar


def dibujar(gen):
    fig = _dibujar_v1(gen)
    a1, a2, a3, a4 = fig.axes[:4]
    r = reg.loc[gen]
    h = mod[mod.gen == gen].hebra.iloc[0]
    # a: rehacer la fila del dominio
    for p in list(a1.patches):
        if abs(p.get_y() + 0.2) < 1e-9:
            p.remove()
    for t in list(a1.texts):
        if t.get_color() == C_DOM or t.get_text().startswith("hebra"):
            t.remove()
    grp, nm = grupos_dominio(gen)
    for k, o in enumerate(grp):
        a1.add_patch(Rectangle((o['g_ini'] / 1e6, -0.2), (o['g_fin'] - o['g_ini'] + 1) / 1e6, 0.4, color=C_DOM, lw=0))
        lab = " · ".join(o['lab']) + (f" (marco {o['marco']})" if nm > 1 else "")
        y, va = (-0.3, "top") if k % 2 == 0 else (0.3, "bottom")
        a1.text((o['g_ini'] + o['g_fin']) / 2e6, y, lab, ha="center", va=va, fontsize=6, color=C_DOM)
    a1.set_ylim(-0.85, 2.5)
    # b: unión no fiable que completa un dominio (si la hay)
    v = V.loc[gen, 'dominio_completo_via_union_no_fiable']
    if isinstance(v, str) and v:
        for u in json.loads(v):
            d0, a0 = map(int, u['union'].split("-"))
            w = (a0 - d0) / 1e6
            cx = (a0 + d0) / 2e6
            a2.add_patch(Arc((cx, 0), w, 0.7, theta1=0, theta2=180, color=C_ANOT, lw=0.8, ls=(0, (2, 1.5))))
            a2.text(cx, 0.38, str(u['lecturas']), ha="center", va="bottom", fontsize=6, color=C_ANOT)
        a2.plot([], [], color=C_ANOT, ls=(0, (2, 1.5)), label="unión débil que completa un dominio")
        a2.legend(loc="upper right", fontsize=6, ncol=3, handlelength=1.6, borderaxespad=0)
    a3.set_xlabel(f"{t01.loc[gen, 'cromosoma']} (Mb), gen en la hebra {'+' if h == '+' else '−'}")
    a3.set_ylabel("Cobertura\n(lecturas únicas)")
    return fig


def dibujar2(gen):
    fig = dibujar(gen)
    a1, a2 = fig.axes[0], fig.axes[1]
    # a: un rótulo por grupo de tramos cercanos (≤600 pb), debajo
    for t in list(a1.texts):
        if t.get_color() == C_DOM:
            t.remove()
    grp, nm = grupos_dominio(gen)
    cl = []
    for o in grp:
        if cl and o['g_ini'] - cl[-1]['g_fin'] <= 600:
            cl[-1]['g_fin'] = max(cl[-1]['g_fin'], o['g_fin'])
            cl[-1]['m'].append(o)
        else:
            cl.append(dict(g_ini=o['g_ini'], g_fin=o['g_fin'], m=[o]))
    for c_ in cl:
        items = sorted(((int(l.split("–")[0]), l, o['marco']) for o in c_['m'] for l in o['lab']))
        lab = " · ".join(l + (f" (marco {m})" if nm > 1 else "") for _, l, m in items)
        a1.text((c_['g_ini'] + c_['g_fin']) / 2e6, -0.3, lab, ha="center", va="top", fontsize=6, color=C_DOM)
    # b: escalonar números de arcos cercanos y dar aire a la leyenda
    span = np.diff(a2.get_xlim())[0]
    prev = []
    for t in sorted([t for t in a2.texts if t.get_text().isdigit()], key=lambda t: t.get_position()[0]):
        x_, y_ = t.get_position()
        while any(abs(x_ - px) < 0.035 * span and abs(y_ - py) < 0.16 for px, py in prev):
            y_ += 0.17
        t.set_position((x_, y_))
        prev.append((x_, y_))
    top = max([p[1] for p in prev] + [1])
    a2.set_ylim(0, top + 0.75)
    return fig


def dibujar3(gen):
    fig = dibujar2(gen)
    a1, a2 = fig.axes[0], fig.axes[1]
    labs = [t for t in a1.texts if t.get_color() == C_DOM]
    for k, t in enumerate(sorted(labs, key=lambda t: t.get_position()[0])):
        if k % 2 == 1:
            t.set_position((t.get_position()[0], 0.3))
            t.set_va("bottom")
    span = np.diff(a2.get_xlim())[0]
    prev = []
    for t in sorted([t for t in a2.texts if t.get_text().isdigit()], key=lambda t: t.get_position()[0]):
        x_, y_ = t.get_position()
        y_ = min(y_, 1.2)
        while any(abs(x_ - px) < 0.045 * span and abs(y_ - py) < 0.16 for px, py in prev):
            y_ += 0.17
        t.set_position((x_, y_))
        prev.append((x_, y_))
    a2.set_ylim(0, max(p[1] for p in prev + [(0, 1)]) + 1.05)
    return fig


def dibujar4(gen):
    fig = dibujar3(gen)
    a2 = fig.axes[1]
    nums = [t for t in a2.texts if t.get_text().isdigit()]
    if nums:
        mx = max(int(t.get_text()) for t in nums)
        for t in nums:
            if int(t.get_text()) < 0.10 * mx and t.get_color() != C_ANOT or (int(t.get_text()) < 0.10 * mx and t.get_text() != "2"):
                t.set_visible(False)
        # la unión débil que completa un dominio conserva su número (gris, discontinua)
        for t in nums:
            if t.get_color() == C_ANOT and V.loc[gen, 'dominio_completo_via_union_no_fiable'] == V.loc[gen, 'dominio_completo_via_union_no_fiable'] and t.get_text() == "2":
                t.set_visible(True)
        fig.canvas.draw()
        rr = fig.canvas.get_renderer()
        placed = []
        for t in sorted([t for t in nums if t.get_visible()], key=lambda t: t.get_position()[0]):
            for _ in range(10):
                bb = t.get_window_extent(rr)
                if not any(bb.overlaps(p) for p in placed):
                    break
                x_, y_ = t.get_position()
                t.set_position((x_, y_ + 0.14))
            placed.append(t.get_window_extent(rr))
        top = max(t.get_position()[1] for t in nums if t.get_visible())
        a2.set_ylim(0, top + 1.05)
    return fig


