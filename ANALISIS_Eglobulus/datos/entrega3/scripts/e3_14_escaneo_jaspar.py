"""Entrega 3, análisis cis local (JASPAR 2026 + FIMO), paso 1: escaneo de todos los promotores.

Se buscan las matrices de JASPAR 2026 CORE plants non-redundant en los 35.943 promotores de 2 kb
(189 NAC + fondo genómico) con FIMO 5.5.9: --text --thresh 1e-4 --max-strand y fondo de orden 0
calculado de todos los promotores con fasta-get-markov (ambas hebras combinadas). Los promotores se
reparten en N_TROZOS trozos; cada trozo se escanea desde una copia en /tmp (RAM) y sus sitios se
guardan en Z en formato compacto (npz). Para los 189 NAC se guarda además la línea completa de FIMO
con la secuencia del sitio. Nada grande se escribe en el disco de WSL.

Uso (desde la raíz «Revisión NAC», con el python de nac3):
  python e3_14_escaneo_jaspar.py preparar              tabla de motivos, trozos, fondo de Markov
  python e3_14_escaneo_jaspar.py lote [--n 8] [--procesos 8]   escanea los siguientes n trozos pendientes
  python e3_14_escaneo_jaspar.py cerrar                matrices de recuentos, sitios NAC, controles

Controles:
  preparar: MEME y TRANSFAC con los mismos motivos y anchos; filas de probabilidad que suman 1 y
            coinciden con los recuentos de TRANSFAC; trozos = FASTA original (mismas secuencias, en
            el mismo orden); 189 NAC de T40 = es_NAC de la tabla de promotores.
  lote:     cada sitio con ancho = ancho del motivo, dentro del promotor y con p <= umbral.
  cerrar:   los 64 trozos presentes, cada sitio en su trozo; los sitios NAC del escaneo por trozos
            son idénticos a un escaneo aparte del FASTA de los 189 NAC con los mismos parámetros.
Entradas:
  datos/entrega3/cis_jaspar/matrices/JASPAR2026_CORE_plants_non-redundant_pfms_{meme,transfac}.txt
  datos/entrega3/promotores/promotores_2kb_todos.{fa,tsv}, NAC189_promotores_2kb.fa
  tablas/T40_promotores_NAC_2kb.csv (símbolo <-> locus)
Salidas (datos/entrega3/cis_jaspar/):
  motivos_JASPAR2026.tsv, loci_orden.txt, fondo_markov0.txt, trozos/trozo_KK.fa,
  sitios/trozo_KK.npz (seq, motivo, inicio, hebra, score, p), sitios/trozo_KK_NAC.tsv,
  conteos_JASPAR2026.npz (sitios por promotor y motivo con p<1e-4 y p<1e-5),
  sitios_NAC189_JASPAR2026.tsv, resumen_escaneo_JASPAR2026.json
"""
import argparse, hashlib, json, multiprocessing as mp, os, platform, shutil, subprocess, sys, tempfile, time
from array import array
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np
import pandas as pd

BASE = (os.environ.get("NAC_ROOT", ".") + "")
A = os.path.join(BASE, "ANALISIS_Eglobulus")
E3 = os.path.join(A, "datos", "entrega3")
PROM = os.path.join(E3, "promotores")
CJ = os.path.join(E3, "cis_jaspar")
MAT = os.path.join(CJ, "matrices")
TROZOS = os.path.join(CJ, "trozos")
SITIOS = os.path.join(CJ, "sitios")
MEME = os.path.join(MAT, "JASPAR2026_CORE_plants_non-redundant_pfms_meme.txt")
TRANSFAC = os.path.join(MAT, "JASPAR2026_CORE_plants_non-redundant_pfms_transfac.txt")
FA_TODOS = os.path.join(PROM, "promotores_2kb_todos.fa")
FA_NAC = os.path.join(PROM, "NAC189_promotores_2kb.fa")
BG = os.path.join(CJ, "fondo_markov0.txt")
TAB_MOT = os.path.join(CJ, "motivos_JASPAR2026.tsv")
LOCI = os.path.join(CJ, "loci_orden.txt")
MEMEBIN = (os.environ.get("CONDA_ENVS", os.path.expanduser("~/conda/envs")) + "/rnaseq/bin")
FIMO = os.path.join(MEMEBIN, "fimo")
FGM = os.path.join(MEMEBIN, "fasta-get-markov")
N_TROZOS = 64
UMBRAL = 1e-4
UMBRAL_ESTRICTO = 1e-5
OPC_FIMO = ["--text", "--thresh", str(UMBRAL), "--max-strand", "--verbosity", "1"]


def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def leer_fasta(p):
    seqs, nombre, partes = {}, None, []
    with open(p) as fh:
        for linea in fh:
            linea = linea.strip()
            if linea.startswith(">"):
                if nombre is not None:
                    seqs[nombre] = "".join(partes)
                nombre, partes = linea[1:].split()[0], []
            elif linea:
                partes.append(linea)
    if nombre is not None:
        seqs[nombre] = "".join(partes)
    return seqs


def leer_meme(p):
    mot, actual = [], None
    with open(p) as fh:
        for linea in fh:
            t = linea.split()
            if not t:
                continue
            if t[0] == "MOTIF":
                actual = {"motivo_id": t[1], "nombre": " ".join(t[2:]), "filas": []}
                mot.append(actual)
            elif t[0] == "letter-probability":
                actual["ancho"] = int(t[t.index("w=") + 1])
            elif actual is not None and len(t) == 4 and all(x.replace(".", "", 1).replace("e-", "", 1).isdigit() for x in t):
                actual["filas"].append([float(x) for x in t])
    return mot


def leer_transfac(p):
    mot, actual = [], None
    with open(p) as fh:
        for linea in fh:
            linea = linea.rstrip("\n")
            if linea.startswith("AC "):
                actual = {"motivo_id": linea[3:].strip(), "filas": [], "cc": {}}
            elif linea.startswith("ID "):
                actual["nombre"] = linea[3:].strip()
            elif linea[:2].isdigit() and actual is not None:
                actual["filas"].append([float(x) for x in linea.split()[1:5]])
            elif linea.startswith("CC "):
                k, _, v = linea[3:].partition(":")
                actual["cc"][k.strip()] = v.strip()
            elif linea.startswith("//") and actual is not None:
                mot.append(actual)
                actual = None
    return mot


def preparar():
    t0 = time.time()
    os.makedirs(TROZOS, exist_ok=True)
    os.makedirs(SITIOS, exist_ok=True)
    # Motivos: MEME (lo que lee FIMO) frente a TRANSFAC (recuentos y anotación)
    mm, tf = leer_meme(MEME), leer_transfac(TRANSFAC)
    ids = [m["motivo_id"] for m in mm]
    assert len(ids) == len(set(ids)) == len(tf), f"MEME {len(ids)} / TRANSFAC {len(tf)}"
    tfd = {m["motivo_id"]: m for m in tf}
    assert set(ids) == set(tfd), "identificadores distintos entre MEME y TRANSFAC"
    dif_max, suma_max = 0.0, 0.0
    filas = []
    for i, m in enumerate(mm):
        t = tfd[m["motivo_id"]]
        P = np.array(m["filas"])
        C = np.array(t["filas"])
        assert P.shape == (m["ancho"], 4) == C.shape, f"{m['motivo_id']}: anchos {P.shape} {C.shape}"
        suma_max = max(suma_max, np.abs(P.sum(1) - 1).max())
        dif_max = max(dif_max, np.abs(P - C / C.sum(1, keepdims=True)).max())
        assert m["nombre"] == t["nombre"], f"{m['motivo_id']}: nombre {m['nombre']} / {t['nombre']}"
        cc = t["cc"]
        fam = cc.get("tf_family", "")
        filas.append({"indice": i, "motivo_id": m["motivo_id"], "nombre": m["nombre"], "ancho": m["ancho"],
                      "familia_JASPAR": fam, "clase_JASPAR": cc.get("tf_class", ""),
                      "familia": fam if fam else cc.get("tf_class", ""),
                      "tipo_dato": cc.get("data_type", ""), "uniprot": cc.get("uniprot_ids", ""),
                      "pubmed": cc.get("pubmed_ids", ""), "grupo_tax": cc.get("tax_group", "")})
    assert suma_max < 1e-3 and dif_max < 1e-3, f"probabilidades: suma {suma_max:.2e}, dif. con TRANSFAC {dif_max:.2e}"
    mot = pd.DataFrame(filas)
    assert (mot.grupo_tax == "plants").all()
    mot.to_csv(TAB_MOT, sep="\t", index=False)
    print(f"motivos: {len(mot)} (MEME = TRANSFAC; |suma-1| máx {suma_max:.1e}; dif. máx {dif_max:.1e}); "
          f"ancho {mot.ancho.min()}-{mot.ancho.max()}; sin familia JASPAR {int((mot.familia_JASPAR == '').sum())} "
          f"(se usa la clase)")

    # Promotores y NAC
    gen = leer_fasta(FA_TODOS)
    loci = list(gen)
    tab = pd.read_csv(os.path.join(PROM, "promotores_2kb_todos.tsv"), sep="\t")
    t40 = pd.read_csv(os.path.join(A, "tablas", "T40_promotores_NAC_2kb.csv"))
    assert set(tab.locus) == set(loci) and tab.locus.is_unique
    assert set(t40.locus) == set(tab.locus[tab.es_NAC]) and len(t40) == 189, "NAC de T40 != es_NAC"
    assert set("".join(gen.values())) <= set("ACGT"), "bases distintas de ACGT"
    with open(LOCI, "w") as fh:
        fh.write("\n".join(loci) + "\n")
    # Trozos contiguos
    partes = np.array_split(np.arange(len(loci)), N_TROZOS)
    n_esc = 0
    for k, P in enumerate(partes, 1):
        with open(os.path.join(TROZOS, f"trozo_{k:02d}.fa"), "w") as fh:
            for j in P:
                fh.write(f">{loci[j]}\n{gen[loci[j]]}\n")
                n_esc += 1
    rel = [leer_fasta(os.path.join(TROZOS, f"trozo_{k:02d}.fa")) for k in range(1, N_TROZOS + 1)]
    juntos = [x for d in rel for x in d.items()]
    assert [x[0] for x in juntos] == loci and all(gen[n] == s for n, s in juntos), "trozos != FASTA original"
    print(f"promotores: {len(loci)} ({len(t40)} NAC); trozos: {N_TROZOS} de {min(map(len, partes))}-"
          f"{max(map(len, partes))} promotores; releídos = FASTA original")

    # Fondo de orden 0 de todos los promotores (hebras combinadas)
    r = subprocess.run([FGM, "-m", "0", "-dna", FA_TODOS, BG], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr[-500:]
    fr = {}
    for linea in open(BG):
        t = linea.split()
        if len(t) == 2 and t[0] in "ACGT":
            fr[t[0]] = float(t[1])
    print(f"fondo de orden 0: {fr}")
    v = subprocess.run([FIMO, "--version"], capture_output=True, text=True).stdout.strip()
    res = {"fecha": time.strftime("%Y-%m-%d %H:%M:%S"), "fimo": v, "opciones_fimo": OPC_FIMO + ["--bfile", "fondo_markov0.txt"],
           "jaspar": "JASPAR 2026 (release 11) CORE plants non-redundant", "n_motivos": len(mot),
           "md5": {os.path.basename(p): md5(p) for p in (MEME, TRANSFAC, FA_TODOS, FA_NAC, BG)},
           "fecha_descarga_jaspar": open(os.path.join(MAT, "fecha_descarga.txt")).read().strip(),
           "n_promotores": len(loci), "n_NAC": len(t40), "n_trozos": N_TROZOS, "fondo_markov0": fr,
           "python": platform.python_version(), "numpy": np.__version__, "pandas": pd.__version__}
    json.dump(res, open(os.path.join(CJ, "preparacion_escaneo.json"), "w"), ensure_ascii=False, indent=1)
    print(f"preparación lista ({time.time() - t0:.0f} s)")


# ---------- escaneo ----------
_M = _S = _ANCHO = _NAC = None


def _cargar_indices():
    global _M, _S, _ANCHO, _NAC
    mot = pd.read_csv(TAB_MOT, sep="\t")
    _M = dict(zip(mot.motivo_id, mot.indice))
    _ANCHO = mot.ancho.to_numpy()
    loci = open(LOCI).read().split()
    _S = {x: i for i, x in enumerate(loci)}
    t40 = pd.read_csv(os.path.join(A, "tablas", "T40_promotores_NAC_2kb.csv"))
    _NAC = {_S[x] for x in t40.locus}


def escanear(k, fa_origen=None, salida=True):
    """Escanea un trozo (o un FASTA dado). Devuelve estadísticas; si salida, guarda npz y sitios NAC."""
    t0 = time.time()
    fa_z = fa_origen or os.path.join(TROZOS, f"trozo_{k:02d}.fa")
    a_seq, a_mot, a_ini, a_heb, a_sc, a_p = array("i"), array("h"), array("h"), array("b"), array("f"), array("f")
    nac, malos, nombres_raros = [], 0, set()
    with tempfile.TemporaryDirectory(dir="/tmp") as td:
        fa = shutil.copy(fa_z, os.path.join(td, "s.fa"))
        meme = shutil.copy(MEME, os.path.join(td, "m.meme"))
        bg = shutil.copy(BG, os.path.join(td, "bg.txt"))
        err = os.path.join(td, "err.txt")
        with open(err, "w") as ferr:
            pr = subprocess.Popen([FIMO] + OPC_FIMO + ["--bfile", bg, meme, fa], stdout=subprocess.PIPE,
                                  stderr=ferr, text=True, bufsize=1 << 20)
            for linea in pr.stdout:
                if linea.startswith("motif_id") or linea.startswith("#"):
                    continue
                f = linea.rstrip("\n").split("\t")
                mi = _M[f[0]]
                si = _S.get(f[2])
                if si is None:
                    nombres_raros.add(f[2])
                    continue
                ini, fin, pv = int(f[3]), int(f[4]), float(f[7])
                if fin - ini + 1 != _ANCHO[mi] or ini < 1 or fin > 2000 or pv > UMBRAL:
                    malos += 1
                a_seq.append(si); a_mot.append(mi); a_ini.append(ini)
                a_heb.append(1 if f[5] == "+" else -1); a_sc.append(float(f[6])); a_p.append(pv)
                if si in _NAC:
                    nac.append(f)
            rc = pr.wait()
        txt_err = open(err).read()
    if rc != 0:
        raise RuntimeError(f"trozo {k}: FIMO terminó con código {rc}: {txt_err[-400:]}")
    d = {"seq": np.frombuffer(a_seq, dtype=np.int32), "motivo": np.frombuffer(a_mot, dtype=np.int16),
         "inicio": np.frombuffer(a_ini, dtype=np.int16), "hebra": np.frombuffer(a_heb, dtype=np.int8),
         "score": np.frombuffer(a_sc, dtype=np.float32), "p": np.frombuffer(a_p, dtype=np.float32)}
    est = {"trozo": k, "sitios": len(a_seq), "sitios_NAC": len(nac), "malos": malos,
           "nombres_desconocidos": len(nombres_raros), "segundos": round(time.time() - t0, 1),
           "stderr": txt_err.strip()[-300:]}
    if salida:
        tmp = os.path.join(SITIOS, f"trozo_{k:02d}.tmp.npz")
        np.savez_compressed(tmp, **d)
        col = ["motif_id", "motif_alt_id", "sequence_name", "start", "stop", "strand", "score", "p-value",
               "q-value", "matched_sequence"]
        pd.DataFrame(nac, columns=col).to_csv(os.path.join(SITIOS, f"trozo_{k:02d}_NAC.tsv"), sep="\t", index=False)
        os.replace(tmp, os.path.join(SITIOS, f"trozo_{k:02d}.npz"))
        return est
    return est, d, nac


def lote(n, procesos):
    pend = [k for k in range(1, N_TROZOS + 1) if not os.path.exists(os.path.join(SITIOS, f"trozo_{k:02d}.npz"))]
    hechos = N_TROZOS - len(pend)
    sel = pend[:n]
    if not sel:
        print("no quedan trozos pendientes")
        return
    print(f"{time.strftime('%H:%M:%S')} lote: trozos {sel[0]}-{sel[-1]} ({len(sel)}), {procesos} procesos; "
          f"hechos antes {hechos}/{N_TROZOS}", flush=True)
    t0 = time.time()
    _cargar_indices()
    with ProcessPoolExecutor(max_workers=procesos, mp_context=mp.get_context("fork")) as ex:
        fut = {ex.submit(escanear, k): k for k in sel}
        for f in as_completed(fut):
            e = f.result()
            assert e["malos"] == 0 and e["nombres_desconocidos"] == 0, e
            print(f"  trozo {e['trozo']:02d}: {e['sitios']:,} sitios ({e['sitios_NAC']:,} NAC) en {e['segundos']:.0f} s"
                  + (f" | stderr: {e['stderr']}" if e["stderr"] else ""), flush=True)
    dt = time.time() - t0
    quedan = len(pend) - len(sel)
    print(f"{time.strftime('%H:%M:%S')} lote terminado en {dt / 60:.1f} min; hechos {N_TROZOS - quedan}/{N_TROZOS}; "
          f"quedan {quedan} trozos (~{quedan / len(sel) * dt / 60:.0f} min a este ritmo)", flush=True)


def cerrar():
    t0 = time.time()
    _cargar_indices()
    mot = pd.read_csv(TAB_MOT, sep="\t")
    loci = open(LOCI).read().split()
    nS, nM = len(loci), len(mot)
    partes = np.array_split(np.arange(nS), N_TROZOS)
    falta = [k for k in range(1, N_TROZOS + 1) if not os.path.exists(os.path.join(SITIOS, f"trozo_{k:02d}.npz"))]
    assert not falta, f"faltan trozos: {falta}"
    C4 = np.zeros(nS * nM, dtype=np.int64)
    C5 = np.zeros(nS * nM, dtype=np.int64)
    total = 0
    for k, P in enumerate(partes, 1):
        z = np.load(os.path.join(SITIOS, f"trozo_{k:02d}.npz"))
        s, m, p = z["seq"].astype(np.int64), z["motivo"].astype(np.int64), z["p"]
        assert s.min() >= P[0] and s.max() <= P[-1], f"trozo {k}: sitios fuera de su trozo"
        C4 += np.bincount(s * nM + m, minlength=nS * nM)
        e = p <= UMBRAL_ESTRICTO
        C5 += np.bincount(s[e] * nM + m[e], minlength=nS * nM)
        total += len(s)
    C4, C5 = C4.reshape(nS, nM), C5.reshape(nS, nM)
    assert C4.max() < 32767
    np.savez_compressed(os.path.join(CJ, "conteos_JASPAR2026.npz"), p1e4=C4.astype(np.int16), p1e5=C5.astype(np.int16),
                        loci=np.array(loci), motivos=mot.motivo_id.to_numpy().astype(str))
    print(f"sitios totales {total:,}; matriz {nS} x {nM}; p<1e-5: {int(C5.sum()):,}")

    # Sitios NAC del escaneo por trozos
    t40 = pd.read_csv(os.path.join(A, "tablas", "T40_promotores_NAC_2kb.csv"))
    loc2sim = dict(zip(t40.locus, t40.simbolo))
    nac = pd.concat([pd.read_csv(os.path.join(SITIOS, f"trozo_{k:02d}_NAC.tsv"), sep="\t", dtype=str, keep_default_na=False)
                     for k in range(1, N_TROZOS + 1)], ignore_index=True)
    nac.insert(3, "simbolo", nac.sequence_name.map(loc2sim))
    assert nac.simbolo.notna().all()

    # Control: escaneo aparte del FASTA de los 189 NAC (cabeceras = símbolos)
    nacfa = leer_fasta(FA_NAC)
    sim2loc = dict(zip(t40.simbolo, t40.locus))
    assert set(nacfa) == set(sim2loc)
    with tempfile.TemporaryDirectory(dir="/tmp") as td:
        fa = os.path.join(td, "nac_loci.fa")
        with open(fa, "w") as fh:
            for sname, seq in nacfa.items():
                fh.write(f">{sim2loc[sname]}\n{seq}\n")
        est, d, filas = escanear(0, fa_origen=fa, salida=False)
    # Comparación exacta de las líneas de FIMO como texto (multiconjunto). No se comparan p-valores como
    # float: el parser C de pandas y float() de Python pueden diferir en el último bit al leer el mismo texto.
    col = ["motif_id", "motif_alt_id", "sequence_name", "start", "stop", "strand", "score", "p-value", "q-value", "matched_sequence"]
    ref = pd.DataFrame(filas, columns=col).astype(str)
    campos = ["motif_id", "sequence_name", "start", "stop", "strand", "score", "p-value", "matched_sequence"]
    # Empates entre hebras: con --max-strand, si las dos hebras de una posición puntúan igual, FIMO elige la
    # hebra al azar (documentado), así que esa hebra cambia entre corridas. Pasa en sitios palindrómicos y en
    # sitios no palindrómicos de matrices simétricas (entonces la secuencia informada es el complemento
    # inverso). Se compara primero el texto exacto; las diferencias se aceptan solo si, con la misma posición,
    # puntuación y p, la secuencia de una línea es el complemento inverso de la otra (hebra opuesta).
    comp = str.maketrans("ACGTacgt", "TGCAtgca")
    ci = lambda s: s.upper().translate(comp)[::-1]
    for df in (nac, ref):
        df["sitio_palindromico"] = df.matched_sequence.map(lambda s: s.upper() == ci(s))
    def multiconj(df, canonico):
        x = df[campos].copy()
        if canonico:   # la hebra fuera; la secuencia, en su forma canónica (mínima entre ella y su c. inverso)
            x["strand"] = "±"
            x["matched_sequence"] = [min(s.upper(), ci(s)) for s in x.matched_sequence]
        return Counter(map(tuple, x.to_numpy()))
    a, b = multiconj(nac, False), multiconj(ref, False)
    solo_hebra, solo_hebra_palin = 0, 0
    if a != b:
        sa, sb = a - b, b - a
        print(f"diferencias exactas: {sum(sa.values())} líneas solo en los trozos, {sum(sb.values())} solo en el escaneo aparte")
        assert multiconj(nac, True) == multiconj(ref, True), "hay diferencias que no son empates de hebra"
        # Emparejamiento uno a uno de las diferencias: misma posición/puntuación/p, hebra opuesta, c. inverso
        clave_pos = lambda x: (x[0], x[1], x[2], x[3], x[5], x[6])
        pa, pb = {}, {}
        for x, n in sa.items():
            pa.setdefault(clave_pos(x), []).extend([x] * n)
        for x, n in sb.items():
            pb.setdefault(clave_pos(x), []).extend([x] * n)
        assert pa.keys() == pb.keys()
        for kpos in pa:
            assert len(pa[kpos]) == len(pb[kpos])
            for x, y in zip(sorted(pa[kpos]), sorted(pb[kpos], key=lambda y: ci(y[7]))):
                assert x[4] != y[4] and x[7].upper() == ci(y[7]), (x, y)
        solo_hebra = sum(sa.values())
        solo_hebra_palin = sum(n for x, n in sa.items() if x[7].upper() == ci(x[7]))
        print(f"  todas son empates de hebra: {solo_hebra} ({solo_hebra_palin} en sitios palindrómicos, "
              f"{solo_hebra - solo_hebra_palin} en matrices simétricas con secuencia en c. inverso)")
    iguales = True
    assert len(nac) == len(ref), f"sitios NAC: trozos {len(nac)} / escaneo aparte {len(ref)}"
    for c in ("start", "stop"):
        nac[c] = nac[c].astype(int)
    for c in ("score", "p-value"):
        nac[c] = nac[c].astype(float)
    print(f"control: sitios NAC por trozos = escaneo aparte de NAC189 ({len(nac):,} sitios, idénticos)")

    # Tabla de sitios NAC con posición relativa al ATG (el promotor acaba justo antes del ATG)
    w = nac.stop - nac.start + 1
    nac["pos_ini_rel_ATG"] = nac.start - 2001
    nac["pos_fin_rel_ATG"] = nac.stop - 2001
    nac = nac.merge(mot[["motivo_id", "familia", "clase_JASPAR"]], left_on="motif_id", right_on="motivo_id", how="left").drop(columns="motivo_id")
    assert (w == nac.matched_sequence.str.len()).all()
    nac.drop(columns=["q-value"]).sort_values(["simbolo", "start", "motif_id"]).to_csv(
        os.path.join(CJ, "sitios_NAC189_JASPAR2026.tsv"), sep="\t", index=False)

    es_nac = np.zeros(nS, bool)
    es_nac[[loci.index(x) for x in t40.locus]] = True
    por_prom = C4.sum(1)
    res = json.load(open(os.path.join(CJ, "preparacion_escaneo.json")))
    res.update({"cierre": time.strftime("%Y-%m-%d %H:%M:%S"), "sitios_totales": int(total),
                "sitios_p1e5": int(C5.sum()), "sitios_NAC": int(len(nac)),
                "sitios_por_promotor": {"NAC_mediana": float(np.median(por_prom[es_nac])),
                                        "fondo_mediana": float(np.median(por_prom[~es_nac])),
                                        "min": int(por_prom.min()), "max": int(por_prom.max())},
                "promotores_sin_sitios": int((por_prom == 0).sum()),
                "motivos_sin_sitios": mot.motivo_id[C4.sum(0) == 0].tolist(),
                "control_NAC_identico": bool(iguales),
                "control_NAC_nota": ("idéntico línea a línea salvo la hebra de posiciones con empate entre hebras "
                                     "(FIMO --max-strand la elige al azar); posición, puntuación y p idénticas"),
                "NAC_sitios_palindromicos": int(nac.sitio_palindromico.sum()),
                "NAC_lineas_con_hebra_distinta_entre_corridas": int(solo_hebra),
                "NAC_de_ellas_palindromicas": int(solo_hebra_palin)})
    json.dump(res, open(os.path.join(CJ, "resumen_escaneo_JASPAR2026.json"), "w"), ensure_ascii=False, indent=1)
    print(f"sitios por promotor: mediana NAC {np.median(por_prom[es_nac]):.0f}, fondo {np.median(por_prom[~es_nac]):.0f}; "
          f"motivos sin sitios {int((C4.sum(0) == 0).sum())}; cierre en {time.time() - t0:.0f} s")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("accion", choices=["preparar", "lote", "cerrar"])
    ap.add_argument("--n", type=int, default=8)
    ap.add_argument("--procesos", type=int, default=8)
    a = ap.parse_args()
    {"preparar": preparar, "cerrar": cerrar}.get(a.accion, lambda: lote(a.n, a.procesos))()
