"""Entrega 3, paso 3.2 (control): recálculo hipergeométrico independiente del GO.

Rehace en Python, sin clusterProfiler, las pruebas de e3_05_go_enriquecimiento.R a partir de
la anotación con ancestros que exportó el script de R (anotacion_GOALL_universo.tsv.gz).

  N = genes del universo con alguna anotación en la ontología
  M = genes del universo anotados al término (se prueban los términos con 10 <= M <= 500
      y al menos un gen de la lista; en el contraste por genes de E. globulus se usan los
      mismos términos)
  n = genes de la lista con alguna anotación en la ontología
  k = genes de la lista anotados al término
  p = P(X >= k), X ~ hipergeométrica(N, M, n); BH sobre los términos probados

Criterio: mismos términos y diferencia relativa < 1e-6 en p y en p.adjust.
Salida: datos/entrega3/go/control_recalculo_go.json
"""
import json, os
import numpy as np
import pandas as pd
from scipy.stats import hypergeom

A = (os.environ.get("NAC_ROOT", ".") + "/ANALISIS_Eglobulus")
GO = f"{A}/datos/entrega3/go"
MIN_GS, MAX_GS, TOL = 10, 500, 1e-6


def bh(p):
    p = np.asarray(p, float)
    o = np.argsort(p)[::-1]
    r = np.arange(len(p), 0, -1)
    q = np.minimum.accumulate(p[o] * len(p) / r)
    out = np.empty_like(q)
    out[o] = np.minimum(q, 1)
    return out


def pruebas(ann, lista, terminos=None):
    """ann: DataFrame (unidad, GOALL) de una ontología; lista: unidades de la lista."""
    N = ann.unidad.nunique()
    lista = set(lista) & set(ann.unidad)
    n = len(lista)
    M = ann.groupby("GOALL").unidad.nunique()
    k = ann[ann.unidad.isin(lista)].groupby("GOALL").unidad.nunique()
    if terminos is None:
        terminos = sorted(set(k.index) & set(M.index[(M >= MIN_GS) & (M <= MAX_GS)]))
    k = k.reindex(terminos).fillna(0).astype(int)
    M = M.reindex(terminos).astype(int)
    p = hypergeom.sf(k.values - 1, N, M.values, n)
    return pd.DataFrame(dict(ID=terminos, k=k.values, M=M.values, n=n, N=N, pvalue=p, p_adj=bh(p)))


def comparar(py, r):
    m = py.merge(r[["ID", "pvalue", "p.adjust", "Count"]], on="ID", how="outer", suffixes=("_py", "_R"), indicator=True)
    mismos = bool((m._merge == "both").all())
    dp = (m.pvalue_py - m.pvalue_R).abs() / m.pvalue_R.clip(lower=1e-300)
    da = (m.p_adj - m["p.adjust"]).abs() / m["p.adjust"].clip(lower=1e-300)
    return dict(terminos_py=int(len(py)), terminos_R=int(len(r)), mismos_terminos=mismos,
                k_igual=bool((m.k == m.Count).all()), dif_rel_max_p=float(dp.max()), dif_rel_max_padj=float(da.max()),
                pasa=bool(mismos and (m.k == m.Count).all() and dp.max() < TOL and da.max() < TOL))


def main():
    ann = pd.read_csv(f"{GO}/anotacion_GOALL_universo.tsv.gz", sep="\t")
    t41 = pd.read_csv(f"{A}/tablas/T41_homologos_Ath_NAC189.csv", keep_default_na=False, na_values=[""])
    mej = pd.read_csv(f"{GO}/Egl_mejor_homologo_Ath.tsv", sep="\t")
    ok = t41.pasa_filtros.astype(str).eq("True")
    lista_agi, lista_egl = set(t41.loc[ok, "AGI"]), set(t41.loc[ok, "locus"])
    res = {}
    for ont in ("BP", "MF", "CC"):
        a = ann[ann.ONTOLOGIA == ont][["TAIR", "GOALL"]].drop_duplicates().rename(columns={"TAIR": "unidad"})
        r_agi = pd.read_csv(f"{GO}/enrichGO_{ont}_todos.tsv", sep="\t")
        py_agi = pruebas(a, lista_agi)
        a_egl = a.merge(mej[["gen_Egl", "AGI"]], left_on="unidad", right_on="AGI")[["gen_Egl", "GOALL"]]
        a_egl = a_egl.rename(columns={"gen_Egl": "unidad"}).drop_duplicates()
        r_egl = pd.read_csv(f"{GO}/enricherEgl_{ont}_todos.tsv", sep="\t")
        py_egl = pruebas(a_egl, lista_egl, terminos=sorted(py_agi.ID))
        res[ont] = dict(AGI=comparar(py_agi, r_agi), Egl=comparar(py_egl, r_egl))
        print(ont, json.dumps(res[ont]))
    res["pasa_todo"] = all(v[u]["pasa"] for v in res.values() if isinstance(v, dict) for u in ("AGI", "Egl"))
    json.dump(res, open(f"{GO}/control_recalculo_go.json", "w"), indent=1)
    print("pasa_todo:", res["pasa_todo"])


if __name__ == "__main__":
    main()
