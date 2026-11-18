"""Entrega 3, paso 2.2: envío de los 4 lotes de promotores a PlantCARE.

PlantCARE (https://bioinformatics.psb.ugent.be/webtools/plantcare/html/search_CARE.html)
devuelve los resultados por correo. Su página pide ficheros FASTA de texto de hasta 100 Kb
y, si se usan scripts, 60 s de espera entre envíos.

El correo se pasa como argumento en memoria y no se escribe en ningún fichero: el registro
y las respuestas guardadas lo sustituyen por <correo>.

Salidas en datos/entrega3/cis_plantcare/envios/:
  registro_envios.tsv          hora, fichero, bytes, md5, estado HTTP y longitud de la respuesta
  respuesta_lote_<n>.html      respuesta del servidor (sin el correo)
"""
import datetime, hashlib, os, re, time
import requests

URL_FORM = "https://bioinformatics.psb.ugent.be/webtools/plantcare/html/search_CARE.html"
URL_ENVIO = "https://bioinformatics.psb.ugent.be/webtools/plantcare/cgi-bin/CallMat_onCluster.htpl"
ENV = (os.environ.get("NAC_ROOT", ".") + "/ANALISIS_Eglobulus/datos/entrega3/cis_plantcare/envios")
ESPERA_S = 60
LIMITE_BYTES = 100_000


def enviar_lotes(correo, lotes=(1, 2, 3, 4), espera=ESPERA_S):
    assert correo and "@" in correo
    reg = os.path.join(ENV, "registro_envios.tsv")
    nuevo = not os.path.exists(reg)
    filas = []
    with open(reg, "a", encoding="utf8") as log:
        if nuevo:
            log.write("lote\tfichero\tbytes\tmd5\tsecuencias\thora_envio\testado_http\tbytes_respuesta\tmensaje\n")
        for k, n in enumerate(lotes):
            if k:
                time.sleep(espera)
            p = os.path.join(ENV, f"lote_{n}.fa")
            datos = open(p, "rb").read()
            assert len(datos) < LIMITE_BYTES, (p, len(datos))
            nseq = datos.count(b">")
            md5 = hashlib.md5(datos).hexdigest()
            hora = datetime.datetime.now().astimezone().isoformat(timespec="seconds")
            r = requests.post(URL_ENVIO, timeout=180, headers={"Referer": URL_FORM},
                              data={"Field_UserEmail": correo, "Field_REF": f"EglNAC_promotores_2kb_lote_{n}",
                                    "Field_Sequence": "", "MatInspector": "Search"},
                              files={"Field_File": (f"lote_{n}.fa", datos, "text/plain")})
            html = r.text.replace(correo, "<correo>")
            open(os.path.join(ENV, f"respuesta_lote_{n}.html"), "w", encoding="utf8").write(html)
            texto = re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", html)).strip()
            msg = texto[:300].replace("\t", " ")
            log.write(f"{n}\tlote_{n}.fa\t{len(datos)}\t{md5}\t{nseq}\t{hora}\t{r.status_code}\t{len(r.content)}\t{msg}\n")
            log.flush()
            filas.append(dict(lote=n, bytes=len(datos), md5=md5, secuencias=nseq, hora=hora, http=r.status_code,
                              bytes_respuesta=len(r.content), mensaje=msg))
            print(filas[-1], flush=True)
    return filas
