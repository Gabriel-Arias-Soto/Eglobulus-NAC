# 10_auditoria.R -> T21_auditoria_figuras_R.csv
# Comprueba que el juego regenerado cumple las reglas medibles: ancho canonico,
# alto dentro de la pagina y existencia de la version vectorial.
if (!exists("ANCHO_COLUMNA")) source("scripts/00_estilo.R")
suppressPackageStartupMessages(library(png))

# El patron anterior pedia digitos seguidos de guion bajo, asi que F02S_ no
# entraba y la suplementaria no se auditaba: el recuento decia 12 figuras
# omitiendo una. La S opcional la incluye, y la columna suplementaria evita que
# su desborde de alto (declarado, ver 07_expresion.R) se lea como defecto.
pngs <- sort(list.files(FIGURAS, pattern = "^[FS][0-9]+S?_.*\\.png$", full.names = TRUE))

# No todas las figuras se guardan a la misma resolucion: las de R salen a DPI y
# la de estructura (matplotlib sobre renders de PyMOL) a 400. Dividir todo por
# DPI daba F10 por sobreancha cuando en realidad esta al ancho canonico, asi que
# la resolucion se lee del propio PNG (trozo pHYs, pixeles por metro).
dpi_png <- function(f) {
  r   <- readBin(f, "raw", n = 4096)
  sig <- charToRaw("pHYs")
  p   <- which(r[1:(length(r) - 3)] == sig[1])
  p   <- p[vapply(p, function(k) all(r[k:(k + 3)] == sig), logical(1))]
  if (!length(p) || r[p[1] + 12] != as.raw(1)) return(DPI)   # sin unidad fisica
  ppm <- sum(as.integer(r[(p[1] + 4):(p[1] + 7)]) * 256^(3:0))
  if (ppm <= 0) DPI else round(ppm * 0.0254)
}

a <- lapply(pngs, function(f) {
  d <- dim(png::readPNG(f))
  dp <- dpi_png(f)
  data.frame(figura = sub("\\.png$", "", basename(f)),
             px_ancho = d[2], px_alto = d[1], dpi = dp,
             ancho_pulg = round(d[2] / dp, 2), alto_pulg = round(d[1] / dp, 2),
             pdf = file.exists(sub("\\.png$", ".pdf", f)))
}) |> bind_rows() |>
  mutate(suplementaria = grepl("^F[0-9]+S_|^S[0-9]+_", figura),
         ancho_canonico = abs(ancho_pulg - ANCHO_COLUMNA) < 0.06 |
                          abs(ancho_pulg - ANCHO_MEDIA) < 0.06,
         cabe_en_pagina = alto_pulg <= ALTO_MAXIMO + 0.05,
         escala_al_cuerpo = round(ANCHO_COLUMNA / ancho_pulg, 2))
write_csv(a, file.path(TABLAS, "T21_auditoria_figuras_R.csv"))
print(a, row.names = FALSE)
cuerpo <- a[!a$suplementaria, ]
message(sprintf("\n%d figuras (%d de cuerpo + %d suplementarias) | ancho canonico: %d | con PDF: %d",
                nrow(a), nrow(cuerpo), sum(a$suplementaria),
                sum(a$ancho_canonico), sum(a$pdf)))
message(sprintf("cuerpo que cabe en pagina: %d de %d",
                sum(cuerpo$cabe_en_pagina), nrow(cuerpo)))
if (any(!cuerpo$cabe_en_pagina))
  message("desbordan: ", paste(cuerpo$figura[!cuerpo$cabe_en_pagina], collapse = ", "))

# --- conformidad con la norma de Trees -> T35 --------------------------------
# Comprueba requisito a requisito lo que se puede leer del propio fichero. Lo
# que NO se comprueba aqui, porque no esta en el PNG, queda declarado en el
# Metodos: tipografia Arial/Helvetica y cuerpo de 8-12 pt (los fija 00_estilo.R)
# y ausencia de titulo dentro de la imagen (regla 2 de la cabecera de estilo).
ihdr <- function(f) {                     # profundidad y tipo de color del PNG
  r <- readBin(f, "raw", n = 26)
  list(bits = as.integer(r[25]), tipo = as.integer(r[26]))
}
TIPO_PNG <- c("0" = "gris", "2" = "RGB", "3" = "paleta",
              "4" = "gris+alfa", "6" = "RGBA")

tr <- a |>
  rowwise() |>
  mutate(h = list(ihdr(file.path(FIGURAS, paste0(figura, ".png"))))) |>
  ungroup() |>
  # Los mm se calculan desde los pixeles, no desde ancho_pulg: esa columna ya
  # viene redondeada a dos decimales y encadenar redondeos daba 84,1 mm para
  # figuras que miden 84,0. La tolerancia lo absorbia, pero la cifra publicada
  # tiene que ser la de verdad.
  mutate(ancho_mm = round(px_ancho / dpi * 25.4, 1),
         alto_mm  = round(px_alto  / dpi * 25.4, 1),
         bits     = vapply(h, `[[`, integer(1), "bits"),
         color    = TIPO_PNG[as.character(vapply(h, `[[`, integer(1), "tipo"))],
         # Trees admite exactamente dos anchos: 84 mm (columna) y 174 mm (caja).
         ancho_ok = abs(ancho_mm - 84) <= 0.5 | abs(ancho_mm - 174) <= 0.5,
         alto_ok  = alto_mm <= 234.5,
         # arte combinado (trazo + rotulacion): minimo 600 ppp
         dpi_ok   = dpi >= 600,
         # RGB de 8 bits por canal; el alfa sobra y estorba en preimpresion
         color_ok = bits == 8 & color == "RGB",
         vector_ok = pdf,
         conforme = ancho_ok & alto_ok & dpi_ok & color_ok) |>
  select(figura, suplementaria, ancho_mm, alto_mm, dpi, bits, color,
         ancho_ok, alto_ok, dpi_ok, color_ok, vector_ok, conforme)

write_csv(tr, file.path(TABLAS, "T35_conformidad_Trees.csv"))
print(as.data.frame(tr[, c("figura", "ancho_mm", "alto_mm", "dpi", "color", "conforme")]),
      row.names = FALSE)

fallan <- tr[!tr$conforme, ]
message(sprintf("\nconformidad Trees: %d de %d figuras", sum(tr$conforme), nrow(tr)))
if (nrow(fallan)) {
  for (i in seq_len(nrow(fallan))) {
    r <- fallan[i, ]
    motivo <- c("ancho", "alto", "ppp", "color")[!c(r$ancho_ok, r$alto_ok,
                                                    r$dpi_ok, r$color_ok)]
    message(sprintf("  %-34s incumple: %s", r$figura, paste(motivo, collapse = ", ")))
  }
  # Una suplementaria alta es admisible (se imprime aparte); una de cuerpo no.
  if (any(!fallan$suplementaria))
    stop("figuras de cuerpo fuera de la norma de Trees: ",
         paste(fallan$figura[!fallan$suplementaria], collapse = ", "), call. = FALSE)
}
