# 09_estructura.R -> figura de estructura del monomero (y del dimero si existe)
# Las tres cifras de la figura anterior se verificaron exactas contra T17. Lo que
# cambia es la forma: el trazo pasa a ir coloreado con la misma escala pLDDT que
# el render, para que los dos paneles se lean como uno, y los limites del dominio
# salen de T16 en vez de estar teclados.
if (!exists("ANCHO_COLUMNA")) source("scripts/00_estilo.R")
suppressPackageStartupMessages({library(patchwork); library(png); library(grid)})

UMBRAL_DESORDEN <- 50
# escala pLDDT de AlphaFold/Boltz, la misma que usa el render
PAL_PLDDT <- c("#FF7D45", "#FFDB13", "#65CBF3", "#0053D6")

t17 <- suppressMessages(read_csv(file.path(TABLAS, "T17_pLDDT_EglNAC66.csv"), show_col_types = FALSE))
t16 <- suppressMessages(read_csv(file.path(TABLAS, "T16_arquitectura_dominios.csv"), show_col_types = FALSE))
t4  <- suppressMessages(read_csv(file.path(TABLAS, "T04_motivos_MEME.csv"), show_col_types = FALSE))
t5  <- suppressMessages(read_csv(file.path(TABLAS, "T05_motivos_sitios.csv"), show_col_types = FALSE))

foc <- t16 |> filter(simbolo == GEN_FOCAL)
stopifnot(nrow(foc) == 1, nrow(t17) == foc$largo_aa)
DOM_INI <- foc$ali_ini; DOM_FIN <- foc$ali_fin

subE <- t5 |> dplyr::rename(simbolo = simbolo_EglNAC) |>
  left_join(t4 |> select(motivo, subdominio), by = "motivo") |>
  filter(simbolo == GEN_FOCAL, subdominio == "E") |>
  summarise(ini = min(inicio), fin = max(inicio + ancho - 1))

# --- residuos que leen bases, transferidos del co-cristal ANAC019-ADN -------
# Esto es lo que convierte el panel de "aqui esta el modelo" en un resultado: el
# perfil solo dice que el dominio esta bien resuelto; T18 dice ademas DONDE lee
# el ADN. Las posiciones son del co-cristal 3SWP y su equivalente en el gen
# focal quedo verificado contra la secuencia, aminoacido a aminoacido, en
# 09_verificacion_equivalencias.py. Se comprueba aqui otra vez porque el script
# tiene que fallar solo si alguien cambia T18 o la proteina de referencia.
t18 <- suppressMessages(read_csv(file.path(TABLAS, "T18_lectura_bases_3SWP.csv"), show_col_types = FALSE))
qr  <- suppressMessages(read_csv(file.path(RAIZ, "qc_rmsd_modelo_vs_plantillas.csv"), show_col_types = FALSE))

lect <- t18 |>
  mutate(pos = as.integer(gsub("\\D", "", residuo_EglNAC)),
         aa  = substr(residuo_EglNAC, 1, 1)) |>
  left_join(t17, by = c("pos" = "res")) |>
  arrange(pos)
stopifnot(nrow(lect) == 8, all(lect$conservado),
          all(lect$pos >= DOM_INI), all(lect$pos <= DOM_FIN))

# Los ocho caen en dos grupos separados; se etiquetan como rangos y no uno a uno
# porque a esta escala ocho rotulos se pisan entre si. Los rangos se calculan,
# no se teclean, para que sigan al dato si T18 cambia.
rangos <- lect |> mutate(g = cumsum(c(1, as.integer(diff(pos) > 5)))) |>
  group_by(g) |>
  summarise(t = if (min(pos) == max(pos)) as.character(min(pos))
                else sprintf("%d\u2013%d", min(pos), max(pos)), .groups = "drop") |>
  # separador corto: con " and " la segunda linea se salia del panel por la derecha
  pull(t) |> paste(collapse = ", ")
etq_lect <- sprintf("%d base-reading residues (3SWP)\npositions %s\nconserved %d/%d",
                    nrow(lect), rangos, sum(lect$conservado), nrow(lect))
q3 <- qr |> filter(grepl("^3SWP", referencia))

p_dom <- t17 |> filter(res >= DOM_INI, res <= DOM_FIN) |> pull(plddt)
p_ct  <- t17 |> filter(res > DOM_FIN) |> pull(plddt)

registrar_cifra("f08_largo", nrow(t17), "T17_pLDDT_EglNAC66.csv")
registrar_cifra("f08_dominio", sprintf("%d-%d", DOM_INI, DOM_FIN),
                "T16: alineamiento del perfil PF02365")
registrar_cifra("f08_plddt_dominio", sprintf("%.0f", mean(p_dom)), "T17")
registrar_cifra("f08_plddt_cterm", sprintf("%.0f", mean(p_ct)), "T17")
registrar_cifra("f08_largo_cterm", length(p_ct), "T17")
registrar_cifra("f08_bajo_umbral", sum(p_ct < UMBRAL_DESORDEN),
                sprintf("T17: residuos C-terminales con pLDDT < %d", UMBRAL_DESORDEN))
registrar_cifra("f08_umbral_desorden", UMBRAL_DESORDEN, "09_estructura.R")
registrar_cifra("f08_escala_compartida",
                "la barra de pLDDT se aplica a los dos paneles: render y perfil usan la misma escala",
                "09_estructura.R")
registrar_cifra("f08_plddt_min_max",
                sprintf("%.0f a %.0f", min(t17$plddt), max(t17$plddt)), "T17")
registrar_cifra("f08_subdominio_E", sprintf("%d-%d", subE$ini, subE$fin),
                "T05 x T04: sitios del motivo asignado al subdominio E")
registrar_cifra("f08_lectores_n", nrow(lect), "T18_lectura_bases_3SWP.csv")
registrar_cifra("f08_lectores_lista", paste(lect$residuo_EglNAC, collapse = ", "), "T18")
registrar_cifra("f08_lectores_conservados", sprintf("%d de %d", sum(lect$conservado), nrow(lect)),
                "T18: identidad con ANAC019 en las posiciones que contactan bases")
registrar_cifra("f08_lectores_contactos", sum(lect$n_contactos_base),
                "T18: contactos con bases nitrogenadas en el co-cristal 3SWP")
registrar_cifra("f08_lectores_plddt_min", sprintf("%.0f", min(lect$plddt)),
                "T17 evaluada en las posiciones de T18: el peor de los ocho")
registrar_cifra("f08_rmsd_3swp",
                sprintf("%.2f-%.2f \u00c5 sobre %d-%d C\u03b1", min(q3$rmsd_A), max(q3$rmsd_A),
                        min(q3$n_pares), max(q3$n_pares)),
                "qc_rmsd_modelo_vs_plantillas.csv: cadenas A y D de 3SWP")
registrar_cifra("f08_identidad_ANAC019", sprintf("%.1f %%", max(q3$identidad_pct)),
                "qc_rmsd_modelo_vs_plantillas.csv")
# Metodos contrasta la superposicion contra la forma libre y la unida a ADN; la
# libre (1UT7-B) no estaba registrada y el documento no podia citarla.
q1 <- qr |> filter(grepl("^1UT7", referencia))
registrar_cifra("f10_rmsd_1ut7",
                sprintf("%.2f \u00c5 sobre %d C\u03b1", q1$rmsd_A[1], q1$n_pares[1]),
                "qc_rmsd_modelo_vs_plantillas.csv: 1UT7-B, ANAC019 libre")

# El render sale del programa de 3D con margenes blancos anchos: el contenido
# ocupaba el 4 % de los pixeles y el monomero quedaba diminuto dentro del panel.
# Se recorta a la caja del dibujo antes de rasterizar.
recortar_blanco <- function(im, margen = 0.01) {
  lum <- if (length(dim(im)) == 3) apply(im[, , 1:3, drop = FALSE], c(1, 2), min) else im
  fil <- which(apply(lum < 0.96, 1, any)); col <- which(apply(lum < 0.96, 2, any))
  if (!length(fil) || !length(col)) return(im)
  mf <- max(1, round(margen * length(fil))); mc <- max(1, round(margen * length(col)))
  f1 <- max(1, min(fil) - mf); f2 <- min(nrow(lum), max(fil) + mf)
  c1 <- max(1, min(col) - mc); c2 <- min(ncol(lum), max(col) + mc)
  if (length(dim(im)) == 3) im[f1:f2, c1:c2, , drop = FALSE] else im[f1:f2, c1:c2]
}

# --- panel a: render 3D ------------------------------------------------------
ruta_render <- file.path(FIGURAS, "_render_EglNAC66.png")
if (!file.exists(ruta_render)) stop("falta ", basename(ruta_render))
im <- recortar_blanco(readPNG(ruta_render))
pa <- ggplot() +
  annotation_custom(rasterGrob(im, interpolate = TRUE)) +
  coord_fixed(ratio = dim(im)[1] / dim(im)[2], xlim = c(0, 1), ylim = c(0, 1),
              expand = FALSE) +
  theme_void()

# --- panel b: perfil de pLDDT por residuo, coloreado con la misma escala ----
seg <- t17 |> mutate(res2 = dplyr::lead(res), p2 = dplyr::lead(plddt)) |>
  filter(!is.na(res2)) |> mutate(pm = (plddt + p2) / 2)

pb <- ggplot() +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0, ymax = 50, fill = "#FDECEA") +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 50, ymax = 70, fill = "#FEF7DF") +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 70, ymax = 90, fill = "#E6F4FB") +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 90, ymax = 100, fill = "#DDE7FA") +
  annotate("rect", xmin = subE$ini, xmax = subE$fin, ymin = -Inf, ymax = Inf,
           fill = "grey45", alpha = 0.22) +
  geom_hline(yintercept = UMBRAL_DESORDEN, linetype = "dotted",
             colour = "grey35", linewidth = 0.3) +
  geom_segment(data = seg, aes(x = res, xend = res2, y = plddt, yend = p2, colour = pm),
               linewidth = 0.45, lineend = "round") +
  annotate("point", x = lect$pos, y = 100.5, shape = 25, size = 1.15,
           fill = COL_DESTACADO, colour = COL_DESTACADO, stroke = 0.15) +
  annotate("segment", x = nrow(t17) * 0.505, xend = max(lect$pos) + 2,
           y = 86, yend = 98, linewidth = 0.25, colour = COL_DESTACADO) +
  annotate("text", x = nrow(t17) * 0.51, y = 84, hjust = 0, vjust = 1, size = 2.1,
           colour = COL_DESTACADO, label = etq_lect) +
  annotate("segment", x = DOM_INI, xend = DOM_FIN, y = 104, yend = 104,
           linewidth = 0.9, colour = "grey20") +
  annotate("text", x = (DOM_INI + DOM_FIN) / 2, y = 108, size = 2.2, colour = "grey20",
           label = sprintf("NAM core (PF02365), %d-%d", DOM_INI, DOM_FIN)) +
  annotate("text", x = subE$fin + 6, y = 96, hjust = 0, size = 2.1, colour = "grey30",
           label = "subdomain E\n(outside PF02365)") +
  annotate("text", x = nrow(t17) * 0.98, y = 8, hjust = 1, size = 2.1, colour = "grey40",
           label = sprintf("disorder threshold (pLDDT %d)", UMBRAL_DESORDEN)) +
  scale_colour_stepsn(colours = PAL_PLDDT, breaks = c(50, 70, 90),
                      limits = c(0, 100), name = "pLDDT",
                      guide = guide_coloursteps(barheight = unit(2.2, "mm"),
                                                barwidth = unit(26, "mm"),
                                                direction = "horizontal",
                                                title.position = "left")) +
  scale_x_continuous(breaks = c(1, 50, 100, 150, 200, 250, 300, nrow(t17)),
                     expand = expansion(mult = c(0.01, 0.02))) +
  scale_y_continuous(limits = c(0, 112), breaks = seq(0, 100, 20), expand = c(0, 0)) +
  labs(x = "residue", y = "pLDDT") +
  theme_nac(8) +
  theme(legend.position = "bottom", legend.box.margin = margin(-4, 0, 0, 0),
        panel.grid = element_blank())

p <- (pa | pb) + plot_layout(widths = c(1, 1.45)) +
  plot_annotation(tag_levels = "a")

# --- panel c: dimero con ADN, solo si el co-plegamiento ya termino -----------
pdb_dim <- file.path(FIGURAS, "_render_dimero_ADN.png")
if (file.exists(pdb_dim)) {
  imd <- recortar_blanco(readPNG(pdb_dim))
  pc <- ggplot() + annotation_custom(rasterGrob(imd, interpolate = TRUE)) +
    coord_fixed(ratio = dim(imd)[1] / dim(imd)[2], xlim = c(0, 1), ylim = c(0, 1),
                expand = FALSE) + theme_void()
  p <- (pa | pb) / pc + plot_layout(heights = c(1, 0.95)) +
    plot_annotation(tag_levels = "a")
  registrar_cifra("f08_dimero", "incluido", basename(pdb_dim))
  guardar_objeto(p, "P_estructura")
  guardar_figura(p, "F08_monomero_EglNAC66", alto = 5.40)
} else {
  registrar_cifra("f08_dimero", "pendiente",
                  "el co-plegamiento con Boltz-2 no ha producido estructura todavia")
  message("  nota: sin _render_dimero_ADN.png, F08 sale con dos paneles")
  guardar_objeto(p, "P_estructura")
  guardar_figura(p, "F08_monomero_EglNAC66", alto = 2.70)
}

volcar_cifras("f08")
