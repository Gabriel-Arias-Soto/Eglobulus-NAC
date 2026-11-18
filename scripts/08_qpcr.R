# 08_qpcr.R -> panel de validacion qRT-PCR frente a RNA-seq (ggpubr)
# Lee la T10 que reconstruye 18_correlacion_qpcr.R: qPCR desde los RQ de T08 y
# RNA-seq desde el DE contra E. globulus, unidos por locus y no por nombre de
# simbolo. La version anterior de la tabla enfrentaba el qPCR con cifras de
# otros loci y con un 0,000 que marcaba "no evaluado", asi que las cifras de
# esta figura NO son comparables con las de la version sometida; todas se
# recalculan aqui y se registran para la leyenda.
# Las barras de error solo pueden ir en el eje x: la salida cuasi-verosimil de
# edgeR no da error estandar del log2FC.
if (!exists("ANCHO_COLUMNA")) source("scripts/00_estilo.R")
suppressPackageStartupMessages({library(ggpubr); library(ggrepel)})

t10 <- suppressMessages(read_csv(file.path(TABLAS, "T10_correlacion_qpcr_rnaseq.csv"),
                                 show_col_types = FALSE))
censo <- leer_censo()

d <- t10 |>
  dplyr::rename(simbolo = simbolo_unico_EglNAC) |>
  filter(!is.na(qpcr_log2FC)) |>
  left_join(censo |> select(simbolo, subfamilia_final), by = "simbolo") |>
  mutate(subfamilia_final = factor(subfamilia_final, levels = c(SUBFAM_NIVELES, SUBFAM_NO_EVAL)),
         discordante = sign(qpcr_log2FC) * sign(rnaseq_log2FC) < 0,
         etiqueta = sprintf("%s | %s", simbolo, contraste))

# Un par es evaluable si el locus tiene estimacion de RNA-seq. Con la tabla
# anterior "no evaluado" se codificaba con un 0,000 exacto, indistinguible de
# una estimacion nula; ahora es un ausente explicito, que solo aparece si el
# locus no pasa el filtro de expresion de edgeR.
# Ejes simetricos y comunes a los tres paneles. La escala tiene que ser la
# misma en x y en y para que la diagonal 1:1 salga a 45 grados y la
# concordancia se lea como distancia a esa recta; por eso mas abajo coord_fixed.
# Lo que coord_fixed necesita para que la diagonal salga a 45 grados es que las
# dos escalas midan lo mismo, NO que esten centradas en cero. La version
# anterior imponia ademas simetria respecto a cero y redondeo a la unidad: con
# un qRT-PCR que baja hasta -6,9 y un RNA-seq que no pasa de +5,6, eso obligaba
# a dibujar un marco de 14 unidades de lado para un dato que ocupa 12,5, y el
# exceso caia entero en esquinas donde no hay ningun gen.
VALS <- c(d$qpcr_log2FC + d$ee, d$qpcr_log2FC - d$ee, d$rnaseq_log2FC)
RG   <- range(VALS, na.rm = TRUE)
LIMS <- RG + c(-1, 1) * 0.04 * diff(RG)
registrar_cifra("f03_limite_eje", sprintf("%.2f a %.2f", LIMS[1], LIMS[2]),
                "T10: rango comun a los dos ejes, no simetrico")

d$evaluable <- !is.na(d$rnaseq_log2FC)
ev <- filter(d, evaluable)
# n = 11 sin empates: cor.test usa por defecto la p exacta de permutacion, que
# es la que corresponde a este tamano; el memo cito la aproximacion asintotica.
sp_ev    <- suppressWarnings(cor.test(ev$qpcr_log2FC, ev$rnaseq_log2FC, method = "spearman"))
sp_todos <- suppressWarnings(cor.test(d$qpcr_log2FC,  d$rnaseq_log2FC,  method = "spearman"))
sesgo    <- median(ev$rnaseq_log2FC - ev$qpcr_log2FC)
concord  <- sum(sign(ev$qpcr_log2FC) == sign(ev$rnaseq_log2FC))
ev$discordante <- sign(ev$qpcr_log2FC) != sign(ev$rnaseq_log2FC)
d$clase <- ifelse(!d$evaluable, "no evaluado", ifelse(
  sign(d$qpcr_log2FC) != sign(d$rnaseq_log2FC), "signo opuesto", "mismo signo"))
d$clase <- factor(d$clase, levels = c("mismo signo", "signo opuesto", "no evaluado"))
# los paneles siguen la secuencia del ensayo (sequia, frio antes y tras la
# helada) y no el orden alfabetico del codigo
ORDEN_CONTR <- c("C_DS", "NA_CABF", "NA_CAAF")
d$panel <- factor(etq_contraste(d$contraste),
                  levels = unname(etq_contraste(intersect(ORDEN_CONTR, unique(d$contraste)))))
disc <- filter(ev, discordante)
misma_dir <- nrow(disc) > 0 && all(disc$qpcr_log2FC < 0 & disc$rnaseq_log2FC > 0)

# --- rho por contraste, dentro de cada viñeta --------------------------------
# El guion ya calculaba una rho por contraste y solo la volcaba al registro de
# cifras; en la figura no aparecia. El cuadrante inferior derecho esta vacio en
# los tres paneles (ningun gen sube por qRT-PCR y baja por RNA-seq), asi que
# ese hueco es donde cabe sin tapar nada. Las coordenadas viajan como columnas
# del propio data: aes() no puede mirar a `lim`, porque el objeto se guarda en
# .rds y se recompone en 12_paneles.R, donde esa variable ya no existe.
rho_panel <- d |>
  filter(evaluable) |>
  group_by(panel) |>
  summarise(rho = suppressWarnings(cor(qpcr_log2FC, rnaseq_log2FC, method = "spearman")),
            k = dplyr::n(), .groups = "drop") |>
  mutate(x = LIMS[2] - 0.02 * diff(LIMS), y = LIMS[1] + 0.02 * diff(LIMS),
         txt = ifelse(k >= 3, sprintf("\u03c1 = %.2f\nn = %d", rho, k),
                      sprintf("n = %d", k)))

# tres lineas cortas: a 6,6 pt sobre 6,3 in no caben mas de ~120 caracteres
# El respaldo del eje de RNA-seq se lee de T33, no se afirma de memoria.
dis33 <- suppressMessages(read_csv(file.path(TABLAS, "T33_diseno_rnaseq_hoja.csv"),
                                   na = character(), show_col_types = FALSE))
nbio <- tapply(dis33$biomuestra, dis33$condicion, function(z) length(unique(z)))
ETQ_REPLICAS <- if (all(nbio == 1)) {
  sprintf(paste("qRT-PCR measures biological replicates; the public RNA-seq contributes %d",
                "libraries from a single biosample per condition, so its log\u2082FC is descriptive."),
          round(mean(table(dis33$condicion))))
} else {
  sprintf("RNA-seq: %d libraries over %d biosamples.",
          nrow(dis33), length(unique(dis33$biomuestra)))
}
registrar_cifra("f03_respaldo_rnaseq", ETQ_REPLICAS, "T33.biomuestra")

ETQ_STAT <- paste(
  sprintf("Spearman \u03c1 = %.2f (p = %s exact, n = %d pairs%s).",
          sp_ev$estimate, format.pval(sp_ev$p.value, digits = 2), nrow(ev),
          if (nrow(ev) < nrow(d))
            sprintf(" evaluable; including the %d without an estimate \u03c1 = %.2f, n = %d",
                    nrow(d) - nrow(ev), sp_todos$estimate, nrow(d)) else ""),
  sprintf("Sign agrees in %d of %d%s.", concord, nrow(ev),
          if (nrow(disc) == 0) "" else
            sprintf("; the %d discordant %s", nrow(disc),
                    if (misma_dir) "are all qRT-PCR negative against RNA-seq positive"
                    else "do not share direction")),
  sprintf("Median bias %+.2f log\u2082 units (RNA-seq \u2212 qRT-PCR).", sesgo),
  # Las dos plataformas no tienen el mismo respaldo: el qRT-PCR mide replicas
  # biologicas y el RNA-seq publico cuelga las tres librerias de cada condicion
  # de una sola biomuestra (T33.biomuestra). Decirlo aqui evita que el lector
  # lea la correlacion como un contraste entre dos estimaciones equivalentes.
  ETQ_REPLICAS,
  sep = "\n")

# ggplot no ajusta el subtitulo al ancho del lienzo: la linea que sobra se
# recorta por la derecha sin aviso, y eso es justo lo que le pasaba a la
# advertencia sobre las replicas. Se parte a mano, respetando los saltos que ya
# separan cada afirmacion. 104 caracteres es lo que entra a 6,5 pt en
# ANCHO_COLUMNA; medido, no supuesto (ver comprobacion al final del guion).
ANCHO_SUBTITULO <- 104
ETQ_STAT <- paste(vapply(strsplit(ETQ_STAT, "\n", fixed = TRUE)[[1]],
                         function(l) paste(strwrap(l, width = ANCHO_SUBTITULO),
                                           collapse = "\n"),
                         character(1), USE.NAMES = FALSE),
                  collapse = "\n")

registrar_cifra("f03_rho_evaluables", round(sp_ev$estimate, 2), "T10")
registrar_cifra("f03_p_evaluables", signif(sp_ev$p.value, 3), "T10")
registrar_cifra("f03_n_evaluables", nrow(ev), "T10")
registrar_cifra("f03_concordancia_signo", sprintf("%d de %d", concord, nrow(ev)), "T10")
registrar_cifra("f03_sesgo_mediano", round(sesgo, 2), "T10")
registrar_cifra("f03_discordantes_misma_direccion", misma_dir, "T10")
for (cc in unique(d$contraste)) {
  x <- filter(ev, contraste == cc)
  if (nrow(x) >= 3) registrar_cifra(paste0("f03_rho_", cc),
      round(suppressWarnings(cor(x$qpcr_log2FC, x$rnaseq_log2FC, method = "spearman")), 2), "T10")
}

p <- ggplot(d, aes(x = qpcr_log2FC, y = rnaseq_log2FC)) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = 0, ymax = Inf, fill = "grey96") +
  annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = 0, fill = "grey96") +
  geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.25) +
  geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.25) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              colour = "grey45", linewidth = 0.3) +
  geom_errorbar(aes(xmin = qpcr_log2FC - ee, xmax = qpcr_log2FC + ee),
                orientation = "y", width = 0, colour = "grey50", linewidth = 0.35) +
  geom_text(data = rho_panel, aes(x = x, y = y, label = txt), inherit.aes = FALSE,
            hjust = 1, vjust = 0, size = 2.3, colour = "grey35", lineheight = 1.05) +
  geom_point(aes(fill = subfamilia_final, shape = clase),
             size = 2.6, colour = "grey20", stroke = 0.35) +
  geom_text_repel(aes(label = simbolo), size = 2.1, colour = "grey25",
                  min.segment.length = 0, segment.size = 0.2,
                  segment.colour = "grey65", max.overlaps = Inf,
                  box.padding = 0.35, force = 3, seed = 7) +
  facet_wrap(~panel, nrow = 1) +
  scale_fill_manual(values = PAL_SUBFAM, name = "subfamily", drop = TRUE) +
  scale_shape_manual(values = c("mismo signo" = 21, "signo opuesto" = 24,
                                "no evaluado" = 22), name = NULL, drop = TRUE,
                     labels = c("mismo signo" = "same sign on both platforms",
                                "signo opuesto" = "opposite sign",
                                "no evaluado" = "no RNA-seq estimate")) +
  coord_fixed(xlim = LIMS, ylim = LIMS) +
  labs(x = expression(log[2]~"FC by qRT-PCR"),
       y = expression(atop(log[2]~"FC by RNA-seq",
                           scriptstyle("(3 libraries per condition, ref. "*italic("E. globulus")*")"))),
       subtitle = ETQ_STAT) +
  guides(fill  = guide_legend(override.aes = list(shape = 21, size = 2), order = 1, nrow = 1),
         shape = guide_legend(override.aes = list(fill = "grey75", size = 2), order = 2, nrow = 1)) +
  theme_nac(8) +
  theme(legend.position = "bottom", legend.box = "vertical",
        legend.box.spacing = unit(2, "pt"), legend.spacing.y = unit(0, "pt"),
        legend.margin = margin(0, 0, 0, 0),
        strip.text = element_text(size = 6.8, lineheight = 1.05),
        plot.subtitle = element_text(size = 6.5, colour = "grey25", lineheight = 1.15),
        panel.grid.major.x = element_line(colour = "grey94", linewidth = 0.2),
        panel.spacing = unit(4, "pt"))

guardar_objeto(p, "P_qpcr")
guardar_figura(p, "F03_correlacion_qpcr_rnaseq", alto = 3.70)

# Comprobacion del ancho del subtitulo: se mide la linea mas larga con la
# tipografia y el cuerpo reales y se compara con el lienzo. Si alguna vez las
# cifras crecen y una linea vuelve a desbordar, el aviso salta aqui y no en la
# revision de la revista.
local({
  tf <- tempfile(fileext = ".png")
  grDevices::png(tf, width = ANCHO_COLUMNA, height = 1, units = "in", res = 300)
  ancho_in <- max(vapply(strsplit(ETQ_STAT, "\n", fixed = TRUE)[[1]],
                         function(l) as.numeric(grid::convertWidth(
                           grid::stringWidth(l), "in")) * 6.5 / 12,
                         numeric(1), USE.NAMES = FALSE))
  grDevices::dev.off(); unlink(tf)
  if (ancho_in > ANCHO_COLUMNA)
    warning(sprintf("subtitulo de F03: %.2f in sobre un lienzo de %.2f in; baja ANCHO_SUBTITULO",
                    ancho_in, ANCHO_COLUMNA), call. = FALSE)
})

volcar_cifras("f03")
