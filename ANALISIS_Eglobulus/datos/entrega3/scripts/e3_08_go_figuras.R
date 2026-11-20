# e3_08_go_figuras.R -> figuras/F13_GO_dotplot, F13S_GO_arbol, F14_GO_subfamilias
#                       tablas/T42_anotacion_GO_NAC189.csv, T43_enriquecimiento_GO.csv,
#                       tablas/T44_compareCluster_GO.csv, tablas/_cifras/e3_go.csv
#
# Entrega 3, paso 3.4 (figuras y tablas de GO). No recalcula nada: lee lo que dejan
# e3_05_go_enriquecimiento.R y e3_07_go_comparecluster.R en datos/entrega3/go y lo dibuja con
# el estilo del juego (scripts/00_estilo.R). Corre en r-nac, que tiene ggtree y patchwork
# pero no clusterProfiler; por eso las figuras se hacen con ggplot sobre los TSV en vez de
# enrichplot::dotplot/treeplot, cuyo aspecto reproducen.
#
# Antes de dibujar se comprueba (si algo no cuadra, el script para):
#   - que la regla de significacion (p y p.adjust <= 0,05; q <= 0,2) da los recuentos de
#     resumen_go.json y resumen_go_figuras.json, en las tres ontologias y en los contrastes;
#   - que los 17 representativos de simplify estan entre los 35 significativos y que las
#     hojas del arbol son esos 35;
#   - que las listas de compareCluster (cc_listas.tsv) salen de T41 + subfamilia;
#   - que, para cada termino, los EglNAC que se deducen de los AGI (unidad del analisis
#     principal) son los mismos que da el contraste con genes de E. globulus como unidades.
#
# Decisiones de dibujo
#   F13  dotplot de los 17 terminos BP representativos, como en el articulo 2: x = GeneRatio,
#        tamano = Count, color = -log10 p.adjust. Se ordenan por los 5 grupos del arbol (F13S)
#        y cada grupo lleva el mismo nombre en las dos figuras. Un dagger marca los terminos
#        que dejan de ser significativos con los genes de E. globulus como unidades.
#   F13S arbol de los 35 significativos (similitud de Wang, ward.D, 5 grupos), como
#        enrichplot::treeplot. Negrita = representativo (el que sale en F13). Misma escala de
#        color y de tamano que F13, para que un punto se lea igual en las dos.
#   F14  compareCluster por subfamilia. Filas = union de los representativos de cada
#        subfamilia. Hay punto en toda subfamilia donde el termino es significativo con >= 2
#        genes, sea o no su representativo alli: si solo se pintara el representativo, el
#        hueco se leeria como "no enriquecido" cuando el termino si lo esta y simplify lo quito
#        por redundante con otro. Tamano = GeneRatio, porque las listas tienen de 3 a 6 AGI y
#        Count no es comparable entre columnas. Las subfamilias probadas sin ningun punto
#        llevan rotulo en su columna (un hueco sin marca se leeria como cero sin explicar por
#        que).
#
# Se ejecuta desde la RAIZ de "Revision NAC":
#   Rscript ANALISIS_Eglobulus/datos/entrega3/scripts/e3_08_go_figuras.R \
#     > ANALISIS_Eglobulus/logs/entrega3_09_go_figuras.log 2>&1
if (!exists("ANCHO_COLUMNA")) source("scripts/00_estilo.R")
suppressPackageStartupMessages({
  library(ggtree); library(ape); library(patchwork); library(jsonlite); library(ggforce)
})

A   <- dirname(FIGURAS)                         # ANALISIS_Eglobulus
GOD <- file.path(A, "datos/entrega3/go")
P_CORTE <- 0.05; Q_CORTE <- 0.2; MIN_COUNT <- 2
PARCIAL <- "no evaluada (dominio parcial)"     # rotulo de T01/T41 (e3_07 lo excluye igual)
TAM_TXT <- 8 / .pt                              # 8 pt, el minimo de la norma de figuras de Springer, para geom_text
DAGA    <- "\u2020"

rd      <- function(f) read.delim(file.path(GOD, f), stringsAsFactors = FALSE, check.names = FALSE)
signif_ <- function(d) d$pvalue <= P_CORTE & d$p.adjust <= P_CORTE & (is.na(d$qvalue) | d$qvalue <= Q_CORTE)
frac    <- function(s) vapply(strsplit(s, "/", fixed = TRUE),
                              function(x) as.numeric(x[1]) / as.numeric(x[2]), numeric(1))
mlp     <- function(p) -log10(p)
num_sim <- function(s) as.integer(sub("^EglNAC", "", s))
mm2in   <- function(x) x / 25.4

ONT <- c("BP", "MF", "CC")
rg  <- fromJSON(file.path(GOD, "resumen_go.json"))
rf  <- fromJSON(file.path(GOD, "resumen_go_figuras.json"))
eg  <- setNames(lapply(ONT, function(o) rd(sprintf("enrichGO_%s_todos.tsv", o))), ONT)
ee  <- setNames(lapply(ONT, function(o) rd(sprintf("enricherEgl_%s_todos.tsv", o))), ONT)
se  <- setNames(lapply(ONT, function(o) rd(sprintf("sens_maxInf_%s_todos.tsv", o))), ONT)
bs  <- rd("enrichGO_BP_simplify.tsv")
arb <- rd("arbol_BP_grupos.tsv")
cst <- rd("cc_subfamilia_BP_todos.tsv"); css <- rd("cc_subfamilia_BP_simplify.tsv")
cgt <- rd("cc_grupo_BP_todos.tsv");      cgs <- rd("cc_grupo_BP_simplify.tsv")
cl  <- rd("cc_listas.tsv")
ad  <- rd("anotacion_directa_lista.tsv")
tr  <- read.tree(file.path(GOD, "arbol_BP_signif.nwk"))
t41 <- read.csv(file.path(TABLAS, "T41_homologos_Ath_NAC189.csv"), stringsAsFactors = FALSE,
                check.names = FALSE)

# --- 1. comprobaciones contra los resumenes ----------------------------------------------
for (o in ONT) {
  r <- rg[[o]]; s <- rf[[paste0("sensibilidad_", o)]]
  ok_o <- c(probados = nrow(eg[[o]]) == r$terminos_probados,
            signif   = sum(signif_(eg[[o]])) == r$significativos,
            egl_prob = nrow(ee[[o]]) == r$probados_contraste_Egl,
            egl      = sum(signif_(ee[[o]])) == r$significativos_contraste_Egl,
            compart  = length(intersect(eg[[o]]$ID[signif_(eg[[o]])],
                                        ee[[o]]$ID[signif_(ee[[o]])])) == r$compartidos,
            sens     = nrow(se[[o]]) == s$probados && sum(signif_(se[[o]])) == s$significativos,
            sens_col = identical(signif_(se[[o]]), as.logical(se[[o]]$significativo)))
  if (!all(ok_o)) stop(o, ": no cuadra con los resumenes -> ", paste(names(ok_o)[!ok_o], collapse = ", "))
}
sig35 <- eg$BP[signif_(eg$BP), ]
tr$tip.label <- sub("^GO-", "GO:", tr$tip.label)   # write.tree cambia ':' por '-' en las hojas
stopifnot(nrow(bs) == rg$BP$tras_simplify, all(bs$ID %in% sig35$ID),
          setequal(tr$tip.label, sig35$ID), setequal(arb$ID, sig35$ID),
          all(as.integer(table(arb$grupo_arbol)[names(rf$arbol$grupos)]) == unlist(rf$arbol$grupos)),
          MIN_COUNT == rf$parametros$MIN_COUNT_figura)
sig_egl <- setNames(lapply(ONT, function(o) ee[[o]]$ID[signif_(ee[[o]])]), ONT)

# genes de la lista: los 141 EglNAC con homologo que pasa los filtros de T41 (46 AGI)
ok41 <- t41$pasa_filtros %in% c(TRUE, "True")
z <- t41[ok41, c("simbolo", "locus", "AGI", "subfamilia")]
z$grupo <- ifelse(z$subfamilia == PARCIAL, PARCIAL, sub("^([IVX]+).*$", "\\1", z$subfamilia))
stopifnot(all(nzchar(z$AGI)), !any(nzchar(t41$AGI[!ok41]) & !is.na(t41$AGI[!ok41])),
          nrow(z) == rg$parametros$lista_Egl, length(unique(z$AGI)) == rg$parametros$lista_AGI)
for (i in seq_len(nrow(cl))) {
  col <- if (cl$nivel[i] == "grupo") "grupo" else "subfamilia"
  if (!setequal(strsplit(cl$AGI[i], ";")[[1]], unique(z$AGI[z[[col]] == cl$grupo[i]])))
    stop("cc_listas no sale de T41 para ", cl$nivel[i], " ", cl$grupo[i])
}
sim_de <- function(ids, clave = "AGI", pool = z) {
  s <- unique(pool$simbolo[pool[[clave]] %in% ids])
  paste(s[order(num_sim(s))], collapse = ";")
}
partir <- function(x) strsplit(x, "/", fixed = TRUE)

# mismos EglNAC por las dos vias en todo termino probado en ambas
dif_via <- 0L
for (o in ONT) {
  m <- merge(eg[[o]][, c("ID", "geneID")], ee[[o]][, c("ID", "geneID")], by = "ID", suffixes = c("_agi", "_egl"))
  dif_via <- dif_via + sum(vapply(partir(m$geneID_agi), sim_de, "") !=
                           vapply(partir(m$geneID_egl), sim_de, "", clave = "locus"))
}
if (dif_via) stop(dif_via, " terminos con EglNAC distintos segun la unidad (AGI o gen de E. globulus)")

# --- 2. nombres de los grupos del arbol ----------------------------------------------------
# Cada nombre resume los terminos de su grupo y va anclado a uno de ellos: si el arbol cambia
# de numeracion el nombre sigue a su termino, y si dos anclas caen en el mismo grupo, para.
ANCLAS <- c("positive regulation of biosynthetic process"        = "Positive regulation of biosynthesis",
            "positive regulation of DNA-templated transcription" = "Regulation of transcription",
            "plant-type secondary cell wall biogenesis"          = "Cell wall biogenesis",
            "xylem development"                                  = "Vascular tissue development",
            "plant organ senescence"                             = "Organ development and senescence")
g_ancla <- arb$grupo_arbol[match(names(ANCLAS), arb$Description)]
stopifnot(!anyNA(g_ancla), !anyDuplicated(g_ancla), setequal(g_ancla, arb$grupo_arbol))
arb$nombre_grupo <- unname(ANCLAS)[match(arb$grupo_arbol, g_ancla)]

pt <- ggtree(tr, ladderize = FALSE, linewidth = 0.3, colour = "grey35")
yt <- pt$data[pt$data$isTip, c("label", "y")]
arb$y <- yt$y[match(arb$ID, yt$label)]
# cada grupo ocupa hojas seguidas (cutree sobre el mismo hclust); si no, los corchetes mentirian
stopifnot(!anyNA(arb$y), all(tapply(arb$y, arb$grupo_arbol, function(y) diff(range(y)) + 1 == length(y))))
ORDEN_GRUPOS <- names(sort(tapply(arb$y, arb$nombre_grupo, mean), decreasing = TRUE))  # de arriba abajo

# --- 3. escalas compartidas F13 / F13S ------------------------------------------------------
LIM_MLP   <- range(mlp(sig35$p.adjust))
LIM_COUNT <- c(0, max(sig35$Count))
# negrita, como legend.title. El signo menos va como caracter de texto (U+2212) y no como
# operador de plotmath: cairo_pdf() dibuja los operadores con la fuente de simbolos, que
# fontconfig resuelve a DejaVu Sans, y el PDF acababa con una segunda fuente para un solo glifo.
TIT_P     <- expression(bold("\u2212log"[10]*"(adj. "*bolditalic(P)*")"))
esc_fill <- function(lim, horizontal = FALSE)
  scale_fill_viridis_c(option = "mako", direction = -1, begin = 0.08, end = 0.9, limits = lim,
                       breaks = scales::breaks_pretty(4), name = TIT_P,
                       guide = guide_colourbar(order = 1, theme = if (horizontal)
                         theme(legend.key.width = unit(20, "mm"), legend.key.height = unit(3, "mm")) else
                         theme(legend.key.width = unit(3, "mm"), legend.key.height = unit(16, "mm"))))
esc_count <- scale_size_area(max_size = 4, limits = LIM_COUNT, breaks = c(2, 5, 10, 15),
                             name = expression(bold("Genes ("*bolditalic(k)*")")),
                             guide = guide_legend(order = 2))

# --- 4. F13: dotplot de los representativos ----------------------------------------------
d13 <- bs
d13$gr <- frac(d13$GeneRatio); d13$mlp <- mlp(d13$p.adjust)
d13$nombre_grupo <- arb$nombre_grupo[match(d13$ID, arb$ID)]
d13$replica  <- d13$ID %in% sig_egl$BP
d13$etiqueta <- ifelse(d13$replica, d13$Description, paste(d13$Description, DAGA))
d13 <- d13[order(match(d13$nombre_grupo, ORDEN_GRUPOS), -d13$gr, d13$p.adjust), ]   # orden de lectura
d13$etiqueta     <- factor(d13$etiqueta, levels = rev(unique(d13$etiqueta)))   # ggplot pone el 1.o abajo
d13$nombre_grupo <- factor(d13$nombre_grupo, levels = intersect(ORDEN_GRUPOS, d13$nombre_grupo))
n_lista <- rg$BP$n_lista_anotada
stopifnot(all(as.integer(sub("/.*", "", d13$GeneRatio)) == d13$Count),
          all(as.integer(sub(".*/", "", d13$GeneRatio)) == n_lista))

p13 <- ggplot(d13, aes(gr, etiqueta)) +
  geom_point(aes(size = Count, fill = mlp), shape = 21, colour = "grey15", stroke = 0.3) +
  facet_col(vars(nombre_grupo), scales = "free_y", space = "free") +
  scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.06)),
                     breaks = scales::breaks_pretty(4)) +
  esc_count + esc_fill(LIM_MLP) +
  labs(x = bquote("Gene ratio ("*italic(k)*" of "*.(n_lista)*" "*italic("A. thaliana")*" homologs)"),
       y = NULL,
       caption = if (any(!d13$replica))
         bquote(.(DAGA)*" not significant with "*italic("E. globulus")*" genes as units")) +
  theme_nac() +
  theme(panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.25),
        strip.text = element_text(hjust = 0, margin = margin(1.5, 2, 1.5, 2)),
        plot.caption = element_text(size = 8, hjust = 0, colour = "grey25", margin = margin(t = 5)),
        plot.caption.position = "plot")
alto13 <- mm2in(nrow(d13) * 4.6 + nlevels(d13$nombre_grupo) * 5.2 + 22)

# --- 5. F13S: arbol de los 35 significativos -------------------------------------------------
# Los rotulos (de termino y de grupo) van como texto de eje a la derecha para que ggplot
# reserve su ancho: medidos a mano con systemfonts salian un 6 % cortos y se recortaban.
dt <- arb[order(arb$y), ]
dt$mlp <- mlp(dt$p.adjust); dt$rep <- dt$ID %in% bs$ID; dt$replica <- dt$ID %in% sig_egl$BP
dt$etiqueta <- ifelse(dt$replica, dt$Description, paste(dt$Description, DAGA))
lab_dt <- as.expression(lapply(seq_len(nrow(dt)), function(i)
  if (dt$rep[i]) bquote(bold(.(dt$etiqueta[i]))) else dt$etiqueta[i]))
YLIM <- c(0.4, nrow(dt) + 0.6)
br <- do.call(rbind, lapply(split(dt, dt$nombre_grupo), function(d)
  data.frame(nombre_grupo = d$nombre_grupo[1], y0 = min(d$y) - 0.32, y1 = max(d$y) + 0.32,
             ym = mean(range(d$y)))))
br$rotulo <- vapply(br$nombre_grupo, function(s) paste(strwrap(s, width = 17), collapse = "\n"), "")
eje_dcha <- function(...) theme_void(base_family = FAMILIA_FIG) +
  theme(axis.text.y.right = element_text(size = 8, colour = "grey30", hjust = 0, ...))

p_arbol <- suppressMessages(pt + scale_y_continuous(limits = YLIM, expand = c(0, 0))) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0))) + theme(plot.margin = margin(2, 0, 2, 2))
p_lab <- ggplot(dt, aes(0, y)) +
  geom_point(aes(size = Count, fill = mlp), shape = 21, colour = "grey15", stroke = 0.3) +
  scale_x_continuous(limits = c(-1, 1), expand = c(0, 0)) +
  scale_y_continuous(limits = YLIM, expand = c(0, 0), breaks = dt$y, labels = lab_dt,
                     position = "right") +
  esc_count + esc_fill(LIM_MLP, horizontal = TRUE) +
  coord_cartesian(clip = "off") +
  eje_dcha(margin = margin(l = 1.6, unit = "mm")) + theme(plot.margin = margin(2, 6, 2, 0))
p_br <- ggplot(br) +
  geom_segment(aes(x = 0, xend = 0, y = y0, yend = y1), linewidth = 0.4, colour = "grey30") +
  geom_segment(aes(x = 0, xend = -1, y = y0, yend = y0), linewidth = 0.4, colour = "grey30") +
  geom_segment(aes(x = 0, xend = -1, y = y1, yend = y1), linewidth = 0.4, colour = "grey30") +
  scale_x_continuous(limits = c(-1, 0), expand = c(0, 0)) +
  scale_y_continuous(limits = YLIM, expand = c(0, 0), breaks = br$ym, labels = br$rotulo,
                     position = "right") +
  eje_dcha(lineheight = 0.9, margin = margin(l = 1.4, unit = "mm")) + theme(plot.margin = margin(2, 2, 2, 0))
f13s <- (p_arbol | p_lab | p_br) +
  plot_layout(widths = unit(c(1, 5, 1.2), c("null", "mm", "mm")), guides = "collect") +
  plot_annotation(caption = bquote(bold("Bold")*": representative term after redundancy reduction;   "*
                                   .(DAGA)*" not significant with "*italic("E. globulus")*" genes as units"),
                  theme = theme(plot.caption = element_text(family = FAMILIA_FIG, size = 8, hjust = 0,
                                                            colour = "grey25"))) &
  theme(legend.position = "bottom", legend.box = "horizontal",
        legend.title = element_text(family = FAMILIA_FIG, size = 8, face = "bold"),
        legend.text = element_text(family = FAMILIA_FIG, size = 8))
alto13s <- mm2in(nrow(dt) * 4.15 + 30)

# --- 6. F14: compareCluster por subfamilia -----------------------------------------------
cc1  <- rf$compareCluster_subfamilia
niv  <- intersect(SUBFAM_NIVELES, names(cc1$entran))
stopifnot(length(niv) == length(cc1$entran), setequal(unique(cst$Cluster), niv))
nagi <- unlist(cc1$entran)[niv]
for (k in niv) stopifnot(sum(cst$Cluster == k) == cc1$probados[[k]],
                         sum(cst$significativo[cst$Cluster == k]) == cc1$significativos[[k]],
                         sum(css$Cluster == k) == cc1$tras_simplify[[k]])
ids_u <- unique(css$ID)
d14 <- cst[cst$significativo & cst$Count >= MIN_COUNT & cst$ID %in% ids_u, ]
stopifnot(all(paste(css$Cluster, css$ID) %in% paste(d14$Cluster, d14$ID)))
d14$x  <- factor(d14$Cluster, levels = niv)
d14$gr <- frac(d14$GeneRatio); d14$mlp <- mlp(d14$p.adjust)
d14$rep_aqui <- paste(d14$Cluster, d14$ID) %in% paste(css$Cluster, css$ID)
# orden en escalera: por la primera subfamilia (orden de la casa) donde sale el termino y,
# dentro de ella, por p.adjust
prim  <- tapply(as.integer(d14$x), d14$ID, min)
p_pri <- vapply(names(prim), function(i) min(d14$p.adjust[d14$ID == i & as.integer(d14$x) == prim[[i]]]), 0)
orden14 <- names(prim)[order(prim, p_pri)]
d14$y <- factor(d14$ID, levels = rev(orden14))
desc14 <- setNames(cst$Description[match(levels(d14$y), cst$ID)], levels(d14$y))
n_sig  <- tapply(cst$significativo, factor(cst$Cluster, levels = niv), sum)
vacias <- setdiff(niv, as.character(unique(d14$x)))
scm <- unlist(cc1$significativos_count_min)
stopifnot(setequal(vacias, names(scm)[scm == 0]))
tv <- data.frame(x = factor(vacias, levels = niv), y = (length(orden14) + 1) / 2,
                 lab = ifelse(n_sig[vacias] == 0, "no enriched term", "single-gene terms only"))
p14_tiras <- ggplot(data.frame(x = factor(niv, levels = niv)), aes(x, 1, fill = x)) +
  geom_tile(width = 0.82, height = 1, colour = NA) +
  scale_fill_manual(values = PAL_SUBFAM[niv], guide = "none") +
  scale_x_discrete(position = "top", labels = function(k) sprintf("%s\n(%d)", k, nagi[k])) +
  scale_y_continuous(expand = c(0, 0), breaks = NULL) +
  labs(x = expression("Subfamily (number of "*italic("A. thaliana")*" homologs)"), y = NULL) +
  theme_nac() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.ticks = element_blank(),
        axis.text.y = element_blank(), axis.text.x.top = element_text(lineheight = 0.9),
        plot.margin = margin(2, 3, 0, 2))
p14 <- ggplot(d14, aes(x, y)) +
  geom_point(aes(size = gr, fill = mlp), shape = 21, colour = "grey15", stroke = 0.3) +
  geom_text(data = tv, aes(x, y, label = lab), angle = 90, size = TAM_TXT, colour = "grey40",
            family = FAMILIA_FIG, inherit.aes = FALSE) +
  scale_x_discrete(drop = FALSE) +
  scale_y_discrete(labels = desc14) +
  scale_size_area(max_size = 3.7, limits = c(0, 1), breaks = c(0.25, 0.5, 0.75, 1),
                  name = "Gene ratio", guide = guide_legend(order = 2)) +
  esc_fill(range(d14$mlp)) +
  labs(x = NULL, y = NULL) +
  theme_nac() +
  theme(panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.25),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        plot.margin = margin(0, 3, 2, 2))
f14 <- (p14_tiras / p14) + plot_layout(heights = unit(c(2.4, 1), c("mm", "null")), guides = "collect")
alto14 <- mm2in(length(orden14) * 3.95 + 30)

# --- 7. guardar figuras -------------------------------------------------------------------
for (h in c(alto13, alto13s, alto14)) stopifnot(h <= ALTO_MAXIMO)
guardar_figura(p13,  "F13_GO_dotplot",     alto13,  ANCHO_COLUMNA)
guardar_figura(f13s, "F13S_GO_arbol",      alto13s, ANCHO_COLUMNA)
guardar_figura(f14,  "F14_GO_subfamilias", alto14,  ANCHO_COLUMNA)

# --- 8. tablas --------------------------------------------------------------------------
# T42: una fila por EglNAC (los 189). GO directa (sin ancestros) del homologo, con codigos de
# evidencia; los 48 sin homologo que pase los filtros quedan sin anotacion.
ad1 <- unique(ad[, c("TAIR", "GO", "ONTOLOGY", "termino", "EVIDENCE")])
pares <- aggregate(EVIDENCE ~ TAIR + GO + ONTOLOGY + termino, data = ad1,
                   FUN = function(e) paste(sort(unique(e)), collapse = ","))
stopifnot(nrow(pares) == rf$anotacion_directa$pares_AGI_GO,
          length(unique(pares$TAIR)) == rf$anotacion_directa$AGI)
pares <- pares[order(pares$TAIR, pares$ONTOLOGY, pares$GO), ]
pares$txt <- sprintf("%s %s [%s]", pares$GO, pares$termino, pares$EVIDENCE)
t42 <- t41[, c("simbolo", "locus", "subfamilia", "AGI", "ANAC", "nombre_Ath", "es_NAC_Ath", "identidad",
               "cobertura_consulta", "cobertura_diana", "evalue", "bitscore", "pasa_filtros")]
t42$pasa_filtros <- t42$pasa_filtros %in% c(TRUE, "True")
for (o in ONT) {
  po <- pares[pares$ONTOLOGY == o, ]
  t42[[paste0("n_GO_", o)]] <- ifelse(t42$pasa_filtros, vapply(t42$AGI, function(a) sum(po$TAIR == a), 0L), NA)
  t42[[paste0("GO_", o)]] <- ifelse(t42$pasa_filtros,
                                    vapply(t42$AGI, function(a) paste(po$txt[po$TAIR == a], collapse = "; "), ""), NA)
}
rep_de <- function(a, d) paste(d$Description[vapply(partir(d$geneID), function(g) a %in% g, NA)], collapse = "; ")
t42$terminos_BP_representativos <- ifelse(t42$pasa_filtros, vapply(t42$AGI, rep_de, "", d = bs), NA)
t42$n_terminos_BP_significativos <- ifelse(t42$pasa_filtros, vapply(t42$AGI, function(a)
  sum(vapply(partir(sig35$geneID), function(g) a %in% g, NA)), 0L), NA)
write_csv(t42, file.path(TABLAS, "T42_anotacion_GO_NAC189.csv"), na = "")

# T43: todos los terminos probados en BP, MF y CC: analisis principal (AGI) + contraste con
# genes de E. globulus + sensibilidad sin limite superior de tamano. Una fila por termino.
t43 <- do.call(rbind, lapply(ONT, function(o) {
  p <- eg[[o]]
  p <- data.frame(ontologia = o, ID = p$ID, termino = p$Description, GeneRatio = p$GeneRatio,
                  BgRatio = p$BgRatio, FoldEnrichment = p$FoldEnrichment, pvalue = p$pvalue,
                  p.adjust = p$p.adjust, qvalue = p$qvalue, Count = p$Count, significativo = signif_(p),
                  AGI = gsub("/", ";", p$geneID), EglNAC = vapply(partir(p$geneID), sim_de, ""))
  p$representativo <- if (o == "BP") p$ID %in% bs$ID else NA
  p$grupo_arbol    <- if (o == "BP") arb$nombre_grupo[match(p$ID, arb$ID)] else NA
  e <- ee[[o]]
  e <- data.frame(ID = e$ID, GeneRatio_Egl = e$GeneRatio, BgRatio_Egl = e$BgRatio,
                  p.adjust_Egl = e$p.adjust, Count_Egl = e$Count, significativo_Egl = signif_(e))
  s <- se[[o]]
  s <- data.frame(ID = s$ID, termino_s = s$Description, GeneRatio_sinMax = s$GeneRatio,
                  BgRatio_sinMax = s$BgRatio, p.adjust_sinMax = s$p.adjust,
                  significativo_sinMax = as.logical(s$significativo),
                  supera_maxGSSize = as.logical(s$supera_maxGSSize))
  m <- merge(merge(p, e, by = "ID", all = TRUE), s, by = "ID", all = TRUE)
  m$ontologia <- o
  m$termino[is.na(m$termino)] <- m$termino_s[is.na(m$termino)]
  m$termino_s <- NULL
  m[order(!m$significativo %in% TRUE, m$p.adjust, m$p.adjust_sinMax), ]
}))
stopifnot(!anyNA(t43$termino), !anyDuplicated(paste(t43$ontologia, t43$ID)))
write_csv(t43, file.path(TABLAS, "T43_enriquecimiento_GO.csv"), na = "")

# T44: compareCluster (BP) por subfamilia y por grupo, todos los terminos probados.
cc_tabla <- function(d, ds, nivel, col) {
  n <- cl[cl$nivel == nivel, ]
  data.frame(nivel = nivel, grupo = d$Cluster, n_AGI = n$n_AGI[match(d$Cluster, n$grupo)],
             n_Egl = n$n_Egl[match(d$Cluster, n$grupo)], ID = d$ID, termino = d$Description,
             GeneRatio = d$GeneRatio, BgRatio = d$BgRatio, FoldEnrichment = d$FoldEnrichment,
             pvalue = d$pvalue, p.adjust = d$p.adjust, qvalue = d$qvalue, Count = d$Count,
             significativo = as.logical(d$significativo),
             representativo = paste(d$Cluster, d$ID) %in% paste(ds$Cluster, ds$ID),
             en_figura = nivel == "subfamilia" & paste(d$Cluster, d$ID) %in% paste(d14$Cluster, d14$ID),
             AGI = gsub("/", ";", d$geneID),
             EglNAC = mapply(function(g, k) sim_de(g, pool = z[z[[col]] == k, ]), partir(d$geneID), d$Cluster))
}
t44 <- rbind(cc_tabla(cst, css, "subfamilia", "subfamilia"), cc_tabla(cgt, cgs, "grupo", "grupo"))
t44 <- t44[order(t44$nivel != "subfamilia", match(t44$grupo, c(niv, unique(cgt$Cluster))),
                 !t44$significativo, t44$p.adjust), ]
stopifnot(sum(t44$en_figura) == nrow(d14), !any(t44$EglNAC[t44$significativo] == ""))
write_csv(t44, file.path(TABLAS, "T44_compareCluster_GO.csv"), na = "")

# --- 9. cifras ------------------------------------------------------------------------------
F_ <- "datos/entrega3/go/"
registrar_cifra("go_lista_agi", rg$parametros$lista_AGI, paste0(F_, "resumen_go.json"),
                "AGI distintos entre los homologos de los 189 EglNAC que pasan los filtros")
registrar_cifra("go_lista_egl", rg$parametros$lista_Egl, "T41.pasa_filtros", "EglNAC con homologo")
registrar_cifra("go_universo_agi", rg$parametros$universo_AGI, paste0(F_, "resumen_go.json"),
                "AGI con homologo en el proteoma de E. globulus")
registrar_cifra("go_universo_egl", rg$parametros$universo_Egl, paste0(F_, "resumen_go.json"),
                "genes de E. globulus con homologo (contraste)")
for (o in ONT) {
  lo <- tolower(o)
  registrar_cifra(paste0("go_", lo, "_lista_anotada"), rg[[o]]$n_lista_anotada, paste0(F_, "resumen_go.json"))
  registrar_cifra(paste0("go_", lo, "_universo_anotado"), rg[[o]]$N_universo_anotado, paste0(F_, "resumen_go.json"))
  registrar_cifra(paste0("go_", lo, "_probados"), rg[[o]]$terminos_probados, sprintf("%senrichGO_%s_todos.tsv", F_, o))
  registrar_cifra(paste0("go_", lo, "_signif"), rg[[o]]$significativos, sprintf("%senrichGO_%s_todos.tsv", F_, o))
  registrar_cifra(paste0("go_", lo, "_signif_egl"), rg[[o]]$significativos_contraste_Egl,
                  sprintf("%senricherEgl_%s_todos.tsv", F_, o))
  registrar_cifra(paste0("go_", lo, "_compartidos_egl"), rg[[o]]$compartidos, paste0(F_, "resumen_go.json"),
                  "significativos con AGI que tambien lo son con genes de E. globulus")
  s <- rf[[paste0("sensibilidad_", o)]]
  registrar_cifra(paste0("go_", lo, "_sinmax_signif"), s$significativos, sprintf("%ssens_maxInf_%s_todos.tsv", F_, o))
  registrar_cifra(paste0("go_", lo, "_sinmax_mayores500"), s$significativos_mayores_500,
                  sprintf("%ssens_maxInf_%s_todos.tsv", F_, o))
  top <- se[[o]][signif_(se[[o]]) & se[[o]]$supera_maxGSSize %in% c(TRUE, "TRUE"), ]
  top <- top[order(top$p.adjust), ][1, ]
  registrar_cifra(paste0("go_", lo, "_sinmax_primero"),
                  sprintf("%s (%s), %s, adj. P = %s", top$Description, top$ID, top$GeneRatio,
                          formatC(top$p.adjust, format = "e", digits = 1)),
                  sprintf("%ssens_maxInf_%s_todos.tsv", F_, o))
}
registrar_cifra("go_bp_representativos", nrow(bs), paste0(F_, "enrichGO_BP_simplify.tsv"), "simplify Wang 0,7")
registrar_cifra("go_bp_representativos_egl", sum(d13$replica), "F13", "representativos que se mantienen con genes de E. globulus")
registrar_cifra("go_bp_no_replicados", paste(sig35$Description[!sig35$ID %in% sig_egl$BP], collapse = "; "),
                "F13S", "significativos con AGI que no lo son con genes de E. globulus")
registrar_cifra("go_bp_padj_rango", sprintf("%s-%s", formatC(min(sig35$p.adjust), format = "e", digits = 1),
                                            formatC(max(sig35$p.adjust), format = "e", digits = 1)), "F13S")
registrar_cifra("go_bp_fold_rango", sprintf("%.1f-%.1f", min(bs$FoldEnrichment), max(bs$FoldEnrichment)), "F13")
registrar_cifra("go_arbol_grupos", length(unique(arb$grupo_arbol)), paste0(F_, "arbol_BP_grupos.tsv"))
registrar_cifra("go_arbol_nombres", paste(sprintf("%s (%d)", ORDEN_GRUPOS, as.integer(table(arb$nombre_grupo)[ORDEN_GRUPOS])),
                                          collapse = "; "), "F13S")
registrar_cifra("go_arbol_similitud_mediana", rf$arbol$similitud_mediana_fuera_diagonal, paste0(F_, "resumen_go_figuras.json"))
registrar_cifra("go_cc_min_agi", rf$parametros$MIN_AGI, paste0(F_, "resumen_go_figuras.json"))
registrar_cifra("go_cc_min_count", MIN_COUNT, paste0(F_, "resumen_go_figuras.json"))
registrar_cifra("go_cc_subfam_probadas", length(niv), "F14", paste(niv, collapse = ", "))
registrar_cifra("go_cc_subfam_excluidas", length(cc1$fuera), paste0(F_, "cc_listas.tsv"),
                paste(names(cc1$fuera), collapse = ", "))
registrar_cifra("go_cc_subfam_con_terminos", length(niv) - length(vacias), "F14",
                paste(setdiff(niv, vacias), collapse = ", "))
registrar_cifra("go_cc_subfam_sin_terminos", paste(sprintf("%s (%s)", vacias, tv$lab), collapse = "; "), "F14")
registrar_cifra("go_cc_terminos_figura", length(orden14), "F14", "union de representativos por subfamilia")
registrar_cifra("go_cc_puntos_figura", nrow(d14), "F14", "significativos con Count >= 2 entre esos terminos")
registrar_cifra("go_cc_puntos_no_representativos", sum(!d14$rep_aqui), "F14",
                "puntos de terminos significativos que simplify no eligio en esa subfamilia")
registrar_cifra("go_anot_pares", rf$anotacion_directa$pares_AGI_GO, "T42", "pares AGI-GO directos")
registrar_cifra("go_anot_egl_con_go", sum(rowSums(t42[, paste0("n_GO_", ONT)], na.rm = TRUE) > 0), "T42")
for (k in names(rg$fuentes)) registrar_cifra(paste0("go_version_", gsub("[^A-Za-z0-9]+", "_", k)), rg$fuentes[[k]],
                                             paste0(F_, "resumen_go.json"))
registrar_cifra("go_version_GOSemSim", rf$versiones$GOSemSim, paste0(F_, "resumen_go_figuras.json"))
volcar_cifras("e3_go")

message(sprintf("  F13: %d terminos en %d grupos | F13S: %d hojas | F14: %d terminos x %d subfamilias, %d puntos",
                nrow(d13), nlevels(d13$nombre_grupo), nrow(dt), length(orden14), length(niv), nrow(d14)))
message(sprintf("  T42: %d filas | T43: %d filas | T44: %d filas", nrow(t42), nrow(t43), nrow(t44)))
sessionInfo()
