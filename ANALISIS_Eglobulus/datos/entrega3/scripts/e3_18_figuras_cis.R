# e3_18_figuras_cis.R -> figuras/F18_cis_PlantCARE_por_gen, F19_cis_enriquecimiento_fondo_GC
#                        tablas/_cifras/e3_cis.csv
#
# Entrega 3, figuras de elementos cis. No recalcula nada: dibuja las tablas que dejan
# e3_15_contraste_jaspar.py (T51, T52), e3_16_lectura_plantcare.py (T54-T56, T59) y
# e3_17_contraste_plantcare.py (T57, T58), mas las sumas del fondo emparejado por GC
# (fondo_sumas_*.npz, 10 000 conjuntos de 189 promotores) para el intervalo del azar.
#
# Antes de dibujar se comprueba (si algo no cuadra, el script para):
#   - que T56 tiene los 189 genes de T01, con la misma subfamilia, y que sus totales por
#     clase son la suma de las 70 unidades de cada clase (y coinciden con T54);
#   - que el orden de filas sale del arbol de F06 (ladderize, como 01_filogenia.R):
#     183 puntas GL_, y los 6 genes sin arbol son los de dominio parcial;
#   - que las medias de las matrices del fondo coinciden con media_fondo de T52/T57/T58
#     (asi se sabe que las columnas del .npz se emparejan bien con las filas de la tabla).
#
# Decisiones de dibujo
#   F18  replica la Fig. 6 del articulo 1 con 189 genes: mapa de recuentos por gen y
#        unidad, columnas agrupadas por clase (A-E, mismo orden que T58) y barras apiladas
#        por clase a la derecha. Filas: facetas de subfamilia como F07 y, dentro, orden de
#        las puntas de F06. Columnas: las 39 unidades con sitio en >= 10 % de los promotores
#        (19 de 189); con las 70 las columnas no admiten rotulos de 8 pt en 174 mm. Las
#        barras suman las 70 unidades, asi que no se pierde ningun sitio.
#   F19  observado / esperado frente al fondo emparejado por GC, en escala log2. La banda
#        gris es el intervalo central del 95 % de las 10 000 replicas del fondo dividido por
#        su media: lo que da el azar con ese GC. Relleno = q < 0,05 y si se mantiene en todas
#        las variantes del contraste (T57: 5, T58: 3, T52: 2).
#
# Se ejecuta desde la RAIZ de "Revision NAC":
#   Rscript ANALISIS_Eglobulus/datos/entrega3/scripts/e3_18_figuras_cis.R \
#     > ANALISIS_Eglobulus/logs/entrega3_20_figuras_cis.log 2>&1
if (!exists("ANCHO_COLUMNA")) source("scripts/00_estilo.R")
suppressPackageStartupMessages(library(patchwork))

E3   <- file.path(ANALISIS, "datos", "entrega3")
# Rotulacion: theme_nac() con su base de 9 deja la figura en 8-12 pt, el rango de Trees
# (8 ejes y leyendas, 8,5 facetas, 9 titulos de eje, 12 letra de panel), como F13/F14.
# Con theme_nac(8) los ejes x y las leyendas salian a 7 pt (auditado en el PDF). La unica
# excepcion son las facetas de subfamilia de F18, a 5 pt como en F07 (189 filas).
CLASES <- c("Luz" = "A", "Desarrollo" = "B", "Fitohormonas" = "C", "Estr\u00e9s" = "D", "Otros" = "E")
PAL_CLASE <- c(A = "#E69F00", B = "#009E73", C = "#0072B2", D = "#D55E00", E = "#999999")  # Okabe-Ito
MIN_GENES <- 19                                   # 10 % de 189
CORTES    <- c(-Inf, 0, 1, 2, 3, 5, Inf)
ETQ_CORTE <- c("0", "1", "2", "3", "4\u20135", "\u2265 6")
PAL_CORTE <- setNames(c("#F2F2F2", "#C7E9B4", "#7FCDBB", "#41B6C4", "#2C7FB8", "#253494"), ETQ_CORTE)
SIG <- c(rob = "q < 0.05 in all variants", pri = "q < 0.05, not in all variants", ns = "q \u2265 0.05")
PAL_SIG <- setNames(c(COL_DESTACADO, "#F4A582", "white"), SIG)
leer <- function(n) suppressMessages(read_csv(file.path(TABLAS, n), show_col_types = FALSE))

# --- genes y orden de filas ---------------------------------------------------
cen <- leer_censo()
t56 <- leer("T56_recuentos_PlantCARE_por_gen_NAC189.csv")
t54 <- leer("T54_elementos_cis_PlantCARE_clases.csv")
stopifnot(nrow(t56) == 189, setequal(t56$locus, cen$locus))
chk <- left_join(t56[, c("locus", "simbolo", "subfamilia")], cen, by = "locus", suffix = c("", ".t01"))
stopifnot(all(chk$simbolo == chk$simbolo.t01),
          all(subfam_facetas(chk$subfamilia) == subfam_facetas(chk$subfamilia_final)))

arb   <- ape::ladderize(ape::read.tree(file.path(RAIZ, "arbol_NAC_189.treefile")))
punta <- arb$tip.label[arb$edge[arb$edge[, 2] <= ape::Ntip(arb), 2]]
orden <- sub("^GL_", "", grep("^GL_", punta, value = TRUE))
# 18-11-2026: 185 puntas (EglNAC37/87 con modelo corregido); el control sale del censo
stopifnot(length(orden) == sum(cen$subfamilia_final != SUBFAM_NO_EVAL), !anyDuplicated(orden))

gen <- cen |>
  transmute(locus, simbolo, sub = subfam_facetas(subfamilia_final), pos = match(locus, orden)) |>
  arrange(sub, pos)
stopifnot(sum(is.na(gen$pos)) == sum(cen$subfamilia_final == SUBFAM_NO_EVAL), all(gen$sub[is.na(gen$pos)] == ETQ_PARCIAL))
gen <- mutate(gen, fila = factor(simbolo, levels = rev(simbolo)))   # primera punta arriba

# --- unidades -----------------------------------------------------------------
u5 <- t54 |> filter(clase %in% names(CLASES)) |> mutate(letra = unname(CLASES[clase]))
stopifnot(nrow(u5) == 70, all(u5$unidad %in% names(t56)))
stopifnot(all(colSums(t56[, u5$unidad]) == u5$sitios_NAC189))
for (cl in names(CLASES))
  stopifnot(all(t56[[cl]] == rowSums(t56[, u5$unidad[u5$clase == cl]])))
stopifnot(all(t56$total_5_clases == rowSums(t56[, names(CLASES)])))
pres  <- colSums(t56[, u5$unidad] > 0)
stopifnot(all(pres == u5$genes_con_sitio))
u_fig <- u5 |> mutate(genes = unname(pres[unidad])) |> filter(genes >= MIN_GENES) |>
  arrange(letra)                                    # dentro de clase, orden de T54
message(sprintf("unidades en F18: %d de %d (con sitio en >= %d promotores)",
                nrow(u_fig), nrow(u5), MIN_GENES))
message("  por clase: ", paste(names(table(u_fig$letra)), table(u_fig$letra), sep = "=", collapse = " "))
message("  columnas de T56 que no son unidad ni clase: ",
        paste(setdiff(names(t56), c(t54$unidad, "locus", "simbolo", "subfamilia",
                                    names(CLASES), "total_5_clases")), collapse = ", "))
etq_clase <- u5 |> distinct(letra, clase_en) |> arrange(letra)
stopifnot(nrow(etq_clase) == 5)
ETQ_CLASE <- setNames(paste(etq_clase$letra, etq_clase$clase_en), etq_clase$letra)
# Rotulo de las unidades en F18 y F19: el guion bajo de PlantCARE (GCN4_motif) se escribe como
# espacio. En F18 las etiquetas van giradas con 7,5 pt de paso y el guion bajo de GCN4_motif
# caia junto a la O de la etiqueta vecina, que se leia "Ō2-site". Las tablas guardan el nombre
# de PlantCARE.
etq_unidad <- function(x) gsub("_", " ", x, fixed = TRUE)

# --- F18: mapa de recuentos + barras por clase ----------------------------------
m <- t56 |> select(locus, all_of(u_fig$unidad)) |>
  pivot_longer(-locus, names_to = "unidad", values_to = "n") |>
  left_join(gen, by = "locus") |>
  left_join(select(u_fig, unidad, letra), by = "unidad") |>
  mutate(unidad = factor(unidad, levels = u_fig$unidad),
         nivel  = cut(n, CORTES, labels = ETQ_CORTE))
stopifnot(nrow(m) == 189 * nrow(u_fig), !anyNA(m$nivel), !anyNA(m$fila))
focal <- filter(gen, simbolo == GEN_FOCAL)
stopifnot(nrow(focal) == 1)
marco <- u_fig |> count(letra, name = "k") |>
  mutate(sub = focal$sub, fila = focal$fila, x = (k + 1) / 2)

tema_filas <- theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                    panel.grid = element_blank(), panel.spacing.y = unit(1.8, "pt"),
                    strip.background = element_blank())

p_mapa <- ggplot(m, aes(x = unidad, y = fila, fill = nivel)) +
  geom_tile() +
  geom_tile(data = marco, aes(x = x, y = fila, width = k), height = 1, inherit.aes = FALSE,
            fill = NA, colour = COL_DESTACADO, linewidth = 0.35) +
  facet_grid(sub ~ letra, scales = "free", space = "free", switch = "y") +
  scale_fill_manual(values = PAL_CORTE, name = "Sites per unit", drop = FALSE) +
  scale_x_discrete(expand = c(0, 0), labels = etq_unidad) + scale_y_discrete(expand = c(0, 0)) +
  labs(x = NULL, y = NULL) +
  guides(fill = guide_legend(nrow = 1, keywidth = unit(3, "mm"), keyheight = unit(2.6, "mm"))) +
  theme_nac() + tema_filas +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
        panel.spacing.x = unit(2, "pt"), strip.placement = "outside",
        strip.text.x = element_text(size = 8, face = "bold"),
        # las subfamilias de 1-2 genes (VIb) tienen menos alto que la letra: sin esto el
        # rotulo sale recortado ("VIh")
        strip.clip = "off",
        strip.text.y.left = element_text(angle = 0, size = 5, hjust = 1, face = "plain"))

b <- t56 |> select(locus, all_of(names(CLASES))) |>
  pivot_longer(-locus, names_to = "clase", values_to = "n") |>
  mutate(letra = factor(unname(CLASES[clase]), levels = names(PAL_CLASE))) |>
  left_join(gen, by = "locus")
tot <- t56 |> select(locus, total = total_5_clases) |> left_join(gen, by = "locus")
# EglNAC66 tiene la barra mas larga (68 sitios), asi que su rotulo no cabe dentro del panel.
# Va como texto del eje y derecho (vacio en las demas filas): queda junto al final de su barra,
# fuera del panel, y el panel conserva el recorte. Con coord_cartesian(clip = "off") el marco
# del panel se dibuja entero y sale con el doble de grosor que el de los paneles del mapa.
stopifnot(tot$simbolo[which.max(tot$total)] == GEN_FOCAL)

p_barras <- ggplot(b, aes(x = n, y = fila, fill = letra)) +
  geom_col(width = 0.85, position = position_stack(reverse = TRUE)) +
  facet_grid(sub ~ ., scales = "free_y", space = "free_y") +
  scale_fill_manual(values = PAL_CLASE, labels = ETQ_CLASE, name = "Class") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.03)), breaks = c(0, 20, 40, 60)) +
  scale_y_discrete(position = "right", labels = function(x) ifelse(x == GEN_FOCAL, x, "")) +
  labs(x = "Sites per promoter", y = NULL) +
  guides(fill = guide_legend(ncol = 2, keywidth = unit(3, "mm"), keyheight = unit(2.6, "mm"))) +
  theme_nac() + tema_filas +
  theme(strip.text.y = element_blank(), panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.25),
        axis.text.y.right = element_text(size = 8, face = "bold", colour = COL_DESTACADO, hjust = 0,
                                         margin = margin(l = 1.5)))

f18 <- (p_mapa | p_barras) + plot_layout(widths = c(nrow(u_fig), 12), guides = "collect") &
  theme(legend.position = "bottom", legend.box = "horizontal", legend.margin = margin(0, 4, 0, 4))
guardar_figura(f18, "F18_cis_PlantCARE_por_gen", alto = ALTO_PANELES, ancho = ANCHO_COLUMNA)

registrar_cifra("f18_unidades_mostradas", nrow(u_fig), "T56", sprintf("con sitio en >= %d promotores", MIN_GENES))
registrar_cifra("f18_unidades_5_clases", nrow(u5), "T54")
registrar_cifra("f18_sitios_5_clases", sum(t56$total_5_clases), "T56")
registrar_cifra("f18_sitios_por_promotor_rango", sprintf("%d-%d", min(t56$total_5_clases), max(t56$total_5_clases)), "T56")
registrar_cifra("f18_sitios_por_promotor_mediana", median(t56$total_5_clases), "T56")
registrar_cifra("f18_sitios_focal", tot$total[tot$simbolo == GEN_FOCAL], "T56", GEN_FOCAL)
registrar_cifra("f18_rango_focal", sum(t56$total_5_clases > tot$total[tot$simbolo == GEN_FOCAL]) + 1, "T56",
                "puesto por sitios en las 5 clases (1 = mas sitios)")

# --- F19: observado / esperado frente al fondo ----------------------------------
leer_npy <- function(npz, nombre) {
  con <- unz(npz, paste0(nombre, ".npy"), open = "rb"); on.exit(close(con))
  stopifnot(identical(readBin(con, "raw", 6), as.raw(c(0x93, 0x4e, 0x55, 0x4d, 0x50, 0x59))))
  ver <- as.integer(readBin(con, "raw", 2))
  hl  <- if (ver[1] == 1L) readBin(con, "integer", 1, size = 2, signed = FALSE, endian = "little") else
    readBin(con, "integer", 1, size = 4, endian = "little")
  cab   <- rawToChar(readBin(con, "raw", hl))
  descr <- sub(".*'descr': *'([^']+)'.*", "\\1", cab)
  stopifnot(!grepl("'fortran_order': *True", cab))
  forma <- as.integer(strsplit(gsub("[^0-9,]", "", sub(".*'shape': *\\(([^)]*)\\).*", "\\1", cab)), ",")[[1]])
  n <- prod(forma)
  if (grepl("^<U[0-9]+$", descr)) {
    w <- as.integer(sub("<U", "", descr))
    v <- readBin(con, "integer", n * w, size = 4, endian = "little"); stopifnot(length(v) == n * w)
    return(apply(matrix(v, nrow = w), 2, function(z) intToUtf8(z[z > 0])))
  }
  stopifnot(descr == "<i4")
  v <- readBin(con, "integer", n, size = 4, endian = "little"); stopifnot(length(v) == n)
  if (length(forma) == 2) matrix(v, nrow = forma[1], ncol = forma[2], byrow = TRUE) else v
}
intervalo <- function(npz, matriz, nombres, tabla, clave, obs) {
  M <- leer_npy(npz, matriz); colnames(M) <- leer_npy(npz, nombres)
  stopifnot(nrow(M) == 10000, setequal(colnames(M), tabla[[clave]]))
  # sin esperado no hay razon (T52: ARID, VOZ y "sin anotacion" no tienen ninguna region
  # ni en los NAC ni en el fondo); se quedan fuera y se dice en el log
  fuera <- tabla[[clave]][!(tabla$media_fondo > 0)]
  if (length(fuera)) {
    stopifnot(all(M[, fuera] == 0), all(tabla[[obs]][tabla[[clave]] %in% fuera] == 0))
    message(sprintf("  %s: sin esperado (0 en NAC y en las 10 000 replicas), fuera de la figura: %s",
                    matriz, paste(fuera, collapse = ", ")))
  }
  tabla <- tabla[tabla$media_fondo > 0, ]
  M  <- M[, tabla[[clave]], drop = FALSE]
  mu <- colMeans(M)
  stopifnot(all(abs(mu - tabla$media_fondo) <= 0.01 + 1e-5 * mu),
            all(abs(tabla[[obs]] / tabla$media_fondo - tabla$razon) <= 1e-4 * tabla$razon))
  qq <- apply(M, 2, quantile, probs = c(0.025, 0.975), names = FALSE)   # no "q": es columna de la tabla
  n_var <- sum(grepl("_q$", names(tabla)))
  tabla |> mutate(lo = qq[1, ] / mu, hi = qq[2, ] / mu,
                  sig = factor(ifelse(q < 0.05 & !is.na(variantes_concordantes) & variantes_concordantes == n_var,
                                      SIG["rob"], ifelse(q < 0.05, SIG["pri"], SIG["ns"])), levels = SIG))
}
NPZ_P <- file.path(E3, "cis_plantcare", "fondo_sumas_PlantCARE.npz")
NPZ_J <- file.path(E3, "cis_jaspar", "fondo_sumas_JASPAR2026.npz")
t58 <- leer("T58_enriquecimiento_PlantCARE_clases.csv")
t57 <- leer("T57_enriquecimiento_PlantCARE_elementos.csv")
t52 <- leer("T52_enriquecimiento_JASPAR_familias.csv")

dc <- intervalo(NPZ_P, "c_sitios", "grupos", t58, "clase", "sitios_NAC") |>
  mutate(letra = unname(CLASES[clase]),
         etq = ifelse(is.na(letra), clase_en, paste(letra, clase_en)),
         etq = factor(etq, levels = rev(etq[order(is.na(letra), letra)])))
du <- intervalo(NPZ_P, "u_sitios", "unidades", t57, "unidad", "sitios_NAC") |>
  filter(unidad %in% u_fig$unidad) |>
  mutate(letra = unname(CLASES[clase]), etq = factor(etq_unidad(unidad), levels = rev(etq_unidad(u_fig$unidad))))
stopifnot(nrow(du) == nrow(u_fig), all(du$letra == u_fig$letra[match(du$unidad, u_fig$unidad)]))
# T52 guarda los nombres armonizados en espanol ("bZIP grupo D", "AP2/EREBP sin subfamilia");
# la figura va en ingles
etq_familia <- function(x) {
  x <- sub("^bZIP grupo ([A-Z])$", "bZIP group \\1", x)
  x <- sub("^(.+) sin (grupo|subfamilia)$", "\\1 (unassigned)", x)
  x <- gsub("_", " ", x)
  stopifnot(!any(grepl("grupo|subfamilia|sin |anotaci", x)), !anyDuplicated(x))
  x
}
dj <- intervalo(NPZ_J, "fam_regiones", "familias", t52, "familia", "regiones_NAC") |>
  mutate(etq = etq_familia(familia), etq = factor(etq, levels = etq[order(razon)]))
# el eje es log2: un 0 en la razon o en el borde de la banda no se puede dibujar. Se mira solo
# en lo que sale en la figura (21 de las 123 unidades tienen 0 en el percentil 2,5 del fondo,
# todas raras y fuera de F18/F19)
for (d in list(dc, du, dj)) stopifnot(all(d$lo > 0), all(d$razon > 0), all(is.finite(d$hi)))

# show.legend = TRUE: desde ggplot2 3.5 cada panel solo dibuja la clave de los niveles que
# tiene (a: ninguno significativo; b: el robusto; c: los no robustos), las tres leyendas
# difieren y patchwork no las funde. Forzando todas las claves son identicas y sale una.
# Cortes por panel: a abarca 0,81-1,15, b 0,45-4,9 y c 0,64-1,41 (razon y banda juntas); con
# los de b en a y c solo se rotularia el 1. En c se extiende el eje a 0,6 para que la banda
# mas ancha (LEAFY, bZIP grupo K) tenga un corte a su izquierda.
panel_enr <- function(d, xlab, cortes) {
  ggplot(d, aes(y = etq)) +
    geom_vline(xintercept = 1, colour = "grey55", linewidth = 0.3) +
    geom_linerange(aes(xmin = lo, xmax = hi), colour = "grey80", linewidth = 1.5) +
    geom_point(aes(x = razon, fill = sig), shape = 21, size = 1.7, stroke = 0.35, colour = "grey15",
               show.legend = TRUE) +
    scale_x_continuous(transform = "log2", breaks = cortes, labels = function(x) sprintf("%g", x)) +
    scale_fill_manual(values = PAL_SIG, drop = FALSE, name = NULL) +
    labs(x = xlab, y = NULL) +
    theme_nac() +
    theme(panel.grid.major.y = element_blank(), axis.text.y = element_text(size = 8))
}
ETQ_X_SIT <- "Observed / expected sites"
pa <- panel_enr(dc, ETQ_X_SIT, c(0.8, 0.9, 1, 1.1))
pb <- panel_enr(du, ETQ_X_SIT, c(0.5, 1, 2, 4)) +
  facet_grid(letra ~ ., scales = "free_y", space = "free_y") +
  theme(strip.text.y = element_text(angle = 0, size = 8, face = "bold"), panel.spacing.y = unit(2, "pt"))
# el titulo del eje de c es mas ancho que el panel y, centrado, llegaba a 0,4 pt del borde de
# la figura; el margen derecho de 12 pt (3 en theme_nac) lo deja a unos 7 pt
pc <- panel_enr(dj, "Observed / expected motif regions", c(0.6, 0.8, 1, 1.2, 1.4)) +
  expand_limits(x = 0.6) +
  theme(plot.margin = margin(2, 12, 2, 2))
izq <- pa / pb + plot_layout(heights = c(nrow(dc), nrow(du) + 3))
f19 <- (izq | pc) + plot_layout(widths = c(1, 1), guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(legend.position = "bottom", legend.margin = margin(0, 0, 0, 0))
guardar_figura(f19, "F19_cis_enriquecimiento_fondo_GC", alto = 7.8, ancho = ANCHO_COLUMNA)

for (nm in c("dc", "du", "dj")) {
  d <- get(nm)
  message(sprintf("%s: n=%d | razon %.3f-%.3f | banda %.3f-%.3f | %s", nm, nrow(d), min(d$razon), max(d$razon),
                  min(d$lo), max(d$hi), paste(names(table(d$sig)), table(d$sig), sep = "=", collapse = "; ")))
  registrar_cifra(sprintf("f19_%s_q005", nm), sum(d$q < 0.05), c(dc = "T58", du = "T57", dj = "T52")[[nm]],
                  "filas con q < 0,05 en el contraste principal")
  # la banda es el 95 % central del fondo sin corregir (p empirica bilateral < 0,05 fuera de
  # ella); q es BH dentro de cada nivel. Lo que sale fuera de la banda con q >= 0,05 va a la
  # leyenda para que el punto blanco fuera de la banda no se lea como contradiccion
  fb <- d$razon < d$lo | d$razon > d$hi
  registrar_cifra(sprintf("f19_%s_fuera_banda", nm), sprintf("%d de %d", sum(fb), nrow(d)),
                  c(dc = "T58", du = "T57", dj = "T52")[[nm]],
                  paste("con q >= 0,05:", paste(as.character(d$etq)[fb & d$q >= 0.05], collapse = ", ")))
}
for (k in which(dj$q < 0.05))
  registrar_cifra(paste0("f19_familia_", dj$familia[k]), sprintf("%.2f (q = %.3f)", dj$razon[k], dj$q[k]), "T52",
                  sprintf("variantes concordantes: %d de 2", dj$variantes_concordantes[k]))
registrar_cifra("f19_familias_mostradas", sprintf("%d de %d", nrow(dj), nrow(t52)), "T52",
                paste("sin esperado:", paste(setdiff(t52$familia, dj$familia), collapse = ", ")))
k <- which(du$unidad == "RY-element")
registrar_cifra("f19_RY_element", sprintf("%.2f (q = %.3g; banda %.2f-%.2f)", du$razon[k], du$q[k], du$lo[k], du$hi[k]),
                "T57 + fondo_sumas_PlantCARE.npz", sprintf("variantes concordantes: %d de 5", du$variantes_concordantes[k]))
volcar_cifras("e3_cis")
