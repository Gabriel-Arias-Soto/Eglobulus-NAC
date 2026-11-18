# 01b_filogenia_cuerpo.R -> Fig. 2 del cuerpo (árbol de la familia), rediseño del 02-11-2026.
# Pedido de Gabriel: la versión anterior quedaba centrada con mucho margen, con las puntas apretadas y
# líneas finas. Cambios: el abanico ocupa todo el ancho (sin leyenda lateral de 19 colores), cada
# subfamilia se rotula con un arco exterior y su nombre (rotulado directo), ramas y puntas más gruesas.
# Reutiliza árbol, censo y color de rama por consenso de 01_filogenia.R (se ejecuta primero).
# Uso: Rscript scripts/01b_filogenia_cuerpo.R [longitudes|cladograma]
args <- commandArgs(trailingOnly = TRUE)
MODO <- if (length(args)) args[1] else "longitudes"
stopifnot(MODO %in% c("longitudes", "cladograma"))
source("scripts/00_estilo.R")
options(nac.f06_v2 = TRUE)
source("scripts/01_filogenia.R")          # define arb, meta, cons, cn, MIN_SOPORTE, COL_DESTACADO

LW_RAMA <- 0.42
q <- ggtree(arb, layout = "circular", linewidth = LW_RAMA, colour = "grey80",
            branch.length = if (MODO == "cladograma") "none" else "branch.length")
q <- q %<+% meta
q$data$soporte <- suppressWarnings(as.numeric(sub(".*/", "", q$data$label)))
q$data$soporte[q$data$isTip] <- NA_real_
q$data$sop_alto <- !is.na(q$data$soporte) & q$data$soporte >= MIN_SOPORTE
q$data$rama <- cons[q$data$node]
q$data$rama[q$data$isTip & q$data$especie == "grandis"] <- ETQ_ANCLA
q$data$rama <- factor(q$data$rama, levels = c(SUBFAM_NIVELES, SUBFAM_NO_EVAL, ETQ_ANCLA))
stopifnot(sum(q$data$sop_alto) == sum(p$data$sop_alto))      # mismos nodos que la versión anterior

R_MAX <- max(q$data$x, na.rm = TRUE)
# --- arcos de subfamilia: tramos consecutivos de puntas (en el orden del abanico) de una misma
#     subfamilia; las anclas de E. grandis intercaladas no cortan el tramo ------------------------
tp <- q$data[q$data$isTip, ] |> arrange(y)
tp$sf <- ifelse(tp$especie == "globulus", as.character(tp$subfamilia_final), NA)
tp$sf_rel <- tp$sf
for (i in seq_len(nrow(tp))) if (is.na(tp$sf_rel[i]) && i > 1) tp$sf_rel[i] <- tp$sf_rel[i - 1]
run <- cumsum(c(TRUE, tp$sf_rel[-1] != tp$sf_rel[-nrow(tp)]))
arcos <- tp |> mutate(run = run) |> filter(!is.na(sf_rel)) |>
  group_by(run, sf = sf_rel) |>
  summarise(y0 = min(y[!is.na(sf)]), y1 = max(y[!is.na(sf)]), n = sum(!is.na(sf)), .groups = "drop")
# cada subfamilia se rotula una vez, en su tramo mayor; los tramos sueltos quedan con arco y sin texto
etq <- arcos |> group_by(sf) |> slice_max(n, n = 1, with_ties = FALSE) |> ungroup() |>
  mutate(n_tot = unname(setNames(cn$n, as.character(cn$subfamilia_final))[sf]),
         texto = sprintf("%s (%d)", sf, n_tot))
stopifnot(sum(etq$n_tot) == n_glob, nrow(etq) == length(unique(arcos$sf)))
N_TIP <- nrow(tp)
R_ARC0 <- R_MAX * 1.025; R_ARC1 <- R_MAX * 1.06; R_TXT <- R_MAX * 1.08
ang_de <- function(y) approx(tp$y, tp$angle, y, rule = 2)$y
foc <- q$data |> filter(!is.na(es_focal), es_focal)
# todos los tramos de 3 o más genes llevan rótulo; el mayor, con el total de la subfamilia
etq <- arcos |> filter(n >= 3 | run %in% (arcos |> group_by(sf) |> slice_max(n, n = 1, with_ties = FALSE))$run) |>
  group_by(sf) |> mutate(mayor = n == max(n) & !duplicated(n == max(n) & TRUE)) |> ungroup() |>
  mutate(n_tot = unname(setNames(cn$n, as.character(cn$subfamilia_final))[sf]),
         texto = ifelse(mayor, sprintf("%s (%d)", sf, n_tot), sf),
         ym = (y0 + y1) / 2)
stopifnot(sum(etq$n_tot[etq$mayor]) == n_glob)   # 18-11-2026: 185 con EglNAC37/87 corregidos
# el rótulo de la subfamilia del gen focal se aparta del rótulo de EglNAC66
i_f <- which(etq$sf == GEN_FOCAL_SUBFAM & etq$y0 <= foc$y[1] & etq$y1 >= foc$y[1])
if (length(i_f) && abs(etq$ym[i_f] - foc$y[1]) < 3)
  etq$ym[i_f] <- if (foc$y[1] - etq$y0[i_f] > etq$y1[i_f] - foc$y[1]) etq$y0[i_f] + 1.5 else etq$y1[i_f] - 1.5
rad <- function(d) d |> mutate(a = ang_de(ym), ang_txt = ifelse(a > 90 & a < 270, a + 180, a),
                               hj = ifelse(a > 90 & a < 270, 1, 0))
etq <- rad(etq)
rot_f <- rad(tibble(ym = foc$y[1]))

q <- q +
  geom_tree(aes(colour = rama), linewidth = LW_RAMA) +
  geom_rect(data = arcos, aes(xmin = R_ARC0, xmax = R_ARC1, ymin = y0 - 0.45, ymax = y1 + 0.45, fill = sf),
            inherit.aes = FALSE, colour = NA) +
  geom_text(data = etq, aes(x = R_TXT, y = ym, label = texto, angle = ang_txt, hjust = hj),
            inherit.aes = FALSE, size = 2.6, family = FAMILIA_FIG, colour = "grey15") +
  geom_point2(aes(subset = sop_alto), size = 0.75, colour = "grey20", stroke = 0) +
  geom_tippoint(aes(colour = rama), size = 1.05, stroke = 0) +
  geom_point2(aes(subset = es_focal), shape = 21, size = 2.6, fill = COL_DESTACADO, colour = "white",
              stroke = 0.4) +

  scale_colour_manual(values = c(PAL_SUBFAM, setNames("grey62", ETQ_ANCLA)), na.value = "grey80",
                      guide = "none") +
  scale_fill_manual(values = PAL_SUBFAM, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0)), limits = c(0, R_MAX * 1.24)) +
  theme_void(base_family = FAMILIA_FIG) + theme(plot.margin = margin(-27, -24, -22, -24, "mm"))
if (MODO == "longitudes") q <- q + geom_treescale(x = R_MAX * 0.02, y = N_TIP * 0.5, width = 0.5,
                                                  fontsize = 2.5, linesize = 0.4, offset = 3,
                                                  family = FAMILIA_FIG)

# rótulo del gen focal: radial, en el anillo de rótulos, en rojo
q <- q + geom_text(data = rot_f, aes(x = R_TXT, y = ym, label = GEN_FOCAL, angle = ang_txt, hjust = hj),
                   inherit.aes = FALSE, fontface = "bold", size = 2.6, colour = COL_DESTACADO,
                   family = FAMILIA_FIG)

# clave mínima (lo que no es subfamilia) en la esquina inferior izquierda, como texto aparte
clave <- ggplot() + theme_void(base_family = FAMILIA_FIG) +
  annotate("point", x = 0, y = c(3, 2, 1), size = c(1.4, 1.05, 2.6), shape = c(16, 16, 21),
           colour = c("grey62", "grey20", "white"), fill = c(NA, NA, COL_DESTACADO), stroke = 0.4) +
  annotate("text", x = 0.6, y = c(3, 2, 1), hjust = 0, size = 2.5, family = FAMILIA_FIG,
           label = c("E. grandis anchor", sprintf("node with UFBoot \u2265 %d", MIN_SOPORTE),
                     "EglNAC66"), fontface = c("italic", "plain", "plain")) +
  coord_cartesian(xlim = c(-0.5, 9), ylim = c(0.4, 3.6))
fig <- q + patchwork::inset_element(clave, left = 0, bottom = 0, right = 0.30, top = 0.09, align_to = "full")

nombre <- if (MODO == "longitudes") "F06_filogenia_NAC" else "F06_filogenia_NAC_cladograma"
guardar_figura(fig, nombre, alto = 150 * MM, ancho = 150 * MM)
