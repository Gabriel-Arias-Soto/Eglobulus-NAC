# 05_cromosomas.R -> panel de distribucion cromosomica
# Correcciones de la revision:
#  - el titulo del eje situaba los conteos "abajo" cuando estan arriba de cada barra.
#  - el rotulo del sombreado daba 125 genes en tandem, que es la cifra de
#    T01.bloque_tandem, mientras el sombreado se dibuja desde T15, que contiene
#    119. Son dos criterios distintos de tandem; la figura declara el que usa y
#    el aviso de abajo recuerda la diferencia en cada regeneracion.
#  - se rotulaban 189 genes por usar la columna 'cuantificable': ahora solo los
#    que tienen qRT-PCR en T10.
# Ninguna cifra se teclea: todas salen de registrar_cifra().
if (!exists("ANCHO_COLUMNA")) source("scripts/00_estilo.R")
suppressPackageStartupMessages({library(ggrepel); library(GenomicRanges)
                                library(patchwork)})

censo <- leer_censo()
gm <- mapa_cromosomas()
t15 <- suppressMessages(read_csv(file.path(TABLAS, "T15_bloques_tandem.csv"), show_col_types = FALSE))
# los genes rotulados son los que tienen qRT-PCR en T10, no los "cuantificables":
# 'cuantificable' es TRUE para casi todo el censo y rotularlos a todos hacia
# ilegible la figura -189 etiquetas solapadas-.
# VALIDADOS lo define 00_estilo.R a partir de T02.medido_por_qPCR

CROMS <- paste0("chr", 1:11)
# Dos filas de cromosomas, no once columnas en una. Con once en el ancho de
# pagina cada columna mide 0,57 in y el rotulo de un bloque necesita 0,45 in,
# de modo que el texto se montaba sobre el ideograma vecino. Repartidos en dos
# filas la columna pasa a 0,87 in y los rotulos caben en el hueco.
POR_FILA <- 6L
IDX  <- setNames(seq_along(CROMS), CROMS)
FILA <- setNames(((IDX - 1L) %/% POR_FILA) + 1L, CROMS)
XR   <- setNames(((IDX - 1L) %%  POR_FILA) + 1L, CROMS)

largos <- gm$largos |>
  filter(cromosoma %in% CROMS) |>
  mutate(x = XR[cromosoma], fila = FILA[cromosoma], largo_mb = largo / 1e6)

# los genes NAC se situan por su punto medio, en Mb
nac <- censo |>
  filter(cromosoma %in% CROMS) |>
  mutate(x = XR[cromosoma], fila = FILA[cromosoma],
         medio_mb = ((inicio + fin) / 2) / 1e6,
         valida = simbolo %in% VALIDADOS)

gr_bl <- GRanges(as.character(t15$cromosoma), IRanges(t15$inicio, t15$fin))
gr_nac <- GRanges(as.character(censo$cromosoma), IRanges(censo$inicio, censo$fin),
                  simbolo = censo$simbolo)
ov <- findOverlaps(gr_bl, gr_nac)
en_bloque <- unique(subjectHits(ov))

# --- que EglNAC forma cada bloque ------------------------------------------
# La figura nombra los miembros de cada bloque, no solo su posicion. Los miembros salen del mismo solape de coordenada que usa
# F09, asi que las dos figuras nombran exactamente el mismo conjunto de genes.
# Junto al ideograma va SOLO el codigo; la pertenencia entera va en la clave al
# pie. Al lado del cromosoma hay 0,87 in y el listado de B01 mide 40 caracteres,
# de modo que ahi solo caben dos rachas: 8 de los 21 bloques saldrian
# truncados. En la clave caben los 21
# completos, y el codigo del ideograma es el que enlaza con ella.
miembros <- tibble(bloque  = t15$bloque[queryHits(ov)],
                   simbolo = gr_nac$simbolo[subjectHits(ov)]) |>
  mutate(num = as.integer(sub("^EglNAC", "", simbolo))) |>
  group_by(bloque) |>
  summarise(n_nac = n(), ids = rango_compacto(num, 99L), .groups = "drop")
stopifnot(nrow(miembros) == nrow(t15))

# bloques en tandem como intervalos reales (GenomicRanges) proyectados a Mb
bl <- t15 |>
  filter(cromosoma %in% CROMS) |>
  left_join(miembros, by = "bloque") |>
  mutate(x = XR[cromosoma], fila = FILA[cromosoma],
         ini_mb = inicio / 1e6, fin_mb = fin / 1e6,
         medio_mb = (ini_mb + fin_mb) / 2,
         etiqueta = bloque)
stopifnot(!any(is.na(bl$ids)))

conteo <- nac |> count(cromosoma, x, fila, name = "n") |>
  left_join(largos |> select(cromosoma, largo_mb), by = "cromosoma")

stopifnot(sum(conteo$n) == 189)
registrar_cifra("f01_total", sum(conteo$n), "T01")
registrar_cifra("f01_cromosomas", nrow(largos), "genes_todos.tsv: scaffolds CM*")
registrar_cifra("f01_conteos",
                paste(sprintf("%s %d", conteo$cromosoma, conteo$n), collapse = ", "), "T01")
registrar_cifra("f01_max",
                sprintf("%s (%d)", conteo$cromosoma[which.max(conteo$n)], max(conteo$n)), "T01")
registrar_cifra("f01_min",
                sprintf("%s (%d)", conteo$cromosoma[which.min(conteo$n)], min(conteo$n)), "T01")
registrar_cifra("f01_bloques", nrow(t15), "T15_bloques_tandem.csv")
# La misma cadena que dibuja la clave al pie, para que el texto del manuscrito
# y la figura no puedan divergir. Ya no hay bloques abreviados.
registrar_cifra("f01_bloques_miembros",
                paste(sprintf("%s (%s): EglNAC %s", bl$bloque, bl$cromosoma,
                              bl$ids), collapse = "; "),
                "solape de coordenada de T15 con T01")
stopifnot(nrow(bl) == nrow(t15))  # ningun bloque se pierde al filtrar cromosomas
registrar_cifra("f01_genes_tandem", sum(t15$n_genes), "T15.n_genes")
registrar_cifra("f01_genes_tandem_gr", length(en_bloque),
                "GenomicRanges: NAC que solapan un intervalo de bloque")
registrar_cifra("f01_validados", sum(nac$valida), "T10: genes con qRT-PCR")
registrar_cifra("f01_validados_cuales", paste(sort(VALIDADOS), collapse = ", "), "T10")
n_t01 <- sum(!is.na(censo$bloque_tandem))
registrar_cifra("f01_genes_tandem_T01", n_t01, "T01.bloque_tandem (criterio mas laxo)")
if (length(en_bloque) != n_t01)
  message(sprintf("  AVISO: los intervalos de T15 contienen %d genes NAC; T01.bloque_tandem marca %d (criterio distinto)",
                  length(en_bloque), n_t01))

# Ideogramas verticales: los rotulos de bloque necesitan crecer a lo ancho y en
# la disposicion horizontal solo cabian sobre la propia barra. En vertical, el
# hueco entre cromosomas es espacio libre y es donde va el texto.
ETQ_QPCR   <- "gene measured by qRT-PCR"
ETQ_TANDEM <- "tandem block (T15)"
LARGO_MAX <- max(largos$largo_mb)
Y_NOMBRE  <- -0.115 * LARGO_MAX   # nombre del cromosoma, sobre el ideograma
Y_ARRIBA  <- -0.055 * LARGO_MAX   # y bajo el, cuantos NAC lleva

# El nombre del cromosoma va dentro del panel y no en el eje: con dos filas, un
# eje x compartido rotularia las dos con Chr1-Chr6.
nombres <- largos |> mutate(etq = sub("^chr", "Chr", cromosoma))

# Una sola capa de puntos para las dos marcas: asi la leyenda sale con la forma
# y el color correctos en cada clave. Con una capa por marca, ggplot dibuja el
# glifo de todas las capas en todas las claves.
marcas <- bind_rows(
  nac |> filter(valida) |> transmute(x = x - 0.17, y = medio_mb, fila, tipo = ETQ_QPCR),
  bl  |>                   transmute(x = x + 0.17, y = medio_mb, fila, tipo = ETQ_TANDEM))

p <- ggplot() +
  geom_segment(data = largos, aes(x = x, xend = x, y = 0, yend = largo_mb),
               colour = COL_IDEOGRAMA, linewidth = 3.4, lineend = "round") +
  geom_segment(data = nac, aes(x = x - 0.105, xend = x + 0.105,
                               y = medio_mb, yend = medio_mb),
               colour = COL_NAC, linewidth = 0.30) +
  # el gen rotulado lleva su marca en el mismo color que el rotulo
  geom_segment(data = filter(nac, valida),
               aes(x = x - 0.125, xend = x + 0.125, y = medio_mb, yend = medio_mb),
               colour = COL_DESTACADO, linewidth = 0.60) +
  # Los bloques se marcan con un triangulo al costado del ideograma y no con un
  # rectangulo sobre el: a escala de 70 Mb un bloque de decenas de kb colapsa a
  # una linea, de modo que el rectangulo no mostraba extension alguna y no
  # coincidia con su propia clave de leyenda. La extension real de cada bloque
  # esta en T15 y se dibuja a escala en F09.
  geom_point(data = marcas, aes(x = x, y = y, colour = tipo, shape = tipo), size = 1.2) +
  geom_text_repel(data = bl, aes(x = x + 0.24, y = medio_mb, label = etiqueta),
                  size = 1.85, lineheight = 0.86, colour = COL_BLOQUE, fontface = "bold",
                  hjust = 0, nudge_x = 0.13, direction = "both",
                  segment.size = 0.18, segment.colour = COL_BLOQUE,
                  min.segment.length = 0, box.padding = 0.10, point.padding = 0,
                  max.overlaps = Inf, seed = 1) +
  geom_text_repel(data = filter(nac, valida),
                  aes(x = x - 0.24, y = medio_mb, label = simbolo),
                  size = 1.95, colour = COL_DESTACADO, fontface = "bold",
                  hjust = 1, nudge_x = -0.13, direction = "y",
                  segment.size = 0.18, segment.colour = COL_DESTACADO,
                  min.segment.length = 0, box.padding = 0.10, point.padding = 0,
                  max.overlaps = Inf, seed = 1) +
  geom_text(data = nombres, aes(x = x, y = Y_NOMBRE, label = etq),
            size = 2.85, fontface = "bold", colour = "grey15") +
  geom_text(data = conteo, aes(x = x, y = Y_ARRIBA, label = paste0("n = ", n)),
            size = 2.05, colour = "grey30") +
  scale_colour_manual(values = setNames(c(COL_DESTACADO, COL_BLOQUE),
                                        c(ETQ_QPCR, ETQ_TANDEM)), name = NULL) +
  scale_shape_manual(values = setNames(c(16, 17), c(ETQ_QPCR, ETQ_TANDEM)), name = NULL) +
  # Margenes asimetricos: ggrepel confina el rotulo al panel, de modo que sin
  # hueco a la izquierda el nombre del gen se empuja sobre el primer ideograma
  # de la fila. El margen derecho aloja el rotulo de bloque del ultimo.
  scale_x_continuous(limits = c(0.02, 7.08), expand = c(0, 0)) +
  scale_y_reverse(breaks = seq(0, 200, 20)) +
  expand_limits(y = 1.16 * Y_NOMBRE) +
  facet_wrap(~ fila, ncol = 1) +
  labs(x = NULL, y = "position on chromosome (Mb)") +
  tema_genomico(9) + leyenda_interna() +
  theme(panel.grid.major.y = element_line(colour = "grey93", linewidth = 0.25),
        panel.spacing.y = unit(7, "pt"),
        strip.text = element_blank(), strip.background = element_blank(),
        axis.ticks.x = element_blank(), axis.text.x = element_blank())

# --- clave de bloques al pie -------------------------------------------------
# Tres columnas de siete entradas, ordenadas por codigo y recorridas por
# columna, para que buscar un B## sea bajar en vertical. Va como panel propio y
# no como caption porque asi las columnas se alinean por coordenada y no
# dependen de la anchura de la fuente.
CLAVE_COLS <- 3L
clave <- bl |>
  arrange(bloque) |>
  mutate(i = row_number() - 1L, n_filas = ceiling(n() / CLAVE_COLS),
         cx = i %/% n_filas, cy = i %% n_filas,
         txt = paste0(bloque, "  ", ids))

p_clave <- ggplot(clave, aes(x = cx, y = cy, label = txt)) +
  geom_text(hjust = 0, size = 1.85, colour = COL_BLOQUE) +
  scale_x_continuous(limits = c(-0.02, CLAVE_COLS), expand = c(0, 0)) +
  scale_y_reverse(expand = expansion(add = 0.7)) +
  labs(title = "EglNAC members of each tandem block (T15)") +
  theme_void() +
  theme(plot.title = element_text(size = 6.2, face = "bold", colour = "grey30",
                                  margin = margin(b = 2)),
        plot.margin = margin(2, 2, 0, 2))

# Alto fijo para la clave; el mapa se queda con el resto del lienzo.
ALTO_CLAVE <- 0.95
fig <- p / p_clave +
  plot_layout(heights = grid::unit(c(1, ALTO_CLAVE), c("null", "in")))

# El objeto que se guarda es el mapa sin clave: es lo componible. Si alguna vez
# se paneliza F01 hay que llevarse la clave, o los codigos B## quedan huerfanos.
guardar_objeto(p, "P_cromosomas")
guardar_figura(fig, "F01_cromosomas_NAC", alto = 7.10 + ALTO_CLAVE)

volcar_cifras("f01")
