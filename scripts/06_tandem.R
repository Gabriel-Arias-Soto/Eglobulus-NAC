# 06_tandem.R -> panel de vecindades de los bloques en tandem (gggenes)
# Tres correcciones de la revision:
#  - el eje x era rango genico: B01 (266 kb) y B18 (22 kb) se dibujaban con el
#    mismo paso entre flechas. Ahora son coordenadas fisicas en kb, con escala
#    libre por bloque, de modo que el espaciado es el real.
#  - se dibujan los 21 bloques, no 9: la figura anterior mostraba 79 de 119
#    genes y saltaba de B08 a B18 sin advertirlo.
#  - la pertenencia a bloque se resuelve por SOLAPE DE COORDENADA con los
#    intervalos de T15 y no por la etiqueta T01.bloque_tandem: las dos tablas
#    implementan criterios distintos de tandem (T01 marca 125 genes en bloques
#    T01-T21; T15 da 119 en B01-B21, subconjunto estricto de los anteriores) y
#    el eje fisico de esta figura necesita los intervalos de T15. Los avisos de
#    abajo dejan constancia de la diferencia cada vez que se regenera.
if (!exists("ANCHO_COLUMNA")) source("scripts/00_estilo.R")
suppressPackageStartupMessages({library(gggenes); library(GenomicRanges)})

censo <- leer_censo()
gm <- mapa_cromosomas()
t14 <- suppressMessages(read_csv(file.path(TABLAS, "T14_pares_NAC_adyacentes.csv"), show_col_types = FALSE))
t15 <- suppressMessages(read_csv(file.path(TABLAS, "T15_bloques_tandem.csv"), show_col_types = FALSE))
t14$misma_subfamilia <- as.logical(t14$misma_subfamilia)   # el CSV trae "True"/"False"
stopifnot(!any(is.na(t14$misma_subfamilia)))

inv <- setNames(names(gm$mapa), gm$mapa)          # chrN -> accesion del ensamblado
gr_bl <- GRanges(unname(inv[t15$cromosoma]), IRanges(t15$inicio, t15$fin),
                 bloque = t15$bloque)

# --- genes NAC de cada bloque, por solape con el intervalo del bloque -------
nac_todos <- censo |>
  transmute(cromosoma, locus, simbolo, inicio, fin, hebra,
            clase = as.character(subfamilia_final), es_focal)
gr_nac <- GRanges(unname(inv[nac_todos$cromosoma]), IRanges(nac_todos$inicio, nac_todos$fin))
ovn <- findOverlaps(gr_bl, gr_nac)
nac_bl <- nac_todos[subjectHits(ovn), ] |>
  mutate(bloque = gr_bl$bloque[queryHits(ovn)]) |>
  select(bloque, cromosoma, simbolo, inicio, fin, hebra, clase, es_focal)

n_t01 <- sum(!is.na(censo$bloque_tandem))
if (nrow(nac_bl) != n_t01)
  message(sprintf("  AVISO: los intervalos de T15 contienen %d genes NAC; T01.bloque_tandem marca %d (criterio distinto)",
                  nrow(nac_bl), n_t01))
if (nrow(nac_bl) != sum(t15$n_genes))
  message(sprintf("  AVISO: el solape da %d genes y T15.n_genes suma %d", nrow(nac_bl), sum(t15$n_genes)))

# --- genes no-NAC intercalados, por solape real de intervalos ---------------
g <- gm$genes |> filter(scaf %in% unname(inv))
gr_g <- GRanges(g$scaf, IRanges(g$ini, g$fin), locus = g$locus, hebra = g$hebra)
ovo <- findOverlaps(gr_bl, gr_g)
otros <- tibble(bloque   = gr_bl$bloque[queryHits(ovo)],
                cromosoma = t15$cromosoma[match(gr_bl$bloque[queryHits(ovo)], t15$bloque)],
                simbolo  = gr_g$locus[subjectHits(ovo)],
                inicio   = start(gr_g)[subjectHits(ovo)],
                fin      = end(gr_g)[subjectHits(ovo)],
                hebra    = as.character(gr_g$hebra[subjectHits(ovo)])) |>
  filter(!simbolo %in% censo$locus) |>
  mutate(clase = ETQ_NO_NAC, es_focal = FALSE)

# --- coordenadas relativas al inicio de cada bloque y rotulo del panel ------
dat <- bind_rows(nac_bl, otros) |>
  left_join(t15 |> select(bloque, ini_bloque = inicio, n_genes, extension_kb), by = "bloque") |>
  mutate(x1 = (inicio - ini_bloque) / 1000, x2 = (fin - ini_bloque) / 1000,
         clase = factor(clase, levels = c(SUBFAM_NIVELES, SUBFAM_NO_EVAL,
                                          ETQ_NO_NAC)),
         # mismo separador y mismo orden -codigo primero- que el rotulo de F01,
         # que es por donde el lector salta de una figura a la otra
         etiqueta = paste(bloque, sub("^chr", "Chr", cromosoma),
                          sprintf("%d genes", n_genes),
                          sprintf("%.0f kb", extension_kb), sep = SEP_BLOQUE))
stopifnot(!any(is.na(dat$x1)), !any(is.na(dat$clase)))
dat$etiqueta <- factor(dat$etiqueta, levels = unique(dat$etiqueta[order(dat$bloque)]))
# el gris de "no evaluada (dominio parcial)" y el de los genes intercalados eran
# el mismo tono: en B02, que son 12 NAC de dominio parcial, el panel parecia
# vacio de NAC. Los intercalados pasan a blanco con borde.
PAL_BL <- c(PAL_SUBFAM, setNames("white", ETQ_NO_NAC))

registrar_cifra("f09_bloques", nrow(t15), "T15_bloques_tandem.csv")
registrar_cifra("f09_genes_nac", nrow(nac_bl), "solape de coordenada con T15")
registrar_cifra("f09_genes_nac_T01", n_t01, "T01.bloque_tandem (criterio mas laxo)")
registrar_cifra("f09_no_nac", nrow(otros), "GenomicRanges: solape con el intervalo del bloque")
registrar_cifra("f09_nac_dominio_parcial",
                sum(nac_bl$clase == SUBFAM_NO_EVAL),
                "genes de bloque sin subfamilia asignada (dominio parcial)")
bl_parcial <- nac_bl |> group_by(bloque) |>
  summarise(todos_parciales = all(clase == SUBFAM_NO_EVAL), n = n(), .groups = "drop") |>
  filter(todos_parciales, n >= 2)
registrar_cifra("f09_bloques_solo_parciales",
                if (nrow(bl_parcial)) paste(sprintf("%s (%d genes)", bl_parcial$bloque,
                                                    bl_parcial$n), collapse = "; ") else "ninguno",
                "bloques cuyos NAC son todos de dominio parcial")
registrar_cifra("f09_pares", nrow(t14), "T14_pares_NAC_adyacentes.csv")
registrar_cifra("f09_misma_hebra", sum(grepl("misma hebra", t14$orientacion)), "T14.orientacion")
registrar_cifra("f09_misma_hebra_pct",
                sprintf("%.0f %%", 100 * mean(grepl("misma hebra", t14$orientacion))), "T14")
registrar_cifra("f09_misma_subfam", sum(t14$misma_subfamilia), "T14.misma_subfamilia")
registrar_cifra("f09_misma_subfam_pct",
                sprintf("%.0f %%", 100 * mean(t14$misma_subfamilia)), "T14")
registrar_cifra("f09_extension_rango",
                sprintf("%.0f a %.0f kb", min(t15$extension_kb), max(t15$extension_kb)), "T15")

# B02 abarca 419 kb: un gen de 2 kb ocupa 2 px y el borde gris tapaba el relleno,
# de modo que sus 12 NAC parecian genes sin clasificar. El borde se mapea a la
# misma clase, asi que solo se ve en los intercalados (blancos).
PAL_BORDE <- c(PAL_SUBFAM, setNames(COL_NAC, ETQ_NO_NAC))
p <- ggplot(dat, aes(xmin = x1, xmax = x2, y = 1, fill = clase, colour = clase,
                     forward = hebra == "+")) +
  geom_gene_arrow(linewidth = 0.2,
                  arrowhead_height = unit(2.4, "mm"), arrowhead_width = unit(1.1, "mm"),
                  arrow_body_height = unit(2.0, "mm")) +
  geom_text(data = filter(dat, es_focal), aes(x = (x1 + x2) / 2, y = 1.42, label = simbolo),
            inherit.aes = FALSE, size = 1.9, fontface = "bold", colour = COL_DESTACADO) +
  facet_wrap(~ etiqueta, scales = "free_x", ncol = 2, strip.position = "top") +
  scale_fill_manual(values = PAL_BL, name = NULL,
                    breaks = c(SUBFAM_NIVELES, SUBFAM_NO_EVAL, ETQ_NO_NAC)) +
  scale_colour_manual(values = PAL_BORDE, guide = "none") +
  scale_y_continuous(limits = c(0.5, 1.7), expand = c(0, 0)) +
  labs(x = "position within block (kb from block start)", y = NULL) +
  guides(fill = guide_legend(ncol = 3, keyheight = unit(2.6, "mm"))) +
  tema_genomico(7) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        strip.text = element_text(size = 5.4, face = "bold", colour = "grey15"),
        legend.position = "bottom", legend.box.margin = margin(-3, 0, 0, 0))

guardar_objeto(p, "P_tandem")
guardar_figura(p, "F09_vecindades_tandem", alto = 6.90)

volcar_cifras("f09")
