# 15_arbol_SNAC.R -> F11: de quien es ortologo el gen focal
#
# Insumos de 14_ortologia_SNAC.py (entorno conda `filo`); este script solo dibuja.
# Dos paneles porque la pregunta tiene dos respuestas que no coinciden:
#   a  el arbol de genes del dominio NAM, que es la evidencia buena
#   b  las identidades por pares, que es la evidencia que usa el manuscrito y que
#      aqui se muestra para que se vea por que no decide: las dos regiones de la
#      misma proteina ordenan el trio al reves.
if (!exists("ANCHO_COLUMNA")) source("scripts/00_estilo.R")
suppressPackageStartupMessages({
  library(ggtree); library(treeio); library(ape); library(patchwork); library(ggrepel)
})

SNAC <- file.path(DATOS, "snac")
etq  <- suppressMessages(read_csv(file.path(SNAC, "etiquetas.csv"), show_col_types = FALSE))
t30  <- suppressMessages(read_csv(file.path(TABLAS, sprintf("T30_identidad_%s_vs_anclas.csv", GEN_FOCAL)),
                                  show_col_types = FALSE))
t32  <- suppressMessages(read_csv(file.path(TABLAS, "T32_arboles_SNAC.csv"), show_col_types = FALSE))
cif  <- suppressMessages(read_csv(file.path(TABLAS, "_cifras", "ortologia.csv"), show_col_types = FALSE))
val  <- function(k) cif$valor[match(k, cif$clave)]

# --- panel a: arbol del dominio -------------------------------------------
arb <- read.iqtree(file.path(SNAC, "arbol_dominio.treefile"))
phy <- as.phylo(arb)
stopifnot(is.rooted(phy) || all(c("ORE1", "NAP") %in% phy$tip.label))

# Grupos de color. El trio de Arabidopsis y el gen focal son el resultado; el
# resto es contexto y va en gris, siguiendo la regla de un solo destacado.
# Las claves son el valor guardado en T30.grupo; los valores, el rotulo visible.
GRUPOS <- c("trio de estres" = "ANAC019 / ANAC055 / RD26",
            "hermano ATAF"   = "ATAF clade",
            "SNAC-A de arroz" = "ATAF clade",
            "grupo externo"  = "outgroup")
etq <- etq |>
  mutate(clase = ifelse(focal, sprintf("%s (focal gene)", GEN_FOCAL),
                 ifelse(grupo %in% names(GRUPOS), GRUPOS[grupo], grupo)),
         clase = ifelse(startsWith(grupo, "subfamilia") & !focal,
                        sprintf("E. globulus subfamily %s", GEN_FOCAL_SUBFAM), clase))
# Vb lleva el mismo color que en el árbol de la familia (Fig. 2); el trío, el azul de los dos paneles
PAL <- c("#B2182B", "#2166AC", unname(SUBFAM_COLORES["Vb"]), "#9E9E9E") |>
  setNames(c(sprintf("%s (focal gene)", GEN_FOCAL), "ANAC019 / ANAC055 / RD26",
             sprintf("E. globulus subfamily %s", GEN_FOCAL_SUBFAM), "ATAF clade"))
PAL <- c(PAL, "outgroup" = "#C9C9C9")
etq$clase <- factor(etq$clase, levels = names(PAL))
stopifnot(!any(is.na(etq$clase)))

# El clado que contiene al gen focal y al trio: es la unidad de ortologia.
nodo_estres <- MRCA(phy, c("ANAC019", "RD26", GEN_FOCAL))
sop <- function(x) ifelse(is.na(x$SH_aLRT) | is.na(x$UFboot), NA_character_,
                          sprintf("%.0f/%.0f", x$SH_aLRT, x$UFboot))

base <- ggtree(arb, linewidth = 0.45, colour = "grey25") %<+% etq
# Los soportes se reparten con repel y no con un desplazamiento fijo: hay nodos
# a menos de una decima de rama y con hjust constante se pisan las cifras.
nodos <- base$data |> filter(!isTip, !is.na(SH_aLRT), !is.na(UFboot))
stopifnot(nrow(nodos) == Nnode(phy) - 1)   # la raiz no lleva soporte

p_arb <- base +
  geom_hilight(node = nodo_estres, fill = "#B2182B", alpha = 0.07, extend = 0.42) +
  geom_tiplab(aes(colour = clase, fontface = ifelse(focal %in% TRUE, "bold", "plain")),
              size = 2.25, offset = 0.012, show.legend = FALSE) +
  geom_text_repel(data = nodos, aes(x = x, y = y, label = sprintf("%.0f/%.0f", SH_aLRT, UFboot)),
                  size = 2.0, colour = "grey35", inherit.aes = FALSE,
                  nudge_x = -0.012, nudge_y = 0.34, direction = "both",
                  box.padding = 0.06, point.padding = 0, min.segment.length = 0.5,
                  segment.size = 0.15, segment.colour = "grey75", max.overlaps = Inf,
                  seed = 1) +
  geom_tippoint(aes(colour = clase), size = 1.25, stroke = 0) +
  geom_treescale(width = 0.1, fontsize = 2.1, linesize = 0.3, offset = 0.5, x = 0, y = 17.6) +
  scale_colour_manual(values = PAL, name = NULL, drop = FALSE) +
  guides(colour = guide_legend(override.aes = list(size = 1.7), nrow = 2, byrow = TRUE)) +
  theme_nac(8) +
  theme(panel.border = element_blank(), panel.grid = element_blank(), panel.grid.major.y = element_blank(),
        axis.text = element_blank(), axis.ticks = element_blank(),
        axis.title = element_blank(),
        legend.position = "bottom", legend.margin = margin(0, 0, 0, 0),
        legend.key.height = unit(0.28, "cm")) +
  coord_cartesian(clip = "off") +
  xlim(0, max(node.depth.edgelength(phy)) * 1.30)

# --- panel b: las dos identidades ------------------------------------------
# Un mismo par de proteinas da dos ordenes segun que region se compare. Las
# lineas que se cruzan entre las dos columnas son el argumento.
largo <- t30 |>
  select(ancla, grupo, dominio = ident_dominio_pct, completa = ident_completa_pct) |>
  pivot_longer(c(dominio, completa), names_to = "region", values_to = "pct") |>
  mutate(region = factor(region, levels = c("dominio", "completa"),
                         labels = c("NAM domain\n(125 aa)", "full protein")),
         trio = grupo == "trio de estres")

p_id <- ggplot(largo, aes(region, pct, group = ancla)) +
  geom_line(aes(colour = trio, linewidth = trio)) +
  geom_point(aes(colour = trio), size = 1.15, stroke = 0) +
  geom_text_repel(data = filter(largo, region == levels(largo$region)[1]),
                  aes(label = sprintf("%s  %.1f", ancla, pct), colour = trio),
                  size = 2.1, hjust = 1, direction = "y", xlim = c(NA, 0.92),
                  segment.size = 0.15, segment.colour = "grey70",
                  box.padding = 0.22, min.segment.length = 0.1, seed = 1, show.legend = FALSE) +
  geom_text_repel(data = filter(largo, region == levels(largo$region)[2]),
                  aes(label = sprintf("%.1f  %s", pct, ancla), colour = trio),
                  size = 2.1, hjust = 0, direction = "y", xlim = c(2.08, NA),
                  segment.size = 0.15, segment.colour = "grey70",
                  box.padding = 0.22, min.segment.length = 0.1, seed = 1, show.legend = FALSE) +
  scale_colour_manual(values = c(`TRUE` = "#2166AC", `FALSE` = "#9E9E9E"), guide = "none") +
  scale_linewidth_manual(values = c(`TRUE` = 0.55, `FALSE` = 0.3), guide = "none") +
  scale_x_discrete(expand = expansion(mult = 0.80)) +
  scale_y_continuous(labels = function(x) paste0(x, " %")) +
  # el título del eje y quedaba cortado por el panel a: va como título del panel
  labs(x = NULL, y = NULL, title = sprintf("Identity with %s", GEN_FOCAL)) +
  theme_nac(8) +
  theme(plot.title = element_text(size = 8, hjust = 0.5, margin = margin(0, 0, 3, 0)),
        panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.25))

p <- (p_arb | p_id) + plot_layout(widths = c(1.45, 1), guides = "collect") +
  plot_annotation(tag_levels = "a",
                  theme = theme(legend.position = "bottom",
                                legend.box.margin = margin(-6, 0, 0, 0)))
guardar_figura(p, "F11_ortologia_SNAC", alto = 4.25)

# --- cifras para la leyenda -------------------------------------------------
trio <- filter(t30, grupo == "trio de estres")
registrar_cifra("f11_clado_estres", val("ort_clado_estres"),
                "arbol_dominio.treefile", "clado minimo que contiene al gen focal y al trio")
registrar_cifra("f11_clado_soporte", val("ort_clado_estres_soporte"),
                "T32", "SH-aLRT/UFboot del clado de estres, arbol del dominio")
registrar_cifra("f11_clado_soporte_completa",
                t32$clado_estres_soporte[t32$etiqueta == "completa"],
                "T32", "el mismo clado en el arbol de proteina completa")
registrar_cifra("f11_hermano_soporte", val("ort_hermano_soporte"),
                "T32", "soporte del nodo que separa a RD26; por debajo de 80/95 = no resuelto")
registrar_cifra("f11_au_p", val("ort_au_p_trio"), "au_test.iqtree",
                "prueba AU: la topologia con el trio monofiletico no se rechaza")
registrar_cifra("f11_modelo", val("ort_modelo"), "IQ-TREE ModelFinder", "arbol del dominio")
registrar_cifra("f11_dominio_mejor", trio$ancla[which.max(trio$ident_dominio_pct)],
                "T30", "ancla mas parecida en el dominio NAM")
registrar_cifra("f11_completa_mejor", trio$ancla[which.max(trio$ident_completa_pct)],
                "T30", "ancla mas parecida en la proteina completa: no es la misma")
registrar_cifra("f11_rango_trio_dominio",
                sprintf("%.1f-%.1f %%", min(trio$ident_dominio_pct), max(trio$ident_dominio_pct)),
                "T30", "las tres identidades del trio en el dominio caben en este rango")
registrar_cifra("f11_lectores_identicos", val("ort_lectores_identicos_trio"),
                "T31", "posiciones que leen bases identicas en los tres de Arabidopsis")
volcar_cifras("f11")
