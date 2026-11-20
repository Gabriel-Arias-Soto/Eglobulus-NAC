# 01_filogenia.R -> panel de la filogenia de maxima verosimilitud
# El color de punta sale SIEMPRE de T01.subfamilia_final. La revision encontro
# que F06 pintaba IVa 47 donde el censo da 53 porque su capa de color venia de
# una asignacion anterior: uniendo censo y arbol por locus eso no puede repetirse.
if (!exists("ANCHO_COLUMNA")) source("scripts/00_estilo.R")
suppressPackageStartupMessages({library(ggtree); library(treeio); library(ape)})

MIN_SOPORTE <- 95    # UFBoot a partir del cual se marca el nodo

censo <- leer_censo()
arb <- ape::read.tree(file.path(RAIZ, "arbol_NAC_189.treefile"))
# ladderize ordena los clados por tamano: sin esto el abanico alterna ramas
# largas y cortas y los clados de subfamilia salen troceados en el circulo.
arb <- ape::ladderize(arb)

# las puntas llevan prefijo de especie: GL_<locus> (globulus) y GR_<id> (anclas)
meta <- tibble(label = arb$tip.label) |>
  mutate(locus   = sub("^G[LR]_", "", label),
         especie = ifelse(grepl("^GL_", label), "globulus", "grandis")) |>
  left_join(censo |> select(locus, simbolo, subfamilia_final, es_focal), by = "locus") |>
  mutate(grupo = ifelse(especie == "grandis", ETQ_ANCLA,
                        as.character(subfamilia_final)),
         grupo = factor(grupo, levels = c(SUBFAM_NIVELES, SUBFAM_NO_EVAL,
                                          ETQ_ANCLA)),
         es_focal = !is.na(es_focal) & es_focal)

n_glob <- sum(meta$especie == "globulus")
# las puntas de globulus son los genes del censo con subfamilia; el control se
# deriva del censo, no se teclea
stopifnot(n_glob == sum(censo$subfamilia_final != SUBFAM_NO_EVAL), sum(meta$es_focal) == 1)
if (any(meta$especie == "globulus" & is.na(meta$subfamilia_final)))
  stop("hay puntas de globulus sin subfamilia en el censo")

# conteos que iran a la leyenda, tomados del censo y no del arbol
cn <- censo |> count(subfamilia_final, name = "n") |> filter(n > 0)
registrar_cifra("f06_puntas_total",   length(arb$tip.label), "arbol_NAC_189.treefile")
registrar_cifra("f06_puntas_globulus", n_glob, "arbol x T01 por locus")
registrar_cifra("f06_anclas", sum(meta$especie == "grandis"), "arbol: puntas GR_")
registrar_cifra("f06_subfamilias",
                sum(cn$subfamilia_final != SUBFAM_NO_EVAL), "T01.subfamilia_final")
registrar_cifra("f06_con_subfamilia",
                sum(censo$subfamilia_final != SUBFAM_NO_EVAL), "T01.subfamilia_final")
registrar_cifra("f06_dominio_parcial",
                sum(censo$subfamilia_final == SUBFAM_NO_EVAL),
                "T01: los 6 loci sin subfamilia son los 6 ausentes del arbol")
registrar_cifra("f06_conteos",
                paste(sprintf("%s %d", cn$subfamilia_final, cn$n), collapse = ", "), "T01")

p <- ggtree(arb, layout = "fan", open.angle = 14, linewidth = 0.28, colour = "grey85")
p <- p %<+% meta
# soportes: IQ-TREE escribe "SH-aLRT/UFBoot" en la etiqueta de nodo
p$data$soporte <- suppressWarnings(as.numeric(sub(".*/", "", p$data$label)))
p$data$soporte[p$data$isTip] <- NA_real_
# la condicion se materializa como COLUMNA y no como expresion dentro de aes():
# el objeto se guarda en .rds y se compone en otro proceso, donde MIN_SOPORTE ya
# no existe. Una aes() que mire constantes del script romperia la panelizacion.
p$data$sop_alto <- !is.na(p$data$soporte) & p$data$soporte >= MIN_SOPORTE

registrar_cifra("f06_nodos_soporte", sum(p$data$sop_alto),
                sprintf("etiquetas de nodo con UFBoot >= %d", MIN_SOPORTE))
registrar_cifra("f06_min_soporte", MIN_SOPORTE, "01_filogenia.R")

# etiquetas y paleta se materializan AQUI, no como closures dentro de la escala:
# el objeto viaja en .rds a otro proceso y una funcion que mire 'cn' romperia.
# --- color de rama por consenso de subfamilia --------------------------------
# Una rama se pinta con el color de una subfamilia solo si TODAS las hojas de
# globulus que descienden de ella pertenecen a esa subfamilia; si se mezclan, va
# en gris. Las 156 anclas de grandis llevan otra nomenclatura de locus y no
# estan en el censo, asi que no cuentan para el consenso: su rama terminal se
# pinta con el gris de ancla, coherente con la clave de la leyenda.
n_tip <- length(arb$tip.label); n_tot <- n_tip + arb$Nnode
ix <- match(arb$tip.label, meta$label)
stopifnot(!any(is.na(ix)))
subf_hoja <- ifelse(meta$especie[ix] == "globulus",
                    as.character(meta$subfamilia_final[ix]), NA_character_)
acum <- vector("list", n_tot)
acum[seq_len(n_tip)] <- as.list(subf_hoja)
# postorden: garantiza que los hijos se acumulan antes que el padre
ed_po <- ape::reorder.phylo(arb, "postorder")$edge
for (i in seq_len(nrow(ed_po))) {
  pa <- ed_po[i, 1]; ch <- ed_po[i, 2]
  acum[[pa]] <- unique(c(acum[[pa]], acum[[ch]]))
}
cons <- vapply(acum, function(v) {
  v <- v[!is.na(v)]
  if (length(v) == 1) v else NA_character_
}, character(1))
p$data$rama <- cons[p$data$node]
p$data$rama[p$data$isTip & p$data$especie == "grandis"] <- ETQ_ANCLA
p$data$rama <- factor(p$data$rama, levels = c(SUBFAM_NIVELES, SUBFAM_NO_EVAL,
                                              ETQ_ANCLA))
registrar_cifra("f06_ramas_puras",
                sprintf("%d de %d", sum(!is.na(cons)), length(cons)),
                "ramas cuyas hojas de globulus son de una sola subfamilia")
registrar_cifra("f06_ramas_puras_internas",
                sprintf("%d de %d", sum(!is.na(cons[(n_tip + 1):n_tot])), arb$Nnode),
                "nodos internos con descendencia de una sola subfamilia")

R_MAX  <- max(p$data$x, na.rm = TRUE)
registrar_cifra("f06_rama_max", sprintf("%.2f", max(arb$edge.length)),
                "arbol_NAC_189.treefile: rama mas larga")
registrar_cifra("f06_profundidad", sprintf("%.2f", R_MAX),
                "arbol_NAC_189.treefile: distancia raiz-punta maxima")

# La clave suma 185, no 189, y es correcto: los 4 genes de dominio NAM parcial
# no entran en el alineamiento y por tanto no tienen punta en el arbol (la
# cifra f06_dominio_parcial lo deja registrado para el pie de figura).
BRK <- c(SUBFAM_NIVELES, ETQ_ANCLA)
mcn <- setNames(cn$n, as.character(cn$subfamilia_final))
ETIQ <- ifelse(BRK %in% names(mcn), sprintf("%s (%d)", BRK, mcn[BRK]), BRK)
# el ancla es contexto y va al gris mas claro; el sin evaluar, al gris medio
PAL_ARBOL <- c(PAL_SUBFAM, setNames("grey78", ETQ_ANCLA))
BRK  <- BRK[BRK %in% unique(as.character(p$data$grupo))]
ETIQ <- ETIQ[c(SUBFAM_NIVELES, ETQ_ANCLA) %in% BRK]

p <- p +
  # las ramas se repintan por consenso de subfamilia sobre el arbol gris de base
  geom_tree(aes(colour = rama), linewidth = 0.28) +
  # a 0,28 pt los soportes no se distinguian del grosor de la rama
  geom_point2(aes(subset = sop_alto),
              size = 0.5, colour = "grey25", stroke = 0) +
  # Sin anillo de puntas alineadas: el color de rama da la misma lectura junto a
  # la rama, y un anillo a radio constante deja una corona vacia entre las puntas
  # reales y el borde.
  # La punta de cada gen conserva su color porque una rama terminal nunca mezcla
  # subfamilias; las anclas quedan en el gris de su clave.
  geom_tippoint(aes(colour = rama), size = 0.7, stroke = 0) +
  geom_point2(aes(subset = es_focal), shape = 8, size = 2.3,
              colour = COL_DESTACADO, stroke = 0.6) +
  geom_treescale(width = 0.5, fontsize = 2.2, linesize = 0.3, offset = 2) +
  scale_colour_manual(values = PAL_ARBOL, name = "subfamily",
                      breaks = BRK, labels = ETIQ, na.value = "grey85") +
  # sin relleno por dentro: la expansion inferior por defecto (5 %) abre un
  # agujero en el centro del abanico y roba ese mismo 5 % al radio util.
  scale_x_continuous(expand = expansion(mult = c(0, 0.01))) +
  # una sola columna: el abanico mantiene relacion de aspecto 1, asi que su
  # diametro es el minimo entre alto y ancho del panel. La leyenda a dos
  # columnas gastaba 1,6 in de ancho y ese lado era el que limitaba. A una
  # columna baja a ~0,95 in y el circulo crece sin tocar el alto de pagina.
  guides(colour = guide_legend(override.aes = list(size = 1.9), ncol = 1,
                               title.position = "top", keyheight = unit(7, "pt"))) +
  theme_nac(8) +
  # La leyenda va a la derecha: el abanico mantiene relacion de aspecto fija, de
  # modo que con la leyenda abajo el circulo solo ocupaba el 60 % del ancho y
  # sobraba papel a los lados. Aqui llena el alto y la leyenda usa ese hueco.
  #
  # Los hijos de panel.grid se anulan uno a uno: theme(panel.grid = element_blank())
  # NO borra panel.grid.major.x/y si el tema base los define explicitamente, y
  # en el abanico esas lineas salen como radios tenues desde el centro.
  theme(legend.position = "right", legend.justification = "centre",
        legend.margin = margin(0, 0, 0, -6),
        panel.border = element_blank(),
        panel.grid = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(), panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.text = element_blank(), axis.ticks = element_blank(),
        axis.title = element_blank())

# --- rotulo del gen protagonista, FUERA del anillo de puntas -----------------
# Antes el rotulo se anclaba en las coordenadas de la propia punta. Como esa
# punta cae muy adentro (radio 2,20 frente a R_MAX 4,29), la caja blanca se
# comia los clados vecinos. Ahora: linea guia radial desde la punta hasta el
# borde del abanico y texto ya fuera del anillo, al estilo align = TRUE.
foc <- p$data |> filter(!is.na(es_focal), es_focal)
ANG_FOCAL <- foc$angle[1]

# Orientacion del texto. Un rotulo radial nunca cruza otro clado, pero cerca de
# los polos (arriba / abajo del circulo) sale casi vertical y se lee mal; ahi
# conviene horizontal, que ademas cabe porque el panel es mas ancho que alto.
# El umbral de 35 grados es la franja donde lo horizontal todavia no invade
# clados vecinos. Con ANG_FOCAL = 265,9 el gen cae en esa franja.
CERCA_POLO <- min(abs(ANG_FOCAL - 90), abs(ANG_FOCAL - 270)) <= 35
INVERTIDO  <- ANG_FOCAL > 90 && ANG_FOCAL < 270   # a la izquierda saldria del reves

# Radio del rotulo. Anclarlo al anillo GLOBAL (R_MAX) lo dejaba a 2 unidades de
# su propia punta, con una guia larguisima cruzando papel vacio, porque el clado
# Vb es corto (sus puntas llegan a 2,41 y el anillo global a 4,29). Se ancla al
# vecindario angular: el punto mas externo de las puntas que el rotulo podria
# tocar, mas un margen. Un rotulo horizontal de ~9 caracteres abarca unos 23
# grados de arco a este radio, de ahi la ventana de +-15 grados.
VECINDAD_ROTULO <- 15
dtip  <- p$data[p$data$isTip, ]
cerca <- abs(((dtip$angle - ANG_FOCAL + 180) %% 360) - 180) <= VECINDAD_ROTULO
R_VECINO <- max(dtip$x[cerca])

# Un rotulo radial si sale por fuera de todo, porque crece a lo largo del radio
# y no barre arco; el horizontal se queda pegado a su clado.
R_ROTULO <- if (CERCA_POLO) R_VECINO + 0.12 * R_MAX else R_MAX * 1.03
# Lo que se reserve de radio encoge el arbol en la misma proporcion: coord_polar
# mete el circulo en el 80 % del lado menor del panel.
R_LIMITE <- if (CERCA_POLO) max(R_MAX * 1.02, R_ROTULO * 1.06) else R_MAX * 1.42

registrar_cifra("f06_radio_rotulo", sprintf("%.2f (vecindario %.2f, anillo %.2f)",
                                            R_ROTULO, R_VECINO, R_MAX),
                "01_filogenia.R: colocacion del rotulo focal")

# Todo viaja como columnas del data de la capa: el objeto se guarda en .rds y
# se recompone en 12_paneles.R, donde estas variables ya no existen.
rot <- tibble(x = foc$x[1], xend = R_ROTULO * 0.94, xr = R_ROTULO,
              y = foc$y[1], yend = foc$y[1],
              etiqueta = GEN_FOCAL,
              ang = if (CERCA_POLO) 0 else if (INVERTIDO) ANG_FOCAL + 180 else ANG_FOCAL,
              hj  = if (CERCA_POLO) 0.5 else if (INVERTIDO) 1 else 0)

p <- p +
  geom_segment(data = rot, aes(x = x, xend = xend, y = y, yend = yend),
               inherit.aes = FALSE, linewidth = 0.2, linetype = "12",
               colour = COL_DESTACADO, lineend = "round") +
  geom_text(data = rot, aes(x = xr, y = y, label = etiqueta,
                            angle = ang, hjust = hj),
            inherit.aes = FALSE, size = 2.5, fontface = "bold",
            colour = COL_DESTACADO) +
  expand_limits(x = R_LIMITE)

registrar_cifra("f06_angulo_focal", sprintf("%.0f grados", ANG_FOCAL),
                "posicion angular de EglNAC66 en el abanico")

guardar_objeto(p, "P_filogenia")
# coord_polar fuerza panel cuadrado: el diametro lo fija el lado MENOR del
# panel. Medido con ggplotGrob, la leyenda y los margenes gastan 1,05 in de los
# 6,30 de ancho, asi que el panel mide 5,25 in de ancho; con alto 5,30 el panel
# queda en 5,24 x 5,25 y el circulo alcanza su maximo. Por encima de ese alto
# solo se anade papel en blanco.
# la Fig. 2 del cuerpo la dibuja 01b_filogenia_cuerpo.R, que ejecuta este guion antes
if (!isTRUE(getOption("nac.f06_v2"))) guardar_figura(p, "F06_filogenia_NAC", alto = 5.30)

volcar_cifras("f06")
