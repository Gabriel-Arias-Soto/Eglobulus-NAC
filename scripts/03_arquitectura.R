# 03_arquitectura.R -> panel de arquitectura de dominios (gggenes)
# Correcciones que exigio la revision:
#  - el criterio de cobertura del perfil es explicito y la cifra sale de aqui:
#    la figura anterior afirmaba 163 loci con "cobertura completa" y ese numero
#    no se reproduce con ningun umbral (cobertura 1,0 -> 94; 0,99 -> 153).
#  - se separan las dos definiciones de helice transmembrana que estaban
#    confundidas: predicha (inicio_TM no nulo) y NTL (a <= UMBRAL_NTL_AA del
#    extremo C). La leyenda declara ambas.
#  - el orden de filas es subfamilia (SUBFAM_NIVELES) y luego largo, el mismo
#    que 04_motivos.R, para que los dos paneles se lean fila a fila.
if (!exists("ANCHO_COLUMNA")) source("scripts/00_estilo.R")
suppressPackageStartupMessages({library(gggenes); library(forcats)})

UMBRAL_COBERTURA <- 0.99   # fraccion del perfil PF02365 cubierta por el alineamiento
UMBRAL_NTL_AA    <- 50     # distancia maxima de la helice TM al extremo C

t16 <- suppressMessages(read_csv(file.path(TABLAS, "T16_arquitectura_dominios.csv"),
                                 show_col_types = FALSE))
stopifnot(nrow(t16) == 189)

arq <- t16 |>
  mutate(subfamilia_final = subfam_facetas(subfamilia_final),
         tm_predicha  = !is.na(inicio_TM),
         ntl_criterio = tm_predicha & (largo_aa - fin_TM) <= UMBRAL_NTL_AA,
         es_focal     = simbolo == GEN_FOCAL) |>
  arrange(subfamilia_final, largo_aa) |>
  mutate(fila = factor(simbolo, levels = simbolo))
stopifnot(sum(arq$es_focal) == 1)

cc <- filter(arq, cobertura_hmm >= UMBRAL_COBERTURA)
registrar_cifra("f07_umbral_cobertura", sprintf("%.0f %%", 100 * UMBRAL_COBERTURA),
                "03_arquitectura.R", "criterio fijado en el script, no heredado")
registrar_cifra("f07_loci_cobertura", nrow(cc),
                sprintf("T16: cobertura_hmm >= %.2f", UMBRAL_COBERTURA))
registrar_cifra("f07_nam_mediana", round(median(cc$largo_NAM)), "T16.largo_NAM")
registrar_cifra("f07_nam_dt",      round(sd(cc$largo_NAM)),     "T16.largo_NAM")
registrar_cifra("f07_nam_rango",
                sprintf("%.0f-%.0f", min(cc$largo_NAM), max(cc$largo_NAM)), "T16.largo_NAM")
registrar_cifra("f07_tr_rango", sprintf("%.0f a %.0f", min(arq$TR), max(arq$TR)), "T16.TR")
registrar_cifra("f07_tr_dt", round(sd(arq$TR)), "T16.TR")
registrar_cifra("f07_tm_predichas", sum(arq$tm_predicha), "T16: inicio_TM no nulo")
registrar_cifra("f07_ntl", sum(arq$ntl_criterio),
                sprintf("T16: helice TM a <= %d aa del extremo C", UMBRAL_NTL_AA))
registrar_cifra("f07_ntl_columna", sum(as.logical(t16$NTL)), "T16.NTL (columna publicada)")
if (sum(arq$ntl_criterio) != sum(as.logical(t16$NTL)))
  message(sprintf("  AVISO F07: el criterio de %d aa da %d NTL y la columna NTL da %d",
                  UMBRAL_NTL_AA, sum(arq$ntl_criterio), sum(as.logical(t16$NTL))))

ALTO_BARRA <- alto_barra(nrow(arq))
registrar_cifra("f07_alto_barra",
                sprintf("%.2f mm sobre un paso de %.2f mm",
                        as.numeric(convertHeight(ALTO_BARRA, "mm")),
                        as.numeric(convertHeight(ALTO_BARRA, "mm")) / 0.62),
                "03_arquitectura.R: alto_barra()",
                "a 1,5 mm las barras se solapaban y el panel se veia continuo")
sin_cabeza <- list(arrowhead_width = unit(0, "mm"), arrowhead_height = unit(0, "mm"),
                   arrow_body_height = ALTO_BARRA)
tm <- filter(arq, tm_predicha) |>
  mutate(clase = ifelse(ntl_criterio, DOM_NTL, DOM_TM))

p <- ggplot(arq, aes(xmin = 1, xmax = largo_aa, y = fila)) +
  do.call(geom_gene_arrow, c(list(aes(fill = DOM_TR), colour = NA), sin_cabeza)) +
  do.call(geom_subgene_arrow, c(list(aes(xsubmin = ali_ini, xsubmax = ali_fin,
                                         fill = DOM_NAM), colour = NA), sin_cabeza)) +
  do.call(geom_subgene_arrow, c(list(data = tm,
                                     aes(xsubmin = inicio_TM, xsubmax = fin_TM, fill = clase),
                                     colour = NA), sin_cabeza)) +
  # sin marca en el margen el rotulo rojo no se asocia a ninguna fila
  # la x del marcador viaja como columna del propio data de la capa: aes() no
  # puede mirar a `arq` ni a una variable suelta, porque el objeto guardado en
  # .rds se recompone en 12_paneles.R, donde ese entorno ya no existe.
  geom_point(data = filter(arq, es_focal) |>
               mutate(x_marca = -0.022 * max(arq$largo_aa)),
             aes(x = x_marca, y = fila),
             inherit.aes = FALSE, shape = 62, size = 1.9, colour = COL_DESTACADO) +
  geom_text(data = filter(arq, es_focal),
            aes(x = largo_aa + 12, y = fila, label = simbolo), inherit.aes = FALSE,
            # vjust: el gen focal es la primera fila de Vb; centrado, el borde del panel cortaba el rótulo
            hjust = 0, vjust = 0.9, size = 2.4, fontface = "bold", colour = COL_DESTACADO) +
  facet_grid(subfamilia_final ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_fill_manual(values = PAL_DOMINIO, name = NULL,
                    breaks = c(DOM_NAM, DOM_TR, DOM_NTL, DOM_TM)) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.06))) +
  labs(x = ETQ_EJE_AA, y = NULL) +
  guides(fill = guide_legend(ncol = 1, keyheight = unit(2.6, "mm"))) +
  theme_nac(8) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        # 0,6 pt no bastaba para distinguir el corte entre subfamilias del
        # aire que ahora separa filas dentro de cada bloque
        panel.grid.major.y = element_blank(), panel.spacing.y = unit(1.8, "pt"),
        strip.placement = "outside", strip.background = element_blank(),
        strip.text.y.left = element_text(angle = 0, size = 6.5, hjust = 1, face = "plain"),
        legend.text = element_text(size = 7), strip.clip = "off",
        legend.position = "bottom", legend.box.margin = margin(-3, 0, 0, 0))

guardar_objeto(p, "P_arquitectura")
# mismo alto que FIG02: la altura de barra se calcula al construir el objeto y
# se recompone tal cual en 12_paneles.R, asi que las dos salidas deben coincidir
guardar_figura(p, "F07_arquitectura_dominios", alto = ALTO_PANELES, ancho = ANCHO_MEDIA)

volcar_cifras("f07")
