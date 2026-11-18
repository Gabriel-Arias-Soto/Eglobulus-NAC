# 07_expresion.R -> mapa de calor de expresion (ComplexHeatmap)
# La revision encontro tres cosas en F02: el orden de las filas de subfamilia era
# ASCII (artefacto de ordenar cadenas), la figura medía 16,8 in de alto a ancho de
# columna (el alto util son 8,70) y las etiquetas de fila quedaban a ~4 pt.
# ComplexHeatmap resuelve las tres: la subfamilia pasa de ir en el texto de cada
# fila a una banda de anotacion, row_split agrupa en el orden de SUBFAM_NIVELES,
# el reparto cuerpo/suplementario recorta el alto.
#
# CAMBIO DE INSUMO. Antes leia DE_NAC_modelo_conjunto.csv, en loci Eucgr de
# E. grandis y con los simbolos del manuscrito como unica clave: el mismo texto
# "EglNAC33" nombra loci distintos en esa nomenclatura y en la del censo, y
# cruzar por el nombre emparejaba filas con el locus equivocado. Ahora lee
# T36_DE_NAC_globulus.csv, que 17_de_globulus.R produce cuantificando contra el
# genoma conspecifico: viene en loci ACJRO7_*, que es la clave primaria del
# censo, y no queda ningun emparejamiento por texto.
if (!exists("ANCHO_COLUMNA")) source("scripts/00_estilo.R")
suppressPackageStartupMessages({library(ComplexHeatmap); library(circlize); library(tidyr)})

UMBRAL_FDR   <- 0.05
UMBRAL_LOGFC <- 1
# Un locus tiene senal propia si al menos la mitad de sus lecturas las asigna
# salmon de forma unica, en vez de repartirlas por el EM entre transcritos
# indistinguibles (T35). Sustituye a la antigua banda "senal atribuible", que
# se definia por el numero de loci del censo que compartian co-ortologo en
# E. grandis: contra el conspecifico esa pregunta ya no se plantea.
UMBRAL_UNICA <- 0.5
# Contrastes de edgeR tal como vienen en la tabla, en el orden en que deben
# leerse. Si el manuscrito cambia de nomenclatura se edita AQUI y nada mas.
# el frio va en la secuencia del ensayo: aclimatacion antes de la helada (CABF)
# y despues (CAAF), no en orden alfabetico
# El mapa de calor va SOLO en hoja y sin la comparacion entre condiciones de
# frio: cada columna enfrenta una condicion al control (NA), igual que los tres
# paneles de F03. Eso deja fuera los cuatro contrastes de raiz y el CABF_CAAF
# (frio contra frio), que si esta en la tabla de DE.
TEJIDO <- "leaf"
CONTRASTES <- c("Leaf_NA_DS", "Leaf_NA_CABF", "Leaf_NA_CAAF")

censo <- leer_censo()
de <- suppressMessages(read_csv(file.path(TABLAS, "T36_DE_NAC_globulus.csv"),
                                show_col_types = FALSE)) |>
  filter(contraste %in% CONTRASTES)
stopifnot(nrow(de) > 0, all(de$locus %in% censo$locus))

ancha <- de |> select(locus, contraste, logFC) |>
  pivot_wider(names_from = contraste, values_from = logFC)
fdr <- de |> select(locus, contraste, FDR) |>
  pivot_wider(names_from = contraste, values_from = FDR)
falta <- setdiff(CONTRASTES, names(ancha))
if (length(falta)) stop("faltan contrastes en la tabla: ", paste(falta, collapse = ", "))

PAL_F02 <- PAL_SUBFAM

info <- tibble(locus = ancha$locus) |>
  left_join(de |> distinct(locus, simbolo, subfamilia_final, es_focal, frac_unica),
            by = "locus") |>
  mutate(subfamilia_final = factor(as.character(subfamilia_final),
                                   levels = c(SUBFAM_NIVELES, SUBFAM_NO_EVAL)),
         # dos propiedades distintas que la version anterior confundia en una
         # sola banda rotulada "qRT-PCR": si la senal del locus es suya (T35) y
         # si el gen se midio por qPCR (los cuatro de T02.medido_por_qPCR)
         unica = !is.na(frac_unica) & frac_unica >= UMBRAL_UNICA,
         qpcr  = simbolo %in% VALIDADOS)
stopifnot(!anyNA(info$simbolo), !anyNA(info$subfamilia_final))

info$etiqueta <- info$simbolo
# dentro de cada subfamilia, por numero y no por texto (antes: 1, 101, 12, 13, 2)
orden <- order(info$subfamilia_final, as.integer(sub("^\\D*(\\d+).*$", "\\1", info$etiqueta)))
info <- info[orden, ]
mat  <- as.matrix(ancha[orden, CONTRASTES]); rownames(mat) <- info$etiqueta
mfdr <- as.matrix(fdr[orden, CONTRASTES])

sig <- apply(mfdr < UMBRAL_FDR & abs(mat) >= UMBRAL_LOGFC, 1,
             function(z) any(z, na.rm = TRUE))
# "NA_CABF" no dice nada al lector: los codigos de Materiales y metodos se
# rotulan como la comparacion que representan (referencia -> condicion)
etiq_col <- sub("^(Leaf|Root)_", "", CONTRASTES)
etiq_col <- sub("_", " \u2192 ", etiq_col)

registrar_cifra("f02_genes", nrow(mat), "T36_DE_NAC_globulus.csv")
registrar_cifra("f02_genes_del_censo", sprintf("%d de %d", nrow(mat), nrow(censo)),
                "loci del censo que pasan el filtro de expresion de edgeR")
registrar_cifra("f02_contrastes", ncol(mat), "DE.contraste")
registrar_cifra("f02_contrastes_cuales", paste(CONTRASTES, collapse = ", "), "DE.contraste")
registrar_cifra("f02_tejidos", TEJIDO, "DE.contraste: solo hoja por decision de figura")
registrar_cifra("f02_contrastes_excluidos",
                "Root_NA_DS, Root_NA_CABF, Root_NA_CAAF, Root_CABF_CAAF, Leaf_CABF_CAAF",
                "estan en la tabla de DE y no se pintan")
registrar_cifra("f02_significativos", sum(sig),
                sprintf("FDR < %.2f y |log2FC| >= %g en al menos un contraste",
                        UMBRAL_FDR, UMBRAL_LOGFC))
registrar_cifra("f02_umbrales", sprintf("FDR < %.2f y |log2FC| >= %g", UMBRAL_FDR, UMBRAL_LOGFC),
                "07_expresion.R")
registrar_cifra("f02_celdas_sin_dato", sum(is.na(mat)), "DE: combinaciones sin estimacion")
registrar_cifra("f02_rango_logfc",
                sprintf("%.1f a %.1f", min(mat, na.rm = TRUE), max(mat, na.rm = TRUE)), "DE.logFC")
registrar_cifra("f02_medidos_qpcr", sum(info$qpcr), "T02.medido_por_qPCR")
# los cuatro genes con qRT-PCR tienen que estar en la figura; si alguno no pasa
# el filtro de expresion hay que declararlo, porque la banda solo puede marcar
# los que tienen estimacion
fuera <- setdiff(VALIDADOS, info$simbolo)
if (length(fuera)) {
  message("  AVISO F02: ", length(fuera), " de ", length(VALIDADOS),
          " genes con qRT-PCR no pasan el filtro de expresion: ",
          paste(fuera, collapse = ", "))
  registrar_cifra("f02_qpcr_sin_estimacion_DE", paste(sort(fuera), collapse = ", "),
                  "VALIDADOS ausentes de T36_DE_NAC_globulus.csv")
}
registrar_cifra("f02_senal_unica", sum(info$unica),
                sprintf("T35: mediana de lecturas asignadas de forma unica >= %.0f %%",
                        100 * UMBRAL_UNICA))
registrar_cifra("f02_regla_senal_unica",
                sprintf(paste("mediana sobre las 12 librerias de UniqueCount / (UniqueCount +",
                              "AmbigCount) agregada por locus (salmon aux_info/ambig_info.tsv);",
                              "umbral %.2f. Mediana del juego %.2f, rango %.2f a %.2f"),
                        UMBRAL_UNICA, median(info$frac_unica, na.rm = TRUE),
                        min(info$frac_unica, na.rm = TRUE), max(info$frac_unica, na.rm = TRUE)),
                "T35_atribuibilidad_globulus.csv")

lim <- ceiling(max(abs(mat), na.rm = TRUE))
col_fun <- colorRamp2(c(-lim, 0, lim), c("#2166AC", "white", "#B2182B"))
# OJO: sin 'at' explicito ComplexHeatmap deja las marcas a pretty(), que con
# lim = 6 estira la barra hasta +-10 y satura el degradado: el lector ve una
# escala que llega al doble del maximo real y casi todas las celdas palidas.
at_lim <- unique(round(c(-lim, -lim / 2, 0, lim / 2, lim)))
gp6 <- gpar(fontsize = 6); gp7b <- gpar(fontsize = 7, fontface = "bold")

pinta <- function(idx, tam_fila, ancho_celda, lado = "right") {
  # al costado las leyendas se apilan en vertical: la barra de log2FC deja de
  # ser una regla horizontal y la clave de subfamilia pasa a una sola columna,
  # que es lo que hace que la pila quepa en el margen sin ensancharla
  vert <- lado %in% c("right", "left")
  m <- mat[idx, , drop = FALSE]; iv <- info[idx, ]
  colnames(m) <- etiq_col
  # a ancho de columna el rotulo largo del grupo sin evaluar se sale del lienzo;
  # la clave completa sigue en la leyenda de subfamilia
  tit <- levels(droplevels(iv$subfamilia_final))
  tit[tit == SUBFAM_NO_EVAL] <- sub(" ", "\n", ETQ_PARCIAL)
  ha <- rowAnnotation(
    subfamily = droplevels(iv$subfamilia_final),
    `unique signal` = ifelse(iv$unica, "yes", "shared"),
    `qRT-PCR` = ifelse(iv$qpcr, "measured", "no"),
    col = list(subfamily = PAL_F02,
               `unique signal` = c("yes" = "grey35", "shared" = "grey85"),
               `qRT-PCR` = c("measured" = COL_DESTACADO, "no" = "grey93")),
    # los rotulos de las bandas van ARRIBA: al pie chocan con las etiquetas de
    # columna giradas 45 grados
    annotation_name_side = "top",
    annotation_name_gp = gpar(fontsize = 5.5), gap = unit(0.5, "mm"),
    simple_anno_size = unit(2.4, "mm"),
    annotation_legend_param = list(
      subfamily = list(ncol = if (vert) 1 else 2, labels_gp = gp6, title_gp = gp7b,
                       grid_height = unit(2.4, "mm"), grid_width = unit(2.4, "mm")),
      `unique signal` = list(labels_gp = gp6, title_gp = gp7b,
                       grid_height = unit(2.4, "mm"), grid_width = unit(2.4, "mm")),
      `qRT-PCR` = list(labels_gp = gp6, title_gp = gp7b,
                       grid_height = unit(2.4, "mm"), grid_width = unit(2.4, "mm"))))
  Heatmap(m, name = "log\u2082FC", col = col_fun, na_col = "grey90",
          cluster_rows = FALSE, cluster_columns = FALSE,
          row_split = droplevels(iv$subfamilia_final), row_gap = unit(0.5, "mm"),
          row_title = tit, row_title_rot = 0, row_title_gp = gpar(fontsize = 5.5),
          column_title = TEJIDO, column_title_gp = gp7b,
          width = unit(ncol(m) * ancho_celda, "mm"),
          show_row_names = TRUE, row_names_gp = gpar(fontsize = tam_fila),
          column_names_gp = gpar(fontsize = 6), column_names_rot = 45,
          left_annotation = ha, border = "grey30",
          heatmap_legend_param = c(
            list(at = at_lim, labels = as.character(at_lim),
                 labels_gp = gp6, title_gp = gp7b),
            if (vert) list(direction = "vertical", legend_height = unit(26, "mm"))
            else      list(direction = "horizontal", legend_width = unit(26, "mm"))))
}

dibuja <- function(idx, tam_fila, ancho_celda, hueco = 3, lado = "right") function() {
  # legend_gap lleva (vertical, horizontal). Al costado las leyendas se apilan
  # y manda el primero; al pie van en una sola fila y manda el segundo, que ahi
  # no lo fija el cuerpo de la leyenda sino el titulo: "unique signal" es mas
  # ancho que su cuadro y sus dos etiquetas cortas y se desborda sobre el
  # titulo siguiente, por eso el hueco al pie iba por figura.
  ht_opt$legend_gap <- if (lado %in% c("right", "left")) unit(c(3, 4), "mm")
                       else unit(c(4, hueco), "mm")
  draw(pinta(idx, tam_fila, ancho_celda, lado), heatmap_legend_side = lado,
       annotation_legend_side = lado, merge_legend = TRUE,
       padding = unit(c(1, 1, 1, 1), "mm"))
}

# el cuerpo lleva los significativos MAS los cuatro medidos por qPCR, aunque
# alguno no pase el umbral: la figura tiene que mostrar lo que se valido
idx_cuerpo <- sort(union(which(sig), which(info$qpcr)))
guardar_base(dibuja(idx_cuerpo, 6, 7, hueco = 6), ancho = ANCHO_MEDIA,
             "F02_heatmap_expresion",
             alto = min(ALTO_MAXIMO, 2.2 + 0.105 * length(idx_cuerpo)))
# La suplementaria existe para poder leer los 114 nombres de gen, y el tope de
# pagina la dejaba en 4,1 pt de alto de fila para una fuente de 5: los nombres
# se tocaban. Al no ir en el cuerpo no la ata ALTO_MAXIMO, asi que usa su alto
# natural con 0,072 in por fila, que deja 5,2 pt por nombre.
guardar_base(dibuja(seq_len(nrow(mat)), 5, 14, hueco = 10), ancho = ANCHO_COLUMNA,
             "F02S_heatmap_expresion_completo",
             alto = 2.2 + 0.072 * nrow(mat), suplementaria = TRUE)

pres <- levels(droplevels(info$subfamilia_final[idx_cuerpo]))
registrar_cifra("f02_filas_cuerpo", length(idx_cuerpo), "significativos union medidos por qPCR")
registrar_cifra("f02_qpcr_en_figura",
                paste(sort(info$etiqueta[idx_cuerpo][info$qpcr[idx_cuerpo]]), collapse = ", "),
                "T02.medido_por_qPCR")
registrar_cifra("f02_qpcr_no_significativos",
                paste(sort(info$etiqueta[info$qpcr & !sig]), collapse = ", ") |>
                  (\(x) if (nzchar(x)) x else "ninguno")(),
                "T02.medido_por_qPCR frente a los umbrales de la figura")
registrar_cifra("f02_dominio_parcial_n", sum(info$subfamilia_final == SUBFAM_NO_EVAL), "T01.subfamilia_final")
registrar_cifra("f02_subfam_en_figura", paste(pres, collapse = ", "), "T01.subfamilia_final")
registrar_cifra("f02_subfam_sin_representacion",
                paste(setdiff(SUBFAM_NIVELES, pres), collapse = ", "),
                "subfamilias del censo sin ninguna fila en la figura")
registrar_cifra("f02_nomenclatura_filas", "censo (T01.simbolo)",
                "misma nomenclatura que el resto del juego de figuras")
registrar_cifra("f02_referencia_expresion",
                "GCA_046226855.1 (E. globulus), PRJNA486291, hoja, 4 condiciones x 3 replicas",
                "17_de_globulus.R")
volcar_cifras("f02")
