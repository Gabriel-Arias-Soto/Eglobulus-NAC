# 04_motivos.R -> panel de subdominios A-E y motivos C-terminales (gggenes)
# La afirmacion central de la figura anterior SI se verifico contra T06 (VII sin
# C 0/21, II sin B 0/10, XI sin A 0/7); se conserva y se completa con las
# perdidas parciales que el titulo omitia. El orden de filas es el de
# 03_arquitectura.R para que los dos paneles se lean fila a fila.
if (!exists("ANCHO_COLUMNA")) source("scripts/00_estilo.R")
suppressPackageStartupMessages({library(gggenes); library(tidyr)})

UMBRAL_PARCIAL <- 0.70   # por debajo de esta fraccion se considera perdida parcial

censo <- leer_censo()
t4 <- suppressMessages(read_csv(file.path(TABLAS, "T04_motivos_MEME.csv"), show_col_types = FALSE))
t5 <- suppressMessages(read_csv(file.path(TABLAS, "T05_motivos_sitios.csv"), show_col_types = FALSE))
t6 <- suppressMessages(read_csv(file.path(TABLAS, "T06_subdominios_189.csv"), show_col_types = FALSE))

orden <- censo |>
  arrange(subfamilia_final, largo_aa) |>
  transmute(simbolo, subfamilia_final = subfam_facetas(subfamilia_final),
            largo_aa, es_focal,
            fila = factor(simbolo, levels = simbolo))

sitios <- t5 |>
  dplyr::rename(simbolo = simbolo_EglNAC) |>
  left_join(t4 |> select(motivo, subdominio), by = "motivo") |>
  mutate(fin = inicio + ancho - 1,
         subdominio = ifelse(subdominio %in% MOTIVO_CTERM_CSV, MOTIVO_CTERM, subdominio),
         subdominio = ifelse(is.na(subdominio) | subdominio == "",
                             SUBDOM_NA, subdominio)) |>
  inner_join(orden, by = "simbolo")

# --- ausencias totales y parciales, desde T06 -------------------------------
cob <- t6 |>
  pivot_longer(A:E, names_to = "subdom", values_to = "presente") |>
  group_by(subfamilia_final, subdom) |>
  summarise(n = n(), k = sum(presente), .groups = "drop") |>
  mutate(frac = k / n)
# subfam_facetas() ya renombro el grupo sin evaluar al rotulo corto de banda
tot <- cob |> filter(k == 0, subfamilia_final != ETQ_PARCIAL)
par <- cob |> filter(k > 0, frac < UMBRAL_PARCIAL, subfamilia_final != ETQ_PARCIAL)

registrar_cifra("f04_ausencias_totales",
                paste(sprintf("%s sin %s (0 de %d)", tot$subfamilia_final, tot$subdom, tot$n),
                      collapse = "; "), "T06_subdominios_189.csv")
registrar_cifra("f04_n_ausencias_totales", nrow(tot), "T06")
registrar_cifra("f04_perdidas_parciales",
                paste(sprintf("%s conserva %s en %d de %d", par$subfamilia_final,
                              par$subdom, par$k, par$n), collapse = "; "),
                sprintf("T06: presencia < %.0f %%", 100 * UMBRAL_PARCIAL))
registrar_cifra("f04_motivos_meme", nrow(t4), "T04_motivos_MEME.csv")
registrar_cifra("f04_sitios", nrow(sitios), "T05 x censo")
registrar_cifra("f04_genes_con_sitio", dplyr::n_distinct(sitios$simbolo), "T05")

ALTO_BARRA <- alto_barra(nrow(orden))
registrar_cifra("f04_alto_barra",
                sprintf("%.2f mm sobre un paso de %.2f mm",
                        as.numeric(convertHeight(ALTO_BARRA, "mm")),
                        as.numeric(convertHeight(ALTO_BARRA, "mm")) / 0.62),
                "04_motivos.R: alto_barra()",
                "identico al de F07 para que los dos paneles se lean fila a fila")
sin_cabeza <- list(arrowhead_width = unit(0, "mm"), arrowhead_height = unit(0, "mm"),
                   arrow_body_height = ALTO_BARRA)

# La barra de fondo era gris sin clave en la leyenda: el lector no podia saber
# que representa el largo de la proteina. Entra como un nivel mas de la escala.
ETQ_LARGO <- "protein length"
# T04 clasifica 8 motivos como "C-terminal / especifico de subfamilia": ese valor
# no estaba en la paleta, se dibujaba con el gris de na.value y se quedaba sin
# clave en la leyenda. El guardian de abajo impide que vuelva a pasar en silencio.
ETQ_CTERM <- MOTIVO_CTERM
PAL_F04 <- c(PAL_SUBDOM, setNames("grey45", ETQ_CTERM), setNames("grey93", ETQ_LARGO))
desconocidos <- setdiff(unique(sitios$subdominio), names(PAL_F04))
if (length(desconocidos))
  stop("clases de motivo sin color definido: ", paste(desconocidos, collapse = ", "))
# la fila del gen focal mide 2,7 pt: sin un triangulo en el margen el rotulo
# rojo no se puede asociar a ninguna fila
x_tri <- -0.022 * max(orden$largo_aa)
orden$clase_fondo <- ETQ_LARGO

p <- ggplot(orden, aes(xmin = 1, xmax = largo_aa, y = fila)) +
  # el relleno de la barra de fondo viaja como columna del data, no como
  # variable suelta: aes() se evalua al recomponer el .rds en 12_paneles.R
  do.call(geom_gene_arrow, c(list(aes(fill = clase_fondo), colour = NA), sin_cabeza)) +
  do.call(geom_subgene_arrow, c(list(data = sitios,
                                     aes(xsubmin = inicio, xsubmax = fin, fill = subdominio),
                                     colour = NA), sin_cabeza)) +
  geom_point(data = filter(orden, es_focal) |> mutate(x_tri = x_tri),
             aes(x = x_tri, y = fila),
             inherit.aes = FALSE, shape = 62, size = 1.9, colour = COL_DESTACADO) +
  geom_text(data = filter(orden, es_focal),
            aes(x = largo_aa + 12, y = fila, label = simbolo), inherit.aes = FALSE,
            # vjust: el gen focal es la primera fila de Vb; centrado, el borde del panel cortaba el rótulo
            hjust = 0, vjust = 0.9, size = 2.4, fontface = "bold", colour = COL_DESTACADO) +
  facet_grid(subfamilia_final ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_fill_manual(values = PAL_F04, name = NULL,
                    breaks = c("A", "B", "C", "D", "E", ETQ_CTERM, ETQ_LARGO),
                    labels = c("A", "B", "C", "D", "E",
                               "C-terminal (subfamily)", ETQ_LARGO)) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.06))) +
  labs(x = ETQ_EJE_AA, y = NULL) +
  # con una sola columna la septima clave se salia del alto fijo de la figura
  guides(fill = guide_legend(ncol = 2, keyheight = unit(2.6, "mm"),
                             keywidth = unit(2.6, "mm"))) +
  theme_nac(8) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        panel.grid.major.y = element_blank(), panel.spacing.y = unit(1.8, "pt"),
        strip.placement = "outside", strip.background = element_blank(),
        strip.text.y.left = element_text(angle = 0, size = 6.5, hjust = 1, face = "plain"),
        legend.text = element_text(size = 7), strip.clip = "off",
        legend.position = "bottom", legend.box.margin = margin(-3, 0, 0, 0))

guardar_objeto(p, "P_motivos")
guardar_figura(p, "F04_motivos_NAC", alto = ALTO_PANELES, ancho = ANCHO_MEDIA)

registrar_cifra("f04_motivos_cterm",
                sum(t4$subdominio %in% MOTIVO_CTERM_CSV, na.rm = TRUE),
                "T04: motivos MEME clasificados como C-terminal / especifico de subfamilia")
volcar_cifras("f04")
