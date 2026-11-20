# 02_subfamilias.R -> panel de genes por subfamilia
# F05 fue la unica figura cuyos conteos reprodujeron el censo gen por gen: se
# conserva la logica y se corrigen dos cosas que la revision marco -- los 6 genes
# de dominio parcial ahora aparecen, y el color pasa a ser el del mapa unico.
if (!exists("ANCHO_COLUMNA")) source("scripts/00_estilo.R")

censo <- leer_censo()
# Una barra por subfamilia.
d <- censo |>
  count(subfamilia_final, name = "n")
tot <- d |> group_by(subfamilia_final) |> summarise(n = sum(n), .groups = "drop") |>
  arrange(n)
# 'no evaluada (dominio parcial)' va al fondo del orden y no intercalada entre
# subfamilias por tamano: si no, la linea de separacion cae en medio de la lista
# y sugiere una agrupacion que no existe.
niveles <- c(SUBFAM_NO_EVAL,
             as.character(tot$subfamilia_final[tot$subfamilia_final != SUBFAM_NO_EVAL]))
d$subfamilia_final   <- factor(d$subfamilia_final, levels = niveles)
tot$subfamilia_final <- factor(tot$subfamilia_final, levels = niveles)

stopifnot(sum(d$n) == 189)
registrar_cifra("f05_total", sum(d$n), "T01")
registrar_cifra("f05_subfamilias", sum(tot$subfamilia_final != SUBFAM_NO_EVAL), "T01")
registrar_cifra("f05_con_subfamilia",
                sum(censo$subfamilia_final != SUBFAM_NO_EVAL), "T01")
registrar_cifra("f05_ya_manuscrito", sum(censo$en_manuscrito), "T01.en_manuscrito")
registrar_cifra("f05_nuevos", sum(!censo$en_manuscrito), "T01.en_manuscrito")
registrar_cifra("f05_mayor",
                sprintf("%s (%d)", tot$subfamilia_final[which.max(tot$n)], max(tot$n)), "T01")
registrar_cifra("f05_no_evaluadas", tot$n[tot$subfamilia_final == SUBFAM_NO_EVAL],
                "T01: dominio NAM parcial, excluidos de la clasificacion")
registrar_cifra("f05_focal_subfam",
                sprintf("%s (%d genes)", GEN_FOCAL_SUBFAM,
                        tot$n[tot$subfamilia_final == GEN_FOCAL_SUBFAM]), "T01")

p <- ggplot(d, aes(x = n, y = subfamilia_final, fill = subfamilia_final)) +
  geom_col(width = 0.74) +
  geom_text(data = tot, aes(x = n, y = subfamilia_final, label = n),
            inherit.aes = FALSE, hjust = -0.35, size = 2.5, colour = "grey25") +
  annotate("segment", x = -Inf, xend = Inf, y = 1.5, yend = 1.5,
           linetype = "dotted", colour = "grey55", linewidth = 0.3) +
  scale_fill_manual(values = PAL_SUBFAM, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.12))) +
  # el rotulo largo del grupo sin evaluar se comia un cuarto del ancho util y
  # dejaba el texto a ras del borde: en ejes y bandas va la forma corta, la
  # larga se reserva para las leyendas, donde hay sitio
  scale_y_discrete(labels = function(x) ifelse(x == SUBFAM_NO_EVAL, ETQ_PARCIAL, x)) +
  labs(x = "number of NAC genes", y = NULL) +
  theme_nac(9) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.25))

guardar_objeto(p, "P_subfamilias")
guardar_figura(p, "F05_subfamilias_NAC", alto = 3.50)

volcar_cifras("f05")
